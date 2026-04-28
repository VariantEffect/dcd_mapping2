"""Align MaveDB target sequences to a human reference genome."""
import functools
import logging
import os
import shlex
import subprocess
import tempfile
from collections.abc import Mapping
from pathlib import Path
from urllib.parse import urlparse

import httpx
from Bio import Align as BioAlign
from Bio.Seq import Seq, UndefinedSequenceError
from cool_seq_tool.schemas import Strand

from dcd_mapping.exceptions import (
    AlignmentError,
    BlatNotFoundError,
    ResourceAcquisitionError,
    ScoresetNotSupportedError,
)
from dcd_mapping.lookup import get_chromosome_identifier, get_gene_location
from dcd_mapping.mavedb_data import LOCAL_STORE_PATH, patch_target_sequence_type
from dcd_mapping.resource_utils import CDOT_URL, http_download
from dcd_mapping.schemas import (
    AlignmentQc,
    AlignmentResult,
    GeneLocation,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
)

__all__ = ["align"]


_logger = logging.getLogger(__name__)

# The assembly name is a derived property of the 2bit file built into the Docker image—
# if the URL is ever changed to another build, update REFERENCE_GENOME_ASSEMBLY here as
# well so the audit field in AlignmentResult.reference_assembly stays accurate.
REFERENCE_GENOME_URL = (
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit"
)
REFERENCE_GENOME_ASSEMBLY = "GRCh38"

# BLAT invocation parameters
BLAT_MIN_SCORE = 20
BLAT_OUT_FORMAT = "pslx"


@functools.lru_cache
def _get_blat_version(bin_name: str) -> str:
    """Return BLAT's self-reported version line by invoking the binary with no arguments.

    BLAT prints usage text to stdout when called without args; the first line
    typically reads: ``blat - Standalone BLAT v. 36 fast sequence search command line tool``
    We return that line verbatim rather than trying to parse a version number out of it,
    so the record is robust to BLAT changing its version string format.
    Result is cached per binary path so we only probe once per process.
    Returns ``"unknown"`` if the binary cannot be invoked or produces no output,
    so the audit field is always a string rather than null.
    """
    try:
        #
        result = subprocess.run([bin_name], capture_output=True, timeout=5)  # noqa: S603
        output = result.stdout.decode("utf-8", errors="replace")
        first_line = output.splitlines()[0].strip() if output.strip() else None
        if first_line:
            return first_line
    except (OSError, IndexError):
        _logger.debug(
            "Could not determine BLAT version for %s", bin_name, exc_info=True
        )
    return "unknown"


def _write_query_file(file: Path, lines: list[str]) -> None:
    """Write lines to query file. This method is broken out to enable easy mocking while
    testing.

    :param file: path to query file
    :param lines: list of lines to write (should be header and then sequence)
    """
    with file.open("w") as f:
        for line in lines:
            f.write(f"{line}\n")


def _build_query_file(scoreset_metadata: ScoresetMetadata, query_file: Path) -> Path:
    """Construct BLAT query file.

    :param scoreset_metadata: MaveDB scoreset metadata object
    :param query_file: path for query file
    :return: Yielded Path to constructed file. Deletes file once complete.
    """
    _logger.debug("Writing BLAT query to %s", query_file)
    lines = []
    for target_gene in scoreset_metadata.target_genes:
        lines.append(f">{target_gene}")
        lines.append(scoreset_metadata.target_genes[target_gene].target_sequence)

    _write_query_file(query_file, lines)
    return query_file


def get_ref_genome_file(
    silent: bool = True, dcd_mapping_dir: Path | None = None
) -> Path:
    """Acquire reference genome file in 2bit format from UCSC.

    :param build: genome build to acquire
    :param silent: if True, suppress console output
    :param dcd_mapping_dir: optionally declare genome file storage location
    :return: path to acquired file
    :raise ResourceAcquisitionError: if unable to acquire file.
    """
    url = REFERENCE_GENOME_URL
    parsed_url = urlparse(url)
    if not dcd_mapping_dir:
        dcd_mapping_dir = LOCAL_STORE_PATH

    # this file shouldn't change, so no need to think about more advanced caching
    genome_file = dcd_mapping_dir / Path(parsed_url.path).name
    if not genome_file.exists():
        try:
            http_download(url, genome_file, silent)
        except httpx.HTTPStatusError as e:
            msg = f"HTTPError when fetching reference genome file from {url}"
            _logger.error(msg)
            raise ResourceAcquisitionError(msg) from e

    return genome_file


def _run_blat(
    target_args: str,
    query_file: Path,
    reference_file: Path,
    out_file: str,
    silent: bool,
    min_score: int = BLAT_MIN_SCORE,
    out_format: str = BLAT_OUT_FORMAT,
) -> tuple[subprocess.CompletedProcess, dict]:
    """Execute BLAT binary with the given parameters.

    Returns both the process result and a parameters dict built from the exact
    values passed in — so the QC record is inseparably paired with what ran.

    Currently, we rely on a system-installed BLAT binary accessible in the containing
    environment's PATH, or under env var ``BLAT_BIN_PATH``. This is sort of awkward and
    it'd be nice to make use of some direct bindings or better packaging if that's possible.

    * Perhaps `gget`? https://pachterlab.github.io/gget/en/blat.html
    * ``PxBlat``? https://github.com/ylab-hi/pxblat

    :param target_args: target params eg ``"-q=prot -t=dnax"`` (can be empty)
    :param query_file: path to query FASTA file
    :param out_file: path-like string to output file (could be "/dev/stdout")
    :param silent: if True, suppress all console output
    :param min_score: value passed as -minScore to BLAT
    :param out_format: value passed as -out to BLAT
    :return: (process result, aligner_parameters dict)
    """
    bin_name = os.environ.get("BLAT_BIN_PATH", "blat")

    # -out=pslx emits PSL plus per-block aligned bases (tseqs/qseqs). Biopython's
    # psl parser surfaces those as aln.target.seq / aln.query.seq, which is
    # what we need to report real per-base mismatches and gaps on the output.
    # target_args is a space-separated string of BLAT flags (e.g. "-q=prot -t=dnax");
    # all values originate from internal constants, so shlex.split is safe here.
    cmd = [bin_name, str(reference_file)]
    if target_args:
        cmd.extend(shlex.split(target_args))

    cmd.extend(
        [f"-minScore={min_score}", f"-out={out_format}", str(query_file), out_file]
    )
    _logger.debug("Running BLAT command: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            shell=False,  # noqa: S603 -- BLAT args are all internally derived constants, not user input, so shell=False is safe here
            capture_output=True,
            timeout=600,
        )
    except subprocess.TimeoutExpired as e:
        msg = f"BLAT timed out after 600 s: {target_args} {query_file} {out_file}"
        raise AlignmentError(msg) from e
    except FileNotFoundError as e:
        raise BlatNotFoundError from e
    _logger.debug("BLAT command finished with result %s", result.returncode)

    if result.returncode == 127:
        raise BlatNotFoundError
    if result.returncode != 0:
        msg = f"BLAT process returned error code {result.returncode}: {target_args} {query_file} {out_file}"
        raise AlignmentError(msg)

    params = {
        "aligner": "blat",
        "aligner_version": _get_blat_version(bin_name),
        "min_score": min_score,
        "out_format": out_format,
        "target_args": target_args,
    }
    return result, params


def _write_blat_output_tempfile(result: subprocess.CompletedProcess) -> str:
    """Write BLAT pslx stdout to a tempfile and return its path.

    BLAT is invoked with ``-out=pslx``, so the file contains the standard 21
    PSL columns plus 2 extra columns (tseqs, qseqs) of aligned bases.
    ``Bio.Align.parse(..., "psl")`` reads this directly.

    :param result: BLAT process result object
    :return: path-like string representing file location
    """
    raw_output = result.stdout.split(b"Loaded")[0]
    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(raw_output)
    return tmp.name


def _count_query_insert_blocks(aln: "BioAlign.Alignment") -> int:
    """Count gap-open events where extra bases appear in the query (query inserts).

    Corresponds to PSL ``qNumInsert`` — the **number of gap-opening blocks**
    in the query, not the total number of inserted bases (that would be
    ``qBaseInsert``).  A single contiguous run of inserted query bases counts
    as one, regardless of how many bases it spans.
    (PSL field reference: https://genome.ucsc.edu/FAQ/FAQformat.html#format2)

    Used for BLAT-style identity and score, which penalize query inserts but
    ignore target inserts (introns).
    """
    coords = aln.coordinates
    n = 0
    for i in range(coords.shape[1] - 1):
        ts, te = int(coords[0][i]), int(coords[0][i + 1])
        qs, qe = int(coords[1][i]), int(coords[1][i + 1])
        if ts == te and qs != qe:
            n += 1

    return n


def _count_target_insert_blocks(aln: "BioAlign.Alignment") -> int:
    """Count gap-open events where the target advances but the query does not.

    Corresponds to PSL ``tNumInsert`` — the **number of gap-opening blocks**
    in the target, not the total bases skipped (that would be
    ``tBaseInsert``).  A single contiguous deleted region in the query counts
    as one block, regardless of its length.
    (PSL field reference: https://genome.ucsc.edu/FAQ/FAQformat.html#format2)

    Together with ``qNumInsert`` these make up the insert penalty in the
    canonical BLAT PSL score formula.
    """
    coords = aln.coordinates
    n = 0
    for i in range(coords.shape[1] - 1):
        ts, te = int(coords[0][i]), int(coords[0][i + 1])
        qs, qe = int(coords[1][i]), int(coords[1][i + 1])
        if qs == qe and ts != te:
            n += 1

    return n


def _blat_score(aln: "BioAlign.Alignment") -> int:
    """Compute the BLAT PSL score for a single HSP.

    Equivalent to the ``score`` column that BLAT writes to PSL output::

        score = matches - misMatches - qNumInsert - tNumInsert

    This formula implements the canonical PSL scoring described in the UCSC
    BLAT documentation (https://genome.ucsc.edu/FAQ/FAQblat.html#blat4) and
    matches the ``score`` field definition in the PSL format spec
    (https://genome.ucsc.edu/FAQ/FAQformat.html#format2):
    ``matches - misMatches - qNumInsert - tNumInsert``.

    **Block count, not base count:** ``qNumInsert`` and ``tNumInsert`` are
    gap-open event counts (one per contiguous inserted run), not total inserted
    bases (those are ``qBaseInsert`` / ``tBaseInsert`` in PSL).  The helpers
    :func:`_count_query_insert_blocks` and :func:`_count_target_insert_blocks`
    implement this correctly by counting coordinate transitions where one
    dimension is zero-width — each such transition is exactly one gap-open
    block in Biopython's alignment representation.

    This was the original basis for best-HSP selection before the BioPython
    port and correctly penalises both mismatches and gap-open events without
    double-counting individual gap bases.  Penalising target inserts is safe
    here because ``_get_best_hsp`` operates on a chromosome-filtered list
    where target gaps do not represent tolerated intron splices.

    :param aln: Bio.Align.Alignment representing a single BLAT HSP
    :return: integer BLAT score (can be negative for very noisy alignments)
    """
    counts = aln.counts()
    q_inserts = _count_query_insert_blocks(aln)
    t_inserts = _count_target_insert_blocks(aln)
    return int(counts.identities) - int(counts.mismatches) - q_inserts - t_inserts


def _blat_style_identity(
    identities: int, mismatches: int, query_insert_blocks: int
) -> float | None:
    """Compute BLAT-style percent identity: ``matches / (matches + mismatches + qNumInsert)``.

    Penalizes real query indels but ignores target gaps (intron-style splits),
    matching what BLAT itself reports for spliced cDNA-to-genome alignments.
    """
    denom = identities + mismatches + query_insert_blocks
    if denom == 0:
        _logger.warning(
            "BLAT identity denominator is zero — alignment has no matches, mismatches, or query inserts."
        )
        return None

    return 100.0 * identities / denom


def _compact_alignment_string(aln_str: str) -> str:
    """Collapse consecutive all-gap display blocks in a Biopython alignment string.

    Biopython renders every 60-bp genomic window as a three-line display group
    (target / match / query).  For cDNA→genome alignments the intronic regions
    produce long runs of ``?`` characters in the target row and ``-`` in the
    match/query rows, ballooning the string to tens of thousands of characters.

    This function replaces each run of consecutive all-gap groups (where the
    target sequence is entirely ``?``) with a single summary line::

        ... [4,560 bp gap: chr14:90400123-90404683] ...

    Groups that contain any real sequence characters are left untouched, so
    exon-intron boundaries and the flanking partial-gap rows remain visible.

    The function is a no-op on protein or purely genomic alignments that
    contain no ``?`` characters.
    """
    groups = aln_str.split("\n\n")
    output: list[str] = []

    # State for the currently-accumulating gap run.
    gap_chrom: str | None = None
    gap_start: int | None = None
    gap_end: int | None = None

    def _flush_gap() -> None:
        nonlocal gap_chrom, gap_start, gap_end
        if gap_start is not None:
            span = gap_end - gap_start  # type: ignore[operator]
            output.append(
                f"... [{span:,} bp gap: {gap_chrom}:{gap_start}-{gap_end}] ..."
            )
        gap_chrom = gap_start = gap_end = None

    for group in groups:
        lines = group.strip("\n").splitlines()
        if len(lines) != 3:
            _flush_gap()
            if group.strip():
                output.append(group)
            continue

        target_parts = lines[0].split()
        # Expect at least: <chrom> <start> <sequence>
        if len(target_parts) < 3:
            _flush_gap()
            output.append(group)
            continue

        try:
            t_start = int(target_parts[1])
        except ValueError:
            _flush_gap()
            output.append(group)
            continue

        seq = target_parts[2]

        if all(c == "?" for c in seq):
            # Pure-gap block — accumulate into the current run.
            t_end = t_start + len(seq)
            if gap_start is None:
                gap_chrom = target_parts[0]
                gap_start = t_start
                gap_end = t_end
            else:
                gap_end = t_end
        else:
            _flush_gap()
            output.append(group)

    _flush_gap()
    return "\n\n".join(output)


def _log_alignment_summary(label: str, result: "AlignmentResult") -> None:
    """Emit a single human-readable INFO log summarizing a finished alignment.

    ``label`` identifies the alignment context (target gene name, or
    ``"<gene> (protein)"`` for protein-to-protein BLAT). The pairwise
    visualization is logged at DEBUG to keep INFO compact.
    """
    qc = result.alignment_qc
    score = result.score
    next_best = result.next_best_score
    pct = result.percent_identity
    chrom = result.chrom or "N/A"
    span = (
        f"{result.hit_range.start}-{result.hit_range.end}"
        if result.hit_range is not None
        else "N/A"
    )
    mismatches = qc.mismatch_count if qc is not None else "N/A"
    gaps = qc.gap_count if qc is not None else "N/A"
    aln_len = qc.alignment_length if qc is not None else "N/A"
    cigar = qc.cigar if qc is not None else None

    _logger.info(
        "Alignment[%s]: chrom=%s span=%s strand=%s "
        "score=%s next_best=%s identity=%s length=%s mismatches=%s gaps=%s%s",
        label,
        chrom,
        span,
        result.strand.value if result.strand is not None else "N/A",
        f"{score:g}" if score is not None else "N/A",
        f"{next_best:g}" if next_best is not None else "N/A",
        f"{pct:.2f}%" if pct is not None else "N/A",
        aln_len,
        mismatches,
        gaps,
        f" cigar={cigar}" if cigar else "",
    )
    if qc is not None and qc.alignment_string:
        _logger.debug("Alignment[%s] visualization:\n%s", label, qc.alignment_string)
    else:
        _logger.debug(
            "Alignment[%s] visualization unavailable (BLAT output may have lacked aligned sequence content).",
            label,
        )


def _cigar_from_coords(aln: "BioAlign.Alignment") -> str | None:
    """Build a CIGAR string directly from alignment coordinates.

    Biopython's SAM formatter raises ``ValueError: Unequal step sizes`` for
    minus-strand PSL alignments because query coordinates are stored in
    descending order (step = -1) while target steps are +1.  Building the
    CIGAR ourselves from absolute step sizes avoids this entirely.
    """
    coords = aln.coordinates
    assert coords.shape[0] == 2, (  # noqa: S101
        f"_cigar_from_coords expects a pairwise alignment (shape[0]==2), got {coords.shape[0]}"
    )
    if coords.shape[1] < 2:
        return None

    parts: list[str] = []

    # Leading soft clip: query bases before the first aligned block.
    # min() handles both strand orientations (descending coords for minus-strand).
    q_aligned_min = int(min(coords[1][0], coords[1][-1]))
    if q_aligned_min > 0:
        parts.append(f"{q_aligned_min}S")

    for i in range(coords.shape[1] - 1):
        t_step = abs(int(coords[0][i + 1]) - int(coords[0][i]))
        q_step = abs(int(coords[1][i + 1]) - int(coords[1][i]))
        # query insert: query advances, target does not
        if t_step == 0 and q_step > 0:
            parts.append(f"{q_step}I")
        # deletion: target advances, query does not
        elif t_step > 0 and q_step == 0:
            parts.append(f"{t_step}D")
        # aligned block: both advance (match or mismatch)
        elif t_step > 0 and q_step > 0:
            parts.append(f"{t_step}M")
        # t_step == 0 and q_step == 0: degenerate zero-width column; skip

    # Trailing soft clip.
    q_aligned_max = int(max(coords[1][0], coords[1][-1]))
    q_len = len(aln.query.seq)  # works even when .seq is UndefinedSequenceData
    if q_len > q_aligned_max:
        parts.append(f"{q_len - q_aligned_max}S")

    return "".join(parts) if parts else None


def _build_alignment_qc(aln: "BioAlign.Alignment") -> "AlignmentQc":
    """Derive the per-(target, alignment_level) QC block from a Bio.Align
    ``Alignment`` parsed out of BLAT's pslx output.

    Aggregate stats come from ``Alignment.counts()``; CIGAR from
    ``_cigar_from_coords()``; the pairwise visualization from ``str(aln)``.  A single
    walk over ``.coordinates`` collects mismatch positions and gap intervals on
    the target (genome) for downstream per-variant flagging; those lists live
    in memory only and are not serialized.

    When per-base sequence content is unavailable (e.g. BLAT pslx omitted
    tseqs/qseqs for cDNA→genome runs), ``mismatch_count`` from
    ``Alignment.counts()`` remains accurate but ``mismatch_positions`` will be
    empty and ``mismatch_positions_unavailable`` will be ``True``.  Callers
    must treat ``at_mismatched_locus`` as ``None`` (not evaluated) in that
    situation to avoid a silent disagreement where ``mismatch_count > 0`` but
    no variant is ever flagged.
    """
    counts = aln.counts()
    coords = aln.coordinates
    mismatch_positions: list[int] = []
    gap_intervals: list[tuple[int, int]] = []
    query_insert_blocks = 0

    # Per-base mismatch detection requires concrete sequence content for both
    # target and query. Biopython's pslx parser sometimes leaves bases outside
    # aligned blocks (or even within, depending on BLAT output) as undefined
    # placeholders. The aggregate mismatch *count* still comes from
    # ``aln.counts()`` below; we treat per-position recording as best-effort
    # and set mismatch_positions_unavailable=True if the sequence content isn't
    # materialized so callers know at_mismatched_locus is unreliable.
    mismatch_positions_unavailable = False
    sequences_available = True
    for i in range(coords.shape[1] - 1):
        ts, te = int(coords[0][i]), int(coords[0][i + 1])
        qs, qe = int(coords[1][i]), int(coords[1][i + 1])

        # gap in reference (extra bases in query): zero-width interval
        # at ts so variants exactly at ts are "near" it within any
        # positive window. This is a query insert in PSL terms.
        if ts == te and qs != qe:
            gap_intervals.append((ts, ts))
            query_insert_blocks += 1

        # gap in query (target advances, query does not): deletion in the
        # query relative to the reference.  This is a target insert in PSL
        # terms (tNumInsert).  Record the full genomic span so variants
        # overlapping or adjacent to the deleted region can be flagged as
        # near_gap.
        elif qs == qe and ts != te:
            gap_intervals.append((min(ts, te), max(ts, te)))

        # Aligned block -- target coords are always ascending (ts < te).
        # For minus-strand query blocks (qe < qs), Biopython stores the
        # query SeqRecord in plus-strand orientation; we need the reverse
        # complement of the [qe:qs] slice to match the target bases.
        elif ts != te and qs != qe and sequences_available:
            try:
                t_bases = str(aln.target.seq[ts:te])
                if qe < qs:
                    q_bases = str(aln.query.seq[qe:qs].reverse_complement())
                else:
                    q_bases = str(aln.query.seq[qs:qe])
            except UndefinedSequenceError:
                _logger.debug(
                    "Skipping per-base mismatch positions: aligned-block "
                    "sequence content unavailable in BLAT output. "
                    "mismatch_count from counts() is still accurate; "
                    "mismatch_positions_unavailable will be set to True."
                )
                mismatch_positions = []
                mismatch_positions_unavailable = True
                sequences_available = False
                continue

            for j, (tb, qb) in enumerate(zip(t_bases, q_bases, strict=True)):
                if qb.upper() != tb.upper():
                    mismatch_positions.append(ts + j)

    cigar = _cigar_from_coords(aln)
    percent_identity = _blat_style_identity(
        counts.identities, counts.mismatches, query_insert_blocks
    )

    # PSL target coordinates are always on the + strand (ascending), so the
    # target span equals the number of reference bases consumed by M+D ops —
    # which is what alignment_length represents.  Using aln.length directly
    # raises "Unequal step sizes" for minus-strand alignments because Biopython
    # derives length from coordinate steps and rejects descending query strides.
    alignment_length = int(coords[0][-1]) - int(coords[0][0])

    return AlignmentQc(
        percent_identity=percent_identity,
        alignment_length=alignment_length,
        mismatch_count=int(counts.mismatches),
        gap_count=len(gap_intervals),
        cigar=cigar,
        alignment_string=_compact_alignment_string(str(aln)),
        mismatch_positions_unavailable=mismatch_positions_unavailable,
        mismatch_positions=mismatch_positions,
        gap_intervals=gap_intervals,
    )


def _get_target_sequence_type(metadata: ScoresetMetadata) -> TargetSequenceType | str:
    """Get overall target sequence type for a score set's target genes.
    Protein if all target sequences are protein sequences, nucleotide if all target
    sequences are nucleotide sequences, and mixed if there is a mix within the score set.
    :param metadata: object containing score set attributes
    :return: TargetSequenceType enum (protein or nucleotide) or string "mixed"
    """
    target_sequence_types = set()
    for target_gene in metadata.target_genes:
        target_sequence_types.add(
            metadata.target_genes[target_gene].target_sequence_type
        )

    if len(target_sequence_types) > 1:
        return "mixed"
    elif len(target_sequence_types) == 1:  # noqa: RET505
        return target_sequence_types.pop()
    else:
        msg = f"Target sequence types not available for score set {metadata.urn}"
        raise ValueError(msg)


def _group_alignments_by_query(
    pslx_path: str,
) -> dict[str, list["BioAlign.Alignment"]]:
    """Parse a pslx file and group Bio.Align Alignments by query ID.

    :param pslx_path: path to a pslx-format file produced by BLAT
    :return: mapping of query ID to the list of all alignments for that query
    :raise ValueError: propagated from Bio.Align.parse on malformed input
    """
    grouped: dict[str, list[BioAlign.Alignment]] = {}
    for aln in BioAlign.parse(pslx_path, "psl"):
        grouped.setdefault(str(aln.query.id), []).append(aln)

    return grouped


def _get_blat_output(
    metadata: ScoresetMetadata, silent: bool
) -> tuple[dict[str, list["BioAlign.Alignment"]], dict]:
    """Run a BLAT query and return a per-gene grouping of Bio.Align Alignments.

    If unable to produce a valid query the first time, retries with
    ``-q=dnax -t=dnax`` flags.

    :param metadata: object containing scoreset attributes
    :param silent: suppress BLAT command output
    :return: ``(grouped, blat_params)`` — ``grouped`` maps each query gene
        label to its alignments; ``blat_params`` is the parameters dict
        returned directly by ``_run_blat`` from the successful invocation.
    :raise AlignmentError: if BLAT subprocess returns error code or both
        parse attempts fail.
    """
    with tempfile.NamedTemporaryFile() as tmp_file:
        query_file = _build_query_file(metadata, Path(tmp_file.name))
        target_sequence_type = _get_target_sequence_type(metadata)

        if target_sequence_type == TargetSequenceType.PROTEIN:
            target_args = "-q=prot -t=dnax"
        elif target_sequence_type == TargetSequenceType.DNA:
            target_args = ""
        else:
            # TODO consider implementing support for mixed types
            msg = "Mapping for score sets with a mix of nucleotide and protein target sequences is not currently supported."
            raise NotImplementedError(msg)

        _logger.info(
            "Running BLAT for %s with mode: %s",
            metadata.urn,
            target_args if target_args else "default (nucleotide DNA)",
        )

        reference_genome_file = get_ref_genome_file(silent=silent)
        process_result, blat_params = _run_blat(
            target_args,
            query_file,
            reference_genome_file,
            "/dev/stdout",
            silent,
        )
        pslx_path = _write_blat_output_tempfile(process_result)

        try:
            grouped = _group_alignments_by_query(pslx_path)
            if not grouped:
                msg = f"BLAT produced no alignments for {metadata.urn} with mode: {target_args if target_args else 'default (nucleotide DNA)'}"
                raise ValueError(msg)

        except ValueError:
            _logger.debug(
                "Initial BLAT parse failed for %s, retrying with -q=dnax -t=dnax",
                metadata.urn,
                exc_info=True,
            )

            process_result, blat_params = _run_blat(
                "-q=dnax -t=dnax",
                query_file,
                reference_genome_file,
                "/dev/stdout",
                silent,
            )
            pslx_path = _write_blat_output_tempfile(process_result)

            try:
                grouped = _group_alignments_by_query(pslx_path)
                if not grouped:
                    msg = f"BLAT produced no alignments for {metadata.urn} after retry with -q=dnax -t=dnax"
                    raise ValueError(msg)

            except ValueError as e:
                msg = f"Unable to run successful BLAT on {metadata.urn}"
                _logger.exception(msg=msg, exc_info=e)
                raise AlignmentError(msg) from e

            _logger.info(
                "BLAT -q=dnax -t=dnax retry succeeded for %s (%d query gene(s)).",
                metadata.urn,
                len(grouped),
            )

    return grouped, blat_params


def _get_best_hit(
    alignments: "list[BioAlign.Alignment]",
    chromosome: str | None,
    scores: dict[int, int] | None = None,
) -> "list[BioAlign.Alignment]":
    """Return the subset of alignments on the best-matching chromosome.

    First, try to find alignments on the expected chromosome taken from scoreset
    metadata.  If no chromosome matches or none is provided, fall back to the
    chromosome whose best-scoring alignment has the highest BLAT score.

    :param alignments: all alignments for one query gene
    :param chromosome: RefSeq chromosome accession, e.g. ``"NC_000001.11"``
    :param scores: optional pre-computed ``{id(aln): score}`` cache; when
        provided avoids redundant ``aln.counts()`` walks inside the fallback.
    :return: alignments filtered to the selected chromosome
    :raise AlignmentError: if the alignment list is empty
    """
    if not alignments:
        msg = f"Empty alignment list passed to _get_best_hit for chromosome {chromosome or 'N/A'}."
        raise AlignmentError(msg)

    if chromosome:
        if chromosome.startswith("refseq:"):
            chromosome = chromosome[7:]

        for aln in alignments:
            hit_chr = str(aln.target.id)
            if hit_chr.startswith("chr"):
                hit_chr = hit_chr[3:]

            hit_chr_ac = get_chromosome_identifier(hit_chr)
            if hit_chr_ac == chromosome:
                matched_chrom = str(aln.target.id)
                return [a for a in alignments if str(a.target.id) == matched_chrom]

        # No hit matched the expected chromosome; log and fall through.
        hit_chrs = list({str(a.target.id) for a in alignments})
        query_id = str(alignments[0].query.id)

        _logger.warning(
            "Failed to match hit chromosomes during alignment for target %s. "
            "Expected chromosome: %s, hit chromosomes: %s",
            query_id,
            chromosome,
            hit_chrs,
        )

    # Fallback: find the chromosome with the highest-scoring alignment.
    chrom_best: dict[str, float] = {}
    for aln in alignments:
        chrom = str(aln.target.id)
        s = float(scores[id(aln)] if scores is not None else _blat_score(aln))
        if chrom not in chrom_best or s > chrom_best[chrom]:
            chrom_best[chrom] = s

    if not chrom_best:
        query_id = str(alignments[0].query.id) if alignments else "unknown"
        msg = f"Couldn't determine best chromosome for target {query_id} — no valid alignments found."
        raise AlignmentError(msg)

    best_chrom = max(chrom_best, key=chrom_best.__getitem__)
    query_id = str(alignments[0].query.id)

    _logger.info(
        "Chromosome fallback selected %s for target %s (best identity score: %.1f, %d candidates).",
        best_chrom,
        query_id,
        chrom_best[best_chrom],
        len(chrom_best),
    )
    return [a for a in alignments if str(a.target.id) == best_chrom]


def _get_best_hsp(
    alignments: "list[BioAlign.Alignment]",
    gene_location: GeneLocation | None = None,
    scores: dict[int, int] | None = None,
) -> "BioAlign.Alignment":
    """Retrieve the preferred alignment from a list filtered to one chromosome.

    If gene location data is available, prefer the alignment with the least
    distance between its target start and the known gene start coordinate.
    Otherwise, take the alignment with the highest BLAT score.

    :param alignments: alignments on a single chromosome for one query gene
    :param gene_location: location data acquired by normalising scoreset metadata
    :param scores: optional pre-computed ``{id(aln): score}`` cache; when
        provided avoids redundant ``aln.counts()`` walks.
    :return: preferred alignment
    :raise AlignmentError: if the list is empty
    """
    if not alignments:
        msg = f"Empty alignment list passed to _get_best_hsp for gene location {gene_location or 'N/A'}."
        raise AlignmentError(msg)

    if gene_location and gene_location.start is not None:
        best = min(
            alignments,
            key=lambda a: abs(int(a.coordinates[0][0]) - gene_location.start),
        )

        _logger.info(
            "Selected HSP by gene-location proximity (gene start: %d, hit start: %d, %d candidate(s)).",
            gene_location.start,
            int(best.coordinates[0][0]),
            len(alignments),
        )
        return best

    best = max(
        alignments,
        key=lambda a: scores[id(a)] if scores is not None else _blat_score(a),
    )

    _logger.info(
        "Selected HSP by best BLAT score (%d, %d candidate(s)).",
        scores[id(best)] if scores is not None else _blat_score(best),
        len(alignments),
    )
    return best


def _get_best_match(
    alignments: "list[BioAlign.Alignment]",
    target_gene: TargetGene,
    seq_len: int,
    blat_params: dict,
) -> AlignmentResult:
    """Obtain the best-matching alignment for one query gene.

    :param alignments: all alignments for this gene across all chromosomes
    :param target_gene: target gene metadata (used for gene-location lookup)
    :param seq_len: query sequence length (for coverage calculation)
    :param blat_params: parameters dict returned by ``_run_blat``
    :return: AlignmentResult
    """
    location = get_gene_location(target_gene)
    chromosome = location.chromosome if location else None

    # Pre-compute BLAT scores once per alignment — each call to _blat_score walks
    # the alignment coordinates, so caching avoids O(n * alignment_length) work
    # when the same alignments are ranked by multiple selection helpers.
    _scores: dict[int, int] = {id(a): _blat_score(a) for a in alignments}
    chrom_alns = _get_best_hit(alignments, chromosome, _scores)
    best_aln = _get_best_hsp(chrom_alns, location, _scores)

    best_counts = best_aln.counts()
    others = [_scores[id(a)] for a in alignments if a is not best_aln]
    next_best = float(max(others)) if others else None

    coords = best_aln.coordinates
    tcoords = coords[0]
    qcoords = coords[1]

    strand = Strand.POSITIVE if int(qcoords[0]) <= int(qcoords[-1]) else Strand.NEGATIVE
    q_start = int(qcoords.min())
    q_end = int(qcoords.max())

    coverage = 100.0 * (q_end - q_start) / seq_len if seq_len else 0.0
    identity = _blat_style_identity(
        best_counts.identities,
        best_counts.mismatches,
        _count_query_insert_blocks(best_aln),
    )
    chrom = str(best_aln.target.id)

    query_subranges: list[SequenceRange] = []
    hit_subranges: list[SequenceRange] = []
    for i in range(coords.shape[1] - 1):
        ts, te = int(tcoords[i]), int(tcoords[i + 1])
        qs, qe = int(qcoords[i]), int(qcoords[i + 1])

        # skip gap-only segments (BLAT-style identity ignores target gaps
        # and penalizes query gaps, so these don't represent real aligned
        # sequence)
        if ts == te or qs == qe:
            continue

        hit_subranges.append(SequenceRange(start=ts, end=te))
        query_subranges.append(SequenceRange(start=min(qs, qe), end=max(qs, qe)))

    alignment_qc = _build_alignment_qc(best_aln)

    return AlignmentResult(
        aligner_parameters=blat_params,
        reference_assembly=REFERENCE_GENOME_ASSEMBLY,
        chrom=chrom,
        strand=strand,
        percent_identity=identity,
        coverage=coverage,
        query_range=SequenceRange(start=q_start, end=q_end),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=int(tcoords[0]), end=int(tcoords[-1])),
        hit_subranges=hit_subranges,
        score=float(_scores[id(best_aln)]),
        next_best_score=next_best,
        alignment_qc=alignment_qc,
    )


def align(
    scoreset_metadata: ScoresetMetadata, silent: bool = True
) -> dict[str, AlignmentResult]:
    """Align target sequence to a reference genome.

    :param scoreset_metadata: object containing scoreset metadata
    :param silent: suppress BLAT process output if true
    :return: dictionary where keys are target gene identifiers and values are alignment result objects
    """
    grouped, blat_params = _get_blat_output(scoreset_metadata, silent)
    alignment_results = {}
    for blat_id, aln_list in grouped.items():
        target_label = blat_id

        # BLAT names the result "query" when there is only one query sequence;
        # replace with the actual target gene name for single-target score sets.
        if target_label == "query" and len(scoreset_metadata.target_genes) == 1:
            target_label = list(scoreset_metadata.target_genes.keys())[0]  # noqa: RUF015

        # BLAT automatically reformats query names; resolve back to metadata key.
        if target_label not in scoreset_metadata.target_genes:
            if len(scoreset_metadata.target_genes) == 1:
                target_label = list(scoreset_metadata.target_genes.keys())[0]  # noqa: RUF015
            else:
                matches = 0
                for target_gene_name in scoreset_metadata.target_genes:
                    blat_target_gene_name = (
                        target_gene_name.split(" ")[0]
                        .replace("(", "")
                        .replace(")", "")
                        .replace(",", "")
                    )

                    if blat_target_gene_name == blat_id:
                        target_label = target_gene_name
                        matches += 1

                if matches == 0:
                    msg = f"BLAT result {blat_id} does not match any target gene names in scoreset {scoreset_metadata.urn}."
                    raise AlignmentError(msg)
                if matches > 1:
                    msg = f"BLAT result {blat_id} matches multiple target gene names in scoreset {scoreset_metadata.urn}"
                    raise AlignmentError(msg)

        target_gene = scoreset_metadata.target_genes[target_label]
        seq_len = len(target_gene.target_sequence or "")
        alignment_results[target_label] = _get_best_match(
            aln_list, target_gene, seq_len, blat_params
        )
        _log_alignment_summary(target_label, alignment_results[target_label])

    # confirm that there is an alignment result for each target gene
    for tgt_gene in scoreset_metadata.target_genes:
        if tgt_gene not in alignment_results:
            msg = f"No BLAT result found for target gene {tgt_gene} in scoreset {scoreset_metadata.urn}"
            raise AlignmentError(msg)

    return alignment_results


def fetch_alignment(
    metadata: ScoresetMetadata, silent: bool
) -> dict[str, AlignmentResult | None]:
    alignment_results: dict[str, AlignmentResult | None] = {}
    for target_gene in metadata.target_genes:
        accession_id = metadata.target_genes[target_gene].target_accession_id
        if not accession_id:
            msg = f"Target gene {target_gene} in scoreset {metadata.urn} is missing an accession ID, which is required for fetching CDOT alignments."
            raise AlignmentError(msg)

        # protein and contig/chromosome accession ids do not need to be aligned to the genome
        if accession_id.startswith(("NP", "ENSP", "NC_")):
            _logger.info(
                "Skipping BLAT for %s (%s): direct protein/contig accession — alignment not required.",
                target_gene,
                accession_id,
            )
            alignment_results[accession_id] = None

        else:
            url = f"{CDOT_URL}/transcript/{accession_id}"
            r = httpx.get(url, timeout=30)

            try:
                r.raise_for_status()
            except httpx.HTTPStatusError as e:
                msg = f"Received HTTPError from {url} for scoreset {metadata.urn}"
                _logger.error(msg)
                raise ResourceAcquisitionError(msg) from e

            cdot_mapping = r.json()
            cdot_alignment = parse_cdot_mapping(cdot_mapping, silent)
            alignment_results[accession_id] = cdot_alignment

            _logger.info(
                "CDOT alignment parsed for %s (%s): %s, strand=%s, %d exon(s)",
                target_gene,
                accession_id,
                cdot_alignment.chrom,
                cdot_alignment.strand,
                len(cdot_alignment.hit_subranges),
            )

    return alignment_results


def parse_cdot_mapping(cdot_mapping: dict, silent: bool) -> AlignmentResult:
    # blat psl & AlignmentResult: 0-based, start inclusive, stop exclusive
    # cdot: 1-based, start inclusive, stop inclusive
    # so, to "translate" cdot ranges to AlignmentResult-style ranges:
    # subtract 1 from start and end to go from 1-based to 0-based coord,
    # and then add 1 to the stop to go from inclusive to exclusive
    # so just subtract 1 from start and do nothing to end

    grch38 = cdot_mapping.get("genome_builds", {}).get("GRCh38")
    grch37 = cdot_mapping.get("genome_builds", {}).get("GRCh37")
    mapping = grch38 if grch38 else grch37
    if mapping is None:
        msg = f"Cdot transcript results for transcript {cdot_mapping.get('id')} do not include GRCh37 or GRCh38 mapping"
        raise AlignmentError(msg)

    chrom = mapping["contig"]
    strand = Strand.POSITIVE if mapping["strand"] == "+" else Strand.NEGATIVE
    query_subranges = []
    hit_subranges = []
    for exon in mapping["exons"]:
        query_subranges.append(SequenceRange(start=exon[3] - 1, end=exon[4]))
        hit_subranges.append(SequenceRange(start=exon[0] - 1, end=exon[1]))

    if strand == Strand.POSITIVE:
        query_range = SequenceRange(
            start=query_subranges[0].start, end=query_subranges[-1].end
        )
        hit_range = SequenceRange(
            start=hit_subranges[0].start, end=hit_subranges[-1].end
        )
    else:
        query_range = SequenceRange(
            start=query_subranges[-1].start, end=query_subranges[0].end
        )
        hit_range = SequenceRange(
            start=hit_subranges[-1].start, end=hit_subranges[0].end
        )

    reference_assembly = "GRCh38" if grch38 else "GRCh37"
    cdot_data_version = cdot_mapping.get("cdot_data_version")
    # Always emit the dict so consumers can distinguish
    # "cdot run, version unknown" (cdot_data_version=None) from "not a cdot run" (None).
    aligner_parameters: dict = {"cdot_data_version": cdot_data_version}

    return AlignmentResult(
        chrom=chrom,
        strand=strand,
        query_range=query_range,
        query_subranges=query_subranges,
        hit_range=hit_range,
        hit_subranges=hit_subranges,
        aligner_parameters=aligner_parameters,
        reference_assembly=reference_assembly,
    )


def build_alignment_result(
    metadata: ScoresetMetadata, silent: bool
) -> Mapping[str, AlignmentResult | None]:
    # NOTE: Score set must contain all accession-based target genes or all sequence-based target genes
    # This decision was made because it is most efficient to run BLAT all together, so the alignment function
    # works on an entire score set rather than per target gene.
    # However, if the need arises, we can allow both types of target genes in a score set.

    # determine whether score set is accession-based or sequence-based
    score_set_type = None
    for target_gene in metadata.target_genes:
        if metadata.target_genes[target_gene].target_accession_id:
            if score_set_type == "sequence":
                msg = "Score set contains both accession-based and sequence-based target genes. This is not currently supported."
                raise ScoresetNotSupportedError(msg)
            score_set_type = "accession"
        else:
            if score_set_type == "accession":
                msg = "Score set contains both accession-based and sequence-based target genes. This is not currently supported."
                raise ScoresetNotSupportedError(msg)
            score_set_type = "sequence"

    _logger.info(
        "Building alignment result for %s (path: %s, targets: %s)",
        metadata.urn,
        score_set_type,
        list(metadata.target_genes.keys()),
    )

    if score_set_type == "sequence":
        try:
            alignment_result = align(metadata, silent)
        except AlignmentError as e:
            failed_at_nucleotide_level = any(
                target_gene.target_sequence_type == TargetSequenceType.DNA
                for target_gene in metadata.target_genes.values()
            )

            if failed_at_nucleotide_level:
                msg = f"BLAT alignment failed for {metadata.urn} at the nucleotide level. This alignment will be retried at the protein level."
                _logger.warning(msg)
            else:
                raise AlignmentError from e

            # So long as force=True, the content of the records dict is irrelevant.
            try:
                alignment_result = align(
                    patch_target_sequence_type(metadata, {}, force=True), silent
                )
                _logger.info(
                    "Protein-level BLAT retry succeeded for %s. "
                    "The nucleotide-level alignment result has been discarded and replaced with the protein-level result.",
                    metadata.urn,
                )

            except AlignmentError as e2:
                msg = f"BLAT alignment failed for {metadata.urn} at the protein level after failing at the nucleotide level."
                _logger.error(msg)
                raise AlignmentError(msg) from e2

    else:
        alignment_result = fetch_alignment(metadata, silent)

    return alignment_result


def align_target_to_protein(
    target_sequence: str, reference_sequence: str, silent: bool
) -> AlignmentResult:
    """Align a protein target sequence against a reference protein using BLAT.

    :param target_sequence: query protein sequence
    :param reference_sequence: reference protein sequence
    :param silent: suppress BLAT process output if true
    :return: AlignmentResult for the best alignment
    :raise AlignmentError: if BLAT produces no usable alignment
    """
    target_args = "-q=prot -t=prot"
    with tempfile.NamedTemporaryFile() as query_tmp_file, tempfile.NamedTemporaryFile() as reference_tmp_file:
        _write_query_file(Path(query_tmp_file.name), [">query", target_sequence])
        _write_query_file(
            Path(reference_tmp_file.name), [">reference", reference_sequence]
        )

        process_result, blat_params = _run_blat(
            target_args,
            query_tmp_file.name,
            reference_tmp_file.name,
            "/dev/stdout",
            silent,
        )
        pslx_path = _write_blat_output_tempfile(process_result)

        try:
            alignments = list(BioAlign.parse(pslx_path, "psl"))
            if not alignments:
                raise ValueError

        except ValueError as e:
            msg = "Unable to run successful BLAT on target sequence against selected reference protein sequence."
            _logger.error(msg)
            raise AlignmentError(msg) from e

    _scores: dict[int, int] = {id(a): _blat_score(a) for a in alignments}
    best_aln = max(alignments, key=lambda a: _scores[id(a)])
    best_counts = best_aln.counts()
    others = [_scores[id(a)] for a in alignments if a is not best_aln]
    next_best = float(max(others)) if others else None

    coords = best_aln.coordinates
    tcoords = coords[0]
    qcoords = coords[1]

    query_subranges: list[SequenceRange] = []
    hit_subranges: list[SequenceRange] = []
    for i in range(coords.shape[1] - 1):
        ts, te = int(tcoords[i]), int(tcoords[i + 1])
        qs, qe = int(qcoords[i]), int(qcoords[i + 1])
        if ts == te or qs == qe:
            continue
        hit_subranges.append(SequenceRange(start=ts, end=te))
        query_subranges.append(SequenceRange(start=min(qs, qe), end=max(qs, qe)))

    # Attach full sequences so _build_alignment_qc can do per-base mismatch
    # detection. BLAT protein pslx does not populate tseqs/qseqs columns, so
    # the parser leaves SeqRecords with undefined content; we supply the same
    # strings we passed to BLAT, which share coordinate space with the PSL output.
    best_aln.target.seq = Seq(reference_sequence)
    best_aln.query.seq = Seq(target_sequence)
    alignment_qc = _build_alignment_qc(best_aln)

    result = AlignmentResult(
        query_range=SequenceRange(start=int(qcoords.min()), end=int(qcoords.max())),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=int(tcoords[0]), end=int(tcoords[-1])),
        hit_subranges=hit_subranges,
        percent_identity=_blat_style_identity(
            best_counts.identities,
            best_counts.mismatches,
            _count_query_insert_blocks(best_aln),
        ),
        score=float(_scores[id(best_aln)]),
        next_best_score=next_best,
        alignment_qc=alignment_qc,
        aligner_parameters=blat_params,
    )

    _log_alignment_summary("protein->protein", result)
    return result
