"""Align MaveDB target sequences to a human reference genome."""
import logging
import os
import subprocess
import tempfile
from collections.abc import Mapping
from pathlib import Path
from urllib.parse import urlparse

import requests
from Bio.SearchIO import HSP
from Bio.SearchIO import parse as parse_blat
from Bio.SearchIO import read as read_blat
from Bio.SearchIO._model import Hit, QueryResult
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
    AlignmentResult,
    GeneLocation,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
)

__all__ = ["align"]

_logger = logging.getLogger(__name__)


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
    url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit"
    parsed_url = urlparse(url)
    if not dcd_mapping_dir:
        dcd_mapping_dir = LOCAL_STORE_PATH
    genome_file = dcd_mapping_dir / Path(parsed_url.path).name
    # this file shouldn't change, so no need to think about more advanced caching
    if not genome_file.exists():
        try:
            http_download(url, genome_file, silent)
        except requests.HTTPError as e:
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
) -> subprocess.CompletedProcess:
    """Execute BLAT binary with relevant params.

    Currently, we rely on a system-installed BLAT binary accessible in the containing
    environment's PATH, or under env var ``BLAT_BIN_PATH``. This is sort of awkward and
    it'd be nice to make use of some direct bindings or better packaging if that's possible.

    * Perhaps `gget`? https://pachterlab.github.io/gget/en/blat.html
    * ``PxBlat``? https://github.com/ylab-hi/pxblat

    :param target_args: target params eg ``"-q=prot -t=dnax"`` (can be empty)
    :param query_file: path to query FASTA file
    :param out_file: path-like string to output fill (could be "/dev/stdout")
    :param silent: if True, suppress all console output
    :return: process result
    """
    bin_name = os.environ["BLAT_BIN_PATH"] if "BLAT_BIN_PATH" in os.environ else "blat"  # noqa: SIM401
    command = f"{bin_name} {reference_file} {target_args} -minScore=20 {query_file} {out_file}"
    _logger.debug("Running BLAT command: %s", command)
    result = subprocess.run(  # noqa: UP022
        command,
        shell=True,  # noqa: S602
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _logger.debug("BLAT command finished with result %s", result.returncode)
    if result.returncode == 127:
        raise BlatNotFoundError
    if result.returncode != 0:
        msg = f"BLAT process returned error code {result.returncode}: {target_args} {query_file} {out_file}"
        raise AlignmentError(msg)
    return result


def _write_blat_output_tempfile(result: subprocess.CompletedProcess) -> str:
    """Create temp BLAT output file. Not immediately deleted, but should eventually
    be cleared by the OS.

    :param result: BLAT process result object
    :return: path-like string representing file location
    """
    raw_output = result.stdout.split(b"Loaded")[0]
    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(raw_output)
    return tmp.name


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


def _get_blat_output(
    metadata: ScoresetMetadata, silent: bool
) -> dict[str, QueryResult]:
    """Run a BLAT query and returns a path to the output object.

    If unable to produce a valid query the first time, then try a query using ``dnax``
    bases.

    :param scoreset_metadata: object containing scoreset attributes
    :param silent: suppress BLAT command output
    :return: dict where keys are target gene identifiers and values are BLAT query result objects
    :raise AlignmentError: if BLAT subprocess returns error code
    """
    with tempfile.NamedTemporaryFile() as tmp_file:
        query_file = _build_query_file(metadata, Path(tmp_file.name))
        target_sequence_type = _get_target_sequence_type(metadata)
        if target_sequence_type == TargetSequenceType.PROTEIN:
            target_args = "-q=prot -t=dnax"
        elif target_sequence_type == TargetSequenceType.DNA:
            target_args = ""
        else:
            # TODO consider implementing support for mixed types, not hard to do - just split blat into two files and run command with each set of arguments.
            msg = "Mapping for score sets with a mix of nucleotide and protein target sequences is not currently supported."
            raise NotImplementedError(msg)
        reference_genome_file = get_ref_genome_file(silent=silent)
        process_result = _run_blat(
            target_args, query_file, reference_genome_file, "/dev/stdout", silent
        )
        out_file = _write_blat_output_tempfile(process_result)

        try:
            output = parse_blat(out_file, "blat-psl")

        except ValueError:
            target_args = "-q=dnax -t=dnax"
            process_result = _run_blat(
                target_args, query_file, reference_genome_file, "/dev/stdout", silent
            )
            out_file = _write_blat_output_tempfile(process_result)
            try:
                output = parse_blat(out_file, "blat-psl")
            except ValueError as e:
                msg = f"Unable to run successful BLAT on {metadata.urn}"
                raise AlignmentError(msg) from e

    return output


def _get_best_hit(output: QueryResult, chromosome: str | None) -> Hit:
    """Get best hit from BLAT output.

    First, try to return hit corresponding to expected chromosome taken from scoreset
    metadata. If chromosome doesn't match any of the outputs or is unavailable, take
    the hit with the single highest-scoring HSP.

    :param output: BLAT output
    :param chromosome: refseq chromosome ID, e.g. ``"NC_000001.11"``
    :return: best Hit
    :raise AlignmentError: if unable to get hits from output
    """
    if chromosome:
        if chromosome.startswith("refseq"):
            chromosome = chromosome[7:]

        for hit in output:
            hit_chr = hit.id
            if hit_chr.startswith("chr"):
                hit_chr = hit_chr[3:]
            hit_chr_ac = get_chromosome_identifier(hit_chr)
            if hit_chr_ac == chromosome:
                return hit
        else:
            if list(output):
                hit_chrs = [h.id for h in output]
                # TODO should this be an error rather than a warning? it seems like a problem if we can't find a hit on the expected chromosome
                _logger.warning(
                    "Failed to match hit chromosomes during alignment for target %s. Expected chromosome: %s, hit chromosomes: %s",
                    output.id,
                    chromosome,
                    hit_chrs,
                )

    best_score = 0
    best_score_hit = None
    for hit in output:
        best_local_score = max(hit, key=lambda i: i.score).score
        if best_local_score > best_score:
            best_score = best_local_score
            best_score_hit = hit

    if best_score_hit is None:
        msg = f"Couldn't get BLAT hits for target {output.id}."
        raise AlignmentError(msg)

    return best_score_hit


def _get_best_hsp(hit: Hit, gene_location: GeneLocation | None = None) -> HSP:
    """Retrieve preferred HSP from BLAT Hit object.

    If gene location data is available, prefer the HSP with the least distance
    between the start of the hit and the start coordinate of the gene. Otherwise,
    take the HSP with the highest score value.

    :param hit: hit object from BLAT result
    :param gene_location: location data acquired by normalizing scoreset metadata
    :return: Preferred HSP object
    :raise AlignmentError: if hit object appears to be empty (should be impossible)
    """
    best_hsp = None
    if gene_location and gene_location.start is not None:
        best_hsp = min(hit, key=lambda hsp: abs(hsp.hit_start - gene_location.start))
    else:
        best_hsp = max(hit, key=lambda hsp: hsp.score)
    if best_hsp is None:
        msg = f"Unable to get best HSP from BLAT hit: {hit}"
        raise AlignmentError(msg)
    return best_hsp


def _get_best_match(output: QueryResult, target_gene: TargetGene) -> AlignmentResult:
    """Obtain best high-scoring pairs (HSP) object for query sequence.

    :param metadata: scoreset metadata
    :param output: BLAT result object
    :return: alignment result ??
    """
    location = get_gene_location(target_gene)
    chromosome = location.chromosome if location else None
    best_hit = _get_best_hit(output, chromosome)
    best_hsp = _get_best_hsp(best_hit, location)

    strand = Strand.POSITIVE if best_hsp[0].query_strand == 1 else Strand.NEGATIVE
    coverage = 100 * (best_hsp.query_end - best_hsp.query_start) / output.seq_len
    identity = best_hsp.ident_pct
    chrom = best_hsp.hit_id

    query_subranges = []
    hit_subranges = []
    for fragment in best_hsp:
        query_subranges.append(
            SequenceRange(start=fragment.query_start, end=fragment.query_end)
        )
        hit_subranges.append(
            SequenceRange(start=fragment.hit_start, end=fragment.hit_end)
        )

    return AlignmentResult(
        chrom=chrom,
        strand=strand,
        ident_pct=identity,
        coverage=coverage,
        query_range=SequenceRange(start=best_hsp.query_start, end=best_hsp.query_end),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=best_hsp.hit_start, end=best_hsp.hit_end),
        hit_subranges=hit_subranges,
    )


def align(
    scoreset_metadata: ScoresetMetadata, silent: bool = True
) -> dict[str, AlignmentResult]:
    """Align target sequence to a reference genome.

    :param scoreset_metadata: object containing scoreset metadata
    :param silent: suppress BLAT process output if true
    :return: dictionary where keys are target gene identifiers and values are alignment result objects
    """
    blat_output = _get_blat_output(scoreset_metadata, silent)
    alignment_results = {}
    for blat_result in blat_output:
        target_label = blat_result.id
        # blat names the result id "query" if there is only one query; replace "query" with the target gene name for single-target score sets
        if target_label == "query" and len(scoreset_metadata.target_genes) == 1:
            target_label = list(scoreset_metadata.target_genes.keys())[0]  # noqa: RUF015
        # blat automatically reformats query names, so sometimes they don't match our metadata
        if target_label not in scoreset_metadata.target_genes:
            # if single-target score set, don't need to match by name
            if len(scoreset_metadata.target_genes) == 1:
                target_label = list(scoreset_metadata.target_genes.keys())[0]  # noqa: RUF015
            else:
                # try to match query name to a target gene in the metadata
                matches = 0
                for target_gene_name in scoreset_metadata.target_genes:
                    blat_target_gene_name = (
                        target_gene_name.split(" ")[0]
                        .replace("(", "")
                        .replace(")", "")
                        .replace(",", "")
                    )
                    if blat_target_gene_name == target_label:
                        target_label = target_gene_name
                        matches += 1
                # we may be missing some blat reformatting rules here - if so, this error will be thrown
                if matches == 0:
                    msg = f"BLAT result {target_label} does not match any target gene names in scoreset {scoreset_metadata.urn}."
                    raise AlignmentError(msg)
                if matches > 1:
                    # could happen if multiple target genes have the same first word in their label (unlikely)
                    msg = f"BLAT result {target_label} matches multiple target gene names in scoreset {scoreset_metadata.urn}"
        target_gene = scoreset_metadata.target_genes[target_label]
        alignment_results[target_label] = _get_best_match(blat_result, target_gene)
    # confirm that there is an alignment result for each target gene
    for target_gene in scoreset_metadata.target_genes:
        if target_gene not in alignment_results:
            msg = f"No BLAT result found for target gene {target_gene} in scoreset {scoreset_metadata.urn}"
            raise AlignmentError(msg)
    return alignment_results


def fetch_alignment(
    metadata: ScoresetMetadata, silent: bool
) -> dict[str, AlignmentResult | None]:
    alignment_results = {}
    for target_gene in metadata.target_genes:
        accession_id = metadata.target_genes[target_gene].target_accession_id
        # protein and contig/chromosome accession ids do not need to be aligned to the genome
        if accession_id.startswith(("NP", "ENSP", "NC_")):
            alignment_results[accession_id] = None
        else:
            url = f"{CDOT_URL}/transcript/{accession_id}"
            r = requests.get(url, timeout=30)

            try:
                r.raise_for_status()
            except requests.HTTPError as e:
                msg = f"Received HTTPError from {url} for scoreset {metadata.urn}"
                _logger.error(msg)
                raise ResourceAcquisitionError(msg) from e

            cdot_mapping = r.json()
            alignment_results[accession_id] = parse_cdot_mapping(cdot_mapping, silent)
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

    return AlignmentResult(
        chrom=chrom,
        strand=strand,
        query_range=query_range,
        query_subranges=query_subranges,
        hit_range=hit_range,
        hit_subranges=hit_subranges,
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
    with tempfile.NamedTemporaryFile() as query_tmp_file, tempfile.NamedTemporaryFile() as reference_tmp_file:
        _write_query_file(Path(query_tmp_file.name), [">query", target_sequence])
        _write_query_file(
            Path(reference_tmp_file.name), [">reference", reference_sequence]
        )
        target_args = "-q=prot -t=prot"

        process_result = _run_blat(
            target_args,
            query_tmp_file.name,
            reference_tmp_file.name,
            "/dev/stdout",
            silent,
        )
        out_file = _write_blat_output_tempfile(process_result)

        try:
            blat_output = read_blat(out_file, "blat-psl")

        except ValueError as e:
            msg = "Unable to run successful BLAT on target sequence against selected reference protein sequence."
            raise AlignmentError(msg) from e

    best_hit = _get_best_hit(blat_output, None)
    best_hsp = _get_best_hsp(best_hit)

    query_subranges = []
    hit_subranges = []
    for fragment in best_hsp:
        query_subranges.append(
            SequenceRange(start=fragment.query_start, end=fragment.query_end)
        )
        hit_subranges.append(
            SequenceRange(start=fragment.hit_start, end=fragment.hit_end)
        )

    return AlignmentResult(
        query_range=SequenceRange(start=best_hsp.query_start, end=best_hsp.query_end),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=best_hsp.hit_start, end=best_hsp.hit_end),
        hit_subranges=hit_subranges,
    )
