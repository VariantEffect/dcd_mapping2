"""Select best reference sequence."""
import logging
import re
from collections.abc import Mapping

from Bio import Align
from Bio.Data.CodonTable import IUPACData
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from cool_seq_tool.schemas import TranscriptPriority

from dcd_mapping.exceptions import TxSelectError
from dcd_mapping.lookup import (
    get_chromosome_identifier,
    get_mane_transcripts,
    get_protein_accession,
    get_seqrepo,
    get_sequence,
    get_transcripts,
    get_uniprot_sequence,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    ManeDescription,
    ScoreRow,
    ScoresetMetadata,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TranscriptDescription,
    TxSelectResult,
)

__all__ = ["select_transcript"]

_logger = logging.getLogger(__name__)


async def _get_compatible_transcripts(
    align_result: AlignmentResult,
) -> set[tuple[str, str]]:
    """Acquire transcripts and their HGNC symbols which overlap with all hit subranges
    of an alignment result.

    :param metadata: metadata for scoreset
    :param align_result: output of ``align()`` method
    :return: Set of compatible transcripts
    """
    aligned_chrom = (
        align_result.chrom[3:]
        if align_result.chrom.startswith("chr")
        else align_result.chrom
    )
    chromosome = get_chromosome_identifier(aligned_chrom)

    transcript_matches: set[tuple[str, str]] = set()
    for hit_range in align_result.hit_subranges:
        matches_list = await get_transcripts(chromosome, hit_range.start, hit_range.end)
        if not transcript_matches:
            transcript_matches = set(matches_list)

        transcript_matches.intersection_update(matches_list)

    return transcript_matches


def _local_alignment_identity(query: str, ref: str) -> float:
    """Compute local alignment percent identity between two protein sequences.

    Uses Smith-Waterman local alignment with BLOSUM62; gap open -10, gap extend -0.5.
    Returns alignment score, or -1000 if either sequence is empty.
    """
    if not query or not ref:
        return -1000

    aligner = Align.PairwiseAligner(
        mode="local",
        substitution_matrix=Align.substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10,
        extend_gap_score=-0.5,
    )

    try:
        alignments = aligner.align(query, ref)
    except Exception as e:
        # Do not fallback to approximate similarity; propagate failure
        msg = "Local alignment failed"
        error_message = f"{msg}: {e!s}"
        raise TxSelectError(error_message) from e

    if not alignments:
        return -1000

    return alignments[0].score


def _choose_most_similar_transcript(
    protein_sequence: str, mane_transcripts: list[TranscriptDescription]
) -> TranscriptDescription | None:
    """Choose the transcript whose protein reference is most similar to the
    provided sequence.

    Selects the highest similarity; ties keep first encountered (stable).
    """
    if not mane_transcripts:
        return None
    if len(mane_transcripts) == 1:
        return mane_transcripts[0]

    best: TranscriptDescription | None = None
    best_score = -1.0
    for tx in mane_transcripts:
        ref_seq = get_sequence(tx.refseq_prot)
        score = _local_alignment_identity(protein_sequence, ref_seq)
        if score > best_score:
            best_score = score
            best = tx

    return best


def _choose_best_mane_transcript(
    mane_transcripts: list[ManeDescription],
) -> ManeDescription | None:
    """Choose best transcript (Select > Plus Clinical) given MANE status. This was
    originally a little longer but I think all we have to worry about is grabbing based
    on MANE status.

    TODO: this shouldn't be necessary anymore, we already sort them

    :param mane_transcripts: list of MANE transcript descriptions
    :return: best transcript
    """
    if not mane_transcripts:
        return None
    for transcript in mane_transcripts:
        if transcript.transcript_priority == TranscriptPriority.MANE_SELECT:
            return transcript
    for transcript in mane_transcripts:
        if transcript.transcript_priority == TranscriptPriority.MANE_PLUS_CLINICAL:
            return transcript
    return None


async def _get_longest_compatible_transcript(
    transcripts: list[str],
) -> TranscriptDescription | None:
    """Get longest transcript from a list of compatible transcripts.

    I think there's a chance of some discord between UTA and Seqrepo and we might get
    KeyErrors here. If so, we should do a filter further up to drop any transcript
    accession IDs not recognized by SeqRepo.

    :param transcripts:
    :return:
    """
    transcripts.sort(key=lambda tx: len(get_seqrepo().sr[tx]))
    nm = transcripts[-1]
    np = await get_protein_accession(nm)
    if not np:
        return None
    return TranscriptDescription(
        refseq_nuc=nm,
        refseq_prot=np,
        transcript_priority=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )


def _get_protein_sequence(target_sequence: str) -> str:
    """Get protein sequence if necessary.

    It'd be nice if there was a more elegant way to check if the sequence was already a
    protein sequence (it should be possible for protein sequences to contain <5 unique
    bases, albeit unlikely with a large enough length).

    :param target_sequence: sequence set as baseline in MAVE experiment (might already
        be set to protein)
    :return: resulting protein sequence
    """
    if len(set(target_sequence)) > 4:
        protein_sequence = target_sequence
    else:
        protein_sequence = str(Seq(target_sequence).translate(table="1")).replace(
            "*", ""
        )
    return protein_sequence


async def _select_protein_reference(
    target_gene: TargetGene, align_result: AlignmentResult
) -> TxSelectResult:
    """Select preferred transcript for protein reference sequence

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: Best transcript and associated metadata
    :raise TxSelectError: if no matching MANE transcripts and unable to get UniProt ID/
    reference sequence
    """
    matching_transcripts = await _get_compatible_transcripts(align_result)

    # Map HGNC symbols to their compatible transcripts
    hgnc_to_transcripts: dict[str, list[str]] = {}
    for tx, hgnc in matching_transcripts:
        hgnc_to_transcripts.setdefault(hgnc, []).append(tx)

    per_gene_best: list[ManeDescription | TranscriptDescription] = []
    best_tx: ManeDescription | TranscriptDescription | None = None

    # Choose one best transcript per gene (based on MANE priority, falling back to longest)
    for _, transcripts in hgnc_to_transcripts.items():
        if not transcripts:
            continue

        mane_transcripts = get_mane_transcripts(transcripts)
        best_tx = _choose_best_mane_transcript(mane_transcripts)

        if not best_tx:
            best_tx = await _get_longest_compatible_transcript(transcripts)

        if best_tx:
            per_gene_best.append(best_tx)

    # If we found any per-gene best candidates, Step 2: choose the most similar among them and
    # select it.
    if per_gene_best:
        if not target_gene.target_sequence:
            msg = f"Unable to find target sequence for target gene {target_gene.target_gene_name}"
            raise TxSelectError(msg)

        protein_sequence = _get_protein_sequence(target_gene.target_sequence)
        best_tx = _choose_most_similar_transcript(protein_sequence, per_gene_best)

        # As a fallback, pick the first candidate
        if not best_tx:
            best_tx = per_gene_best[0]

        ref_sequence = get_sequence(best_tx.refseq_prot)
        is_full_match = ref_sequence.find(protein_sequence) != -1
        start = ref_sequence.find(protein_sequence[:10])

        return TxSelectResult(
            nm=best_tx.refseq_nuc,
            np=best_tx.refseq_prot,
            start=start,
            is_full_match=is_full_match,
            sequence=get_sequence(best_tx.refseq_prot),
            transcript_mode=best_tx.transcript_priority,
            # Only MANE transcripts have symbols
            hgnc_symbol=best_tx.symbol if hasattr(best_tx, "symbol") else None,
        )

    # If we didn't find any suitable transcript, attempt to use a provided UniProt reference
    if not target_gene.target_uniprot_ref:
        msg = f"Unable to find matching transcripts for target gene {target_gene.target_gene_name}"
        raise TxSelectError(msg)

    uniprot_sequence = get_uniprot_sequence(target_gene.target_uniprot_ref.id)
    if not uniprot_sequence:
        msg = f"Unable to grab reference sequence from uniprot.org for target gene {target_gene.target_gene_name}"
        raise TxSelectError(msg)

    is_full_match = uniprot_sequence.find(protein_sequence) != -1
    start = uniprot_sequence.find(protein_sequence[:10])

    return TxSelectResult(
        nm=None,
        np=target_gene.target_uniprot_ref.id,
        start=start,
        is_full_match=is_full_match,
        sequence=protein_sequence,
        transcript_mode=None,
        hgnc_symbol=None,
    )


def _offset_target_sequence(target_gene: TargetGene, records: list[ScoreRow]) -> int:
    """Find start location in target sequence

    :param target_gene: MaveDB metadata for target gene
    :param records: individual score records (including MAVE-HGVS descriptions)
    :return: starting index position (may be 0)
    """
    if not isinstance(records[0].hgvs_pro, str) or records[0].hgvs_pro.startswith("NP"):
        return 0
    protein_change_list = [rec.hgvs_pro.lstrip("p.") for rec in records]

    # build table of parseable amino acids by reference location on target sequence
    amino_acids_by_position = {}
    for protein_change in protein_change_list:
        if protein_change == "_sy" or protein_change == "_wt":
            continue
        if ";" in protein_change:
            protein_changes = protein_change[1:-1].split(";")
        else:
            protein_changes = [protein_change]
        for change in protein_changes:
            aa = change[:3]
            if aa == "=" or change[-3:] not in IUPACData.protein_letters_3to1:
                continue
            loc = change[3:-1] if "=" in change else change[3:-3]
            if loc not in amino_acids_by_position:
                loc = re.sub("[^0-9]", "", loc)
                if loc:
                    amino_acids_by_position[loc] = seq1(aa)

    err_locs = []
    protein_sequence = Seq(target_gene.target_sequence).translate(table="1")
    for i in range(len(protein_sequence)):
        if (
            str(i) in amino_acids_by_position
            and amino_acids_by_position[str(i)] != protein_sequence[i - 1]
        ):
            err_locs.append(i)
    if len(err_locs) == 0:
        return 0

    amino_acids_by_position = {int(k): v for k, v in amino_acids_by_position.items()}
    amino_acids_by_position = sorted(amino_acids_by_position.items())
    amino_acids_by_position = dict(amino_acids_by_position)
    p0, p1, p2, p3, p4 = list(amino_acids_by_position.keys())[0:5]

    seq = ""
    for value in amino_acids_by_position.values():
        seq += value

    protein_sequence = _get_protein_sequence(target_gene.target_sequence)
    offset = 0

    if protein_sequence in seq:
        return offset

    for i, base in enumerate(protein_sequence):
        if all(
            [
                base == amino_acids_by_position[p0],
                protein_sequence[i + p1 - p0] == amino_acids_by_position[p1],
                protein_sequence[i + p2 - p0] == amino_acids_by_position[p2],
                protein_sequence[i + p3 - p0] == amino_acids_by_position[p3],
                protein_sequence[i + p4 - p0] == amino_acids_by_position[p4],
            ]
        ):
            if i + 1 == min(amino_acids_by_position.keys()) or i + 2 == min(
                amino_acids_by_position.keys()
            ):
                offset = 0
            else:
                offset = i
            break
    return offset


def _handle_edge_cases(
    urn: str, transcript_reference: TxSelectResult
) -> TxSelectResult:
    """Handle a few edge case scoresets

    A handful of scoresets have known issues that require minor alterations of
    start position and sequence values. This method performs them if necessary and
    returns a transcript selection object.
    """
    if urn.startswith(("urn:mavedb:00000047", "urn:mavedb:00000048")):
        _logger.warning(
            "Setting transcript start = 0 -- there is discordance between actual and expected amino acid locations in experiments 47 and 48"
        )
        transcript_reference.start = 0
        transcript_reference.sequence = "M" + transcript_reference.sequence
    elif urn.startswith("urn:mavedb:00000058-a-1"):
        _logger.warning(
            "urn:mavedb:00000058-a-1 describes the starting residue as Asp2, but the starting residue is D -- manually reducing offset by 1 to reflect start of Met1."
        )
        transcript_reference.start = 670
        transcript_reference.sequence = "M" + transcript_reference.sequence
    elif urn.startswith("urn:mavedb:00000053"):
        _logger.warning(
            "Experiment 53's target sequence is missing start residue E -- manually reducing offset by 1"
        )
        transcript_reference.start = 309
        transcript_reference.sequence = "E" + transcript_reference.sequence
    return transcript_reference


async def select_transcript(
    scoreset_urn: str,
    target_gene: TargetGene,
    records: list[ScoreRow],
    align_result: AlignmentResult,
) -> TxSelectResult | None:
    """Select appropriate human reference sequence for one target in a score set.

    * Unnecessary for regulatory/other noncoding targets in score sets which report genomic
    variations.
    * For protein-coding targets, identify a matching RefSeq protein reference sequence.

    :param scoreset_urn: MaveDB URN for score set, used for hardcoding for certain score sets
    :param target_gene: Target gene metadata from MaveDB
    :param records:
    :param align_result: alignment results
    :return: Transcript description (accession ID, offset, selected sequence, etc)
    """
    if target_gene.target_gene_category != TargetType.PROTEIN_CODING:
        _logger.debug(
            "%s is regulatory/noncoding -- skipping transcript selection",
            target_gene.target_gene_name,
        )
        return None
    transcript_reference = await _select_protein_reference(target_gene, align_result)
    if (
        transcript_reference
        and target_gene.target_sequence_type == TargetSequenceType.DNA
    ):
        offset = _offset_target_sequence(target_gene, records)
        if offset:
            transcript_reference.start = offset

    return _handle_edge_cases(scoreset_urn, transcript_reference)


async def select_transcripts(
    scoreset_metadata: ScoresetMetadata,
    records: dict[str, list[ScoreRow]],
    align_results: Mapping[str, AlignmentResult | None],
) -> dict[str, TxSelectResult | Exception | None]:
    """Select appropriate human reference sequence for each target in a score set.
    :param scoreset_metadata: Metadata for score set from MaveDB API
    :param records: Variant/score records from MaveDB API
    :param align_results: Alignment results for all targets in score set.
        * Dict where keys are target labels and values are alignment result objects
    :return: dict where keys are target labels and values are objects describing selected transcript (accession ID, offset, selected sequence, etc)
    """
    selected_transcripts: dict[
        str, TxSelectResult | TxSelectError | KeyError | None
    ] = {}
    for target_gene in scoreset_metadata.target_genes:
        if scoreset_metadata.target_genes[target_gene].target_accession_id:
            # for accession-based targets, create tx select objects for protein sequence accessions only
            accession_id = scoreset_metadata.target_genes[
                target_gene
            ].target_accession_id
            # TODO create full list of possible protein accession prefixes
            if accession_id.startswith(("NP_", "ENSP_")):
                # TODO make sequence field optional instead of leaving blank here?
                selected_transcripts[target_gene] = TxSelectResult(
                    np=accession_id,
                    start=0,
                    is_full_match=True,
                    sequence="",
                    transcript_mode=None,
                )
            else:
                selected_transcripts[target_gene] = None
        else:
            try:
                selected_transcripts[target_gene] = await select_transcript(
                    scoreset_urn=scoreset_metadata.urn,
                    target_gene=scoreset_metadata.target_genes[target_gene],
                    records=records[target_gene],
                    align_result=align_results[target_gene],
                )
            except (TxSelectError, KeyError) as e:
                _logger.exception(
                    "Transcript selection failed for %s in %s: %s",
                    target_gene,
                    scoreset_metadata.urn,
                    e,
                )
                selected_transcripts[target_gene] = e

    return selected_transcripts
