"""Select best reference sequence."""
import logging
import re

from Bio.Data.CodonTable import IUPACData
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from cool_seq_tool.schemas import TranscriptPriority

from dcd_mapping.lookup import (
    get_chromosome_identifier,
    get_gene_symbol,
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

__all__ = ["select_transcript", "TxSelectError"]

_logger = logging.getLogger(__name__)


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


async def _get_compatible_transcripts(
    target_gene: TargetGene, align_result: AlignmentResult
) -> list[list[str]]:
    """Acquire matching transcripts

    :param metadata: metadata for scoreset
    :param align_result: output of ``align()`` method
    :return: List of list of compatible transcripts
    """
    if align_result.chrom.startswith("chr"):
        aligned_chrom = align_result.chrom[3:]
    else:
        aligned_chrom = align_result.chrom
    chromosome = get_chromosome_identifier(aligned_chrom)
    gene_symbol = get_gene_symbol(target_gene)
    if not gene_symbol:
        msg = (
            f"Unable to find gene symbol for target gene {target_gene.target_gene_name}"
        )
        raise TxSelectError(msg)
    transcript_matches = []
    for hit_range in align_result.hit_subranges:
        matches_list = await get_transcripts(
            gene_symbol, chromosome, hit_range.start, hit_range.end
        )
        if matches_list:
            transcript_matches.append(matches_list)
    return transcript_matches


def _reduce_compatible_transcripts(matching_transcripts: list[list[str]]) -> list[str]:
    """Reduce list of list of transcripts to a list containing only entries present
    in each sublist

    :param matching_transcripts: list of list of transcript accession IDs
    :return: list of transcripts shared by all sublists
    """
    common_transcripts_set = set(matching_transcripts[0])
    for sublist in matching_transcripts[1:]:
        common_transcripts_set.intersection_update(sublist)
    return list(common_transcripts_set)


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
    matching_transcripts = await _get_compatible_transcripts(target_gene, align_result)
    if not matching_transcripts:
        common_transcripts = None
    else:
        common_transcripts = _reduce_compatible_transcripts(matching_transcripts)
    if not common_transcripts:
        if not target_gene.target_uniprot_ref:
            msg = f"Unable to find matching transcripts for target gene {target_gene.target_gene_name}"
            raise TxSelectError(msg)
        protein_sequence = get_uniprot_sequence(target_gene.target_uniprot_ref.id)
        np_accession = target_gene.target_uniprot_ref.id
        ref_sequence = get_uniprot_sequence(target_gene.target_uniprot_ref.id)
        if not ref_sequence:
            msg = f"Unable to grab reference sequence from uniprot.org for target gene {target_gene.target_gene_name}"
            raise TxSelectError(msg)
        nm_accession = None
        tx_mode = None
    else:
        mane_transcripts = get_mane_transcripts(common_transcripts)
        best_tx = _choose_best_mane_transcript(mane_transcripts)
        if not best_tx:
            best_tx = await _get_longest_compatible_transcript(common_transcripts)
        if not best_tx:
            msg = f"Unable to find matching MANE transcripts for target gene {target_gene.target_gene_name}"
            raise TxSelectError(msg)
        ref_sequence = get_sequence(best_tx.refseq_prot)
        nm_accession = best_tx.refseq_nuc
        np_accession = best_tx.refseq_prot
        tx_mode = best_tx.transcript_priority

    protein_sequence = _get_protein_sequence(target_gene.target_sequence)
    is_full_match = ref_sequence.find(protein_sequence) != -1
    start = ref_sequence.find(protein_sequence[:10])

    return TxSelectResult(
        nm=nm_accession,
        np=np_accession,
        start=start,
        is_full_match=is_full_match,
        sequence=protein_sequence,
        transcript_mode=tx_mode,
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
    align_results: dict[str, AlignmentResult | None],
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
                selected_transcripts[target_gene] = e

    return selected_transcripts
