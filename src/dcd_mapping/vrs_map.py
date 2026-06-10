"""Map transcripts to VRS objects."""

import logging
import os
from collections import Counter
from collections.abc import Iterable
from enum import StrEnum
from itertools import cycle

from Bio.Seq import Seq
from bioutils.accessions import infer_namespace
from cool_seq_tool.schemas import AnnotationLayer, Strand
from ga4gh.core import sha512t24u
from ga4gh.vrs._internal.models import (
    Allele,
    Expression,
    Haplotype,
    LiteralSequenceExpression,
    ReferenceLengthExpression,
    SequenceLocation,
    SequenceString,
    Syntax,
)
from ga4gh.vrs.normalize import normalize
from mavehgvs.util import parse_variant_strings
from mavehgvs.variant import Variant

from dcd_mapping.align import align_target_to_protein
from dcd_mapping.exceptions import (
    MissingSequenceIdError,
    UnsupportedReferenceSequenceNameSpaceError,
    UnsupportedReferenceSequencePrefixError,
)
from dcd_mapping.lookup import (
    build_ref_identical_allele,
    cdot_rest,
    coding_hgvs_is_intronic,
    get_chromosome_identifier,
    get_genomic_accession_for_transcript,
    get_seqrepo,
    project_coding_hgvs_to_genomic,
    project_coding_hgvs_to_protein,
    project_genomic_hgvs_to_coding,
    translate_hgvs_to_vrs,
    translate_ref_identical_to_vrs,
)
from dcd_mapping.resource_utils import is_missing_value, request_with_backoff
from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    MappingOutcome,
    ScoreRow,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
    VrsMapResult,
)
from dcd_mapping.transcripts import TxSelectError
from dcd_mapping.vrs_utils import identify_allele, normalize_and_identify

__all__ = ["vrs_map"]


_logger = logging.getLogger(__name__)

CLINGEN_API_URL = os.environ.get("CLINGEN_API_URL", "https://reg.genome.network/allele")


class ProjectionOutcome(StrEnum):
    """Per-variant projection QC (internal; drives the validation log only, never emitted).

    ``failed`` means the selected transcript could not project its own variant -- a
    mis-selection signal.
    """

    PROJECTED = "projected"  # coding (and usually protein) form constructed cleanly
    INTRONIC = "intronic"  # coding projection lands in an intron -- benign, no protein consequence
    SKIPPED = "skipped"  # not a projectable variant row (_wt/_sy/=/fs)
    FAILED = "failed"  # g.->c. projection or coding allele construction failed -- selection suspect


_PROJECTION_FAILURE_WARN_FRACTION = 0.25
"""Projection-failure fraction above which the per-target summary logs at WARNING, not INFO."""


def _hgvs_variant_is_valid(hgvs_string: str) -> bool:
    return not hgvs_string.endswith((".=", ")", "X"))


def _process_any_aa_code(hgvs_pro_string: str) -> str:
    """Substitute "Xaa" for the "?" wildcard in a protein expression.

    Nucleotide strings with "X"-style wildcards are not adjusted -- they're treated as
    invalid (see ``_hgvs_variant_is_valid``).

    :param hgvs_string: MAVE HGVS expression
    :return: processed variation (equivalent to input if no wildcard code found)
    """
    if "?" in hgvs_pro_string:
        _logger.debug("Substituting Xaa for ? in %s", hgvs_pro_string)
        hgvs_pro_string = hgvs_pro_string.replace("?", "Xaa")
    return hgvs_pro_string


def is_intronic_variant(variant: Variant) -> bool:
    """Return True if given Variant is intronic, otherwise return False.
    Supports single or multi-position variants.
    """
    # ref identical variants with no positions should not be considered intronic
    if variant.positions is None:
        return False

    if isinstance(variant.positions, Iterable):
        if any(position.is_intronic() for position in variant.positions):
            return True
    else:
        if variant.positions.is_intronic():
            return True

    return False


def fetch_clingen_genomic_hgvs(hgvs: str) -> str | None:
    """Fetch the genomic HGVS string from ClinGen.

    :param hgvs: The HGVS string to fetch
    :return: The genomic HGVS string on GRCh38, or None if not found
    """
    if CLINGEN_API_URL is None:
        msg = "CLINGEN_API_URL environment variable is not set and default is unavailable."
        _logger.error(msg)
        raise ValueError(msg)
    response = request_with_backoff(url=f"{CLINGEN_API_URL}?hgvs={hgvs}", timeout=30)
    if response.status_code == 200:
        data = response.json()
        for allele in data.get("genomicAlleles", []):
            if allele.get("referenceGenome") == "GRCh38":
                for coordinates in allele.get("hgvs", []):
                    if coordinates.startswith("NC_"):
                        return coordinates
    return None


def _create_pre_mapped_hgvs_strings(
    raw_description: str,
    layer: AnnotationLayer,
    tx: TxSelectResult | None = None,
    alignment: AlignmentResult | None = None,
    accession_id: str | None = None,
) -> list[str]:
    """Generate pre-mapped HGVS strings from a raw string of one or more HGVS substrings.

    Known limitation: the transcript is used as the reference, but pre-mapped variants
    should be relative to the user-provided target sequence; any transcript/target offset
    is not accounted for here.

    :param raw_description: A string containing valid HGVS sub-strings
    :param layer: An enum denoting the targeted annotation layer of these HGVS strings
    :param tx: A TxSelectResult object defining the transcript we are mapping to (or None).
    :param alignment: An AlignmentResult object defining the alignment we are mapping to (or None).
    :param accession_id: An accession id describing the reference sequence (for accession-based target gene variants)
    :return: A list of HGVS strings prior to being mapped to the `tx` or `alignment`
    """
    if layer is AnnotationLayer.PROTEIN and tx is None:
        msg = f"Transcript result must be provided for {layer} annotations (Transcript was `{tx}`)."
        raise ValueError(msg)
    if layer is AnnotationLayer.GENOMIC and alignment is None and accession_id is None:
        msg = f"Alignment result or accession ID must be provided for {layer} annotations (Alignment was `{alignment}`)."
        raise ValueError(msg)

    raw_variant_strings = _parse_raw_variant_str(raw_description)
    variants, errors = parse_variant_strings(raw_variant_strings)

    hgvs_strings = []
    for variant, error in zip(variants, errors, strict=True):
        if error is not None:
            msg = f"Variant could not be parsed by mavehgvs: {error}"
            raise ValueError(msg)

        # ga4gh hgvs_tools can't translate intronic variants -- reject them here.
        if is_intronic_variant(variant):
            msg = f"Variant is intronic and cannot be processed: {variant}"
            raise ValueError(msg)

        if accession_id:
            hgvs_strings.append(accession_id + ":" + str(variant))
        # Ideally we would create an HGVS string namespaced to GA4GH. The line below
        # creates such a string, but it is not able to be parsed by the GA4GH VRS translator.
        # hgvs_strings.append('ga4gh:' + sequence_id + ':' + str(variant))
        elif layer is AnnotationLayer.PROTEIN:
            assert tx  # noqa: S101  # mypy help
            hgvs_strings.append(tx.np + ":" + str(variant))
        elif layer is AnnotationLayer.GENOMIC:
            assert alignment  # noqa: S101  # mypy help
            hgvs_strings.append(
                get_chromosome_identifier(alignment.chrom) + ":" + str(variant)
            )
        else:
            msg = (
                f"Could not generate HGVS strings for invalid AnnotationLayer: {layer}"
            )
            raise ValueError(msg)

    return hgvs_strings


def _create_post_mapped_hgvs_strings(
    raw_description: str,
    layer: AnnotationLayer,
    tx: TxSelectResult | None = None,
    alignment: AlignmentResult | None = None,
    accession_id: str | None = None,
    protein_alignment: AlignmentResult | None = None,
) -> list[str]:
    """Generate a list of (post-mapped) HGVS strings from one long string containing many valid HGVS substrings.

    For protein annotations, these strings must be adjusted to match the alignment of the target to the selected reference protein sequence.
    For genomic annotations, these strings must be adjusted to match the coordinates of the reference alignment.

    :param raw_description: A string containing valid HGVS sub-strings
    :param layer: An enum denoting the targeted annotation layer of these HGVS strings
    :param tx: A TxSelectResult object defining the transcript we are mapping to (or None)
    :param alignment: An AlignmentResult object defining the alignment we are mapping to (or None)
    :param accession_id: An accession id describing the reference sequence (or None). Only used for accession-based variants.
    :return: A list of HGVS strings relative to the `tx` or `alignment`
    """
    if layer is AnnotationLayer.PROTEIN and tx is None:
        msg = f"Transcript result must be provided for {layer} annotations (Transcript was `{tx}`)."
        raise ValueError(msg)
    if layer is AnnotationLayer.GENOMIC and alignment is None and accession_id is None:
        msg = f"Alignment result or accession ID must be provided for {layer} annotations (Alignment was `{alignment}` and Accession ID was `{accession_id}`)."
        raise ValueError(msg)

    raw_variants = _parse_raw_variant_str(raw_description)
    variants, errors = parse_variant_strings(raw_variants)

    hgvs_strings = []
    for variant, error in zip(variants, errors, strict=True):
        if error is not None:
            msg = f"Variant could not be parsed by mavehgvs: {error}"
            raise ValueError(msg)

        # ga4gh hgvs_tools can't translate intronic variants -- reject them here.
        if is_intronic_variant(variant):
            msg = f"Variant is intronic and cannot be processed: {variant}"
            raise ValueError(msg)

        # Reference-identical variants require no positional adjustment
        if str(variant).endswith(".="):
            if layer is AnnotationLayer.PROTEIN:
                assert tx  # noqa: S101  # mypy help
                hgvs_strings.append(tx.np + ":" + str(variant))

            elif layer is AnnotationLayer.GENOMIC:
                if accession_id:
                    hgvs_strings.append(accession_id + ":" + str(variant))
                else:
                    assert alignment  # noqa: S101  # mypy help
                    hgvs_strings.append(
                        get_chromosome_identifier(alignment.chrom) + ":" + str(variant)
                    )

            else:
                msg = f"Could not generate HGVS strings for invalid AnnotationLayer: {layer}"
                raise ValueError(msg)

            continue

        if layer is AnnotationLayer.PROTEIN:
            assert tx  # noqa: S101  # mypy help

            variant = _adjust_protein_variant_to_ref(variant, protein_alignment)
            hgvs_strings.append(tx.np + ":" + str(variant))
        elif layer is AnnotationLayer.GENOMIC:
            if accession_id:
                pre_mapped_hgvs = accession_id + ":" + str(variant)
                # use ClinGen to align accession-based variants to genomic reference
                genomic_hgvs = fetch_clingen_genomic_hgvs(pre_mapped_hgvs)
                if genomic_hgvs:
                    hgvs_strings.append(genomic_hgvs)
                else:
                    msg = f"Could not fetch genomic HGVS on GRCh38 for accession-based variant: {pre_mapped_hgvs}"
                    raise ValueError(msg)
            else:
                assert alignment  # noqa: S101  # mypy help

                variant = _adjust_genomic_variant_to_ref(variant, alignment)
                hgvs_strings.append(
                    get_chromosome_identifier(alignment.chrom) + ":" + str(variant)
                )
        else:
            msg = (
                f"Could not generate HGVS strings for invalid AnnotationLayer: {layer}"
            )
            raise ValueError(msg)

    return hgvs_strings


def _adjust_protein_variant_to_ref(
    variant: Variant,
    alignment: AlignmentResult,
) -> Variant:
    # adjust starts - hgvs uses 1-based numbering for c. sequences, while blat hits are 0-based
    starts = []
    if isinstance(variant.positions, Iterable):
        is_multi_position = True
        for position in variant.positions:
            starts.append(position.position - 1)
    else:
        is_multi_position = False
        starts.append(variant.positions.position - 1)

    # get hit
    query_subrange_containing_hit = None
    target_subrange_containing_hit = None
    for query_subrange, target_subrange in zip(
        alignment.query_subranges, alignment.hit_subranges, strict=True
    ):
        if all(
            start >= query_subrange.start and start < query_subrange.end
            for start in starts
        ):
            query_subrange_containing_hit = query_subrange
            target_subrange_containing_hit = target_subrange
            break

    if query_subrange_containing_hit is None or target_subrange_containing_hit is None:
        msg = f"Cannot process variant {variant} because it is not fully contained within the aligned portion of the query sequence compared to the selected reference sequence"
        raise ValueError(msg)

    for idx, start in enumerate(starts):
        # get variant start relative to the reference (the "hit")
        # distance from beginning of query to variant start position:
        query_to_start = start - query_subrange_containing_hit.start

        # distance from beginning of ref to the variant start position:
        ref_to_start = target_subrange_containing_hit.start + query_to_start

        # add distance from ref to variant start; hgvs is 1-based, so convert back to 1-based
        if is_multi_position:
            variant.positions[idx].position = ref_to_start + 1
        else:
            variant.positions.position = ref_to_start + 1

    return variant


def _adjust_genomic_variant_to_ref(
    variant: Variant,
    alignment: AlignmentResult,
) -> Variant:
    """Adjust a variant relative to a provided alignment.

    :param variant: A variant object relative to a scoreset's target sequence
    :param alignment: An AlignmentResult object denoting the alignment we are mapping to
    :return: A variant object that describes the variant relative to the provided alignment result
    """
    # adjust starts - hgvs uses 1-based numbering for c. sequences, while blat hits are 0-based
    starts = []
    if isinstance(variant.positions, Iterable):
        is_multi_position = True
        for position in variant.positions:
            starts.append(position.position - 1)
    else:
        is_multi_position = False
        starts.append(variant.positions.position - 1)

    # get hit
    query_subrange_containing_hit = None
    target_subrange_containing_hit = None
    for query_subrange, target_subrange in zip(
        alignment.query_subranges, alignment.hit_subranges, strict=True
    ):
        if all(
            start >= query_subrange.start and start < query_subrange.end
            for start in starts
        ):
            query_subrange_containing_hit = query_subrange
            target_subrange_containing_hit = target_subrange
            break

    if query_subrange_containing_hit is None or target_subrange_containing_hit is None:
        msg = f"Cannot process variant {variant} because it is not fully contained within the aligned portion of the query sequence"
        raise ValueError(msg)

    for idx, start in enumerate(starts):
        if alignment.strand is Strand.POSITIVE:
            # get variant start relative to the reference (the "hit")
            # distance from beginning of query to variant start position:
            query_to_start = start - query_subrange_containing_hit.start

            # distance from beginning of ref to the variant start position:
            ref_to_start = target_subrange_containing_hit.start + query_to_start
        else:
            # picture the rev comp of the query/variant as mapping to the positive strand of the ref
            # the start of the reverse complement of the variant is the end of the "original" variant
            # so we need to know where the end of the original variant is, relative to the query molecule
            end = start

            # subtract 1 from end of hit range, because blat ranges are 0-based [start, end)
            ref_to_start = (target_subrange_containing_hit.end - 1) - (
                end - query_subrange_containing_hit.start
            )

        # add distance from ref to variant start; hgvs is 1-based, so convert back to 1-based
        if is_multi_position:
            variant.positions[idx].position = ref_to_start + 1
        else:
            variant.positions.position = ref_to_start + 1

    # get reverse complement of sequence if the target maps to the negative strand of the reference
    if alignment.strand is Strand.NEGATIVE:
        # variant._sequences can be a string or an iterable
        if isinstance(variant._sequences, str):
            variant._sequences = str(Seq(variant._sequences).reverse_complement())
        elif variant._sequences is not None:
            revcomp_sequences_list = []
            for sequence in variant._sequences:
                revcomp_sequences_list.append(str(Seq(sequence).reverse_complement()))
            variant._sequences = revcomp_sequences_list

        # reverse order of positions tuple
        if is_multi_position:
            variant._positions = tuple(reversed(list(variant.positions)))

    # change prefix from c. to g. since variant is now relative to chr reference
    variant._prefix = "g"

    return variant


def _parse_raw_variant_str(raw_description: str) -> list[str]:
    """Parse a string which may contain many HGVS strings into a list of each one.

    :param raw_description: A string that may contain a list of variant descriptions or a single variant description
    :return: A list of HGVS strings
    """
    # some variant strings follow mavehgvs format, meaning they don't have a reference sequence id and colon preceding the c./g./n./p. prefix
    # the reference sequence information has previously been parsed for score sets with multiple targets,
    # so can discard the reference sequence id and colon if they are present
    # TODO check assumption of no colon unless reference sequence identifier is supplied!
    if ":" in raw_description:
        raw_description = raw_description.split(":")[1]
    if "[" in raw_description:
        prefix = raw_description[0:2]
        return [prefix + var for var in set(raw_description[3:-1].split(";"))]

    return [raw_description]


def _map_protein_coding_pro(
    row: ScoreRow,
    sequence_id: str,
    transcript: TxSelectResult,
    protein_align_result: AlignmentResult,
) -> MappedScore:
    """Construct VRS object mapping for ``hgvs_pro`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param sequence_id: The GA4GH accession for the provided sequence
    :param transcript: The transcript selection information for a target
    :param protein_align_result: The protein-protein alignment result for a target against its selected protein reference
    :return: VRS mapping object if mapping succeeds
    """
    if row.hgvs_pro in {"_wt", "_sy"} or is_missing_value(row.hgvs_pro):
        _logger.warning(
            "Can't process variant syntax %s for %s", row.hgvs_pro, row.accession
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            error_message=f"Can't process variant syntax {row.hgvs_pro}",
        )

    if "fs" in row.hgvs_pro:
        _logger.warning(
            "Can't process variant syntax %s for %s because protein frameshift variants are not supported",
            row.hgvs_pro,
            row.accession,
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            error_message="Protein frameshift variants are not supported",
        )

    # Special case for experiment set urn:mavedb:0000097
    if row.hgvs_pro.startswith("NP_009225.1:p."):
        vrs_variation = translate_hgvs_to_vrs(row.hgvs_pro)
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            alignment_level=AnnotationLayer.PROTEIN,
            pre_mapped=vrs_variation,
            post_mapped=vrs_variation,
        )

    try:
        pre_mapped_hgvs_strings = _create_pre_mapped_hgvs_strings(
            row.hgvs_pro,
            AnnotationLayer.PROTEIN,
            tx=transcript,
        )
        pre_mapped_protein = _construct_vrs_allele(
            pre_mapped_hgvs_strings,
            AnnotationLayer.PROTEIN,
            sequence_id,
            True,
        )
    except Exception as e:
        _logger.warning(
            "An error occurred while generating pre-mapped protein variant for %s, accession %s: %s",
            row.hgvs_pro,
            row.accession,
            e,
            exc_info=True,
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            error_message=f"{type(e).__name__}: {e}",
        )

    try:
        post_mapped_hgvs_strings = _create_post_mapped_hgvs_strings(
            row.hgvs_pro,
            AnnotationLayer.PROTEIN,
            tx=transcript,
            protein_alignment=protein_align_result,
        )
        post_mapped_protein = _construct_vrs_allele(
            post_mapped_hgvs_strings,
            AnnotationLayer.PROTEIN,
            None,
            False,
        )
    except Exception as e:
        _logger.warning(
            "An error occurred while generating post-mapped protein variant for %s, accession %s: %s",
            row.hgvs_pro,
            row.accession,
            e,
            exc_info=True,
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            alignment_level=AnnotationLayer.PROTEIN,
            pre_mapped=pre_mapped_protein,
            error_message=f"{type(e).__name__}: {e}",
        )

    return MappedScore(
        accession_id=row.accession,
        score=row.score,
        alignment_level=AnnotationLayer.PROTEIN,
        pre_mapped=pre_mapped_protein,
        post_mapped=post_mapped_protein,
    )


def _map_genomic(
    row: ScoreRow,
    sequence_id: str,
    align_result: AlignmentResult | None,
) -> MappedScore:
    """Construct VRS object mapping for ``hgvs_nt`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param sequence_id: The GA4GH accession for the provided sequence
    :param align_result: The transcript selection information for a score set
    :return: VRS mapping object if mapping succeeds
    """
    namespace = infer_namespace(sequence_id)
    if namespace is None:
        if sequence_id.startswith("SQ."):
            # if the sequence id starts with SQ, it is a target sequence which is in the ga4gh namespace
            namespace = "ga4gh"
        else:
            return MappedScore(
                accession_id=row.accession,
                score=row.score,
                error_message=f"Namespace could not be inferred from sequence: {sequence_id}",
            )

    if (
        row.hgvs_nt in {"_wt", "_sy", "="}
        or "fs"
        in row.hgvs_nt  # TODO I think this line can be taken out, I don't think fs nomenclature can be used for nt anyway
    ):
        _logger.warning(
            "Can't process variant syntax %s for %s", row.hgvs_nt, row.accession
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            error_message=f"Can't process variant syntax {row.hgvs_nt}",
        )

    if align_result is None:
        # for contig accession based score sets, no mapping is performed,
        # so pre- and post-mapped alleles are the same
        try:
            pre_mapped_hgvs_strings = post_mapped_hgvs_strings = (
                _create_pre_mapped_hgvs_strings(
                    row.hgvs_nt,
                    AnnotationLayer.GENOMIC,
                    accession_id=sequence_id,
                )
            )
            # accession-based pre-mapped alleles should be constructed like post-mapped alleles (sequence id is gathered from hgvs string rather than manually provided)
            pre_mapped_genomic = _construct_vrs_allele(
                pre_mapped_hgvs_strings,
                AnnotationLayer.GENOMIC,
                None,
                False,
            )
            post_mapped_genomic = _construct_vrs_allele(
                post_mapped_hgvs_strings,
                AnnotationLayer.GENOMIC,
                None,
                False,
            )
        except Exception as e:
            _logger.warning(
                "An error occurred while generating genomic variant for %s, accession %s: %s",
                row.hgvs_nt,
                row.accession,
                e,
                exc_info=True,
            )
            return MappedScore(
                accession_id=row.accession,
                score=row.score,
                error_message=f"{type(e).__name__}: {e}",
            )

    elif namespace.lower() in ("refseq", "ncbi", "ensembl"):
        # nm/enst way
        try:
            pre_mapped_hgvs_strings = _create_pre_mapped_hgvs_strings(
                row.hgvs_nt,
                AnnotationLayer.GENOMIC,
                accession_id=sequence_id,
            )
            # accession-based pre-mapped alleles should be constructed like post-mapped alleles (sequence id is gathered from hgvs string rather than manually provided)
            pre_mapped_genomic = _construct_vrs_allele(
                pre_mapped_hgvs_strings,
                AnnotationLayer.GENOMIC,
                None,
                False,
            )
        except Exception as e:
            _logger.warning(
                "An error occurred while generating pre-mapped genomic variant for %s, accession %s: %s",
                row.hgvs_nt,
                row.accession,
                e,
                exc_info=True,
            )
            return MappedScore(
                accession_id=row.accession,
                score=row.score,
                error_message=f"{type(e).__name__}: {e}",
            )
        try:
            post_mapped_hgvs_strings = _create_post_mapped_hgvs_strings(
                row.hgvs_nt,
                AnnotationLayer.GENOMIC,
                accession_id=sequence_id,
            )
            post_mapped_genomic = _construct_vrs_allele(
                post_mapped_hgvs_strings,
                AnnotationLayer.GENOMIC,
                None,
                False,
            )
        except Exception as e:
            _logger.warning(
                "An error occurred while generating post-mapped genomic variant for %s, accession %s: %s",
                row.hgvs_nt,
                row.accession,
                e,
                exc_info=True,
            )
            return MappedScore(
                accession_id=row.accession,
                score=row.score,
                alignment_level=AnnotationLayer.GENOMIC,
                pre_mapped=pre_mapped_genomic,
                error_message=f"{type(e).__name__}: {e}",
            )
    elif namespace.lower() == "ga4gh":
        # target seq way
        try:
            pre_mapped_hgvs_strings = _create_pre_mapped_hgvs_strings(
                row.hgvs_nt,
                AnnotationLayer.GENOMIC,
                alignment=align_result,
            )
            pre_mapped_genomic = _construct_vrs_allele(
                pre_mapped_hgvs_strings,
                AnnotationLayer.GENOMIC,
                sequence_id,
                True,
            )
        except Exception as e:
            _logger.warning(
                "An error occurred while generating pre-mapped genomic variant for %s, accession %s: %s",
                row.hgvs_nt,
                row.accession,
                e,
                exc_info=True,
            )
            return MappedScore(
                accession_id=row.accession,
                score=row.score,
                error_message=f"{type(e).__name__}: {e}",
            )

        try:
            post_mapped_hgvs_strings = _create_post_mapped_hgvs_strings(
                row.hgvs_nt,
                AnnotationLayer.GENOMIC,
                alignment=align_result,
            )
            post_mapped_genomic = _construct_vrs_allele(
                post_mapped_hgvs_strings,
                AnnotationLayer.GENOMIC,
                None,
                False,
            )
        except Exception as e:
            _logger.warning(
                "An error occurred while generating post-mapped genomic variant for %s, accession %s: %s",
                row.hgvs_nt,
                row.accession,
                e,
                exc_info=True,
            )
            return MappedScore(
                accession_id=row.accession,
                score=row.score,
                alignment_level=AnnotationLayer.GENOMIC,
                pre_mapped=pre_mapped_genomic,
                error_message=f"{type(e).__name__}: {e}",
            )
    else:
        msg = f"Unsupported reference sequence namespace: {namespace}"
        raise UnsupportedReferenceSequenceNameSpaceError(msg)

    return MappedScore(
        accession_id=row.accession,
        score=row.score,
        alignment_level=AnnotationLayer.GENOMIC,
        pre_mapped=pre_mapped_genomic,
        post_mapped=post_mapped_genomic,
    )


def _map_coding(row: ScoreRow, sequence_id: str) -> MappedScore:
    """Map a cdna-source (``NM_``/``ENST``) measured variant on its native transcript.

    The measured ``c.`` variant is already the least-derived layer, and the target *is* the
    reference transcript, so pre- and post-mapped are the same coding allele -- no genomic
    round-trip. Genomic/protein forms come separately from :func:`_construct_projected_layers`.

    :param row: a MaveDB score row
    :param sequence_id: the coding transcript accession the variant is described against
    :return: a CDNA-layer mapping, or a failed ``MappedScore`` carrying the error
    """
    if row.hgvs_nt in {"_wt", "_sy", "="} or "fs" in row.hgvs_nt:
        _logger.warning(
            "Can't process variant syntax %s for %s", row.hgvs_nt, row.accession
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            error_message=f"Can't process variant syntax {row.hgvs_nt}",
        )

    try:
        coding_hgvs_strings = _create_pre_mapped_hgvs_strings(
            row.hgvs_nt, AnnotationLayer.CDNA, accession_id=sequence_id
        )
        coding_allele = _construct_vrs_allele(
            coding_hgvs_strings, AnnotationLayer.CDNA, None, False
        )
    except Exception as e:
        _logger.warning(
            "An error occurred while generating coding variant for %s, accession %s: %s",
            row.hgvs_nt,
            row.accession,
            e,
            exc_info=True,
        )
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            error_message=f"{type(e).__name__}: {e}",
        )

    return MappedScore(
        accession_id=row.accession,
        score=row.score,
        alignment_level=AnnotationLayer.CDNA,
        pre_mapped=coding_allele,
        post_mapped=coding_allele,
    )


def _get_allele_sequence(allele: Allele) -> str:
    """Get sequence for Allele

    :param allele: VRS allele
    :return: sequence
    :raise ValueError: if sequence is none
    """
    dp = get_seqrepo()
    start = allele.location.start
    end = allele.location.end
    sequence = dp.get_sequence(
        f"ga4gh:{allele.location.sequenceReference.refgetAccession}", start, end
    )
    if sequence is None:
        raise ValueError
    return sequence


def store_sequence(sequence: str) -> str:
    """Store sequence in SeqRepo.

    :param sequence: raw sequence (ie nucleotides or amino acids)
    :return: sequence ID (sans prefix, which is ``"ga4gh"``)
    """
    sequence_id = f"SQ.{sha512t24u(sequence.encode('ascii'))}"
    alias_dict_list = [{"namespace": "ga4gh", "alias": sequence_id}]
    sr = get_seqrepo()
    sr.sr.store(sequence, alias_dict_list)
    return sequence_id


def _hgvs_nt_is_valid(hgvs_nt: str) -> bool:
    """Check for invalid or unavailable nucleotide MAVE-HGVS variation

    :param hgvs_nt: MAVE_HGVS nucleotide expression
    :return: True if expression appears populated and valid
    """
    return (not is_missing_value(hgvs_nt)) and (hgvs_nt not in {"_wt", "_sy", "="})


def _hgvs_pro_is_valid(hgvs_pro: str) -> bool:
    """Check for invalid or unavailable protein MAVE-HGVS variation

    :param hgvs_nt: MAVE_HGVS protein expression
    :return: True if expression appears populated and valid
    """
    return (
        (hgvs_pro not in {"_wt", "_sy"})
        and (not is_missing_value(hgvs_pro))
        and ("fs" not in hgvs_pro)
    )


def _map_protein_coding(
    metadata: TargetGene,
    records: list[ScoreRow],
    transcript: TxSelectResult | TxSelectError,
    align_result: AlignmentResult,
    silent: bool,
) -> tuple[list[MappedScore], AlignmentResult | None]:
    """Perform mapping on protein coding experiment results.

    :param metadata: Target gene metadata from MaveDB API
    :param records: The list of MAVE variants in a given score set
    :param transcript: The transcript data for a score set, or an error message describing why an expected transcript is missing
    :param align_result: The alignment data for a score set
    :return: ``(mappings, protein_align_result)``. The target-protein-to-reference-protein
        alignment is surfaced so downstream annotation can flag protein-layer
        variants at mismatched/near-gap loci using its QC block; pipelines that
        don't need it can ignore the second element.
    """
    if metadata.target_sequence_type == TargetSequenceType.DNA:
        sequence = str(
            Seq(metadata.target_sequence).translate(table="1", stop_symbol="")
        )
        psequence_id = store_sequence(sequence)
        gsequence_id = store_sequence(metadata.target_sequence)
    else:
        sequence = metadata.target_sequence
        psequence_id = gsequence_id = store_sequence(sequence)

    # align protein sequence to selected reference protein sequence to get offsets and gaps
    protein_align_result: AlignmentResult | None = None
    if isinstance(transcript, TxSelectResult):
        protein_align_result = align_target_to_protein(
            sequence, transcript.sequence, silent
        )
        _logger.info(
            "Protein-to-protein alignment produced for %s (ref protein: %s). "
            "This alignment will be used for pro-layer variants in place of the genomic alignment.",
            metadata.target_gene_name,
            transcript.np,
        )
    else:
        _logger.warning(
            "Skipping protein-to-protein alignment for %s: no valid transcript selected (%s). "
            "Pro-layer variants will not be mapped.",
            metadata.target_gene_name,
            transcript,
        )

    # Sequence-based coding target: measures a genomic variant, carries a BLAT-selected
    # coding transcript. Project onto the unmeasured layers (c., and p. if no protein was
    # measured); build_scoreset_mapping routes them via preferred_layer_only.
    project_nm = transcript.nm if isinstance(transcript, TxSelectResult) else None
    projection_outcomes: list[ProjectionOutcome] = []

    variations: list[MappedScore] = []
    for row in records:
        hgvs_nt_mappings = None
        hgvs_pro_mappings = None
        projected: list[MappedScore] = []
        if _hgvs_nt_is_valid(row.hgvs_nt):
            hgvs_nt_mappings = _map_genomic(row, gsequence_id, align_result)
            if project_nm and align_result is not None:
                projected, outcome = _construct_projected_layers(
                    row,
                    AnnotationLayer.GENOMIC,
                    project_nm,
                    alignment=align_result,
                    # Protein is measured here only if a valid hgvs_pro is present; when
                    # it is, _map_protein_coding_pro maps it directly, so do not project
                    # a redundant protein layer.
                    project_protein=not _hgvs_pro_is_valid(row.hgvs_pro),
                )
                projection_outcomes.append(outcome)

        if (
            isinstance(transcript, TxSelectError) and not hgvs_nt_mappings
        ):  # only create error message if there is not an hgvs nt mapping
            # TODO create pre mapped allele
            hgvs_pro_mappings = MappedScore(
                accession_id=row.accession,
                score=row.score,
                error_message=str(transcript).strip("'"),
            )
        else:
            if _hgvs_pro_is_valid(row.hgvs_pro):
                if protein_align_result is not None:
                    hgvs_pro_mappings = _map_protein_coding_pro(
                        row, psequence_id, transcript, protein_align_result
                    )
                # Skip this error when an nt mapping exists: protein alignment is then
                # expected to fail, so the message would be redundant.
                elif protein_align_result is None and not hgvs_nt_mappings:
                    hgvs_pro_mappings = MappedScore(
                        accession_id=row.accession,
                        score=row.score,
                        error_message="Could not perform mapping for protein variant because transcript sequence is missing or could not be aligned to reference sequence",
                    )
            elif (
                not hgvs_nt_mappings
            ):  # only create error message if there is not an hgvs nt mapping
                hgvs_pro_mappings = MappedScore(
                    accession_id=row.accession,
                    score=row.score,
                    error_message="Invalid protein variant syntax",
                )

        # append both pro and nt mappings if both available, plus the deterministic
        # projected layers (suppressed as variants by the API, emitted by the CLI).
        if hgvs_pro_mappings:
            variations.append(hgvs_pro_mappings)
        if hgvs_nt_mappings:
            variations.append(hgvs_nt_mappings)
        variations.extend(projected)

    if project_nm and projection_outcomes:
        _log_projection_validation(
            metadata.target_accession_id or metadata.target_gene_name,
            project_nm,
            projection_outcomes,
        )

    return variations, protein_align_result


def _map_regulatory_noncoding(
    metadata: TargetGene,
    records: list[ScoreRow],
    align_result: AlignmentResult,
) -> list[MappedScore]:
    """Perform mapping on noncoding/regulatory experiment results

    :param metadata: Target gene metadata from MaveDB API
    :param records: list of MAVE experiment result rows
    :param align_result: An AlignmentResult object for a score set
    :return: A list of VRS mappings
    """
    variations: list[MappedScore] = []
    sequence_id = store_sequence(metadata.target_sequence)

    for row in records:
        hgvs_nt_mappings = _map_genomic(row, sequence_id, align_result)
        variations.append(hgvs_nt_mappings)

    return variations


def store_accession(
    accession_id: str,
) -> None:
    namespace = infer_namespace(accession_id)
    alias_dict_list = [{"namespace": namespace, "alias": accession_id}]
    cd = cdot_rest()
    sequence = cd.get_seq(accession_id)
    sr = get_seqrepo()
    sr.sr.store(sequence, alias_dict_list)


def _coding_pivot_hgvs_strings(
    row: ScoreRow,
    source_layer: AnnotationLayer,
    transcript_nm: str,
    accession_id: str | None,
    alignment: AlignmentResult | None,
) -> list[str]:
    """Build the coding (``c.``) HGVS form(s) of a measured variant -- the pivot every
    deterministic projection runs through.

    The coding form is reached differently per assay level:

    * **genomic source** (``NC_`` accession or sequence-based BLAT alignment) -- build
      the genomic HGVS first (pre-mapped on the accession, or post-mapped onto the
      reference contig via the alignment), then project ``g. -> c.`` onto the transcript.
    * **cdna source** (``NM_``/``ENST`` accession) -- the measured variant is already
      coding on ``transcript_nm``; prefix the accession onto the parsed variant string(s).

    :raise Exception: if the genomic build or ``g. -> c.`` projection fails (the caller
        classifies this as ``FAILED``).
    """
    if source_layer is AnnotationLayer.CDNA:
        # Measured variant is already coding on transcript_nm; just prefix the accession.
        return [f"{accession_id}:{v}" for v in _parse_raw_variant_str(row.hgvs_nt)]

    if accession_id is not None:
        genomic_hgvs_strings = _create_pre_mapped_hgvs_strings(
            row.hgvs_nt, AnnotationLayer.GENOMIC, accession_id=accession_id
        )
    else:
        genomic_hgvs_strings = _create_post_mapped_hgvs_strings(
            row.hgvs_nt, AnnotationLayer.GENOMIC, alignment=alignment
        )
    return [
        project_genomic_hgvs_to_coding(g, transcript_nm) for g in genomic_hgvs_strings
    ]


def _projection_record(
    row: ScoreRow,
    layer: AnnotationLayer,
    *,
    outcome: MappingOutcome,
    allele: Allele | Haplotype | None = None,
    error_message: str | None = None,
) -> MappedScore:
    """Build one projected-layer record carrying its typed :class:`MappingOutcome`.

    Benign outcomes (``INTRONIC`` / ``NO_PROTEIN_CONSEQUENCE``) leave ``error_message``
    ``None`` so a populated ``error_message`` always means a genuine failure.
    """
    return MappedScore(
        accession_id=row.accession,
        score=row.score,
        alignment_level=layer,
        pre_mapped=allele,
        post_mapped=allele,
        error_message=error_message,
        outcome=outcome,
    )


def _construct_projected_layers(
    row: ScoreRow,
    source_layer: AnnotationLayer,
    transcript_nm: str,
    *,
    accession_id: str | None = None,
    alignment: AlignmentResult | None = None,
    genomic_accession: str | None = None,
    project_protein: bool = True,
) -> tuple[list[MappedScore], ProjectionOutcome]:
    """Project a measured variant onto its own deterministic non-measured forms.

    A 1-to-1 transform of one measured variant onto the other layers -- distinct from the
    equivalence-class expansion (synonymous codons) the reverse-translation job owns. The
    non-measured forms depend on ``source_layer``:

    * **genomic source** (``NC_`` accession, or sequence-based via ``alignment``) -- cdna
      (``g. -> c.``) and, unless ``project_protein`` is False, protein (``c. -> p.``).
    * **cdna source** (``NM_``/``ENST``) -- genomic (``c. -> g.`` onto ``genomic_accession``,
      from :func:`lookup.get_genomic_accession_for_transcript`) and protein (``c. -> p.``).

    Every expected level yields exactly one record (never a silent omission), each with a
    typed :class:`MappingOutcome`: ``MAPPED``, benign ``INTRONIC``/``NO_PROTEIN_CONSEQUENCE``
    (no ``error_message``), or ``FAILED`` (``error_message`` populated). Non-variant rows
    (``_wt``/``_sy``/``=``/``fs``) yield none. The returned :class:`ProjectionOutcome` (for
    the per-target validation log) tracks the load-bearing nucleotide form -- ``FAILED``
    there is a mis-selection signal; a protein-stage miss never fails the aggregate.

    :param row: a MaveDB score row
    :param source_layer: the assay (measured) level -- ``GENOMIC`` or ``CDNA``
    :param transcript_nm: the coding transcript the cdna form is expressed against
    :param accession_id: accession the measured variant is described against; ``None`` for
        a sequence-based genomic source, which uses ``alignment``
    :param alignment: BLAT alignment for a sequence-based genomic source
    :param genomic_accession: reference contig for the ``c. -> g.`` projection (cdna source)
    :param project_protein: construct the protein form (False when protein was measured)
    :return: one typed ``MappedScore`` per expected non-measured level, and the outcome
    """
    if row.hgvs_nt in {"_wt", "_sy", "="} or "fs" in row.hgvs_nt:
        # Not an input variant (wildtype/synonymous/reference/frameshift) -- the measured
        # path already emits the row's record; there is nothing to project.
        return [], ProjectionOutcome.SKIPPED

    # The deterministic non-measured layer set, with the load-bearing nucleotide
    # re-expression first: cdna for a genomic source, genomic for a cdna source.
    load_bearing_layer = (
        AnnotationLayer.CDNA
        if source_layer is AnnotationLayer.GENOMIC
        else AnnotationLayer.GENOMIC
    )
    expected_layers = [load_bearing_layer]
    if project_protein:
        expected_layers.append(AnnotationLayer.PROTEIN)

    try:
        coding_hgvs_strings = _coding_pivot_hgvs_strings(
            row, source_layer, transcript_nm, accession_id, alignment
        )
    except Exception as e:
        # The coding pivot underpins every projection; if it cannot be built, none of the
        # expected layers can. Emit a FAILED record for each rather than omitting them.
        _logger.debug(
            "Projection for %s (accession %s) failed: coding pivot failed: %s",
            row.hgvs_nt,
            row.accession,
            e,
        )
        records = [
            _projection_record(
                row,
                layer,
                outcome=MappingOutcome.FAILED,
                error_message=f"coding pivot ({source_layer} -> c.) failed: {e}",
            )
            for layer in expected_layers
        ]
        return records, ProjectionOutcome.FAILED

    if any(coding_hgvs_is_intronic(c) for c in coding_hgvs_strings):
        # Intronic: the coding form is not VRS-representable and there is no protein
        # consequence. A benign absence for every expected layer (no error_message).
        _logger.debug(
            "Projection for %s (accession %s): coding projection is intronic.",
            row.hgvs_nt,
            row.accession,
        )
        records = [
            _projection_record(row, layer, outcome=MappingOutcome.INTRONIC)
            for layer in expected_layers
        ]
        return records, ProjectionOutcome.INTRONIC

    # 1. Resolve this layer's HGVS strings. The cdna form is the pivot itself; the
    #    genomic form is the c.->g. re-expression (cdna source); the protein form is
    #    the c.->p. consequence. A protein string-projection miss is a *benign*
    #    no-consequence; a load-bearing nucleotide miss is a FAILED record.
    records: list[MappedScore] = []
    aggregate = ProjectionOutcome.PROJECTED
    for layer in expected_layers:
        if layer is AnnotationLayer.CDNA:
            hgvs_strings = coding_hgvs_strings

        elif layer is AnnotationLayer.GENOMIC:
            if genomic_accession is None:
                records.append(
                    _projection_record(
                        row,
                        layer,
                        outcome=MappingOutcome.FAILED,
                        error_message="no resolvable reference contig for c. -> g. projection",
                    )
                )
                aggregate = ProjectionOutcome.FAILED
                continue

            try:
                hgvs_strings = [
                    project_coding_hgvs_to_genomic(c, genomic_accession)
                    for c in coding_hgvs_strings
                ]
            except Exception as e:
                records.append(
                    _projection_record(
                        row,
                        layer,
                        outcome=MappingOutcome.FAILED,
                        error_message=f"c. -> g. projection failed: {e}",
                    )
                )
                aggregate = ProjectionOutcome.FAILED
                continue

        elif layer is AnnotationLayer.PROTEIN:
            try:
                hgvs_strings = [
                    project_coding_hgvs_to_protein(c) for c in coding_hgvs_strings
                ]
            except Exception as e:
                # Non-projectable protein is a benign no-consequence (e.g. UTR); the
                # load-bearing nucleotide layer carries any mis-selection signal.
                _logger.debug(
                    "No protein consequence for %s (accession %s): %s",
                    row.hgvs_nt,
                    row.accession,
                    e,
                )
                records.append(
                    _projection_record(
                        row, layer, outcome=MappingOutcome.NO_PROTEIN_CONSEQUENCE
                    )
                )
                continue

        else:
            msg = f"Unexpected annotation layer: {layer}"
            raise ValueError(msg)

        # 2. Construct the VRS allele for this layer.
        try:
            allele = _construct_vrs_allele(hgvs_strings, layer, None, False)
            records.append(
                _projection_record(
                    row, layer, outcome=MappingOutcome.MAPPED, allele=allele
                )
            )
        except Exception as e:
            records.append(
                _projection_record(
                    row,
                    layer,
                    outcome=MappingOutcome.FAILED,
                    error_message=f"{layer} allele construction failed: {e}",
                )
            )
            if layer is load_bearing_layer:
                aggregate = ProjectionOutcome.FAILED

    return records, aggregate


def _log_projection_validation(
    accession_id: str,
    transcript_nm: str,
    outcomes: list[ProjectionOutcome],
) -> None:
    """Log a per-target summary of how cleanly a target's variants projected.

    For a selected transcript (genomic/sequence source) a high failure fraction signals a
    mis-selection; for a cdna source (``accession_id == transcript_nm``) nothing was
    selected, so this just reports the c. -> g./p. projection. WARNING above
    ``_PROJECTION_FAILURE_WARN_FRACTION`` of attempted (non-skipped) variants, else INFO.
    """
    counts = Counter(outcomes)
    projected = counts[ProjectionOutcome.PROJECTED]
    intronic = counts[ProjectionOutcome.INTRONIC]
    failed = counts[ProjectionOutcome.FAILED]
    skipped = counts[ProjectionOutcome.SKIPPED]
    attempted = projected + intronic + failed  # non-variant skipped rows excluded
    if attempted == 0:
        return

    failure_fraction = failed / attempted
    log = (
        _logger.warning
        if failure_fraction > _PROJECTION_FAILURE_WARN_FRACTION
        else _logger.info
    )

    if accession_id == transcript_nm:
        scope = f"Projecting cdna target {accession_id} to its genomic/protein forms"
    else:
        scope = f"Projection of {accession_id} onto selected transcript {transcript_nm}"

    log(
        "%s: %d projected, %d intronic, %d failed, %d skipped (%.0f%% of %d attempted failed).",
        scope,
        projected,
        intronic,
        failed,
        skipped,
        100 * failure_fraction,
        attempted,
    )


def _map_accession(
    metadata: TargetGene,
    records: list[ScoreRow],
    align_result: AlignmentResult
    | None,  # NP and NC accessions won't have alignment results
    transcript: TxSelectResult | None,
) -> list[MappedScore]:
    variations: list[MappedScore] = []
    sequence_id = metadata.target_accession_id
    if sequence_id is None:
        msg = " No target_accession_id was provided by target gene metadata. Target gene metadata must have a target_accession_id to map to VRS."
        raise MissingSequenceIdError(msg)

    store_accession(sequence_id)

    if metadata.target_accession_id.startswith(
        (
            "NP",
            "ENSP",
        )
    ):  # TODO full list of protein accession id prefixes
        for row in records:
            hgvs_pro_mappings = _map_protein_coding_pro(
                row,
                sequence_id,
                transcript,
            )
            variations.append(hgvs_pro_mappings)

    # Coding accession targets project their measured variant onto the unmeasured
    # deterministic layers; build_scoreset_mapping routes them via preferred_layer_only.
    #   NC_ (genomic source): measured variant is genomic; project g.->c.->p. onto the
    #     gene's selected MANE transcript (transcript.nm).
    #   NM_/ENST (cdna source): measured variant is already coding; project c.->g. (onto
    #     the transcript's reference contig) and c.->p.
    elif metadata.target_accession_id.startswith(
        (
            "NM",
            "ENST",
            "NC",
        )
    ):  # TODO full list of transcript and contig accession id prefixes
        is_genomic_source = metadata.target_accession_id.startswith("NC")
        project_nm = transcript.nm if isinstance(transcript, TxSelectResult) else None

        # NM_/ENST accession is itself the coding transcript even without a
        # TxSelectResult (none is produced for transcript accessions).
        if not is_genomic_source and project_nm is None:
            project_nm = sequence_id

        # The c. -> g. projection (cdna source only) needs the transcript's reference
        # contig; resolve it once per target, not per variant.
        genomic_accession = (
            get_genomic_accession_for_transcript(project_nm)
            if project_nm and not is_genomic_source
            else None
        )

        # Map the measured variant in its native layer: genomic for an NC_ source,
        # coding (directly on the transcript) for an NM_/ENST source.
        projection_outcomes: list[ProjectionOutcome] = []
        for row in records:
            if is_genomic_source:
                hgvs_nt_mappings = _map_genomic(row, sequence_id, align_result)
            else:
                hgvs_nt_mappings = _map_coding(row, sequence_id)
            variations.append(hgvs_nt_mappings)

            # cdna source: map a measured hgvs_pro directly rather than projecting it --
            # it's already in reference coordinates (pre_mapped == post_mapped). Needs the
            # NP_ from a TxSelectResult; if unavailable, skip and log rather than fall back
            # to a projection that could silently diverge from the measured value.
            has_measured_protein = False
            if not is_genomic_source and _hgvs_pro_is_valid(row.hgvs_pro):
                has_measured_protein = True
                np_accession = (
                    transcript.np
                    if isinstance(transcript, TxSelectResult) and transcript.np
                    else None
                )
                if np_accession:
                    try:
                        pro_hgvs = _create_pre_mapped_hgvs_strings(
                            row.hgvs_pro, AnnotationLayer.PROTEIN, tx=transcript
                        )
                        pro_allele = _construct_vrs_allele(
                            pro_hgvs, AnnotationLayer.PROTEIN, None, False
                        )
                        variations.append(
                            MappedScore(
                                accession_id=row.accession,
                                score=row.score,
                                alignment_level=AnnotationLayer.PROTEIN,
                                pre_mapped=pro_allele,
                                post_mapped=pro_allele,
                            )
                        )
                    except Exception as e:
                        _logger.warning(
                            "Could not map measured hgvs_pro %s for %s: %s",
                            row.hgvs_pro,
                            row.accession,
                            e,
                        )
                        variations.append(
                            MappedScore(
                                accession_id=row.accession,
                                score=row.score,
                                alignment_level=AnnotationLayer.PROTEIN,
                                error_message=f"measured protein mapping failed: {e}",
                            )
                        )
                else:
                    _logger.warning(
                        "hgvs_pro %s present for %s but no NP_ accession resolvable "
                        "from transcript; skipping protein layer",
                        row.hgvs_pro,
                        row.accession,
                    )

            if project_nm:
                projected, outcome = _construct_projected_layers(
                    row,
                    AnnotationLayer.GENOMIC
                    if is_genomic_source
                    else AnnotationLayer.CDNA,
                    project_nm,
                    accession_id=sequence_id,
                    genomic_accession=genomic_accession,
                    project_protein=not has_measured_protein,
                )
                variations.extend(projected)
                projection_outcomes.append(outcome)
        if project_nm:
            _log_projection_validation(sequence_id, project_nm, projection_outcomes)

    else:
        msg = f"Unrecognized accession prefix for accession id: {metadata.target_accession_id}"
        raise UnsupportedReferenceSequencePrefixError(msg)

    return variations


def _rle_to_lse(
    rle: ReferenceLengthExpression, location: SequenceLocation
) -> LiteralSequenceExpression:
    """Coerce ReferenceLengthExpression to LiteralSequenceExpression.

    RLEs are helpful for long repeating sequences but a) unnecessary here and b)
    create incompatibilities with some data extraction further down so to simplify,
    we'll just turn them into equivalent LiteralSequenceExpressions.
    """
    sr = get_seqrepo()
    sequence_id = location.sequenceReference.refgetAccession
    start: int = location.start
    end = start + rle.repeatSubunitLength
    subsequence = sr.get_sequence(f"ga4gh:{sequence_id}", start, end)
    c = cycle(subsequence)
    derived_sequence = ""
    for _ in range(rle.length):
        derived_sequence += next(c)
    return LiteralSequenceExpression(sequence=derived_sequence)


def _construct_vrs_allele(
    hgvs_strings: list[str],
    layer: AnnotationLayer,
    sequence_id: str | None,
    pre_map: bool,
) -> Allele | Haplotype:
    alleles: list[Allele] = []
    for hgvs_string in hgvs_strings:
        _logger.debug("Processing HGVS string: %s", hgvs_string)

        # Reference-identical variants must be built as Alleles with a
        # ReferenceLengthExpression state; ga4gh hgvs_tools can't translate them from HGVS.
        if hgvs_string.endswith(".="):
            if pre_map:
                if sequence_id is None:
                    msg = "Must provide sequence id to construct pre-mapped ref-identical VRS allele"
                    raise ValueError(msg)
                allele = build_ref_identical_allele(sequence_id)
            else:
                allele = translate_ref_identical_to_vrs(hgvs_string)

            alleles.append(normalize_and_identify(allele))
            continue

        allele = translate_hgvs_to_vrs(hgvs_string)

        if pre_map:
            if sequence_id is None:
                msg = "Must provide sequence id to construct pre-mapped VRS allele"
                raise ValueError(msg)
            allele.location.sequenceReference.refgetAccession = sequence_id
        else:
            allele.location.sequenceReference.label = hgvs_string.split(":")[0]

        if "dup" in hgvs_string:
            allele.state.sequence = SequenceString(2 * _get_allele_sequence(allele))

        if allele.state.sequence.root == "N" and layer == AnnotationLayer.GENOMIC:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele))

        if "=" in hgvs_string and layer == AnnotationLayer.PROTEIN:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele))

        allele = normalize(allele, data_proxy=get_seqrepo())

        if isinstance(allele.state, ReferenceLengthExpression):
            _logger.debug(
                "Coercing state for %s into LSE: %s",
                hgvs_string,
                allele.state.model_dump_json(),
            )
            allele.state = _rle_to_lse(allele.state, allele.location)

        # Carry the verbatim c. HGVS so annotate need not reconstruct it (it can't: its
        # reconstructor only handles g./p. frames).
        if not pre_map and layer is AnnotationLayer.CDNA:
            allele.expressions = [Expression(syntax=Syntax.HGVS_C, value=hgvs_string)]

        allele.id = identify_allele(allele)
        alleles.append(allele)

    if not alleles:
        msg = f"Input variant hgvs_string(s) could not be translated to an allele: {hgvs_strings}."
        raise ValueError(msg)

    if len(alleles) > 1:
        return Haplotype(members=alleles)

    return alleles[0]


def vrs_map(
    metadata: TargetGene,
    align_result: AlignmentResult | None,
    records: list[ScoreRow],
    transcript: TxSelectResult | TxSelectError | None = None,
    silent: bool = True,
) -> VrsMapResult:
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    :param metadata: target gene metadata from MaveDB API
    :param align_result: output from the sequence alignment process
    :param records: scoreset records
    :param transcript: output of transcript selection process
    :param silent: If true, suppress console output
    :return: :class:`~dcd_mapping.schemas.VrsMapResult` with ``mappings`` and
        ``protein_align_result``.  ``protein_align_result`` is the
        target-protein-to-reference-protein alignment used for protein-level
        mapping, surfaced so downstream annotation can flag protein-layer
        variants at mismatched/near-gap loci. It is ``None`` for accession-based
        and regulatory/non-coding targets.
    """
    if metadata.target_accession_id:
        return VrsMapResult(
            mappings=_map_accession(metadata, records, align_result, transcript),
            protein_align_result=None,
        )

    if metadata.target_gene_category == TargetType.PROTEIN_CODING and transcript:
        mappings, protein_align_result = _map_protein_coding(
            metadata,
            records,
            transcript,
            align_result,
            silent,
        )
        return VrsMapResult(
            mappings=mappings, protein_align_result=protein_align_result
        )

    return VrsMapResult(
        mappings=_map_regulatory_noncoding(metadata, records, align_result),
        protein_align_result=None,
    )
