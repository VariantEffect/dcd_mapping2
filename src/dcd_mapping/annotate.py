"""Annotate MaveDB score set metadata with mapped scores."""

import datetime
import json
import logging
from pathlib import Path

import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant
from Bio.SeqUtils import seq3
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u
from ga4gh.core._internal.models import Extension
from ga4gh.vrs._internal.models import (
    Allele,
    Expression,
    Haplotype,
    LiteralSequenceExpression,
    ReferenceLengthExpression,
    Syntax,
)

from dcd_mapping import vrs_v1_schemas
from dcd_mapping.lookup import (
    _get_hgnc_symbol,
    get_chromosome_identifier,
    get_chromosome_identifier_from_vrs_id,
    get_gene_symbol,
    get_overlapping_features_for_region,
    get_seqrepo,
    get_ucsc_chromosome_name,
    get_vrs_id_from_identifier,
)
from dcd_mapping.resource_utils import LOCAL_STORE_PATH
from dcd_mapping.schemas import (
    AlignmentResult,
    ComputedReferenceSequence,
    GeneInfo,
    MappedReferenceSequence,
    MappedScore,
    ScoreAnnotation,
    ScoreAnnotationWithLayer,
    ScoresetMapping,
    ScoresetMetadata,
    TargetAnnotation,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
    VrsVersion,
)
from dcd_mapping.transcripts import TxSelectError

_logger = logging.getLogger(__name__)


async def compute_target_gene_info(
    target_key: str,
    transcripts: dict[str, TxSelectResult | TxSelectError | None],
    alignment_results: dict[str, AlignmentResult | None],
    metadata: ScoresetMetadata,
    mapped_scores: list[MappedScore] | None = None,
) -> GeneInfo | None:
    """Determine a single gene symbol per target with provenance.

    Priority:
    1. HGNC symbol from selected transcript when present
    2. Overlap-based inference from alignment hit subranges
    3. Overlap-based inference from variant spans
    4. Fallback to normalized target metadata symbol
    """
    # If a target is not coding, we can't select it's gene.
    # TODO#66: Handle regulatory/non-coding targets more intelligently. Our current method of querying
    #          Ensembl for overlap isn't robust for these cases.
    if (
        metadata.target_genes[target_key].target_gene_category
        != TargetType.PROTEIN_CODING
    ):
        _logger.info(
            "Target %s is not protein coding. Skipped computing target gene.",
            target_key,
        )
        return GeneInfo(hgnc_symbol=None, selection_method="target_category")

    # Prefer returning the HGNC symbol directly from a selected transcript.
    try:
        tx = transcripts.get(target_key)
        if tx and isinstance(tx, TxSelectResult):
            _logger.info(
                "Using selected transcript for gene info for target %s. Computed HGNC: %s",
                target_key,
                tx.hgnc_symbol,
            )
            return GeneInfo(hgnc_symbol=tx.hgnc_symbol, selection_method="tx_selection")

    except Exception:
        _logger.exception("Error computing target gene info for target %s", target_key)
        return None

    # If we cannot compute gene info via a selected transcript, try to infer a gene symbol from alignment results.
    try:
        gene_info = _compute_target_gene_info_from_alignment(
            target_key, alignment_results.get(target_key)
        )

        if gene_info is not None:
            return gene_info

    except Exception:
        _logger.exception("Error computing target gene info for target %s", target_key)
        return None

    # If we cannot infer a transcript from transcript selection or alignment results, try to compute gene info
    # from mapped variant spans.
    try:
        gene_info = _compute_target_gene_info_from_mapped_variant_spans(
            target_key, mapped_scores
        )

        if gene_info is not None:
            return gene_info

    except Exception:
        _logger.exception(
            "Error computing gene info from mapped variant spans for target %s",
            target_key,
        )
        return None

    # Fallback to target metadata normalization in cases where inference fails.
    try:
        symbol = get_gene_symbol(metadata.target_genes[target_key])
        if symbol:
            _logger.warning(
                "Using target metadata for gene info for target %s. Computed HGNC: %s",
                target_key,
                symbol,
            )
            return GeneInfo(hgnc_symbol=symbol, selection_method="target_metadata")

    except Exception:
        _logger.exception(
            "Error computing gene info from target metadata for target %s",
            target_key,
        )
        return None

    _logger.warning(
        "No gene info could be determined for target %s",
        target_key,
    )
    return None


def _compute_target_gene_info_from_mapped_variant_spans(
    target_key: str, mapped_scores: list[MappedScore] | None
) -> GeneInfo | None:
    if mapped_scores is None:
        return None

    spans = _iter_genomic_spans_from_mapped_scores(mapped_scores)

    # Although multiple chromosomes being present in the same target seems exceedingly unlikely, we handle it just in case
    # given it isn't much additional complexity or extra work.
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in spans:
        by_chrom.setdefault(chrom, []).append((start, end))

    # interval merging helper. Sorts intervals, then merges any which overlap.
    # Merging intervals will reduce the number of overlap queries we need to perform.
    def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
        if not intervals:
            return []

        intervals.sort()
        m = [intervals[0]]
        for s, e in intervals[1:]:
            ls, le = m[-1]
            if s <= le:
                m[-1] = (ls, max(le, e))
            else:
                m.append((s, e))

        return m

    covered_bases: dict[str, int] = {}
    for chrom, intervals in by_chrom.items():
        chromosomal_coverage = (
            _covered_bases_from_overlapping_genes_of_chromosomal_intervals(
                chrom, merge(intervals)
            )
        )

        # Since HGNC terms are unambiguous and will not exist on multiple chromosomes,
        # we can safely merge these dictionaries.
        covered_bases = {**covered_bases, **chromosomal_coverage}

    if not covered_bases:
        return None

    # Find the gene with the maximum coverage.
    max_cov = max(covered_bases.values())
    candidates = sorted([g for g, cov in covered_bases.items() if cov == max_cov])

    # It doesn't seem guaranteed that the first label in the candidates list will be stable, nor
    # that it will be the correct choice if there are multiple candidates.
    if len(candidates) > 1:
        _logger.warning(
            "Multiple genes with maximum coverage for target %s: %s. No gene info will be computed for this target.",
            target_key,
            candidates,
        )
        return None

    _logger.info(
        "Using mapped variant spans for gene info for target %s. Computed HGNC: %s",
        target_key,
        candidates[0],
    )
    return GeneInfo(
        hgnc_symbol=_get_hgnc_symbol(candidates[0]),
        selection_method="variants_max_covered_bases",
    )


def _compute_target_gene_info_from_alignment(
    target_key: str, alignment_result: AlignmentResult | None
) -> GeneInfo | None:
    if alignment_result is None:
        return None

    chrom = get_chromosome_identifier(alignment_result.chrom)

    # Find all 'gene' features from Ensembl which at least partially overlap each hit range.
    aligned_intervals = [(sub.start, sub.end) for sub in alignment_result.hit_subranges]
    covered_bases = _covered_bases_from_overlapping_genes_of_chromosomal_intervals(
        chrom, aligned_intervals
    )

    if not covered_bases:
        return None

    # Find the gene with the maximum coverage.
    max_cov = max(covered_bases.values())
    candidates = sorted([g for g, cov in covered_bases.items() if cov == max_cov])

    # It doesn't seem guaranteed that the first label in the candidates list will be stable, nor
    # that it will be the correct choice if there are multiple candidates.
    if len(candidates) > 1:
        _logger.warning(
            "Multiple genes with maximum coverage for target %s: %s. No gene info will be computed for this target.",
            target_key,
            candidates,
        )
        return None

    _logger.info(
        "Using alignment results for gene info for target %s. Computed HGNC: %s",
        target_key,
        candidates[0],
    )
    return GeneInfo(
        hgnc_symbol=_get_hgnc_symbol(candidates[0]),
        selection_method="alignment_max_covered_bases",
    )


def _covered_bases_from_overlapping_genes_of_chromosomal_intervals(
    chromosome: str,
    intervals: list[tuple[int, int]],
) -> dict[str, int]:
    """Compute the number of bases from gene features that overlap given chromosomal intervals.

    This function iterates over a list of genomic intervals, queries overlapping
    gene features, and sums the number of overlapping bases per
    HGNC symbol across all intervals. Intervals that cause errors during feature
    lookup are skipped, and features missing start or end coordinates are ignored.

    Parameters
    ----------
    chromosome : str
        Chromosome identifier (e.g., "1", "chr1", "X").
    intervals : list[tuple[int, int]]
        A list of (start, end) tuples describing 0-based chromosomal intervals,
        where start is inclusive and end is exclusive.

    Returns
    -------
    dict[str, int]
        A mapping from HGNC gene symbol to the total number of bases that overlap
        the provided intervals.

    Notes
    -----
    - Overlap for each feature is computed as: max(0, min(interval_end, feature_end) - max(interval_start, feature_start)).
    - If feature lookup raises an exception, that interval contributes nothing.
    - Features without HGNC symbol or valid start/end are skipped.

    """
    covered_bases: dict[str, int] = {}
    for s, e in intervals:
        try:
            matches = get_overlapping_features_for_region(
                chromosome, s, e, features=["gene"]
            )
        except Exception:
            _logger.exception(
                "Error fetching overlapping gene features for region %s:%d-%d",
                chromosome,
                s,
                e,
            )
            matches = []

        # The overlapping bases of each feature within our aligned region, 8 in the example below:
        # feature:  ----------------
        # interval:         --------------
        #                   ^^^^^^^^
        for feature in matches:
            hgnc = feature.get("external_name")
            feature_start = feature.get("start")
            feature_end = feature.get("end")

            if hgnc and feature_start is not None and feature_end is not None:
                overlapping_bases = max(0, min(e, feature_end) - max(s, feature_start))
                covered_bases[hgnc] = covered_bases.get(hgnc, 0) + overlapping_bases
                _logger.debug(
                    "Feature %s overlaps region %s:%d-%d by %d bases",
                    hgnc,
                    chromosome,
                    s,
                    e,
                    overlapping_bases,
                )
            else:
                _logger.warning(
                    "Skipping feature with missing HGNC symbol or invalid coordinates: %s",
                    feature,
                )

    return covered_bases


def _iter_genomic_spans_from_mapped_scores(
    mapped_scores: list[MappedScore],
) -> list[tuple[str, int, int]]:
    """Extract (chrom, start, end) spans from post-mapped VRS structures.

    Only considers genomic annotations. Returns an empty list if none.
    """
    spans: list[tuple[str, int, int]] = []
    for ms in mapped_scores:
        if ms.annotation_layer != AnnotationLayer.GENOMIC or ms.post_mapped is None:
            continue
        if isinstance(ms.post_mapped, Allele):
            loc = ms.post_mapped.location
            refget_chrom = get_chromosome_identifier_from_vrs_id(
                f"ga4gh:{loc.sequenceReference.refgetAccession}"
            )
            if refget_chrom:
                spans.append(
                    (get_ucsc_chromosome_name(refget_chrom), loc.start, loc.end)
                )

        elif isinstance(ms.post_mapped, Haplotype):
            for allele in ms.post_mapped.members:
                loc = allele.location
                refget_chrom = get_chromosome_identifier_from_vrs_id(
                    f"ga4gh:{loc.sequenceReference.refgetAccession}"
                )
                if refget_chrom:
                    spans.append(
                        (get_ucsc_chromosome_name(refget_chrom), loc.start, loc.end)
                    )

    return spans


def _allele_to_v1_allele(allele: Allele) -> vrs_v1_schemas.Allele:
    """Convert VRS 2.0 allele to VRS 1.3 allele.

    :param allele: VRS 2.0a allele
    :return: equivalent VRS 1.3 allele
    """
    start = allele.location.start
    end = allele.location.end
    sequence_id = f"ga4gh:{allele.location.sequenceReference.refgetAccession}"  # type: ignore
    location_raw = f'{{"end":{{"type":"Number","value":{end}}},"sequence_id":"{sequence_id.split(".")[1]}","start":{{"type":"Number","value":{start}}},"type":"SequenceLocation"}}'
    location_id = sha512t24u(location_raw.encode("ascii"))
    sequence = "" if not allele.state.sequence else allele.state.sequence.root
    allele_raw = f'{{"location":"{location_id}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'
    allele_id = sha512t24u(allele_raw.encode("ascii"))

    return vrs_v1_schemas.Allele(
        id=f"ga4gh:VA.{allele_id}",
        location=vrs_v1_schemas.SequenceLocation(
            id=location_id,
            sequence_id=sequence_id,
            interval=vrs_v1_schemas.SequenceInterval(
                start=vrs_v1_schemas.Number(value=start, type="Number"),
                end=vrs_v1_schemas.Number(value=end, type="Number"),
            ),
        ),
        state=vrs_v1_schemas.LiteralSequenceExpression(sequence=sequence),
    )


def _allele_to_vod(allele: Allele) -> vrs_v1_schemas.VariationDescriptor:
    """Convert VRS 2.0 allele to comparable VRSATILE VariationDescriptor.

    Some allele properties aren't available in the 1.3 allele, so we have to lift them
    up to the VariationDescriptor.
    """
    allele_v1 = _allele_to_v1_allele(allele)
    if allele.expressions:
        original_expression = allele.expressions[0]
        expressions = [
            vrs_v1_schemas.Expression(
                syntax=original_expression.syntax,
                value=original_expression.value,
                syntax_version=None,
            )
        ]
    else:
        expressions = []
    return vrs_v1_schemas.VariationDescriptor(
        id=allele_v1.id,
        variation=allele_v1,
        type="VariationDescriptor",
        expressions=expressions,
        vrs_ref_allele_seq=allele.extensions[0].value,
        extensions=[],
    )


def _haplotype_to_haplotype_1_3(haplotype: Haplotype) -> vrs_v1_schemas.Haplotype:
    """Convert VRS 2.0 Haplotype to VRS 1.3 Haplotype.

    :param haplotype: VRS 2.0 haplotype
    :return: VRS 1.3 haplotype that contains VRSATILE variation descriptors (not alleles)
    """
    members = []
    allele: Allele
    for allele in haplotype.members:  # type: ignore
        members.append(_allele_to_vod(allele))
    return vrs_v1_schemas.Haplotype(members=members)


def _offset_allele_ref_seq(ss: str, start: int, end: int) -> tuple[int, int]:
    """Handle known edge cases in start and end coordinates for vrs_ref_allele_seq."""
    if ss.startswith("urn:mavedb:00000060-a-1"):
        _logger.warning(
            "urn:mavedb:00000060-a-1 reports the entire human reference sequence as the target sequence, but the start and end positions need to be manually offset by 289"
        )
        return (start + 289, end + 289)
    if ss.startswith("urn:mavedb:00000060-a-2"):
        _logger.warning(
            "urn:mavedb:00000060-a-2 reports the entire human reference sequence as the target sequence, but the start and end positions need to be manually offset by 331"
        )
        return (start + 331, end + 331)
    return (start, end)


def _get_vrs_ref_allele_seq(
    allele: Allele,
    metadata: TargetGene,
    urn: str,
    tx_select_results: TxSelectResult | None,
) -> Extension | None:
    """Create `vrs_ref_allele_seq` property."""
    start, end = _offset_allele_ref_seq(urn, allele.location.start, allele.location.end)
    if (
        urn.startswith(
            (
                "urn:mavedb:00000047",
                "urn:mavedb:00000048",
                "urn:mavedb:00000053",
                "urn:mavedb:00000058-a-1",
            )
        )
        and tx_select_results is not None
        and isinstance(tx_select_results, TxSelectResult)
    ):
        seq = tx_select_results.sequence
        ref = seq[start:end]
    else:
        seq = f"ga4gh:{allele.location.sequenceReference.refgetAccession}"  # type: ignore
        sr = get_seqrepo()
        ref = sr.get_sequence(seq, start, end)

    if not ref:
        msg = f"Could not retrieve reference sequence for allele {allele.id} in urn {urn} with start {start} and end {end}"
        _logger.warning(msg)
        return None

    return Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)


def _get_hgvs_string(allele: Allele, accession: str) -> tuple[str, Syntax]:
    """Return an HGVS string for a given VRS allele

    :param allele: A post-mapped VRS allele
    :param accession: A RefSeq accession
    :return An HGVS string and the Syntax value
    """
    if accession.startswith("NP"):
        syntax = Syntax.HGVS_P
        syntax_value = "p"
    else:
        syntax = Syntax.HGVS_G
        syntax_value = "g"

    # Reference-identical variants are represented as simple Alleles with a ReferenceLengthExpression state rather than being translated from HGVS.
    # We can generate a simplified HGVS string for these variants without needing to perform a full translation.
    if isinstance(allele.state, ReferenceLengthExpression):
        return f"{accession}:{syntax_value}.=", syntax

    start: int = allele.location.start
    end: int = allele.location.end

    dp = get_seqrepo()
    if start == end:
        ref = None
        aas = dp.get_sequence(accession, start - 1, start)
        aae = dp.get_sequence(accession, end, end + 1)
        end += 1
    else:
        ref = dp.get_sequence(accession, start, end)
        aas = dp.get_sequence(accession, start, start + 1)
        aae = dp.get_sequence(accession, end - 1, end)
        start += 1

    if syntax == Syntax.HGVS_P:
        ival = hgvs.location.Interval(
            start=hgvs.location.AAPosition(base=start, aa=aas),
            end=hgvs.location.AAPosition(base=end, aa=aae),
        )
    else:
        ival = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=start),
            end=hgvs.location.SimplePosition(base=end),
        )
    if isinstance(allele.state, LiteralSequenceExpression):
        alt = allele.state.sequence.root
    else:
        msg = (
            f"Unable to handle string for non-LSE based allele in {allele.model_dump()}"
        )
        raise NotImplementedError(msg)

    edit = ""  # empty by default
    if alt == ref:
        edit = "="
    if ref and (2 * ref == alt or len(ref) == 1 and set(ref) == set(alt)):
        edit = "dup"
    if alt == "":
        edit = "del"

    if edit != "dup" or edit != "del" or edit != "=":
        edit = (
            hgvs.edit.AARefAlt(ref=ref, alt=alt)
            if syntax == Syntax.HGVS_P
            else hgvs.edit.NARefAlt(ref=ref, alt=alt)
        )

    if alt != ref:
        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
    else:
        posedit = (
            f"{seq3(ref)}{start!s}=" if syntax == Syntax.HGVS_P else f"{end!s}{ref}="
        )

    var = str(
        hgvs.sequencevariant.SequenceVariant(
            ac=accession, type=syntax_value, posedit=posedit
        )
    )
    if var.endswith("delins"):
        var = var.replace("delins", "del")
    return var, syntax


def _annotate_allele_mapping(
    mapped_score: MappedScore,
    tx_results: TxSelectResult | TxSelectError | None,
    metadata: TargetGene,
    urn: str,
    vrs_version: VrsVersion = VrsVersion.V_2,
) -> ScoreAnnotationWithLayer:
    """Perform annotations and, if necessary, create VRS 1.3 equivalents for allele mappings."""
    pre_mapped: Allele = mapped_score.pre_mapped
    post_mapped: Allele = mapped_score.post_mapped

    # get vrs_ref_allele_seq for pre-mapped variants if they aren't reference-identical variants, which have a ReferenceLengthExpression state
    # and for which the vrs_ref_allele_seq would be redundant with the length and sequence reference information already present in the allele.
    # We also want to avoid fetching the reference sequence for long reference-identical variants, as this can cause performance issues and the
    # vrs_ref_allele_seq doesn't add much value in these cases.
    if not isinstance(pre_mapped.state, ReferenceLengthExpression):
        ref_allele_seq_extension = _get_vrs_ref_allele_seq(
            pre_mapped, metadata, urn, tx_results
        )
        if ref_allele_seq_extension is not None:
            pre_mapped.extensions = [ref_allele_seq_extension]

    if post_mapped:
        # Determine reference sequence
        if mapped_score.annotation_layer == AnnotationLayer.GENOMIC:
            sequence_id = f"ga4gh:{mapped_score.post_mapped.location.sequenceReference.refgetAccession}"
            accession = get_chromosome_identifier_from_vrs_id(sequence_id)
            if accession is None:
                accession = None
                mapped_score.error_message = "Could not determine accession for this annotation. No allele expression is available."
            elif accession.startswith("refseq:"):
                accession = accession[7:]
        else:
            if tx_results is None or isinstance(tx_results, TxSelectError):
                accession = None
                mapped_score.error_message = "Could not determine accession for this annotation. No allele expression is available."
            else:
                accession = tx_results.np

        sr = get_seqrepo()
        loc = mapped_score.post_mapped.location
        sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"

        # Skip getting refereence sequence for RLE Alleles, see above for pre-mapped alleles.
        if not isinstance(post_mapped.state, ReferenceLengthExpression):
            ref = sr.get_sequence(sequence_id, loc.start, loc.end)
            post_mapped.extensions = [
                Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)
            ]

        if accession:
            hgvs_string, syntax = _get_hgvs_string(post_mapped, accession)
            post_mapped.expressions = [Expression(syntax=syntax, value=hgvs_string)]

    if vrs_version == VrsVersion.V_1_3:
        pre_mapped = _allele_to_vod(pre_mapped)
        post_mapped = _allele_to_vod(post_mapped) if post_mapped else None

    return ScoreAnnotationWithLayer(
        pre_mapped=pre_mapped,
        post_mapped=post_mapped,
        vrs_version=vrs_version,
        mavedb_id=mapped_score.accession_id,
        score=float(mapped_score.score) if mapped_score.score is not None else None,
        annotation_layer=mapped_score.annotation_layer,
        error_message=mapped_score.error_message,
    )


def _annotate_haplotype_mapping(
    mapped_score: MappedScore,
    tx_results: TxSelectResult | TxSelectError | None,
    metadata: TargetGene,
    urn: str,
    vrs_version: VrsVersion = VrsVersion.V_2,
) -> ScoreAnnotationWithLayer:
    """Perform annotations and, if necessary, create VRS 1.3 equivalents for haplotype mappings."""
    pre_mapped: Haplotype = mapped_score.pre_mapped  # type: ignore
    post_mapped: Haplotype = mapped_score.post_mapped  # type: ignore

    # see comment in _annotate_allele_mapping regarding why we skip getting vrs_ref_allele_seq for reference-identical variants.
    for allele in pre_mapped.members:
        if not isinstance(allele.state, ReferenceLengthExpression):
            ref_allele_seq_extension = _get_vrs_ref_allele_seq(
                allele, metadata, urn, tx_results
            )
            if ref_allele_seq_extension is not None:
                allele.extensions = [ref_allele_seq_extension]

    if post_mapped:
        # Determine reference sequence
        if mapped_score.annotation_layer == AnnotationLayer.GENOMIC:
            sequence_id = f"ga4gh:{post_mapped.members[0].location.sequenceReference.refgetAccession}"
            accession = get_chromosome_identifier_from_vrs_id(sequence_id)
            if accession is None:
                accession = None
                mapped_score.error_message = "Could not determine accession for this annotation. No allele expression is available."
            elif accession.startswith("refseq:"):
                accession = accession[7:]
        else:
            if tx_results is None or isinstance(tx_results, TxSelectError):
                # impossible by definition
                accession = None
                mapped_score.error_message = "Could not determine accession for this annotation. No allele expression is available."
            else:
                accession = tx_results.np

        sr = get_seqrepo()
        for allele in post_mapped.members:
            loc = allele.location
            sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"

            # Again, skip getting reference sequence for RLE Alleles.
            if not isinstance(allele.state, ReferenceLengthExpression):
                ref = sr.get_sequence(
                    sequence_id, loc.start, loc.end
                )  # TODO type issues??
                allele.extensions = [
                    Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)
                ]

            if accession:
                hgvs, syntax = _get_hgvs_string(allele, accession)
                allele.expressions = [Expression(syntax=syntax, value=hgvs)]

    if vrs_version == VrsVersion.V_1_3:
        pre_mapped = _haplotype_to_haplotype_1_3(pre_mapped)
        post_mapped = _haplotype_to_haplotype_1_3(post_mapped) if post_mapped else None

    return ScoreAnnotationWithLayer(
        pre_mapped=pre_mapped,
        post_mapped=post_mapped,
        vrs_version=vrs_version,
        mavedb_id=mapped_score.accession_id,
        score=float(mapped_score.score) if mapped_score.score is not None else None,
        annotation_layer=mapped_score.annotation_layer,
        error_message=mapped_score.error_message,
    )


def annotate(
    mapped_scores: list[MappedScore],
    tx_results: TxSelectResult | TxSelectError | None,
    metadata: TargetGene,
    urn: str,
    vrs_version: VrsVersion = VrsVersion.V_2,
) -> list[ScoreAnnotationWithLayer]:
    """Given a list of mappings, add additional contextual data:

    1. ``vrs_ref_allele_seq``: The sequence between the start and end positions
        indicated in the variant
    2. ``hgvs``: An HGVS string describing the variant (only included for post-mapped
        variants)
    3. ``transcript_accession``: A description of the MANE annotation of the transcript,
        if any  # < -- TODO I think we took this out

    ...and provide VRS 1.3-converted equivalents, too.

    :param vrs_results: in-progress variant mappings
    :param tx_select_results: transcript selection if available
    :param metadata: Target gene metadata from MaveDB API
    :param urn: Score set URN
    :return: annotated mappings objects
    """
    score_annotations = []
    for mapped_score in mapped_scores:
        if mapped_score.pre_mapped is None:
            score_annotations.append(
                ScoreAnnotationWithLayer(
                    mavedb_id=mapped_score.accession_id,
                    score=float(mapped_score.score)
                    if mapped_score.score is not None
                    else None,
                    vrs_version=vrs_version,
                    error_message=mapped_score.error_message,
                )
            )
        elif isinstance(mapped_score.pre_mapped, Haplotype) and (
            isinstance(mapped_score.post_mapped, Haplotype)
            or mapped_score.post_mapped is None
        ):
            score_annotations.append(
                _annotate_haplotype_mapping(
                    mapped_score, tx_results, metadata, urn, vrs_version
                )
            )
        elif isinstance(mapped_score.pre_mapped, Allele) and (
            isinstance(mapped_score.post_mapped, Allele)
            or mapped_score.post_mapped is None
        ):
            score_annotations.append(
                _annotate_allele_mapping(
                    mapped_score, tx_results, metadata, urn, vrs_version
                )
            )
        else:
            score_annotations.append(
                ScoreAnnotationWithLayer(
                    pre_mapped=mapped_score.pre_mapped,
                    post_mapped=mapped_score.post_mapped,
                    vrs_version=vrs_version,
                    mavedb_id=mapped_score.accession_id,
                    score=float(mapped_score.score)
                    if mapped_score.score is not None
                    else None,
                    error_message=f"Multiple issues with annotation: Inconsistent variant structure (Allele and Haplotype mix).{' ' + mapped_score.error_message if mapped_score.error_message else ''}",
                )
            )

    return score_annotations


def _get_computed_reference_sequence(
    metadata: TargetGene,
    layer: AnnotationLayer,
    tx_output: TxSelectResult | TxSelectError | None = None,
) -> ComputedReferenceSequence | MappedReferenceSequence | None:
    """Report the computed reference sequence for a score set

    :param metadata: Target gene metadata from MaveDB API
    :param layer: AnnotationLayer
    :param tx_output: Transcript data for a score set
    :return A ComputedReferenceSequence object,
    or if the target gene is accession-based, a mapped reference sequence describing the pre-mapped reference
    """
    # accession-based target genes always use accession id as pre-mapped reference sequence
    if metadata.target_accession_id:
        seq_id = get_vrs_id_from_identifier(metadata.target_accession_id)
        # use MappedReferenceSequence type because there should be an accession id but no sequence.
        # for accession-based target genes, the object returned by this function describes the provided reference accession
        # whereas the object returned by _get_mapped_reference_sequence describes the mapped reference accession, which could be a chromosome for ex.
        seq_type: TargetSequenceType
        if metadata.target_accession_id.startswith(("NP", "ENSP")):
            seq_type = TargetSequenceType.PROTEIN
        elif metadata.target_accession_id.startswith(("NM", "ENST", "NC")):
            seq_type = TargetSequenceType.DNA
        else:
            msg = f"Unrecognized accession prefix for accession id {metadata.target_accession_id}"
            raise ValueError(msg)
        return MappedReferenceSequence(
            sequence_type=seq_type,
            sequence_id=seq_id,
            sequence_accessions=[metadata.target_accession_id],
        )
    if layer == AnnotationLayer.PROTEIN:
        if tx_output is None or isinstance(tx_output, TxSelectError):
            return None
        seq_id = f"ga4gh:SQ.{sha512t24u(tx_output.sequence.encode('ascii'))}"
        return ComputedReferenceSequence(
            sequence=tx_output.sequence,
            sequence_type=TargetSequenceType.PROTEIN,
            sequence_id=seq_id,
        )
    seq_id = f"ga4gh:SQ.{sha512t24u(metadata.target_sequence.encode('ascii'))}"
    return ComputedReferenceSequence(
        sequence=metadata.target_sequence,
        sequence_type=TargetSequenceType.DNA,
        sequence_id=seq_id,
    )


def _get_mapped_reference_sequence(
    metadata: TargetGene,
    layer: AnnotationLayer,
    tx_output: TxSelectResult | TxSelectError | None = None,
    align_result: AlignmentResult | None = None,
) -> MappedReferenceSequence | None:
    """Report the mapped reference sequence for a score set

    :param metadata: Target gene metadata from MaveDB API
    :param layer: AnnotationLayer
    :param tx_output: Transcript data for a score set
    :return A MappedReferenceSequence object
    """
    if (
        layer == AnnotationLayer.PROTEIN
        and tx_output is not None
        and isinstance(tx_output, TxSelectResult)
    ):
        if tx_output.np is None:
            return None
        vrs_id = get_vrs_id_from_identifier(tx_output.np)
        if vrs_id is None:
            return None
        return MappedReferenceSequence(
            sequence_type=TargetSequenceType.PROTEIN,
            sequence_id=vrs_id,
            sequence_accessions=[tx_output.np],
        )
    # accession-based score sets with genomic accession do not have alignment results
    if (
        align_result is None
        and metadata.target_accession_id
        and metadata.target_accession_id.startswith("NC")
    ):
        seq_id = metadata.target_accession_id
    else:
        seq_id = get_chromosome_identifier(align_result.chrom)
    vrs_id = get_vrs_id_from_identifier(seq_id)
    if vrs_id is None:
        return None
    return MappedReferenceSequence(
        sequence_type=TargetSequenceType.DNA,
        sequence_id=vrs_id,
        sequence_accessions=[seq_id],
    )


def _set_scoreset_layer(
    urn: str, mappings: list[ScoreAnnotationWithLayer]
) -> AnnotationLayer:
    """Many individual score results provide both genomic and protein variant
    expressions. If genomic expressions are available, that's what we'd like to use.
    This function tells us how to filter the final annotation objects.
    """
    for mapping in mappings:
        if mapping.annotation_layer == AnnotationLayer.GENOMIC:
            return AnnotationLayer.GENOMIC
    return AnnotationLayer.PROTEIN


def write_scoreset_mapping_to_json(
    urn: str,
    scoreset_mapping: ScoresetMapping,
    output_path: Path | None = None,
) -> Path:
    """Write the given ScoresetMapping as a JSON at the specified
    or default ouput path.
    """
    if not output_path:
        now = datetime.datetime.now(tz=datetime.UTC).isoformat()
        output_path = LOCAL_STORE_PATH / f"{urn}_mapping_{now}.json"

    _logger.info("Saving mapping output to %s", output_path)
    with output_path.open("w") as file:
        json.dump(
            scoreset_mapping.model_dump(exclude_unset=False, exclude_none=True),
            file,
            indent=4,
        )
    return output_path


def save_mapped_output_json(
    metadata: ScoresetMetadata,
    mappings: dict[str, list[ScoreAnnotationWithLayer]],
    align_results: dict[str, AlignmentResult],
    tx_output: dict[str, TxSelectResult | TxSelectError | None],
    gene_info: dict[str, GeneInfo | None],
    preferred_layer_only: bool = False,
    output_path: Path | None = None,
) -> Path:
    """Save mapping output for a score set in a JSON file

    :param urn: Score set accession
    :param mappings: A dictionary of annotated VRS mappings for each target
    :param align_result: A dictionary of alignment information for each target
    :param tx_output: A dictionary of transcript output for each target
    :param output_path: specific location to save output to. Default to
        <dcd_mapping_data_dir>/urn:mavedb:00000XXX-X-X_mapping_<ISO8601 datetime>.json
    :return: output location
    """
    # set preferred layers for each target, to allow a mix of coding and noncoding targets
    reference_sequences: dict[str, TargetAnnotation] = {}
    mapped_scores: list[ScoreAnnotation] = []
    for target_gene in mappings:
        if preferred_layer_only:
            preferred_layers = {
                _set_scoreset_layer(metadata.urn, mappings[target_gene]),
            }
        else:
            preferred_layers = {
                mapping.annotation_layer for mapping in mappings[target_gene]
            }

        # use target gene name in reference sequence dictionary, rather than the label, which differs between score sets
        target_gene_name = metadata.target_genes[target_gene].target_gene_name

        reference_sequences[target_gene_name] = TargetAnnotation(
            gene_info=gene_info.get(target_gene),
            layers={
                layer: {
                    "computed_reference_sequence": None,
                    "mapped_reference_sequence": None,
                }
                for layer in preferred_layers
            },
        )
        # sometimes Nonetype layers show up in preferred layers dict; remove these
        preferred_layers.discard(None)
        for layer in preferred_layers:
            reference_sequences[target_gene_name].layers[layer][
                "computed_reference_sequence"
            ] = _get_computed_reference_sequence(
                metadata.target_genes[target_gene], layer, tx_output[target_gene]
            )
            reference_sequences[target_gene_name].layers[layer][
                "mapped_reference_sequence"
            ] = _get_mapped_reference_sequence(
                metadata.target_genes[target_gene],
                layer,
                tx_output[target_gene],
                align_results[target_gene],
            )

        # if genomic layer, not accession-based, and target gene type is coding, add cdna entry (just the sequence accession) to reference_sequences dict
        if (
            AnnotationLayer.GENOMIC in reference_sequences[target_gene_name].layers
            and metadata.target_genes[target_gene].target_gene_category
            == TargetType.PROTEIN_CODING
            and metadata.target_genes[target_gene].target_accession_id is None
            and tx_output[target_gene] is not None
            and isinstance(tx_output[target_gene], TxSelectResult)
            and tx_output[target_gene].nm is not None
        ):
            reference_sequences[target_gene_name].layers[AnnotationLayer.CDNA] = {
                "computed_reference_sequence": None,
                "mapped_reference_sequence": {
                    "sequence_accessions": [tx_output[target_gene].nm]
                },
            }

        for m in mappings[target_gene]:
            if m.pre_mapped is None:
                mapped_scores.append(ScoreAnnotation(**m.model_dump()))
            elif m.annotation_layer in preferred_layers:
                # drop annotation layer from mapping object
                mapped_scores.append(ScoreAnnotation(**m.model_dump()))

        # drop Nonetype reference sequences
        for target_gene in reference_sequences:
            for layer in list(reference_sequences[target_gene].layers.keys()):
                if (
                    reference_sequences[target_gene].layers[layer][
                        "mapped_reference_sequence"
                    ]
                    is None
                    and reference_sequences[target_gene].layers[layer][
                        "computed_reference_sequence"
                    ]
                    is None
                ) or layer is None:
                    del reference_sequences[target_gene].layers[layer]

    output = ScoresetMapping(
        metadata=metadata.model_dump(),
        reference_sequences=reference_sequences,
        mapped_scores=mapped_scores,
    )

    return write_scoreset_mapping_to_json(metadata.urn, output, output_path)
