"""Annotate MaveDB score set metadata with mapped scores."""

import bisect
import datetime
import json
import logging
from pathlib import Path
from typing import Any

import cdot
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
from dcd_mapping.resource_utils import CDOT_URL, LOCAL_STORE_PATH
from dcd_mapping.schemas import (
    AlignmentResult,
    ComputedReferenceSequence,
    GeneInfo,
    MappedReferenceSequence,
    MappedScore,
    MappingOutcome,
    ScoreAnnotation,
    ScoresetMapping,
    ScoresetMetadata,
    TargetAnnotation,
    TargetGene,
    TargetMapping,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
    VrsVersion,
)
from dcd_mapping.transcripts import TxSelectError
from dcd_mapping.version import dcd_mapping_version

_logger = logging.getLogger(__name__)

# How close (in reference bases) a variant must be to an alignment gap before
# we flag it as ``near_gap``. This is somewhat arbitrary, but seems like a reasonable window to
# capture potential edge effects without being so large that it flags a large
# fraction of variants as near-gap. We might consider re-evaluating this window with more
# rigor in the future, or making it configurable to support situations where callers
# want to be more or less conservative in their definition of "near".
NEAR_GAP_WINDOW = 5


def _allele_ref_positions(post_mapped: object) -> list[tuple[int, int]]:
    """Return the reference-frame half-open intervals covered by an Allele or
    Haplotype-of-Alleles.

    Returns an empty list for:
    - Unrecognised shapes (e.g. pre-VRS-2.0 representations) — caller leaves
      locus flags as None.
    - ReferenceLengthExpression alleles (reference-identical variants) — their
      span covers the entire reference allele unchanged, so per-variant locus
      flags carry no meaningful audit value (and could span entire chromosomes).
    """
    members = (
        post_mapped.members if isinstance(post_mapped, Haplotype) else [post_mapped]
    )
    intervals: list[tuple[int, int]] = []
    for member in members:
        if (
            not isinstance(member, Allele)
            or member.location is None
            or isinstance(member.state, ReferenceLengthExpression)
        ):
            continue

        try:
            start = int(member.location.start)
            end = int(member.location.end)
        except (TypeError, AttributeError, ValueError):
            continue

        # SequenceLocation is half-open [start, end). Ensure we don't produce
        # a zero-length interval when start == end.
        intervals.append((start, max(start + 1, end)))

    return intervals


def _stamp_alignment_locus_flags(
    annotations: list[ScoreAnnotation],
    align_result: AlignmentResult | None,
    target_layer: AnnotationLayer,
    near_gap_window: int = NEAR_GAP_WINDOW,
) -> None:
    """Stamp ``at_mismatched_locus`` and ``near_gap`` on each annotation whose
    ``annotation_layer`` matches ``target_layer``, using positions from the
    given alignment's QC block.

    The caller is responsible for pairing each layer with its own alignment:
    GENOMIC variants get checked against the BLAT genomic alignment; PROTEIN
    variants get checked against the target-protein-to-reference-protein
    alignment. CDNA flags are not currently evaluated -- we'd need either a
    dedicated cDNA alignment or to translate genomic mismatch/gap positions
    through exon structure -- cdna annotations leave the flags as None to
    signal "not evaluated".

    When ``AlignmentQc.mismatch_positions_unavailable`` is ``True`` and
    ``mismatch_count > 0``, the per-base positions are absent but the aggregate
    count is non-zero.  In that case ``at_mismatched_locus`` is left as ``None``
    (not evaluated) for all variants to prevent a silent disagreement between
    ``mismatch_count`` and the per-variant flag.  ``near_gap`` is still evaluated
    from gap intervals, which are always coordinate-derived.
    """
    if align_result is None or align_result.alignment_qc is None:
        return

    qc = align_result.alignment_qc
    mismatches = sorted(qc.mismatch_positions)
    gap_intervals = sorted(qc.gap_intervals)

    # When sequence content was unavailable during QC computation, mismatch_count
    # may be non-zero but mismatch_positions is empty. Marking at_mismatched_locus
    # as False would silently disagree with mismatch_count; leave it as None instead.
    positions_reliable = not (
        qc.mismatch_positions_unavailable and qc.mismatch_count > 0
    )

    layer_anns = [
        ann
        for ann in annotations
        if ann.alignment_level == target_layer and ann.post_mapped is not None
    ]
    if not layer_anns:
        return

    # Explicitly record that this target/layer had no mismatches or gaps. Use
    # mismatch_count (not the positions list) as the source of truth: when
    # mismatch_positions_unavailable is True the list is empty even though
    # mismatches exist, so testing `not mismatches` would be a faulty premise.
    # mismatch_count comes from aln.counts() which is purely coordinate-based
    # and is always accurate regardless of sequence availability.
    if qc.mismatch_count == 0 and not gap_intervals:
        for ann in layer_anns:
            ann.at_mismatched_locus = False
            ann.near_gap = False
        return

    # Pre-compute sorted lists for O(log G) gap overlap queries.
    gap_starts = [gs for gs, _ge in gap_intervals]
    gap_ends = [ge for _gs, ge in gap_intervals]

    for ann in layer_anns:
        intervals = _allele_ref_positions(ann.post_mapped)
        if not intervals:
            # Empty means either RLE allele or unrecognised shape; leave flags None.
            continue

        at_mismatch = False
        near_g = False
        for int_start, int_end in intervals:
            # -- Mismatch check: any mismatch position in [int_start, int_end)?  O(log M)
            if mismatches and not at_mismatch:
                idx = bisect.bisect_left(mismatches, int_start)
                if idx < len(mismatches) and mismatches[idx] < int_end:
                    at_mismatch = True

            # -- Gap proximity check: any gap overlaps [int_start-W, int_end+W)?  O(log G)
            # A gap (gs, ge) overlaps window [L, R) iff gs < R AND ge > L.
            # With sorted (non-overlapping) gap intervals:
            #   Case 1: gap start in [L, R) — definitively overlaps.
            #   Case 2: gap start < L — overlaps only if its end > L.
            if not near_g:
                win_l = int_start - near_gap_window
                win_r = int_end + near_gap_window
                first_at_or_after_l = bisect.bisect_left(gap_starts, win_l)
                first_at_or_after_r = bisect.bisect_left(gap_starts, win_r)
                # At least one gap starts in [win_l, win_r) OR the last gap starting before win_l extends into [win_l, win_r).
                if first_at_or_after_l < first_at_or_after_r or (
                    first_at_or_after_l > 0
                    and gap_ends[first_at_or_after_l - 1] > win_l
                ):
                    near_g = True

            if at_mismatch and near_g:
                break

        # Stamp flags.
        # Note: when positions_reliable is False (mismatch_positions_unavailable
        # AND count>0), at_mismatch will always be False here because mismatches
        # is an empty list — we cannot trust it. Leaving at_mismatched_locus as
        # None (not evaluated) rather than stamping a misleading False. The
        # near_gap branch is unaffected: gap_intervals are coordinate-derived and
        # always reliable, so near_g is evaluated regardless.
        if positions_reliable:
            ann.at_mismatched_locus = at_mismatch
        ann.near_gap = near_g


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
        if ms.alignment_level != AnnotationLayer.GENOMIC or ms.post_mapped is None:
            continue
        if isinstance(ms.post_mapped, Allele):
            loc = ms.post_mapped.location
            refget_chrom = get_chromosome_identifier_from_vrs_id(
                f"ga4gh:{loc.sequenceReference.refgetAccession}"
            )
            if refget_chrom:
                spans.append(
                    (
                        get_ucsc_chromosome_name(refget_chrom),
                        loc.start,
                        loc.end,
                    )
                )

        elif isinstance(ms.post_mapped, Haplotype):
            for allele in ms.post_mapped.members:
                loc = allele.location
                refget_chrom = get_chromosome_identifier_from_vrs_id(
                    f"ga4gh:{loc.sequenceReference.refgetAccession}"
                )
                if refget_chrom:
                    spans.append(
                        (
                            get_ucsc_chromosome_name(refget_chrom),
                            loc.start,
                            loc.end,
                        )
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
    if ref and (2 * ref == alt or (len(ref) == 1 and set(ref) == set(alt))):
        edit = "dup"
    if alt == "":
        edit = "del"

    if edit not in ("dup", "del", "="):
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


def _resolve_outcome(mapped_score: MappedScore) -> MappingOutcome:
    """Resolve the typed :class:`MappingOutcome` for an emitted annotation.

    Projected records arrive with an explicit outcome (``MAPPED`` / ``INTRONIC`` /
    ``NO_PROTEIN_CONSEQUENCE`` / ``FAILED``) -- keep it. Measured records (and any legacy
    record) carry none, so derive it so *every* emitted annotation is typed uniformly:
    ``MAPPED`` when a post-mapped allele was produced, else ``FAILED``.
    """
    if mapped_score.outcome is not None:
        return mapped_score.outcome

    return (
        MappingOutcome.MAPPED
        if mapped_score.post_mapped is not None
        else MappingOutcome.FAILED
    )


def _resolve_postmapped_accession(
    alignment_level: AnnotationLayer | None,
    representative_allele: Allele,
    metadata: TargetGene,
    tx_results: TxSelectResult | TxSelectError | None,
) -> str | None:
    """Resolve the accession to reconstruct a post-mapped HGVS against: the contig for
    genomic, the transcript's protein for protein or coding for cdna.
    ``representative_allele`` is any member for a haplotype. ``None`` if unresolvable
    (the caller surfaces that as an annotation error).
    """
    if alignment_level == AnnotationLayer.GENOMIC:
        sequence_id = (
            f"ga4gh:{representative_allele.location.sequenceReference.refgetAccession}"
        )
        accession = get_chromosome_identifier_from_vrs_id(sequence_id)
        if accession is not None and accession.startswith("refseq:"):
            return accession[len("refseq:") :]
        return accession
    if alignment_level == AnnotationLayer.CDNA:
        return _resolve_cdna_accession(metadata, tx_results, None)

    if tx_results is None or isinstance(tx_results, TxSelectError):
        return None

    return tx_results.np


def _annotate_allele_mapping(
    mapped_score: MappedScore,
    tx_results: TxSelectResult | TxSelectError | None,
    metadata: TargetGene,
    urn: str,
    vrs_version: VrsVersion = VrsVersion.V_2,
) -> ScoreAnnotation:
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
        sr = get_seqrepo()
        loc = mapped_score.post_mapped.location
        sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"

        # Skip getting refereence sequence for RLE Alleles, see above for pre-mapped alleles.
        if not isinstance(post_mapped.state, ReferenceLengthExpression):
            ref = sr.get_sequence(sequence_id, loc.start, loc.end)
            post_mapped.extensions = [
                Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)
            ]

        # Trust a carried HGVS (coding records carry their c.); else reconstruct (g./p. only).
        if not post_mapped.expressions:
            accession = _resolve_postmapped_accession(
                mapped_score.alignment_level, post_mapped, metadata, tx_results
            )
            if accession is None:
                mapped_score.error_message = "Could not determine accession for this annotation. No allele expression is available."
            else:
                hgvs_string, syntax = _get_hgvs_string(post_mapped, accession)
                post_mapped.expressions = [Expression(syntax=syntax, value=hgvs_string)]

    if vrs_version == VrsVersion.V_1_3:
        pre_mapped = _allele_to_vod(pre_mapped)
        post_mapped = _allele_to_vod(post_mapped) if post_mapped else None

    return ScoreAnnotation(
        pre_mapped=pre_mapped,
        post_mapped=post_mapped,
        vrs_version=vrs_version,
        mavedb_id=mapped_score.accession_id,
        score=float(mapped_score.score) if mapped_score.score is not None else None,
        error_message=mapped_score.error_message,
        alignment_level=mapped_score.alignment_level,
        outcome=_resolve_outcome(mapped_score),
    )


def _annotate_haplotype_mapping(
    mapped_score: MappedScore,
    tx_results: TxSelectResult | TxSelectError | None,
    metadata: TargetGene,
    urn: str,
    vrs_version: VrsVersion = VrsVersion.V_2,
) -> ScoreAnnotation:
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
        # Members share one reference; resolve the reconstruction accession once.
        accession = _resolve_postmapped_accession(
            mapped_score.alignment_level, post_mapped.members[0], metadata, tx_results
        )

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

            if allele.expressions:
                continue
            if accession is None:
                mapped_score.error_message = "Could not determine accession for this annotation. No allele expression is available."
            else:
                hgvs, syntax = _get_hgvs_string(allele, accession)
                allele.expressions = [Expression(syntax=syntax, value=hgvs)]

    if vrs_version == VrsVersion.V_1_3:
        pre_mapped = _haplotype_to_haplotype_1_3(pre_mapped)
        post_mapped = _haplotype_to_haplotype_1_3(post_mapped) if post_mapped else None

    return ScoreAnnotation(
        pre_mapped=pre_mapped,
        post_mapped=post_mapped,
        vrs_version=vrs_version,
        mavedb_id=mapped_score.accession_id,
        score=float(mapped_score.score) if mapped_score.score is not None else None,
        error_message=mapped_score.error_message,
        alignment_level=mapped_score.alignment_level,
        outcome=_resolve_outcome(mapped_score),
    )


def annotate(
    mapped_scores: list[MappedScore],
    tx_results: TxSelectResult | TxSelectError | None,
    metadata: TargetGene,
    urn: str,
    vrs_version: VrsVersion = VrsVersion.V_2,
) -> list[ScoreAnnotation]:
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
                ScoreAnnotation(
                    mavedb_id=mapped_score.accession_id,
                    score=float(mapped_score.score)
                    if mapped_score.score is not None
                    else None,
                    vrs_version=vrs_version,
                    error_message=mapped_score.error_message,
                    alignment_level=mapped_score.alignment_level,
                    outcome=_resolve_outcome(mapped_score),
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
                ScoreAnnotation(
                    pre_mapped=mapped_score.pre_mapped,
                    post_mapped=mapped_score.post_mapped,
                    vrs_version=vrs_version,
                    mavedb_id=mapped_score.accession_id,
                    score=float(mapped_score.score)
                    if mapped_score.score is not None
                    else None,
                    error_message=f"Multiple issues with annotation: Inconsistent variant structure (Allele and Haplotype mix).{' ' + mapped_score.error_message if mapped_score.error_message else ''}",
                    alignment_level=mapped_score.alignment_level,
                    outcome=_resolve_outcome(mapped_score),
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
    if layer == AnnotationLayer.CDNA:
        nm_accession: str | None = None
        if isinstance(tx_output, TxSelectResult) and tx_output.nm:
            nm_accession = tx_output.nm
        elif metadata.target_accession_id and metadata.target_accession_id.startswith(
            ("NM", "ENST")
        ):
            nm_accession = metadata.target_accession_id
        if nm_accession is None:
            return None
        vrs_id = get_vrs_id_from_identifier(nm_accession)
        if vrs_id is None:
            return None
        return MappedReferenceSequence(
            sequence_type=TargetSequenceType.DNA,
            sequence_id=vrs_id,
            sequence_accessions=[nm_accession],
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


def _align_result_for_target(
    target_gene: str,
    metadata: ScoresetMetadata,
    align_results: dict[str, "AlignmentResult"],
) -> "AlignmentResult | None":
    """Return the AlignmentResult for a target, handling both sequence-based and
    accession-based targets.

    Sequence-based targets key ``align_results`` by target gene name.
    Accession-based targets key it by ``target_accession_id`` (set by
    ``parse_cdot_mapping`` in ``align.py``). This helper does the fallback so
    callers don't have to repeat the lookup logic.
    """
    result = align_results.get(target_gene)
    if result is None:
        accession_id = metadata.target_genes[target_gene].target_accession_id
        if accession_id:
            result = align_results.get(accession_id)

    return result


def _pick_preferred_layer(
    target_meta: TargetGene,
    mappings: list[ScoreAnnotation],
) -> AnnotationLayer:
    """Return the preferred (assay) layer: the least-derived standard frame the target's
    input lets us assert.

    A variant is most faithfully represented in the coordinate system its assay natively
    *described* it in. Every other layer is a projection -- a transform that can lose or
    distort information -- so the assay layer is the one record we did not have to derive,
    and the one carried as preferred (``preferred_layer_only=True`` keeps exactly this
    layer and suppresses the projected forms). Two qualifiers sharpen "least-derived":

    * *standard frame* -- a target's bespoke sequence coordinates are not a shareable
      reference, so they never count; we express against g./c./p. references only.
    * *we can assert* -- the frame must actually be reachable for this target; an
      unreachable native frame falls through to the next-best reachable one.

    Accession-based targets name their frame in the accession, and it is always assertable
    (the accession *is* the reference), so selection is deterministic:

    * ``NP_``/``ENSP`` -> PROTEIN  (variants described as protein consequences)
    * ``NM_``/``ENST`` -> CDNA     (variants already coding on the transcript)
    * ``NC_``          -> GENOMIC  (variants described against the contig)

    Sequence-based targets have no standard native frame, so we fall to the nearest one
    reachable by alignment: genomic -- the universal anchor -- when the genome alignment
    succeeded (the common case, hence the reachability check against ``mappings``), else
    protein (e.g. a construct that does not place on the genome but aligns to a reference
    protein). This is why ``mappings`` is consulted only here, never for accessions.
    """
    accession = target_meta.target_accession_id
    if accession is not None:
        if accession.startswith(("NP", "ENSP")):
            return AnnotationLayer.PROTEIN
        if accession.startswith(("NM", "ENST")):
            return AnnotationLayer.CDNA
        if accession.startswith("NC"):
            return AnnotationLayer.GENOMIC

    # Sequence-based: prefer the genomic anchor when it was reachable, else protein.
    for mapping in mappings:
        if mapping.alignment_level == AnnotationLayer.GENOMIC:
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
            scoreset_mapping.model_dump(exclude_none=True),
            file,
            indent=4,
        )
    return output_path


def _resolve_protein_accession(
    target_metadata: TargetGene,
    tx_result: TxSelectResult | TxSelectError | None,
    _align_result: AlignmentResult
    | None,  # protein accession comes from tx_result only
) -> str | None:
    if isinstance(tx_result, TxSelectResult) and tx_result.np:
        return tx_result.np

    accession_id = target_metadata.target_accession_id
    if accession_id and accession_id.startswith(("NP", "ENSP")):
        return accession_id

    return None


def _resolve_cdna_accession(
    target_metadata: TargetGene,
    tx_result: TxSelectResult | TxSelectError | None,
    _align_result: AlignmentResult | None,  # cDNA accession comes from tx_result only
) -> str | None:
    if isinstance(tx_result, TxSelectResult) and tx_result.nm:
        return tx_result.nm

    accession_id = target_metadata.target_accession_id
    if accession_id and accession_id.startswith(("NM", "ENST")):
        return accession_id

    return None


def _resolve_genomic_accession(
    target_metadata: TargetGene,
    _tx_result: TxSelectResult
    | TxSelectError
    | None,  # genomic accession comes from align_result only
    align_result: AlignmentResult | None,
) -> str | None:
    if align_result is not None and align_result.chrom:
        try:
            chrom_id = get_chromosome_identifier(align_result.chrom)
        except Exception:
            _logger.exception(
                "Could not resolve chromosome identifier for %s", align_result.chrom
            )
            return None

        if chrom_id and chrom_id.startswith("refseq:"):
            return chrom_id[len("refseq:") :]

        return chrom_id

    accession_id = target_metadata.target_accession_id
    if accession_id and accession_id.startswith("NC"):
        return accession_id

    return None


_LAYER_ACCESSION_RESOLVERS = {
    AnnotationLayer.PROTEIN: _resolve_protein_accession,
    AnnotationLayer.CDNA: _resolve_cdna_accession,
    AnnotationLayer.GENOMIC: _resolve_genomic_accession,
}


def _reference_accession_for_target_level(
    alignment_level: AnnotationLayer,
    target_metadata: TargetGene,
    tx_result: TxSelectResult | TxSelectError | None,
    align_result: AlignmentResult | None,
) -> str | None:
    """Pick the most informative RefSeq-style accession describing the reference
    sequence that variants at ``alignment_level`` were mapped against.

    Falls back to None when nothing sensible is available; the API column is nullable.
    """
    resolver = _LAYER_ACCESSION_RESOLVERS.get(alignment_level)
    if resolver is None:
        return None

    return resolver(target_metadata, tx_result, align_result)


def _build_target_mapping(
    target_gene_identifier: str,
    target_meta: TargetGene,
    alignment_level: AnnotationLayer,
    preferred: bool,
    tx_result: TxSelectResult | TxSelectError | None,
    align_result: AlignmentResult | None,
    vrs_version: VrsVersion,
    annotations: list[ScoreAnnotation],
    protein_align_result: AlignmentResult | None = None,
    near_gap_window: int = NEAR_GAP_WINDOW,
) -> TargetMapping:
    """Assemble one target_mappings[] row for a (target, alignment_level) pair.

    Provenance and aggregate counts are populated from the data the pipeline
    already has. Alignment QC columns are populated from
    ``align_result.alignment_qc`` for GENOMIC rows or
    ``protein_align_result.alignment_qc`` for PROTEIN rows.

    CDNA rows are scored here for cdna-source (``NM_``/``ENST``) targets and, in
    all-layers mode, for projected cdna. They carry no alignment QC -- there is no
    cdna-frame alignment, so ``qc_source`` stays ``None``. Reference-only cdna/protein
    layers (no variants) come from ``_build_identity_target_mapping`` instead.
    """
    reference_accession = _reference_accession_for_target_level(
        alignment_level, target_meta, tx_result, align_result
    )
    reference_sequence_id: str | None = None
    if reference_accession is not None:
        try:
            reference_sequence_id = get_vrs_id_from_identifier(reference_accession)
        except Exception:
            _logger.exception(
                "Could not resolve VRS sequence id for %s", reference_accession
            )

    # Pick the alignment that lives in this row's coordinate frame:
    # GENOMIC -> BLAT genomic alignment; PROTEIN -> target-protein-to-reference
    # alignment from vrs_map. CDNA has no own-frame alignment, so qc_source stays None.
    qc_source: AlignmentResult | None = None
    if alignment_level == AnnotationLayer.GENOMIC:
        qc_source = align_result
    elif alignment_level == AnnotationLayer.PROTEIN:
        qc_source = protein_align_result

    tool_parameters: dict[str, Any] | None
    # Accession-based targets split into two cases:
    #   NC_ (chromosome/contig) — no alignment tool was run; no tool_parameters.
    #   NM_/ENST (transcript)   — cdot was used for transcript placement; record its parameters.
    # For non-accession-based targets we pull aligner parameters from the BLAT alignment as usual.
    accession_id = target_meta.target_accession_id
    if accession_id is not None and accession_id.startswith("NC"):
        # Direct chromosome/contig accession: the coordinate frame is defined by
        # the accession itself — no placement tool was invoked.
        tool_parameters = {
            "aligner": "reference_accession_passthrough",
        }
    elif accession_id is not None:
        # Transcript accession (NM_/ENST): cdot was used for transcript placement.
        # cdot_data_version is present in aligner_parameters (possibly None) for cdot
        # paths after parse_cdot_mapping, so we can read it unconditionally.
        cdot_data_version = (
            align_result.aligner_parameters.get("cdot_data_version")
            if align_result is not None and align_result.aligner_parameters is not None
            else None
        )
        tool_parameters = {
            "aligner": "cdot_transcript_placement",
            "aligner_version": cdot.__version__,
            "cdot_url": CDOT_URL,
            # Always emit the key; None means "cdot run, version unknown".
            "cdot_data_version": cdot_data_version,
        }
    else:
        tool_parameters = (
            qc_source.aligner_parameters if qc_source is not None else None
        )

    qc = qc_source.alignment_qc if qc_source is not None else None
    percent_identity = qc.percent_identity if qc is not None else None
    alignment_length = qc.alignment_length if qc is not None else None
    mismatch_count = qc.mismatch_count if qc is not None else None
    gap_count = qc.gap_count if qc is not None else None
    alignment_string = qc.alignment_string if qc is not None else None

    # Record the near-gap window alongside the CIGAR so consumers can interpret
    # the per-variant ``near_gap`` flag. Stamped on every layer that had its
    # flags evaluated against an alignment (any layer with a non-None qc).
    alignment_metadata: dict[str, Any] | None = None
    if qc is not None and qc.cigar:
        alignment_metadata = {"cigar": qc.cigar}
    if qc is not None:
        alignment_metadata = alignment_metadata or {}
        alignment_metadata["near_gap_window"] = near_gap_window
        # Always emit this bit so consumers never have to infer it from absence.
        alignment_metadata[
            "at_mismatched_locus_evaluated"
        ] = not qc.mismatch_positions_unavailable

    alignment_score = qc_source.score if qc_source is not None else None
    next_best_alignment_score = (
        qc_source.next_best_score if qc_source is not None else None
    )

    total = len(annotations)
    failed = sum(1 for a in annotations if a.post_mapped is None)
    clean = total - failed

    # near_gap is always evaluable (derived from gap positions in the alignment,
    # not per-base sequence content), so count it unconditionally.
    # at_mismatched_locus is only reliable when mismatch_positions_unavailable=False;
    # when pslx tseq/qseq data was absent we leave it as None (falsy) and omit it
    # from the tally so we don't silently under-count. The schema's
    # at_mismatched_locus_evaluated flag in alignment_metadata lets consumers
    # know they may be seeing only the near_gap portion of this count.
    mismatch_unavailable = qc is not None and qc.mismatch_positions_unavailable
    alignment_warnings = sum(
        1
        for a in annotations
        if a.near_gap or (not mismatch_unavailable and a.at_mismatched_locus)
    )

    return TargetMapping(
        target_gene_identifier=target_gene_identifier,
        alignment_level=alignment_level,
        preferred=preferred,
        tool_name="dcd-mapping",
        tool_version=dcd_mapping_version,
        tool_parameters=tool_parameters,
        reference_assembly=align_result.reference_assembly
        if align_result is not None
        else None,
        reference_accession=reference_accession,
        reference_sequence_id=reference_sequence_id,
        vrs_version=vrs_version,
        percent_identity=percent_identity,
        alignment_score=alignment_score,
        next_best_alignment_score=next_best_alignment_score,
        alignment_length=alignment_length,
        mismatch_count=mismatch_count,
        gap_count=gap_count,
        alignment_string=alignment_string,
        alignment_metadata=alignment_metadata,
        total_variants=total,
        variants_mapped_cleanly=clean,
        variants_with_alignment_warnings=alignment_warnings,
        variants_failed=failed,
    )


def _build_identity_target_mapping(
    target_gene_identifier: str,
    alignment_level: AnnotationLayer,
    reference_accession: str,
    vrs_version: VrsVersion,
) -> TargetMapping:
    """Assemble an identity ``target_mappings[]`` row for a coding layer that
    produced no per-variant mappings this run.

    A coding target's coding transcript (and protein reference) is known/selected
    by the mapper even when it does not emit per-variant cdna/protein mappings --
    the genomic-accession (``NC_``) path projects onto a gene-selected MANE
    transcript, and sequence-based / ``NM_`` targets carry the transcript trivially.
    This row surfaces that reference identity as target metadata so the API/RT can
    resolve the projection transcript without re-deriving it. It carries **null QC
    and null counts**: no ``mapped_score`` joins it (the
    ``mapped_score -> TargetGeneMapping`` join guarantee applies only to scored
    layers), so there is nothing to tally.
    """
    reference_sequence_id: str | None = None
    try:
        reference_sequence_id = get_vrs_id_from_identifier(reference_accession)
    except Exception:
        _logger.exception(
            "Could not resolve VRS sequence id for %s", reference_accession
        )

    return TargetMapping(
        target_gene_identifier=target_gene_identifier,
        alignment_level=alignment_level,
        preferred=False,
        tool_name="dcd-mapping",
        tool_version=dcd_mapping_version,
        tool_parameters={"aligner": "transcript_identity"},
        reference_accession=reference_accession,
        reference_sequence_id=reference_sequence_id,
        vrs_version=vrs_version,
    )


def build_scoreset_mapping(
    metadata: ScoresetMetadata,
    raw_metadata: dict,
    mappings: dict[str, list[ScoreAnnotation]],
    align_results: dict[str, AlignmentResult],
    tx_output: dict[str, TxSelectResult | TxSelectError | None],
    gene_info: dict[str, GeneInfo | None],
    preferred_layer_only: bool = False,
    vrs_version: VrsVersion = VrsVersion.V_2,
    protein_align_results: dict[str, AlignmentResult | None] | None = None,
    near_gap_window: int = NEAR_GAP_WINDOW,
) -> ScoresetMapping:
    """Assemble and return a :class:`ScoresetMapping` without writing to disk.

    :param metadata: parsed scoreset metadata (used for pipeline logic)
    :param raw_metadata: the full MaveDB API response dict, stored verbatim in
        the output ``metadata`` field so no fields are lost
    :param mappings: annotated VRS mappings per target. **Mutated in place** —
        ``_stamp_alignment_locus_flags`` sets ``at_mismatched_locus`` and
        ``near_gap`` on each ``ScoreAnnotation``. Do not pass the same dict to
        two successive calls if you need reproducible per-call results; deep-copy
        it first.
    :param align_results: alignment results per target (genomic frame, from BLAT)
    :param tx_output: transcript selection results per target
    :param gene_info: gene info per target
    :param preferred_layer_only: if True, only emit the preferred annotation layer
    :param vrs_version: VRS schema version for this run
    :param protein_align_results: per-target target-protein-to-reference-protein
        alignments (from vrs_map). Used to flag protein-layer variants at
        mismatched/near-gap loci and as the QC source for PROTEIN target_mappings
        rows. Pass ``None`` if you don't have them.
    :param near_gap_window: half-width (in reference bases) of the window used to
        flag variants as ``near_gap``. Defaults to :data:`NEAR_GAP_WINDOW` (5).
        Increase for a more conservative call set; decrease to restrict to
        immediately adjacent variants.
    :return: fully assembled ScoresetMapping
    """
    run_mapped_date_utc = datetime.datetime.now(tz=datetime.UTC)
    # set preferred layers for each target, to allow a mix of coding and noncoding targets
    reference_sequences: dict[str, TargetAnnotation] = {}
    mapped_scores: list[ScoreAnnotation] = []
    target_mappings: list[TargetMapping] = []
    for target_gene in mappings:
        # Stamp per-variant locus flags first so the warning counts in
        # TargetMapping and the emitted ScoreAnnotations are consistent.
        protein_align_for_target = (protein_align_results or {}).get(target_gene)
        genomic_align_for_target = _align_result_for_target(
            target_gene, metadata, align_results
        )
        _stamp_alignment_locus_flags(
            mappings[target_gene],
            genomic_align_for_target,
            AnnotationLayer.GENOMIC,
            near_gap_window,
        )
        _stamp_alignment_locus_flags(
            mappings[target_gene],
            protein_align_for_target,
            AnnotationLayer.PROTEIN,
            near_gap_window,
        )

        preferred_layer_for_target = _pick_preferred_layer(
            metadata.target_genes[target_gene], mappings[target_gene]
        )
        if preferred_layer_only:
            preferred_layers = {preferred_layer_for_target}
        else:
            preferred_layers = {
                mapping.alignment_level for mapping in mappings[target_gene]
            }

        # use target gene name in reference sequence dictionary, rather than the label, which differs between score sets
        target_gene_name = metadata.target_genes[target_gene].target_gene_name

        # sometimes Nonetype layers show up in preferred layers dict; remove these
        # before constructing TargetAnnotation (which requires valid AnnotationLayer keys).
        preferred_layers.discard(None)

        _logger.info(
            "For target %s, preferred layer is %s and layers seen are %s",
            target_gene_name,
            preferred_layer_for_target,
            preferred_layers,
        )

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
                genomic_align_for_target,
            )

        # If genomic layer and coding target, add a cdna entry (just the sequence
        # accession) to the reference_sequences dict. Covers both sequence-based
        # targets and genomic-accession (NC_) coding targets -- both carry the
        # selected coding transcript on tx_output[...].nm. NM_/ENST and NP_ accession
        # targets do not qualify (no TxSelectResult with nm), so the tx checks below
        # are the gate rather than an explicit accession-prefix test.
        if (
            AnnotationLayer.GENOMIC in reference_sequences[target_gene_name].layers
            and metadata.target_genes[target_gene].target_gene_category
            == TargetType.PROTEIN_CODING
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

        # Every input variant must yield exactly one record at the preferred layer.
        # Records already at a preferred layer are emitted directly. A variant with none
        # gets a single re-attributed failure there: its own null-layer failure if it has
        # one, else a synthesized failure (it mapped only at non-preferred layers -- e.g. a
        # wild-type p.= on a genomic-preferred target). A variant already represented is
        # never also re-attributed, so a dead attempt can't duplicate its real record.
        preferred_mappings = [
            m for m in mappings[target_gene] if m.alignment_level in preferred_layers
        ]
        represented_ids = {m.mavedb_id for m in preferred_mappings}

        records_by_id: dict[str, list[ScoreAnnotation]] = {}
        for m in mappings[target_gene]:
            records_by_id.setdefault(m.mavedb_id, []).append(m)

        reattributed: list[ScoreAnnotation] = []
        for variant_id, recs in records_by_id.items():
            if variant_id in represented_ids:
                continue
            null_failure = next(
                (m for m in recs if m.alignment_level is None and m.pre_mapped is None),
                None,
            )
            if null_failure is not None:
                score_dict = null_failure.model_dump()
            else:
                # Mapped only at non-preferred layers; synthesize a failure there.
                mapped_layers = sorted(
                    {
                        m.alignment_level.value
                        for m in recs
                        if m.alignment_level is not None
                    }
                )
                score_dict = recs[0].model_dump()
                score_dict["pre_mapped"] = None
                score_dict["post_mapped"] = None
                score_dict["outcome"] = MappingOutcome.FAILED
                score_dict["error_message"] = (
                    f"No representation at preferred layer {preferred_layer_for_target.value}"
                    + (
                        f"; mapped only at: {', '.join(mapped_layers)}"
                        if mapped_layers
                        else ""
                    )
                )
            score_dict["alignment_level"] = preferred_layer_for_target
            score_dict["target_gene_identifier"] = target_gene_name
            reattributed.append(ScoreAnnotation(**score_dict))

        for m in preferred_mappings:
            score_dict = m.model_dump()
            score_dict["target_gene_identifier"] = target_gene_name
            mapped_scores.append(ScoreAnnotation(**score_dict))
        mapped_scores.extend(reattributed)

        # One TargetMapping per preferred layer that produced records. The
        # (target_gene_identifier, alignment_level) key must be unique per run and
        # cover every mapped_score, so the preferred layer also gets a row whenever
        # failures were re-attributed to it. Re-attributed failures count toward its totals.
        layers_seen: dict[AnnotationLayer, list[ScoreAnnotation]] = {}
        for m in preferred_mappings:
            layers_seen.setdefault(m.alignment_level, []).append(m)
        if reattributed:
            layers_seen.setdefault(preferred_layer_for_target, [])

        for layer, layer_annotations in layers_seen.items():
            annotations_for_tm = (
                layer_annotations + reattributed
                if layer == preferred_layer_for_target
                else layer_annotations
            )
            target_mappings.append(
                _build_target_mapping(
                    target_gene_identifier=target_gene_name,
                    target_meta=metadata.target_genes[target_gene],
                    alignment_level=layer,
                    preferred=(layer == preferred_layer_for_target),
                    tx_result=tx_output.get(target_gene),
                    align_result=genomic_align_for_target,
                    vrs_version=vrs_version,
                    annotations=annotations_for_tm,
                    protein_align_result=protein_align_for_target,
                    near_gap_window=near_gap_window,
                )
            )

        # Identity target_mappings for coding layers that produced no per-variant
        # mappings this run. A coding target's coding transcript (and protein
        # reference) may be known to the mapper even when it emits no per-variant
        # cdna/protein mappings: the genomic-accession (NC_) path projects onto a
        # gene-selected MANE transcript; sequence-based and NM_/ENST targets carry
        # it trivially. Surfacing it as a cdna (and protein) identity row lets consumers
        # prefer the mapper-selected transcript over the NP_->NM_ fallback.
        # Decoupled from layers_seen on purpose -- the cdna/protein per-variant
        # forms are filtered, so those layers never appear in layers_seen.
        target_meta = metadata.target_genes[target_gene]
        if target_meta.target_gene_category == TargetType.PROTEIN_CODING:
            scored_levels = {m.alignment_level for m in preferred_mappings}
            scored_levels.discard(None)
            for identity_level in (AnnotationLayer.CDNA, AnnotationLayer.PROTEIN):
                if identity_level in scored_levels:
                    continue
                identity_accession = _reference_accession_for_target_level(
                    identity_level,
                    target_meta,
                    tx_output.get(target_gene),
                    genomic_align_for_target,
                )
                if identity_accession is None:
                    continue
                target_mappings.append(
                    _build_identity_target_mapping(
                        target_gene_identifier=target_gene_name,
                        alignment_level=identity_level,
                        reference_accession=identity_accession,
                        vrs_version=vrs_version,
                    )
                )

    # Drop layers where both reference sequence entries are None, and any None-keyed
    # layers. Moved outside the per-target loop to avoid O(n²) scans and eliminate
    # the variable-shadowing hazard (the inner loop previously reused `target_gene`).
    for tgt_name in reference_sequences:
        for layer in list(reference_sequences[tgt_name].layers.keys()):
            if (
                reference_sequences[tgt_name].layers[layer]["mapped_reference_sequence"]
                is None
                and reference_sequences[tgt_name].layers[layer][
                    "computed_reference_sequence"
                ]
                is None
            ) or layer is None:
                del reference_sequences[tgt_name].layers[layer]

    # Runtime invariant: every mapped_score must be attributable to a TargetMapping
    # row via (target_gene_identifier, alignment_level). If any orphaned keys exist
    # the API join would silently drop those scores, undermining the audit guarantee.
    # Using an explicit RuntimeError rather than `assert` so this is never compiled
    # away with `-O` optimizations in production deployments.
    tm_keys = {
        (tm.target_gene_identifier, tm.alignment_level) for tm in target_mappings
    }
    orphaned = {
        (ms.target_gene_identifier, ms.alignment_level)
        for ms in mapped_scores
        if ms.alignment_level is not None
    } - tm_keys
    if orphaned:
        raise RuntimeError(
            f"build_scoreset_mapping produced {len(orphaned)} orphaned mapped_score "  # noqa: EM102
            f"key(s) with no corresponding TargetMapping row: {orphaned!r}. "
            "This is a pipeline invariant violation — every mapped_score must have "
            "a parent TargetMapping so the API join is unambiguous."
        )

    return ScoresetMapping(
        metadata=raw_metadata,
        mapped_date=run_mapped_date_utc,
        reference_sequences=reference_sequences,
        mapped_scores=mapped_scores,
        target_mappings=target_mappings,
    )


def save_mapped_output_json(
    metadata: ScoresetMetadata,
    raw_metadata: dict,
    mappings: dict[str, list[ScoreAnnotation]],
    align_results: dict[str, AlignmentResult],
    tx_output: dict[str, TxSelectResult | TxSelectError | None],
    gene_info: dict[str, GeneInfo | None],
    preferred_layer_only: bool = False,
    output_path: Path | None = None,
    vrs_version: VrsVersion = VrsVersion.V_2,
    protein_align_results: dict[str, AlignmentResult | None] | None = None,
) -> Path:
    """Build and write mapping output for a score set to a JSON file.

    :return: output file path
    """
    output = build_scoreset_mapping(
        metadata=metadata,
        raw_metadata=raw_metadata,
        mappings=mappings,
        align_results=align_results,
        tx_output=tx_output,
        gene_info=gene_info,
        preferred_layer_only=preferred_layer_only,
        vrs_version=vrs_version,
        protein_align_results=protein_align_results,
    )
    return write_scoreset_mapping_to_json(metadata.urn, output, output_path)
