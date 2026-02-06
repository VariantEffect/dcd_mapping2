"""Tests for dcd_mapping.annotate"""
from unittest import mock

import pytest
from ga4gh.vrs._internal.models import (
    Allele,
    LiteralSequenceExpression,
    SequenceLocation,
    SequenceReference,
)

from dcd_mapping.annotate import (
    _compute_target_gene_info_from_alignment,
    _compute_target_gene_info_from_mapped_variant_spans,
    _covered_bases_from_overlapping_genes_of_chromosomal_intervals,
    compute_target_gene_info,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    AnnotationLayer,
    GeneInfo,
    MappedScore,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)


@pytest.fixture()
def target_dna_pc():
    return TargetGene(
        target_gene_name="BRAF",
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence="ATGGCG...",
        target_sequence_type=TargetSequenceType.DNA,
        target_accession_id=None,
        target_uniprot_ref=None,
    )


@pytest.fixture()
def scoreset_metadata(target_dna_pc):
    return ScoresetMetadata(
        urn="urn:mavedb:TEST",
        score_count=1,
        target_genes={"label": target_dna_pc},
        mapped=False,
    )


def make_align(hit_intervals):
    return AlignmentResult(
        chrom="NC_000001.11",
        strand=1,
        coverage=None,
        ident_pct=None,
        query_range=SequenceRange(start=1, end=10),
        query_subranges=[SequenceRange(start=1, end=10)],
        hit_range=SequenceRange(start=1, end=10),
        hit_subranges=[SequenceRange(start=s, end=e) for s, e in hit_intervals],
    )


@pytest.mark.asyncio()
async def test_compute_target_gene_info_non_coding_category():
    meta = ScoresetMetadata(
        urn="urn:mavedb:TEST",
        score_count=0,
        target_genes={
            "t": TargetGene(
                target_gene_name="REG",
                target_gene_category=TargetType.OTHER_NC,
                target_sequence="ACGT",
                target_sequence_type=TargetSequenceType.DNA,
            )
        },
        mapped=False,
    )

    res = await compute_target_gene_info("t", {}, {}, meta, None)
    assert isinstance(res, GeneInfo)
    assert res.hgnc_symbol is None
    assert res.selection_method == "target_category"


@pytest.mark.asyncio()
async def test_compute_target_gene_info_tx_selection(scoreset_metadata):
    tx = TxSelectResult(
        nm="NM_000001.1",
        np="NP_000001.1",
        start=0,
        is_full_match=True,
        sequence="MSEQUENCE",
        transcript_mode=None,
        hgnc_symbol="BRAF",
    )
    res = await compute_target_gene_info(
        "label", {"label": tx}, {"label": None}, scoreset_metadata, None
    )
    assert isinstance(res, GeneInfo)
    assert res.hgnc_symbol == "BRAF"
    assert res.selection_method == "tx_selection"


def test_compute_target_gene_info_alignment_path(scoreset_metadata):
    align = make_align([(100, 120)])
    with (
        mock.patch(
            "dcd_mapping.annotate.get_chromosome_identifier",
            side_effect=lambda c: c,
        ),
        mock.patch(
            "dcd_mapping.annotate.get_overlapping_features_for_region",
            return_value=[{"external_name": "GENE1", "start": 95, "end": 130}],
        ),
        mock.patch("dcd_mapping.annotate._get_hgnc_symbol", side_effect=lambda s: s),
    ):
        res = _compute_target_gene_info_from_alignment("label", align)
        assert isinstance(res, GeneInfo)
        assert res.hgnc_symbol == "GENE1"
        assert res.selection_method == "alignment_max_covered_bases"


def test_compute_target_gene_info_alignment_tie_returns_none(scoreset_metadata):
    align = make_align([(100, 120)])
    with (
        mock.patch(
            "dcd_mapping.annotate.get_chromosome_identifier",
            side_effect=lambda c: c,
        ),
        mock.patch(
            "dcd_mapping.annotate.get_overlapping_features_for_region",
            return_value=[
                {"external_name": "GENE1", "start": 100, "end": 120},
                {"external_name": "GENE2", "start": 100, "end": 120},
            ],
        ),
    ):
        res = _compute_target_gene_info_from_alignment("label", align)
        assert res is None


def test_compute_target_gene_info_mapped_variants_path(scoreset_metadata):
    # Build mapped scores to yield spans that overlap single gene best
    allele = Allele(
        location=SequenceLocation(
            sequenceReference=SequenceReference(
                refgetAccession="SQ.1234567890abcdef1234567890abcdef"
            ),
            start=100,
            end=110,
        ),
        state=LiteralSequenceExpression(sequence="A"),
        expressions=[],
    )
    ms = MappedScore(
        mavedb_id="id",
        accession_id="id",
        pre_mapped=None,
        post_mapped=allele,
        annotation_layer=AnnotationLayer.GENOMIC,
        score=None,
        error_message=None,
    )
    with (
        mock.patch(
            "dcd_mapping.annotate.get_chromosome_identifier_from_vrs_id",
            return_value="refseq:NC_000001.11",
        ),
        mock.patch(
            "dcd_mapping.annotate.get_ucsc_chromosome_name",
            return_value="NC_000001.11",
        ),
        mock.patch("dcd_mapping.annotate._get_hgnc_symbol", side_effect=lambda s: s),
        mock.patch(
            "dcd_mapping.annotate.get_overlapping_features_for_region",
            return_value=[{"external_name": "GENE3", "start": 50, "end": 200}],
        ),
    ):
        res = _compute_target_gene_info_from_mapped_variant_spans("label", [ms])
        assert isinstance(res, GeneInfo)
        assert res.hgnc_symbol == "GENE3"
        assert res.selection_method == "variants_max_covered_bases"


@pytest.mark.asyncio()
async def test_compute_target_gene_info_fallback_metadata(scoreset_metadata):
    # No tx, no alignment, no mapped scores -> fallback
    with mock.patch("dcd_mapping.annotate.get_gene_symbol", return_value="META"):
        res = await compute_target_gene_info(
            "label", {"label": None}, {"label": None}, scoreset_metadata, None
        )
        assert isinstance(res, GeneInfo)
        assert res.hgnc_symbol == "META"
        assert res.selection_method == "target_metadata"


@pytest.mark.asyncio()
async def test_compute_target_gene_info_fallback_unavailable(scoreset_metadata):
    # No tx, no alignment, no mapped scores -> fallback
    with mock.patch("dcd_mapping.annotate.get_gene_symbol", return_value=None):
        res = await compute_target_gene_info(
            "label", {"label": None}, {"label": None}, scoreset_metadata, None
        )
        assert res is None


def test_covered_bases_sums_and_skips_invalid():
    with (
        mock.patch(
            "dcd_mapping.annotate.get_overlapping_features_for_region",
            return_value=[
                {"external_name": "A", "start": 10, "end": 30},
                {"external_name": None, "start": 10, "end": 30},
                {"external_name": "B", "start": None, "end": 100},
            ],
        ),
    ):
        cov = _covered_bases_from_overlapping_genes_of_chromosomal_intervals(
            "NC_000001.11", [(15, 25), (20, 35)]
        )
        # Overlaps: A: (15-25)=10 and (20-30)=10 => total 20
        assert cov == {"A": 20}


def test_interval_merging_overlapping_and_adjacent_merge(scoreset_metadata):
    # Two intervals on same chrom: overlapping and adjacent; both should merge before coverage

    def make_ms(start, end):
        allele = Allele(
            location=SequenceLocation(
                sequenceReference=SequenceReference(
                    refgetAccession="SQ.1234567890abcdef1234567890abcdef"
                ),
                start=start,
                end=end,
            ),
            state=LiteralSequenceExpression(sequence="A"),
            expressions=[],
        )
        return MappedScore(
            mavedb_id=f"id_{start}",
            accession_id=f"id_{start}",
            pre_mapped=None,
            post_mapped=allele,
            annotation_layer=AnnotationLayer.GENOMIC,
            score=None,
            error_message=None,
        )

    ms_list = [
        make_ms(100, 110),
        make_ms(108, 120),
        make_ms(120, 130),
    ]  # overlap then adjacent

    # After merging, intervals become [(100, 130)], expect single coverage call and length 30 contributing to GENEZ
    def fake_overlap(chrom, s, e, features=None):
        assert (s, e) == (100, 130)
        return [{"external_name": "GENEZ", "start": 90, "end": 200}]

    with (
        mock.patch(
            "dcd_mapping.annotate.get_chromosome_identifier_from_vrs_id",
            return_value="refseq:NC_000001.11",
        ),
        mock.patch(
            "dcd_mapping.annotate.get_ucsc_chromosome_name",
            return_value="NC_000001.11",
        ),
        mock.patch("dcd_mapping.annotate._get_hgnc_symbol", side_effect=lambda s: s),
        mock.patch(
            "dcd_mapping.annotate.get_overlapping_features_for_region",
            side_effect=fake_overlap,
        ),
    ):
        res = _compute_target_gene_info_from_mapped_variant_spans("label", ms_list)
        assert isinstance(res, GeneInfo)
        assert res.hgnc_symbol == "GENEZ"
        assert res.selection_method == "variants_max_covered_bases"
