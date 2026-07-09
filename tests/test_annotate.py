"""Tests for dcd_mapping.annotate"""

from unittest import mock

import pytest
from ga4gh.vrs._internal.models import (
    Allele,
    LiteralSequenceExpression,
    SequenceLocation,
    SequenceReference,
)

from dcd_mapping import vrs_v1_schemas
from dcd_mapping.annotate import (
    _compute_target_gene_info_from_alignment,
    _compute_target_gene_info_from_mapped_variant_spans,
    _covered_bases_from_overlapping_genes_of_chromosomal_intervals,
    _get_mapped_reference_sequence,
    _stamp_alignment_locus_flags,
    compute_target_gene_info,
)
from dcd_mapping.schemas import (
    AlignmentQc,
    AlignmentResult,
    AnnotationLayer,
    GeneInfo,
    MappedReferenceSequence,
    MappedScore,
    ScoreAnnotation,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)


@pytest.fixture
def target_dna_pc():
    return TargetGene(
        target_gene_name="BRAF",
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence="ATGGCG...",
        target_sequence_type=TargetSequenceType.DNA,
        target_accession_id=None,
        target_uniprot_ref=None,
    )


@pytest.fixture
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
        percent_identity=None,
        query_range=SequenceRange(start=1, end=10),
        query_subranges=[SequenceRange(start=1, end=10)],
        hit_range=SequenceRange(start=1, end=10),
        hit_subranges=[SequenceRange(start=s, end=e) for s, e in hit_intervals],
    )


@pytest.mark.asyncio
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


@pytest.mark.asyncio
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
        alignment_level=AnnotationLayer.GENOMIC,
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


@pytest.mark.asyncio
async def test_compute_target_gene_info_fallback_metadata(scoreset_metadata):
    # No tx, no alignment, no mapped scores -> fallback
    with mock.patch("dcd_mapping.annotate.get_gene_symbol", return_value="META"):
        res = await compute_target_gene_info(
            "label", {"label": None}, {"label": None}, scoreset_metadata, None
        )
        assert isinstance(res, GeneInfo)
        assert res.hgnc_symbol == "META"
        assert res.selection_method == "target_metadata"


@pytest.mark.asyncio
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
            alignment_level=AnnotationLayer.GENOMIC,
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


# ---------------------------------------------------------------------------
# Tests for _apply_alignment_locus_flags / mismatch_positions_unavailable
# ---------------------------------------------------------------------------


def _make_genomic_annotation(start: int, end: int) -> "ScoreAnnotation":
    """Return a ScoreAnnotation with a simple genomic post_mapped Allele."""
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
    return ScoreAnnotation(
        mavedb_id=f"id_{start}",
        pre_mapped=None,
        post_mapped=allele,
        alignment_level=AnnotationLayer.GENOMIC,
        score=None,
    )


def _make_align_result_with_qc(**qc_kwargs) -> "AlignmentResult":
    """Build a minimal AlignmentResult with an AlignmentQc built from qc_kwargs."""
    qc = AlignmentQc(
        alignment_length=100,
        mismatch_count=qc_kwargs.pop("mismatch_count", 0),
        gap_count=qc_kwargs.pop("gap_count", 0),
        **qc_kwargs,
    )
    return AlignmentResult(
        chrom="NC_000001.11",
        strand=1,
        coverage=None,
        percent_identity=None,
        query_range=SequenceRange(start=0, end=100),
        query_subranges=[SequenceRange(start=0, end=100)],
        hit_range=SequenceRange(start=0, end=100),
        hit_subranges=[SequenceRange(start=0, end=100)],
        alignment_qc=qc,
    )


def test_apply_alignment_locus_flags_mismatch_positions_unavailable_leaves_at_mismatched_locus_none():
    """When mismatch_positions_unavailable=True and mismatch_count>0, at_mismatched_locus
    must stay None (not evaluated) to prevent silent disagreement with mismatch_count.
    """
    ann = _make_genomic_annotation(10, 11)
    assert ann.at_mismatched_locus is None  # pre-condition

    align_result = _make_align_result_with_qc(
        mismatch_count=5,
        mismatch_positions_unavailable=True,
    )
    _stamp_alignment_locus_flags([ann], align_result, AnnotationLayer.GENOMIC)

    assert ann.at_mismatched_locus is None, (
        "at_mismatched_locus must remain None when mismatch_positions_unavailable=True "
        "and mismatch_count>0 -- it cannot claim False when positions are unknown"
    )
    assert ann.near_gap is False  # gap_intervals empty → near_gap still evaluated


def test_apply_alignment_locus_flags_mismatch_positions_unavailable_zero_count_stamps_false():
    """When mismatch_positions_unavailable=True but mismatch_count==0, there are no
    mismatches to miss, so at_mismatched_locus=False is safe.
    """
    ann = _make_genomic_annotation(10, 11)
    align_result = _make_align_result_with_qc(
        mismatch_count=0,
        mismatch_positions_unavailable=True,
    )
    _stamp_alignment_locus_flags([ann], align_result, AnnotationLayer.GENOMIC)

    assert ann.at_mismatched_locus is False
    assert ann.near_gap is False


def test_apply_alignment_locus_flags_early_exit_does_not_fire_when_count_nonzero():
    """The early-exit 'no mismatches and no gaps → stamp False' must NOT fire when
    mismatch_positions_unavailable=True and mismatch_count>0, even if gap_intervals
    is also empty.  Previously the check was ``not mismatch_positions``, which is
    always True under unavailability, causing the early-exit to stamp near_gap=False
    on annotations whose positions could not be extracted (and thus were never evaluated).
    The fix uses ``mismatch_count == 0`` as the ground-truth condition.
    """
    # A VRS v1 Allele is not an instance of ga4gh.vrs.Allele (v2), so
    # _allele_ref_positions returns [] and the main-loop body is skipped via `continue`.
    v1_allele = vrs_v1_schemas.Allele(
        id="ga4gh:VA.test",
        location=vrs_v1_schemas.SequenceLocation(
            id="loc1",
            sequence_id="ga4gh:SQ.test",
            interval=vrs_v1_schemas.SequenceInterval(
                start=vrs_v1_schemas.Number(value=10),
                end=vrs_v1_schemas.Number(value=11),
            ),
        ),
        state=vrs_v1_schemas.LiteralSequenceExpression(sequence="A"),
    )
    ann = ScoreAnnotation.model_construct(
        mavedb_id="id_v1",
        post_mapped=v1_allele,
        alignment_level=AnnotationLayer.GENOMIC,
        at_mismatched_locus=None,
        near_gap=None,
    )

    align_result = _make_align_result_with_qc(
        mismatch_count=3,
        mismatch_positions_unavailable=True,
    )
    _stamp_alignment_locus_flags([ann], align_result, AnnotationLayer.GENOMIC)

    # Both flags must remain None: positions could not be extracted, and the early-exit
    # must not have fired (it would have incorrectly stamped near_gap=False).
    assert ann.at_mismatched_locus is None
    assert ann.near_gap is None, (
        "near_gap must be None when positions are unextractable and mismatch_count>0 "
        "prevents the early-exit from firing -- old code stamped False via faulty early-exit"
    )


def _make_nm_target(accession_id: str | None = "NM_007294.3") -> TargetGene:
    return TargetGene(
        target_gene_name="BRCA1",
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence="ATGG",
        target_sequence_type=TargetSequenceType.DNA,
        target_accession_id=accession_id,
        target_uniprot_ref=None,
    )


def _make_tx_result(nm: str = "NM_007294.3", np: str = "NP_009225.1") -> TxSelectResult:
    return TxSelectResult(
        nm=nm,
        np=np,
        start=0,
        is_full_match=True,
        sequence="MAST",
    )


def _make_genomic_align(chrom: str = "NC_000017.11") -> AlignmentResult:
    return AlignmentResult(
        chrom=chrom,
        strand=1,
        coverage=None,
        percent_identity=None,
        query_range=SequenceRange(start=1, end=10),
        query_subranges=[SequenceRange(start=1, end=10)],
        hit_range=SequenceRange(start=1, end=10),
        hit_subranges=[SequenceRange(start=1, end=10)],
    )


class TestGetMappedReferenceSequence:
    """_get_mapped_reference_sequence must return the NM/ENST transcript for the CDNA
    layer -- never the genomic chromosome -- regardless of whether an align_result with
    a chromosome is also available.
    """

    def test_cdna_layer_nm_accession_target_uses_nm(self):
        """NM_ accession target: the NM_ accession itself is the cdna mapped reference."""
        target = _make_nm_target("NM_007294.3")
        align = _make_genomic_align("NC_000017.11")
        vrs_id = "ga4gh:SQ.fake_nm"
        with mock.patch(
            "dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=vrs_id
        ):
            result = _get_mapped_reference_sequence(
                target, AnnotationLayer.CDNA, None, align
            )
        assert isinstance(result, MappedReferenceSequence)
        assert result.sequence_accessions == ["NM_007294.3"]
        assert result.sequence_id == vrs_id
        assert result.sequence_type == TargetSequenceType.DNA

    def test_cdna_layer_tx_result_nm_takes_priority_over_accession(self):
        """When a TxSelectResult carries an nm, it takes precedence over the target accession."""
        target = _make_nm_target("NM_007294.3")
        tx = _make_tx_result(nm="NM_007294.4")  # versioned upgrade
        align = _make_genomic_align("NC_000017.11")
        vrs_id = "ga4gh:SQ.fake_nm_v4"
        with mock.patch(
            "dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=vrs_id
        ):
            result = _get_mapped_reference_sequence(
                target, AnnotationLayer.CDNA, tx, align
            )
        assert isinstance(result, MappedReferenceSequence)
        assert result.sequence_accessions == ["NM_007294.4"]

    def test_cdna_layer_nc_accession_target_uses_tx_nm(self):
        """NC_ accession target: the cdna mapped reference is the MANE transcript from tx_result."""
        target = _make_nm_target("NC_000017.11")
        target = TargetGene(
            target_gene_name="BRCA1",
            target_gene_category=TargetType.PROTEIN_CODING,
            target_sequence="ATGG",
            target_sequence_type=TargetSequenceType.DNA,
            target_accession_id="NC_000017.11",
            target_uniprot_ref=None,
        )
        tx = _make_tx_result(nm="NM_007294.3")
        align = _make_genomic_align("NC_000017.11")
        vrs_id = "ga4gh:SQ.fake_nm"
        with mock.patch(
            "dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=vrs_id
        ):
            result = _get_mapped_reference_sequence(
                target, AnnotationLayer.CDNA, tx, align
            )
        assert isinstance(result, MappedReferenceSequence)
        assert result.sequence_accessions == ["NM_007294.3"]

    def test_cdna_layer_no_nm_returns_none(self):
        """When no NM is resolvable for the CDNA layer, return None rather than a chromosome."""
        target = TargetGene(
            target_gene_name="BRCA1",
            target_gene_category=TargetType.PROTEIN_CODING,
            target_sequence="ATGG",
            target_sequence_type=TargetSequenceType.DNA,
            target_accession_id=None,
            target_uniprot_ref=None,
        )
        align = _make_genomic_align("NC_000017.11")
        result = _get_mapped_reference_sequence(
            target, AnnotationLayer.CDNA, None, align
        )
        assert result is None

    def test_genomic_layer_returns_chromosome(self):
        """GENOMIC layer: existing behaviour -- returns the chromosome from align_result.chrom."""
        target = _make_nm_target(None)
        align = _make_genomic_align("NC_000017.11")
        vrs_id = "ga4gh:SQ.fake_nc"
        with (
            mock.patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000017.11",
            ),
            mock.patch(
                "dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=vrs_id
            ),
        ):
            result = _get_mapped_reference_sequence(
                target, AnnotationLayer.GENOMIC, None, align
            )
        assert isinstance(result, MappedReferenceSequence)
        assert result.sequence_accessions == ["NC_000017.11"]

    def test_protein_layer_returns_np(self):
        """PROTEIN layer: existing behaviour -- returns the NP_ accession from tx_result."""
        target = _make_nm_target(None)
        tx = _make_tx_result(np="NP_009225.1")
        vrs_id = "ga4gh:SQ.fake_np"
        with mock.patch(
            "dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=vrs_id
        ):
            result = _get_mapped_reference_sequence(
                target, AnnotationLayer.PROTEIN, tx, None
            )
        assert isinstance(result, MappedReferenceSequence)
        assert result.sequence_accessions == ["NP_009225.1"]
        assert result.sequence_type == TargetSequenceType.PROTEIN
