"""Tests for annotate._reference_accession_for_target_level and build_scoreset_mapping."""

from unittest.mock import patch

from cool_seq_tool.schemas import AnnotationLayer

from dcd_mapping.annotate import (
    _align_result_for_target,
    _reference_accession_for_target_level,
    _stamp_alignment_locus_flags,
    build_scoreset_mapping,
)
from dcd_mapping.schemas import (
    AlignmentQc,
    AlignmentResult,
    GeneInfo,
    ScoreAnnotation,
    ScoresetMapping,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
    VrsVersion,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_seq_target(name: str = "GENE1") -> TargetGene:
    return TargetGene(
        target_gene_name=name,
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence="ATGCATGC",
        target_sequence_type=TargetSequenceType.DNA,
        target_accession_id=None,
    )


def _make_acc_target(accession_id: str, name: str = "GENE1") -> TargetGene:
    return TargetGene(
        target_gene_name=name,
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence=None,
        target_sequence_type=None,
        target_accession_id=accession_id,
    )


def _make_tx(nm: str = "NM_000001.1", np: str = "NP_000001.1") -> TxSelectResult:
    return TxSelectResult(
        nm=nm,
        np=np,
        start=0,
        is_full_match=True,
        sequence="MAAAA",
        hgnc_symbol="GENE1",
    )


def _make_align(chrom: str = "chr1") -> AlignmentResult:
    return AlignmentResult(
        chrom=chrom,
        query_range=SequenceRange(start=0, end=100),
        query_subranges=[SequenceRange(start=0, end=100)],
        hit_range=SequenceRange(start=1000, end=1100),
        hit_subranges=[SequenceRange(start=1000, end=1100)],
        aligner_parameters={"aligner": "blat", "min_score": 20, "out_format": "pslx"},
    )


# ---------------------------------------------------------------------------
# _reference_accession_for_target_level
# ---------------------------------------------------------------------------


class TestReferenceAccessionForTargetLevel:
    def test_protein_layer_from_tx(self):
        target = _make_seq_target()
        tx = _make_tx(np="NP_001234.1")
        result = _reference_accession_for_target_level(
            AnnotationLayer.PROTEIN, target, tx, None
        )
        assert result == "NP_001234.1"

    def test_protein_layer_from_accession(self):
        target = _make_acc_target("NP_001234.1")
        result = _reference_accession_for_target_level(
            AnnotationLayer.PROTEIN, target, None, None
        )
        assert result == "NP_001234.1"

    def test_protein_layer_ensp_accession(self):
        target = _make_acc_target("ENSP00000001.1")
        result = _reference_accession_for_target_level(
            AnnotationLayer.PROTEIN, target, None, None
        )
        assert result == "ENSP00000001.1"

    def test_protein_layer_no_info(self):
        target = _make_seq_target()
        result = _reference_accession_for_target_level(
            AnnotationLayer.PROTEIN, target, None, None
        )
        assert result is None

    def test_cdna_layer_from_tx(self):
        target = _make_seq_target()
        tx = _make_tx(nm="NM_000001.1")
        result = _reference_accession_for_target_level(
            AnnotationLayer.CDNA, target, tx, None
        )
        assert result == "NM_000001.1"

    def test_cdna_layer_from_accession(self):
        target = _make_acc_target("NM_000001.1")
        result = _reference_accession_for_target_level(
            AnnotationLayer.CDNA, target, None, None
        )
        assert result == "NM_000001.1"

    def test_cdna_layer_no_info(self):
        target = _make_seq_target()
        result = _reference_accession_for_target_level(
            AnnotationLayer.CDNA, target, None, None
        )
        assert result is None

    def test_genomic_layer_from_alignment(self):
        target = _make_seq_target()
        align = _make_align(chrom="chr1")
        with patch(
            "dcd_mapping.annotate.get_chromosome_identifier",
            return_value="refseq:NC_000001.11",
        ):
            result = _reference_accession_for_target_level(
                AnnotationLayer.GENOMIC, target, None, align
            )
        assert result == "NC_000001.11"

    def test_genomic_layer_nc_accession(self):
        target = _make_acc_target("NC_000001.11")
        result = _reference_accession_for_target_level(
            AnnotationLayer.GENOMIC, target, None, None
        )
        assert result == "NC_000001.11"

    def test_genomic_layer_no_info(self):
        target = _make_seq_target()
        result = _reference_accession_for_target_level(
            AnnotationLayer.GENOMIC, target, None, None
        )
        assert result is None


# ---------------------------------------------------------------------------
# build_scoreset_mapping golden test
# ---------------------------------------------------------------------------


def _make_annotation(
    layer: AnnotationLayer | None,
    post_mapped=None,
    error_message: str | None = None,
) -> ScoreAnnotation:
    """Minimal ScoreAnnotation for golden test."""
    return ScoreAnnotation(
        mavedb_id="urn:mavedb:00000001-a-1#1",
        alignment_level=layer,
        pre_mapped=None,
        post_mapped=post_mapped,
        error_message=error_message,
        score=None,
    )


class TestBuildScoresetMapping:
    def test_emits_one_target_mapping_per_layer(self):
        """Each (target, layer) pair appearing in mappings generates exactly one TargetMapping."""
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={
                "GENE1": _make_seq_target("GENE1"),
            },
        )
        # Two annotations: one genomic, one protein
        g_ann = _make_annotation(AnnotationLayer.GENOMIC)
        p_ann = _make_annotation(AnnotationLayer.PROTEIN)
        mappings = {"GENE1": [g_ann, p_ann]}

        with (
            patch(
                "dcd_mapping.annotate.get_vrs_id_from_identifier",
                return_value="ga4gh:SQ.test",
            ),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings=mappings,
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": GeneInfo(hgnc_symbol="GENE1")},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        assert isinstance(result, ScoresetMapping)
        assert result.target_mappings is not None
        layers_seen = {tm.alignment_level for tm in result.target_mappings}
        assert AnnotationLayer.GENOMIC in layers_seen
        assert AnnotationLayer.PROTEIN in layers_seen
        # Exactly one mapping per layer
        assert len(result.target_mappings) == 2

    def test_preferred_flag_on_exactly_one_row_per_target(self):
        """preferred=True must appear on exactly one TargetMapping per target."""
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        g_ann = _make_annotation(AnnotationLayer.GENOMIC)
        p_ann = _make_annotation(AnnotationLayer.PROTEIN)

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [g_ann, p_ann]},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        preferred = [tm for tm in result.target_mappings if tm.preferred]
        assert len(preferred) == 1
        assert preferred[0].alignment_level == AnnotationLayer.GENOMIC

    def test_annotation_qc_counts_sum_correctly(self):
        """total_variants == clean + warnings + failed."""
        from ga4gh.vrs._internal.models import (
            Allele,
            LiteralSequenceExpression,
            SequenceLocation,
            SequenceReference,
        )

        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )

        allele = Allele(
            location=SequenceLocation(
                sequenceReference=SequenceReference(refgetAccession="SQ." + "A" * 32),
                start=0,
                end=1,
            ),
            state=LiteralSequenceExpression(sequence="A"),
        )

        annotations = [
            _make_annotation(AnnotationLayer.GENOMIC, post_mapped=allele),  # clean
            _make_annotation(
                AnnotationLayer.GENOMIC, post_mapped=allele, error_message="warn"
            ),  # warning
            _make_annotation(AnnotationLayer.GENOMIC),  # failed
        ]

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": annotations},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        assert result.target_mappings
        tm = result.target_mappings[0]
        assert tm.total_variants == 3
        assert tm.variants_mapped_cleanly == 1
        assert tm.variants_with_mapping_warnings == 1
        assert tm.variants_failed == 1
        assert (
            (tm.variants_mapped_cleanly or 0)
            + (tm.variants_with_mapping_warnings or 0)
            + (tm.variants_failed or 0)
        ) == tm.total_variants

    def test_tool_parameters_for_sequence_based_target(self):
        """tool_parameters should contain BLAT aligner key for sequence-based targets."""
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        ann = _make_annotation(AnnotationLayer.GENOMIC)

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [ann]},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        tm = result.target_mappings[0]
        assert tm.tool_parameters is not None
        assert tm.tool_parameters.get("aligner") == "blat"
        assert "min_score" in tm.tool_parameters

    def test_tool_parameters_for_accession_based_target(self):
        """tool_parameters should reference cdot for accession-based targets."""
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_acc_target("NM_000001.1", "GENE1")},
        )
        ann = _make_annotation(AnnotationLayer.CDNA)

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.CDNA,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [ann]},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        tm = result.target_mappings[0]
        assert tm.tool_parameters is not None
        assert tm.tool_parameters.get("aligner") == "cdot_transcript_placement"
        assert "cdot_url" in tm.tool_parameters

    def test_tool_parameters_for_nc_accession_target(self):
        """NC_ (chromosome/contig) targets should use 'reference_accession_passthrough', not cdot."""
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_acc_target("NC_000001.11", "GENE1")},
        )
        ann = _make_annotation(AnnotationLayer.GENOMIC)

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [ann]},
                # NC_ targets produce align_result=None in fetch_alignment
                align_results={"NC_000001.11": None},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        tm = result.target_mappings[0]
        assert tm.tool_parameters is not None
        assert tm.tool_parameters.get("aligner") == "reference_accession_passthrough"
        # Must not bleed cdot fields into a non-cdot path
        assert "cdot_data_version" not in tm.tool_parameters
        assert "cdot_url" not in tm.tool_parameters
        # Accession is already in reference_accession; must not duplicate here
        assert "accession" not in tm.tool_parameters

    def test_mapped_scores_alignment_levels_subset_of_target_mappings(self):
        """Invariant: every alignment_level in mapped_scores must have a parent TargetMapping row.

        No mapped_score should be orphaned (i.e. have an alignment_level that
        does not correspond to any TargetMapping.alignment_level for this run).
        This mirrors the join the MaveDB API performs via target_gene_mapping_id.
        """
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        # Genomic success, protein success, and a completely-failed (NULL-layer) variant
        g_ann = _make_annotation(AnnotationLayer.GENOMIC)
        p_ann = _make_annotation(AnnotationLayer.PROTEIN)
        null_ann = _make_annotation(None)  # completely failed - no layer attribution

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [g_ann, p_ann, null_ann]},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        assert result.mapped_scores is not None
        assert result.target_mappings is not None

        tm_levels = {tm.alignment_level for tm in result.target_mappings}
        ms_levels = {ms.alignment_level for ms in result.mapped_scores}

        # Core invariant: every alignment_level in mapped_scores has a parent row
        assert ms_levels <= tm_levels, (
            f"Orphaned alignment levels in mapped_scores: {ms_levels - tm_levels}"
        )

    def test_null_layer_failures_attributed_to_preferred_layer(self):
        """NULL-layer (completely-failed) variants must be re-attributed to the preferred
        layer so the MaveDB API can join them via alignment_level.
        """
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        g_ann = _make_annotation(AnnotationLayer.GENOMIC)
        null_ann = _make_annotation(None)  # completely failed

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [g_ann, null_ann]},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=True,
                vrs_version=VrsVersion.V_2,
            )

        assert result.mapped_scores is not None
        assert result.target_mappings is not None

        # NULL-layer failure should appear in mapped_scores as the preferred layer
        for ms in result.mapped_scores:
            assert ms.alignment_level == AnnotationLayer.GENOMIC, (
                f"Expected alignment_level=GENOMIC for all mapped_scores; got {ms.alignment_level}"
            )

        # The preferred layer's TargetMapping must count the null failure
        tm = result.target_mappings[0]
        assert tm.alignment_level == AnnotationLayer.GENOMIC
        assert tm.total_variants == 2  # one genomic + one null failure

    def test_preferred_layer_only_total_variants_includes_null_failures(self):
        """With preferred_layer_only=True, total_variants in TargetMapping must
        account for null-layer failures so len(mapped_scores) == total_variants.
        """
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        # 3 genomic successes + 2 completely-failed variants
        g_anns = [_make_annotation(AnnotationLayer.GENOMIC) for _ in range(3)]
        null_anns = [_make_annotation(None) for _ in range(2)]

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate.get_chromosome_identifier",
                return_value="NC_000001.11",
            ),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.GENOMIC,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": g_anns + null_anns},
                align_results={"GENE1": _make_align()},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=True,
                vrs_version=VrsVersion.V_2,
            )

        assert result.mapped_scores is not None
        assert result.target_mappings is not None

        tm = result.target_mappings[0]
        assert len(result.mapped_scores) == 5
        assert tm.total_variants == 5
        assert tm.variants_failed == 5  # all have pre_mapped=None (failed)


# ---------------------------------------------------------------------------
# Locus flag stamping
# ---------------------------------------------------------------------------


def _make_allele_annotation(
    start: int, end: int, layer: AnnotationLayer
) -> ScoreAnnotation:
    """Build a ScoreAnnotation with a minimal VRS 2 Allele post_mapped."""
    from ga4gh.vrs._internal.models import (
        Allele,
        LiteralSequenceExpression,
        SequenceLocation,
        SequenceReference,
    )

    allele = Allele(
        location=SequenceLocation(
            sequenceReference=SequenceReference(refgetAccession="SQ." + "A" * 32),
            start=start,
            end=end,
        ),
        state=LiteralSequenceExpression(sequence="T"),
    )
    return ScoreAnnotation(
        mavedb_id=f"urn:mavedb:00000001-a-1#{start}",
        alignment_level=layer,
        post_mapped=allele,
        score=None,
    )


def _make_rle_annotation(start: int, end: int) -> ScoreAnnotation:
    """Build a ScoreAnnotation with a ReferenceLengthExpression (reference-identical)."""
    from ga4gh.vrs._internal.models import (
        Allele,
        ReferenceLengthExpression,
        SequenceLocation,
        SequenceReference,
    )

    allele = Allele(
        location=SequenceLocation(
            sequenceReference=SequenceReference(refgetAccession="SQ." + "A" * 32),
            start=start,
            end=end,
        ),
        state=ReferenceLengthExpression(length=end - start, sequence=None),
    )
    return ScoreAnnotation(
        mavedb_id=f"urn:mavedb:00000001-a-1#rle{start}",
        alignment_level=AnnotationLayer.GENOMIC,
        post_mapped=allele,
        score=None,
    )


def _make_align_result_with_qc(
    mismatches: list[int],
    gap_intervals: list[tuple[int, int]],
    mismatch_count: int | None = None,
) -> AlignmentResult:
    qc = AlignmentQc(
        alignment_length=1000,
        mismatch_count=mismatch_count
        if mismatch_count is not None
        else len(mismatches),
        gap_count=len(gap_intervals),
        mismatch_positions=mismatches,
        gap_intervals=gap_intervals,
    )
    return AlignmentResult(
        chrom="chr1",
        query_range=SequenceRange(start=0, end=100),
        query_subranges=[SequenceRange(start=0, end=100)],
        hit_range=SequenceRange(start=0, end=1000),
        hit_subranges=[SequenceRange(start=0, end=1000)],
        alignment_qc=qc,
    )


class TestStampAlignmentLocusFlags:
    def test_variant_at_mismatch(self):
        ann = _make_allele_annotation(10, 11, AnnotationLayer.GENOMIC)
        align = _make_align_result_with_qc(mismatches=[10], gap_intervals=[])
        _stamp_alignment_locus_flags([ann], align, AnnotationLayer.GENOMIC)
        assert ann.at_mismatched_locus is True
        assert ann.near_gap is False

    def test_variant_near_gap(self):
        ann = _make_allele_annotation(10, 11, AnnotationLayer.GENOMIC)
        # Gap at [20, 30); within 5-base window of position 10+5=15 < 20, so NOT near
        # Place gap closer: [13, 20)
        align = _make_align_result_with_qc(mismatches=[], gap_intervals=[(13, 20)])
        _stamp_alignment_locus_flags(
            [ann], align, AnnotationLayer.GENOMIC, near_gap_window=5
        )
        assert ann.at_mismatched_locus is False
        assert ann.near_gap is True  # pos 10, window end 15, gap start 13 < 15

    def test_variant_not_near_gap(self):
        ann = _make_allele_annotation(0, 1, AnnotationLayer.GENOMIC)
        # Gap far away
        align = _make_align_result_with_qc(mismatches=[], gap_intervals=[(1000, 1100)])
        _stamp_alignment_locus_flags([ann], align, AnnotationLayer.GENOMIC)
        assert ann.near_gap is False

    def test_rle_allele_flags_left_as_none(self):
        """RLE alleles (reference-identical) must leave both flags as None."""
        ann = _make_rle_annotation(0, 198295559)  # whole-chromosome span
        # Use a real mismatch/gap so the function would flag a normal allele
        align = _make_align_result_with_qc(
            mismatches=list(range(100)),
            gap_intervals=[(200, 300)],
        )
        _stamp_alignment_locus_flags([ann], align, AnnotationLayer.GENOMIC)
        assert ann.at_mismatched_locus is None
        assert ann.near_gap is None

    def test_idempotent_double_call(self):
        """Calling _stamp_alignment_locus_flags twice must not change the result."""
        ann = _make_allele_annotation(10, 11, AnnotationLayer.GENOMIC)
        align = _make_align_result_with_qc(mismatches=[10], gap_intervals=[(15, 25)])
        _stamp_alignment_locus_flags([ann], align, AnnotationLayer.GENOMIC)
        first_at = ann.at_mismatched_locus
        first_gap = ann.near_gap
        _stamp_alignment_locus_flags([ann], align, AnnotationLayer.GENOMIC)
        assert ann.at_mismatched_locus == first_at
        assert ann.near_gap == first_gap

    def test_no_align_result_leaves_flags_none(self):
        ann = _make_allele_annotation(10, 11, AnnotationLayer.GENOMIC)
        _stamp_alignment_locus_flags([ann], None, AnnotationLayer.GENOMIC)
        assert ann.at_mismatched_locus is None
        assert ann.near_gap is None


# ---------------------------------------------------------------------------
# _align_result_for_target
# ---------------------------------------------------------------------------


class TestAlignResultForTarget:
    def test_sequence_based_found_by_gene_name(self):
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        align = _make_align()
        result = _align_result_for_target("GENE1", metadata, {"GENE1": align})
        assert result is align

    def test_accession_based_found_by_accession(self):
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_acc_target("NM_000001.1", "GENE1")},
        )
        align = _make_align()
        # align_results keyed by accession (as produced by parse_cdot_mapping)
        result = _align_result_for_target("GENE1", metadata, {"NM_000001.1": align})
        assert result is align

    def test_returns_none_when_not_found(self):
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_seq_target("GENE1")},
        )
        result = _align_result_for_target("GENE1", metadata, {})
        assert result is None


# ---------------------------------------------------------------------------
# accession-based tool_parameters: cdot_data_version propagation
# ---------------------------------------------------------------------------


class TestCdotDataVersionPropagation:
    def test_cdot_data_version_propagates_when_present(self):
        """tool_parameters.cdot_data_version is emitted even when explicitly None."""
        metadata = ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1",
            target_genes={"GENE1": _make_acc_target("NM_000001.1", "GENE1")},
        )
        ann = _make_annotation(AnnotationLayer.CDNA)

        align_with_version = AlignmentResult(
            chrom="chr1",
            query_range=SequenceRange(start=0, end=100),
            query_subranges=[SequenceRange(start=0, end=100)],
            hit_range=SequenceRange(start=1000, end=1100),
            hit_subranges=[SequenceRange(start=1000, end=1100)],
            aligner_parameters={
                "aligner": "cdot_transcript_placement",
                "cdot_data_version": "0.2.26",
            },
        )

        with (
            patch("dcd_mapping.annotate.get_vrs_id_from_identifier", return_value=None),
            patch(
                "dcd_mapping.annotate._pick_preferred_genomic_or_protein_layer",
                return_value=AnnotationLayer.CDNA,
            ),
            patch(
                "dcd_mapping.annotate._get_computed_reference_sequence",
                return_value=None,
            ),
            patch(
                "dcd_mapping.annotate._get_mapped_reference_sequence", return_value=None
            ),
        ):
            result = build_scoreset_mapping(
                metadata=metadata,
                raw_metadata={},
                mappings={"GENE1": [ann]},
                align_results={"NM_000001.1": align_with_version},
                tx_output={"GENE1": _make_tx()},
                gene_info={"GENE1": None},
                preferred_layer_only=False,
                vrs_version=VrsVersion.V_2,
            )

        tm = result.target_mappings[0]
        assert tm.tool_parameters is not None
        assert "cdot_data_version" in tm.tool_parameters
        assert tm.tool_parameters["cdot_data_version"] == "0.2.26"
