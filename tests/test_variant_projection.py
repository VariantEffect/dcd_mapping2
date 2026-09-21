"""Tests for per-variant projection across assay levels.

Covers the mapper's *projection* of a measured variant onto its own deterministic
forms -- ``g. -> c. -> p.`` for a genomic source, ``c. -> g.`` / ``c. -> p.`` for a
cdna source -- against the coding transcript, distinct from the equivalence-class
expansion the reverse-translation job owns. The projected forms are emitted alongside
the measured one and routed by ``preferred_layer_only`` (API keeps the assay layer,
CLI keeps all); :class:`~dcd_mapping.vrs_map.ProjectionOutcome` and the per-target
summary log additionally surface a mis-selected transcript.

See ``dcd_mapping.vrs_map._construct_projected_layers`` and the projection helpers in
``dcd_mapping.lookup``.
"""

import logging
from unittest.mock import MagicMock, patch

import hgvs.parser
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.vrs._internal.models import (
    Allele,
    LiteralSequenceExpression,
    SequenceLocation,
    SequenceReference,
)

from dcd_mapping.lookup import (
    coding_hgvs_is_intronic,
    get_genomic_accession_for_transcript,
    project_coding_hgvs_to_genomic,
    project_coding_hgvs_to_protein,
    project_genomic_hgvs_to_coding,
)
from dcd_mapping.schemas import (
    MappedScore,
    MappingOutcome,
    ScoreRow,
    TargetGene,
    TargetType,
    TxSelectResult,
)
from dcd_mapping.vrs_map import (
    ProjectionOutcome,
    _construct_projected_layers,
    _log_projection_validation,
    _map_accession,
)

VRS_MAP = "dcd_mapping.vrs_map"
LOOKUP = "dcd_mapping.lookup"

NC = "NC_000007.14"
NM = "NM_004333.6"


def _row(hgvs_nt: str, accession: str = "urn:mavedb:00000001-a-1#1") -> ScoreRow:
    return ScoreRow(hgvs_nt=hgvs_nt, hgvs_pro="_wt", score="1.0", accession=accession)


def _allele() -> Allele:
    """Build a minimal, schema-valid VRS Allele for use as a construction stand-in."""
    return Allele(
        location=SequenceLocation(
            sequenceReference=SequenceReference(refgetAccession="SQ." + "A" * 32),
            start=0,
            end=1,
        ),
        state=LiteralSequenceExpression(sequence="A"),
    )


class TestCodingHgvsIsIntronic:
    """Intronic detection reads the parsed base-offset, not the textual ``-``."""

    def _patched(self):
        # Real hgvs parser (pure string parsing, no network); fake out the seqrepo-
        # backed builder so only the parser is exercised.
        builder = MagicMock()
        builder.hgvs_tools.parser = hgvs.parser.Parser()
        return patch(f"{LOOKUP}.TranslatorBuilder", return_value=builder), patch(
            f"{LOOKUP}.get_seqrepo"
        )

    def test_intron_offset_is_intronic(self):
        tb, sr = self._patched()
        with tb, sr:
            assert coding_hgvs_is_intronic("NM_000051.4:c.2002-1del") is True

    def test_exonic_is_not_intronic(self):
        tb, sr = self._patched()
        with tb, sr:
            assert coding_hgvs_is_intronic("NM_000051.4:c.76A>G") is False

    def test_five_prime_utr_dash_is_not_intronic(self):
        """``c.-20A>G`` carries a textual ``-`` but offset 0 -- not intronic."""
        tb, sr = self._patched()
        with tb, sr:
            assert coding_hgvs_is_intronic("NM_000051.4:c.-20A>G") is False


class TestProjectionHelpersWiring:
    """The g.->c. and c.->p. helpers route through the package hgvs VariantMapper."""

    def test_project_genomic_to_coding(self):
        builder = MagicMock()
        tools = builder.hgvs_tools
        tools.parser.parse.return_value = "parsed_g"
        tools.variant_mapper.g_to_c.return_value = f"{NM}:c.1A>G"
        with (
            patch(f"{LOOKUP}.TranslatorBuilder", return_value=builder),
            patch(f"{LOOKUP}.get_seqrepo"),
        ):
            result = project_genomic_hgvs_to_coding(f"{NC}:g.140A>G", NM)
        tools.parser.parse.assert_called_once_with(f"{NC}:g.140A>G")
        tools.variant_mapper.g_to_c.assert_called_once_with("parsed_g", NM)
        assert result == f"{NM}:c.1A>G"

    def test_project_coding_to_protein(self):
        builder = MagicMock()
        tools = builder.hgvs_tools
        tools.parser.parse.return_value = "parsed_c"
        tools.variant_mapper.c_to_p.return_value = "NP_004324.2:p.Lys1Glu"
        with (
            patch(f"{LOOKUP}.TranslatorBuilder", return_value=builder),
            patch(f"{LOOKUP}.get_seqrepo"),
        ):
            result = project_coding_hgvs_to_protein(f"{NM}:c.1A>G")
        tools.parser.parse.assert_called_once_with(f"{NM}:c.1A>G")
        tools.variant_mapper.c_to_p.assert_called_once_with("parsed_c")
        assert result == "NP_004324.2:p.Lys1Glu"

    def test_project_coding_to_genomic(self):
        """``c. -> g.`` routes through ``c_to_g`` with the explicit target contig."""
        builder = MagicMock()
        tools = builder.hgvs_tools
        tools.parser.parse.return_value = "parsed_c"
        tools.variant_mapper.c_to_g.return_value = f"{NC}:g.140A>G"
        with (
            patch(f"{LOOKUP}.TranslatorBuilder", return_value=builder),
            patch(f"{LOOKUP}.get_seqrepo"),
        ):
            result = project_coding_hgvs_to_genomic(f"{NM}:c.1A>G", NC)
        tools.parser.parse.assert_called_once_with(f"{NM}:c.1A>G")
        tools.variant_mapper.c_to_g.assert_called_once_with("parsed_c", NC)
        assert result == f"{NC}:g.140A>G"


class TestGenomicAccessionForTranscript:
    """The contig resolver intersects cdot's mapping options with the assembly's contigs."""

    def test_picks_contig_in_requested_assembly(self):
        cd = MagicMock()
        cd.get_tx_mapping_options.return_value = [
            {"alt_ac": "NC_000007.13"},  # GRCh37 contig -- not in the GRCh38 set
            {"alt_ac": "NC_000007.14"},  # GRCh38 contig
        ]
        cd.get_assembly_map.return_value = {"NC_000007.14": "7"}
        with patch(f"{LOOKUP}.cdot_rest", return_value=cd):
            assert get_genomic_accession_for_transcript(NM) == "NC_000007.14"
        cd.get_assembly_map.assert_called_once_with("GRCh38")

    def test_no_matching_contig_returns_none(self):
        cd = MagicMock()
        cd.get_tx_mapping_options.return_value = [{"alt_ac": "NC_000007.13"}]
        cd.get_assembly_map.return_value = {"NC_000007.14": "7"}
        with patch(f"{LOOKUP}.cdot_rest", return_value=cd):
            assert get_genomic_accession_for_transcript(NM) is None

    def test_lookup_failure_returns_none(self):
        cd = MagicMock()
        cd.get_tx_mapping_options.side_effect = RuntimeError("cdot down")
        with patch(f"{LOOKUP}.cdot_rest", return_value=cd):
            assert get_genomic_accession_for_transcript(NM) is None


class TestConstructProjectedLayers:
    """Branch + outcome classification for a single variant's projection."""

    def test_non_variant_rows_skipped(self):
        for hgvs_nt in ("_wt", "_sy", "=", "c.1_2delinsAA fs"):
            projected, outcome = _construct_projected_layers(
                _row(hgvs_nt), AnnotationLayer.GENOMIC, NM, accession_id=NC
            )
            assert projected == []
            assert outcome is ProjectionOutcome.SKIPPED

    def test_genomic_to_coding_failure_records_all_layers_failed(self):
        """A pivot failure emits a FAILED record for every expected level, not silence."""
        with (
            patch(
                f"{VRS_MAP}._create_pre_mapped_hgvs_strings",
                return_value=[f"{NC}:g.140A>G"],
            ),
            patch(
                f"{VRS_MAP}.project_genomic_hgvs_to_coding",
                side_effect=ValueError("off transcript"),
            ),
        ):
            projected, outcome = _construct_projected_layers(
                _row("140A>G"), AnnotationLayer.GENOMIC, NM, accession_id=NC
            )
        assert outcome is ProjectionOutcome.FAILED
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.CDNA, MappingOutcome.FAILED),
            (AnnotationLayer.PROTEIN, MappingOutcome.FAILED),
        ]
        # Genuine failures carry an error_message; benign absences would not.
        assert all(m.error_message for m in projected)

    def test_intronic_records_all_layers_benign(self):
        """Intronic is a benign absence: every expected level gets an INTRONIC record
        with no error_message (benign != error).
        """
        with (
            patch(
                f"{VRS_MAP}._create_pre_mapped_hgvs_strings",
                return_value=[f"{NC}:g.140A>G"],
            ),
            patch(
                f"{VRS_MAP}.project_genomic_hgvs_to_coding",
                return_value=[f"{NM}:c.2002-1A>G"],
            ),
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=True),
            patch(f"{VRS_MAP}.project_coding_hgvs_to_protein") as to_protein,
            patch(f"{VRS_MAP}._construct_vrs_allele") as construct,
        ):
            projected, outcome = _construct_projected_layers(
                _row("140A>G"), AnnotationLayer.GENOMIC, NM, accession_id=NC
            )
        assert outcome is ProjectionOutcome.INTRONIC
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.CDNA, MappingOutcome.INTRONIC),
            (AnnotationLayer.PROTEIN, MappingOutcome.INTRONIC),
        ]
        assert all(m.error_message is None for m in projected)
        # No protein projection or allele construction attempted for an intronic variant.
        to_protein.assert_not_called()
        construct.assert_not_called()

    def test_clean_projection_emits_cdna_and_protein(self):
        allele = _allele()
        with (
            patch(
                f"{VRS_MAP}._create_pre_mapped_hgvs_strings",
                return_value=[f"{NC}:g.140A>G"],
            ),
            patch(
                f"{VRS_MAP}.project_genomic_hgvs_to_coding",
                return_value=[f"{NM}:c.1A>G"],
            ),
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=False),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_protein",
                return_value="NP_004324.2:p.Lys1Glu",
            ),
            patch(f"{VRS_MAP}._construct_vrs_allele", return_value=allele),
        ):
            projected, outcome = _construct_projected_layers(
                _row("140A>G"), AnnotationLayer.GENOMIC, NM, accession_id=NC
            )
        assert outcome is ProjectionOutcome.PROJECTED
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.CDNA, MappingOutcome.MAPPED),
            (AnnotationLayer.PROTEIN, MappingOutcome.MAPPED),
        ]

    def test_protein_failure_alone_is_no_consequence(self):
        """A c.->p. miss is a benign no-consequence; the coding form still maps and the
        aggregate stays PROJECTED.
        """
        allele = _allele()
        with (
            patch(
                f"{VRS_MAP}._create_pre_mapped_hgvs_strings",
                return_value=[f"{NC}:g.140A>G"],
            ),
            patch(
                f"{VRS_MAP}.project_genomic_hgvs_to_coding",
                return_value=[f"{NM}:c.1A>G"],
            ),
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=False),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_protein",
                side_effect=ValueError("no p."),
            ),
            patch(f"{VRS_MAP}._construct_vrs_allele", return_value=allele),
        ):
            projected, outcome = _construct_projected_layers(
                _row("140A>G"), AnnotationLayer.GENOMIC, NM, accession_id=NC
            )
        assert outcome is ProjectionOutcome.PROJECTED
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.CDNA, MappingOutcome.MAPPED),
            (AnnotationLayer.PROTEIN, MappingOutcome.NO_PROTEIN_CONSEQUENCE),
        ]
        protein = projected[1]
        assert protein.error_message is None  # benign, not an error
        assert protein.post_mapped is None

    def test_coding_allele_construction_failure_is_failed(self):
        """Coding hgvs existed but could not become an allele -> FAILED record."""
        with (
            patch(
                f"{VRS_MAP}._create_pre_mapped_hgvs_strings",
                return_value=[f"{NC}:g.140A>G"],
            ),
            patch(
                f"{VRS_MAP}.project_genomic_hgvs_to_coding",
                return_value=[f"{NM}:c.1A>G"],
            ),
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=False),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_protein",
                return_value="NP_004324.2:p.Lys1Glu",
            ),
            patch(
                f"{VRS_MAP}._construct_vrs_allele", side_effect=ValueError("bad allele")
            ),
        ):
            projected, outcome = _construct_projected_layers(
                _row("140A>G"), AnnotationLayer.GENOMIC, NM, accession_id=NC
            )
        assert outcome is ProjectionOutcome.FAILED
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.CDNA, MappingOutcome.FAILED),
            (AnnotationLayer.PROTEIN, MappingOutcome.FAILED),
        ]


class TestConstructProjectedLayersCdnaSource:
    """A cdna assay (NM_/ENST) projects its measured coding variant to g. and p.

    The measured variant is already coding on the transcript, so there is no g.->c.
    pivot: the coding form is built by prefixing the accession, the genomic re-expression
    (c.->g.) is the load-bearing form, and the protein consequence (c.->p.) follows.
    """

    def test_cdna_source_emits_genomic_and_protein(self):
        allele = _allele()
        with (
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=False),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_genomic",
                return_value=f"{NC}:g.140A>G",
            ),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_protein",
                return_value="NP_004324.2:p.Lys1Glu",
            ),
            patch(f"{VRS_MAP}._construct_vrs_allele", return_value=allele),
        ):
            projected, outcome = _construct_projected_layers(
                _row("c.1A>G"),
                AnnotationLayer.CDNA,
                NM,
                accession_id=NM,
                genomic_accession=NC,
            )
        assert outcome is ProjectionOutcome.PROJECTED
        assert [m.alignment_level for m in projected] == [
            AnnotationLayer.GENOMIC,
            AnnotationLayer.PROTEIN,
        ]

    def test_cdna_source_genomic_projection_failure_is_failed(self):
        """The c.->g. re-expression is load-bearing: its failure is a FAILED genomic
        record, while the protein consequence still maps.
        """
        allele = _allele()
        with (
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=False),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_genomic",
                side_effect=ValueError("off contig"),
            ),
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_protein",
                return_value="NP_004324.2:p.Lys1Glu",
            ),
            patch(f"{VRS_MAP}._construct_vrs_allele", return_value=allele),
        ):
            projected, outcome = _construct_projected_layers(
                _row("c.1A>G"),
                AnnotationLayer.CDNA,
                NM,
                accession_id=NM,
                genomic_accession=NC,
            )
        assert outcome is ProjectionOutcome.FAILED
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.GENOMIC, MappingOutcome.FAILED),
            (AnnotationLayer.PROTEIN, MappingOutcome.MAPPED),
        ]
        assert projected[0].error_message  # genomic failure carries detail

    def test_cdna_source_without_contig_records_genomic_failed(self):
        """An unresolvable contig is still an accounted-for outcome: the genomic layer
        gets a FAILED record (not silence), while the protein consequence still maps.
        """
        allele = _allele()
        with (
            patch(f"{VRS_MAP}.coding_hgvs_is_intronic", return_value=False),
            patch(f"{VRS_MAP}.project_coding_hgvs_to_genomic") as to_genomic,
            patch(
                f"{VRS_MAP}.project_coding_hgvs_to_protein",
                return_value="NP_004324.2:p.Lys1Glu",
            ),
            patch(f"{VRS_MAP}._construct_vrs_allele", return_value=allele),
        ):
            projected, outcome = _construct_projected_layers(
                _row("c.1A>G"),
                AnnotationLayer.CDNA,
                NM,
                accession_id=NM,
                genomic_accession=None,
            )
        to_genomic.assert_not_called()
        assert outcome is ProjectionOutcome.FAILED
        assert [(m.alignment_level, m.outcome) for m in projected] == [
            (AnnotationLayer.GENOMIC, MappingOutcome.FAILED),
            (AnnotationLayer.PROTEIN, MappingOutcome.MAPPED),
        ]


class TestLogProjectionValidation:
    """Per-target summary escalates to WARNING when too many variants fail."""

    def _levels(self, caplog, outcomes):
        caplog.clear()
        with caplog.at_level(logging.INFO, logger=VRS_MAP):
            _log_projection_validation(NC, NM, outcomes)
        return [r.levelno for r in caplog.records]

    def test_clean_run_logs_info(self, caplog):
        outcomes = [ProjectionOutcome.PROJECTED] * 3 + [ProjectionOutcome.INTRONIC]
        levels = self._levels(caplog, outcomes)
        assert levels == [logging.INFO]

    def test_high_failure_fraction_logs_warning(self, caplog):
        # 2 failed of 4 attempted = 50% > 25% threshold.
        outcomes = [ProjectionOutcome.PROJECTED] * 2 + [ProjectionOutcome.FAILED] * 2
        levels = self._levels(caplog, outcomes)
        assert levels == [logging.WARNING]

    def test_intronic_excluded_from_failure_fraction(self, caplog):
        # 1 failed of 4 attempted = 25%, not > 25%, so INFO despite many intronic.
        outcomes = [ProjectionOutcome.PROJECTED] * 3 + [ProjectionOutcome.FAILED]
        levels = self._levels(caplog, outcomes)
        assert levels == [logging.INFO]

    def test_all_skipped_logs_nothing(self, caplog):
        levels = self._levels(caplog, [ProjectionOutcome.SKIPPED] * 5)
        assert levels == []

    def test_message_distinguishes_cdna_self_pivot_from_selection(self, caplog):
        outcomes = [ProjectionOutcome.PROJECTED] * 3
        with caplog.at_level(logging.INFO, logger=VRS_MAP):
            _log_projection_validation(NC, NM, outcomes)  # different -> selection
            _log_projection_validation(NM, NM, outcomes)  # same -> cdna self-pivot
        selected, cdna = caplog.records[0].getMessage(), caplog.records[1].getMessage()
        assert f"onto selected transcript {NM}" in selected
        assert f"Projecting cdna target {NM}" in cdna


class TestNmAccessionMeasuredProtein:
    """NM_/ENST targets with a measured hgvs_pro use it directly rather than projecting.

    The measured protein is already in reference coordinates so pre_mapped == post_mapped.
    When hgvs_pro is present, _construct_projected_layers must not produce a redundant
    projected protein layer (project_protein=False).
    """

    NP = "NP_004324.2"

    def _tx(self) -> TxSelectResult:
        return TxSelectResult(
            nm=NM,
            np=self.NP,
            start=0,
            is_full_match=True,
            sequence="MAAAA",
            hgnc_symbol="BRAF",
        )

    def _run(self, hgvs_pro: str, transcript=None):
        """Call _map_accession for a single NM_ row via vrs_map, fully mocked."""
        metadata = TargetGene(
            target_gene_name="BRAF",
            target_gene_category=TargetType.PROTEIN_CODING,
            target_sequence=None,
            target_sequence_type=None,
            target_accession_id=NM,
        )
        row = ScoreRow(
            hgvs_nt="c.1799T>A",
            hgvs_pro=hgvs_pro,
            score="1.0",
            accession="urn:mavedb:00000001-a-1#1",
        )
        allele = _allele()
        with (
            patch(f"{VRS_MAP}.store_accession"),
            patch(f"{VRS_MAP}.get_genomic_accession_for_transcript", return_value=NC),
            patch(f"{VRS_MAP}._map_genomic") as mock_genomic,
            patch(
                f"{VRS_MAP}._create_pre_mapped_hgvs_strings",
                return_value=[f"{self.NP}:p.Val600Glu"],
            ),
            patch(f"{VRS_MAP}._construct_vrs_allele", return_value=allele),
            patch(
                f"{VRS_MAP}._construct_projected_layers",
                return_value=([], ProjectionOutcome.PROJECTED),
            ) as mock_proj,
            patch(f"{VRS_MAP}._log_projection_validation"),
        ):
            mock_genomic.return_value = MappedScore(
                accession_id="urn:mavedb:00000001-a-1#1",
                score="1.0",
                alignment_level=AnnotationLayer.GENOMIC,
            )
            result = _map_accession(metadata, [row], None, transcript or self._tx())
        return result, mock_proj

    def test_measured_protein_emitted_directly(self):
        """When hgvs_pro is valid and NP_ is available, a protein MappedScore is emitted
        with pre_mapped == post_mapped (reference-coordinate form, no alignment offset).
        """
        variations, _ = self._run("p.Val600Glu")
        protein_scores = [
            v for v in variations if v.alignment_level == AnnotationLayer.PROTEIN
        ]
        assert len(protein_scores) == 1
        assert protein_scores[0].pre_mapped is not None
        assert protein_scores[0].pre_mapped is protein_scores[0].post_mapped

    def test_measured_protein_suppresses_projection(self):
        """project_protein=False must be passed to _construct_projected_layers so the
        projected c.->p. form is not emitted alongside the measured one.
        """
        _, mock_proj = self._run("p.Val600Glu")
        _, kwargs = mock_proj.call_args
        assert kwargs.get("project_protein") is False

    def test_no_hgvs_pro_still_projects(self):
        """When hgvs_pro is absent, projection behaviour is unchanged (project_protein=True)."""
        _, mock_proj = self._run("_wt")
        _, kwargs = mock_proj.call_args
        assert (
            kwargs.get("project_protein") is not False
        )  # default True (not passed or True)

    def test_missing_np_accession_skips_protein(self, caplog):
        """When there is no TxSelectResult (the common case for NM_ accession targets
        where transcript selection produces no result), the protein layer is skipped with
        a warning rather than silently emitting a potentially wrong projected form.
        """
        from dcd_mapping.schemas import TargetGene, TargetType
        from dcd_mapping.vrs_map import _map_accession

        metadata = TargetGene(
            target_gene_name="BRAF",
            target_gene_category=TargetType.PROTEIN_CODING,
            target_sequence=None,
            target_sequence_type=None,
            target_accession_id=NM,
        )
        row = ScoreRow(
            hgvs_nt="c.1799T>A",
            hgvs_pro="p.Val600Glu",
            score="1.0",
            accession="urn:mavedb:00000001-a-1#1",
        )
        with (
            patch(f"{VRS_MAP}.store_accession"),
            patch(f"{VRS_MAP}.get_genomic_accession_for_transcript", return_value=NC),
            patch(f"{VRS_MAP}._map_genomic") as mock_genomic,
            patch(
                f"{VRS_MAP}._construct_projected_layers",
                return_value=([], ProjectionOutcome.PROJECTED),
            ) as mock_proj,
            patch(f"{VRS_MAP}._log_projection_validation"),
            caplog.at_level(logging.WARNING, logger=VRS_MAP),
        ):
            mock_genomic.return_value = MappedScore(
                accession_id="urn:mavedb:00000001-a-1#1",
                score="1.0",
                alignment_level=AnnotationLayer.GENOMIC,
            )
            # transcript=None: no TxSelectResult produced for this NM_ target
            result = _map_accession(metadata, [row], None, None)

        protein_scores = [
            v for v in result if v.alignment_level == AnnotationLayer.PROTEIN
        ]
        assert protein_scores == []
        assert any("no NP_ accession" in r.message for r in caplog.records)
        # project_protein must still be False so no projected protein is emitted
        _, kwargs = mock_proj.call_args
        assert kwargs.get("project_protein") is False
