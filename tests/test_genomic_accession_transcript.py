"""Tests for genomic-accession (NC_) coding-target transcript selection.

Covers the path added so that genomic-accession coding targets resolve a coding
(MANE) transcript -- gene inferred from the variants' genomic loci (preferred) or
the declared target metadata (fallback) -- which the mapper surfaces as cdna target
metadata for reverse translation. See ``dcd_mapping.transcripts``.
"""

from unittest.mock import patch

import pytest
from cool_seq_tool.schemas import TranscriptPriority

from dcd_mapping.exceptions import NoCodingTranscriptError
from dcd_mapping.lookup import infer_hgnc_symbol_from_genomic_loci
from dcd_mapping.schemas import (
    ManeDescription,
    ScoreRow,
    ScoresetMetadata,
    TargetGene,
    TargetType,
    TxSelectResult,
)
from dcd_mapping.transcripts import (
    _genomic_positions_from_records,
    _select_genomic_accession_reference,
    _select_refseq_cdna_counterpart,
    _select_refseq_protein_counterpart,
    select_transcripts,
)

MODULE = "dcd_mapping.transcripts"
LOOKUP = "dcd_mapping.lookup"


def _mane(
    nm: str = "NM_004333.6",
    np: str = "NP_004324.2",
    symbol: str = "BRAF",
    priority: TranscriptPriority = TranscriptPriority.MANE_SELECT,
) -> ManeDescription:
    return ManeDescription(
        refseq_nuc=nm,
        refseq_prot=np,
        transcript_priority=priority,
        ncbi_gene_id="673",
        ensembl_gene_id="ENSG00000157764",
        hgnc_gene_id="HGNC:1097",
        symbol=symbol,
        name="B-Raf proto-oncogene",
        ensembl_nuc="ENST00000646891.2",
        ensembl_prot="ENSP00000493543.1",
        grch38_chr="7",
        chr_start=140730665,
        chr_end=140924764,
        chr_strand="-",
    )


def _nc_coding_target(
    accession: str = "NC_000007.14", name: str = "BRAF"
) -> TargetGene:
    return TargetGene(
        target_gene_name=name,
        target_gene_category=TargetType.PROTEIN_CODING,
        target_accession_id=accession,
        target_accession_assembly="GRCh38",
    )


def _row(hgvs_nt: str) -> ScoreRow:
    return ScoreRow(hgvs_nt=hgvs_nt, hgvs_pro="_wt", score="1.0", accession="urn#1")


def _accession_target(accession: str) -> TargetGene:
    return TargetGene(
        target_gene_name="BRAF",
        target_gene_category=TargetType.PROTEIN_CODING,
        target_accession_id=accession,
    )


class TestSelectRefseqProteinCounterpart:
    def test_maps_ensembl_protein_to_mane_refseq(self):
        with (
            patch(
                f"{MODULE}.get_gene_symbol_from_ensembl_protein", return_value="BRAF"
            ) as symbol,
            patch(
                f"{MODULE}.get_mane_transcripts_for_gene", return_value=[_mane()]
            ) as mane,
        ):
            result = _select_refseq_protein_counterpart("ENSP00000493543.1")

        assert isinstance(result, TxSelectResult)
        assert result.nm == "NM_004333.6"
        assert result.np == "NP_004324.2"
        assert result.hgnc_symbol == "BRAF"
        assert result.transcript_mode == TranscriptPriority.MANE_SELECT
        symbol.assert_called_once_with("ENSP00000493543.1")
        mane.assert_called_once_with("BRAF")

    def test_returns_none_when_gene_symbol_unresolved(self):
        with patch(f"{MODULE}.get_gene_symbol_from_ensembl_protein", return_value=None):
            assert _select_refseq_protein_counterpart("ENSP00000000000.1") is None

    def test_returns_none_when_no_mane_transcript(self):
        with (
            patch(
                f"{MODULE}.get_gene_symbol_from_ensembl_protein", return_value="BRAF"
            ),
            patch(f"{MODULE}.get_mane_transcripts_for_gene", return_value=[]),
        ):
            assert _select_refseq_protein_counterpart("ENSP00000493543.1") is None


class TestSelectRefseqCdnaCounterpart:
    def test_maps_ensembl_transcript_to_mane_refseq(self):
        with (
            patch(
                f"{MODULE}.get_gene_symbol_from_ensembl_transcript", return_value="BRAF"
            ) as symbol,
            patch(
                f"{MODULE}.get_mane_transcripts_for_gene", return_value=[_mane()]
            ) as mane,
        ):
            result = _select_refseq_cdna_counterpart("ENST00000646891.2")

        assert isinstance(result, TxSelectResult)
        assert result.nm == "NM_004333.6"
        assert result.np == "NP_004324.2"
        assert result.hgnc_symbol == "BRAF"
        assert result.transcript_mode == TranscriptPriority.MANE_SELECT
        symbol.assert_called_once_with("ENST00000646891.2")
        mane.assert_called_once_with("BRAF")

    def test_returns_none_when_gene_symbol_unresolved(self):
        with patch(
            f"{MODULE}.get_gene_symbol_from_ensembl_transcript", return_value=None
        ):
            assert _select_refseq_cdna_counterpart("ENST00000000000.1") is None

    def test_returns_none_when_no_mane_transcript(self):
        with (
            patch(
                f"{MODULE}.get_gene_symbol_from_ensembl_transcript", return_value="BRAF"
            ),
            patch(f"{MODULE}.get_mane_transcripts_for_gene", return_value=[]),
        ):
            assert _select_refseq_cdna_counterpart("ENST00000646891.2") is None


class TestSelectGenomicAccessionReference:
    def test_inferred_gene_preferred(self):
        """The gene inferred from genomic loci is preferred and drives MANE selection."""
        target = _nc_coding_target()
        with (
            patch(
                f"{MODULE}.infer_hgnc_symbol_from_genomic_loci", return_value="BRAF"
            ) as infer,
            patch(f"{MODULE}.get_gene_symbol") as declared,
            patch(
                f"{MODULE}.get_mane_transcripts_for_gene", return_value=[_mane()]
            ) as mane,
        ):
            result = _select_genomic_accession_reference(target, [_row("g.123A>G")])

        assert isinstance(result, TxSelectResult)
        assert result.nm == "NM_004333.6"
        assert result.np == "NP_004324.2"
        assert result.hgnc_symbol == "BRAF"
        assert result.transcript_mode == TranscriptPriority.MANE_SELECT
        infer.assert_called_once()
        mane.assert_called_once_with("BRAF")
        # Declared-metadata normalization is not consulted when inference succeeds.
        declared.assert_not_called()

    def test_falls_back_to_declared_gene(self):
        """When locus inference yields nothing, fall back to target-metadata gene."""
        target = _nc_coding_target(name="Wildtype BRAF")
        with (
            patch(f"{MODULE}.infer_hgnc_symbol_from_genomic_loci", return_value=None),
            patch(f"{MODULE}.get_gene_symbol", return_value="BRAF") as declared,
            patch(f"{MODULE}.get_mane_transcripts_for_gene", return_value=[_mane()]),
        ):
            result = _select_genomic_accession_reference(target, [_row("g.123A>G")])

        assert result.nm == "NM_004333.6"
        declared.assert_called_once()

    def test_no_gene_raises_no_coding_transcript(self):
        """No resolvable gene -> typed NoCodingTranscriptError (recoverable skip)."""
        target = _nc_coding_target(name="mystery element")
        with (
            patch(f"{MODULE}.infer_hgnc_symbol_from_genomic_loci", return_value=None),
            patch(f"{MODULE}.get_gene_symbol", return_value=None),
            pytest.raises(NoCodingTranscriptError),
        ):
            _select_genomic_accession_reference(target, [_row("g.123A>G")])

    def test_no_mane_raises_no_coding_transcript(self):
        """Gene resolves but has no MANE transcript -> NoCodingTranscriptError."""
        target = _nc_coding_target()
        with (
            patch(f"{MODULE}.infer_hgnc_symbol_from_genomic_loci", return_value="BRAF"),
            patch(f"{MODULE}.get_mane_transcripts_for_gene", return_value=[]),
            pytest.raises(NoCodingTranscriptError),
        ):
            _select_genomic_accession_reference(target, [_row("g.123A>G")])

    def test_mane_plus_clinical_secondary(self):
        """MANE Plus Clinical is selected when no MANE Select is present."""
        target = _nc_coding_target()
        plus = _mane(priority=TranscriptPriority.MANE_PLUS_CLINICAL)
        with (
            patch(f"{MODULE}.infer_hgnc_symbol_from_genomic_loci", return_value="BRAF"),
            patch(f"{MODULE}.get_mane_transcripts_for_gene", return_value=[plus]),
        ):
            result = _select_genomic_accession_reference(target, [_row("g.123A>G")])
        assert result.transcript_mode == TranscriptPriority.MANE_PLUS_CLINICAL


class TestInferHgncSymbolFromGenomicLoci:
    def _feature(
        self,
        name: str,
        start: int = 0,
        end: int = 10**9,
        biotype: str = "protein_coding",
    ) -> dict:
        return {
            "external_name": name,
            "feature_type": "gene",
            "start": start,
            "end": end,
            "biotype": biotype,
        }

    def test_returns_gene_containing_loci_with_one_request(self):
        """A single Ensembl request over the bounding span resolves the gene -- not one
        request per position.
        """
        with (
            patch(
                f"{LOOKUP}.get_overlapping_features_for_region",
                return_value=[self._feature("BARD1", 214725646, 214809707)],
            ) as overlap,
            patch(f"{LOOKUP}._get_hgnc_symbol", return_value="BARD1"),
        ):
            result = infer_hgnc_symbol_from_genomic_loci(
                "NC_000002.12", [214728639, 214728788, 214728789, 214728789]
            )
        assert result == "BARD1"
        # One request for the whole target, regardless of variant/duplicate count.
        overlap.assert_called_once()

    def test_picks_gene_containing_most_loci(self):
        """When the bounding span overlaps several genes, the one containing the most
        query loci wins -- resolved locally from coordinates, no extra requests.
        """
        with (
            patch(
                f"{LOOKUP}.get_overlapping_features_for_region",
                return_value=[
                    self._feature("MAINGENE", 100, 200),
                    self._feature(
                        "EDGEGENE", 199, 1000
                    ),  # contains only the last locus
                ],
            ),
            patch(f"{LOOKUP}._get_hgnc_symbol", side_effect=lambda s: s),
        ):
            result = infer_hgnc_symbol_from_genomic_loci(
                "NC_000002.12", [110, 120, 130, 500]
            )
        assert result == "MAINGENE"

    def test_ignores_non_coding_overlapping_genes(self):
        """A non-coding gene overlapping the same loci (which has no MANE transcript) is
        excluded, so it neither wins nor forces a tie against the coding gene.
        """
        with (
            patch(
                f"{LOOKUP}.get_overlapping_features_for_region",
                return_value=[
                    self._feature("CODINGGENE", 0, 1000, biotype="protein_coding"),
                    self._feature("AS-LNCRNA", 0, 1000, biotype="lncRNA"),
                ],
            ),
            patch(f"{LOOKUP}._get_hgnc_symbol", side_effect=lambda s: s),
        ):
            result = infer_hgnc_symbol_from_genomic_loci("NC_000002.12", [10, 20, 30])
        assert result == "CODINGGENE"

    def test_no_protein_coding_overlap_returns_none(self):
        """Only non-coding genes overlap -> no coding gene to infer."""
        with patch(
            f"{LOOKUP}.get_overlapping_features_for_region",
            return_value=[self._feature("AS-LNCRNA", 0, 1000, biotype="lncRNA")],
        ):
            assert infer_hgnc_symbol_from_genomic_loci("NC_000002.12", [10]) is None

    def test_falls_back_to_ensembl_symbol_when_normalizer_misses(self):
        """A locus-confirmed gene is not dropped just because the gene normalizer
        returns nothing -- the Ensembl external_name is used directly.
        """
        with (
            patch(
                f"{LOOKUP}.get_overlapping_features_for_region",
                return_value=[self._feature("BARD1", 214725646, 214809707)],
            ),
            patch(f"{LOOKUP}._get_hgnc_symbol", return_value=None),
        ):
            assert (
                infer_hgnc_symbol_from_genomic_loci("NC_000002.12", [214728639])
                == "BARD1"
            )

    def test_no_overlap_returns_none(self):
        with patch(f"{LOOKUP}.get_overlapping_features_for_region", return_value=[]):
            assert infer_hgnc_symbol_from_genomic_loci("NC_000002.12", [1]) is None

    def test_tie_returns_none(self):
        # Two genes each contain the same loci equally -> no single dominant gene.
        with patch(
            f"{LOOKUP}.get_overlapping_features_for_region",
            return_value=[
                self._feature("GENEA", 0, 100),
                self._feature("GENEB", 0, 100),
            ],
        ):
            assert infer_hgnc_symbol_from_genomic_loci("NC_000002.12", [10, 20]) is None

    def test_empty_positions_returns_none(self):
        assert infer_hgnc_symbol_from_genomic_loci("NC_000002.12", []) is None

    def test_span_too_wide_skips_query(self):
        """A bounding span too wide to be one gene declines inference without firing a
        doomed oversized Ensembl request.
        """
        with patch(f"{LOOKUP}.get_overlapping_features_for_region") as overlap:
            result = infer_hgnc_symbol_from_genomic_loci(
                "NC_000002.12", [1, 1 + 6_000_000]
            )
        assert result is None
        overlap.assert_not_called()


class TestGenomicPositionsFromRecords:
    def test_extracts_and_skips(self):
        rows = [
            _row("g.140753336A>T"),
            _row("_wt"),
            _row("_sy"),
            _row("="),
        ]
        positions = _genomic_positions_from_records(rows)
        # 1-based 140753336 -> 0-based 140753335; reference/special rows skipped.
        assert positions == [140753335]

    def test_unparseable_rows_skipped(self):
        positions = _genomic_positions_from_records([_row("not-a-variant")])
        assert positions == []


class TestSelectTranscriptsRouting:
    def _metadata(self, target: TargetGene) -> ScoresetMetadata:
        return ScoresetMetadata(
            urn="urn:mavedb:00000001-a-1", target_genes={"T": target}
        )

    @pytest.mark.asyncio
    async def test_nc_coding_routes_to_genomic_accession_path(self):
        target = _nc_coding_target()
        metadata = self._metadata(target)
        with patch(
            f"{MODULE}._select_genomic_accession_reference",
            return_value=TxSelectResult(
                nm="NM_004333.6",
                np="NP_004324.2",
                start=0,
                is_full_match=True,
                sequence="",
                transcript_mode=TranscriptPriority.MANE_SELECT,
                hgnc_symbol="BRAF",
            ),
        ) as sel:
            result = await select_transcripts(
                metadata, {"T": [_row("g.123A>G")]}, {"T": None}
            )
        assert isinstance(result["T"], TxSelectResult)
        assert result["T"].nm == "NM_004333.6"
        sel.assert_called_once()

    @pytest.mark.asyncio
    async def test_nc_regulatory_returns_none(self):
        target = _nc_coding_target()
        target.target_gene_category = TargetType.REGULATORY
        metadata = self._metadata(target)
        with patch(f"{MODULE}._select_genomic_accession_reference") as sel:
            result = await select_transcripts(
                metadata, {"T": [_row("g.123A>G")]}, {"T": None}
            )
        # Regulatory NC_ target: no coding transcript expected, no typed error.
        assert result["T"] is None
        sel.assert_not_called()

    @pytest.mark.asyncio
    async def test_nc_coding_no_transcript_stores_typed_error(self):
        target = _nc_coding_target()
        metadata = self._metadata(target)
        with patch(
            f"{MODULE}._select_genomic_accession_reference",
            side_effect=NoCodingTranscriptError("no gene"),
        ):
            result = await select_transcripts(
                metadata, {"T": [_row("g.123A>G")]}, {"T": None}
            )
        assert isinstance(result["T"], NoCodingTranscriptError)

    @pytest.mark.asyncio
    async def test_ensembl_protein_target_maps_to_refseq_counterpart(self):
        target = _accession_target("ENSP00000493543.1")
        metadata = self._metadata(target)
        with patch(
            f"{MODULE}._select_refseq_protein_counterpart",
            return_value=TxSelectResult(
                nm="NM_004333.6",
                np="NP_004324.2",
                start=0,
                is_full_match=True,
                sequence="",
                transcript_mode=TranscriptPriority.MANE_SELECT,
                hgnc_symbol="BRAF",
            ),
        ) as sel:
            result = await select_transcripts(metadata, {"T": []}, {"T": None})
        assert result["T"].nm == "NM_004333.6"
        assert result["T"].np == "NP_004324.2"
        sel.assert_called_once_with("ENSP00000493543.1")

    @pytest.mark.asyncio
    async def test_ensembl_protein_target_falls_back_when_no_counterpart(self):
        target = _accession_target("ENSP00000493543.1")
        metadata = self._metadata(target)
        with patch(f"{MODULE}._select_refseq_protein_counterpart", return_value=None):
            result = await select_transcripts(metadata, {"T": []}, {"T": None})
        # No RefSeq counterpart found: falls through to the bare accession passthrough.
        assert result["T"].nm is None
        assert result["T"].np == "ENSP00000493543.1"

    @pytest.mark.asyncio
    async def test_refseq_protein_target_skips_counterpart_lookup(self):
        target = _accession_target("NP_004324.2")
        metadata = self._metadata(target)
        with patch(f"{MODULE}._select_refseq_protein_counterpart") as sel:
            result = await select_transcripts(metadata, {"T": []}, {"T": None})
        assert result["T"].np == "NP_004324.2"
        sel.assert_not_called()

    @pytest.mark.asyncio
    async def test_ensembl_transcript_target_maps_to_refseq_counterpart(self):
        target = _accession_target("ENST00000646891.2")
        metadata = self._metadata(target)
        with patch(
            f"{MODULE}._select_refseq_cdna_counterpart",
            return_value=TxSelectResult(
                nm="NM_004333.6",
                np="NP_004324.2",
                start=0,
                is_full_match=True,
                sequence="",
                transcript_mode=TranscriptPriority.MANE_SELECT,
                hgnc_symbol="BRAF",
            ),
        ) as sel:
            result = await select_transcripts(metadata, {"T": []}, {"T": None})
        assert result["T"].nm == "NM_004333.6"
        assert result["T"].np == "NP_004324.2"
        sel.assert_called_once_with("ENST00000646891.2")

    @pytest.mark.asyncio
    async def test_ensembl_transcript_target_none_when_no_counterpart(self):
        target = _accession_target("ENST00000646891.2")
        metadata = self._metadata(target)
        with patch(f"{MODULE}._select_refseq_cdna_counterpart", return_value=None):
            result = await select_transcripts(metadata, {"T": []}, {"T": None})
        # No RefSeq counterpart found: np is a required TxSelectResult field, so
        # there's no bare-accession object to fall back to -- left None, same as a
        # declared NM_ accession, letting annotation use the declared accession.
        assert result["T"] is None

    @pytest.mark.asyncio
    async def test_refseq_cdna_target_skips_counterpart_lookup(self):
        target = _accession_target("NM_004333.6")
        metadata = self._metadata(target)
        with patch(f"{MODULE}._select_refseq_cdna_counterpart") as sel:
            result = await select_transcripts(metadata, {"T": []}, {"T": None})
        assert result["T"] is None
        sel.assert_not_called()
