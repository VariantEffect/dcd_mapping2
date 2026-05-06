"""Test ``transcripts`` module.

Todo:
----
* Get a test case where there are no common transcripts and you have to use the UniProt ref.


"""
import re
from collections.abc import Coroutine
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pytest
from cool_seq_tool.schemas import TranscriptPriority

from dcd_mapping.mavedb_data import _load_scoreset_records, get_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)
from dcd_mapping.transcripts import (
    TranscriptDescription,
    _choose_most_similar_transcript,
    _percent_similarity,
    _select_protein_reference,
    select_transcript,
)


@pytest.fixture()
def mock_cst(mocker: MagicMock, mock_seqrepo_access):
    """Mock CoolSeqTool instance."""

    async def _execute_query(query: str) -> Coroutine[Any, Any, Any]:
        query = query.strip().replace("\n", "")
        query = re.sub(r"\s+", " ", query)
        calls = {
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202207252 BETWEEN alt_start_i AND alt_end_i OR 202207321 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005782.1"},
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_001005782.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371393.1"},
                {"tx_ac": "NM_001371394.1"},
                {"tx_ac": "NM_003352.8"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202210734 BETWEEN alt_start_i AND alt_end_i OR 202210806 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_001005782.1"},
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005782.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371394.1"},
                {"tx_ac": "NM_003352.8"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202214356 BETWEEN alt_start_i AND alt_end_i OR 202214434 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005782.1"},
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_001005782.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371393.1"},
                {"tx_ac": "NM_001371394.1"},
                {"tx_ac": "NM_003352.8"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202220031 BETWEEN alt_start_i AND alt_end_i OR 202220109 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_003352.8"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371393.1"},
                {"tx_ac": "NM_001371394.1"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37397802 BETWEEN alt_start_i AND alt_end_i OR 37397854 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37400114 BETWEEN alt_start_i AND alt_end_i OR 37400294 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37401601 BETWEEN alt_start_i AND alt_end_i OR 37401678 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%'; ": [
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37402434 BETWEEN alt_start_i AND alt_end_i OR 37402588 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37402748 BETWEEN alt_start_i AND alt_end_i OR 37402880 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37403170 BETWEEN alt_start_i AND alt_end_i OR 37403325 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SCN5A' AND (38551475 BETWEEN alt_start_i AND alt_end_i OR 38551511 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000003.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_198056.2"},
                {"tx_ac": "NM_001160161.1"},
                {"tx_ac": "NM_001354701.1"},
                {"tx_ac": "NM_001099405.1"},
                {"tx_ac": "NM_000335.4"},
                {"tx_ac": "NM_001160160.1"},
                {"tx_ac": "NM_001099404.1"},
                {"tx_ac": "NM_001160160.2"},
                {"tx_ac": "NM_001160161.2"},
                {"tx_ac": "NM_001354701.2"},
                {"tx_ac": "NM_001099404.2"},
                {"tx_ac": "NM_001099405.2"},
                {"tx_ac": "NM_000335.5"},
                {"tx_ac": "NM_198056.3"},
            ],
        }
        return calls[query]

    mock_uta_instance = mocker.MagicMock()
    mock_uta_instance.schema = "20210129b"
    mock_uta_instance.execute_query.side_effect = _execute_query

    mock_cst_instance = mocker.MagicMock()
    mock_cst_instance.uta = mock_uta_instance
    mock_cst_instance.seqrepo_access = mock_seqrepo_access
    mocker.patch(
        "dcd_mapping.lookup.CoolSeqToolBuilder.__new__", return_value=mock_cst_instance
    )

    return mock_cst_instance


def check_transcript_results_equality(actual: TxSelectResult, expected: TxSelectResult):
    """Check equality of transcript selection result vs fixture"""
    assert actual.np == expected.np
    assert actual.start == expected.start
    assert actual.is_full_match is expected.is_full_match
    assert actual.nm == expected.nm
    assert actual.transcript_mode == expected.transcript_mode


@pytest.mark.asyncio(scope="module")
async def test_1_b_2(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    urn = "urn:mavedb:00000001-b-2"
    metadata = scoreset_metadata_fixture[urn]
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    align_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, align_result)
    assert actual, "`select_transcript()` should return a transcript selection result"
    check_transcript_results_equality(actual, expected)


# --- Similarity helper tests ---


def make_mane(nm: str, np: str, priority: TranscriptPriority):
    # Use generic TranscriptDescription for testing similarity logic
    return TranscriptDescription(
        refseq_nuc=nm, refseq_prot=np, transcript_priority=priority
    )


def test_percent_similarity_basic():
    assert _percent_similarity("AAAA", "AAAA") == 1.0
    assert _percent_similarity("AAAA", "BAAAA") == 1.0  # substring fast path
    assert 0.0 <= _percent_similarity("ABCD", "WXYZ") <= 1.0


def test_choose_most_similar_transcript_simple(monkeypatch):
    # Query is most similar to NP_2
    query = "MKTFFV"
    seqs = {
        "NP_1": "MKAAAA",
        "NP_2": "MKTFFV",
        "NP_3": "MKTYFV",
    }

    def fake_get_sequence(ac):
        return seqs[ac]

    monkeypatch.setattr("dcd_mapping.transcripts.get_sequence", fake_get_sequence)

    mane_list = [
        make_mane("NM_1", "NP_1", TranscriptPriority.MANE_PLUS_CLINICAL),
        make_mane("NM_2", "NP_2", TranscriptPriority.MANE_SELECT),
        make_mane("NM_3", "NP_3", TranscriptPriority.MANE_PLUS_CLINICAL),
    ]

    best = _choose_most_similar_transcript(query, mane_list)
    assert best.refseq_prot == "NP_2"


def test_choose_most_similar_transcript_tie_keeps_first(monkeypatch):
    # NP_1 and NP_2 have identical sequences vs query; NP_1 should win (stable tie)
    query = "ABCDE"
    seqs = {
        "NP_1": "ABCDE",
        "NP_2": "ABCDE",
        "NP_3": "ABCXX",
    }

    def fake_get_sequence(ac):
        return seqs[ac]

    monkeypatch.setattr("dcd_mapping.transcripts.get_sequence", fake_get_sequence)

    mane_list = [
        make_mane("NM_1", "NP_1", TranscriptPriority.MANE_PLUS_CLINICAL),
        make_mane("NM_2", "NP_2", TranscriptPriority.MANE_SELECT),
        make_mane("NM_3", "NP_3", TranscriptPriority.MANE_PLUS_CLINICAL),
    ]

    best = _choose_most_similar_transcript(query, mane_list)
    assert best.refseq_prot == "NP_1"


@pytest.mark.asyncio(scope="module")
async def test_end_to_end_per_gene_then_similarity(monkeypatch):
    """E2E: per-gene best via MANE priority, then global similarity among winners."""
    # Mock compatible transcripts grouped by HGNC symbol
    compatible = {
        ("NM_G1_A", "G1"),
        ("NM_G1_B", "G1"),
        ("NM_G2_A", "G2"),
        ("NM_G2_B", "G2"),
    }

    async def fake_get_compatible(_align_result):
        return compatible

    monkeypatch.setattr(
        "dcd_mapping.transcripts._get_compatible_transcripts", fake_get_compatible
    )

    # MANE per gene: G1 has a MANE_SELECT; G2 no select -> plus clinical wins
    def fake_get_mane_transcripts(tx_list):
        s = set(tx_list)
        if s == {"NM_G1_A", "NM_G1_B"}:
            return [
                make_mane("NM_G1_A", "NP_G1_A", TranscriptPriority.MANE_PLUS_CLINICAL),
                make_mane("NM_G1_B", "NP_G1_B", TranscriptPriority.MANE_SELECT),
            ]
        if s == {"NM_G2_A", "NM_G2_B"}:
            return [
                make_mane("NM_G2_A", "NP_G2_A", TranscriptPriority.MANE_PLUS_CLINICAL),
                make_mane(
                    "NM_G2_B",
                    "NP_G2_B",
                    TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
                ),
            ]
        return []

    monkeypatch.setattr(
        "dcd_mapping.transcripts.get_mane_transcripts", fake_get_mane_transcripts
    )

    # Sequence database: make query closest to NP_G2_A (G2 winner)
    seqs = {
        "NP_G1_A": "MKAAAA",
        "NP_G1_B": "MKBBBB",
        "NP_G2_A": "MKTFFV",
        "NP_G2_B": "MKTYFV",
    }

    def fake_get_sequence(ac):
        return seqs[ac]

    monkeypatch.setattr("dcd_mapping.transcripts.get_sequence", fake_get_sequence)

    # Build target gene and minimal align result
    query = "MKTFFV"
    target_gene = TargetGene(
        target_gene_name="Dummy",
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence=query,
        target_sequence_type=TargetSequenceType.PROTEIN,
    )

    align_result = AlignmentResult(
        chrom="chr1",
        strand=1,
        coverage=None,
        percent_identity=None,
        query_range=SequenceRange(start=1, end=6),
        query_subranges=[SequenceRange(start=1, end=6)],
        hit_range=SequenceRange(start=1, end=6),
        hit_subranges=[SequenceRange(start=1, end=6)],
    )

    tx = await _select_protein_reference(target_gene, align_result)
    assert tx.np == "NP_G2_A"


@pytest.mark.asyncio(scope="module")
async def test_tx_src(
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_scn5a(
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000098-a-1"""
    urn = "urn:mavedb:00000098-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_hbb(
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000018-a-1"""
    urn = "urn:mavedb:00000018-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual is None


@pytest.mark.asyncio(scope="module")
async def test_tx_raf(
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000061-h-1"""
    urn = "urn:mavedb:00000061-h-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_tp53(
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000068-a-1"""
    urn = "urn:mavedb:00000068-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)
