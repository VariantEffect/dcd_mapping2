"""Test ``vrs_map.py``

Todo:
----
* Sample a second mapping row for each test case (also, order shouldn't matter)
* Add a few more test cases


"""
from pathlib import Path
from typing import Dict

import pytest

from dcd_mapping.lookup import get_seqrepo
from dcd_mapping.mavedb_data import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    TxSelectResult,
)

# VrsMapping,
from dcd_mapping.vrs_map import vrs_map


def test_41_a_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    urn = "urn:mavedb:00000041-a-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_select = transcript_results_fixture[urn]

    _ = vrs_map(get_seqrepo(), metadata, records, align_result, tx_select)
    pytest.fail("test")
