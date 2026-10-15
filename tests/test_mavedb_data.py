"""Test `dcd_mapping.resources` module."""

import copy
import json
import shutil
from pathlib import Path

import httpx
import pytest
import respx

from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records


@pytest.fixture
def resources_data_dir():
    """Temporarily store data resources"""
    path = Path(__file__).parent / "tmp"
    if path.exists():  # make sure it's empty
        shutil.rmtree(str(path.absolute()))
    else:
        path.mkdir()
    yield path
    shutil.rmtree(str(path.absolute()))  # clean up afterward


@pytest.fixture
def scoreset_metadata_response(fixture_data_dir: Path):
    """Provide response that client receives from MaveDB API"""
    with (fixture_data_dir / "scoreset_metadata_response.json").open() as f:
        return json.load(f)


def test_get_scoreset_metadata(
    resources_data_dir: Path, scoreset_metadata_response: dict
):
    urn = "urn:mavedb:00000093-a-1"
    with respx.mock:
        respx.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}").mock(
            return_value=httpx.Response(200, json=scoreset_metadata_response[urn])
        )
        scoreset_metadata = get_scoreset_metadata(
            urn, dcd_mapping_dir=resources_data_dir
        )
        assert scoreset_metadata.urn == urn
        assert scoreset_metadata.target_genes
        assert (
            scoreset_metadata.target_genes[
                "BRCA1 translation start through RING domain"
            ].target_uniprot_ref.id
            == "uniprot:P38398"
        )
        assert (
            scoreset_metadata.target_genes[
                "BRCA1 translation start through RING domain"
            ].target_uniprot_ref.offset
            == 0
        )


def test_get_scoreset_metadata_multi_target_uniprot_refs(
    resources_data_dir: Path, scoreset_metadata_response: dict
):
    """Each target in a multi-target score set keeps its own UniProt reference."""
    urn = "urn:mavedb:00000093-a-1"
    response = copy.deepcopy(scoreset_metadata_response[urn])
    first = response["targetGenes"][0]
    first["targetSequence"]["label"] = "first"
    second = copy.deepcopy(first)
    second["targetSequence"]["label"] = "second"
    for ext_id in second["externalIdentifiers"]:
        if ext_id["identifier"]["dbName"] == "UniProt":
            ext_id["identifier"]["identifier"] = "P15056"
            ext_id["offset"] = 5
    response["targetGenes"].append(second)

    with respx.mock:
        respx.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}").mock(
            return_value=httpx.Response(200, json=response)
        )
        scoreset_metadata = get_scoreset_metadata(
            urn, dcd_mapping_dir=resources_data_dir
        )

    first_ref = scoreset_metadata.target_genes["first"].target_uniprot_ref
    second_ref = scoreset_metadata.target_genes["second"].target_uniprot_ref
    assert (first_ref.id, first_ref.offset) == ("uniprot:P38398", 0)
    assert (second_ref.id, second_ref.offset) == ("uniprot:P15056", 5)


def test_get_scoreset_records(
    resources_data_dir: Path, fixture_data_dir: Path, scoreset_metadata_response: dict
):
    urn = "urn:mavedb:00000093-a-1"
    with (fixture_data_dir / f"{urn}_scores.csv").open() as f:
        scores_csv_text = f.read()
    with respx.mock:
        respx.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}").mock(
            return_value=httpx.Response(200, json=scoreset_metadata_response[urn])
        )
        scoreset_metadata = get_scoreset_metadata(
            urn, dcd_mapping_dir=resources_data_dir
        )
        respx.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores").mock(
            return_value=httpx.Response(200, text=scores_csv_text)
        )
        scoreset_records = get_scoreset_records(
            scoreset_metadata, dcd_mapping_dir=resources_data_dir
        )
        assert (
            len(scoreset_records["BRCA1 translation start through RING domain"]) == 853
        )
