"""Test `dcd_mapping.resources` module."""

import copy
import json
import shutil
from pathlib import Path

import httpx
import pytest
import respx

from dcd_mapping.mavedb_data import (
    get_scoreset_metadata,
    get_scoreset_records,
    patch_target_sequence_type,
)
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMetadata,
    TargetGene,
    TargetSequenceType,
    TargetType,
)


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


def _dna_target(name: str, sequence: str = "ATGGCCAAATAA") -> TargetGene:
    return TargetGene(
        target_gene_name=name,
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence=sequence,
        target_sequence_type=TargetSequenceType.DNA,
    )


def _row(hgvs_nt: str = "NA", hgvs_pro: str = "NA") -> ScoreRow:
    return ScoreRow(
        hgvs_nt=hgvs_nt, hgvs_pro=hgvs_pro, score="1.0", accession="urn:x#1"
    )


# The targets of urn:mavedb:00001237-a-1, in order. Every variant is an n. variant;
# PTPN11, KRAS, NF2 and MSH2 have none.
_ORF_TARGETS = [
    "PTPN11",
    "KRAS",
    "NF2",
    "BCL10",
    "MSH2",
    "RHO",
    "BRAF",
    "FH",
    "CTCF",
    "RPS19",
    "SOS2",
    "LMNA",
    "RET",
    "F9",
    "RAD51D",
    "SHOC2",
    "SDHD",
    "ZC4H2",
    "CCM2",
    "FLCN",
    "FANCC",
    "BRIP1",
    "G6PD",
    "BAP1",
    "BRCA1",
    "ABCD1",
    "RAF1",
    "FARS2",
]
_TARGETS_WITHOUT_VARIANTS = {"PTPN11", "KRAS", "NF2", "MSH2"}


def test_patch_target_sequence_type_ignores_targets_without_variants():
    metadata = ScoresetMetadata(
        urn="urn:mavedb:00001237-a-1",
        target_genes={name: _dna_target(name) for name in _ORF_TARGETS},
    )
    records = {
        name: [_row(hgvs_nt=f"{name}:n.4G>A")]
        for name in _ORF_TARGETS
        if name not in _TARGETS_WITHOUT_VARIANTS
    }

    patch_target_sequence_type(metadata, records)

    assert {t.target_sequence_type for t in metadata.target_genes.values()} == {
        TargetSequenceType.DNA
    }
    assert {t.target_sequence for t in metadata.target_genes.values()} == {
        "ATGGCCAAATAA"
    }


def test_patch_target_sequence_type_converts_only_targets_without_nucleotide_variants():
    metadata = ScoresetMetadata(
        urn="urn:mavedb:00000001-a-1",
        target_genes={
            "PRO": _dna_target("PRO"),
            "NT": _dna_target("NT"),
            "EMPTY": _dna_target("EMPTY"),
        },
    )
    records = {
        "PRO": [
            _row(hgvs_pro="PRO:p.Ala2Val"),
            _row(hgvs_pro="PRO:p.Lys3Arg"),
            _row(hgvs_nt="_wt", hgvs_pro="_wt"),
            _row(hgvs_nt="_sy", hgvs_pro="_sy"),
        ],
        "NT": [_row(hgvs_nt="NT:c.4G>A", hgvs_pro="NT:p.Ala2Thr")],
    }

    patch_target_sequence_type(metadata, records)

    converted = metadata.target_genes["PRO"]
    assert converted.target_sequence_type == TargetSequenceType.PROTEIN
    assert converted.target_sequence == "MAK"
    for label in ("NT", "EMPTY"):
        target = metadata.target_genes[label]
        assert target.target_sequence_type == TargetSequenceType.DNA
        assert target.target_sequence == "ATGGCCAAATAA"
