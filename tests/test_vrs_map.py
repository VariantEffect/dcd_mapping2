"""Test ``vrs_map.py``

* Use 2.0a VA IDs rather than 1.3 IDs. Currently converting to 1.3 and testing b/c
we're focused on remaining consistent w/ previous results.
* Move expected data into a separate JSON file or something?
"""

from collections.abc import Generator
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from cdot.hgvs.dataproviders import ChainedSeqFetcher
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.vrs._internal.models import (
    Allele,
    Haplotype,
    ReferenceLengthExpression,
    SequenceLocation,
    SequenceReference,
)

from dcd_mapping.annotate import _allele_to_v1_allele
from dcd_mapping.mavedb_data import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    ScoreRow,
    ScoresetMetadata,
    SequenceRange,
    TargetGene,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)
from dcd_mapping.transcripts import TxSelectError
from dcd_mapping.vrs_map import (
    NO_TRANSCRIPT_PROTEIN_ERROR,
    PROTEIN_ONLY_ALIGNMENT_ERROR,
    _construct_vrs_allele,
    _map_protein_layer,
    vrs_map,
)


def _assert_correct_vrs_map(
    mapping: MappedScore,
    expected_mappings_data: dict[tuple[str, AnnotationLayer], dict],
):
    """Note that we're testing against VRS 1.3 VA IDs (temporary?)."""
    key = (mapping.accession_id, mapping.alignment_level)
    assert key in expected_mappings_data, (
        "Score row/layer combination is not in expected mappings"
    )
    expected = expected_mappings_data[key]
    vrs_2_to_1 = lambda var: _allele_to_v1_allele(var)  # noqa: E731
    if isinstance(mapping.pre_mapped, Haplotype) and isinstance(
        mapping.post_mapped, Haplotype
    ):
        assert all(
            len(x) == len(mapping.pre_mapped.members)
            for x in (
                mapping.post_mapped.members,
                expected["pre_mapped"],
                expected["post_mapped"],
            )
        ), "mappings are different lengths"
        for va_id in expected["pre_mapped"]:
            for variant in mapping.pre_mapped.members:
                if vrs_2_to_1(variant).id == va_id:
                    break
            else:
                pytest.fail(f"Failed to find {va_id} in pre-mapped variants.")
        for va_id in expected["post_mapped"]:
            for variant in mapping.post_mapped.members:
                if vrs_2_to_1(variant).id == va_id:
                    break
            else:
                pytest.fail(f"Failed to find {va_id} in post-mapped variants.")
    elif isinstance(mapping.pre_mapped, Allele) and isinstance(
        mapping.post_mapped, Allele
    ):
        assert vrs_2_to_1(mapping.pre_mapped).id == expected["pre_mapped"]
        assert vrs_2_to_1(mapping.post_mapped).id == expected["post_mapped"]
    else:
        pytest.fail("mapping format appears to be broken")


@pytest.fixture
def get_fixtures(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    def _get_fixtures(urn: str):
        return (
            _load_scoreset_records(
                fixture_data_dir / f"{urn}_scores.csv", scoreset_metadata_fixture[urn]
            ),
            scoreset_metadata_fixture[urn],
            align_result_fixture[urn],
            transcript_results_fixture[urn],
        )

    return _get_fixtures


def test_2_a_2(
    get_fixtures,
    mock_seqrepo_access: MagicMock,
    data_provider: Generator[ChainedSeqFetcher, None, None],
):
    urn = "urn:mavedb:00000002-a-2"
    records, metadata, align_result, tx_result = get_fixtures(urn)
    expected_mappings_data = {
        ("urn:mavedb:00000002-a-2#1", AnnotationLayer.PROTEIN): {
            "pre_mapped": [
                "ga4gh:VA.jvd3wir-9AwP7Ay9FWnqqGVxETG9Dl0M",
                "ga4gh:VA.aKrgTa26FUdF8b4wMk7mwFCZnTQmIe5i",
            ],
            "post_mapped": [
                "ga4gh:VA.-IuFaR_wBHzEQJSdDAoke-r2LRe9jtMQ",
                "ga4gh:VA.aF9h1d9DvWWGlkhRAbdz1Ni9DQUOXIhL",
            ],
        },
    }

    mappings = {}
    with patch("dcd_mapping.lookup.seqfetcher") as mock_cdot_seqfetcher:
        mock_cdot_seqfetcher.return_value = data_provider
        for target_gene in metadata.target_genes:
            mappings[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=align_result[target_gene],
                records=records[target_gene],
                transcript=tx_result[target_gene],
            ).mappings
    assert mappings["hYAP65 WW domain"] is not None
    assert len(mappings["hYAP65 WW domain"]) == 1

    for m in mappings["hYAP65 WW domain"]:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        [
            "DVPLPAGWEMAKTSSGQRYFLNHIDQTTTWQDPR",
            [{"namespace": "ga4gh", "alias": "SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R"}],
        ],
        [
            "GACGTTCCACTGCCGGCTGGTTGGGAAATGGCTAAAACTAGTTCTGGTCAGCGTTACTTCCTGAACCACATCGACCAGACCACCACGTGGCAGGACCCGCGT",
            [{"namespace": "ga4gh", "alias": "SQ.CR4vhHED4FbR2mgu0k14tfJ0ldnMoi2A"}],
        ],
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_41_a_1(
    get_fixtures,
    mock_seqrepo_access: MagicMock,
    data_provider: Generator[ChainedSeqFetcher, None, None],
):
    urn = "urn:mavedb:00000041-a-1"
    records, metadata, align_result, tx_result = get_fixtures(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000041-a-1#548", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.NJgaCF0JPFERdw9Y7fW4bXkaP1tSa5fv",
            "post_mapped": "ga4gh:VA.csCB31gWoiiD38TlR35dZnyAI156YWgW",
        },
        ("urn:mavedb:00000041-a-1#50", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.RfNyaPcZg8o9f3YK0qa4Vu3LwBEQdmVf",
            "post_mapped": "ga4gh:VA.QP9KLxvN6_b7sWeY7L8FBs5XKMxsCiLE",
        },
        ("urn:mavedb:00000041-a-1#51", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.sZNa3SNPlv_gU2JSiH7Q03nNT7oFy1NX",
            "post_mapped": "ga4gh:VA.2kmNV4T4Bp_UAnI002QOfrzd_yqb21vs",
        },
        ("urn:mavedb:00000041-a-1#977", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.h7QtWm0WzlOr0zbB9y4tJOIEUet6VLcB",
            "post_mapped": "ga4gh:VA.bakDMkeUFIa46_HAwXb7gUjhywggVNIN",
        },
        ("urn:mavedb:00000041-a-1#52", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.gl5xiWNmwUfEMZe5Aub15HIiztUaizay",
            "post_mapped": "ga4gh:VA.3Pp5-tRnYkmm8f6qxk06GvTpn81DqiQV",
        },
    }

    mappings = {}
    with patch("dcd_mapping.lookup.seqfetcher") as mock_cdot_seqfetcher:
        mock_cdot_seqfetcher.return_value = data_provider
        for target_gene in metadata.target_genes:
            mappings[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=align_result[target_gene],
                records=records[target_gene],
                transcript=tx_result[target_gene],
            ).mappings
    assert mappings["Src catalytic domain"] is not None
    assert len(mappings["Src catalytic domain"]) == 5

    for m in mappings["Src catalytic domain"]:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        [
            "LRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKEPEERPTFEYLQAFL",
            [{"namespace": "ga4gh", "alias": "SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB"}],
        ],
        [
            "CTGCGGCTGGAGGTCAAGCTGGGCCAGGGCTGCTTTGGCGAGGTGTGGATGGGGACCTGGAACGGTACCACCAGGGTGGCCATCAAAACCCTGAAGCCTGGCACGATGTCTCCAGAGGCCTTCCTGCAGGAGGCCCAGGTCATGAAGAAGCTGAGGCATGAGAAGCTGGTGCAGTTGTATGCTGTGGTTTCAGAGGAGCCCATTTACATCGTCACGGAGTACATGAGCAAGGGGAGTTTGCTGGACTTTCTCAAGGGGGAGACAGGCAAGTACCTGCGGCTGCCTCAGCTGGTGGACATGGCTGCTCAGATCGCCTCAGGCATGGCGTACGTGGAGCGGATGAACTACGTCCACCGGGACCTTCGTGCAGCCAACATCCTGGTGGGAGAGAACCTGGTGTGCAAAGTGGCCGACTTTGGGCTGGCTCGGCTCATTGAAGACAATGAGTACACGGCGCGGCAAGGTGCCAAATTCCCCATCAAGTGGACGGCTCCAGAAGCTGCCCTCTATGGCCGCTTCACCATCAAGTCGGACGTGTGGTCCTTCGGGATCCTGCTGACTGAGCTCACCACAAAGGGACGGGTGCCCTACCCTGGGATGGTGAACCGCGAGGTGCTGGACCAGGTGGAGCGGGGCTACCGGATGCCCTGCCCGCCGGAGTGTCCCGAGTCCCTGCACGACCTCATGTGCCAGTGCTGGCGGAAGGAGCCTGAGGAGCGGCCCACCTTCGAGTACCTGCAGGCCTTCCTG",
            [{"namespace": "ga4gh", "alias": "SQ.SBlyhKAQ5u14huYar7I7UmXZk4eRPexx"}],
        ],
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_99_a_1(
    get_fixtures,
    mock_seqrepo_access: MagicMock,
    data_provider: Generator[ChainedSeqFetcher, None, None],
):
    urn = "urn:mavedb:00000099-a-1"
    records, metadata, align_result, tx_result = get_fixtures(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000099-a-1#8", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.ES3p2PRZ6Lhial9xnlVxV8uXUbl03ESx",
            "post_mapped": "ga4gh:VA.ES3p2PRZ6Lhial9xnlVxV8uXUbl03ESx",
        },
        ("urn:mavedb:00000099-a-1#8", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.dw7KwkqdRsRL5UU55KznGYFqGCK6sWgt",
            "post_mapped": "ga4gh:VA.T8SkVZHuAF9J23etkAXG2Sz2w2yJMSO4",
        },
        ("urn:mavedb:00000099-a-1#96", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.bu9RmB2kvJiTFKfThhtpYam50_GrIkbx",
            "post_mapped": "ga4gh:VA.bu9RmB2kvJiTFKfThhtpYam50_GrIkbx",
        },
        ("urn:mavedb:00000099-a-1#96", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.4cQNoG_u6yGDk_GMmVvbYX-8EniN6RKd",
            "post_mapped": "ga4gh:VA.4-NsddlZX70Jzi7OxSL6KYRe8T-h4rwT",
        },
        ("urn:mavedb:00000099-a-1#194", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.eT6tN7HdCmDQnTHuetRHtRZRHPxpO53x",
            "post_mapped": "ga4gh:VA.eT6tN7HdCmDQnTHuetRHtRZRHPxpO53x",
        },
        ("urn:mavedb:00000099-a-1#194", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.gRO5ANz9hJ6y6NuyzHOzl94b39SAbbIE",
            "post_mapped": "ga4gh:VA.JSmTABIybCChEE6mmjnPuubOt9eRWj_4",
        },
        ("urn:mavedb:00000099-a-1#211", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.MdeuhDWe0mccL5lPkb5J3d_I39m-Zg1R",
            "post_mapped": "ga4gh:VA.MdeuhDWe0mccL5lPkb5J3d_I39m-Zg1R",
        },
        ("urn:mavedb:00000099-a-1#211", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.K4bnyLQ5M3WCmli8lvAqHfqkZuk4KvwP",
            "post_mapped": "ga4gh:VA.HSfipwsg28LbqwITCawzumz_OWZYu_jM",
        },
    }
    mappings = {}
    with patch("dcd_mapping.lookup.seqfetcher") as mock_cdot_seqfetcher:
        mock_cdot_seqfetcher.return_value = data_provider
        for target_gene in metadata.target_genes:
            mappings[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=align_result[target_gene],
                records=records[target_gene],
                transcript=tx_result[target_gene],
            ).mappings
    assert mappings["RHO"] is not None
    assert len(mappings["RHO"]) == 8  # includes protein and genomic for all 4 rows

    for m in mappings["RHO"]:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA",
            [{"namespace": "ga4gh", "alias": "SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9"}],
        ),
        (
            "ATGAATGGCACAGAAGGCCCTAACTTCTACGTGCCCTTCTCCAATGCGACGGGTGTGGTACGCAGCCCCTTCGAGTACCCACAGTACTACCTGGCTGAGCCATGGCAGTTCTCCATGCTGGCCGCCTACATGTTTCTGCTGATCGTGCTGGGCTTCCCCATCAACTTCCTCACGCTCTACGTCACCGTCCAGCACAAGAAGCTGCGCACGCCTCTCAACTACATCCTGCTCAACCTAGCCGTGGCTGACCTCTTCATGGTCCTAGGTGGCTTCACCAGCACCCTCTACACCTCTCTGCATGGATACTTCGTCTTCGGGCCCACAGGATGCAATTTGGAGGGCTTCTTTGCCACCCTGGGCGGTGAAATTGCCCTGTGGTCCTTGGTGGTCCTGGCCATCGAGCGGTACGTGGTGGTGTGTAAGCCCATGAGCAACTTCCGCTTCGGGGAGAACCATGCCATCATGGGCGTTGCCTTCACCTGGGTCATGGCGCTGGCCTGCGCCGCACCCCCACTCGCCGGCTGGTCCAGGTACATCCCCGAGGGCCTGCAGTGCTCGTGTGGAATCGACTACTACACGCTCAAGCCGGAGGTCAACAACGAGTCTTTTGTCATCTACATGTTCGTGGTCCACTTCACCATCCCCATGATTATCATCTTTTTCTGCTATGGGCAGCTCGTCTTCACCGTCAAGGAGGCCGCTGCCCAGCAGCAGGAGTCAGCCACCACACAGAAGGCAGAGAAGGAGGTCACCCGCATGGTCATCATCATGGTCATCGCTTTCCTGATCTGCTGGGTGCCCTACGCCAGCGTGGCATTCTACATCTTCACCCACCAGGGCTCCAACTTCGGTCCCATCTTCATGACCATCCCAGCGTTCTTTGCCAAGAGCGCCGCCATCTACAACCCTGTCATCTATATCATGATGAACAAGCAGTTCCGGAACTGCATGCTCACCACCATCTGCTGCGGCAAGAACCCACTGGGTGACGATGAGGCCTCTGCTACCGTGTCCAAGACGGAGACGAGCCAGGTGGCCCCGGCCTAA",
            [{"namespace": "ga4gh", "alias": "SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_103_c_1(
    get_fixtures,
    mock_seqrepo_access: MagicMock,
    data_provider: Generator[ChainedSeqFetcher, None, None],
):
    urn = "urn:mavedb:00000103-c-1"
    records, metadata, align_result, tx_result = get_fixtures(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000103-c-1#376", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.KPfJzFOpDl49yIPqqeF07BsnCsOPbf5Z",
            "post_mapped": "ga4gh:VA.KPfJzFOpDl49yIPqqeF07BsnCsOPbf5Z",
        },
        ("urn:mavedb:00000103-c-1#55", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.AVg1O4zA7DP71z-QZD7-5864A52Cp4RO",
            "post_mapped": "ga4gh:VA.AVg1O4zA7DP71z-QZD7-5864A52Cp4RO",
        },
        ("urn:mavedb:00000103-c-1#2548", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.4dd3ml7ZuqyoMhwJhMOWi2n1729MtE-6",
            "post_mapped": "ga4gh:VA.4dd3ml7ZuqyoMhwJhMOWi2n1729MtE-6",
        },
        ("urn:mavedb:00000103-c-1#6810", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.H_7xCjw--slTAAKFwP42WCUo1tITw39d",
            "post_mapped": "ga4gh:VA.H_7xCjw--slTAAKFwP42WCUo1tITw39d",
        },
    }

    mappings = {}
    with patch("dcd_mapping.lookup.seqfetcher") as mock_cdot_seqfetcher:
        mock_cdot_seqfetcher.return_value = data_provider
        for target_gene in metadata.target_genes:
            mappings[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=align_result[target_gene],
                records=records[target_gene],
                transcript=tx_result[target_gene],
            ).mappings
    assert mappings["MAPK1"] is not None
    assert len(mappings["MAPK1"]) == 4
    for m in mappings["MAPK1"]:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQHLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQEDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLEQYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS",
            [{"namespace": "ga4gh", "alias": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_1_b_2(
    get_fixtures,
    mock_seqrepo_access: MagicMock,
    data_provider: Generator[ChainedSeqFetcher, None, None],
):
    urn = "urn:mavedb:00000001-b-2"
    records, metadata, align_result, tx_result = get_fixtures(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000001-b-2#444", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.ojIs4GEPbxiMt5E6InnF6k6m9ix_z3SH",
            "post_mapped": "ga4gh:VA.ojIs4GEPbxiMt5E6InnF6k6m9ix_z3SH",
        },
        ("urn:mavedb:00000001-b-2#444", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.1XfnARSJMWLoHTQDJU-m0ZYLUsj-qRNM",
            "post_mapped": "ga4gh:VA.ifSwfAlXaZWIcTQaXrZv7LCsa6sywRvw",
        },
        ("urn:mavedb:00000001-b-2#57", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.EwrANE7HrTJ3OidbdJcL4DAfiYiAJ0Zn",
            "post_mapped": "ga4gh:VA.EwrANE7HrTJ3OidbdJcL4DAfiYiAJ0Zn",
        },
        ("urn:mavedb:00000001-b-2#57", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.-zYenjM1Wsuu5Ia06nn856cn4JKEjMwR",
            "post_mapped": "ga4gh:VA.LODyOWFdnBsGn7dciC-MCSLHCNtoyjaf",
        },
        ("urn:mavedb:00000001-b-2#2311", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.KQrLWh1WBJNOO7tEkJ7ujkkHCEGGTSCo",
            "post_mapped": "ga4gh:VA.KQrLWh1WBJNOO7tEkJ7ujkkHCEGGTSCo",
        },
        ("urn:mavedb:00000001-b-2#2311", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.kdSdynoQvgoau0GqsUzhduF7riQVvx5-",
            "post_mapped": "ga4gh:VA.LdKB-BHsNMueerA0u4ngGU1oxK2oBDHs",
        },
        ("urn:mavedb:00000001-b-2#2312", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.ZdHxEMecv2kdYzRJz5iaqiS5qUd05Wim",
            "post_mapped": "ga4gh:VA.ZdHxEMecv2kdYzRJz5iaqiS5qUd05Wim",
        },
        ("urn:mavedb:00000001-b-2#2312", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.pkEhKe6fG7eWaUTA7DJG2v-zggeP5N1m",
            "post_mapped": "ga4gh:VA.R4iEL0X_2Mr4o4cmuLE4AW_UrTotjZ8M",
        },
    }

    mappings = {}
    with patch("dcd_mapping.lookup.seqfetcher") as mock_cdot_seqfetcher:
        mock_cdot_seqfetcher.return_value = data_provider
        for target_gene in metadata.target_genes:
            mappings[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=align_result[target_gene],
                records=records[target_gene],
                transcript=tx_result[target_gene],
            ).mappings
    assert mappings["SUMO1"] is not None
    assert len(mappings["SUMO1"]) == 8
    for m in mappings["SUMO1"]:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MSDQEAKPSTEDLGDKKEGEYIKLKVIGQDSSEIHFKVKMTTHLKKLKESYCQRQGVPMNSLRFLFEGQRIADNHTPKELGMEEEDVIEVYQEQTGGHSTV",
            [{"namespace": "ga4gh", "alias": "SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip"}],
        ),
        (
            "ATGTCTGACCAGGAGGCAAAACCTTCAACTGAGGACTTGGGGGATAAGAAGGAAGGTGAATATATTAAACTCAAAGTCATTGGACAGGATAGCAGTGAGATTCACTTCAAAGTGAAAATGACAACACATCTCAAGAAACTCAAAGAATCATACTGTCAAAGACAGGGTGTTCCAATGAATTCACTCAGGTTTCTCTTTGAGGGTCAGAGAATTGCTGATAATCATACTCCAAAAGAACTGGGAATGGAGGAAGAAGATGTGATTGAAGTTTATCAGGAACAAACGGGGGGTCATTCAACAGTTTAG",
            [{"namespace": "ga4gh", "alias": "SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


class TestMapProteinLayerReason:
    """``_map_protein_layer`` returns ``(mapping, reason)``: no mapping plus a specific
    reason when there is no protein layer to map. A row that maps at no layer is failed
    once, layer-agnostically, carrying this detailed reason rather than a failure
    arbitrarily tagged to the protein layer.
    """

    def _row(self, hgvs_pro: str) -> ScoreRow:
        return ScoreRow(
            hgvs_nt="_wt",
            hgvs_pro=hgvs_pro,
            score="1.0",
            accession="urn:mavedb:00000001-a-2#1",
        )

    def _tx(self) -> TxSelectResult:
        return TxSelectResult(
            nm="NM_000001.1",
            np="NP_000001.1",
            start=0,
            is_full_match=True,
            sequence="MAA",
            hgnc_symbol="GENE1",
        )

    def _align(self) -> AlignmentResult:
        return AlignmentResult(
            chrom="chr1",
            query_range=SequenceRange(start=0, end=9),
            query_subranges=[SequenceRange(start=0, end=9)],
            hit_range=SequenceRange(start=1000, end=1009),
            hit_subranges=[SequenceRange(start=1000, end=1009)],
            aligner_parameters={"aligner": "blat"},
            query_sequence_type=TargetSequenceType.DNA,
        )

    def test_transcript_error_surfaces_its_message(self):
        err = TxSelectError("no transcript selected for GENE1")
        mapping, reason = _map_protein_layer(self._row("p.Ala2Val"), "SQ.x", err, None)
        assert mapping is None
        assert reason == "no transcript selected for GENE1"

    @pytest.mark.parametrize("hgvs_pro", ["_wt", "_sy"])
    def test_invalid_protein_variant_names_inputs(self, hgvs_pro):
        mapping, reason = _map_protein_layer(
            self._row(hgvs_pro), "SQ.x", self._tx(), self._align()
        )
        assert mapping is None
        assert "Can't process variant syntax" in reason
        assert hgvs_pro in reason  # names the actual input, not a layer

    def test_valid_protein_without_alignment_reports_alignment_failure(self):
        mapping, reason = _map_protein_layer(
            self._row("p.Ala2Val"), "SQ.x", self._tx(), None
        )
        assert mapping is None
        assert "could not be aligned" in reason


def _whole_sequence_allele(length: int = 648) -> Allele:
    return Allele(
        location=SequenceLocation(
            sequenceReference=SequenceReference(
                refgetAccession="SQ._YJMuvFaxcDIHC_IXWkSjTtNDx4CIWLU"
            ),
            start=0,
            end=length,
        ),
        state=ReferenceLengthExpression(length=length, repeatSubunitLength=length),
    )


@pytest.mark.parametrize("pre_map", [True, False])
def test_reference_identical_allele_keeps_the_whole_sequence(pre_map):
    """A ``p.=`` allele must not go through VRS normalization, which rewrites an
    unchanged sequence as ``length=0``: a deletion of the whole sequence.
    """
    with (
        patch(
            "dcd_mapping.vrs_map.build_ref_identical_allele",
            return_value=_whole_sequence_allele(),
        ),
        patch(
            "dcd_mapping.vrs_map.translate_ref_identical_to_vrs",
            return_value=_whole_sequence_allele(),
        ),
        patch("dcd_mapping.vrs_map.normalize") as normalize,
    ):
        allele = _construct_vrs_allele(
            ["NP_002871.1:p.="],
            AnnotationLayer.PROTEIN,
            "SQ._YJMuvFaxcDIHC_IXWkSjTtNDx4CIWLU",
            pre_map=pre_map,
        )

    normalize.assert_not_called()
    assert allele.state.length == allele.location.end - allele.location.start == 648
    # The digest MaveDB's backfilled whole-sequence allele for NP_002871.1:p.= carries.
    assert allele.id == "ga4gh:VA.JRYlMMc2SCGXI6dQ1yGyONWLkL1wMa4d"


class TestProteinOnlyAlignment:
    """Modeled on urn:mavedb:00000639-0-1: a codon-optimized DNA target that aligned to
    the genome only as protein. Nucleotide positions can't be placed through an alignment
    whose query blocks count amino acids, so no row gets a genomic record.
    """

    target = TargetGene(
        target_gene_name="ABLIM2 (R14S;R16S)",
        target_gene_category=TargetType.PROTEIN_CODING,
        target_sequence="YKIYPYDSLIVTNSISVKLPKDVDRTRLERHLSPEEFQEVFGMSIEEFDRLALWKRNDLKKKALLF",
        target_sequence_type=TargetSequenceType.PROTEIN,
    )
    # Blocks from the stored production alignment, CIGAR 2S81M16160D111M.
    align_result = AlignmentResult(
        chrom="chr4",
        query_range=SequenceRange(start=2, end=66),
        query_subranges=[
            SequenceRange(start=2, end=29),
            SequenceRange(start=29, end=66),
        ],
        hit_range=SequenceRange(start=7966992, end=7983344),
        hit_subranges=[
            SequenceRange(start=7983263, end=7983344),
            SequenceRange(start=7966992, end=7967103),
        ],
        aligner_parameters={"aligner": "blat", "target_args": "-q=prot -t=dnax"},
        query_sequence_type=TargetSequenceType.PROTEIN,
    )
    transcript = TxSelectResult(
        nm="NM_032432.6",
        np="NP_115808.1",
        start=0,
        is_full_match=False,
        sequence="MSAVSQPQAAPSPLEKSPSTAILCNTCGNVCKGEVLRVQDKYFHIKCFVCKACGCDLAEGGFFVRQGEYICTLDYQRLYGTRCFSCDQFIEGEVVSALGKTYHPDCFVCAVCRLPFPPGDRVTFNGKECMCQKCSLPVSVGSSAHLSQGLRSCGGCGTEIKNGQALVALDKHWHLGCFKCKSCGKLLNAEYISKDGLPYCEADYHAKFGIRCDSCEKYITGRVLEAGEKHYHPSCALCVRCGQMFAEGEEMYLQGSSIWHPACRQAARTEDRNKETRTSSESIISVPASSTSGSPSRVIYAKLGGEILDYRDLAALPKSKAIYDIDRPDMISYSPYISHSAGDRQSYGEGDQDDRSYKQCRTSSPSSTGSVSLGRYTPTSRSPQHYSRPGSESGRSTPSLSVLSDSKPPPSTYQQAPRHFHVPDTGVKDNIYRKPPIYRQHAARRSDGEDGSLDQDNRKKSSWLMLKGDADTRTNSPDLDTQSLSHSSGTDRDPLQRMAGDSFHSRFPYSKSDPLPGHGKNGLDQRNANLAPCGADPDASWGMREYKIYPYDSLIVTNRIRVKLPKDVDRTRLERHLSPEEFQEVFGMSIEEFDRLALWKRNDLKKKALLF",
        hgnc_symbol="ABLIM2",
    )

    def _map(self, rows: list[ScoreRow]):
        with (
            patch("dcd_mapping.vrs_map.store_sequence", return_value="SQ.x"),
            patch("dcd_mapping.vrs_map.align_target_to_protein") as align_protein,
            patch(
                "dcd_mapping.vrs_map._map_protein_coding_pro",
                side_effect=lambda row, *_: MappedScore(
                    accession_id=row.accession,
                    score=row.score,
                    alignment_level=AnnotationLayer.PROTEIN,
                ),
            ),
            patch("dcd_mapping.vrs_map._map_genomic") as map_genomic,
            patch("dcd_mapping.vrs_map._construct_projected_layers") as project,
        ):
            align_protein.return_value = self.align_result
            result = vrs_map(self.target, self.align_result, rows, self.transcript)
        map_genomic.assert_not_called()
        project.assert_not_called()
        return result.mappings

    def test_rows_with_a_protein_variant_map_only_at_the_protein_layer(self):
        rows = [
            ScoreRow(hgvs_nt="c.=", hgvs_pro="p.=", score="0.1", accession="urn:x#1"),
            ScoreRow(
                hgvs_nt="c.[1T>C;3C>G;81T>C]",
                hgvs_pro="p.Tyr1Gln",
                score="0.2",
                accession="urn:x#2",
            ),
        ]

        mappings = self._map(rows)

        assert [(m.accession_id, m.alignment_level) for m in mappings] == [
            ("urn:x#1", AnnotationLayer.PROTEIN),
            ("urn:x#2", AnnotationLayer.PROTEIN),
        ]

    def test_nucleotide_only_row_fails_once_naming_the_protein_alignment(self):
        rows = [
            ScoreRow(
                hgvs_nt="c.150C>T", hgvs_pro="NA", score="0.3", accession="urn:x#3"
            )
        ]

        mappings = self._map(rows)

        assert len(mappings) == 1
        assert mappings[0].alignment_level is None
        assert mappings[0].error_message == PROTEIN_ONLY_ALIGNMENT_ERROR

    def test_rows_without_a_transcript_fail_naming_their_own_layer(self):
        rows = [
            ScoreRow(
                hgvs_nt="NA", hgvs_pro="p.Tyr1Gln", score="0.1", accession="urn:x#1"
            ),
            ScoreRow(
                hgvs_nt="c.150C>T", hgvs_pro="NA", score="0.3", accession="urn:x#2"
            ),
        ]

        with patch("dcd_mapping.vrs_map._map_genomic") as map_genomic:
            mappings = vrs_map(self.target, self.align_result, rows, None).mappings
        map_genomic.assert_not_called()

        assert [m.error_message for m in mappings] == [
            NO_TRANSCRIPT_PROTEIN_ERROR,
            PROTEIN_ONLY_ALIGNMENT_ERROR,
        ]
