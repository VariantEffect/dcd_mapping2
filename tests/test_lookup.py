"""Tests for dcd_mapping.lookup"""

from unittest.mock import patch

import requests

from dcd_mapping.lookup import get_overlapping_features_for_region

RAW_OVERLAP_RESPONSE = [
    {
        "seq_region_name": "22",
        "version": 1,
        "biotype": "protein_coding",
        "feature_type": "gene",
        "description": "novel transcript",
        "logic_name": "havana_homo_sapiens",
        "start": 19717220,
        "id": "ENSG00000284874",
        "source": "havana",
        "canonical_transcript": "ENST00000455843.5",
        "assembly_name": "GRCh38",
        "end": 19724772,
        "gene_id": "ENSG00000284874",
        "strand": 1,
    },
    {
        "source": "ensembl_havana",
        "canonical_transcript": "ENST00000366425.4",
        "assembly_name": "GRCh38",
        "end": 19724776,
        "gene_id": "ENSG00000203618",
        "strand": 1,
        "external_name": "GP1BB",
        "seq_region_name": "22",
        "version": 7,
        "biotype": "protein_coding",
        "logic_name": "ensembl_havana_gene_homo_sapiens",
        "feature_type": "gene",
        "start": 19723539,
        "description": "glycoprotein Ib platelet subunit beta [Source:HGNC Symbol;Acc:HGNC:4440]",
        "id": "ENSG00000203618",
    },
    {
        "end": 19724224,
        "gene_id": "ENSG00000184702",
        "strand": 1,
        "canonical_transcript": "ENST00000455784.7",
        "source": "ensembl_havana",
        "assembly_name": "GRCh38",
        "seq_region_name": "22",
        "version": 21,
        "biotype": "protein_coding",
        "description": "septin 5 [Source:HGNC Symbol;Acc:HGNC:9164]",
        "feature_type": "gene",
        "logic_name": "ensembl_havana_gene_homo_sapiens",
        "start": 19714467,
        "id": "ENSG00000184702",
        "external_name": "SEPTIN5",
    },
]


class _FakeResponse:
    def __init__(self, data):
        self._data = data
        self.status_code = 200

    def json(self):
        return self._data

    def raise_for_status(self):
        return None


def test_get_overlapping_features_for_region_success():
    with (
        patch(
            "dcd_mapping.lookup.request_with_backoff",
            return_value=_FakeResponse(RAW_OVERLAP_RESPONSE),
        ),
        patch("dcd_mapping.lookup.get_chromosome_identifier", side_effect=lambda c: c),
    ):
        result = get_overlapping_features_for_region(
            "NC_000022.11", 19714000, 19725000, features=["gene"]
        )
        assert isinstance(result, list)
        assert result == RAW_OVERLAP_RESPONSE


def test_get_overlapping_features_for_region_error():
    class ErrorResponse(_FakeResponse):
        def __init__(self):
            super().__init__(None)
            self.status_code = 500

        def raise_for_status(self):
            msg = f"HTTP {self.status_code} Error"
            raise requests.RequestException(msg)

    with (
        patch("dcd_mapping.lookup.request_with_backoff", return_value=ErrorResponse()),
        patch("dcd_mapping.lookup.get_chromosome_identifier", side_effect=lambda c: c),
    ):
        result = get_overlapping_features_for_region(
            "NC_000022.11", 19714000, 19725000, features=["gene"]
        )
        assert result == []
