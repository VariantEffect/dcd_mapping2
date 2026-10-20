"""Check that the mapper Docker image ships the UCSC liftover chain files.

This is an image check, not a unit test: the Dockerfile runs it during
``docker build``, and it skips everywhere else. Outside the image agct is free
to download the chain files on demand. Inside it, that download goes through
certifi, which fails whenever certifi drops a root UCSC's TLS chain relies on
(VariantEffect/dcd_mapping2#117), so the Dockerfile bakes the files in and this
check fails the build if agct can't find them.
"""

import os
from unittest import mock

import pytest
from agct import Converter

# Only the Dockerfile sets DCD_MAPPING_RESOURCES_DIR, so it marks the image.
pytestmark = pytest.mark.skipif(
    "DCD_MAPPING_RESOURCES_DIR" not in os.environ,
    reason="image check; runs during docker build",
)


@pytest.mark.parametrize(("from_db", "to_db"), [("hg19", "hg38"), ("hg38", "hg19")])
def test_converter_loads_baked_chainfile(from_db: str, to_db: str):
    with mock.patch(
        "agct.converter.download_http",
        side_effect=AssertionError("agct tried to download a chain file"),
    ):
        Converter(from_db, to_db)
