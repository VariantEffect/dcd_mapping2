"""Handle requests for MaveDB data, such as scoresets or scoreset metadata.

Much of this can/should be replaced by the ``mavetools`` library? (and/or ``wags-tails``.)
"""

import csv
import json
import logging
import tempfile
import zipfile
from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import Any

import requests
from fastapi import HTTPException
from pydantic import ValidationError

from dcd_mapping.lookup import DataLookupError
from dcd_mapping.resource_utils import (
    LOCAL_STORE_PATH,
    MAVEDB_BASE_URL,
    ResourceAcquisitionError,
    authentication_header,
    http_download,
)
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMapping,
    ScoresetMetadata,
    TargetGene,
    UniProtRef,
)

__all__ = [
    "get_scoreset_urns",
    "get_human_urns",
    "get_raw_scoreset_metadata",
    "get_scoreset_records",
    "get_scoreset_metadata",
    "get_human_urns",
]

_logger = logging.getLogger(__name__)


class ScoresetNotSupportedError(Exception):
    """Raise when a score set cannot be mapped because it has characteristics that are not currently supported."""


def get_scoreset_urns() -> set[str]:
    """Fetch all scoreset URNs. Since species is annotated at the scoreset target level,
    we can't yet filter on anything like `homo sapien` -- meaning this is fairly slow.

    :return: set of URN strings
    """
    r = requests.get(
        f"{MAVEDB_BASE_URL}/api/v1/experiments/",
        timeout=30,
        headers=authentication_header(),
    )
    r.raise_for_status()
    scoreset_urn_lists = [
        experiment["scoreSetUrns"]
        for experiment in r.json()
        if experiment.get("scoreSetUrns")
    ]
    return {urn for urns in scoreset_urn_lists for urn in urns}


def _metadata_response_is_human(json_response: dict) -> bool:
    """Check that response from scoreset metadata API refers to a human genome target.

    :param json_response: response from scoreset metadata API
    :return: True if contains a target tagged as ``"Homo sapiens"``
    """
    for target_gene in json_response.get("targetGenes", []):
        organism = (
            target_gene.get("targetSequence", {})
            .get("taxonomy", {})
            .get("organismName")
        )
        if organism == "Homo sapiens":
            return True
    return False


def get_human_urns() -> list[str]:
    """Fetch all human scoreset URNs. Pretty slow, shouldn't be used frequently because
    it requires requesting every single scoreset.

    :return: list of human scoreset URNs
    """
    scoreset_urns = get_scoreset_urns()
    human_scoresets: list[str] = []
    for urn in scoreset_urns:
        r = requests.get(
            f"{MAVEDB_BASE_URL}/api/v1/score-sets/{urn}",
            timeout=30,
            headers=authentication_header(),
        )
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError:
            _logger.info("Unable to retrieve scoreset data for URN %s", urn)
            continue
        data = r.json()

        if _metadata_response_is_human(data):
            human_scoresets.append(urn)
    return human_scoresets


def _get_uniprot_ref(scoreset_json: dict[str, Any]) -> UniProtRef | None:
    """Extract UniProt reference from scoreset metadata if available.

    :param scoreset_json: parsed JSON from scoresets API
    :return: UniProt ID if available
    """
    ext_ids = scoreset_json["targetGenes"][0].get("externalIdentifiers")
    if not ext_ids:
        return None
    for ext_id in ext_ids:
        if ext_id.get("identifier", {}).get("dbName") == "UniProt":
            return UniProtRef(
                id=f"uniprot:{ext_id['identifier']['identifier']}",
                offset=ext_id["offset"],
            )
    return None


def get_raw_scoreset_metadata(
    scoreset_urn: str, dcd_mapping_dir: Path | None = None
) -> dict:
    """Get raw (original JSON) metadata for scoreset.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :param dcd_mapping_dir: optionally declare location to save metadata to.
    :return: Complete JSON response for object
    :raise ResourceAcquisitionError: if HTTP request fails
    """
    if not dcd_mapping_dir:
        dcd_mapping_dir = LOCAL_STORE_PATH

    metadata_file = dcd_mapping_dir / f"{scoreset_urn}_metadata.json"
    if not metadata_file.exists():
        url = f"{MAVEDB_BASE_URL}/api/v1/score-sets/{scoreset_urn}"
        r = requests.get(url, timeout=30, headers=authentication_header())
        try:
            r.raise_for_status()
        except requests.HTTPError as e:
            msg = f"Received HTTPError from {url} for scoreset {scoreset_urn}"
            _logger.error(msg)
            raise ResourceAcquisitionError(msg) from e
        metadata = r.json()
        with metadata_file.open("w") as f:
            json.dump(metadata, f)
    else:
        with metadata_file.open() as f:
            metadata = json.load(f)
    return metadata


def get_scoreset_metadata(
    scoreset_urn: str, dcd_mapping_dir: Path | None = None
) -> ScoresetMetadata:
    """Acquire metadata for scoreset.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :param dcd_mapping_dir: optionally declare location to save metadata to.
    :return: Object containing salient metadata
    :raise ResourceAcquisitionError: if unable to acquire metadata
    """
    metadata = get_raw_scoreset_metadata(scoreset_urn, dcd_mapping_dir)
    target_genes = {}
    multi_target = len(metadata["targetGenes"]) > 1

    for gene in metadata["targetGenes"]:
        target_sequence_gene = gene.get("targetSequence")
        if target_sequence_gene is None:
            msg = f"No target sequence available for {scoreset_urn}. Accession-based score sets are not currently supported."
            raise ScoresetNotSupportedError(msg)
        if not _metadata_response_is_human(metadata):
            msg = f"Experiment for {scoreset_urn} contains no human targets"
            raise ScoresetNotSupportedError(msg)
        try:
            target_sequence_label = gene["targetSequence"].get("label")
            if target_sequence_label is None:
                # if there are not multiple targets, label is not required by mavedb,
                # so use target gene name as the label.
                if not multi_target:
                    target_sequence_label = gene["name"]
                else:
                    msg = f"No target label provided for target in multi-target score set {scoreset_urn}."
                    raise DataLookupError(msg)
            target_genes[target_sequence_label] = TargetGene(
                target_gene_name=gene["name"],
                target_gene_category=gene["category"],
                target_sequence=gene["targetSequence"]["sequence"],
                target_sequence_type=gene["targetSequence"]["sequenceType"],
                target_sequence_label=target_sequence_label,
                target_uniprot_ref=_get_uniprot_ref(metadata),
            )
        except (KeyError, ValidationError) as e:
            msg = f"Unable to extract metadata from API response for scoreset {scoreset_urn}: {e}"
            _logger.error(msg)
            raise ScoresetNotSupportedError(msg) from e

    return ScoresetMetadata(urn=scoreset_urn, target_genes=target_genes)


def _load_scoreset_records(
    path: Path, metadata: ScoresetMetadata
) -> dict[str, list[ScoreRow]]:
    """Load scoreset records from CSV file.
    Organize scoreset records by reference sequence prefix / target gene label.
    If no reference sequence prefix is provided, the score set should only have one
    target, so use the one target's label.

    This method is intentionally identified as "private", but is refactored out for
    use during testing.
    """
    scores_data: dict[str, list[ScoreRow]] = {}
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["score"] == "NA":
                row["score"] = None
            else:
                row["score"] = row["score"]
            if row["hgvs_nt"] != "NA":
                # TODO check assumption of no colon in hgvs unless reference sequence identifier present
                prefix = row["hgvs_nt"].split(":")[0] if ":" in row["hgvs_nt"] else None
            elif row["hgvs_pro"] != "NA":
                # TODO check assumption of no colon in hgvs unless reference sequence identifier present
                prefix = (
                    row["hgvs_pro"].split(":")[0] if ":" in row["hgvs_pro"] else None
                )
            else:
                # Should we quit the whole mapping job if this comes up, or just skip this row and only quit if none contain hgvs_nt or hgvs_pro?
                msg = f"Each score row in {metadata.urn} must contain hgvs_nt or hgvs_pro variant description "
                raise ScoresetNotSupportedError(msg)
            # If no reference sequence prefix is provided, the score set should only have one
            # target, so use the one target's label.
            if prefix is None:
                if len(metadata.target_genes) == 1:
                    # don't know name of key, so loop through dictionary
                    # and use the key of the single-entry metadata.target_genes dict as the key for this new dict
                    for target_gene in metadata.target_genes:
                        prefix = target_gene
                else:
                    msg = f"Multi-target score set {metadata.urn} does not contain target sequence labels in variant HGVS strings."
                    raise ScoresetNotSupportedError(msg)
            if prefix in scores_data:
                scores_data[prefix].append(ScoreRow(**row))
            else:
                scores_data[prefix] = [ScoreRow(**row)]
    return scores_data


def _get_experiment_53_scores(outfile: Path, silent: bool) -> None:
    """Scores for `urn:mavedb:00000053-a-1` can be hard to acquire from the server
    on account of their considerable size. This method uses a basic workaround to fetch
    a copy hosted on a GitHub issue until we have a final resolution.
    """
    url = "https://github.com/VariantEffect/mavedb-api/files/13746791/00000053-a-1.zip"
    with tempfile.NamedTemporaryFile() as temp_file:
        path = Path(temp_file.name)
        http_download(url, path, silent)
        with zipfile.ZipFile(path, "r") as zip_ref:
            with zip_ref.open("00000053-a-1/00000053-a-1.scores.csv") as file:
                contents = file.read()
            with outfile.open("wb") as f:
                f.write(contents)


def get_scoreset_records(
    metadata: ScoresetMetadata, silent: bool = True, dcd_mapping_dir: Path | None = None
) -> dict[str, list[ScoreRow]]:
    """Get scoreset records.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param urn: URN for scoreset
    :param silent: if true, suppress console output
    :param dcd_mapping_dir: optionally declare save location for records
    :return: Array of individual ScoreRow objects, containing information like protein
        changes and scores
    :raise ResourceAcquisitionError: if unable to fetch file
    """
    if not dcd_mapping_dir:
        dcd_mapping_dir = LOCAL_STORE_PATH
    scores_csv = dcd_mapping_dir / f"{metadata.urn}_scores.csv"
    # TODO use smarter/more flexible caching methods
    if not scores_csv.exists():
        if metadata.urn == "urn:mavedb:00000053-a-1":
            _get_experiment_53_scores(scores_csv, silent)
        else:
            url = f"{MAVEDB_BASE_URL}/api/v1/score-sets/{metadata.urn}/scores"
            try:
                http_download(url, scores_csv, silent)
            except requests.HTTPError as e:
                msg = f"HTTPError when fetching scores CSV from {url}"
                _logger.error(msg)
                raise ResourceAcquisitionError(msg) from e

    return _load_scoreset_records(scores_csv, metadata)


def with_mavedb_score_set(fn: Callable) -> Callable:
    @wraps(fn)
    async def wrapper(*args, **kwargs) -> ScoresetMapping:  # noqa: ANN002
        urn = args[0] if args else kwargs["urn"]
        silent = kwargs.get("silent", False)

        with tempfile.TemporaryDirectory(
            prefix=f"{LOCAL_STORE_PATH.as_posix()}/"
        ) as temp_dir:
            # Set up metadata and scores for the current run. Now they will be accessible by these functions
            # without the need to download the data again.
            temp_dir_as_path = Path(temp_dir)
            try:
                metadata = get_scoreset_metadata(urn, temp_dir_as_path)
                get_scoreset_records(metadata, silent, temp_dir_as_path)
            except ScoresetNotSupportedError as e:
                return ScoresetMapping(
                    metadata=None,
                    error_message=str(e).strip("'"),
                )
            except ResourceAcquisitionError as e:
                msg = f"Unable to acquire resource from MaveDB: {e}"
                raise HTTPException(status_code=500, detail=msg) from e

            # Pass the storage path of the temp directory to the wrapped function as a kwarg.
            kwargs["store_path"] = temp_dir_as_path
            v: ScoresetMapping = await fn(*args, **kwargs)

            return v

    return wrapper
