"""Provide core MaveDB mapping methods."""
import json
import logging
from typing import List

import click

from dcd_mapping.align import AlignmentError, align
from dcd_mapping.resources import (
    LOCAL_STORE_PATH,
    ResourceAcquisitionError,
    get_scoreset_metadata,
    get_scoreset_records,
)
from dcd_mapping.schemas import ScoreRow, ScoresetMetadata, VrsMappingResult
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


def _save_results(
    metadata: ScoresetMetadata, mapping_results: VrsMappingResult
) -> None:
    """Save results to file.

    Todo:
    ----
    * Embed in original metadata JSON
    * Option to save as VRS 1.x

    :param metadata: scoreset metadata
    :param mapping results: mapped objects
    """
    outfile = LOCAL_STORE_PATH / f"{metadata.urn}_mapping_results.json"
    with open(outfile, "w") as f:
        json.dump(mapping_results.model_dump_json(indent=2), f)


async def map_scoreset(
    metadata: ScoresetMetadata, records: List[ScoreRow], silent: bool = True
) -> None:
    """Given information about a MAVE experiment, map to VRS.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param silent:
    :return: something (TODO)
    """
    try:
        alignment_result = align(metadata, silent, True)
    except AlignmentError:
        _logger.error(f"Alignment failed for scoreset {metadata.urn}")
        return None

    try:
        transcript = await select_transcript(
            metadata, records, alignment_result, silent
        )
    except TxSelectError:
        _logger.error(f"Transcript selection failed for scoreset {metadata.urn}")
        return None

    try:
        vrs_results = vrs_map(metadata, alignment_result, transcript, records)
    except VrsMapError:
        _logger.error(f"VRS mapping failed for scoreset {metadata.urn}")
        return None

    if vrs_results:
        _save_results(metadata, vrs_results)


async def map_scoreset_urn(urn: str, silent: bool = True) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :return: something (TODO)
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        return None
    mapped = await map_scoreset(metadata, records, silent)
    return mapped  # TODO not sure what this will be