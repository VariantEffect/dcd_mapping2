"""Provide core MaveDB mapping methods."""

import logging
import os
import subprocess
from pathlib import Path

import click

from dcd_mapping.align import build_alignment_result
from dcd_mapping.annotate import annotate, save_mapped_output_json
from dcd_mapping.log_utilities import _emit_info
from dcd_mapping.lookup import check_gene_normalizer, check_seqrepo, check_uta
from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMetadata,
)
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


async def _check_data_prereqs(silent: bool) -> None:
    """Non-exhaustive check that data prereqs are properly configured and available."""
    _emit_info("Checking data prereqs....", silent)
    success = True
    try:
        await check_uta()
    except Exception:
        success = False
        _emit_info(
            "* UTA appears to be unavailable. Check the logs for more information. For troubleshooting, we recommend checking the UTA readme (https://github.com/biocommons/uta?tab=readme-ov-file#installing-uta-locally) and the Cool-Seq-Tool installation instructions (https://coolseqtool.readthedocs.io/0.4.0-dev3/install.html#set-up-uta). Remember that the UTA connection is configurable via a libpq URI provided under the environment variable UTA_DB_URL (see Cool-Seq-Tool docs: https://coolseqtool.readthedocs.io/0.4.0-dev3/usage.html#environment-configuration) -- otherwise, by default it attempts a connection to `postgresql://uta_admin:uta@localhost:5433/uta/uta_20210129b`.",
            silent,
            logging.ERROR,
        )
    try:
        check_seqrepo()
    except Exception:
        success = False
        _emit_info(
            "* SeqRepo appears inaccessible or unusable. Check the logs for more information. Ensure that a local SeqRepo snapshot has been downloaded (it should've taken a while -- see https://github.com/biocommons/biocommons.seqrepo?tab=readme-ov-file#requirements), that it's located either at `/usr/local/share/seqrepo/latest` or at the location designated by the `SEQREPO_ROOT_DIR` environment variable, and that it's writeable (see https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/store.rst).",
            silent,
            logging.ERROR,
        )
    try:
        check_gene_normalizer()
    except Exception:
        success = False
        _emit_info(
            "* Gene Normalizer appears to be unavailable. Check the logs for more information. Note that a data snapshot needs to be acquired, or the data update routine must be routine (this should've taken at least a few seconds, if not several minutes). For troubleshooting, review the Gene Normalizer installation instructions and documentation: https://gene-normalizer.readthedocs.io/0.3.0-dev1/install.html",
            silent,
            logging.ERROR,
        )
    try:
        configured_blat_bin = os.environ.get("BLAT_BIN_PATH")
        if configured_blat_bin:
            result = subprocess.run(  # noqa: ASYNC101
                configured_blat_bin,  # noqa: S603
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if result.returncode == 127:
                success = False
                _emit_info(
                    f"* BLAT binary at location pointed to by BLAT_BIN_PATH env var, {configured_blat_bin}, appears to be missing. Please check that a BLAT executable is at that location.",
                    silent,
                    logging.ERROR,
                )
            elif result.returncode != 0 and result.returncode != 255:
                success = False
                _emit_info(
                    f"* BLAT binary at location pointed to by BLAT_BIN_PATH env var, {configured_blat_bin}, doesn't appear to be properly executable. Please double-check that the executable is at the proper location and has correct permissions.",
                    silent,
                    logging.ERROR,
                )
        else:
            result = subprocess.run(  # noqa: ASYNC101
                "blat",  # noqa: S603 S607
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if result.returncode == 127:
                success = False
                _emit_info(
                    "* Unable to run BLAT. The BLAT binary must be acquired separately by the user and either accessible in the $PATH or at the location pointed to by the BLAT_BIN_PATH environment variable.",
                    silent,
                    logging.ERROR,
                )
            elif result.returncode != 0 and result.returncode != 255:
                success = False
                _emit_info(
                    "* BLAT binary appear to be properly executable. Please double-check that the executable is at the proper location and has correct permissions.",
                    silent,
                    logging.ERROR,
                )
    except Exception:
        success = False
        _emit_info(
            "Encountered unexpected error while testing availability of BLAT. Check logs for more information. See README for more information about BLAT setup.",
            silent,
            logging.ERROR,
        )
    if not success:
        raise LookupError
    _emit_info("Data prereqs checks pass.", silent)


async def map_scoreset(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    output_path: Path | None = None,
    include_vrs_2: bool = False,
    prefer_genomic: bool = False,
    silent: bool = True,
) -> None:
    """Given information about a MAVE experiment, map to VRS and save output as JSON.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param output_path: optional path to save output at
    :param include_vrs_2: if true, include VRS 2.0 mappings in output JSON
    :param silent: if True, suppress console information output
    """
    await _check_data_prereqs(silent)

    alignment_result = build_alignment_result(metadata, silent)

    # transcript =

    # accession based
    # pre mapped - take hgvs string literally from hgvs_pro and hgvs_nt columns
    # protein accession (NP, ENSP):
    # if they gave us np or ensp, they at least have the hgvs_pro filled out
    # create p. hgvs string (this is just what they gave us) and vrs object
    # need to pass the np accession id (which is actually in the metadata)
    # if hgvs_nt column is available, then we can also create genomic object
    # TODO do noncoding accessions that are not contigs exist? NR
    # transcript accession (NM, ENST):
    # assume they  have the ngvs_nt filled out?
    # no need to create p. hgvs string or p. vrs object
    # post mapped: NC (from alignment output) + adjusted variant coordinates (from hgvs_nt column)

    # if accession is NP - basically do nothing unless they filled out the hgvs_nt column
    # but if they give us an NP and then fill out the hgvs_nt column, we would need to find the mapping
    # to an hg38 contig, but cdot api doesn't seem to support np/ensp accessions, only transcripts
    # but will we ever get an np accession?

    # np - skip for now
    # nm - map (already done)
    # assuming that if nm/enst provided, hgvs_nt column is filled out
    # pre map - hgvs string is exactly what was provided by user (NM_xxx.c.xxx)
    # post map - same flow as for target sequence. NC_xxx.g.<adjusted_variant>
    # nc
    # pre map and post map are both exactly the hgvs string provided by user
    # always make sure

    _emit_info("Selecting reference sequence...", silent)
    try:
        transcript = await select_transcript(metadata, records, alignment_result)
    except TxSelectError as e:
        _emit_info(
            f"Transcript selection failed for scoreset {metadata.urn}",
            silent,
            logging.ERROR,
        )
        raise e
    _emit_info("Reference selection complete.", silent)

    _emit_info("Mapping to VRS...", silent)
    try:
        vrs_results = vrs_map(metadata, alignment_result, records, transcript, silent)
    except VrsMapError as e:
        _emit_info(
            f"VRS mapping failed for scoreset {metadata.urn}", silent, logging.ERROR
        )
        raise e
    if vrs_results is None:
        _emit_info(f"No mapping available for {metadata.urn}", silent)
        return
    _emit_info("VRS mapping complete.", silent)

    _emit_info("Annotating metadata and saving to file...", silent)
    vrs_results = annotate(vrs_results, transcript, metadata)
    final_output = save_mapped_output_json(
        metadata.urn,
        vrs_results,
        alignment_result,
        transcript,
        include_vrs_2,
        prefer_genomic,
        output_path,
    )
    _emit_info(f"Annotated scores saved to: {final_output}.", silent)


async def map_scoreset_urn(
    urn: str,
    output_path: Path | None = None,
    include_vrs_2: bool = False,
    prefer_genomic: bool = False,
    silent: bool = True,
) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param output_path: optional path to save output at
    :param include_vrs_2: if true, include VRS 2.0 mappings in output JSON
    :param silent: if True, suppress console information output
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        raise e
    await map_scoreset(
        metadata, records, output_path, include_vrs_2, prefer_genomic, silent
    )
