"""Provide core MaveDB mapping methods."""

import logging
import os
import subprocess
from pathlib import Path

import click
from requests import HTTPError

from dcd_mapping.align import AlignmentError, BlatNotFoundError, build_alignment_result
from dcd_mapping.annotate import (
    annotate,
    save_mapped_output_json,
    write_scoreset_mapping_to_json,
)
from dcd_mapping.lookup import (
    DataLookupError,
    check_gene_normalizer,
    check_seqrepo,
    check_uta,
)
from dcd_mapping.mavedb_data import (
    ScoresetNotSupportedError,
    get_scoreset_metadata,
    get_scoreset_records,
    with_mavedb_score_set,
)
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMapping,
    ScoresetMetadata,
    VrsVersion,
)
from dcd_mapping.transcripts import select_transcripts
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


def _emit_info(msg: str, silent: bool, log_level: int = logging.INFO) -> None:
    if not silent:
        click.echo(msg)
    if log_level == logging.INFO:
        _logger.info(msg)
    elif log_level == logging.ERROR:
        _logger.error(msg)
    else:
        msg = f"Unexpected log level requested: {log_level}"
        raise ValueError(msg)


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
    records: dict[str, list[ScoreRow]],
    output_path: Path | None = None,
    vrs_version: VrsVersion = VrsVersion.V_2,
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

    _emit_info(f"Performing alignment for {metadata.urn}...", silent)
    try:
        # dictionary where keys are target gene labels or accession ids, and values are alignment result objects
        alignment_results = build_alignment_result(metadata, silent)
    except BlatNotFoundError as e:
        msg = "BLAT command appears missing. Ensure it is available on the $PATH or use the environment variable BLAT_BIN_PATH to point to it. See instructions in the README prerequisites section for more."
        _emit_info(msg, silent, logging.ERROR)
        raise e
    except ResourceAcquisitionError as e:
        _emit_info(f"BLAT resource could not be acquired: {e}", silent, logging.ERROR)
        raise e
    except AlignmentError as e:
        _emit_info(
            f"Alignment failed for scoreset  {metadata.urn} {e}", silent, logging.ERROR
        )
        final_output = write_scoreset_mapping_to_json(
            metadata.urn,
            ScoresetMapping(metadata=metadata, error_message=str(e).strip("'")),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    except ScoresetNotSupportedError as e:
        _emit_info(f"Score set not supported: {e}", silent, logging.ERROR)
        final_output = write_scoreset_mapping_to_json(
            metadata.urn,
            ScoresetMapping(metadata=metadata, error_message=str(e).strip("'")),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    _emit_info("Alignment complete.", silent)

    _emit_info("Selecting reference sequence...", silent)
    try:
        transcripts = await select_transcripts(metadata, records, alignment_results)
    # NOTE: transcript selection errors are handled in select_transcripts,
    # and they do not cause the entire mapping process to exit; instead, an error will be reported
    # on the target level and on the variant level for variants relative to that target
    # HTTPErrors and DataLookupErrors cause the mapping process to exit because these indicate
    # underlying issues with data providers.
    except HTTPError as e:
        _emit_info(
            f"HTTP error occurred during transcript selection: {e}",
            silent,
            logging.ERROR,
        )
        raise e
    except DataLookupError as e:
        _emit_info(
            f"Data lookup error occurred during transcript selection: {e}",
            silent,
            logging.ERROR,
        )
        raise e
    _emit_info("Reference selection complete.", silent)

    _emit_info("Mapping to VRS...", silent)
    vrs_results = {}
    try:
        for target_gene in metadata.target_genes:
            vrs_results[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=alignment_results[target_gene],
                records=records[target_gene],
                transcript=transcripts[target_gene],
                silent=silent,
            )
    except VrsMapError as e:
        _emit_info(
            f"VRS mapping failed for scoreset {metadata.urn}", silent, logging.ERROR
        )
        final_output = write_scoreset_mapping_to_json(
            metadata.urn,
            ScoresetMapping(metadata=metadata, error_message=str(e).strip("'")),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    if not vrs_results or all(
        mapping_result is None for mapping_result in vrs_results.values()
    ):
        _emit_info(f"No mapping available for {metadata.urn}", silent, logging.ERROR)
        final_output = write_scoreset_mapping_to_json(
            metadata.urn,
            ScoresetMapping(
                metadata=metadata,
                error_message="No variant mappings available for this score set",
            ),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    _emit_info("VRS mapping complete.", silent)

    _emit_info("Annotating metadata and saving to file...", silent)
    # annotate each target's variants separately, since annotation is target specific (e.g. obtaining reference sequence)
    annotated_vrs_results = {}
    for target_gene in vrs_results:
        try:
            annotated_vrs_results[target_gene] = annotate(
                vrs_results[target_gene],
                transcripts[target_gene],
                metadata.target_genes[target_gene],
                metadata.urn,
                vrs_version,
            )
        except Exception as e:
            _emit_info(
                f"VRS annotation failed for scoreset {metadata.urn}",
                silent,
                logging.ERROR,
            )
            final_output = write_scoreset_mapping_to_json(
                metadata.urn,
                ScoresetMapping(metadata=metadata, error_message=str(e).strip("'")),
                output_path,
            )
            _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
            return
    if not annotated_vrs_results or all(
        mapping_result is None for mapping_result in annotated_vrs_results.values()
    ):
        _emit_info(f"No annotation available for {metadata.urn}", silent, logging.ERROR)
        final_output = write_scoreset_mapping_to_json(
            metadata.urn,
            ScoresetMapping(
                metadata=metadata,
                error_message="No annotated variant mappings available for this score set",
            ),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    try:
        final_output = save_mapped_output_json(
            metadata,
            annotated_vrs_results,
            alignment_results,
            transcripts,
            prefer_genomic,
            output_path,
        )
    except Exception as e:
        _emit_info(
            f"Error in creating or saving final score set mapping for {metadata.urn} {e}",
            silent,
            logging.ERROR,
        )
        final_output = write_scoreset_mapping_to_json(
            metadata.urn,
            ScoresetMapping(metadata=metadata, error_message=str(e).strip("'")),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    _emit_info(f"Annotated scores saved to: {final_output}.", silent)


@with_mavedb_score_set
async def map_scoreset_urn(
    urn: str,
    output_path: Path | None = None,
    vrs_version: VrsVersion = VrsVersion.V_2,
    prefer_genomic: bool = False,
    silent: bool = True,
    store_path: Path | None = None,
) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param output_path: optional path to save output at
    :param vrs_version: version of VRS objects to output (1.3 or 2)
    :param silent: if True, suppress console information output
    """
    try:
        metadata = get_scoreset_metadata(urn, store_path)
        records = get_scoreset_records(metadata, silent, store_path)
    except ScoresetNotSupportedError as e:
        _emit_info(f"Score set not supported: {e}", silent, logging.ERROR)
        final_output = write_scoreset_mapping_to_json(
            urn,
            ScoresetMapping(
                metadata=None,
                error_message=str(e).strip("'"),
            ),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        raise e

    if not records:
        _emit_info("Score set contains no variants to map", silent, logging.ERROR)
        final_output = write_scoreset_mapping_to_json(
            urn,
            ScoresetMapping(
                metadata=metadata,
                error_message="Score set contains no variants to map",
            ),
            output_path,
        )
        _emit_info(f"Score set mapping output saved to: {final_output}.", silent)
        return

    await map_scoreset(
        metadata, records, output_path, vrs_version, prefer_genomic, silent
    )
