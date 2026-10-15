"""Provide mapping router"""

import logging
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from httpx import HTTPStatusError

from dcd_mapping.align import build_alignment_result
from dcd_mapping.annotate import (
    annotate,
    build_scoreset_mapping,
    compute_target_gene_info,
)
from dcd_mapping.exceptions import (
    AlignmentError,
    BlatNotFoundError,
    DataLookupError,
    MissingSequenceIdError,
    ResourceAcquisitionError,
    ScoresetNotSupportedError,
    UnsupportedReferenceSequenceNameSpaceError,
    UnsupportedReferenceSequencePrefixError,
    VrsMapError,
)
from dcd_mapping.mavedb_data import (
    get_raw_scoreset_metadata,
    get_scoreset_metadata,
    get_scoreset_records,
    patch_target_sequence_type,
    with_mavedb_score_set,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMapping,
    VrsVersion,
)
from dcd_mapping.transcripts import select_transcripts
from dcd_mapping.vrs_map import vrs_map

router = APIRouter(
    prefix="/api/v1", tags=["mappings"], responses={404: {"description": "Not found"}}
)

_logger = logging.getLogger(__name__)


@router.post(path="/map/{urn}", status_code=200, response_model=ScoresetMapping)
@with_mavedb_score_set
async def map_scoreset(urn: str, store_path: Path | None = None) -> JSONResponse:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param store_path: optional path to save output at
    """
    try:
        raw_metadata = get_raw_scoreset_metadata(urn, store_path)
        metadata = get_scoreset_metadata(urn, store_path)
        records = get_scoreset_records(metadata, True, store_path)
        metadata = patch_target_sequence_type(metadata, records, force=False)
    except ScoresetNotSupportedError as e:
        _logger.error("Scoreset not supported for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=None,
                error_message=str(e).strip("'"),
            ).model_dump(exclude_none=True)
        )
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.error(msg)
        raise HTTPException(status_code=500, detail=msg) from e

    if not records:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata,
                error_message="Score set contains no variants to map",
            ).model_dump(exclude_none=True)
        )
    total_score_records = sum(len(v) for v in records.values())

    try:
        alignment_results = build_alignment_result(metadata, True)
    except BlatNotFoundError as e:
        msg = "BLAT command appears missing. Ensure it is available on the $PATH or use the environment variable BLAT_BIN_PATH to point to it. See instructions in the README prerequisites section for more."
        _logger.error("BLAT not found for %s: %s", urn, e)
        raise HTTPException(status_code=500, detail=msg) from e
    except ResourceAcquisitionError as e:
        msg = f"BLAT resource could not be acquired: {e}"
        _logger.error(msg)
        raise HTTPException(status_code=500, detail=msg) from e
    except AlignmentError as e:
        _logger.error("Alignment error for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    except ScoresetNotSupportedError as e:
        _logger.error("Scoreset not supported during alignment for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    try:
        transcripts = await select_transcripts(metadata, records, alignment_results)
    # NOTE: transcript selection errors are handled in select_transcripts,
    # and they do not cause the entire mapping process to exit; instead, an error will be reported
    # on the target level and on the variant level for variants relative to that target
    # HTTPErrors and DataLookupErrors cause the mapping process to exit because these indicate
    # underlying issues with data providers.
    except HTTPStatusError as e:
        msg = f"HTTP error occurred during transcript selection: {e}"
        _logger.error(msg)
        raise HTTPException(status_code=500, detail=msg) from e
    except DataLookupError as e:
        msg = f"Data lookup error occurred during transcript selection: {e}"
        _logger.error(msg)
        raise HTTPException(status_code=500, detail=msg) from e

    vrs_results = {}
    protein_align_results: dict[str, AlignmentResult | None] = {}
    try:
        for target_gene in metadata.target_genes:
            target_records = records.get(target_gene)

            # e.g. base-editor score sets that declare separate protein and
            # cDNA accession targets for the same variant: every row's hgvs_nt
            # prefix groups under the cDNA target, so the protein target has no
            # record group of its own and contributes nothing independently.
            if target_records is None:
                _logger.info(
                    "No score records reference target %s directly; skipping standalone VRS mapping for this target.",
                    target_gene,
                )
                continue

            vrs_map_result = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=alignment_results[target_gene],
                records=target_records,
                transcript=transcripts[target_gene],
                silent=True,
            )
            vrs_results[target_gene] = vrs_map_result.mappings
            protein_align_results[target_gene] = vrs_map_result.protein_align_result
    except (
        UnsupportedReferenceSequenceNameSpaceError,
        VrsMapError,
        UnsupportedReferenceSequencePrefixError,
        MissingSequenceIdError,
    ) as e:
        _logger.error("VRS mapping error for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    nonetype_vrs_results = [
        result is None
        for target_gene in vrs_results
        for result in vrs_results[target_gene]
    ]

    if not vrs_results or all(nonetype_vrs_results):
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata,
                error_message="No variant mappings available for this score set",
            ).model_dump(exclude_none=True)
        )
    if any(nonetype_vrs_results):
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata,
                error_message="Some variants generated vrs results, but not all. If any variants were mapped, all should have been.",
            ).model_dump(exclude_none=True)
        )

    annotated_vrs_results = {}
    try:
        for target_gene in vrs_results:
            annotated_vrs_results[target_gene] = annotate(
                vrs_results[target_gene],
                transcripts[target_gene],
                metadata.target_genes[target_gene],
                metadata.urn,
                VrsVersion.V_2,
            )
    except Exception as e:
        _logger.error("Unexpected error during annotation for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    nonetype_annotated_vrs_results = [
        result is None
        for target_gene in annotated_vrs_results
        for result in annotated_vrs_results[target_gene]
    ]

    if not annotated_vrs_results or all(nonetype_annotated_vrs_results):
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata,
                error_message="No annotated variant mappings available for this score set",
            ).model_dump(exclude_none=True)
        )
    if any(nonetype_annotated_vrs_results):
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata,
                error_message="Some variants generated annotated vrs results, but not all. If any variants were annotated, all should have been.",
            ).model_dump(exclude_none=True)
        )

    try:
        gene_info = {}
        for target_gene in annotated_vrs_results:
            gene_info[target_gene] = await compute_target_gene_info(
                target_key=target_gene,
                transcripts=transcripts,
                alignment_results=alignment_results,
                metadata=metadata,
                mapped_scores=annotated_vrs_results[target_gene],
            )
        output = build_scoreset_mapping(
            metadata=metadata,
            mappings=annotated_vrs_results,
            align_results=alignment_results,
            tx_output=transcripts,
            gene_info=gene_info,
            preferred_layer_only=True,
            vrs_version=VrsVersion.V_2,
            raw_metadata=raw_metadata,
            protein_align_results=protein_align_results,
        )
    except Exception as e:
        _logger.error("Unexpected error during result assembly for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    # With preferred_layer_only=True, build_scoreset_mapping emits exactly one
    # mapped_score per input variant: preferred-layer successes/failures plus
    # completely-failed variants (annotation_layer=None) re-attributed to the
    # preferred layer.  The count must equal total_score_records.
    if len(output.mapped_scores or []) != total_score_records:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=raw_metadata,
                error_message=f"Mismatch between number of mapped scores ({len(output.mapped_scores or [])}) and total score records ({total_score_records}). This is unexpected and indicates an issue with the mapping process.",
            ).model_dump(exclude_none=True)
        )

    return JSONResponse(content=output.model_dump(exclude_none=True))
