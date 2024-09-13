""""Provide mapping router"""
from cool_seq_tool.schemas import AnnotationLayer
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from requests import HTTPError

from dcd_mapping.align import AlignmentError, BlatNotFoundError, align
from dcd_mapping.annotate import (
    _get_computed_reference_sequence,
    _get_mapped_reference_sequence,
    _set_scoreset_layer,
    annotate,
)
from dcd_mapping.lookup import DataLookupError
from dcd_mapping.mavedb_data import (
    ScoresetNotSupportedError,
    get_raw_scoreset_metadata,
    get_scoreset_metadata,
    get_scoreset_records,
)
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import ScoreAnnotation, ScoresetMapping, VrsVersion
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.vrs_map import VrsMapError, vrs_map

router = APIRouter(
    prefix="/api/v1", tags=["mappings"], responses={404: {"description": "Not found"}}
)


@router.post(path="/map/{urn}", status_code=200, response_model=ScoresetMapping)
async def map_scoreset(urn: str) -> ScoresetMapping:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param output_path: optional path to save output at
    :param vrs_version: version of VRS objects to output (1.3 or 2)
    :param silent: if True, suppress console information output
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, True)
    except ScoresetNotSupportedError as e:
        return ScoresetMapping(
            metadata=None,
            error_message=str(e).strip("'"),
        )
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        raise HTTPException(status_code=500, detail=msg) from e

    try:
        alignment_result = align(metadata, True)
    except BlatNotFoundError as e:
        msg = "BLAT command appears missing. Ensure it is available on the $PATH or use the environment variable BLAT_BIN_PATH to point to it. See instructions in the README prerequisites section for more."
        raise HTTPException(status_code=500, detail=msg) from e
    except ResourceAcquisitionError as e:
        msg = f"BLAT resource could not be acquired: {e}"
        raise HTTPException(status_code=500, detail=msg) from e
    except AlignmentError as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    try:
        transcript = await select_transcript(metadata, records, alignment_result)
    except (TxSelectError, KeyError, ValueError) as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    except HTTPError as e:
        msg = f"HTTP error occurred during transcript selection: {e}"
        raise HTTPException(status_code=500, detail=msg) from e
    except DataLookupError as e:
        msg = f"Data lookup error occurred during transcript selection: {e}"
        raise HTTPException(status_code=500, detail=msg) from e

    try:
        vrs_results = vrs_map(metadata, alignment_result, records, transcript, True)
    except VrsMapError as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    if vrs_results is None:
        return ScoresetMapping(
            metadata=metadata,
            error_message="No variant mappings available for this score set",
        )

    try:
        vrs_results = annotate(vrs_results, transcript, metadata, VrsVersion.V_2)
    except Exception as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    if vrs_results is None:
        return ScoresetMapping(
            metadata=metadata,
            error_message="No annotated variant mappings available for this score set",
        )

    try:
        raw_metadata = get_raw_scoreset_metadata(urn)
        preferred_layers = {
            _set_scoreset_layer(urn, vrs_results),
        }

        reference_sequences = {
            layer: {
                "computed_reference_sequence": None,
                "mapped_reference_sequence": None,
            }
            for layer in AnnotationLayer
        }

        for layer in preferred_layers:
            reference_sequences[layer][
                "computed_reference_sequence"
            ] = _get_computed_reference_sequence(urn, layer, transcript)
            reference_sequences[layer][
                "mapped_reference_sequence"
            ] = _get_mapped_reference_sequence(layer, transcript, alignment_result)

        mapped_scores: list[ScoreAnnotation] = []
        for m in vrs_results:
            if m.annotation_layer in preferred_layers:
                # drop annotation layer from mapping object
                mapped_scores.append(ScoreAnnotation(**m.model_dump()))
    except Exception as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    return JSONResponse(
        content=ScoresetMapping(
            metadata=raw_metadata,
            computed_protein_reference_sequence=reference_sequences[
                AnnotationLayer.PROTEIN
            ]["computed_reference_sequence"],
            mapped_protein_reference_sequence=reference_sequences[
                AnnotationLayer.PROTEIN
            ]["mapped_reference_sequence"],
            computed_genomic_reference_sequence=reference_sequences[
                AnnotationLayer.GENOMIC
            ]["computed_reference_sequence"],
            mapped_genomic_reference_sequence=reference_sequences[
                AnnotationLayer.GENOMIC
            ]["mapped_reference_sequence"],
            mapped_scores=mapped_scores,
        ).model_dump(exclude_none=True)
    )
