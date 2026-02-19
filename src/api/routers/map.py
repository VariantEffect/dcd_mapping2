""""Provide mapping router"""
import logging
from pathlib import Path

from cool_seq_tool.schemas import AnnotationLayer
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from requests import HTTPError

from dcd_mapping.align import build_alignment_result
from dcd_mapping.annotate import (
    _get_computed_reference_sequence,
    _get_mapped_reference_sequence,
    _set_scoreset_layer,
    annotate,
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
    ScoreAnnotation,
    ScoresetMapping,
    TargetAnnotation,
    TargetType,
    TxSelectResult,
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
                metadata=metadata,
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
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    except ScoresetNotSupportedError as e:
        _logger.error("Scoreset not supported during alignment for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    try:
        transcripts = await select_transcripts(metadata, records, alignment_results)
    # NOTE: transcript selection errors are handled in select_transcripts,
    # and they do not cause the entire mapping process to exit; instead, an error will be reported
    # on the target level and on the variant level for variants relative to that target
    # HTTPErrors and DataLookupErrors cause the mapping process to exit because these indicate
    # underlying issues with data providers.
    except HTTPError as e:
        msg = f"HTTP error occurred during transcript selection: {e}"
        _logger.error(msg)
        raise HTTPException(status_code=500, detail=msg) from e
    except DataLookupError as e:
        msg = f"Data lookup error occurred during transcript selection: {e}"
        _logger.error(msg)
        raise HTTPException(status_code=500, detail=msg) from e

    vrs_results = {}
    try:
        for target_gene in metadata.target_genes:
            vrs_results[target_gene] = vrs_map(
                metadata=metadata.target_genes[target_gene],
                align_result=alignment_results[target_gene],
                records=records[target_gene],
                transcript=transcripts[target_gene],
                silent=True,
            )
    except (
        UnsupportedReferenceSequenceNameSpaceError,
        VrsMapError,
        UnsupportedReferenceSequencePrefixError,
        MissingSequenceIdError,
    ) as e:
        _logger.error("VRS mapping error for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
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
                metadata=metadata,
                error_message="No variant mappings available for this score set",
            ).model_dump(exclude_none=True)
        )
    if any(nonetype_vrs_results):
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata,
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
                metadata=metadata, error_message=str(e).strip("'")
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
                metadata=metadata,
                error_message="No annotated variant mappings available for this score set",
            ).model_dump(exclude_none=True)
        )
    if any(nonetype_annotated_vrs_results):
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata,
                error_message="Some variants generated annotated vrs results, but not all. If any variants were annotated, all should have been.",
            ).model_dump(exclude_none=True)
        )

    try:
        raw_metadata = get_raw_scoreset_metadata(urn, store_path)
        reference_sequences: dict[str, TargetAnnotation] = {}
        mapped_scores: list[ScoreAnnotation] = []
        for target_gene in annotated_vrs_results:
            preferred_layers = {
                _set_scoreset_layer(urn, annotated_vrs_results[target_gene]),
            }
            target_gene_name = metadata.target_genes[target_gene].target_gene_name
            reference_sequences[target_gene_name] = TargetAnnotation()
            reference_sequences[target_gene_name].layers = {
                layer: {
                    "computed_reference_sequence": None,
                    "mapped_reference_sequence": None,
                }
                for layer in preferred_layers
            }

            # sometimes Nonetype layers show up in preferred layers dict; remove these
            preferred_layers.discard(None)

            # Determine one gene symbol per target and its selection method
            gene_info = await compute_target_gene_info(
                target_key=target_gene,
                transcripts=transcripts,
                alignment_results=alignment_results,
                metadata=metadata,
                mapped_scores=annotated_vrs_results[target_gene],
            )

            for layer in preferred_layers:
                reference_sequences[target_gene_name].layers[layer][
                    "computed_reference_sequence"
                ] = _get_computed_reference_sequence(
                    metadata.target_genes[target_gene], layer, transcripts[target_gene]
                )
                reference_sequences[target_gene_name].layers[layer][
                    "mapped_reference_sequence"
                ] = _get_mapped_reference_sequence(
                    metadata.target_genes[target_gene],
                    layer,
                    transcripts[target_gene],
                    alignment_results[target_gene],
                )

            if gene_info is not None:
                reference_sequences[target_gene_name].gene_info = gene_info

            for m in annotated_vrs_results[target_gene]:
                if m.pre_mapped is None:
                    mapped_scores.append(ScoreAnnotation(**m.model_dump()))
                elif m.annotation_layer in preferred_layers:
                    # drop annotation layer from mapping object
                    mapped_scores.append(ScoreAnnotation(**m.model_dump()))

            # if genomic layer, not accession-based, and target gene type is coding, add cdna entry (just the sequence accession) to reference_sequences dict
            if (
                AnnotationLayer.GENOMIC in reference_sequences[target_gene_name].layers
                and metadata.target_genes[target_gene].target_gene_category
                == TargetType.PROTEIN_CODING
                and metadata.target_genes[target_gene].target_accession_id is None
                and transcripts[target_gene] is not None
                and isinstance(transcripts[target_gene], TxSelectResult)
                and transcripts[target_gene].nm is not None
            ):
                reference_sequences[target_gene_name].layers[AnnotationLayer.CDNA] = {
                    "computed_reference_sequence": None,
                    "mapped_reference_sequence": {
                        "sequence_accessions": [transcripts[target_gene].nm]
                    },
                }

            # drop Nonetype reference sequences
            for target_gene in reference_sequences:
                for layer in list(reference_sequences[target_gene].layers.keys()):
                    if (
                        reference_sequences[target_gene].layers[layer][
                            "mapped_reference_sequence"
                        ]
                        is None
                        and reference_sequences[target_gene].layers[layer][
                            "computed_reference_sequence"
                        ]
                        is None
                    ) or layer is None:
                        del reference_sequences[target_gene].layers[layer]

    except Exception as e:
        _logger.error("Unexpected error during result assembly for %s: %s", urn, e)
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    if len(mapped_scores) != total_score_records:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata,
                error_message=f"Mismatch between number of mapped scores ({len(mapped_scores)}) and total score records ({total_score_records}). This is unexpected and indicates an issue with the mapping process.",
            ).model_dump(exclude_none=True)
        )

    return JSONResponse(
        content=ScoresetMapping(
            metadata=raw_metadata,
            reference_sequences=reference_sequences,
            mapped_scores=mapped_scores,
        ).model_dump(exclude_none=True)
    )
