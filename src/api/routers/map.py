""""Provide mapping router"""
from pathlib import Path

from cool_seq_tool.schemas import AnnotationLayer
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from requests import HTTPError

from dcd_mapping.align import AlignmentError, BlatNotFoundError, build_alignment_result
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
    with_mavedb_score_set,
)
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreAnnotation,
    ScoresetMapping,
    TargetType,
    TxSelectResult,
    VrsVersion,
)
from dcd_mapping.transcripts import select_transcripts
from dcd_mapping.vrs_map import VrsMapError, vrs_map

router = APIRouter(
    prefix="/api/v1", tags=["mappings"], responses={404: {"description": "Not found"}}
)


@router.post(path="/map/{urn}", status_code=200, response_model=ScoresetMapping)
@with_mavedb_score_set
async def map_scoreset(urn: str, store_path: Path | None = None) -> ScoresetMapping:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param store_path: optional path to save output at
    """
    try:
        metadata = get_scoreset_metadata(urn, store_path)
        records = get_scoreset_records(metadata, True, store_path)
    except ScoresetNotSupportedError as e:
        return ScoresetMapping(
            metadata=None,
            error_message=str(e).strip("'"),
        )
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        raise HTTPException(status_code=500, detail=msg) from e

    if not records:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata,
                error_message="Score set contains no variants to map",
            ).model_dump(exclude_none=True)
        )

    try:
        alignment_results = build_alignment_result(metadata, True)
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
    except ScoresetNotSupportedError as e:
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
        raise HTTPException(status_code=500, detail=msg) from e
    except DataLookupError as e:
        msg = f"Data lookup error occurred during transcript selection: {e}"
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
    except VrsMapError as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    # TODO this should instead check if all values in dict are none. or might not need this at all.
    if vrs_results is None or len(vrs_results) == 0:
        return ScoresetMapping(
            metadata=metadata,
            error_message="No variant mappings available for this score set",
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
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )
    # TODO this should instead check if all values in dict are none. or might not need this at all.
    if annotated_vrs_results is None or len(annotated_vrs_results) == 0:
        return ScoresetMapping(
            metadata=metadata,
            error_message="No annotated variant mappings available for this score set",
        )

    try:
        raw_metadata = get_raw_scoreset_metadata(urn, store_path)
        reference_sequences: dict[str, dict] = {}
        mapped_scores: list[ScoreAnnotation] = []
        for target_gene in annotated_vrs_results:
            preferred_layers = {
                _set_scoreset_layer(urn, annotated_vrs_results[target_gene]),
            }
            target_gene_name = metadata.target_genes[target_gene].target_gene_name
            reference_sequences[target_gene_name] = {
                layer: {
                    "computed_reference_sequence": None,
                    "mapped_reference_sequence": None,
                }
                for layer in preferred_layers
            }
            # sometimes Nonetype layers show up in preferred layers dict; remove these
            preferred_layers.discard(None)
            for layer in preferred_layers:
                reference_sequences[target_gene_name][layer][
                    "computed_reference_sequence"
                ] = _get_computed_reference_sequence(
                    metadata.target_genes[target_gene], layer, transcripts[target_gene]
                )
                reference_sequences[target_gene_name][layer][
                    "mapped_reference_sequence"
                ] = _get_mapped_reference_sequence(
                    metadata.target_genes[target_gene],
                    layer,
                    transcripts[target_gene],
                    alignment_results[target_gene],
                )

            for m in annotated_vrs_results[target_gene]:
                if m.pre_mapped is None:
                    mapped_scores.append(ScoreAnnotation(**m.model_dump()))
                elif m.annotation_layer in preferred_layers:
                    # drop annotation layer from mapping object
                    mapped_scores.append(ScoreAnnotation(**m.model_dump()))

            # drop Nonetype reference sequences
            for target_gene in reference_sequences:
                for layer in list(reference_sequences[target_gene].keys()):
                    if (
                        reference_sequences[target_gene][layer][
                            "mapped_reference_sequence"
                        ]
                        is None
                        and reference_sequences[target_gene][layer][
                            "computed_reference_sequence"
                        ]
                        is None
                    ) or layer is None:
                        del reference_sequences[target_gene][layer]

            # if genomic layer, not accession-based, and target gene type is coding, add cdna entry (just the sequence accession) to reference_sequences dict
            if (
                AnnotationLayer.GENOMIC in reference_sequences[target_gene_name]
                and metadata.target_genes[target_gene].target_gene_category
                == TargetType.PROTEIN_CODING
                and metadata.target_genes[target_gene].target_accession_id is None
                and transcripts[target_gene] is not None
                and isinstance(transcripts[target_gene], TxSelectResult)
                and transcripts[target_gene].nm is not None
            ):
                reference_sequences[target_gene_name][AnnotationLayer.CDNA] = {
                    "computed_reference_sequence": None,
                    "mapped_reference_sequence": {
                        "sequence_accessions": [transcripts[target_gene].nm]
                    },
                }

    except Exception as e:
        return JSONResponse(
            content=ScoresetMapping(
                metadata=metadata, error_message=str(e).strip("'")
            ).model_dump(exclude_none=True)
        )

    return JSONResponse(
        content=ScoresetMapping(
            metadata=raw_metadata,
            reference_sequences=reference_sequences,
            mapped_scores=mapped_scores,
        ).model_dump(exclude_none=True)
    )
