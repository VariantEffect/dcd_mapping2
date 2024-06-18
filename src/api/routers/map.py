from cool_seq_tool.schemas import AnnotationLayer
from fastapi import APIRouter

from dcd_mapping.align import align
from dcd_mapping.annotate import (
    _get_computed_reference_sequence,
    _get_mapped_reference_sequence,
    _set_scoreset_layer,
    annotate,
)
from dcd_mapping.mavedb_data import (
    get_raw_scoreset_metadata,
    get_scoreset_metadata,
    get_scoreset_records,
)
from dcd_mapping.schemas import ScoreAnnotation, ScoresetMapping
from dcd_mapping.transcripts import select_transcript
from dcd_mapping.vrs_map import vrs_map

router = APIRouter(prefix="/api/v1", tags=["mappings"], responses={404: {"description": "Not found"}})

@router.post(path="/map/{urn}", status_code=200, response_model=ScoresetMapping)
async def map_scoreset(urn: str) -> ScoresetMapping:
    metadata = get_scoreset_metadata(urn)
    records = get_scoreset_records(urn, True)

    alignment_result = align(metadata, True)

    transcript = await select_transcript(metadata, records, alignment_result)

    vrs_results = vrs_map(metadata, alignment_result, records, transcript, True)

    # TODO raise server error if vrs_results is None
    if vrs_results is None:
        return None

    vrs_results = annotate(vrs_results, transcript, metadata)

    raw_metadata = get_raw_scoreset_metadata(urn)
    # TODO change vrs map back to always use only the preferred layer
    #preferred_layers = {mapping.annotation_layer for mapping in vrs_results}
    preferred_layers = {
        _set_scoreset_layer(urn, vrs_results),
    }

    reference_sequences = {
        layer: {"computed_reference_sequence": None, "mapped_reference_sequence": None}
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

    output = ScoresetMapping(
        metadata=raw_metadata,
        computed_protein_reference_sequence=reference_sequences[
            AnnotationLayer.PROTEIN
        ]["computed_reference_sequence"],
        mapped_protein_reference_sequence=reference_sequences[AnnotationLayer.PROTEIN][
            "mapped_reference_sequence"
        ],
        computed_genomic_reference_sequence=reference_sequences[
            AnnotationLayer.GENOMIC
        ]["computed_reference_sequence"],
        mapped_genomic_reference_sequence=reference_sequences[AnnotationLayer.GENOMIC][
            "mapped_reference_sequence"
        ],
        mapped_scores=mapped_scores,
    )

    return output
