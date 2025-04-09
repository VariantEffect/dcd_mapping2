"""Provide class definitions for commonly-used information objects."""
import datetime
from enum import Enum
from typing import Any, Literal

from cool_seq_tool.schemas import AnnotationLayer, Strand, TranscriptPriority
from ga4gh.vrs._internal.models import Allele, Haplotype
from pydantic import BaseModel, ConfigDict, Field, StrictBool, StrictInt, StrictStr

from dcd_mapping import vrs_v1_schemas
from dcd_mapping.version import dcd_mapping_version


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class TargetType(str, Enum):
    """Define target gene types."""

    PROTEIN_CODING = "protein_coding"
    REGULATORY = "regulatory"
    OTHER_NC = "other_noncoding"


class VrsVersion(str, Enum):
    """Define VRS versions"""

    V_1_3 = "1.3"
    V_2 = "2"


class UniProtRef(BaseModel):
    """Store metadata associated with MaveDB UniProt reference"""

    id: str
    offset: int


class TargetGene(BaseModel):
    """Store metadata for a target gene from a MaveDB score set"""

    target_gene_name: str
    target_gene_category: TargetType
    target_sequence: str | None = None
    target_sequence_type: TargetSequenceType | None = None
    target_sequence_label: str | None = None
    target_uniprot_ref: UniProtRef | None = None
    target_accession_id: str | None = None
    target_accession_assembly: str | None = None


class ScoresetMetadata(BaseModel):
    """Store all relevant metadata from metadata reported for scoreset by MaveDB"""

    urn: str
    target_genes: dict[str, TargetGene]


class ScoreRow(BaseModel):
    """Row from a MAVE score result"""

    hgvs_pro: str
    hgvs_nt: str
    score: str | None
    accession: str


class SequenceRange(BaseModel):
    """Define range over a sequence. Useful for expressing alignment query and hit results."""

    start: int
    end: int


class GeneLocation(BaseModel):
    """Gene location info, gathered from normalizer result. Likely to be incomplete."""

    chromosome: str | None = None
    start: int | None = None
    end: int | None = None


class ReferenceSequence(BaseModel):
    """Base reference sequence class."""

    sequence_type: TargetSequenceType
    sequence_id: StrictStr


class ComputedReferenceSequence(ReferenceSequence):
    """Define metadata describing a computed reference sequence"""

    sequence: StrictStr


class MappedReferenceSequence(ReferenceSequence):
    """Define metadata describing a mapped, human reference sequence"""

    sequence_accessions: list[StrictStr]


class AlignmentResult(BaseModel):
    """Define BLAT alignment output."""

    chrom: str
    strand: Strand
    coverage: float | None = None
    ident_pct: float | None = None
    query_range: SequenceRange
    query_subranges: list[SequenceRange]
    hit_range: SequenceRange
    hit_subranges: list[SequenceRange]


class TranscriptDescription(BaseModel):
    """Structured transcript description.

    Provides less information than the MANE results, but should convey what we need.
    """

    refseq_nuc: str
    refseq_prot: str
    transcript_priority: TranscriptPriority


class ManeDescription(TranscriptDescription):
    """Structured MANE data retrieval result."""

    ncbi_gene_id: str
    ensembl_gene_id: str
    hgnc_gene_id: str
    symbol: str
    name: str
    ensembl_nuc: str
    ensembl_prot: str
    grch38_chr: str
    chr_start: int
    chr_end: int
    chr_strand: str


class TxSelectResult(BaseModel):
    """Define response object from transcript selection process."""

    nm: str | None = None
    np: str
    start: StrictInt
    is_full_match: StrictBool
    transcript_mode: TranscriptPriority | None = None
    sequence: str


class MappedScore(BaseModel):
    """Provide mappings for an individual experiment score.

    This model defines the output of the VRS mapping phase of the pipeline.
    """

    model_config = ConfigDict(use_enum_values=True)

    accession_id: StrictStr
    score: str | None
    annotation_layer: AnnotationLayer | None = None
    pre_mapped: Allele | Haplotype | None = None
    post_mapped: Allele | Haplotype | None = None
    error_message: str | None = None


class ScoreAnnotation(BaseModel):
    """Provide extra annotations on top of mappings for an individual experiment score.

    This model defines what an individual mapping instance looks like in the final JSON.
    """

    mavedb_id: StrictStr
    relation: Literal[
        "SO:is_homologous_to"
    ] = "SO:is_homologous_to"  # TODO this should probably be None if pre_mapped is false?
    pre_mapped: vrs_v1_schemas.VariationDescriptor | vrs_v1_schemas.Haplotype | Allele | Haplotype | None = None
    post_mapped: vrs_v1_schemas.VariationDescriptor | vrs_v1_schemas.Haplotype | Allele | Haplotype | None = None
    vrs_version: VrsVersion | None = None
    score: float | None = None
    error_message: str | None = None


class ScoreAnnotationWithLayer(ScoreAnnotation):
    """Couple annotations with an easily-computable definition of the annotation layer
    from which they originate.

    Used for filtering individual annotations just before saving the final JSON product.
    """

    annotation_layer: AnnotationLayer | None = None


class ScoresetMapping(BaseModel):
    """Provide all mapped scores for a scoreset."""

    metadata: Any  # TODO get exact MaveDB metadata structure?
    dcd_mapping_version: str = Field(default=dcd_mapping_version)
    mapped_date_utc: str = Field(
        default=datetime.datetime.now(tz=datetime.UTC).isoformat()
    )
    # TODO re-implement metadata change later to support multi-target score sets. will require corresponding changes in mavedb-api
    # reference_sequences: dict[
    #     str,
    #     dict[
    #         AnnotationLayer,
    #         dict[str, ComputedReferenceSequence | MappedReferenceSequence | None],
    #     ],
    # ] | None = None
    computed_protein_reference_sequence: ComputedReferenceSequence | None = None
    mapped_protein_reference_sequence: MappedReferenceSequence | None = None
    computed_genomic_reference_sequence: ComputedReferenceSequence | None = None
    mapped_genomic_reference_sequence: MappedReferenceSequence | None = None
    mapped_scores: list[ScoreAnnotation] | None = None
    error_message: str | None = None
