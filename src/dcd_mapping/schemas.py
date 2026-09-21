"""Provide class definitions for commonly-used information objects."""

import datetime
from enum import StrEnum
from typing import Any, Literal, NamedTuple

from cool_seq_tool.schemas import AnnotationLayer, Strand, TranscriptPriority
from ga4gh.vrs._internal.models import Allele, Haplotype
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StrictBool,
    StrictInt,
    StrictStr,
    field_serializer,
)

from dcd_mapping import vrs_v1_schemas
from dcd_mapping.version import dcd_mapping_version


class TargetSequenceType(StrEnum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class TargetType(StrEnum):
    """Define target gene types."""

    PROTEIN_CODING = "protein_coding"
    REGULATORY = "regulatory"
    OTHER_NC = "other_noncoding"


class VrsVersion(StrEnum):
    """Define VRS versions"""

    V_1_3 = "1.3"
    V_2 = "2"


class MappingOutcome(StrEnum):
    """Per-record outcome for one (variant, annotation level) pair.

    The mapper's output is a complete accounting: for every variant and every
    annotation level in that variant's deterministically-reachable set, there is one
    record carrying its outcome -- never a silent omission. This field is uniform across
    measured (assay-level) and projected (deterministic non-assay) records so the two can
    be treated identically by consumers; it distinguishes a benign absence from a genuine
    failure, which a populated ``error_message`` alone cannot.

    - ``MAPPED`` -- a VRS allele was produced (``pre_mapped``/``post_mapped`` populated).
    - ``INTRONIC`` -- the variant's coding projection is intronic: no VRS-representable
      coding form and no protein consequence. Benign (``error_message`` is ``None``).
    - ``NO_PROTEIN_CONSEQUENCE`` -- the protein layer was reachable but yields no
      projectable protein change (e.g. UTR). Benign (``error_message`` is ``None``).
    - ``FAILED`` -- the mapping/projection genuinely failed (mis-selected transcript,
      projection error, unresolvable reference contig). ``error_message`` carries detail.
    """

    MAPPED = "mapped"
    INTRONIC = "intronic"
    NO_PROTEIN_CONSEQUENCE = "no_protein_consequence"
    FAILED = "failed"


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


class AlignmentQc(BaseModel):
    """Aggregate QC for a BLAT pslx alignment.

    Mismatch positions and gap intervals are kept on the in-memory model so
    downstream code can flag individual variants that fall on a mismatched
    base or near a gap, but they are excluded from JSON serialization to keep
    payloads compact -- consumers reconstruct per-base detail from the CIGAR
    or use the per-variant flags emitted on each :class:`ScoreAnnotation`.

    **cDNA→genome BLAT limitation:** BLAT's pslx format is supposed to include
    per-block target and query sequences (tseqs/qseqs), but in practice
    Biopython's PSL parser may leave aligned-block bases as undefined when
    the pslx output is incomplete or when blocks fall outside the portion of
    the database sequence that BLAT includes in its output.  When this occurs,
    ``mismatch_count`` (derived from ``Alignment.counts()``, which is purely
    coordinate-based) remains accurate, but per-position detail cannot be
    extracted.  The ``mismatch_positions_unavailable`` flag is set to ``True``
    in that case; consumers should treat ``at_mismatched_locus`` as ``None``
    (not evaluated) for affected rows rather than trusting a ``False`` value.

    A future improvement would be to fetch target slices directly via SeqRepo
    when per-base detail is required, decoupling per-position accuracy from
    what pslx happens to include.

    **Downstream effect** — ``annotate._stamp_alignment_locus_flags`` is the
    sole consumer of ``mismatch_positions`` and ``gap_intervals``.  When
    ``mismatch_positions_unavailable`` is ``True`` it skips the
    ``at_mismatched_locus`` assignment entirely (leaving it ``None``) so
    that a ``False`` value is never incorrectly written for a variant whose
    position could not actually be checked against the alignment.
    """

    percent_identity: float | None = None
    alignment_length: int
    mismatch_count: int
    gap_count: int
    cigar: StrictStr | None = None
    # Compact Biopython alignment visualization. Consecutive all-gap display
    # blocks (intronic regions rendered as '?' characters in cDNA→genome
    # alignments) are collapsed to a single summary line of the form
    # ``... [N bp gap: chrom:start-end] ...`` by ``_compact_alignment_string``
    # in align.py. This keeps the string proportional to exon count rather
    # than genomic span.
    alignment_string: StrictStr | None = None
    # True when per-base sequence content was unavailable during QC computation
    # (e.g. BLAT pslx did not populate tseqs/qseqs). mismatch_count is still
    # correct in that case, but mismatch_positions is empty and at_mismatched_locus
    # should be treated as None (not evaluated) for all variants in this target/level.
    mismatch_positions_unavailable: bool = False
    # In-memory only -- not serialized. Used to compute per-variant flags for variants that
    # fall on a mismatched base or near a gap.
    mismatch_positions: list[int] = Field(default_factory=list, exclude=True)
    gap_intervals: list[tuple[int, int]] = Field(default_factory=list, exclude=True)


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

    chrom: str | None = None
    strand: Strand | None = None
    coverage: float | None = None
    percent_identity: float | None = None
    query_range: SequenceRange
    query_subranges: list[SequenceRange]
    hit_range: SequenceRange
    hit_subranges: list[SequenceRange]
    # BLAT PSL score (matches - misMatches - qNumInsert - tNumInsert) for the
    # winning HSP. Only set for BLAT-derived results.
    score: float | None = None
    # BLAT PSL score (same units as ``score``) of the next-best HSP considered
    # during selection. Presence signals that an ambiguity check was performed;
    # None means no ambiguity check was run (not that there was no ambiguity).
    # For instance, when returning results for protein->protein alignments, BLAT
    # returns only the highest scoring HSP.
    next_best_score: float | None = None
    # Structured per-base QC derived from the pslx alignment (populated for
    # BLAT-sourced alignments; None for cdot-fetched results).
    alignment_qc: AlignmentQc | None = None
    # Aligner invocation parameters recorded for reproducibility.
    aligner_parameters: dict[str, Any] | None = None
    # Human reference genome assembly used for this alignment (e.g. "GRCh38").
    # None for alignments that have no genomic coordinate frame (protein-vs-protein).
    reference_assembly: StrictStr | None = None


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
    hgnc_symbol: str | None = None


class MappedScore(BaseModel):
    """Provide mappings for an individual experiment score.

    This model defines the output of the VRS mapping phase of the pipeline.
    """

    model_config = ConfigDict(use_enum_values=True)

    accession_id: StrictStr
    score: str | None
    alignment_level: AnnotationLayer | None = None
    pre_mapped: Allele | Haplotype | None = None
    post_mapped: Allele | Haplotype | None = None
    error_message: str | None = None
    # Typed outcome for this (variant, level) record. None until stamped (legacy /
    # pre-annotation); ``annotate`` resolves it for every emitted record.
    outcome: MappingOutcome | None = None


class ScoreAnnotation(BaseModel):
    """Provide extra annotations on top of mappings for an individual experiment score.

    This model defines what an individual mapping instance looks like in the final JSON.
    """

    mavedb_id: StrictStr
    relation: Literal["SO:is_homologous_to"] = (
        "SO:is_homologous_to"  # TODO this should probably be None if pre_mapped is false?
    )
    # Identifies which target gene this score belongs to. Required for the API
    # to join mapped_scores → target_mappings unambiguously when a score set has
    # more than one target gene (joining on alignment_level alone is ambiguous in
    # that case). Nullable for backwards compatibility with older pipeline outputs.
    target_gene_identifier: StrictStr | None = None
    pre_mapped: (
        vrs_v1_schemas.VariationDescriptor
        | vrs_v1_schemas.Haplotype
        | Allele
        | Haplotype
        | None
    ) = None
    post_mapped: (
        vrs_v1_schemas.VariationDescriptor
        | vrs_v1_schemas.Haplotype
        | Allele
        | Haplotype
        | None
    ) = None
    vrs_version: VrsVersion | None = None
    score: float | None = None
    error_message: str | None = None
    alignment_level: AnnotationLayer | None = None
    # Typed outcome for this (variant, level) record -- see MappingOutcome. Always
    # populated on emitted annotations; distinguishes a benign absence (intronic, no
    # protein consequence) from a genuine failure that error_message alone cannot.
    outcome: MappingOutcome | None = None
    # Per-variant alignment-locus flags. None means "not evaluated" (e.g. no
    # genomic alignment available, or non-genomic annotation layer); True/False
    # mean the flag was computed against this run's alignment.
    at_mismatched_locus: bool | None = None
    near_gap: bool | None = None


class GeneInfo(BaseModel):
    """Basic gene metadata for a target, including symbol and selection method."""

    hgnc_symbol: str | None = None
    selection_method: str | None = None


class TargetAnnotation(BaseModel):
    """Represents annotations associated with a biological target, including optional gene metadata
    and structured annotation layers.

    Attributes
    ----------
    gene_info : GeneInfo | None
        Optional metadata describing the gene associated with the target,
        including identifiers and descriptive information where available.

    layers : dict[AnnotationLayer, dict[str, ComputedReferenceSequence | MappedReferenceSequence | dict | None]]
        A mapping of annotation layers to keyed layer data. Each layer is identified by an
        AnnotationLayer key and contains a dictionary where:
          - keys are string identifiers for items within the layer (e.g., feature names),
          - values are one of:
              - ComputedReferenceSequence: a computed sequence representation for the item,
              - MappedReferenceSequence: a sequence mapped to a reference coordinate system,
              - dict: a generic dictionary for custom layer-specific payloads,
              - None: indicating missing or intentionally omitted data.

    Notes
    -----
    - The default value for 'layers' is an empty dictionary.
    - This model is intended to standardize layer-based annotations for downstream processing
      and validation, allowing both computed and mapped sequence data to coexist within the same
      structure.

    """

    gene_info: GeneInfo | None = None
    layers: dict[
        AnnotationLayer,
        dict[str, ComputedReferenceSequence | MappedReferenceSequence | dict | None],
    ] = {}


class ScoresetMapping(BaseModel):
    """Provide all mapped scores for a scoreset."""

    metadata: Any  # TODO get exact MaveDB metadata structure?
    mapped_date: datetime.datetime | None = None
    reference_sequences: dict[str, TargetAnnotation] | None = None
    mapped_scores: list[ScoreAnnotation] | None = None
    target_mappings: list["TargetMapping"] | None = None
    error_message: str | None = None

    @field_serializer("mapped_date")
    def _serialize_mapped_date(self, value: datetime.datetime | None) -> str | None:
        """Serialize mapped_date as an ISO 8601 string so it is always JSON-safe."""
        return value.isoformat() if value is not None else None


class TargetMapping(BaseModel):
    """Per-(target, alignment_level) provenance and QC block.

    Field names mirror the corresponding columns on `target_gene_mappings` in the
    MaveDB API so the API worker can deserialize directly with minimal transformation.
    Any aligner-specific structured details go in ``tool_parameters`` /
    ``alignment_metadata`` (JSONB on the API side).

    ``reference_assembly`` is a top-level column (not nested in ``tool_parameters``)
    because it describes the coordinate frame of the mapping result, not aligner
    configuration.  It is ``None`` for alignments with no genomic frame (e.g.
    protein-vs-protein).

    ``tool_parameters`` shape per aligner
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLAT genomic (sequence-based, nucleotide or protein-vs-genome)::

        {
          "aligner":         "blat",
          "aligner_version": "<BLAT version string | null>",
          "min_score":       <int>,
          "out_format":      "<pslx|...>",
          "target_args":     "<e.g. '-q=prot -t=dnax' | ''>",
        }

    BLAT protein-vs-protein (sequence-based, protein annotation layer)::

        {
          "aligner":         "blat",
          "aligner_version": "<BLAT version string | null>",
          "min_score":       <int>,
          "out_format":      "<pslx|...>",
          "target_args":     "-q=prot -t=prot",
        }

    Passthrough for fully qualified reference accessions (accession-based)::

        {
          "aligner":         "reference_accession_passthrough",
        }

    cdot transcript placement (accession-based)::

        {
          "aligner":           "cdot_transcript_placement",
          "aligner_version":   "<cdot library semver, e.g. '0.2.26'>",
          "cdot_url":          "<CDOT REST endpoint>",
          "cdot_data_version": "<data version string | null>",
        }

    ``cdot_data_version`` is ``null`` when the cdot REST response did not
    include a ``cdot_data_version`` field; its presence (even as ``null``)
    distinguishes "cdot run, version unknown" from "not a cdot run".
    """

    # Identity
    target_gene_identifier: StrictStr
    # Single-char cool-seq-tool AnnotationLayer value ("p", "c", "g").
    alignment_level: AnnotationLayer
    preferred: StrictBool = False

    # Provenance (required on the API side: tool_name, tool_version)
    tool_name: StrictStr = "dcd-mapping"
    tool_version: StrictStr = dcd_mapping_version
    tool_parameters: dict[str, Any] | None = None
    # Human reference genome assembly (e.g. "GRCh38"). None for protein-vs-protein alignments.
    reference_assembly: StrictStr | None = None
    reference_accession: StrictStr | None = None
    reference_sequence_id: StrictStr | None = None
    vrs_version: VrsVersion | None = None

    # Alignment QC -- all optional; omit (leave None) fields the aligner doesn't produce.
    percent_identity: float | None = None
    alignment_score: float | None = None
    next_best_alignment_score: float | None = None
    alignment_length: int | None = None
    mismatch_count: int | None = None
    gap_count: int | None = None
    # Compact alignment visualization (see AlignmentQc.alignment_string).
    # Intronic gap runs are collapsed to ``... [N bp gap: chrom:start-end] ...``
    # summaries so the string length is proportional to exon count, not
    # genomic span. None for protein-layer or cdot-derived alignments.
    alignment_string: StrictStr | None = None
    # Structured per-alignment metadata payload (JSONB on the API side).
    # Keys present when available:
    #   "cigar"                     — CIGAR string for the winning HSP.
    #   "near_gap_window"           — half-width (in ref bases) of the window
    #                                 used to flag ``near_gap`` on each variant.
    #   "at_mismatched_locus_evaluated" — False when per-base sequence content
    #                                 was unavailable (BLAT pslx omitted
    #                                 tseqs/qseqs); treat ``at_mismatched_locus``
    #                                 as None (not evaluated) for all variants in
    #                                 this target/level.  True otherwise.
    alignment_metadata: dict[str, Any] | None = None

    # Annotation QC -- totals for this target x level.
    # variants_with_alignment_warnings: the variant's reference position fell on
    # a mismatched locus or near a gap in the underlying alignment
    # (``at_mismatched_locus`` or ``near_gap`` was True). Mapping itself
    # succeeded cleanly; this is purely about alignment context.
    total_variants: int | None = None
    # Count of variants where post_mapped is not None.
    # Alignment warnings (at_mismatched_locus, near_gap) are counted separately in
    # variants_with_alignment_warnings and do NOT reduce this count; a variant can
    # be "cleanly mapped" and still sit at a mismatched locus or near a gap.
    variants_mapped_cleanly: int | None = None
    variants_with_alignment_warnings: int | None = None
    variants_failed: int | None = None


class VrsMapResult(NamedTuple):
    """Provide mappings for an individual experiment score."""

    mappings: list["MappedScore"]
    protein_align_result: "AlignmentResult | None"
