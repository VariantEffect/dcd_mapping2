# dcd-map: Map MaveDB data to computable and interoperable variant objects

[![image](https://img.shields.io/pypi/v/dcd_mapping.svg)](https://pypi.python.org/pypi/dcd_mapping)
[![image](https://img.shields.io/pypi/l/dcd_mapping.svg)](https://pypi.python.org/pypi/dcd_mapping)
[![image](https://img.shields.io/pypi/pyversions/dcd_mapping.svg)](https://pypi.python.org/pypi/dcd_mapping)
[![Actions status](https://github.com/ave-dcd/dcd_mapping/actions/workflows/checks.yaml/badge.svg)](https://github.com/ave-dcd/dcd_mapping/actions/checks.yaml)
[![DOI](https://zenodo.org/badge/472473437.svg)](https://zenodo.org/doi/10.5281/zenodo.11406657)

<!-- description -->

This library implements a novel method for mapping [MaveDB scoreset data](https://mavedb.org/) to [GA4GH Variation Representation Specification (VRS)](https://vrs.ga4gh.org/en/stable/) objects, enhancing interoperability for genomic medicine applications. See [Arbesfeld et. al. (2023)](https://www.biorxiv.org/content/10.1101/2023.06.20.545702v1) for a preprint edition of the mapping manuscript, or [download the resulting mappings directly](https://mavedb-mapping.s3.us-east-2.amazonaws.com/mappings.tar.gz).

<!-- /description -->

## Prerequisites

* Universal Transcript Archive (UTA): see [README](https://github.com/biocommons/uta?tab=readme-ov-file#installing-uta-locally) for setup instructions. Users with access to Docker on their local devices can use the available Docker image; otherwise, start a relatively recent (version 14+) PostgreSQL instance and add data from the available database dump.
* SeqRepo: see [README](https://github.com/biocommons/biocommons.seqrepo?tab=readme-ov-file#requirements) for setup instructions. The SeqRepo data directory must be writeable; see specific instructions [here](https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/store.rst) for more.
* Gene Normalizer: see [documentation](https://gene-normalizer.readthedocs.io/0.3.0-dev1/install.html) for data setup instructions.
* blat: Must be available on the local PATH and executable by the user. Otherwise, its location can be set manually with the `BLAT_BIN_PATH` env var. See the [UCSC Genome Browser FAQ](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3) for download instructions.


## Installation

Install from [PyPI](https://pypi.python.org/pypi/dcd_mapping):

```
python3 -m pip install dcd-mapping
```

## Usage

Use the `dcd-map` command with a scoreset URN, eg

```shell
$ dcd-map urn:mavedb:00000083-c-1
```

Output is saved in the format `<URN>_mapping_results_<ISO datetime>.json` in the directory specified by the environment variable `MAVEDB_STORAGE_DIR`, or `~/.local/share/dcd-mapping` by default.

Use `dcd-map --help` to see other available options.

## Mapping output

Each mapping run produces a single JSON document conforming to [`schema.json`](schema.json) (the JSON Schema serialization of `ScoresetMapping`). Top-level keys:

* `metadata` — the verbatim MaveDB API scoreset response, stored unchanged so no upstream fields are lost.
* `mapped_date` — ISO 8601 UTC timestamp of when this run completed.
* `reference_sequences` — per-target reference sequence info per annotation layer.
* `mapped_scores` — flat list of per-variant `ScoreAnnotation` records (see below).
* `target_mappings` — per-`(target, alignment_level)` provenance and alignment QC rows. The MaveDB API consumes these as `target_gene_mappings` and uses them to attribute every `mapped_score` back to the alignment that produced it.
* `error_message` — populated only when the run failed before producing scores.

### `metadata`

The verbatim MaveDB API scoreset response. Stored unchanged so downstream consumers retain access to every upstream field (URN, title, description, target gene definitions, score-column metadata, etc.) without having to query MaveDB again.

### `reference_sequences`

A `dict[target_gene_name, TargetAnnotation]` describing the reference sequences each target was mapped against, organized by annotation layer. Each `TargetAnnotation` carries:

* `gene_info` — `hgnc_symbol` plus the `selection_method` that picked it (transcript-derived, alignment-overlap-derived, variant-overlap-derived, or metadata fallback).
* `layers` — a `dict[AnnotationLayer, {computed_reference_sequence, mapped_reference_sequence}]` populated only for layers that actually produced mappings. `computed_reference_sequence` is the in-pipeline sequence (e.g. translated protein); `mapped_reference_sequence` lists the canonical accession(s) the variants were ultimately grounded in. Layers with no usable reference are pruned, not emitted as `null`.

This block is the human-readable "what was used as reference" view; programmatic auditing should use `target_mappings` instead.

### `mapped_scores`

A flat list of per-variant `ScoreAnnotation` records. One entry per `(score_record, emitted annotation_layer)` pair. Key fields:

* `mavedb_id`, `score` — identifier and numeric score copied from the MaveDB record.
* `relation` — fixed at `"SO:is_homologous_to"` while `pre_mapped` is populated.
* `target_gene_identifier`, `alignment_level` — composite key linking back to a `target_mappings` row (see below).
* `pre_mapped`, `post_mapped` — VRS variant objects in the target's coordinate frame and in the reference frame, respectively. Either may be `null` for failed mappings.
* `vrs_version` — VRS schema version used for this record.
* `error_message` — populated when `post_mapped` is `null` *or* when mapping succeeded with a caveat (e.g. RLE fallback, ambiguous reference allele).
* `at_mismatched_locus`, `near_gap` — per-variant audit flags, described below.

### `target_mappings`

Per-`(target, alignment_level)` provenance and alignment QC rows. The MaveDB API consumes these as `target_gene_mappings` and uses them to attribute every `mapped_score` back to the alignment that produced it. (See [`schema.json`](schema.json) `TargetMapping` for the wire format.)

### `error_message`

Populated only when the run failed before producing any scores; otherwise omitted. Per-variant errors live on `mapped_scores[].error_message`, not here.

## Audit and provenance details

### `target_mappings` fields

Each row describes the alignment that one set of mapped variants is grounded in:

| Field | Notes |
|---|---|
| `target_gene_identifier`, `alignment_level`, `preferred` | Composite key. `(target_gene_identifier, alignment_level)` is unique per run. Exactly one row per target has `preferred=True`. |
| `tool_name`, `tool_version`, `tool_parameters` | Aligner provenance. `tool_parameters.aligner` is `"blat"` for sequence-based targets and `"cdot_transcript_placement"` for accession-based targets. |
| `reference_accession`, `reference_sequence_id`, `vrs_version` | Coordinate-frame and run provenance. |
| `percent_identity`, `alignment_score`, `next_best_alignment_score`, `alignment_length`, `mismatch_count`, `gap_count` | Aggregate QC for the winning HSP. `alignment_score` is the canonical PSL score (`identities − mismatches − qNumInsert − tNumInsert`). |
| `alignment_string`, `alignment_metadata` | Pairwise visualization plus a small structured payload (CIGAR, `near_gap_window`, `at_mismatched_locus_evaluated`). |
| `total_variants`, `variants_mapped_cleanly`, `variants_with_mapping_warnings`, `variants_with_alignment_warnings`, `variants_failed` | Per-row variant counts. `variants_with_alignment_warnings` counts variants whose reference position fell on a mismatched base or near a gap. |

### Per-variant audit flags

Each `ScoreAnnotation` is attributable to exactly one `target_mappings` row via the composite key `(target_gene_identifier, alignment_level)`. The pipeline enforces this as a runtime invariant — orphaned scores raise `RuntimeError` rather than silently corrupting downstream joins.

Per-variant locus flags:

* `at_mismatched_locus` — `True` when any base in the variant's reference span mismatches between the target sequence and the reference; `False` when evaluated and no mismatch was found; `None` when per-base sequence content was unavailable for that layer (see `alignment_metadata.at_mismatched_locus_evaluated`), or when the variant is a `ReferenceLengthExpression` allele (large deletions/duplications, always `None`/`None`).
* `near_gap` — `True` when the variant lies within `alignment_metadata.near_gap_window` reference bases of any alignment gap; `None` for layers without an alignment (e.g. `cdna`).

Completely-failed variants (`pre_mapped is None` and no annotation layer was determined) are attributed to the target's preferred layer so the join invariant holds.

### Regenerating `schema.json`

`schema.json` is checked in and consumed by downstream services (notably the MaveDB API). After any change to `src/dcd_mapping/schemas.py` that alters the public output contract, regenerate it:

```shell
python scripts/generate_schema.py
```

Commit the regenerated `schema.json` in the same change.

## Notebooks

Notebooks for manuscript data analysis and figure generation are provided within `notebooks/analysis`. See [`notebooks/analysis/README.md`](notebooks/analysis/README.md) for more information.

Following installation instructions for [CoolSeqTool](https://coolseqtool.readthedocs.io/latest/install.html) and [Gene Normalizer](https://gene-normalizer.readthedocs.io/latest/install.html) should take care of the external data dependencies.

Note that Gene Normalizer's `pg` dependency group must be installed to make use of the PostgreSQL-based backend:

```shell
python3 -m pip install 'gene-normalizer[pg]'
```

## Development

Clone the repo

```
git clone https://github.com/ave-dcd/dcd_mapping
cd dcd_mapping
```

Create and activate a virtual environment

```
python3 -m virtualenv venv
source venv/bin/activate
```

Install as editable and with developer dependencies

```
python3 -m pip install -e '.[dev,tests]'
```

Add pre-commit hooks

```
pre-commit install
```

Run tests with `pytest`

```
pytest
```
