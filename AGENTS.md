# AGENTS.md — dcd_mapping

Facts about this repo that you can't work out from the code. Keep this file short. Rules shared with the other MaveDB repos are in the workspace `AGENTS.md` one folder up. Update both in the same change that makes any of it untrue.

The mapper turns a MaveDB score set's variants into GA4GH VRS objects. MaveDB's worker calls it over HTTP (`POST /api/v1/map/{urn}`, served from `src/api/`); the `dcd-map` CLI does the same job locally. This repo is MaveDB's fork of `ave-dcd/dcd_mapping`: work goes to `VariantEffect/dcd_mapping2`, branched from and merged into `mavedb-dev`, which is released to `mavedb-main`.

## Setup

- `python3 -m venv .venv && .venv/bin/pip install -e '.[dev,tests]'`, on Python 3.11. The `bin/` scripts use `.venv` when it exists.
- A real mapping run needs UTA, SeqRepo, the Gene Normalizer database and BLAT. The README's Prerequisites section covers each.

## Commands

| Task | Command |
|---|---|
| Lint, format and `schema.json` checks | `bin/check` |
| What CI tests | `bin/t tests/test_mavedb_data.py` |
| Other tests | `bin/t tests/test_variant_projection.py` |
| Regenerate `schema.json` | `.venv/bin/python scripts/generate_schema.py` |

- Run `bin/check` after every edit, and the tests for what you touched before handing work back.
- Only `tests/test_mavedb_data.py` runs in CI. The align and `vrs_map` tests need BLAT, and others need SeqRepo or UTA, so some fail on a machine without them. Check a failure against `mavedb-dev` before blaming your change, and report it rather than skipping or loosening the test.

## Rules

- Open pull requests with `gh pr create -R VariantEffect/dcd_mapping2 --base mavedb-dev`. Without `-R`, `gh` may target the upstream `ave-dcd` repo.
- `schema.json` is the output contract the MaveDB api consumes. A change to `schemas.py` that alters it ships with the regenerated file in the same commit; `bin/check` fails otherwise.
- `ga4gh.vrs` is pinned to 2.0.0a6 and `cool-seq-tool` to 0.4.0.dev3, because upgrading either means moving the whole GA4GH stack (gene-normalizer, seqrepo, cdot) at once. Don't bump any of them as part of other work.
- Never store a RefSeq sequence in SeqRepo from another source, such as cdot's genome-assembled transcript. It changes the refget, and with it every allele digest MaveDB has published.
- Docstrings use the reST `:param:` style the code already has, not Google style.
