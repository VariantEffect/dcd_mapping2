#!/usr/bin/env python3
"""Regenerate schema.json from the current Pydantic models.

Run from the repository root after any change to src/dcd_mapping/schemas.py
that alters the public output contract:

    python scripts/generate_schema.py

The script writes the JSON Schema for :class:`~dcd_mapping.schemas.ScoresetMapping`
to ``schema.json`` at the repository root, overwriting the previous file.

If you change schemas.py, re-run this script and commit the result.
"""
import json
import pathlib

from dcd_mapping.schemas import ScoresetMapping

REPO_ROOT = pathlib.Path(__file__).parent.parent
SCHEMA_PATH = REPO_ROOT / "schema.json"


def main() -> None:
    """Generate JSON Schema for ScoresetMapping and write to file."""
    schema = ScoresetMapping.model_json_schema()
    SCHEMA_PATH.write_text(json.dumps(schema, indent=2, sort_keys=True) + "\n")
    print(f"Written: {SCHEMA_PATH}")  # noqa: T201


if __name__ == "__main__":
    main()
