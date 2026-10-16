#!/usr/bin/env python3
"""Regenerate schema.json from the current Pydantic models.

Run from the repository root after any change to src/dcd_mapping/schemas.py
that alters the public output contract:

    python scripts/generate_schema.py

The script writes the JSON Schema for :class:`~dcd_mapping.schemas.ScoresetMapping`
to ``schema.json`` at the repository root, overwriting the previous file. With
``--check`` it writes nothing and exits 1 if ``schema.json`` is out of date.

If you change schemas.py, re-run this script and commit the result.
"""

import json
import pathlib
import sys

from dcd_mapping.schemas import ScoresetMapping

REPO_ROOT = pathlib.Path(__file__).parent.parent
SCHEMA_PATH = REPO_ROOT / "schema.json"


def main() -> None:
    """Generate JSON Schema for ScoresetMapping and write it, or check it with --check."""
    schema = (
        json.dumps(ScoresetMapping.model_json_schema(), indent=2, sort_keys=True) + "\n"
    )
    if "--check" in sys.argv[1:]:
        if SCHEMA_PATH.read_text() != schema:
            msg = f"{SCHEMA_PATH} is out of date; run python scripts/generate_schema.py"
            print(msg)  # noqa: T201
            sys.exit(1)
        return
    SCHEMA_PATH.write_text(schema)
    print(f"Written: {SCHEMA_PATH}")  # noqa: T201


if __name__ == "__main__":
    main()
