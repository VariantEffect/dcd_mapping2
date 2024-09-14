"""Provide basic utilities for fetching and storing external data."""
import os
from pathlib import Path

import click
import requests
from tqdm import tqdm

MAVEDB_API_KEY = os.environ.get("MAVEDB_API_KEY")
MAVEDB_BASE_URL = os.environ.get("MAVEDB_BASE_URL")

LOCAL_STORE_PATH = Path(
    os.environ.get(
        "DCD_MAPPING_RESOURCES_DIR", Path.home() / ".local" / "share" / "dcd_mapping"
    )
)
if not LOCAL_STORE_PATH.exists():
    LOCAL_STORE_PATH.mkdir(exist_ok=True, parents=True)


class ResourceAcquisitionError(Exception):
    """Raise when resource acquisition fails."""


def authentication_header() -> dict | None:
    """Fetch with api key envvar, if available."""
    return {"X-API-key": MAVEDB_API_KEY} if MAVEDB_API_KEY is not None else None


def http_download(url: str, out_path: Path, silent: bool = True) -> Path:
    """Download a file via HTTP.

    :param url: location of file to retrieve
    :param out_path: location to save file to
    :param silent: show TQDM progress bar if true
    :return: Path if download successful
    :raise requests.HTTPError: if request is unsuccessful
    """
    if not silent:
        click.echo(f"Downloading {out_path.name} to {out_path.parents[0].absolute()}")
    with requests.get(
        url, stream=True, timeout=30, headers=authentication_header()
    ) as r:
        r.raise_for_status()
        total_size = int(r.headers.get("content-length", 0))
        with out_path.open("wb") as h:
            if not silent:
                with tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    desc=out_path.name,
                    ncols=80,
                ) as progress_bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            h.write(chunk)
                            progress_bar.update(len(chunk))
            else:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        h.write(chunk)
    return out_path
