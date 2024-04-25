"""Provide basic utilities for fetching and storing external data."""
import os
from pathlib import Path

import click
import requests
from tqdm import tqdm

LOCAL_STORE_PATH = Path(
    os.environ.get(
        "DCD_MAPPING_RESOURCES_DIR", Path.home() / ".local" / "share" / "dcd_mapping"
    )
)
if not LOCAL_STORE_PATH.exists():
    LOCAL_STORE_PATH.mkdir(exist_ok=True, parents=True)


class ResourceAcquisitionError(Exception):
    """Raise when resource acquisition fails."""


def paginated_http_download(
    url: str, out_path: Path, silent: bool = False, page_limit: int = 5000
) -> Path:
    """Use pagination params to handle large files. Currently only needed for experiment 53.

    :param url: url to fetch
    :parm out_path: location to save to
    :param silent: if true, suppress console output
    :param page_limit: number of rows to fetch at once
    """
    if not silent:
        click.echo(
            f"Downloading {out_path.name} with pagination to {out_path.parents[0].absolute()}"
        )
    start = 0
    all_rows = []
    while True:
        params = {"start": start, "limit": page_limit}
        with requests.get(url, params=params, timeout=30) as r:
            r.raise_for_status()
            rows = r.text.split("\r\n")[:-1]
            header, results = rows[0], rows[1:]
            if not all_rows:
                all_rows += [header]
            num_rows = len(results)
            if num_rows == 0:
                break
            start += num_rows
            all_rows += results
    with out_path.open("w") as f:
        for line in all_rows:
            f.write(line)
    return out_path


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
    with requests.get(url, stream=True, timeout=30) as r:
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


def get_mapping_tmp_dir() -> Path:
    """Acquire app-specific "tmp" directory. It's not actually temporary because it's
    manually maintained, but we need a slightly more durable file location than what the
    system tmp directory can provide. Used for storing small, consistently-named files
    like the BLAT query and results files.

    :return: path to temporary file directory
    """
    tmp = LOCAL_STORE_PATH / "tmp"
    tmp.mkdir(exist_ok=True)
    return tmp
