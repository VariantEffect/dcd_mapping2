"""Provide basic utilities for fetching and storing external data."""
import logging
import os
import time
from pathlib import Path

import click
import requests
from tqdm import tqdm

_logger = logging.getLogger(__name__)

MAVEDB_API_KEY = os.environ.get("MAVEDB_API_KEY")
MAVEDB_BASE_URL = os.environ.get("MAVEDB_BASE_URL")

LOCAL_STORE_PATH = Path(
    os.environ.get(
        "DCD_MAPPING_RESOURCES_DIR", Path.home() / ".local" / "share" / "dcd_mapping"
    )
)
if not LOCAL_STORE_PATH.exists():
    LOCAL_STORE_PATH.mkdir(exist_ok=True, parents=True)


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
        url, stream=True, timeout=60, headers=authentication_header()
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


def request_with_backoff(
    method: str, url: str, backoff_limit: int = 5, backoff_wait: int = 10, **kwargs
) -> requests.Response:
    """Make HTTP request with exponential backoff on failure.
    This is a duplicate of the function with same name in MaveDB API codebase.

    :param method: HTTP method (e.g., 'GET', 'POST')
    :param url: URL to make request to
    :param backoff_limit: number of retry attempts
    :param backoff_wait: initial wait time between retries (in seconds)
    :param kwargs: additional keyword arguments to pass to the request
    :return: Response object from the successful request
    """
    attempt = 0
    while attempt <= backoff_limit:
        msg = f"Attempt {attempt+1} of {backoff_limit} for {method} {url}"
        _logger.debug(msg)
        try:
            response = requests.request(method=method, url=url, **kwargs)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as exc:
            msg = f"Request to {url} failed on attempt {attempt+1}."
            _logger.warning(msg, exc_info=exc)
            backoff_time = backoff_wait * (2**attempt)
            attempt += 1
            msg = f"Waiting {backoff_time} seconds before retrying."
            _logger.info(msg)
            time.sleep(backoff_time)
    msg = f"Request to {url} failed after {backoff_limit} attempts."
    raise requests.exceptions.RequestException(msg)
