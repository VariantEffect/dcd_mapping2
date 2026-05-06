"""Provide basic utilities for fetching and storing external data."""

import logging
import os
import time
from pathlib import Path

import click
import httpx
from tqdm import tqdm

from dcd_mapping.exceptions import ResourceAcquisitionError

_logger = logging.getLogger(__name__)

# Common representations of missing/null data in CSV files
MISSING_VALUE_REPRESENTATIONS = frozenset(
    {
        "NA",
        "N/A",
        "na",
        "n/a",
        "NaN",
        "nan",
        "null",
        "NULL",
        "None",
        "none",
        "",
        "-",
        ".",
    }
)

MAVEDB_API_KEY = os.environ.get("MAVEDB_API_KEY")
MAVEDB_BASE_URL = os.environ.get("MAVEDB_BASE_URL")
ENSEMBL_API_URL = os.environ.get("ENSEMBL_API_URL", "https://rest.ensembl.org")  # TODO
CDOT_URL = os.environ.get("CDOT_URL", "cdot-rest:8000")

LOCAL_STORE_PATH = Path(
    os.environ.get(
        "DCD_MAPPING_RESOURCES_DIR", Path.home() / ".local" / "share" / "dcd_mapping"
    )
)
if not LOCAL_STORE_PATH.exists():
    LOCAL_STORE_PATH.mkdir(exist_ok=True, parents=True)


def is_missing_value(value: str | None) -> bool:
    """Check if a value represents missing/null data.

    This function recognizes multiple common representations of missing data
    that may appear in CSV files from external sources, making the codebase
    more resilient to upstream changes in NA representation.

    :param value: The value to check
    :return: True if the value represents missing data, False otherwise
    """
    if value is None:
        return True
    # Strip whitespace and check against known missing value representations
    return value.strip() in MISSING_VALUE_REPRESENTATIONS


def authentication_header() -> dict | None:
    """Fetch with api key envvar, if available."""
    return {"X-API-key": MAVEDB_API_KEY} if MAVEDB_API_KEY is not None else None


def http_download(url: str, out_path: Path, silent: bool = True) -> Path:
    """Download a file via HTTP.

    :param url: location of file to retrieve
    :param out_path: location to save file to
    :param silent: show TQDM progress bar if true
    :return: Path if download successful
    :raise httpx.HTTPStatusError: if request is unsuccessful
    """
    if not silent:
        click.echo(f"Downloading {out_path.name} to {out_path.parents[0].absolute()}")
    with httpx.stream("GET", url, timeout=60, headers=authentication_header()) as r:
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
                    for chunk in r.iter_bytes(chunk_size=8192):
                        if chunk:
                            h.write(chunk)
                            progress_bar.update(len(chunk))
            else:
                for chunk in r.iter_bytes(chunk_size=8192):
                    if chunk:
                        h.write(chunk)
    return out_path


def request_with_backoff(
    url: str, max_retries: int = 5, backoff_factor: float = 0.3, **kwargs
) -> httpx.Response:
    """HTTP GET with exponential backoff only for retryable errors.

    Retries on:
    - Connection timeout or connection errors
    - HTTP 5xx server errors
    - HTTP 429 rate limiting (respecting Retry-After when present)

    Immediately raises on other HTTP errors (e.g., 4xx client errors).
    """
    attempt = 0
    last_exception: Exception | None = None
    while attempt < max_retries:
        try:
            kwargs.setdefault("timeout", 60)
            response = httpx.get(url, **kwargs)
        except (httpx.TimeoutException, httpx.ConnectError) as e:
            last_exception = e
            # Retry on transient network failures
            if attempt == max_retries - 1:
                raise
            _logger.debug(
                "Transient network error fetching %s (attempt %d/%d): %s",
                url,
                attempt + 1,
                max_retries,
                e,
            )
            sleep_time = backoff_factor * (2**attempt)
            time.sleep(sleep_time)
            attempt += 1
            continue

        # If we have a response, decide retry based on status code
        status = response.status_code
        if 200 <= status < 300:
            return response

        # 429: Too Many Requests — optionally use Retry-After
        if status == 429:
            if attempt == max_retries - 1:
                response.raise_for_status()
            retry_after = response.headers.get("Retry-After")
            try:
                sleep_time = (
                    float(retry_after)
                    if retry_after is not None
                    else backoff_factor * (2**attempt)
                )
            except ValueError:
                _logger.debug(
                    "Invalid Retry-After header value: %s, using exponential backoff",
                    retry_after,
                )
                sleep_time = backoff_factor * (2**attempt)
            time.sleep(sleep_time)
            attempt += 1
            continue

        # 5xx: server errors — retry
        if 500 <= status < 600:
            if attempt == max_retries - 1:
                response.raise_for_status()
            _logger.debug(
                "Server error %d fetching %s (attempt %d/%d)",
                status,
                url,
                attempt + 1,
                max_retries,
            )
            sleep_time = backoff_factor * (2**attempt)
            time.sleep(sleep_time)
            attempt += 1
            continue

        # Non-retryable (e.g., 4xx other than 429): raise immediately
        response.raise_for_status()

    # Exhausted retries without success
    msg = f"Failed to fetch {url} after {max_retries} attempts"
    raise ResourceAcquisitionError(msg) from last_exception
