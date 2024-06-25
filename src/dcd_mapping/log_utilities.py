"""Shared logging functions for dcd_mapping"""

import logging

import click

_logger = logging.getLogger(__name__)


def _emit_info(msg: str, silent: bool, log_level: int = logging.INFO) -> None:
    if not silent:
        click.echo(msg)
    if log_level == logging.INFO:
        _logger.info(msg)
    elif log_level == logging.ERROR:
        _logger.error(msg)
    else:
        msg = f"Unexpected log level requested: {log_level}"
        raise ValueError(msg)
