"""Application logging initialization"""

import logging
import os

LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()


def init_logging() -> None:
    """Initialize application-wide logging with configured log level and format.

    This sets the root logger's level based on the LOG_LEVEL environment variable
    and applies a consistent log message format across the application.
    """
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=getattr(logging, LOG_LEVEL, logging.INFO),
        force=True,
    )
