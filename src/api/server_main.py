"""FastAPI server file"""
import logging

import uvicorn
from fastapi import FastAPI

from api.routers import map
from dcd_mapping import dcd_mapping_version

logging.basicConfig(
    format="%(levelname)s:%(name)s:%(message)s",
    level=logging.INFO,
    force=True,
)
_logger = logging.getLogger(__name__)

app = FastAPI()

app.include_router(map.router)

msg = f"Starting DCD Mapping server v{dcd_mapping_version})"
_logger.info(msg)


# If the application is not already being run within a uvicorn server, start uvicorn here.
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)  # noqa: S104
