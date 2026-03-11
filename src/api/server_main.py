"""FastAPI server file"""
import logging

import uvicorn
from fastapi import FastAPI
from fastapi.requests import Request
from fastapi.responses import JSONResponse

from api.routers import map
from application_logging import init_logging
from dcd_mapping import dcd_mapping_version

init_logging()
_logger = logging.getLogger(__name__)
_logger.info("dcd-mapping API: %s", dcd_mapping_version)

app = FastAPI()

app.include_router(map.router)


@app.exception_handler(Exception)
async def unhandled_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Log unexpected exceptions and return a generic 500 response."""
    _logger.exception(
        "Unhandled exception for %s %s: %s",
        request.method,
        request.url.path,
        exc,
    )
    return JSONResponse(
        status_code=500,
        content={"detail": "Internal server error"},
    )


msg = f"Starting DCD Mapping server v{dcd_mapping_version})"
_logger.info(msg)


# If the application is not already being run within a uvicorn server, start uvicorn here.
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)  # noqa: S104
