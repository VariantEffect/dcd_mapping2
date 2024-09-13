"""FastAPI server file"""
import uvicorn
from fastapi import FastAPI

from api.routers import map

app = FastAPI()

app.include_router(map.router)


# If the application is not already being run within a uvicorn server, start uvicorn here.
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)  # noqa: S104
