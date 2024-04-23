FROM --platform=x86_64 python:3.11

RUN apt update
RUN apt install -y postgresql-client

# download and install blat executable
WORKDIR /usr/bin
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
RUN chmod +x blat

RUN mkdir /usr/src/app
WORKDIR /usr/src/app
COPY ./pyproject.toml .
RUN pip install -e '.[dev,tests]'
# use polars-lts-cpu to avoid issues with x86 emulation on arm machine
RUN pip install -U polars-lts-cpu
# install gene normalizer with pg dependencies. TODO: can the pg dependencies be specified in pyproject.toml?
#RUN pip install 'gene-normalizer[pg]'
ENV PYTHONUNBUFFERED 1
COPY . .
ENV PYTHONPATH "${PYTHONPATH}:/usr/src/app/src"

