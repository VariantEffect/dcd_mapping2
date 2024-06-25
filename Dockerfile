FROM --platform=x86_64 python:3.11

RUN apt update
# Install tools necessary used to install samtools and htslib so we can configure fasta files for genomic assembly.
RUN apt-get clean && apt-get update && apt-get install -y \
    postgresql-client \
	build-essential \
	curl \
	git \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
	libperl-dev \
	libssl-dev \
	zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# download and install blat executable
WORKDIR /usr/bin
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
RUN chmod +x blat

# set dcd_mapping resources directory and download reference file
WORKDIR /home/.local/share/dcd_mapping
ENV DCD_MAPPING_RESOURCES_DIR=/home/.local/share/dcd_mapping
RUN curl -LJO https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit

# Install samtools and htslib.
ARG htsversion=1.19
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig && \
    curl -L https://github.com/samtools/samtools/releases/download/${htsversion}/samtools-${htsversion}.tar.bz2 | tar xj && \
    (cd samtools-${htsversion} && ./configure --with-htslib=system && make install) && \
    curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && make install)

# Fetch and index GRCh37 and GRCh38 assemblies for cdot
RUN wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz | gzip -d | bgzip >  GCF_000001405.25_GRCh37.p13_genomic.fna.gz
RUN wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz | gzip -d | bgzip > GCF_000001405.39_GRCh38.p13_genomic.fna.gz
RUN samtools faidx GCF_000001405.25_GRCh37.p13_genomic.fna.gz
RUN samtools faidx GCF_000001405.39_GRCh38.p13_genomic.fna.gz

RUN mkdir /usr/src/app
WORKDIR /usr/src/app
COPY . .

RUN pip install -e '.[dev,tests]'
# use polars-lts-cpu to avoid issues with x86 emulation on arm machine
RUN pip install -U polars-lts-cpu
# install gene normalizer with pg dependencies. TODO: can the pg dependencies be specified in pyproject.toml?
#RUN pip install 'gene-normalizer[pg]'
ENV PYTHONUNBUFFERED 1

ENV PYTHONPATH "${PYTHONPATH}:/usr/src/app/src"
