FROM python:3.11

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
RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
RUN chmod +x blat

# set dcd_mapping resources directory and download reference file
WORKDIR /home/.local/share/dcd_mapping
ENV DCD_MAPPING_RESOURCES_DIR=/home/.local/share/dcd_mapping
# hgdownload.cse.ucsc.edu is a legacy alias whose TLS certificate only covers
# hgdownload.soe.ucsc.edu. Checksum from the same directory's md5sum.txt.
RUN curl -fLJO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit && \
    echo "dcc3ea27079aa6dc3f9deccd7275e0f8  hg38.2bit" | md5sum -c -

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

# Bake in the liftover chain files agct would otherwise download on first use.
# That download goes through certifi, which can drop roots UCSC's TLS chain
# relies on (VariantEffect/dcd_mapping2#117); curl uses the system trust store.
# The paths are what agct's wags-tails lookup expects for unversioned files.
# Checksums from each directory's md5sum.txt.
ENV WAGS_TAILS_DIR=/home/.local/share/wags_tails
RUN mkdir -p $WAGS_TAILS_DIR/ucsc-chainfile && cd $WAGS_TAILS_DIR/ucsc-chainfile && \
    curl -fLO https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz && \
    curl -fLO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz && \
    printf '%s\n' \
        "35887f73fe5e2231656504d1f6430900  hg19ToHg38.over.chain.gz" \
        "ff3031d93792f4cbb86af44055efd903  hg38ToHg19.over.chain.gz" | md5sum -c - && \
    gunzip -c hg19ToHg38.over.chain.gz > chainfile_hg19_to_hg38_.chain && \
    gunzip -c hg38ToHg19.over.chain.gz > chainfile_hg38_to_hg19_.chain && \
    rm hg19ToHg38.over.chain.gz hg38ToHg19.over.chain.gz

RUN mkdir /usr/src/app
WORKDIR /usr/src/app
COPY . .

RUN pip install -e '.[dev,tests]'
# use polars-lts-cpu to avoid issues with x86 emulation on arm machine
RUN pip install -U polars-lts-cpu
RUN python -m pytest -p no:cacheprovider --no-cov tests/test_liftover_chainfiles.py
# install gene normalizer with pg dependencies. TODO: can the pg dependencies be specified in pyproject.toml?
#RUN pip install 'gene-normalizer[pg]'

# not working, needs to happen after db volume is mounted
# ENV GENE_NORM_DB_URL=postgres://postgres:postgres@db:5432/gene_normalizer
# RUN echo "y" | gene_norm_update_remote

ENV PYTHONUNBUFFERED 1

ENV PYTHONPATH "${PYTHONPATH}:/usr/src/app/src"

# Tell Docker that we will listen on port 8000.
EXPOSE 8000

# At container startup, run the application using uvicorn.
CMD ["uvicorn", "api.server_main:app", "--host", "0.0.0.0", "--port", "8000"]
