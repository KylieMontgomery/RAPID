FROM continuumio/miniconda3:23.10.0-1
LABEL maintainer="Kylie Montgomery <kam219@cam.ac.uk>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml .

RUN conda env create -f environment.yml && \
    conda clean -afy

ENV PATH=/opt/conda/envs/rapid_env/bin:$PATH
