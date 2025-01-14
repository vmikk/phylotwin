FROM rocker/r-ver:4.4.2

## Metadata
LABEL maintainer="vladimir.mikryukov@ut.ee" \
      version="v0.5" \
      description="PhyloNext/PhyloTwin diversity-analysis container"

## Accept GITHUB_PAT as a build argument
ARG GITHUB_PAT

## Set environment variables
ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    DEBIAN_FRONTEND=noninteractive \
    GITHUB_PAT=$GITHUB_PAT

## Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl wget git less \
        default-jre parallel \
    && rm -rf /var/lib/apt/lists/*

## Install R packages
RUN install2.r --error --skipmissing --skipinstalled --ncpus -1 \
        optparse \
        R.utils \
        glue \
        data.table \
        ape \
        phangorn \
        phytools \
        plyr \
        rotl \
        rgbif \
        sf \
        terra \
        leaflet \
        leaflet.extras \
        leaflet.extras2 \
        leafsync \
        webshot \
        canaper \
        arrow \
        dplyr \
        duckdb \
        qs \
        future \
        remotes \
    && R -e 'remotes::install_github("crazycapivara/h3-r")' \
    && R -e 'remotes::install_github("cran/PhyloMeasures")'

## Install phyloregion from GitHub
RUN git clone --depth 1 https://github.com/darunabas/phyloregion.git \
    && cd phyloregion \
    && rm build/vignette.rds \
    && R -e 'remotes::install_local(build_manual = FALSE, build_vignettes = FALSE)' \
    && cd .. \
    && rm -rf phyloregion

## Download and setup DuckDB with extensions
RUN curl -L https://github.com/duckdb/duckdb/releases/download/v1.1.3/duckdb_cli-linux-amd64.zip -o duckdb_cli-linux-amd64.zip \
    && unzip duckdb_cli-linux-amd64.zip -d /usr/local/bin \
    && rm duckdb_cli-linux-amd64.zip \
    && echo "SET GLOBAL extension_directory='/usr/local/bin/duckdb_ext'; INSTALL arrow; INSTALL spatial; INSTALL h3 FROM community;" | /usr/local/bin/duckdb

# Set the default command
CMD ["bash"] 
