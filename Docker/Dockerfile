# Base image with R and Shiny Server
FROM rocker/shiny:latest

ENV DEBIAN_FRONTEND=noninteractive

## Install system dependencies and Gemmi CLI
## pkg-config is required because configure.ac uses PKG_CHECK_MODULES(...)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        autotools-dev \
        autoconf \
        automake \
        libtool \
        pkg-config \
        zlib1g-dev \
        libgsl-dev \
        libxml2-dev \
        libglu1-mesa \
        libssl-dev \
        vim-tiny \
        python3-pip \
        gemmi && \
    gemmi --version && \
    rm -rf /var/lib/apt/lists/*

# C code
WORKDIR /build/POPSC
COPY POPSC/ .

## Run the Autotools bootstrap + build
RUN ./bootstrap && \
    ./configure && \
    make -j$(nproc) && \
    make install DESTDIR=/build/install

# R packages
RUN R -e "install.packages(c('shiny', 'bio3d', 'data.table', 'DT', 'rgl', 'remotes', \
                             'digest', 'ggplot2', 'optparse', 'markdown'), \
                             repos='https://cloud.r-project.org')"

RUN R -e "remotes::install_github('AnalytixWare/ShinySky')"

## Compile and install POPSR
COPY POPSR /build/POPSR
WORKDIR /build/POPSR
RUN R CMD build . && \
    R CMD INSTALL *.tar.gz

# Shiny App
WORKDIR /app

COPY POPSR/inst/popsr/ /app/

## Expose Shiny port
EXPOSE 3838

## Run Shiny app
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
