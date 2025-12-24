
# Base image with R and Shiny Server
FROM rocker/shiny:latest

## Install system dependencies
RUN apt-get update && \
	apt-get -y upgrade

RUN	apt-get install -y --no-install-recommends build-essential && \
	apt-get install -y --no-install-recommends autotools-dev && \
	apt-get install -y --no-install-recommends autoconf && \
    apt-get install -y --no-install-recommends automake && \
	apt-get install -y --no-install-recommends libtool

RUN	apt-get install -y --no-install-recommends zlib1g-dev && \
	apt-get install -y --no-install-recommends libgsl-dev && \
	apt-get install -y --no-install-recommends libxml2-dev
	
RUN apt-get install -y --no-install-recommends libglu1-mesa && \
	apt-get install -y --no-install-recommends libssl-dev && \
	apt-get install -y --no-install-recommends vim-tiny

RUN rm -rf /var/lib/apt/lists/*


# C code
WORKDIR /build/POPSC
COPY POPSC/ .

## Run the Autotools bootstrap + build
RUN ./bootstrap && \
    ./configure && \
    make && \
    make install DESTDIR=/build/install


# R code
RUN R -e "install.packages(c('shiny', 'bio3d', 'data.table', 'DT', 'rgl', \
							 'digest', 'devtools', 'ggplot2', 'optparse'))"

RUN R -e "devtools::install_github('AnalytixWare/ShinySky')"

## Compile and install POPSR
COPY POPSR /build/POPSR
WORKDIR /build/POPSR
RUN R CMD build .
RUN R CMD INSTALL *.tar.gz


# Shiny App 
## Create app directory
WORKDIR /app
COPY POPSR/inst/popsr/app.R /app/

WORKDIR /app/www
COPY POPSR/inst/popsr/www/POPScomp.png /app/www

## Expose Shiny port
EXPOSE 3838

## Run Shiny server
#CMD ["/usr/bin/shiny-server"]
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]

