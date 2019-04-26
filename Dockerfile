FROM rocker/shiny:latest

MAINTAINER Kristian Peters <kpeters@ipb-halle.de>

LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

RUN apt-get -y update && apt-get -y install \
  netcdf-bin libnetcdf-dev libdigest-sha-perl \
  xorg-dev libglu1-mesa-dev freeglut3-dev libgomp1 libxml2-dev gcc g++ libcurl4-gnutls-dev libssl-dev gdebi-core

ADD binder/install.R /tmp

ENV PACK_R="BiocManager cba colourpicker devtools DT FactoMineR htmltools Matrix matrixStats plotrix rCharts rmarkdown shiny shinyBS shinyjs squash stringi tools"
ENV PACK_BIOC="mzR pcaMethods xcms"
ENV PACK_GITHUB=""

# Install R packages
RUN for PACK in $PACK_R; do R -e "install.packages(\"$PACK\", repos='https://cran.r-project.org/')"; done

# Install Bioconductor packages
RUN for PACK in $PACK_BIOC; do R -e "BiocManager::install("\"$PACK\"", ask=FALSE)"; done

# Install other R packages from source
#RUN for PACK in $PACK_GITHUB; do R -e "library('devtools'); install_github(\"$PACK\")"; done


#RUN R -e "source('/tmp/install.R')"

ADD . /srv/shiny-server/
