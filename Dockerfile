FROM rocker/shiny:latest

MAINTAINER Kristian Peters <kpeters@ipb-halle.de>

LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

RUN apt-get -y update && apt-get -y install \
  netcdf-bin libnetcdf-dev libdigest-sha-perl \
  xorg-dev libglu1-mesa-dev freeglut3-dev libgomp1 libxml2-dev gcc g++ libcurl4-gnutls-dev libssl-dev gdebi-core

ENV NETCDF_INCLUDE=/usr/include
RUN echo 'sanitize_errors off;disable_protocols xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf

ADD biner /tmp
RUN R -e "source('/tmp/install.R')"

WORKDIR /srv/shiny-server

RUN rm -rf *

ADD MetFamily /srv/shiny-server/


