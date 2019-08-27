FROM rocker/shiny:latest

MAINTAINER Kristian Peters <kpeters@ipb-halle.de>

LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

RUN apt-get -y update && apt-get -y install \
  netcdf-bin libnetcdf-dev libdigest-sha-perl \
  xorg-dev libglu1-mesa-dev freeglut3-dev libgomp1 libxml2-dev gcc g++ libcurl4-gnutls-dev libssl-dev gdebi-core

ENV NETCDF_INCLUDE=/usr/include

WORKDIR /srv/shiny-server

ADD app /srv/shiny-server/
RUN R -e "source('binder/install.R')"

