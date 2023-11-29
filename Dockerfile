FROM sneumann/metfamily-base:latest

MAINTAINER Kristian Peters <kpeters@ipb-halle.de>

LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

ADD . /tmp/MetFamily

RUN R CMD INSTALL /tmp/MetFamily

WORKDIR /srv/shiny-server
RUN rm -rf *
ADD inst/MetFamily /srv/shiny-server/
