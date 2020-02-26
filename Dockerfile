FROM sneumann/metfamily-base:latest

MAINTAINER Kristian Peters <kpeters@ipb-halle.de>

LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

WORKDIR /srv/shiny-server

RUN rm -rf *

ADD MetFamily /srv/shiny-server/


