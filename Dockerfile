#FROM sneumann/metfamily-base:latest
FROM sneumann/metfamily-base:4.4
#FROM sneumann/metfamily-base:4.3.2
#FROM sneumann/metfamily-base:4.0.5
#FROM sneumann/metfamily-base:3.6.3

LABEL maintainer="Steffen Neumann <sneumann@ipb-halle.de>"
LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

ADD . /tmp/MetFamily

RUN R CMD INSTALL /tmp/MetFamily

WORKDIR /srv/shiny-server
RUN rm -rf *
ADD inst/MetFamily /srv/shiny-server/
