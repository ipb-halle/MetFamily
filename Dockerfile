#FROM sneumann/metfamily-base:latest
FROM sneumann/metfamily-base:4.5
#FROM sneumann/metfamily-base:4.3.2
#FROM sneumann/metfamily-base:4.0.5
#FROM sneumann/metfamily-base:3.6.3

LABEL maintainer="Steffen Neumann <sneumann@ipb-halle.de>"
LABEL Description="MetFamily helps identifying metabolites and groups them into metabolite clusters (a.k.a. families)."

## Additional steps for Galaxy integration
## pip for installation, netstat as bioblend dependency
RUN apt-get update && apt-get install -y python3-pip net-tools\
    && apt-get clean

# Install Python dependencies using pip
# Use Steffen's patch https://github.com/bgruening/galaxy_ie_helpers/pull/15
# until it is properly merged and available through pypi
RUN pip3 install --break-system-packages "galaxy-ie-helpers>=0.2.8" || pip3 install --break-system-packages "https://github.com/sneumann/galaxy_ie_helpers/archive/refs/heads/master.zip"

# Install R package
RUN R -e 'devtools::install_github("hexylena/rGalaxyConnector")'
# Needed for gx_get:
RUN mkdir -p /import; chmod 777 /import

## Now install real MetFamily
ADD . /tmp/MetFamily
RUN R CMD INSTALL /tmp/MetFamily

WORKDIR /srv/shiny-server
RUN rm -rf *
ADD inst/MetFamily /srv/shiny-server/
