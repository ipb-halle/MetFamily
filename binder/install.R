# Install R packages

RUN for PACK in $PACK_R; do R -e "install.packages(\"$PACK\", repos='https://cran.r-project.org/')"; done

# Install Bioconductor packages
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite(\"BiocInstaller\", dep=TRUE, ask=FALSE)"
RUN for PACK in $PACK_BIOC; do R -e "library(BiocInstaller); biocLite(\"$PACK\", ask=FALSE)"; done

# Install other R packages from source
#RUN for PACK in $PACK_GITHUB; do R -e "library('devtools'); install_github(\"$PACK\")"; done

