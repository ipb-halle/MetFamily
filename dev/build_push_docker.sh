#!/usr/bin/bash

echo -n "Now building images at " && date && \
test -f Dockerfile-base && \
docker build -t sneumann/metfamily-base:`grep ^FROM Dockerfile-base | cut -d: -f 2 | tr -d " "` -f Dockerfile-base . && \
test -f Dockerfile && \
docker build -t sneumann/metfamily:`grep ^FROM Dockerfile | cut -d: -f 2 | tr -d " "`-`grep Version DESCRIPTION | cut -d: -f 2 | tr -d " "`-`grep metFamilyAppVersion inst/MetFamily/version.R | cut -d'"' -f2` . && \
echo -n "Now pushing to Docker Hub at " && date && \
docker push sneumann/metfamily-base:`grep ^FROM Dockerfile-base | cut -d: -f 2 | tr -d " "` && \
docker push sneumann/metfamily:`grep ^FROM Dockerfile | cut -d: -f 2 | tr -d " "`-`grep Version DESCRIPTION | cut -d: -f 2 | tr -d " "`-`grep metFamilyAppVersion inst/MetFamily/version.R | cut -d'"' -f2` && \
echo -n "Done at " && date && \
xdg-open https://hub.docker.com/r/sneumann/metfamily-base/tags && \
xdg-open https://hub.docker.com/r/sneumann/metfamily/tags 




