# Notes for developpers

## Docker images
The image `sneumann/metfamily-base` contains all dependencies
for the MetFamily web application.

### Building the container(s)
Build via `docker build -t sneumann/metfamily-base -f Dockerfile-base . `or with the correct tagging:
```
echo docker build -t sneumann/metfamily-base:`grep ^FROM Dockerfile-base | cut -d: -f 2 | tr -d " "` -f Dockerfile-base .
```


The image `sneumann/metfamily` is built on top and contains
the actual MetFamily code and web application.

Build via `docker build -t sneumann/metfamily . ` or with the correct tagging:
```
docker build -t sneumann/metfamily:`grep ^FROM Dockerfile | cut -d: -f 2 | tr -d " "`-`grep Version DESCRIPTION | cut -d: -f 2 | tr -d " "`-`grep metFamilyAppVersion inst/MetFamily/version.R | cut -d'"' -f2` .
```

The build of the metfamily-base image https://hub.docker.com/r/sneumann/metfamily-base is triggerd whenever in the `master` branch a tag `basechange-<date>` e.g. `basechange-20190804` is specified.

### Running from a container

To run the resulting container, start with
`docker run --rm -p 3838:3838 sneumann/metfamily:latest`

and point your browser to http://localhost:3838/

### Developing and debugging in a container

If you want to develop and debug stuff, you can build a container
on top of `metfamily:latest` that has an added rstudio server.
First build using `docker build -t metfamily-rstudio -f Dockerfile-rstudio .`
and then run via `docker run -it --rm -p 8787:8787 metfamily-rstudio:latest`.
CAVEAT: the `Dockerfile-rstudio` specifies a fixed user/password combo
of `rstudio:rstudio`. Do not use in Production !

You can also pass a local directory with checked out MetFamily git tree
via the `docker run -v` argument.


