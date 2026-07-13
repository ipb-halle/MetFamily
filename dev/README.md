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

## MetFamily as Galaxy Interactive tool

### Preparing a local Galaxy 

Checkout the Galaxy project:

```
git clone https://github.com/galaxyproject/galaxy.git
git checkout release_26.1	# just to be on a stable side ...
cd galaxy
```

Use config examples that work on Ubuntu 26.04:
```
cp config/galaxy.yml.interactivetools config/galaxy.yml
cp config/job_conf.yml.interactivetools config/job_conf.yml
cp config/tool_conf.xml.sample config/tool_conf.xml
sed -i -e '/--/d' config/tool_conf.xml	                # un-comment interactive tools
```

The run local instance. On first launch, this will take a while because dependencies are being installed:
```
./run.sh
```
The visit your development Galaxy at http://localhost:8080/ , 
register yourself as user at http://localhost:8080/register/start

Add yourself as admin:
```
echo "  admin_users: sneumann@ipb-halle.de" >>config/galaxy.yml
```

Now you should already be able to launch the Rstudio interactive tool.

### Installing MetFamily into Galaxy

There is a MetFamily Galaxy (interactive) tool under development. 
You can copy the tool.xml and Logo into the Galaxy folder:

```
# Get MetFamily files from GitHub:
wget -O tools/interactive/interactivetool_metfamily.xml https://raw.githubusercontent.com/ipb-halle/MetFamily/refs/heads/feature/galaxify/dev/interactivetool_metfamily.xml
wget -O tools/interactive/MetFamily.png https://raw.githubusercontent.com/ipb-halle/MetFamily/refs/heads/feature/galaxify/inst/MetFamily/www/img/MetFamily.png

sed -i -e 's/askomics/metfamily/' config/tool_conf.xml	# add MetFamily
```
(Alternatively to overwriting the `askomics` IE, you can add an own entry for `interactivetool_metfamily.xml`)

### Using MetFamily

- You will need a saved Project file, e.g. the showcase
from https://raw.githubusercontent.com/ipb-halle/MetFamily/refs/heads/devel/inst/extdata/showcase/Project_file_showcase_annotated.csv.gz

- Upload the project file to Galaxy via its upload tool

- You can launch MetFamily within Galaxy via 
http://localhost:8080/?tool_id=interactive_tool_metfamily&version=latest

- Switch the Input Mode to `Preprocessed MetFamily Project File`, and make sure the project file is selected as input.

- Hit `Run Tool` button, and after the initialisation, you can open MetFamily (ideally using the link out ☐↗) 

- After finishing your analysis, you can export the project. EIther you download from MetFamily, or you export to the Galaxy history.

