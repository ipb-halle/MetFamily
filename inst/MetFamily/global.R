# The global.R file is sourced before the app runs, so it’s a good place for library loading.
# Objects created in global.R are loaded in the global environment, and will remain after the app is closed.
library(MetFamily)

## running as Galaxy Interactive Environment ?
## This variable is either set directly by Galaxy, 
## and/or written to /usr/local/lib/R/etc/Renviron.site
## by the interactivetool_metfamily.xml tool wrapper

isGalaxyIE <- !is.na(Sys.getenv("_GALAXY_JOB_HOME_DIR", unset = NA))
