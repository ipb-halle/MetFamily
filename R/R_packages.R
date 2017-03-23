
##############################################################################################################
## GUI
#install.packages("shiny")
library("shiny")
#devtools::install_github("rstudio/htmltools")
library("htmltools")
#install.packages("shinyjs")
library("shinyjs")
#install.packages("DT")
library("DT")
#install.packages("colourpicker")
#library("colourpicker") ## TODO
#install.packages("shinyBS")
library("shinyBS")

##############################################################################################################
## mass spectrometry
#source("https://bioconductor.org/biocLite.R")
#biocLite("mzR")
library("mzR")
#biocLite("xmcs")
library("xcms")


##############################################################################################################
## MS1 analyses
#install.packages("FactoMineR")
library("FactoMineR")
#install.packages("mixOmics")
library("mixOmics")
#source("http://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")
library("pcaMethods")

##############################################################################################################
## tools
#install.packages("matrixStats")
library("matrixStats")
library("Matrix")
library("tools")
#install.packages("stringi")
library("stringr")
#install_github('rCharts', 'ramnathv')
#library("rCharts")

##############################################################################################################
## plot
#install.packages("cba")
library("cba")
#install.packages("squash")
library("squash")
#install.packages("plotrix")
library("plotrix")
