
load_metfamily_dependencies <- function() 
{
  ##############################################################################################################
  ## GUI
  library("shiny")
  library("htmltools")
  library("shinyjs")
  library("DT")
  library("colourpicker")
  library("shinyBS")
  library("shinybusy")
  library(egg)
  
  ##############################################################################################################
  ## mass spectrometry
  library("mzR")
  library("xcms")
  
  ##############################################################################################################
  ## MS1 analyses
  library("FactoMineR")
  library("mixOmics")
  library("pcaMethods")
  library(searchable)
  library(gdata)
  ##############################################################################################################
  ## tools
  library("matrixStats")
  library("Matrix")
  library("tools")
  library("stringr")
  library("slam")
  
  ##############################################################################################################
  ## pdf report
  library("knitr")
  
  ##############################################################################################################
  ## plot
  library("cba")
  library("squash")
  library("plotrix")
  library("plotly")
  library("RColorBrewer")
}
