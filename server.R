
#########################################################################################
#########################################################################################
## libraries
#install.packages("shiny")
library(shiny)
#devtools::install_github("rstudio/htmltools")
library(htmltools)
#install.packages("shinyjs")
library(shinyjs)

#library(devtools)
#install_github("shinyTable", "trestletech")
#library(shinyTable)
#install.packages("DT")
library(DT)
library("Matrix")

#library("lattice")
#library("stats")
#library("gridSVG")

source("ClusteringMS2SpectraGUI.R")
source("FragmentMatrixFunctions.R")
#source("/vol/R/shiny/srv/shiny-server/MetFam/ClusteringMS2SpectraGUI.R")
#source("/vol/R/shiny/srv/shiny-server/MetFam/FragmentMatrixFunctions.R")

#########################################################################################
#########################################################################################
## global variables
shinyAppFolder <- "/vol/R/shiny/srv/shiny-server/MetFam/"

shinyServer(
  func = function(input, output, session) {
    #########################################################################################
    #########################################################################################
    ## global variables per user
    
    ##############################################
    ## constants
    
    ## MetFamily properties
    toolName    <- "MetFamily"
    toolVersion <- "1.0"
    
    ## data import
    proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- 0.9
    mzDeviationAbsolute_mapping <- 0.01
    minimumNumberOfMS2PeaksPerGroup <- 1
    ## annotation
    artifact <- "Ignore"
    artifactColor <- "red"
    ## HCA
    minimumNumberOfPrecursorsForHca <- 6
    maximumNumberOfPrecursorsForHca <- 1000000
    ## selections
    selectionAnalysisName <- "Selection by analysis"
    selectionFragmentName <- "Selection by fragment"
    selectionSearchName   <- "Selection by search"
    precursorSelectionTabSelection  <- "Selection"
    precursorSelectionTabAnnotation <- "Annotation"
    precursorSelectionTabTable      <- "Table"
    selectionAnalysisHcaName <- "Analysis_HCA"
    selectionAnalysisPcaName <- "Analysis_PCA"
    selectionFragmentHcaName <- "Fragment_HCA"
    selectionFragmentPcaName <- "Fragment_PCA"
    selectionSearchHcaName   <- "Search_HCA"
    selectionSearchPcaName   <- "Search_PCA"
    ## GUI
    runRightColumnWidthFull <- 12
    legendColumnWidthFull <- 2
    runRightColumnWidthPart <- 8
    legendColumnWidthPart <- 2
    annoLegendEntryHeight <- 20
    scoresGroupsLegendEntryHeight <- 20
    maximumNumberOfTableEntries <- 50
    
    ##############################################
    ## data
    dataList <- NULL
    currentDistanceMatrixObj <- NULL
    clusterDataList <- NULL
    pcaDataList <- NULL
    ## filter
    filterGlobal <- NULL
    filterHca <- NULL
    filterPca <- NULL
    filterSearch <- NULL
    
    ##############################################
    ## program state
    state <- reactiveValues(
      initialGuiUpdatePerformed = FALSE, 
      importedOrLoadedFile_s_ = NULL, 
      globalMS2filterValid = FALSE, 
      hcafilterValid = FALSE, 
      pcafilterValid = FALSE, 
      searchfilterValid = FALSE, 
      filterSearchActive = FALSE, 
      runRightColumnWidth = runRightColumnWidthPart, 
      legendColumnWidth = legendColumnWidthPart,
      showSideBar = TRUE, 
      analysisType = "HCA",
      anyPlotDrawn = FALSE,
      showHCAplotPanel = FALSE, 
      showPCAplotPanel = FALSE, 
      plotHcaShown = FALSE,
      plotPcaShown = FALSE,
      precursorSetSelected = FALSE,
      selectedSelection = NULL,
      showClusterLabels = TRUE,
      showScoresLabels = TRUE,
      showLoadingsLabels = FALSE,
      showLoadingsAbundance = FALSE,
      #annotationLegendHeight = annoLegendEntryHeight * (2 + 1 + 1)
      annotationLegendHeightHca = -1,
      annotationLegendHeightPca = -1,
      ## plot annotations
      annotationsHca = NULL,
      annotationsPca = NULL,
      scoresGroups = NULL,
      scoresGroupsLegendHeight = -1
    )
    changeSelectionCurrentSelection <- selectionAnalysisName
    precursorSelectionTabCurrentTab <- precursorSelectionTabSelection
    plotToShow <- "Display HCA"
    showSideBar <- TRUE
    
    ##############################################
    ## buttons
    updateProjectDescriptionButtonValue <- 0
    applyGlobalMS2filtersButtonValue <- 0
    applySearchButtonValue <- 0
    clearSearchButtonValue <- 0
    applyHcaFiltersButtonValue <- 0
    applyPcaFiltersButtonValue <- 0
    drawHCAButtonValue <- 0
    drawPCAButtonValue <- 0
    downloadMatrixButtonValue <- 0
    removePresentAnnotationValue <- 0
    setPresentAnnotationPrimaryValue <- 0
    submitNewAnnotationValue <- 0
    submitPreviousAnnotationValue <- 0
    updateArtifactsFromCheckboxesButtonValue <- 0
    clearSelectionButtonValue <- 0
    importMs1Ms2DataButtonValue <- 0
    loadProjectDataButtonValue <- 0
    loadExampleDataButtonValue <- 0
    
    ##############################################
    ## MS2 peaks
    fragmentsX <- NULL
    fragmentsY <- NULL
    fragmentsColor <- NULL
    fragmentsDiscriminativity <- NULL
    fragmentsXhovered <- NULL
    fragmentsYhovered <- NULL
    fragmentsColorHovered <- NULL
    fragmentsDiscriminativityHovered <- NULL
    
    ##############################################
    ## plot ranges
    dendrogramPlotRangeY <- NULL
    dendrogramPlotRange  <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
    fragmentPlotRange    <- reactiveValues(xMin = NULL, xMax = NULL)
    ms2PlotRange         <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
    pcaScoresPlotRange   <- reactiveValues(
      xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL, 
      yMin = NULL, yMax = NULL, yInterval = NULL, yIntervalSize = NULL
    )
    pcaLoadingsPlotRange <- reactiveValues(
      xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL, 
      yMin = NULL, yMax = NULL, yInterval = NULL, yIntervalSize = NULL
    )
    
    ##############################################
    ## selections
    
    ## selection MS2
    selectionFragmentSelectedFragmentIndex <- NULL
    selectionFragmentTreeNodeSet <- NULL
    selectionFragmentPcaLoadingSet <- NULL
    ## selection analysis
    selectionAnalysisTreeNodeSet <- NULL
    selectionAnalysisPcaLoadingSet <- NULL
    ## selection search
    selectionSearchTreeNodeSet <- NULL
    selectionSearchPcaLoadingSet <- NULL
    
    ## table data
    tableCheckboxIdCounter <- 0
    selectedPrecursorSet <- NULL
    
    selectedTable <- NULL
    selectedTable_id <- NULL
    table <- reactiveValues(
      df_Fragment_HCA = NULL,
      df_Fragment_PCA = NULL,
      df_Search_HCA = NULL,
      df_Search_PCA = NULL,
      df_Analysis_HCA = NULL,
      df_Analysis_PCA = NULL
    )
    listForTable_Fragment_HCA = NULL
    listForTable_Fragment_PCA = NULL
    listForTable_Analysis_HCA = NULL
    listForTable_Analysis_PCA = NULL
    listForTable_Search_HCA = NULL
    listForTable_Search_PCA = NULL
    
    table_Fragment_HCA_id = NULL
    table_Fragment_PCA_id = NULL
    table_Analysis_HCA_id = NULL
    table_Analysis_PCA_id = NULL
    table_Search_HCA_id = NULL
    table_Search_PCA_id = NULL
    
    #########################################################################################
    #########################################################################################
    ## functions
    
    ## POI selection
    getSelectedPOI_X <- function(mouseX, poiCoordinatesX, plotWidth, plotRangeX){
      if(any(is.na(c(poiCoordinatesX))))
        return(NULL)
      
      factorX <- plotWidth  / plotRangeX
      
      mouseX <- mouseX * factorX
      poiCoordinatesX <- poiCoordinatesX * factorX
      
      distances <- abs(poiCoordinatesX - mouseX)
      distanceThreshold <- factorX * plotRangeX / 35
      
      minimumIndex <- which.min(distances)
      minimumDistance <- distances[[minimumIndex]]
      
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    getSelectedPOI_XY <- function(mouseX, mouseY, poiCoordinatesX, poiCoordinatesY, plotWidth, plotHeight, plotRangeX, plotRangeY){
      if(any(is.na(c(poiCoordinatesX, poiCoordinatesY))))
        return(NULL)
      
      factorX <- plotWidth  / plotRangeX
      factorY <- plotHeight / plotRangeY
      
      mouseX <- mouseX * factorX
      mouseY <- mouseY * factorY
      poiCoordinatesX <- poiCoordinatesX * factorX
      poiCoordinatesY <- poiCoordinatesY * factorY
      
      distancesX <- poiCoordinatesX - mouseX
      distancesY <- poiCoordinatesY - mouseY
      distances <- sqrt(distancesX * distancesX + distancesY * distancesY)
      distanceThreshold <- factorX * plotRangeX / 35
      
      minimumIndex <- which.min(distances)
      if(is.na(minimumIndex))
        return(NULL)
      minimumDistance <- distances[[minimumIndex]]
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    ## Parse the input file
    resetWorkspace <- function(){
      print(paste("resetWorkspace"))
      
      #########################################################################################
      ## reset
      
      ## reset plots
      doClearPlots()
      
      ## reset variables
      clusterDataList <<- NULL
      pcaDataList <<- NULL
      
      ## selection
      selectionFragmentSelectedFragmentIndex <- NULL
      
      selectionFragmentTreeNodeSet <- NULL
      selectionAnalysisTreeNodeSet <- NULL
      selectionSearchTreeNodeSet <- NULL
      
      selectionFragmentPcaLoadingSet <- NULL
      selectionAnalysisPcaLoadingSet <- NULL
      selectionSearchPcaLoadingSet <- NULL
      
      selectionByFragmentReset()
      selectionByAnalysisReset()
      selectionBySearchReset()
      
      
      ## fragments
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      fragmentsColor <<- NULL
      fragmentsDiscriminativity <<- NULL
      fragmentsXhovered <<- NULL
      fragmentsYhovered <<- NULL
      fragmentsColorHovered <<- NULL
      fragmentsDiscriminativityHovered <<- NULL
      
      ## reset state
      state$importedOrLoadedFile_s_ <<- NULL
      state$analysisType <<- "HCA"
      state$showHCAplotPanel <<- FALSE
      state$showPCAplotPanel <<- FALSE
      state$plotHcaShown <<- FALSE
      state$plotPcaShown <<- FALSE
      state$precursorSetSelected <<- FALSE
      state$anyPlotDrawn <<- FALSE
      selectedPrecursorSet <<- NULL
      selectedTable <<- NULL
      state$annotationsHca <<- NULL
      state$annotationsPca <<- NULL
      state$scoresGroups <<- NULL
      
      ## reset plot range
      dendrogramPlotRangeY <<- NULL
      dendrogramPlotRange$xMin <<- NULL
      dendrogramPlotRange$xMax <<- NULL
      dendrogramPlotRange$xInterval <<- NULL
      dendrogramPlotRange$xIntervalSize <<- NULL
      pcaScoresPlotRange$xMin <<- NULL
      pcaScoresPlotRange$xMax <<- NULL
      pcaScoresPlotRange$xInterval <<- NULL
      pcaScoresPlotRange$xIntervalSize <<- NULL
      pcaScoresPlotRange$yMin <<- NULL
      pcaScoresPlotRange$yMax <<- NULL
      pcaScoresPlotRange$yInterval <<- NULL
      pcaScoresPlotRange$yIntervalSize <<- NULL
      pcaLoadingsPlotRange$xMin <<- NULL
      pcaLoadingsPlotRange$xMax <<- NULL
      pcaLoadingsPlotRange$xInterval <<- NULL
      pcaLoadingsPlotRange$xIntervalSize <<- NULL
      pcaLoadingsPlotRange$yMin <<- NULL
      pcaLoadingsPlotRange$yMax <<- NULL
      pcaLoadingsPlotRange$yInterval <<- NULL
      pcaLoadingsPlotRange$yIntervalSize <<- NULL
      
      #########################################################################################
      ## update fragment plot
      
      min <- min(dataList$masses)
      max <- max(dataList$masses)
      
      fragmentPlotRange$xMin <<- min
      fragmentPlotRange$xMax <<- max
      fragmentPlotRange$xInterval <<- c(min, max)
      fragmentPlotRange$xIntervalSize <<- max - min
      
      output$fragmentPlot <- renderPlot({
        print(paste("### MS2 ### all"))
        plotFragments(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
      #}, bg = "transparent")
      })
      
      #########################################################################################
      ## update filter input values
      
      ## groups
      switch(as.character(length(dataList$groups)), 
        "0"={
          stop("No groups available")
        },
        "1"={
          selectedOne <- dataList$groups[[1]]
          selectedTwo <- dataList$groups[[1]]
        },
        {
          selectedOne <- dataList$groups[[1]]
          selectedTwo <- dataList$groups[[2]]
        }
      )
      updateRadioButtons(session = session, inputId = "hcaFilterGroupOne", choices = dataList$groups, selected = selectedOne)
      updateRadioButtons(session = session, inputId = "hcaFilterGroupTwo", choices = dataList$groups, selected = selectedTwo)
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = dataList$groups)
      
      ## input fields: HCA
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_ppm", value = "20")
      
      ## input fields: HCA
      updateTextInput(session = session, inputId = "hcaFilter_average", value = "0")
      updateTextInput(session = session, inputId = "hcaFilter_lfc", value = "0")
      updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
      ## input fields: PCA
      updateTextInput(session = session, inputId = "pcaFilter_average", value = "0")
      updateTextInput(session = session, inputId = "pcaFilter_lfc", value = "0")
      updateCheckboxInput(session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", value = FALSE)
      ## input fields: search MS1
      updateTextInput(session = session, inputId = "searchMS1mass", value = "")
      updateTextInput(session = session, inputId = "searchMS1massPpm", value = 20)
      ## input fields: search MS2
      updateTextInput(session = session, inputId = "search_ms2_masses1", value = "")
      updateTextInput(session = session, inputId = "search_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "search_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "searchMS2massPpm", value = 20)
      updateCheckboxInput(session = session, inputId = "searchIncludeIgnoredPrecursors", value = FALSE)
      #updateColourInput(session = session, inputId = "newAnnotationColor", allowedCols = colorPalette())
      ## anno
      updateTextInput(session = session, inputId = "newAnnotationValue", value = "")
      updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[init]"))
      updateSelectInput(session = session, inputId = "previousAnnotationValue", choices = c("Artifact"))
      
      #########################################################################################
      ## update filter
      filter <- doPerformFiltering(dataList$groups, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      if(length(dataList$groups) == 1)
        filter2 <- doPerformFiltering(c(dataList$groups[[1]], dataList$groups[[1]]), NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      else
        filter2 <- filter
      
      filterGlobal <<- filter
      filterHca    <<- filter2
      filterPca    <<- filter
      state$filterSearchActive <<- FALSE
      state$searchfilterValid <<- TRUE
      filterSearch    <<- NULL
      
      updateGlobalMS2filterInformation()
      updateHcaFilterInformation()
      updatePcaFilterInformation()
      updateSearchInformation()
      
      state$globalMS2filterValid <<- TRUE
      state$hcaFilterValid <<- TRUE
      state$pcaFilterValid <<- TRUE
      
      checkHcaFilterValidity(filter$numberOfPrecursorsFiltered)
      checkPcaFilterValidity(filter$numberOfPrecursorsFiltered)
      
      state$plotHcaShown <<- FALSE
      state$plotPcaShown <<- FALSE
      
      ## project infos
      updateTextInput    (session = session, inputId = "projectName2",           value = dataList$importParameterSet$projectName)
      shinyjs::toggleState("projectName2", FALSE)
      updateTextInput    (session = session, inputId = "projectDescription2",    value = dataList$importParameterSet$projectDescription)
      #shinyjs::toggleState("projectName2", FALSE)
      
      updateTextInput    (session = session, inputId = "minimumIntensityOfMaximalMS2peak2",            value = dataList$importParameterSet$minimumIntensityOfMaximalMS2peak)
      updateTextInput    (session = session, inputId = "minimumProportionOfMS2peaks2",                 value = dataList$importParameterSet$minimumProportionOfMS2peaks)
      updateTextInput    (session = session, inputId = "mzDeviationAbsolute_grouping2",                value = dataList$importParameterSet$mzDeviationAbsolute_grouping)
      updateTextInput    (session = session, inputId = "mzDeviationInPPM_grouping2",                   value = dataList$importParameterSet$mzDeviationInPPM_grouping)
      updateCheckboxInput(session = session, inputId = "doPrecursorDeisotoping2",                      value = dataList$importParameterSet$doPrecursorDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationAbsolute_precursorDeisotoping2",    value = dataList$importParameterSet$mzDeviationAbsolute_precursorDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationInPPM_precursorDeisotoping2",       value = dataList$importParameterSet$mzDeviationInPPM_precursorDeisotoping)
      updateTextInput    (session = session, inputId = "maximumRtDifference2",                         value = dataList$importParameterSet$maximumRtDifference)
      updateCheckboxInput(session = session, inputId = "doMs2PeakGroupDeisotoping2",                   value = dataList$importParameterSet$doMs2PeakGroupDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping2", value = dataList$importParameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping2",    value = dataList$importParameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping)
      #dataList$importParameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
      #dataList$importParameterSet$mzDeviationAbsolute_mapping
      #dataList$importParameterSet$minimumNumberOfMS2PeaksPerGroup
      updateCheckboxInput(session = session, inputId = "neutralLossesPrecursorToFragments2",           value = dataList$importParameterSet$neutralLossesPrecursorToFragments)
      updateCheckboxInput(session = session, inputId = "neutralLossesFragmentsToFragments2",           value = dataList$importParameterSet$neutralLossesFragmentsToFragments)
      
      shinyjs::toggleState("minimumIntensityOfMaximalMS2peak2", FALSE)
      shinyjs::toggleState("minimumProportionOfMS2peaks2", FALSE)
      shinyjs::toggleState("mzDeviationAbsolute_grouping2", FALSE)
      shinyjs::toggleState("mzDeviationInPPM_grouping2", FALSE)
      shinyjs::toggleState("doPrecursorDeisotoping2", FALSE)
      shinyjs::toggleState("mzDeviationAbsolute_precursorDeisotoping2", FALSE)
      shinyjs::toggleState("mzDeviationInPPM_precursorDeisotoping2", FALSE)
      shinyjs::toggleState("maximumRtDifference2", FALSE)
      shinyjs::toggleState("doMs2PeakGroupDeisotoping2", FALSE)
      shinyjs::toggleState("mzDeviationAbsolute_ms2PeakGroupDeisotoping2", FALSE)
      shinyjs::toggleState("mzDeviationInPPM_ms2PeakGroupDeisotoping2", FALSE)
      shinyjs::toggleState("neutralLossesPrecursorToFragments2", FALSE)
      shinyjs::toggleState("neutralLossesFragmentsToFragments2", FALSE)
      
      ## MS2 plot range
      resetMS2PlotRange()
    }
    ## filter info
    updateGlobalMS2filterInformation <- function(){
      if(is.null(filterGlobal)){
        output$globalMS2filteredPrecursors <- renderText({
          print(paste("update output$globalMS2filteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
      } else {
        output$globalMS2filteredPrecursors <- renderText({
          print(paste("update output$globalMS2filteredPrecursors", sep = ""))
          paste("Number of filtered MS\u00B9 features: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, sep = "")
        })
      }
    }
    updateHcaFilterInformation <- function(){
      if(is.null(filterHca)){
        ## errors
        output$hcaFilteredPrecursors <- renderText({
          print(paste("update output$hcaFilteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values.", sep = "")
        })
      } else {
        ## no errors
        globalMs2Filter <- ifelse(filterGlobal$numberOfPrecursorsFiltered != dataList$numberOfPrecursors, paste("\n(Global MS/MS filter: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, " precursors)", sep = ""), "")
        
        if(filterHca$numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & filterHca$numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
          ## filter valid
          output$hcaFilteredPrecursors <- renderText({
            print(paste("update output$hcaFilteredPrecursors ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
            paste("Number of filtered MS\u00B9 features: ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, globalMs2Filter, sep = "")
          })
        } else {
          ## filter invalid
          
          ## update info
          if(filterHca$numberOfPrecursorsFiltered == 0){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors # = 0", sep = ""))
              paste("There are no MS\u00B9 features which fulfill the given criteria.", globalMs2Filter, sep = "")
            })
          }
          if(filterHca$numberOfPrecursorsFiltered > 0 & filterHca$numberOfPrecursorsFiltered < minimumNumberOfPrecursorsForHca){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors 0 < # < ", minimumNumberOfPrecursorsForHca, sep = ""))
              paste("There are only ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, " MS\u00B9 features which fulfill the given criteria. There must be at least more than five MS\u00B9 features to proceed.", globalMs2Filter, sep = "")
            })
          }
          if(filterHca$numberOfPrecursorsFiltered > maximumNumberOfPrecursorsForHca){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors # > ", maximumNumberOfPrecursorsForHca, sep = ""))
              paste("There are ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, " MS\u00B9 features which fulfill the given criteria. There must be at most ", maximumNumberOfPrecursorsForHca, " MS\u00B9 features to proceed.", globalMs2Filter, sep = "")
            })
          }
        }
      }
    }
    updatePcaFilterInformation <- function(){
      if(is.null(filterPca)){
        output$pcaFilteredPrecursors <- renderText({
          print(paste("update output$pcaFilteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values.", sep = "")
        })
      } else {
        globalMs2Filter <- ifelse(filterGlobal$numberOfPrecursorsFiltered != dataList$numberOfPrecursors, paste("\n(Global MS/MS filter: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, " precursors)", sep = ""), "")
        if(filterPca$numberOfPrecursorsFiltered > 0){
          output$pcaFilteredPrecursors <- renderText({
            print(paste("update output$pcaFilteredPrecursors", sep = ""))
            paste("Number of filtered MS\u00B9 features: ", filterPca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, globalMs2Filter, sep = "")
          })
        } else {
          output$pcaFilteredPrecursors <- renderText({
            print(paste("update output$pcaFilteredPrecursors", sep = ""))
            paste("There are no MS\u00B9 features which fulfill the given criteria.", globalMs2Filter, sep = "")
          })
        }
      }
    }
    updateSearchInformation <- function(){
      if(!state$filterSearchActive)
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo inactive search", sep = ""))
          paste("Please search for MS\u00B9 features", sep = "")
        })
      if(state$filterSearchActive & is.null(filterSearch))
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo invalid search", sep = ""))
          paste("There are invalid or missing search values", sep = "")
        })
      if(state$filterSearchActive & !is.null(filterSearch)){
        
        str1 <- ""
        str2 <- ""
        if(state$showHCAplotPanel & !is.null(listForTable_Search_HCA))
          str1 <- paste(length(listForTable_Search_HCA$precursorSet), " in HCA", sep = "")
        if(state$showPCAplotPanel & !is.null(listForTable_Search_PCA))
          str2 <- paste(length(listForTable_Search_PCA$precursorSet), " in PCA", sep = "")
        
        if(nchar(str1) > 0 & nchar(str2) > 0)
          val <- paste(str1, str2, sep = ", ")
        if(nchar(str1) > 0 & !(nchar(str2) > 0))
          val <- str1
        if(!(nchar(str1) > 0) & nchar(str2) > 0)
          val <- str2
        if(!(nchar(str1) > 0) & !(nchar(str2) > 0))
          val <- "None"
        
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo", sep = ""))
          paste("Number of hits among MS\u00B9 features: ", val, sep = "")
        })
      }
    }
    ## perform filtering
    doPerformFiltering <- function(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, preFilter = NULL){
      print(paste("Observe applyFilters1", "gs", paste(groupSet, collapse = "-"), "a", filter_average, "lfc", filter_lfc, "ms2_1", filter_ms2_masses1, "ms2_2", filter_ms2_masses2, "ms2_3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "ig", includeIgnoredPrecursors))
      
      groupSetOriginal                 <- groupSet
      filter_averageOriginal           <- filter_average
      filter_lfcOriginal               <- filter_lfc
      filter_ms2_masses1Original       <- filter_ms2_masses1
      filter_ms2_masses2Original       <- filter_ms2_masses2
      filter_ms2_masses3Original       <- filter_ms2_masses3
      filter_ms2_ppmOriginal           <- filter_ms2_ppm
      filter_ms1_massesOriginal        <- filter_ms1_masses
      filter_ms1_ppmOriginal           <- filter_ms1_ppm
      includeIgnoredPrecursorsOriginal <- includeIgnoredPrecursors
      
      #################################################
      ## parse inputs
      print(paste("Observe applyFilters2", "gs", paste(groupSet, collapse = "-"), "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "i", includeIgnoredPrecursors))
      print(paste("Observe applyFilters3", "gs", is.null(groupSet), "a", is.null(filter_average), "lfc", is.null(filter_lfc), "ms2 1", is.null(filter_ms2_masses1), "ms2 2", is.null(filter_ms2_masses2), "ms2 3", is.null(filter_ms2_masses3), "ppm", is.null(filter_ms2_ppm), "i", includeIgnoredPrecursors))
      
      #################################################
      ## sanity checks
      if(!is.null(filter_lfc) & length(groupSet) != 2)
        stop("lfc filter for not exactly two groups")
      
      #################################################
      ## check for errors in inputs amd process ms2
      error <- FALSE
      if(any(is.null(groupSet), length(groupSet) == 0, nchar(groupSet) == 0))
        error <- TRUE
      
      if(any(is.null(filter_average), length(filter_average) == 0, nchar(filter_average) == 0))
        filter_average <- NULL
      else{
        filter_average <- as.numeric(filter_average)
        error <- error | is.na(filter_average)
      }
      
      if(any(is.null(filter_lfc), length(filter_lfc) == 0, nchar(filter_lfc) == 0))
        filter_lfc <- NULL
      else{
        filter_lfc <- as.numeric(filter_lfc)
        error <- error | is.na(filter_lfc)
      }
      
      if(any(is.null(filter_ms2_masses1), length(filter_ms2_masses1) == 0, nchar(filter_ms2_masses1) == 0))
        filter_ms2_masses1 <- NULL
      else{
        ms2Masses <- strsplit(x = filter_ms2_masses1, split = "[,; ]+")[[1]]
        filter_ms2_masses1 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses1[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses1))
      }
      if(any(is.null(filter_ms2_masses2), length(filter_ms2_masses2) == 0, nchar(filter_ms2_masses2) == 0))
        filter_ms2_masses2 <- NULL
      else{
        ms2Masses <- strsplit(x = filter_ms2_masses2, split = "[,; ]+")[[1]]
        filter_ms2_masses2 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses2[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses2))
      }
      if(any(is.null(filter_ms2_masses3), length(filter_ms2_masses3) == 0, nchar(filter_ms2_masses3) == 0))
        filter_ms2_masses3 <- NULL
      else{
        ms2Masses <- strsplit(x = filter_ms2_masses3, split = "[,; ]+")[[1]]
        filter_ms2_masses3 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses3[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses3))
      }
      
      if(any(is.null(filter_ms2_ppm), length(filter_ms2_ppm) == 0, nchar(filter_ms2_ppm) == 0))
        filter_ms2_ppm <- NULL
      else{
        filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
        error <- error | is.na(filter_ms2_ppm)
      }
      
      if(any(is.null(filter_ms1_masses), length(filter_ms1_masses) == 0, nchar(filter_ms1_masses) == 0))
        filter_ms1_masses <- NULL
      else{
        ms1Masses <- strsplit(x = filter_ms1_masses, split = "[,; ]+")[[1]]
        filter_ms1_masses <- vector(mode = "numeric", length = length(ms1Masses))
        for(idx in 1:length(ms1Masses))
          filter_ms1_masses[[idx]] <- as.numeric(ms1Masses[[idx]])
        error <- error | any(is.na(filter_ms1_masses))
      }
      
      if(any(is.null(filter_ms1_ppm), length(filter_ms1_ppm) == 0, nchar(filter_ms1_ppm) == 0))
        filter_ms1_ppm <- NULL
      else{
        filter_ms1_ppm <- as.numeric(filter_ms1_ppm)
        error <- error | is.na(filter_ms1_ppm)
      }
      
      ## sanity check
      error <- error | (!is.null(filter_ms1_masses) & is.null(filter_ms1_ppm))
      error <- error | ((!is.null(filter_ms2_masses1) | !is.null(filter_ms2_masses2) | !is.null(filter_ms2_masses3)) & is.null(filter_ms2_ppm))
      
      print(paste("Observe applyFilters4", "e", error, "gs", paste(groupSet, collapse = "-"), "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "i", includeIgnoredPrecursors))
      
      filterList_ms2_masses <- list()
      if(!is.null(filter_ms2_masses1))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses1
      if(!is.null(filter_ms2_masses2))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses2
      if(!is.null(filter_ms2_masses3))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses3
      
      #################################################
      ## do filtering and wrap results
      resultObj <- list()
      resultObj$error  <- error
      
      if(error){
        resultObj$filter <- NULL
      } else {
        filterHere <- filterData(
          dataList = dataList, 
          #groupOne = groupOne, groupTwo = groupTwo, 
          groups = groupSet, filter_average = filter_average, filter_lfc = filter_lfc, 
          filterList_ms2_masses = filterList_ms2_masses, filter_ms2_ppm = filter_ms2_ppm, 
          filter_ms1_masses = filter_ms1_masses, filter_ms1_ppm = filter_ms1_ppm,
          includeIgnoredPrecursors = includeIgnoredPrecursors,
          progress = FALSE
        )
        
        ## set original values
        filterHere$groupSetOriginal                 <- groupSetOriginal                
        filterHere$filter_averageOriginal           <- filter_averageOriginal          
        filterHere$filter_lfcOriginal               <- filter_lfcOriginal              
        filterHere$filter_ms2_masses1Original       <- filter_ms2_masses1Original      
        filterHere$filter_ms2_masses2Original       <- filter_ms2_masses2Original      
        filterHere$filter_ms2_masses3Original       <- filter_ms2_masses3Original      
        filterHere$filter_ms2_ppmOriginal           <- filter_ms2_ppmOriginal          
        filterHere$filter_ms1_massesOriginal        <- filter_ms1_massesOriginal         
        filterHere$filter_ms1_ppmOriginal           <- filter_ms1_ppmOriginal          
        filterHere$includeIgnoredPrecursorsOriginal <- includeIgnoredPrecursorsOriginal
        
        resultObj$error  <- error
        resultObj$filter <- filterHere
        print(paste("Observe applyFilters5", "n", resultObj$filter$numberOfPrecursorsFiltered))
      }
      
      return(resultObj)
    }
    processSearchFilterResult <- function(resultObj){
      state$filterSearchActive <<- TRUE
      #################################################
      ## info / error
      if(resultObj$error){
        filterSearch <<- NULL
        state$searchfilterValid <<- FALSE
        selectionBySearch(NULL)
        updateSearchInformation()
      } else {
        filterSearch <<- resultObj$filter
        state$searchfilterValid <<- TRUE
        selectionBySearch(filterSearch$filter)
        updateSearchInformation()
        
        #################################################
        ## update plots
        if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "update search")
        if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "update search")
      }
    }
    checkHcaFilterValidity <- function(numberOfPrecursorsFiltered){
      if(numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
        ## filter valid
        print(paste("Observe applyFilters ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
        
        shinyjs::enable("drawHCAplots")
        #enableActionButton(session, "drawHCAplots")
        state$hcaFilterValid <<- TRUE
      } else {
        ## filter invalid
        
        shinyjs::disable("drawHCAplots")
        #disableActionButton(session, "drawHCAplots")
        state$hcaFilterValid <<- FALSE
      }
    }
    checkPcaFilterValidity <- function(numberOfPrecursorsFiltered){
      if(numberOfPrecursorsFiltered > 0){
        ## filter valid
        print(paste("Observe applyFilters # > 0", sep = ""))
        
        shinyjs::enable("drawPCAplots")
        state$pcaFilterValid <<- TRUE
      } else {
        ## filter invalid
        print(paste("Observe applyFilters # = 0", sep = ""))
        
        shinyjs::disable("drawPCAplots")
        state$pcaFilterValid <<- FALSE
      }
    }
    ## input
    updateFileInputInfo <- function(){
      fileInputSelection <- input$fileInputSelection
      filePath <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      fileMs1Path <- input$ms1DataFile$datapath
      fileMs1Name <- input$ms1DataFile$name
      fileMs2Path <- input$ms2DataFile$datapath
      fileMs2Name <- input$ms2DataFile$name
      #exampleDataSelection <- input$exampleDataSelection
      
      if(all(fileInputSelection == "Example data"))
        output$fileInfo <- renderText({paste("Please press 'Load example data' to load the full example data set")})
      #if(all(fileInputSelection == "Example data", exampleDataSelection == "Example data set (full)"))
      #  output$fileInfo <- renderText({paste("Please press 'Load example data' to load the full example data set")})
      #if(all(fileInputSelection == "Example data", exampleDataSelection == "Example data set (reduced)"))
      #  output$fileInfo <- renderText({paste("Please press 'Load example data' to load the reduced example data set")})
      if(all(fileInputSelection == "Load project", is.null(filePath)))
        output$fileInfo <- renderText({paste("Please select a project file and press 'Load project data'")})
      if(all(fileInputSelection == "Load project", !is.null(filePath), any(is.null(state$importedOrLoadedFile_s_), fileName != state$importedOrLoadedFile_s_)))
        output$fileInfo <- renderText({paste("Please press 'Load project data'")})
      if(all(fileInputSelection == "Load project", !is.null(filePath), !is.null(state$importedOrLoadedFile_s_), fileName == state$importedOrLoadedFile_s_))
        output$fileInfo <- renderText({paste(fileName)})
      if(all(fileInputSelection == "Import data", is.null(fileMs1Path), is.null(fileMs2Path)))
        output$fileInfo <- renderText({paste("Please select a metabolite profile, a MS/MS library, and press 'Import MS\u00B9 and MS/MS data'")})
      if(all(fileInputSelection == "Import data", is.null(fileMs1Path), !is.null(fileMs2Path)))
        output$fileInfo <- renderText({paste("Please select a metabolite profile and press 'Import MS\u00B9 and MS/MS data'")})
      if(all(fileInputSelection == "Import data", !is.null(fileMs1Path), is.null(fileMs2Path)))
        output$fileInfo <- renderText({paste("Please select a MS/MS library and press 'Import MS\u00B9 and MS/MS data'")})
      if(all(fileInputSelection == "Import data", !is.null(fileMs1Path), !is.null(fileMs2Path), any(is.null(state$importedOrLoadedFile_s_), c(fileMs1Name,fileMs2Name) != state$importedOrLoadedFile_s_)))
        output$fileInfo <- renderText({paste("Please press 'Import MS\u00B9 and MS/MS data'")})
      if(all(fileInputSelection == "Import data", !is.null(fileMs1Path), !is.null(fileMs2Path), !is.null(state$importedOrLoadedFile_s_), c(fileMs1Name,fileMs2Name) == state$importedOrLoadedFile_s_))
        output$fileInfo <- renderText({paste(fileMs1Name, "\n", fileMs2Name, sep = "")})
    }
    ## plots
    doClearPlots <- function(){
      output$plotDendrogram <- renderPlot({
        print(paste("reset output$plotDendrogram"))
        NULL
      })
      output$plotHeatmap <- renderPlot({
        print(paste("reset output$plotHeatmap"))
        NULL
      })
      output$plotHeatmapLegend <- renderPlot({
        print(paste("reset output$plotHeatmapLegend"))
        NULL
      })
      output$plotMS2 <- renderPlot({
        print(paste("reset output$plotMS2"))
        NULL
      })
      output$plotPcaScores <- renderPlot({
        print(paste("reset output$plotPcaScores"))
        NULL
      })
      output$plotPcaLoadings <- renderPlot({
        print(paste("reset output$plotPcaLoadings"))
        NULL
       })
      output$plotAnnoLegendHCA <- renderPlot({
        print(paste("reset output$plotAnnoLegendHCA"))
        NULL
      })
      output$plotAnnoLegendPCA <- renderPlot({
        print(paste("reset output$plotAnnoLegendPCA"))
        NULL
      })
      output$plotMS2Legend <- renderPlot({
        print(paste("reset output$plotMS2Legend"))
        NULL
      })
      output$plotFragmentDiscriminativityLegend <- renderPlot({
        print(paste("reset output$plotFragmentDiscriminativityLegend"))
        NULL
      })
    }
    drawMS2Plot <- function(consoleInfo = NULL){
      output$plotMS2 <- renderPlot({
        print(paste("### MS2 ###", consoleInfo))
        drawMS2PlotImpl()
      })
    }
    drawMS2PlotImpl <- function(){
      calcPlotMS2(
        dataList = dataList, 
        fragmentsX = fragmentsX, 
        fragmentsY = fragmentsY, 
        fragmentsColor = fragmentsColor, 
        fragmentsDiscriminativity = fragmentsDiscriminativity, 
        fragmentsX_02 = fragmentsXhovered, 
        fragmentsY_02 = fragmentsYhovered, 
        fragmentsColor_02 = fragmentsColorHovered, 
        fragmentsDiscriminativity_02 = fragmentsDiscriminativityHovered, 
        xInterval = ms2PlotRange$xInterval, 
        selectedFragmentIndex = selectionFragmentSelectedFragmentIndex  
      )
    }
    drawDendrogramPlot <- function(consoleInfo = NULL, withHeatmap = FALSE){
      output$plotDendrogram <- renderPlot({
        print(paste("### den ###", consoleInfo))
        drawDendrogramPlotImpl()
      })
      if(!withHeatmap)
        return()
      
      output$plotHeatmap <- renderPlot({
        print(paste("### hea ### update range output$plotHeatmap"))
        drawHeatmapPlotImpl()
      })
    }
    drawDendrogramPlotImpl <- function(){
      resultObj <- calcPlotDendrogram(
        dataList = dataList, 
        filter = filterHca$filter, 
        clusterDataList = clusterDataList, 
        annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, 
        annoPresentColorsList = dataList$annoPresentColorsList, 
        distanceMeasure = currentDistanceMatrixObj$distanceMeasure, 
        selectionFragmentTreeNodeSet = selectionFragmentTreeNodeSet,
        selectionAnalysisTreeNodeSet = selectionAnalysisTreeNodeSet,
        selectionSearchTreeNodeSet = selectionSearchTreeNodeSet,
        showClusterLabels = state$showClusterLabels, 
        xInterval = dendrogramPlotRange$xInterval
      )
      
      dendrogramPlotRange <- par("usr")
      dendrogramPlotRangeY <<- dendrogramPlotRange[[4]] - dendrogramPlotRange[[3]]
      
      state$annotationsHca <<- resultObj
      state$annotationLegendHeightHca <<- annoLegendEntryHeight * (length(state$annotationsHca$setOfAnnotations) + 1)
    }
    drawHeatmapPlotImpl <- function(){
      calcPlotHeatmap(dataList = dataList, filterObj = filterHca, clusterDataList = clusterDataList, xInterval = dendrogramPlotRange$xInterval)
    }
    drawPcaPlots <- function(consoleInfo = NULL){
      drawPcaScoresPlot(consoleInfo = consoleInfo)
      drawPcaLoadingsPlot(consoleInfo = consoleInfo)
    }
    drawPcaScoresPlot <- function(consoleInfo = NULL){
      output$plotPcaScores <- renderPlot({
        print(paste("### psc ###", consoleInfo))
        drawPcaScoresPlotImpl()
      })
    }
    drawPcaScoresPlotImpl <- function(){
      resultObj <- calcPlotPCAscores(
        pcaObj = pcaDataList$pcaObj, 
        dataList = dataList, 
        filterObj = filterPca, 
        pcaDimensionOne = pcaDataList$dimensionOne, 
        pcaDimensionTwo = pcaDataList$dimensionTwo, 
        showScoresLabels = state$showScoresLabels, 
        xInterval = pcaScoresPlotRange$xInterval, 
        yInterval = pcaScoresPlotRange$yInterval
      )
      
      state$scoresGroups <<- resultObj
      state$scoresGroupsLegendHeight <<- scoresGroupsLegendEntryHeight * (length(state$scoresGroups$groups) + 1)
    }
    drawPcaLoadingsPlot <- function(consoleInfo = NULL){
      output$plotPcaLoadings <- renderPlot({
        print(paste("### psc ###", consoleInfo))
        drawPcaLoadingsPlotImpl()
      })
    }
    drawPcaLoadingsPlotImpl <- function(){
      resultObj <- calcPlotPCAloadings(
        pcaObj = pcaDataList$pcaObj, 
        dataList = dataList, 
        filter = filterPca$filter, 
        pcaDimensionOne = pcaDataList$dimensionOne, 
        pcaDimensionTwo = pcaDataList$dimensionTwo, 
        selectionFragmentPcaLoadingSet = selectionFragmentPcaLoadingSet,
        selectionAnalysisPcaLoadingSet = selectionAnalysisPcaLoadingSet,
        selectionSearchPcaLoadingSet = selectionSearchPcaLoadingSet,
        showLoadingsLabels = state$showLoadingsLabels, 
        showLoadingsAbundance = state$showLoadingsAbundance, 
        xInterval = pcaLoadingsPlotRange$xInterval, 
        yInterval = pcaLoadingsPlotRange$yInterval
      )
      
      state$annotationsPca <<- resultObj
      state$annotationLegendHeightPca <<- annoLegendEntryHeight * (length(state$annotationsPca$setOfAnnotations) + 1)
    }
    drawHeatmapLegend <- function(consoleInfo = NULL){
      output$plotHeatmapLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawHeatmapLegendImpl()
      })
    }
    drawHeatmapLegendImpl <- function(){
      calcPlotHeatmapLegend(dataList = dataList)
    }
    drawAnnotationLegendHCA <- function(consoleInfo = NULL){
      output$plotAnnoLegendHCA <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcPlotAnnoLegend(state$annotationsHca$setOfAnnotations, state$annotationsHca$setOfColors)
      })
    }
    drawAnnotationLegendPCA <- function(consoleInfo = NULL){
      output$plotAnnoLegendPCA <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcPlotAnnoLegend(state$annotationsPca$setOfAnnotations, state$annotationsPca$setOfColors)
      })
    }
    drawScoresGroupsLegend <- function(consoleInfo = NULL){
      output$plotScoresGroupsLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcPlotScoresGroupsLegend(state$scoresGroups$groups, state$scoresGroups$colors)
      })
    }
    drawMS2Legend <- function(consoleInfo = NULL){
      output$plotMS2Legend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawMS2LegendImpl()
      })
    }
    drawMS2LegendImpl <- function(){
      calcPlotMS2Legend(dataList = dataList)
    }
    drawDendrogramLegend <- function(consoleInfo = NULL){
      output$calcPlotDendrogramLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawDendrogramLegendImpl()
      })
    }
    drawDendrogramLegendImpl <- function(){
      calcPlotDendrogramLegend()
    }
    drawFragmentDiscriminativityLegend <- function(consoleInfo = NULL){
      output$plotFragmentDiscriminativityLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawFragmentDiscriminativityLegendImpl()
      })
    }
    drawFragmentDiscriminativityLegendImpl <- function(){
      calcPlotDiscriminativityLegend()
    }
    drawFragmentPlot <- function(consoleInfo = NULL){
      output$fragmentPlot <- renderPlot({
        print(paste("### MS2 ### all", consoleInfo))
        drawFragmentPlotImpl()
      })
    }
    drawFragmentPlotImpl <- function(){
      plotFragments(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
    }
    ## plot range resets
    resetHcaPlotRange <- function(){
      dendrogramPlotRange$xMin <<- 1
      dendrogramPlotRange$xMax <<- filterHca$numberOfPrecursorsFiltered
      dendrogramPlotRange$xInterval <<- c(1, filterHca$numberOfPrecursorsFiltered)
      dendrogramPlotRange$xIntervalSize <<- filterHca$numberOfPrecursorsFiltered - 1
    }
    resetPcaPlotRange <- function(){
      minX <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
      maxX <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
      minY <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
      maxY <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
      
      if(any(is.na(c(minX, maxX, minY, maxY)))){
        minX <- -1
        maxX <- 1
        minY <- -1
        maxY <- 1
      }
      
      pcaScoresPlotRange$xMin <<- minX
      pcaScoresPlotRange$xMax <<- maxX
      pcaScoresPlotRange$xInterval <<- c(minX, maxX)
      pcaScoresPlotRange$xIntervalSize <<- maxX - minX
      pcaScoresPlotRange$yMin <<- minY
      pcaScoresPlotRange$yMax <<- maxY
      pcaScoresPlotRange$yInterval <<- c(minY, maxY)
      pcaScoresPlotRange$yIntervalSize <<- maxY - minY
      
      minX <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
      maxX <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
      minY <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
      maxY <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
      
      if(any(is.na(c(minX, maxX, minY, maxY)))){
        minX <- -1
        maxX <- 1
        minY <- -1
        maxY <- 1
      }
      
      pcaLoadingsPlotRange$xMin <<- minX
      pcaLoadingsPlotRange$xMax <<- maxX
      pcaLoadingsPlotRange$xInterval <<- c(minX, maxX)
      pcaLoadingsPlotRange$xIntervalSize <<- maxX - minX
      pcaLoadingsPlotRange$yMin <<- minY
      pcaLoadingsPlotRange$yMax <<- maxY
      pcaLoadingsPlotRange$yInterval <<- c(minY, maxY)
      pcaLoadingsPlotRange$yIntervalSize <<- maxY - minY
    }
    resetMS2PlotRange <- function(){
      ms2PlotRange$xMin <<- dataList$minimumMass
      ms2PlotRange$xMax <<- dataList$maximumMass
      ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
    }
    ## annotation stuff
    addAnnotation <- function(precursorSet, annotationValue, annotationColor){
      if(is.na(match(x = annotationValue, table = dataList$annoPresentAnnotationsList))){
        ## new annotation
        annoIdx <- length(dataList$annoPresentAnnotationsList) + 1
        dataList$annoPresentAnnotationsList[[annoIdx]] <<- annotationValue
        dataList$annoPresentColorsList[[annoIdx]] <<- annotationColor
      }
      ## add
      for(precursor in precursorSet){
        if(annotationValue == artifact)
          ## is artifact!
          dataList$annoArrayIsArtifact[[precursor]] <<- TRUE
        else{
          if(is.na(match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])))
            ## new anno
            dataList$annoArrayOfLists[[precursor]] <<- c(dataList$annoArrayOfLists[[precursor]], annotationValue)
        }
      }
      
      ## update gui
      #print("addAnnotation updateAnnoGui")
      updateAnnoGui(precursorSet)
      updateTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    removeAnnotation <- function(precursorSet, annotationValue){
      ## remove
      for(precursor in precursorSet){
        if(annotationValue == artifact)
          dataList$annoArrayIsArtifact[[precursor]] <<- FALSE
        else{
          idx <- match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])
          dataList$annoArrayOfLists[[precursor]] <<- dataList$annoArrayOfLists[[precursor]][-idx]
        }
      }
      
      ## remove anno completely?
      annoThere <- lapply(X = dataList$annoArrayOfLists, FUN = function(x){ match(x = annotationValue, table = x) })
      if(annotationValue != artifact & all(is.na(annoThere))){
        idx <- match(x = annotationValue, table = dataList$annoPresentAnnotationsList)
        dataList$annoPresentAnnotationsList <<- dataList$annoPresentAnnotationsList[-idx]
        dataList$annoPresentColorsList      <<- dataList$annoPresentColorsList[-idx]
      }
      
      ## update gui
      #print("removeAnnotation updateAnnoGui")
      updateAnnoGui(precursorSet)
      updateTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    setArtifactState <- function(precursorSet, isArtifact){
      ## add
      dataList$annoArrayIsArtifact[precursorSet] <<- isArtifact
      
      ## update gui
      updateAnnoGui(precursorSet)
      updateTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    setAnnotationPrimary <- function(precursorSet, annotationValue){
      ## remove
      for(precursor in precursorSet){
        idx <- match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])
        dataList$annoArrayOfLists[[precursor]] <<- c(dataList$annoArrayOfLists[[precursor]][[idx]], dataList$annoArrayOfLists[[precursor]][-idx])
      }
      
      ## update gui
      #print("setAnnotation primary updateAnnoGui")
      updateAnnoGui(precursorSet)
      updateTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    commonAnnotations <- function(precursorSet){
      if(is.null(precursorSet))
        return(NULL)
      if(all(dataList$annoArrayIsArtifact[precursorSet]))
        return(artifact)
      
      ## at least one non-artifact precursor present
      intersection <- unlist(dataList$annoPresentAnnotationsList)
      for(precursor in precursorSet)
        if(!dataList$annoArrayIsArtifact[[precursor]])
          intersection <- intersect(x = intersection, y = unlist(dataList$annoArrayOfLists[[precursor]]))
      return(intersection)
    }
    ## table update
    updateAnnoGui <- function(precursorSet){
      ## annotation
      commonAnnos <- commonAnnotations(precursorSet)
      
      if(!is.null(commonAnnos))
        updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = commonAnnos)
      else
        updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[none]"))
      updateSelectInput(session = session, inputId = "previousAnnotationValue", choices = dataList$annoPresentAnnotationsList)
      updateTabsetPanel(session = session, inputId = "precursorSelectionTabs", selected = precursorSelectionTabCurrentTab)
      
      if(length(commonAnnos) > 0){
        shinyjs::enable("removePresentAnnotation")
        shinyjs::enable("setPresentAnnotationPrimary")
      } else {
        shinyjs::disable("removePresentAnnotation")
        shinyjs::disable("setPresentAnnotationPrimary")
      }
    }
    updateTableGui <- function(precursorSet){
      print("updateTableGui")
      ## table update with new annotations
      if(all(!is.null(selectionFragmentTreeNodeSet), precursorSet %in% listForTable_Fragment_HCA$precursorSet)){ ## HCA
        table$df_Fragment_HCA <<- createTable(listForTable_Fragment_HCA, selectionFragmentHcaName)
        #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
        table_Fragment_HCA_id <<- tableCheckboxIdCounter
      }
      if(all(!is.null(selectionFragmentPcaLoadingSet), precursorSet %in% listForTable_Fragment_PCA$precursorSet)){ ## PCA
        table$df_Fragment_PCA <<- createTable(listForTable_Fragment_PCA, selectionFragmentPcaName)
        #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
        table_Fragment_PCA_id <<- tableCheckboxIdCounter
      }
      
      if(all(!is.null(selectionAnalysisTreeNodeSet), precursorSet %in% listForTable_Analysis_HCA$precursorSet)){ ## HCA
        table$df_Analysis_HCA <<- createTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
        #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
        table_Analysis_HCA_id <<- tableCheckboxIdCounter
      }
      if(all(!is.null(selectionAnalysisPcaLoadingSet), precursorSet %in% listForTable_Analysis_PCA$precursorSet)){ ## PCA
        table$df_Analysis_PCA <<- createTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
        #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
        table_Analysis_PCA_id <<- tableCheckboxIdCounter
      }
      
      if(all(!is.null(selectionSearchTreeNodeSet), precursorSet %in% listForTable_Search_HCA$precursorSet)){ ## HCA
        table$df_Search_HCA <<- createTable(listForTable_Search_HCA, selectionSearchHcaName)
        #output$dt_Search_HCA <- DT::renderDataTable(table$df_Search_HCA)
        table_Search_HCA_id <<- tableCheckboxIdCounter
      }
      if(all(!is.null(selectionSearchPcaLoadingSet), precursorSet %in% listForTable_Search_PCA$precursorSet)){ ## PCA
        table$df_Search_PCA <<- createTable(listForTable_Search_PCA, selectionSearchPcaName)
        #output$dt_Search_PCA <- DT::renderDataTable(table$df_Search_PCA)
        table_Search_PCA_id <<- tableCheckboxIdCounter
      }
      ## update
      updateTableAssignment()
    }
    updateTableAssignment <- function(){
      print(paste("updateTableAssignment '", state$selectedSelection, "'", sep = ""))
      switch(state$selectedSelection, 
             "Analysis_HCA"={ 
               #selectionAnalysisHcaName={ 
               selectedTable_id <<- table_Analysis_HCA_id
               selectedTable <<- table$df_Analysis_HCA
             },"Analysis_PCA"={  
               #},selectionAnalysisPcaName={  
               selectedTable_id <<- table_Analysis_PCA_id
               selectedTable <<- table$df_Analysis_PCA
             },"Fragment_HCA"={  
               #},selectionFragmentHcaName={  
               selectedTable_id <<- table_Fragment_HCA_id
               selectedTable <<- table$df_Fragment_HCA
               state$precursorSetSelected <<- !is.null(listForTable_Fragment_HCA)
             },"Fragment_PCA"={  
               #},selectionFragmentPcaName={  
               selectedTable_id <<- table_Fragment_PCA_id
               selectedTable <<- table$df_Fragment_PCA
             },"Search_HCA"  ={  
               #},selectionSearchHcaName  ={  
               selectedTable_id <<- table_Search_HCA_id
               selectedTable <<- table$df_Search_HCA
             },"Search_PCA"  ={  
               #},selectionSearchPcaName  ={  
               selectedTable_id <<- table_Search_PCA_id
               selectedTable <<- table$df_Search_PCA
               #output$table <- output$dt_Search_PCA
             },{
               print(paste("### unknown state$selectedSelection: '", state$selectedSelection, "'", sep = ""))
             }
      )
      setTable()
    }
    updateSelectedSelection <- function(){
      selection <- input$changeSelection
      if(selection == selectionAnalysisName & state$analysisType == "HCA")
        selectedSelection <- selectionAnalysisHcaName
      if(selection == selectionAnalysisName & state$analysisType == "PCA")
        selectedSelection <- selectionAnalysisPcaName
      if(selection == selectionFragmentName & state$analysisType == "HCA")
        selectedSelection <- selectionFragmentHcaName
      if(selection == selectionFragmentName & state$analysisType == "PCA")
        selectedSelection <- selectionFragmentPcaName
      if(selection == selectionSearchName & state$analysisType == "HCA")
        selectedSelection <- selectionSearchHcaName
      if(selection == selectionSearchName & state$analysisType == "PCA")
        selectedSelection <- selectionSearchPcaName
      
      state$selectedSelection <<- selectedSelection
      
      updateSelectedPrecursorSet()
    }
    updateSelectedPrecursorSet <- function(){
      print(paste("updateSelectionAssignment '", state$selectedSelection, "'", sep = ""))
      switch(state$selectedSelection, 
             "Analysis_HCA"={ 
               state$precursorSetSelected <<- !is.null(listForTable_Analysis_HCA)
               if(!is.null(listForTable_Analysis_HCA)) selectedPrecursorSet <<- listForTable_Analysis_HCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Analysis_PCA"={  
               state$precursorSetSelected <<- !is.null(listForTable_Analysis_PCA)
               if(!is.null(listForTable_Analysis_PCA)) selectedPrecursorSet <<- listForTable_Analysis_PCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Fragment_HCA"={  
               state$precursorSetSelected <<- !is.null(listForTable_Fragment_HCA)
               if(!is.null(listForTable_Fragment_HCA)) selectedPrecursorSet <<- listForTable_Fragment_HCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Fragment_PCA"={  
               state$precursorSetSelected <<- !is.null(listForTable_Fragment_PCA)
               if(!is.null(listForTable_Fragment_PCA)) selectedPrecursorSet <<- listForTable_Fragment_PCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Search_HCA"  ={  
               state$precursorSetSelected <<- !is.null(listForTable_Search_HCA)
               if(!is.null(listForTable_Search_HCA)) selectedPrecursorSet <<- listForTable_Search_HCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Search_PCA"  ={  
               #output$table <- output$dt_Search_PCA
               state$precursorSetSelected <<- !is.null(listForTable_Search_PCA)
               if(!is.null(listForTable_Search_PCA)) selectedPrecursorSet <<- listForTable_Search_PCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
               #selectedPrecursorSet <<- ifelse(!is.null(listForTable_Search_PCA), listForTable_Search_PCA$precursorSet, NULL)
             },{
               print(paste("### unknown state$selectedSelection: '", state$selectedSelection, "'", sep = ""))
             }
      )
      precursorSelectionChanged()
    }
    precursorSelectionChanged <- function(){
      selectionPresent <- !is.null(selectedPrecursorSet)
      
      ####################
      ## anno, table
      updateAnnoGui(selectedPrecursorSet)
      updateTableGui(selectedPrecursorSet)
      
      ####################
      ## MetFrag link
      if(length(selectedPrecursorSet) == 1){
        landingPageUrl <- getMetFragLink(dataList, selectedPrecursorSet)
        if(!is.null(landingPageUrl)){
          output$metFragLink <- renderText({
            print(paste("update output$metFragLink II", landingPageUrl))
            paste("<a href=", gsub(pattern = " ", replacement = "%20", x = landingPageUrl)[[1]], " target=\"_blank\">Send to MetFrag!</a>", sep = "")
          })
        } else {
          output$metFragLink <- renderText({
            print(paste("update output$metFragLink II empty"))
            paste("", sep = "")
          })
        }
      }
      
      ####################
      ## selection info
      selection <- state$selectedSelection
      selectionInfo <- ""
      if(selectionPresent){
        switch(as.character(length(selectedPrecursorSet)), 
          "0"={ selectionInfo <- paste("The set of selected MS\u00B9 features is empty", sep = "")         },
          "1"={ selectionInfo <- paste(length(selectedPrecursorSet), " MS\u00B9 feature selected", sep = "")  },
          {     selectionInfo <- paste(length(selectedPrecursorSet), " MS\u00B9 features selected", sep = "") }
        )
      } else {
        if(selection == selectionAnalysisHcaName)
          selectionInfo <- paste("Please select a cluster or MS\u00B9 feature in the HCA plot", sep = "")
        if(selection == selectionAnalysisPcaName)
          selectionInfo <- paste("Please select a loading in the PCA plot", sep = "")
        if(selection == selectionFragmentHcaName | selection == selectionFragmentPcaName)
          selectionInfo <- paste("Please select a fragment in the Fragment plot above", sep = "")
        if(selection == selectionSearchHcaName | selection == selectionSearchPcaName)
          selectionInfo <- paste("Please select a set of MS\u00B9 features in the 'Search' tab of the sidebar panel", sep = "")
      }
      
      output$selectionInfo <- renderText({
        print(paste("update output$selectionInfo '", selectionInfo, "'", sep = ""))
        selectionInfo
      })
    }
    updatePlotsWithAnnotations <- function(){
      ## plots
      if(state$showHCAplotPanel){
        drawDendrogramPlot(consoleInfo = "updatePlotsWithAnnotations")
        #drawAnnotationLegendHCA(consoleInfo = "updatePlotsWithAnnotations")
      }
      if(state$showPCAplotPanel){
        drawPcaLoadingsPlot(consoleInfo = "updatePlotsWithAnnotations")
        #drawAnnotationLegendPCA(consoleInfo = "updatePlotsWithAnnotations")
      }
    }
    ## create and set table for all six types
    createInputFields <- function(FUN, id, values) {
      ## running id
      tableCheckboxIdCounter <<- tableCheckboxIdCounter + 1
      id <- paste(id, tableCheckboxIdCounter, sep = "_")
      
      ## create a character vector of shiny inputs
      inputs <- character(length(values))
      for (i in 1:length(values)){
        itemId    <- paste(id, "_", i, sep = "")
        inputs[i] <- as.character(FUN(inputId = itemId, label = NULL, value = values[[i]]))
      }
      return(inputs)
    }
    getInputValues <- function(id, counter, len) {
      id <- paste(id, "_", counter, sep = "")
      
      ## obtain the values of inputs
      unlist(lapply(1:len, function(i) {
        itemId <- paste(id, "_", i, sep = "")
        value  <- input[[itemId]]
        if (is.null(value))
          return(NA)
        else
          return(value)
      }))
    }
    createTable <- function(list, type){
      if(is.null(list$precursorSet))
        return(NULL)
      
      if(length(list$precursorSet) > maximumNumberOfTableEntries){
        precursorSet          <- list$precursorSet[1:maximumNumberOfTableEntries]
        ms1abundanceDataFrame <- list$ms1abundanceDataFrame[1:maximumNumberOfTableEntries, ]
        annotationDataFrame   <- list$annotationDataFrame[1:maximumNumberOfTableEntries, ]
        ms2fragmentDataFrame  <- list$ms2fragmentDataFrame[1:maximumNumberOfTableEntries, ]
      } else {
        precursorSet          <- list$precursorSet
        ms1abundanceDataFrame <- list$ms1abundanceDataFrame
        annotationDataFrame   <- list$annotationDataFrame
        ms2fragmentDataFrame  <- list$ms2fragmentDataFrame
      }
        
      isArtifact <- dataList$annoArrayIsArtifact[precursorSet]
      dataFrameIgnore <- data.frame(Ignore = createInputFields(FUN = checkboxInput, id = paste(type, "Ignore", sep = "_"), values = isArtifact))
      
      dataFrame <<- cbind(
        ms1abundanceDataFrame,
        dataFrameIgnore,
        annotationDataFrame,
        ms2fragmentDataFrame
      )
      
      return(dataFrame)
    }
    setTable <- function(){
      output$table <- DT::renderDataTable(
        expr = selectedTable,
        server = FALSE, escape = FALSE, selection = "none",
        options = list(
          preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
          drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }')
        )
      )
    }
    ## selection stuff
    selectionByFragmentReset <- function(){
      selectionFragmentSelectedFragmentIndex <<- NULL
      
      if(!is.null(selectionFragmentTreeNodeSet)){ ## HCA
        selectionFragmentSelectedFragmentIndex <<- NULL
        selectionFragmentTreeNodeSet <<- NULL
        listForTable_Fragment_HCA <<- NULL
        table_Fragment_HCA_id <<- NULL
        table$df_Fragment_HCA <<- NULL
        #output$dt_Fragment_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionFragmentHcaName)
          updateSelectedPrecursorSet()
      }
      if(!is.null(selectionFragmentPcaLoadingSet)){ ## PCA
        selectionFragmentPcaLoadingSet <<- NULL
        listForTable_Fragment_PCA <<- NULL
        table_Fragment_PCA_id <<- NULL
        table$df_Fragment_PCA <<- NULL
        #output$dt_Fragment_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionFragmentPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByFragment <- function(minimumIndex){
      selectionFragmentSelectedFragmentIndex <<- minimumIndex
      
      fragmentMass  <- fragmentsX[[minimumIndex]]
      fragmentIndex <- which(dataList$fragmentMasses == fragmentMass)
      precursorSet  <- which(dataList$featureMatrix[, fragmentIndex] != 0)
      
      if(state$showHCAplotPanel)
        selectionByFragmentInitHca(precursorSet)
      if(state$showPCAplotPanel)
        selectionByFragmentInitPca(precursorSet)
      
      if(input$changeSelection != selectionFragmentName){
        updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionFragmentName)
        updateSelectedSelection()
      }
    }
    selectionByFragmentInitHca <- function(precursorSet){
      ## HCA - fetch subroots of subtrees comprising the selected fragment
      selectionFragmentTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
      precursorSetHca <- intersect(x = precursorSet, y = filterHca$filter)
      
      if(length(precursorSetHca) > 0){
        listForTable_Fragment_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetHca)
        table$df_Fragment_HCA <<- createTable(listForTable_Fragment_HCA, selectionFragmentHcaName)
        table_Fragment_HCA_id <<- tableCheckboxIdCounter
        #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
        if(state$selectedSelection == selectionFragmentHcaName)
          updateSelectedPrecursorSet()
      } else {
        listForTable_Fragment_HCA <<- NULL
        table_Fragment_HCA_id <<- NULL
        table$df_Fragment_HCA <<- NULL
        #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
        if(state$selectedSelection == selectionFragmentHcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByFragmentInitPca <- function(precursorSet){
      ## PCA
      #selectionFragmentPcaLoadingSet <<- which(dataList$featureMatrix[filterPca$filter, fragmentIndex] != 0)
      selectionFragmentPcaLoadingSet <<- which(filterPca$filter %in% precursorSet)
      precursorSetPca <- intersect(x = precursorSet, y = filterPca$filter)
      
      if(length(precursorSetPca) > 0){
        listForTable_Fragment_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetPca)
        table$df_Fragment_PCA <<- createTable(listForTable_Fragment_PCA, selectionFragmentPcaName)
        table_Fragment_PCA_id <<- tableCheckboxIdCounter
        #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
        if(state$selectedSelection == selectionFragmentPcaName)
          updateSelectedPrecursorSet()
      } else {
        listForTable_Fragment_PCA <<- NULL
        table_Fragment_PCA_id <<- NULL
        table$df_Fragment_PCA <<- NULL
        #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
        if(state$selectedSelection == selectionFragmentPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByAnalysisReset <- function(){
      if(!is.null(selectionAnalysisTreeNodeSet)){ ## HCA
        selectionAnalysisTreeNodeSet <<- NULL
        listForTable_Analysis_HCA <<- NULL
        table_Analysis_HCA_id <<- NULL
        table$df_Analysis_HCA <<- NULL
        #output$dt_Analysis_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisHcaName)
          updateSelectedPrecursorSet()
      }
      if(!is.null(selectionAnalysisPcaLoadingSet)){ ## PCA
        selectionAnalysisPcaLoadingSet <<- NULL
        listForTable_Analysis_PCA <<- NULL
        table_Analysis_PCA_id <<- NULL
        table$df_Analysis_PCA <<- NULL
        #output$dt_Analysis_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByHca <- function(minimumLabel){
      selectionAnalysisTreeNodeSet <<- minimumLabel
      precursorSet <- getPrecursorSetFromTreeSelection(clusterDataList = clusterDataList, clusterLabel = minimumLabel)
      listForTable_Analysis_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
      table$df_Analysis_HCA <<- createTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
      #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
      if(state$selectedSelection == selectionAnalysisHcaName)
        updateSelectedPrecursorSet()
      
      ## pca selection
      if(state$showPCAplotPanel)
        selectionByAnalysisInitPca(precursorSet)
      
      if(input$changeSelection != selectionAnalysisName){
        updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionAnalysisName)
        updateSelectedSelection()
      }
    }
    selectionByAnalysisInitPca <- function(precursorSet){
      precursorSetPca <- intersect(x = precursorSet, y = filterPca$filter)
      
      if(length(precursorSetPca) > 0){
        selectionAnalysisPcaLoadingSet <<- which(filterPca$filter %in% precursorSetPca)
        
        listForTable_Analysis_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetPca)
        table$df_Analysis_PCA <<- createTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
        table_Analysis_PCA_id <<- tableCheckboxIdCounter
        #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
        if(state$selectedSelection == selectionAnalysisPcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionAnalysisPcaLoadingSet <<- NULL
        listForTable_Analysis_PCA <<- NULL
        table_Analysis_PCA_id <<- NULL
        table$df_Analysis_PCA <<- NULL
        #output$dt_Analysis_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByPca <- function(minimumIndex){
      precursorIndex <- filterPca$filter[[minimumIndex]]
      
      selectionAnalysisPcaLoadingSet <<- minimumIndex
      listForTable_Analysis_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorIndex)
      table$df_Analysis_PCA <<- createTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
      #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
      if(state$selectedSelection == selectionAnalysisPcaName)
        updateSelectedPrecursorSet()
      
      if(state$showHCAplotPanel)
        selectionByAnalysisInitHca(precursorIndex)
      
      if(input$changeSelection != selectionAnalysisName){
        updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionAnalysisName)
        updateSelectedSelection()
      }
    }
    selectionByAnalysisInitHca <- function(precursorSet){
      if(any(precursorSet %in% filterHca$filter)){
        #selectionAnalysisTreeNode <<- -match(x = precursorIndex, table = filterHca$filter)
        selectionAnalysisTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
        
        listForTable_Analysis_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
        table$df_Analysis_HCA <<- createTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
        table_Analysis_HCA_id <<- tableCheckboxIdCounter
        #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
        if(state$selectedSelection == selectionAnalysisHcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionAnalysisTreeNodeSet <<- NULL
        listForTable_Analysis_HCA <<- NULL
        table_Analysis_HCA_id <<- NULL
        table$df_Analysis_HCA <<- NULL
        #output$dt_Analysis_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisHcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionBySearchReset <- function(){
      if(!is.null(selectionSearchTreeNodeSet)){ ## HCA
        selectionSearchTreeNodeSet <<- NULL
        listForTable_Search_HCA <<- NULL
        table_Search_HCA_id <<- NULL
        table$df_Search_HCA <<- NULL
        if(state$selectedSelection == selectionSearchHcaName)
          updateSelectedPrecursorSet()
      }
      if(!is.null(selectionSearchPcaLoadingSet)){ ## HCA
        selectionSearchPcaLoadingSet <<- NULL
        listForTable_Search_PCA <<- NULL
        table_Search_PCA_id <<- NULL
        table$df_Search_PCA <<- NULL
        if(state$selectedSelection == selectionSearchHcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionBySearch <- function(precursorSet){
      if(state$showHCAplotPanel)
        selectionBySearchInitHca(precursorSet)
      if(state$showPCAplotPanel)
        selectionBySearchInitPca(precursorSet)
      
      if(input$changeSelection != selectionSearchName){
        updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionSearchName)
        updateSelectedSelection()
      }
    }
    selectionBySearchInitHca <- function(precursorSet){
      precursorSetHca <- intersect(x = precursorSet, y = filterHca$filter)
      
      if(length(precursorSetHca) > 0){
        selectionSearchTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
        listForTable_Search_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetHca)
        table$df_Search_HCA <<- createTable(listForTable_Search_HCA, selectionSearchHcaName)
        table_Search_HCA_id <<- tableCheckboxIdCounter
        #output$dt_Search_HCA <- DT::renderDataTable(table$df_Search_HCA)
        if(state$selectedSelection == selectionSearchHcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionSearchTreeNodeSet <<- NULL
        listForTable_Search_HCA <<- NULL
        table_Search_HCA_id <<- NULL
        table$df_Search_HCA <<- NULL
        #output$dt_Search_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionSearchHcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionBySearchInitPca <- function(precursorSet){
      precursorSetPca <- intersect(x = precursorSet, y = filterPca$filter)
      if(length(precursorSetPca) > 0){
        selectionSearchPcaLoadingSet <<- which(filterPca$filter %in% precursorSetPca)
        listForTable_Search_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetPca)
        table$df_Search_PCA <<- createTable(listForTable_Search_PCA, selectionSearchPcaName)
        table_Search_PCA_id <<- tableCheckboxIdCounter
        #output$dt_Search_PCA <- DT::renderDataTable(table$df_Search_PCA)
        if(state$selectedSelection == selectionSearchPcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionSearchPcaLoadingSet <<- NULL
        listForTable_Search_PCA <<- NULL
        table_Search_PCA_id <<- NULL
        table$df_Search_PCA <<- NULL
        #output$dt_Search_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionSearchPcaName)
          updateSelectedPrecursorSet()
      }
    }
    
    #########################################################################################
    #########################################################################################
    ## observer
    
    ### TODO
    #exampleInput <- reactive(function(){
    #  print(input$example) #for debugging
    #  paste(input$example,"example",sep=".")
    #})
    #output$myImage <- reactive(function(){
    #  #from lattice package documentation
    #  Depth <- equal.count(quakes$depth, number=8, overlap=.1)
    #  xyplot.example <- xyplot(lat ~ long | Depth, data = quakes)
    #  dotplot.example <- dotplot(variety ~ yield | site, data = barley, groups = year,
    #                             key = simpleKey(levels(barley$year), space = "right"),
    #                             xlab = "Barley Yield (bushels/acre) ",
    #                             aspect=0.5, layout = c(1,6), ylab=NULL)
    #  barchart.example <- barchart(yield ~ variety | site, data = barley,
    #                               groups = year, layout = c(1,6), stack = TRUE,
    #                               auto.key = list(space = "right"),
    #                               ylab = "Barley Yield (bushels/acre)",
    #                               scales = list(x = list(rot = 45)))
    #  print(get(exampleInput()))
    #  #print((paste(input$example,"example",sep=".")))
    #  tempsvg <- tempfile(fileext=".svg")
    #  on.exit(unlink(tempsvg))
    #  gridToSVG(name=tempsvg)
    #  
    #  svgoutput <- readLines("/home/htreutle/Downloads/MetSWATH/svg/svgTestScoresPlot01.svg", n=-1)
    #  #svgoutput <- readLines(tempsvg, n=-1)
    #  svgoutput
    #})
    
    ## controls
    obsTabs <- observeEvent(input$runTabs, {
      tabId <- input$runTabs
      print(paste("observe tabs", tabId))
      if(tabId == "HCA"){
        state$analysisType <<- "HCA"
        if(state$showHCAplotPanel){
          state$plotHcaShown <<- TRUE
          state$plotPcaShown <<- FALSE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
        }
      }
      if(tabId == "PCA"){
        state$analysisType <<- "PCA"
        if(state$showPCAplotPanel){
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- TRUE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
        }
      }
      if(tabId == "Input" & !state$initialGuiUpdatePerformed){
        ## initial gui update
        print(paste("update GUI initially", tabId))
        shinyjs::disable("importMs1Ms2Data")
        shinyjs::disable("loadProjectData")
        
        state$initialGuiUpdatePerformed <<- TRUE
      }
      
      if(tabId == "HCA" | tabId == "PCA")
        updateSelectedSelection()
    })
    obsChangePlot <- observeEvent(input$changePlot, {
      plot <- input$changePlot
      print(paste("Observe changePlot", plot))
      if(plot == "Display HCA"){
        analysisType <- "HCA"
        if(state$showHCAplotPanel){
          state$plotHcaShown <<- TRUE
          state$plotPcaShown <<- FALSE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
        }
      }
      if(plot == "Display PCA") {
        analysisType <- "PCA"
        if(state$showPCAplotPanel){
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- TRUE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
        }
      }
      state$analysisType <<- analysisType
      
      updateSelectedSelection()
    })
    obsChangeSelection <- observeEvent(input$changeSelection, {
      selection <- input$changeSelection
      print(paste("Observe changeSelection", selection, "for", state$analysisType, ""))
      
      changeSelectionCurrentSelection <<- selection
      
      updateSelectedSelection()
    })
    obsPrecursorSelectionTabs <- observeEvent(input$precursorSelectionTabs, {
      selectedTab <- input$precursorSelectionTabs
      print(paste("Observe selectedTab", selectedTab))
      
      precursorSelectionTabCurrentTab <<- selectedTab
    })
    obsClearSelection <- observeEvent(input$clearSelection, {
      clearSelection  <- as.numeric(input$clearSelection)
      selection       <- input$changeSelection
      
      print(paste("Observe clearSelection", clearSelection))
      
      #################################################
      ## check if button was hit
      if(clearSelection == clearSelectionButtonValue)
        return()
      clearSelectionButtonValue <<- clearSelection
      
      switch(selection, 
          "Selection by analysis"={
          #selectionAnalysisName={
            selectionByAnalysisReset()
        },"Selection by fragment"={
        #},selectionFragmentName={
          selectionByFragmentReset()
        },"Selection by search"={
        #},selectionSearchName={
          selectionBySearchReset()
        },{
          print(paste("### unknown selection '", selection, "'", sep = ""))
        }
      )
      
      #################################################
      ## update plots
      if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "clear selection")
      if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "clear selection")
    })
    obsFile <- observeEvent(input$matrixFile$datapath, {
      filePath <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      print(paste("Observe file for data", fileName))
      if(!is.null(filePath))
        shinyjs::enable("loadProjectData")
      
      updateFileInputInfo()
    })
    obsLoadProjectData <- observeEvent(input$loadProjectData, {
      loadProjectData <- as.numeric(input$loadProjectData)
      
      print(paste("Observe loadProjectData", loadProjectData))
      
      #################################################
      ## check if button was hit
      if(loadProjectData == loadProjectDataButtonValue)
        return()
      loadProjectDataButtonValue <<- loadProjectData
      
      #################################################
      ## files
      filePath <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      
      loadProjectFile(filePath = filePath, fileName = fileName, buttonId = "loadProjectData")
    })
    obsLoadExampleData <- observeEvent(input$loadExampleData, {
      loadExampleData <- as.numeric(input$loadExampleData)
      
      print(paste("Observe loadExampleData", loadExampleData))
      
      #################################################
      ## check if button was hit
      if(loadExampleData == loadExampleDataButtonValue)
        return()
      loadProjectDataButtonValue <<- loadExampleData
      
      #exampleDataSelection <- input$exampleDataSelection
      
      #################################################
      ## files
      ## TODO relative paths or as R package
      #filePath <- paste(shinyAppFolder, "files/Fragment_matrix_showcase.csv", sep = "")
      filePath <- paste(shinyAppFolder, "files/Project_file_showcase_annotated.csv.gz", sep = "")
      fileName <- "Fragment_matrix_showcase.csv"
      #if(exampleDataSelection == "Example data set (full)"){
      #  filePath <- paste(shinyAppFolder, "files/Fragment_matrix_showcase.csv", sep = "")
      #  fileName <- "Fragment_matrix_showcase.csv"
      #}
      #if(exampleDataSelection == "Example data set (reduced)"){
      #  filePath <- paste(shinyAppFolder, "files/Project_file_showcase_reduced.csv.gz", sep = "")
      #  fileName <- "Project_file_showcase.csv.gz"
      #}
      
      loadProjectFile(filePath = filePath, fileName = fileName, buttonId = "loadExampleData")
    })
    #obsExampleDataSelection <- observeEvent(input$exampleDataSelection, {
    #  updateFileInputInfo()
    #})
    loadProjectFile <- function(filePath, fileName, buttonId){
      #########################################################################################
      ## read data
      #output$fileInfo <- renderText({paste("Please wait for the project data to be loaded...")})
      
      session$sendCustomMessage("disableButton", buttonId)
      
      error <- NULL
      withProgress(message = 'Reading file...', value = 0, {
        dataList <<- tryCatch(
        {
          #dataList <<- calcClusterData(file = "/mnt/VOL1/ABT/Alle/Balcke/MetSWATH/data/MS-DIAL/UC Davis/Results/201558139_matrixPrecursorsVersusFragmentsDeisotoped_withoutZerosTest01.txt")
          readClusterDataFromProjectFile(file = filePath, progress = TRUE)
        }, error = function(e) {
          error <- e
        }
        )
      })
      
      if(!is.null(error)){
        output$fileInfo <- renderText({paste("There occurred an error while processing the project file. Please check the file format and content and try again.")})
        session$sendCustomMessage("enableButton", buttonId)
        return()
      }
      
      print(paste("readClusterDataFromProjectFile finished", dataList$minimumMass))
      
      resetWorkspace()
      
      state$importedOrLoadedFile_s_ <<- fileName
      updateFileInputInfo()
      
      session$sendCustomMessage("enableButton", buttonId)
    }
    obsImportMs1DataFile <- observeEvent(input$ms1DataFile$datapath, {
      fileMs1Path <- input$ms1DataFile$datapath
      fileMs1Name <- input$ms1DataFile$name
      fileMs2Path <- input$ms2DataFile$datapath
      fileMs2Name <- input$ms2DataFile$name
      print(paste("Observe import MS1 file", fileMs1Name))
      
      if(all(!is.null(fileMs1Path), !is.null(fileMs2Path)))
        shinyjs::enable("importMs1Ms2Data")
      else
        shinyjs::disable("importMs1Ms2Data")
      
      updateFileInputInfo()
      #processMs1Ms2Files(fileMs1Path = fileMs1Path, fileMs2Path = fileMs2Path)
    })
    obsImportMs2DataFile <- observeEvent(input$ms2DataFile$datapath, {
      fileMs1Path <- input$ms1DataFile$datapath
      fileMs1Name <- input$ms1DataFile$name
      fileMs2Path <- input$ms2DataFile$datapath
      fileMs2Name <- input$ms2DataFile$name
      print(paste("Observe import MS2 file", fileMs2Name))
      
      if(all(!is.null(fileMs1Path), !is.null(fileMs2Path)))
        shinyjs::enable("importMs1Ms2Data")
      else
        shinyjs::disable("importMs1Ms2Data")
      
      updateFileInputInfo()
      #processMs1Ms2Files(fileMs1Path = fileMs1Path, fileMs2Path = fileMs2Path)
    })
    obsImportMs1Ms2Data <- observeEvent(input$importMs1Ms2Data, {
      importMs1Ms2Data <- as.numeric(input$importMs1Ms2Data)
      
      print(paste("Observe importMs1Ms2Data", importMs1Ms2Data))
      
      #################################################
      ## check if button was hit
      if(importMs1Ms2Data == importMs1Ms2DataButtonValue)
        return()
      importMs1Ms2DataButtonValue <<- importMs1Ms2Data
      
      #################################################
      ## files
      fileMs1Path <- input$ms1DataFile$datapath
      fileMs1Name <- input$ms1DataFile$name
      fileMs2Path <- input$ms2DataFile$datapath
      fileMs2Name <- input$ms2DataFile$name
      
      #################################################
      ## params
      
      ## project name
      projectName <- input$projectName
      projectName <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectName)))
      projectDescription <- input$projectDescription
      projectDescription <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectDescription)))
      
      ## minimum MS2 peak intensity
      minimumIntensityOfMaximalMS2peak <- input$minimumIntensityOfMaximalMS2peak
      minimumProportionOfMS2peaks <- input$minimumProportionOfMS2peaks
      ## grouping of MS2 peaks
      mzDeviationAbsolute_grouping <- input$mzDeviationAbsolute_grouping
      mzDeviationInPPM_grouping <- input$mzDeviationInPPM_grouping
      ## precursor deisotoping
      doPrecursorDeisotoping <- input$doPrecursorDeisotoping
      mzDeviationAbsolute_precursorDeisotoping <- input$mzDeviationAbsolute_precursorDeisotoping
      mzDeviationInPPM_precursorDeisotoping <- input$mzDeviationInPPM_precursorDeisotoping
      maximumRtDifference <- input$maximumRtDifference
      ## fragment deisotoping
      doMs2PeakGroupDeisotoping <- input$doMs2PeakGroupDeisotoping
      mzDeviationAbsolute_ms2PeakGroupDeisotoping <- input$mzDeviationAbsolute_ms2PeakGroupDeisotoping
      mzDeviationInPPM_ms2PeakGroupDeisotoping <- input$mzDeviationInPPM_ms2PeakGroupDeisotoping
      ## neutral losses
      neutralLossesPrecursorToFragments <- input$neutralLossesPrecursorToFragments
      neutralLossesFragmentsToFragments <- input$neutralLossesFragmentsToFragments
      #neutralLossesPrecursorToFragments <- TRUE
      #neutralLossesFragmentsToFragments <- FALSE
      
      ## fixed
      proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere <- proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
      mzDeviationAbsolute_mappingHere <- mzDeviationAbsolute_mapping
      minimumNumberOfMS2PeaksPerGroupHere <- minimumNumberOfMS2PeaksPerGroup
      
      #################################################
      ## check params
      error <- FALSE
      if(any(is.null(minimumIntensityOfMaximalMS2peak), length(minimumIntensityOfMaximalMS2peak) == 0, nchar(minimumIntensityOfMaximalMS2peak) == 0))
        error <- TRUE
      else{
        minimumIntensityOfMaximalMS2peak <- as.numeric(minimumIntensityOfMaximalMS2peak)
        error <- error | is.na(minimumIntensityOfMaximalMS2peak)
      }
      if(any(is.null(minimumProportionOfMS2peaks), length(minimumProportionOfMS2peaks) == 0, nchar(minimumProportionOfMS2peaks) == 0))
        error <- TRUE
      else{
        minimumProportionOfMS2peaks <- as.numeric(minimumProportionOfMS2peaks)
        error <- error | is.na(minimumProportionOfMS2peaks)
      }
      if(any(is.null(mzDeviationAbsolute_grouping), length(mzDeviationAbsolute_grouping) == 0, nchar(mzDeviationAbsolute_grouping) == 0))
        error <- TRUE
      else{
        mzDeviationAbsolute_grouping <- as.numeric(mzDeviationAbsolute_grouping)
        error <- error | is.na(mzDeviationAbsolute_grouping)
      }
      if(any(is.null(mzDeviationInPPM_grouping), length(mzDeviationInPPM_grouping) == 0, nchar(mzDeviationInPPM_grouping) == 0))
        error <- TRUE
      else{
        mzDeviationInPPM_grouping <- as.numeric(mzDeviationInPPM_grouping)
        error <- error | is.na(mzDeviationInPPM_grouping)
      }
      if(doPrecursorDeisotoping){
        if(any(is.null(mzDeviationAbsolute_precursorDeisotoping), length(mzDeviationAbsolute_precursorDeisotoping) == 0, nchar(mzDeviationAbsolute_precursorDeisotoping) == 0))
          error <- TRUE
        else{
          mzDeviationAbsolute_precursorDeisotoping <- as.numeric(mzDeviationAbsolute_precursorDeisotoping)
          error <- error | is.na(mzDeviationAbsolute_precursorDeisotoping)
        }
        if(any(is.null(mzDeviationInPPM_precursorDeisotoping), length(mzDeviationInPPM_precursorDeisotoping) == 0, nchar(mzDeviationInPPM_precursorDeisotoping) == 0))
          error <- TRUE
        else{
          mzDeviationInPPM_precursorDeisotoping <- as.numeric(mzDeviationInPPM_precursorDeisotoping)
          error <- error | is.na(mzDeviationInPPM_precursorDeisotoping)
        }
        if(any(is.null(maximumRtDifference), length(maximumRtDifference) == 0, nchar(maximumRtDifference) == 0))
          error <- TRUE
        else{
          maximumRtDifference <- as.numeric(maximumRtDifference)
          error <- error | is.na(maximumRtDifference)
        }
      }
      if(doMs2PeakGroupDeisotoping){
        if(any(is.null(mzDeviationAbsolute_ms2PeakGroupDeisotoping), length(mzDeviationAbsolute_ms2PeakGroupDeisotoping) == 0, nchar(mzDeviationAbsolute_ms2PeakGroupDeisotoping) == 0))
          error <- TRUE
        else{
          mzDeviationAbsolute_ms2PeakGroupDeisotoping <- as.numeric(mzDeviationAbsolute_ms2PeakGroupDeisotoping)
          error <- error | is.na(mzDeviationAbsolute_ms2PeakGroupDeisotoping)
        }
        if(any(is.null(mzDeviationInPPM_ms2PeakGroupDeisotoping), length(mzDeviationInPPM_ms2PeakGroupDeisotoping) == 0, nchar(mzDeviationInPPM_ms2PeakGroupDeisotoping) == 0))
          error <- TRUE
        else{
          mzDeviationInPPM_ms2PeakGroupDeisotoping <- as.numeric(mzDeviationInPPM_ms2PeakGroupDeisotoping)
          error <- error | is.na(mzDeviationInPPM_ms2PeakGroupDeisotoping)
        }
      }
      
      if(any(is.null(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere), length(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere) == 0, nchar(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere) == 0))
        error <- TRUE
      else{
        proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere <- as.numeric(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere)
        error <- error | is.na(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere)
      }
      if(any(is.null(mzDeviationAbsolute_mappingHere), length(mzDeviationAbsolute_mappingHere) == 0, nchar(mzDeviationAbsolute_mappingHere) == 0))
        error <- TRUE
      else{
        mzDeviationAbsolute_mappingHere <- as.numeric(mzDeviationAbsolute_mappingHere)
        error <- error | is.na(mzDeviationAbsolute_mappingHere)
      }
      if(any(is.null(minimumNumberOfMS2PeaksPerGroupHere), length(minimumNumberOfMS2PeaksPerGroupHere) == 0, nchar(minimumNumberOfMS2PeaksPerGroupHere) == 0))
        error <- TRUE
      else{
        minimumNumberOfMS2PeaksPerGroupHere <- as.numeric(minimumNumberOfMS2PeaksPerGroupHere)
        error <- error | is.na(minimumNumberOfMS2PeaksPerGroupHere)
      }
      
      if(error){
        output$fileInfo <- renderText({paste("There are invalid parameter values. Please check the parameters and press 'Import MS\u00B9 and MS/MS data' again.")})
        return()
      }
      
      ## box parameters
      print(paste("Observe importMs1Ms2Data", "e", error, "mi", minimumIntensityOfMaximalMS2peak, "mp", minimumProportionOfMS2peaks, "ga", mzDeviationAbsolute_grouping, "gr", mzDeviationInPPM_grouping, "pd", doPrecursorDeisotoping, "pa", mzDeviationAbsolute_precursorDeisotoping, "pr", mzDeviationInPPM_precursorDeisotoping, "mr", maximumRtDifference, "fd", doMs2PeakGroupDeisotoping, "fa", mzDeviationAbsolute_ms2PeakGroupDeisotoping, "fr", mzDeviationInPPM_ms2PeakGroupDeisotoping, "pm", proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere, "ma", mzDeviationAbsolute_mappingHere, "mn", minimumNumberOfMS2PeaksPerGroupHere))
      parameterSet <- list()
      parameterSet$projectName                                       <- projectName
      parameterSet$projectDescription                                <- projectDescription
      parameterSet$toolVersion                                       <- paste(toolName, toolVersion, sep = " ")
      parameterSet$minimumIntensityOfMaximalMS2peak                  <- minimumIntensityOfMaximalMS2peak
      parameterSet$minimumProportionOfMS2peaks                       <- minimumProportionOfMS2peaks
      parameterSet$mzDeviationAbsolute_grouping                      <- mzDeviationAbsolute_grouping
      parameterSet$mzDeviationInPPM_grouping                         <- mzDeviationInPPM_grouping
      parameterSet$doPrecursorDeisotoping                            <- doPrecursorDeisotoping
      parameterSet$mzDeviationAbsolute_precursorDeisotoping          <- mzDeviationAbsolute_precursorDeisotoping
      parameterSet$mzDeviationInPPM_precursorDeisotoping             <- mzDeviationInPPM_precursorDeisotoping
      parameterSet$maximumRtDifference                               <- maximumRtDifference
      parameterSet$doMs2PeakGroupDeisotoping                         <- doMs2PeakGroupDeisotoping
      parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping       <- mzDeviationAbsolute_ms2PeakGroupDeisotoping
      parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping          <- mzDeviationInPPM_ms2PeakGroupDeisotoping
      parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere
      parameterSet$mzDeviationAbsolute_mapping                       <- mzDeviationAbsolute_mappingHere
      parameterSet$minimumNumberOfMS2PeaksPerGroup                   <- minimumNumberOfMS2PeaksPerGroupHere
      parameterSet$neutralLossesPrecursorToFragments                 <- neutralLossesPrecursorToFragments
      parameterSet$neutralLossesFragmentsToFragments                 <- neutralLossesFragmentsToFragments
      
      #################################################
      ## convert to project file
      session$sendCustomMessage("disableButton", "importMs1Ms2Data")
      
      ## built matrix
      error <- NULL
      withProgress(message = 'Generating matrix...', value = 0, {
        resultObj <- tryCatch(
          {
            convertToProjectFile(
              filePeakMatrix = fileMs1Path, 
              fileSpectra = fileMs2Path, 
              parameterSet = parameterSet, 
              progress = TRUE
            )
          }, error = function(e) {
            error <- e
          }
        )
      })
      
      if(!is.null(error)){
        output$fileInfo <- renderText({paste("There occurred an error while processing the input files. Please check the file format and content and try again.")})
        session$sendCustomMessage("enableButton", "importMs1Ms2Data")
        return()
      }
      if(resultObj == "Number of spectra is zero"){
        output$fileInfo <- renderText({paste("There are no MS/MS spectra which fulfill the given criteria. Please refine parameter 'Spectrum intensity' and try 'Import MS\u00B9 and MS/MS data' again.")})
        session$sendCustomMessage("enableButton", "importMs1Ms2Data")
        return()
      }
      
      matrixRows <- resultObj$matrixRows
      matrixCols <- resultObj$matrixCols
      matrixVals <- resultObj$matrixVals
      
      ## convert matrix to dataframe
      numberOfRows    <- max(matrixRows)
      numberOfColumns <- max(matrixCols)
      
      fragmentMatrix <- matrix(data = rep(x = "", times = numberOfRows * numberOfColumns), nrow = numberOfRows, ncol = numberOfColumns)
      fragmentMatrix[cbind(matrixRows, matrixCols)] <- matrixVals
      
      #dataFrame <- as.data.frame(as.matrix(sparseMatrix(i = matrixRows, j = matrixCols, x = matrixVals)))
      dataFrame <- as.data.frame(fragmentMatrix, stringsAsFactors = FALSE)
      
      ## add import parameters
      dataFrame[[1, 1]] <- serializeParameterSet(parameterSet)
      #################################################
      ## process project file
      withProgress(message = 'Processing matrix...', value = 0, {
        dataList <<- readProjectData(dataFrame = dataFrame, progress = TRUE)
      })
      print(paste("readProjectData do data finished", dataList$minimumMass))
      resetWorkspace()
      
      state$importedOrLoadedFile_s_ <<- c(fileMs1Name, fileMs2Name)
      updateFileInputInfo()
      
      session$sendCustomMessage("enableButton", "importMs1Ms2Data")
    })
    obsFileInputSelection <- observeEvent(input$fileInputSelection, {
      updateFileInputInfo()
    })
    obsApplyImportParameterFile <- observeEvent(input$importParameterFileInput$datapath, {
      filePath <- input$importParameterFileInput$datapath
      fileName <- input$importParameterFileInput$name
      print(paste("Observe importParameterFile", fileName))
      if(is.null(filePath))
        return()
      
      ## read and parse
      fileContent <- readLines(con = filePath)
      parameterSet <- deserializeParameterSetFile(fileContent)
      
      ## apply
      #projectName2 <- parameterSet$projectName
      #projectName2 <- paste(projectName2, " adopted", sep = "")
      #updateTextInput    (session = session, inputId = "projectName",                                 value = projectName2)
      #parameterSet$toolVersion
      updateTextInput    (session = session, inputId = "minimumIntensityOfMaximalMS2peak",            value = parameterSet$minimumIntensityOfMaximalMS2peak)
      updateTextInput    (session = session, inputId = "minimumProportionOfMS2peaks",                 value = parameterSet$minimumProportionOfMS2peaks)
      updateTextInput    (session = session, inputId = "mzDeviationAbsolute_grouping",                value = parameterSet$mzDeviationAbsolute_grouping)
      updateTextInput    (session = session, inputId = "mzDeviationInPPM_grouping",                   value = parameterSet$mzDeviationInPPM_grouping)
      updateCheckboxInput(session = session, inputId = "doPrecursorDeisotoping",                      value = parameterSet$doPrecursorDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationAbsolute_precursorDeisotoping",    value = parameterSet$mzDeviationAbsolute_precursorDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationInPPM_precursorDeisotoping",       value = parameterSet$mzDeviationInPPM_precursorDeisotoping)
      updateTextInput    (session = session, inputId = "maximumRtDifference",                         value = parameterSet$maximumRtDifference)
      updateCheckboxInput(session = session, inputId = "doMs2PeakGroupDeisotoping",                   value = parameterSet$doMs2PeakGroupDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping", value = parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping)
      updateTextInput    (session = session, inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping",    value = parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping)
      #parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
      #parameterSet$mzDeviationAbsolute_mapping
      #parameterSet$minimumNumberOfMS2PeaksPerGroup
      updateCheckboxInput(session = session, inputId = "neutralLossesPrecursorToFragments",           value = parameterSet$neutralLossesPrecursorToFragments)
      updateCheckboxInput(session = session, inputId = "neutralLossesFragmentsToFragments",           value = parameterSet$neutralLossesFragmentsToFragments)
    })
    obsUpdateProjectDescription <- observeEvent(input$updateProjectDescription, {
      updateProjectDescription <- as.numeric(input$updateProjectDescription)
      
      print(paste("Observe updateProjectDescription", updateProjectDescription))
      
      #################################################
      ## check if button was hit
      if(updateProjectDescription == updateProjectDescriptionButtonValue)
        return()
      updateProjectDescriptionButtonValue <<- updateProjectDescription
      
      projectDescription <- input$projectDescription2
      projectDescription <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectDescription)))
      dataList$importParameterSet$projectDescription <<- projectDescription
    })
    obsShowSideBar <- observeEvent(input$showSideBar, {
      showSideBar <- input$showSideBar
      print(paste("Observe showSideBar", showSideBar))
      
      state$showSideBar <<- showSideBar
      showSideBar <<- showSideBar
      #output$showSideBar <- showSideBar
      
      if(showSideBar){
        state$runRightColumnWidth <<- runRightColumnWidthPart
        state$legendColumnWidth <<- legendColumnWidthPart
      }
      else{
        state$runRightColumnWidth <<- runRightColumnWidthFull
        state$legendColumnWidth <<- legendColumnWidthFull
      }
      
      ## restore gui state
      updateCheckboxInput(session = session, inputId = "showClusterLabels", value = state$showClusterLabels)
      updateCheckboxInput(session = session, inputId = "showScoresLabels", value = state$showScoresLabels)
      #values <- c(ifelse(state$showLoadingsLabels, "Show labels", ""), ifelse(state$showLoadingsAbundance, "Show abundance", ""))
      #values <- values[values != ""]
      #if(length(values) > 0)
      #  updateCheckboxGroupInput(session = session, inputId = "pcaLoadingsProperties", value = values)
      updateCheckboxInput(session = session, inputId = "showLoadingsLabels", value = state$showLoadingsLabels)
      updateCheckboxInput(session = session, inputId = "showLoadingsAbundance", value = state$showLoadingsAbundance)
    })
    obsShowClusterLabels <- observeEvent(input$showClusterLabels, {
      showClusterLabels <- input$showClusterLabels
      print(paste("Observe showClusterLabels", showClusterLabels))
      state$showClusterLabels <<- showClusterLabels
      drawDendrogramPlot(consoleInfo = "showClusterLabels")
    })
    obsShowScoresLabels <- observeEvent(input$showScoresLabels, {
      showScoresLabels <- input$showScoresLabels
      print(paste("Observe showScoresLabels", showScoresLabels))
      state$showScoresLabels <<- showScoresLabels
      drawPcaScoresPlot(consoleInfo = "showScoresLabels")
    })
    # observePcaLoadingsProperties <- observeEvent(input$pcaLoadingsProperties, {
    #   pcaLoadingsProperties <- input$pcaLoadingsProperties
    #   print(paste("observe pcaLoadingsProperties", pcaLoadingsProperties))
    #   
    #   switch(as.character(length(pcaLoadingsProperties)),
    #          "0"={## 
    #            showLoadingsLabels <- FALSE
    #            showLoadingsAbundance <- FALSE
    #          },
    #          "1"={## 
    #            if(pcaLoadingsProperties == "Show labels"){
    #              showLoadingsLabels <- TRUE
    #              showLoadingsAbundance <- FALSE
    #            } else {
    #              showLoadingsLabels <- FALSE
    #              showLoadingsAbundance <- TRUE
    #            }
    #          },
    #          "2"={## 
    #            showLoadingsLabels <- TRUE
    #            showLoadingsAbundance <- TRUE
    #          },
    #          {## multiple annotations --> take the one which is primary
    #            stop(paste("unknown pcaLoadingsProperties state: ", paste(pcaLoadingsProperties, collapse = "; ")))
    #          }
    #   )## end switch
    #   
    #   state$showLoadingsLabels <<- showLoadingsLabels
    #   state$showLoadingsAbundance <<- showLoadingsAbundance
    #   drawPcaLoadingsPlot(consoleInfo = "pcaLoadingsProperties")
    # })
    obsShowLoadingsLabels <- observeEvent(input$showLoadingsLabels, {
      showLoadingsLabels <- input$showLoadingsLabels
      print(paste("Observe showLoadingsLabels", showLoadingsLabels))
      state$showLoadingsLabels <<- showLoadingsLabels
      drawPcaLoadingsPlot(consoleInfo = "showLoadingsLabels")
    })
    obsShowLoadingsAbundance <- observeEvent(input$showLoadingsAbundance, {
      showLoadingsAbundance <- input$showLoadingsAbundance
      print(paste("Observe showLoadingsAbundance", showLoadingsAbundance))
      state$showLoadingsAbundance <<- showLoadingsAbundance
      drawPcaLoadingsPlot(consoleInfo = "showLoadingsAbundance")
    })
    ## listen to filter events
    obsClearSearch <- observeEvent(input$clearSearch, {
      clearSearch <- as.numeric(input$clearSearch)
      
      print(paste("Observe clearSearch", clearSearch))
      
      #################################################
      ## check if button was hit
      if(clearSearch == clearSearchButtonValue)
        return()
      clearSearchButtonValue <<- clearSearch
      
      filterSearch <<- NULL
      filterSearchActive <<- FALSE
      state$searchfilterValid <<- TRUE
      
      selectionBySearchReset()
      updateSearchInformation()
      
      #################################################
      ## update plots
      if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "clear search")
      if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "clear search")
    })
    obsApplySearch <- observeEvent(input$applySearch, {
      applySearch <- as.numeric(input$applySearch)
      
      print(paste("Observe applySearch", applySearch))
      
      #################################################
      ## check if button was hit
      if(applySearch == applySearchButtonValue)
        return()
      applySearchButtonValue <<- applySearch
      
      #################################################
      ## MS1 or MS2?
      searchMode <- input$searchMS1orMS2
      if(searchMode == 'MS\u00B9 feature mass'){
        #################################################
        ## get inputs
        filter_ms1_masses <- input$searchMS1mass
        filter_ms1_ppm    <- input$searchMS1massPpm
        
        if(nchar(trimws(filter_ms1_masses)) == 0)
          return()
        
        filter_ms2_masses1  <- NULL
        filter_ms2_masses2  <- NULL
        filter_ms2_masses3  <- NULL
        filter_ms2_ppm      <- NULL
      }
      if(searchMode == 'Fragment mass'){
        #################################################
        ## get inputs
        filter_ms2_masses1 <- input$search_ms2_masses1
        filter_ms2_masses2 <- input$search_ms2_masses2
        filter_ms2_masses3 <- input$search_ms2_masses3
        filter_ms2_ppm     <- input$searchMS2massPpm
        
        filter_ms1_masses <- NULL
        filter_ms1_ppm    <- NULL
      }
      
      filter_lfc      <- NULL
      filter_average  <- NULL
      groupSet        <- dataList$groups
      includeIgnoredPrecursors  <- input$searchIncludeIgnoredPrecursors
      
      #################################################
      ## do filtering
      print(paste("Observe applySearch", "1m", filter_ms1_masses, "1p", filter_ms1_ppm, "i", includeIgnoredPrecursors, "gs", paste(groupSet, collapse = "-")))
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      processSearchFilterResult(resultObj)
    })
    obsApplyGlobalMS2filters <- observeEvent(input$applyGlobalMS2filters, {
      applyGlobalMS2filters <- as.numeric(input$applyGlobalMS2filters)
      
      print(paste("Observe applyGlobalMS2filters", applyGlobalMS2filters))
      
      #################################################
      ## check if button was hit
      if(applyGlobalMS2filters == applyGlobalMS2filtersButtonValue)
        return()
      applyGlobalMS2filtersButtonValue <<- applyGlobalMS2filters
      
      #################################################
      ## get inputs
      filter_ms2_masses1  <- input$globalFilter_ms2_masses1
      filter_ms2_masses2  <- input$globalFilter_ms2_masses2
      filter_ms2_masses3  <- input$globalFilter_ms2_masses3
      filter_ms2_ppm      <- input$globalFilter_ms2_ppm
      
      groupSet        <- dataList$groups
      filter_average  <- NULL
      filter_lfc      <- NULL
      includeIgnoredPrecursors  <- TRUE
      filter_ms1_masses <- NULL
      filter_ms1_ppm <- NULL
      
      #################################################
      ## do filtering
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        filterGlobal <<- NULL
        state$globalMS2filterValid <<- FALSE
      } else {
        filterGlobal <<- resultObj$filter
        state$globalMS2filterValid <<- TRUE
      }
      updateGlobalMS2filterInformation()
    })
    obsApplyHcaFilters <- observeEvent(input$applyHcaFilters, {
      applyHcaFilters <- as.numeric(input$applyHcaFilters)
      
      print(paste("Observe applyHcaFilters", applyHcaFilters))
      
      #################################################
      ## check if button was hit
      if(applyHcaFilters == applyHcaFiltersButtonValue)
        return()
      applyHcaFiltersButtonValue <<- applyHcaFilters
      
      #################################################
      ## get inputs
      #filter_ms2_masses1  <- input$globalFilter_ms2_masses1
      #filter_ms2_masses2  <- input$globalFilter_ms2_masses2
      #filter_ms2_masses3  <- input$globalFilter_ms2_masses3
      #filter_ms2_ppm      <- input$globalFilter_ms2_ppm
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      groupOne        <- input$hcaFilterGroupOne
      groupTwo        <- input$hcaFilterGroupTwo
      filter_average  <- input$hcaFilter_average
      filter_lfc      <- input$hcaFilter_lfc
      includeIgnoredPrecursors  <- input$hcaFilterIncludeIgnoredPrecursors
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      groupSet        <- c(groupOne, groupTwo)
      
      #################################################
      ## do filtering and update
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        shinyjs::disable("drawHCAplots")
        #disableActionButton(session, "drawHCAplots")
        filterHca <<- NULL
        updateHcaFilterInformation()
        state$hcaFilterValid <<- FALSE
        return()
      }
      
      #################################################
      ## check filter validity
      filterHca <<- resultObj$filter
      updateHcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
      checkHcaFilterValidity(numberOfPrecursorsFiltered)
    })
    obsApplyPcaFilters <- observeEvent(input$applyPcaFilters, {
      applyPcaFilters <- as.numeric(input$applyPcaFilters)
      
      print(paste("Observe applyPcaFilters", applyPcaFilters))
      
      #################################################
      ## check if button was hit
      if(applyPcaFilters == applyPcaFiltersButtonValue)
        return()
      applyPcaFiltersButtonValue <<- applyPcaFilters
      
      #################################################
      ## get inputs
      #filter_ms2_masses1  <- input$globalFilter_ms2_masses1
      #filter_ms2_masses2  <- input$globalFilter_ms2_masses2
      #filter_ms2_masses3  <- input$globalFilter_ms2_masses3
      #filter_ms2_ppm      <- input$globalFilter_ms2_ppm
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      groupSet        <- input$pcaGroups
      filter_average  <- input$pcaFilter_average
      filter_lfc      <- input$pcaFilter_lfc
      includeIgnoredPrecursors  <- input$pcaFilterIncludeIgnoredPrecursors
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      
      if(length(groupSet) != 2)
        filter_lfc <- NULL
      
      #################################################
      ## do filtering and update
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        shinyjs::disable("drawPCAplots")
        #disableActionButton(session, "drawPCAplots")
        filterPca <<- NULL
        updatePcaFilterInformation()
        state$pcaFilterValid <<- FALSE
        return()
      }
      
      filterPca <<- resultObj$filter
      updatePcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterPca$numberOfPrecursorsFiltered
      checkPcaFilterValidity(numberOfPrecursorsFiltered)
    })
    observeGroupSet <- observeEvent(input$pcaGroups, {
      print(paste("observe groups change", paste(input$pcaGroups, collapse = "-"), length(input$pcaGroups), length(input$pcaGroups) == 2))
      shinyjs::toggleState("pcaFilter_lfc", length(input$pcaGroups) == 2)
    })
    ## listen to draw button events
    obsDrawHCA <- observeEvent(input$drawHCAplots, {
      #################################################
      ## get input
      drawPlots <- as.numeric(input$drawHCAplots)
      
      print(paste("Observe draw HCA plots", drawPlots))
      
      ## check if button was hit
      if(drawPlots == drawHCAButtonValue)
        return()
      drawHCAButtonValue <<- drawPlots
      
      distanceMeasure <- input$hcaDistanceFunction
      clusterMethod <- input$hcaClusterMethod
      print(paste("Observe draw HCA plots", "D", distanceMeasure, "M", clusterMethod))
      
      ##########################
      ## calc
      
      ## compute distance matrix
      withProgress(message = 'Calculating distances...', value = 0, {
        currentDistanceMatrixObj <<- calculateDistanceMatrix(dataList = dataList, filter = filterHca$filter, distanceMeasure = distanceMeasure, progress = TRUE)
      })
      ## compute cluster
      withProgress(message = 'Calculating cluster...', value = 0, {
        clusterDataList <<- calculateCluster(progress = TRUE, dataList = dataList, filter = filterHca$filter, distanceMatrix = currentDistanceMatrixObj$distanceMatrix, method = clusterMethod)
      })
      
      ##########################
      ## hca selections
      if(!is.null(listForTable_Fragment_PCA)){ ## selection from fragment
        fragmentIndex <- which(dataList$fragmentMasses == fragmentsX[[selectionFragmentSelectedFragmentIndex]])
        precursorSet  <- which(dataList$featureMatrix[, fragmentIndex] != 0)
        selectionByFragmentInitHca(precursorSet)
      } else {
        selectionByFragmentReset()
      }
      if(!is.null(listForTable_Analysis_PCA)){ ## selection from PCA
        precursorSet <- filterPca$filter[selectionAnalysisPcaLoadingSet]
        selectionByAnalysisInitHca(precursorSet)
      } else {
        selectionByAnalysisReset()
      }
      if(!is.null(filterSearch)){ ## selection from search
        selectionBySearchInitHca(filterSearch$filter)
      } else {
        selectionBySearchReset()
      }
      
      ##########################
      ## reset MS2 stuff
      if(!state$showPCAplotPanel){
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
      }
      
      ##########################
      ## draw
      resetHcaPlotRange()
      drawDendrogramPlot(consoleInfo = "init output$plotDendrogram", withHeatmap = TRUE)
      drawMS2Plot(consoleInfo = "init output$plotMS2")
      drawAnnotationLegendHCA(consoleInfo = "init output$plotAnnoLegend")
      
      if(!state$anyPlotDrawn){
        drawMS2Legend(consoleInfo = "init output$ms2LegendPlot")
        drawFragmentDiscriminativityLegend(consoleInfo = "init output$plotFragmentDiscriminativityLegend")
        drawHeatmapLegend(consoleInfo = "init output$plotHeatmapLegend")
        drawDendrogramLegend(consoleInfo = "init output$calcPlotDendrogramLegend")
        state$anyPlotDrawn <<- TRUE
      }
      
      ## state
      state$showHCAplotPanel <<- TRUE
      state$plotHcaShown <<- TRUE
      
      ##########################
      ## update info and tip
      output$information <- renderText({
        print(paste("update output$information clear"))
        paste("", sep = "")
      })
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste(
          "Hover or select a cluster node or leaf node to view information about the corresponding MS\u00B9 feature cluster or MS\u00B9 feature respectively.", 
          "Brush horizontally and double-click to zoom in.", 
          "Double-click to zoom out.", 
          sep = "\n"
        )
      })
    })
    obsDrawPCA <- observeEvent(input$drawPCAplots, {
      #################################################
      ## get input
      drawPlots <- as.numeric(input$drawPCAplots)
      
      print(paste("Observe draw PCA plots", drawPlots))
      
      ## check if button was hit
      if(drawPlots == drawPCAButtonValue)
        return()
      drawPCAButtonValue <<- drawPlots
      
      pcaScaling      <- input$pcaScaling
      pcaLogTransform <- input$pcaLogTransform
      pcaDimensionOne <<- as.numeric(input$pcaDimensionOne)
      pcaDimensionTwo <<- as.numeric(input$pcaDimensionTwo)
      
      #################################################
      ## calc PCA
      pca <- calculatePCA(dataList = dataList, filterObj = filterPca, scaling = pcaScaling, logTransform = pcaLogTransform)
      
      pcaDataList <<- list()
      pcaDataList$pcaObj <<- pca
      pcaDataList$pcaScoresX <<- pca$scores[, pcaDimensionOne]
      pcaDataList$pcaScoresY <<- pca$scores[, pcaDimensionTwo]
      pcaDataList$pcaLoadingsX <<- pca$loadings[, pcaDimensionOne]
      pcaDataList$pcaLoadingsY <<- pca$loadings[, pcaDimensionTwo]
      pcaDataList$dimensionOne <<- pcaDimensionOne
      pcaDataList$dimensionTwo <<- pcaDimensionTwo
      
      ##########################
      ## hca selections
      if(!is.null(listForTable_Fragment_HCA)){ ## selection from fragment
        fragmentIndex <- which(dataList$fragmentMasses == fragmentsX[[selectionFragmentSelectedFragmentIndex]])
        precursorSet  <- which(dataList$featureMatrix[, fragmentIndex] != 0)
        selectionByFragmentInitPca(precursorSet)
      } else {
        selectionByFragmentReset()
      }
      if(!is.null(listForTable_Analysis_HCA)){ ## selection from HCA
        precursorSet <- getPrecursorSetFromTreeSelections(clusterDataList = clusterDataList, clusterLabels = selectionAnalysisTreeNodeSet)
        selectionByAnalysisInitPca(precursorSet)
      } else {
        selectionByAnalysisReset()
      }
      if(!is.null(filterSearch)){ ## selection from search
        selectionBySearchInitPca(filterSearch$filter)
      } else {
        selectionBySearchReset()
      }
      
      ##########################
      ## reset MS2 stuff
      if(!state$showHCAplotPanel){
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
      }
      
      ##########################
      ## draw
      resetPcaPlotRange()
      drawPcaPlots(consoleInfo = "drawPCA output$plotPcaScores")
      drawMS2Plot(consoleInfo = "drawPCA output$plotMS2")
      drawAnnotationLegendPCA(consoleInfo = "init output$plotAnnoLegend")
      drawScoresGroupsLegend(consoleInfo = "init output$plotScoresGroupsLegend")
      
      if(!state$anyPlotDrawn){
        drawMS2Legend(consoleInfo = "init output$ms2LegendPlot")
        drawFragmentDiscriminativityLegend(consoleInfo = "init output$plotFragmentDiscriminativityLegend")
        drawHeatmapLegend(consoleInfo = "init output$plotHeatmapLegend")
        drawDendrogramLegend(consoleInfo = "init output$calcPlotDendrogramLegend")
        state$anyPlotDrawn <<- TRUE
      }
      
      ## state
      state$showPCAplotPanel <<- TRUE
      state$plotPcaShown <<- TRUE
      
      ##########################
      ## update info and tip
      output$information <- renderText({
        print(paste("update output$information clear"))
        paste("", sep = "")
      })
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Hover a score node in the scores plot or a loadings node in the loadings plot to view information about the corresponding sample or MS\u00B9 feature respectively.", "Brush and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    ## listen to dendrogram/heatmap plot mouse events
    obsDendrogramHover <- observeEvent(input$plotDendrogram_hover, {
      hoverX <- input$plotDendrogram_hover$x
      hoverY <- input$plotDendrogram_hover$y
      
      plotWidth  <- session$clientData$output_plotDendrogram_width
      plotHeight <- session$clientData$output_plotDendrogram_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramPlotRangeY
      )
      print(paste("Observe dendrogram hover", minimumIndex))
      if(is.null(minimumIndex)){
        if(!is.null(fragmentsXhovered)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
        }
      } else {
        minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
        
        if(all(!is.null(selectionAnalysisTreeNodeSet), minimumLabel == selectionAnalysisTreeNodeSet)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
          fragmentsColorHovered <<- NULL
          fragmentsDiscriminativityHovered <<- NULL
          #output$information <- renderText({
          #  print(paste("update output$information", resultObj$infoText))
          #  paste("", sep = "")
          #})
        } else {
          #################################################
          ## fetch ms2 spectrum
          resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
          fragmentsXhovered <<- resultObj$fragmentMasses
          fragmentsYhovered <<- resultObj$fragmentAbundances
          fragmentsColorHovered <<- resultObj$fragmentColor
          fragmentsDiscriminativityHovered <<- resultObj$fragmentDiscriminativity
          
          #################################################
          ## output as message
          output$information <- renderText({
            print(paste("update output$information", resultObj$infoText))
            paste(resultObj$infoText, sep = "")
          })
        }
      }
      
      ## MS2 plot
      drawMS2Plot(consoleInfo = "dendrogram hover output$plotMS2")
      
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Hover a cluster node or leaf node to view information about the corresponding MS\u00B9 feature cluster or MS\u00B9 feature respectively.", "Brush horizontally and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    obsDendrogramClick <- observeEvent(input$plotDendrogram_click, {
      clickX <- input$plotDendrogram_click$x
      clickY <- input$plotDendrogram_click$y
      
      brush <- input$plotDendrogram_brush
      
      plotWidth  <- session$clientData$output_plotDendrogram_width
      plotHeight <- session$clientData$output_plotDendrogram_height
      
      if(is.null(clickX) | is.null(clickY))
        return()
      if(!is.null(brush))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramPlotRangeY
      )
      print(paste("Observe dendrogram click", is.null(minimumIndex), minimumIndex))
      
      if(is.null(minimumIndex)){
        ## reset stuff
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        selectionByAnalysisReset()
        selectionByFragmentReset()
      } else {
        ## tree selection
        minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
        #minimumText <- clusterDataList$poiText[[minimumIndex]]
        
        ## fetch ms2 spectrum
        resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
        
        ## keep fragment selection
        selectionFragmentSelectedFragmentIndexNew <- NULL
        if(!is.null(selectionFragmentSelectedFragmentIndex)){
          fragmentMass <- fragmentsX[[selectionFragmentSelectedFragmentIndex]]
          if(fragmentMass %in% resultObj$fragmentMasses)
            selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
        }
        
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsColor <<- resultObj$fragmentColor
        fragmentsDiscriminativity <<- resultObj$fragmentDiscriminativity
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        #################################################
        ## output as message
        selectionByHca(minimumLabel)
        
        ## update the selected fragment in case of overlapping spectra
        if(!is.null(selectionFragmentSelectedFragmentIndexNew))
          selectionByFragment(selectionFragmentSelectedFragmentIndexNew)
        else
          selectionByFragmentReset()
        
        output$information <- renderText({
          print(paste("update output$information", resultObj$infoText))
          paste(resultObj$infoText, sep = "")
        })
      }
      
      #################################################
      ## plots
      
      ## cluster dendrogram
      drawDendrogramPlot(consoleInfo = "dendrogram click output$plotDendrogram")
      
      ## MS2 plot
      resetMS2PlotRange()
      drawMS2Plot(consoleInfo = "dendrogram click output$plotMS2")
      
      if(state$showPCAplotPanel)
        ## update PCA plots
        drawPcaPlots(consoleInfo = "dendrogram click output$plotPcaScores")
    })
    obsDendrogramdblClick <- observeEvent(input$plotDendrogram_dblclick, {
      brush <- input$plotDendrogram_brush
      
      print(paste("observe dendrogram dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        dendrogramPlotRange$xMin <<- brush$xmin
        dendrogramPlotRange$xMax <<- brush$xmax
        dendrogramPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        dendrogramPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
      } else {
        ## reset range
        dendrogramPlotRange$xMin <<- 1
        dendrogramPlotRange$xMax <<- filterHca$numberOfPrecursorsFiltered
        dendrogramPlotRange$xInterval <<- c(1, filterHca$numberOfPrecursorsFiltered)
        dendrogramPlotRange$xIntervalSize <<- filterHca$numberOfPrecursorsFiltered - 1
      }
    })
    obsDendrogramRangeUpdate <- observe({
      xInterval <- dendrogramPlotRange$xInterval
      
      if(!state$showHCAplotPanel)
        return()
      
      print(paste("observe dendrogramPlotRange", xInterval))
      
      #drawDendrogramPlot(consoleInfo = "update range output$plotDendrogram", withHeatmap = TRUE)
    })
    obsHeatmaphover <- observeEvent(input$plotHeatmap_hover, {
      hoverX <- input$plotHeatmap_hover$x
      hoverY <- input$plotHeatmap_hover$y
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(hoverX < 0.5 | hoverX > (clusterDataList$numberOfPrecursorsFiltered + 0.5))
        return()
      if(hoverY < 0 | hoverY > 3)
        return()
      
      print(paste("Observe heatmap hover", hoverX, hoverY))
      
      #################################################
      ## info
      treeLeafIndex2 <- as.numeric(format(x = hoverX, digits = 0))
      treeLeafIndex  <- clusterDataList$cluster$order[[treeLeafIndex2]]
      precursorIndex <- filterHca$filter[[treeLeafIndex]]
      
      msg <- list()
      msg[[length(msg) + 1]] <- "MS\u00B9 feature: "
      msg[[length(msg) + 1]] <- dataList$precursorLabels[[precursorIndex]]
      msg[[length(msg) + 1]] <- "\n"
      
      if(hoverY > 2){
        ## lcf
        msg[[length(msg) + 1]] <- paste("log fold change [ log_2( mean(group ", filterHca$groups[[2]], ") / mean(group ", filterHca$groups[[1]], ") ) ]: ", sep = "")
        lfc <- dataList$dataFrameMeasurements[precursorIndex, dataList$lfcColumnNameFunctionFromName(filterHca$groups[[1]], filterHca$groups[[2]])]
        msg[[length(msg) + 1]] <- as.numeric(format(x = lfc, digits = 2))
      } else {
        if(hoverY > 1){
          ## group 1
          groupHere <- filterHca$groups[[1]]
        } else {
          ## group 2
          groupHere <- filterHca$groups[[2]]
        }
        msg[[length(msg) + 1]] <- paste("Mean abundance of group ", groupHere, ": ", sep = "")
        val <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupHere)]
        msg[[length(msg) + 1]] <- as.numeric(format(x = val, digits = 2))
        msg[[length(msg) + 1]] <- " = mean("
        vals <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataColumnsNameFunctionFromName(groupHere)]
        vals <- as.numeric(format(x = vals, digits = 2))
        msg[[length(msg) + 1]] <- paste(vals, collapse = ", ")
        msg[[length(msg) + 1]] <- ")"
      }
      
      output$information <- renderText({
        print(paste("update output$information heatmap hover", sep = ""))
        paste(msg, collapse = "")
      })
    })
    ## listen to MS2 plot mouse events
    obsMS2hover <- observeEvent(input$plotMS2_hover, {
      hoverX <- input$plotMS2_hover$x
      hoverY <- input$plotMS2_hover$y
      plotWidth  <- session$clientData$output_plotMS2_width
      plotHeight  <- session$clientData$output_plotMS2_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(any(is.null(fragmentsX), length(fragmentsX) == 0))
        return()
      
      ################################################
      ## decide whether the click is close enough to trigger event
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = fragmentsX, poiCoordinatesY = fragmentsY, 
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = ms2PlotRange$xIntervalSize, plotRangeY = 1
      )
      print(paste("Observe MS2 hover", minimumIndex))
      if(is.null(minimumIndex)){
        ## nothing in selection range
        #output$information <- renderText({
        #  print(paste("update output$information"))
        #  paste("", sep = "")
        #})
      } else {
        ## point selected
        fragmentIndex <- which(dataList$fragmentMasses == fragmentsX[[minimumIndex]])
        numberOfPrecursors <- sum(dataList$featureMatrix[clusterDataList$filter, fragmentIndex] != 0)
        output$information <- renderText({
          print(paste("update output$information"))
          paste("Fragment with m/z = ", fragmentsX[[minimumIndex]], " and (average) abundance = ", format(x = fragmentsY[[minimumIndex]], digits = 0, nsmall = 4), " is present in ", numberOfPrecursors, " MS/MS spectra.", sep = "")
        })
      }
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste(
          "Hover or click a fragment node (only 'Fragments from selection') to view information about this fragment.", 
          "Brush horizontally and double-click to zoom in.", 
          "Double-click to zoom out.", 
          sep = "\n"
        )
      })
    })
    obsMS2click <- observeEvent(input$plotMS2_click, {
      clickX <- input$plotMS2_click$x
      clickY <- input$plotMS2_click$y
      
      brush  <- input$plotMS2_brush
      
      plotWidth  <- session$clientData$output_plotMS2_width
      plotHeight <- session$clientData$output_plotMS2_height
      
      if(is.null(clickX) | is.null(clickY))
        return()
      if(!is.null(brush))
        ## ongoing brushing
        return()
      if(any(is.null(fragmentsX), length(fragmentsX) == 0))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = fragmentsX, poiCoordinatesY = fragmentsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = ms2PlotRange$xIntervalSize, plotRangeY = 1
      )
      print(paste("Observe MS2 click", minimumIndex))
      
      if(is.null(minimumIndex)){
        ##########################################
        ## reset click
        selectionByFragmentReset()
      } else {
        ##########################################
        ## peak click
        selectionByFragment(minimumIndex)
      }
      
      ## HCA
      if(state$showHCAplotPanel)
        drawDendrogramPlot(consoleInfo = "MS2 click output$plotDendrogram")
      ## PCA
      if(state$showPCAplotPanel)
        drawPcaLoadingsPlot(consoleInfo = "MS2 click output$plotPcaLoadings")
      
      ## update node selection
      drawMS2Plot(consoleInfo = "MS2 click output$plotMS2")
    })
    obsMS2dblClick <- observeEvent(input$plotMS2_dblclick, {
      brush <- input$plotMS2_brush
      
      if(all(any(is.null(fragmentsX), length(fragmentsX) == 0), any(is.null(fragmentsXhovered), length(fragmentsXhovered) == 0)))
        return()
      
      print(paste("observe MS2 dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        ms2PlotRange$xMin <<- brush$xmin
        ms2PlotRange$xMax <<- brush$xmax
        ms2PlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        ms2PlotRange$xIntervalSize <<- brush$xmax - brush$xmin
      } else {
        ## reset range
        resetMS2PlotRange()
      }
    })
    obsMS2rangeUpdate <- observe({
      xInterval <- ms2PlotRange$xInterval
      
      if(!(state$showHCAplotPanel | state$showPCAplotPanel))
        return()
      
      print(paste("observe ms2PlotRange", paste(xInterval, collapse = ", ")))
      
      #drawMS2Plot(consoleInfo = "range update output$plotMS2")
    })
    ## listen to PCA events
    obsPCAscoresHover <- observeEvent(input$plotPcaScores_hover, {
      hoverX <- input$plotPcaScores_hover$x
      hoverY <- input$plotPcaScores_hover$y
      plotWidth  <- session$clientData$output_plotPcaScores_width
      plotHeight <- session$clientData$output_plotPcaScores_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA scores hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaScoresX, poiCoordinatesY = pcaDataList$pcaScoresY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaScoresPlotRange$xIntervalSize, plotRangeY = pcaScoresPlotRange$yIntervalSize
      )
      print(paste("Observe PCA scores hover", hoverX, hoverY, minimumIndex))
      
      if(is.null(minimumIndex)){
        #output$information <- renderText({
        #  print(paste("update output$information PCA scores hover", sep = ""))
        #  paste("", group, sep = "")
        #})
      }
      else{
        dataColumnName <- dataList$dataColumnsNameFunctionFromNames(filterPca$groups)[[minimumIndex]]
        group <- dataList$groupNameFunctionFromDataColumnName(dataColumnName)
        output$information <- renderText({
          print(paste("update output$information PCA scores hover", sep = ""))
          paste("Sample '", dataColumnName , "' is a replicate of group '", group, "'.", sep = "")
        })
      }
      
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Hover a scores node to view information about the sample.", "Brush and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    obsPCAscoresDblClick <- observeEvent(input$plotPcaScores_dblclick, {
      brush <- input$plotPcaScores_brush
      
      print(paste("observe PCAscores dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        pcaScoresPlotRange$xMin <<- brush$xmin
        pcaScoresPlotRange$xMax <<- brush$xmax
        pcaScoresPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        pcaScoresPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
        pcaScoresPlotRange$yMin <<- brush$ymin
        pcaScoresPlotRange$yMax <<- brush$ymax
        pcaScoresPlotRange$yInterval <<- c(brush$ymin, brush$ymax)
        pcaScoresPlotRange$yIntervalSize <<- brush$ymax - brush$ymin
      } else {
        ## reset range
        minX <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
        maxX <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
        minY <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
        maxY <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
        
        pcaScoresPlotRange$xMin <<- minX
        pcaScoresPlotRange$xMax <<- maxX
        pcaScoresPlotRange$xInterval <<- c(minX, maxX)
        pcaScoresPlotRange$xIntervalSize <<- maxX - minX
        pcaScoresPlotRange$yMin <<- minY
        pcaScoresPlotRange$yMax <<- maxY
        pcaScoresPlotRange$yInterval <<- c(minY, maxY)
        pcaScoresPlotRange$yIntervalSize <<- maxY - minY
      }
    })
    obsPCAscoresRangeUpdate <- observe({
      xInterval <- pcaScoresPlotRange$xInterval
      yInterval <- pcaScoresPlotRange$yInterval
      
      if(!state$showPCAplotPanel)
        return()
      
      print(paste("observe pcaScoresPlotRange", paste(xInterval, collapse = ", "), paste(yInterval, collapse = ", ")))
      
      ## plot PCA
      #drawPcaScoresPlot(consoleInfo = "range update output$plotPcaScores")
    })
    obsPCAloadingsClick <- observeEvent(input$plotPcaLoadings_click, {
      clickX <- input$plotPcaLoadings_click$x
      clickY <- input$plotPcaLoadings_click$y
      plotWidth  <- session$clientData$output_plotPcaLoadings_width
      plotHeight <- session$clientData$output_plotPcaLoadings_height
      brush <- input$plotPcaLoadings_brush
      
      if(!is.null(brush))
        return()
      if(is.null(clickX) | is.null(clickY))
        return()
      print(paste("Observe PCA Loadings click", clickX, clickY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = pcaDataList$pcaLoadingsX, poiCoordinatesY = pcaDataList$pcaLoadingsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaLoadingsPlotRange$xIntervalSize, plotRangeY = pcaLoadingsPlotRange$yIntervalSize
      )
      print(paste("Observe PCA Loadings hover", clickX, clickY, minimumIndex))
      
      if(is.null(minimumIndex)){
        ## reset stuff
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        selectionByAnalysisReset()
        selectionByFragmentReset()
      } else {
        ## loadng selection
        precursorIndex <- filterPca$filter[[minimumIndex]]
        ## fetch ms2 spectrum
        resultObj <- getMS2spectrumOfPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        
        ## keep fragment selection
        selectionFragmentSelectedFragmentIndexNew <- NULL
        if(!is.null(selectionFragmentSelectedFragmentIndex)){
          fragmentMass <- fragmentsX[[selectionFragmentSelectedFragmentIndex]]
          if(fragmentMass %in% resultObj$fragmentMasses)
            selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
        }
        
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsColor <<- resultObj$fragmentColor
        fragmentsDiscriminativity <<- resultObj$fragmentDiscriminativity
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        selectionByPca(minimumIndex)
        
        ## update the selected fragment in case of overlapping spectra
        if(!is.null(selectionFragmentSelectedFragmentIndexNew))
          selectionByFragment(selectionFragmentSelectedFragmentIndexNew)
        else
          selectionByFragmentReset()
      }
      
      ## MS2 plot
      resetMS2PlotRange()
      drawMS2Plot(consoleInfo = "PCA loadings click output$plotMS2")
      drawPcaLoadingsPlot(consoleInfo = "PCA loadings click output$plotPcaLoadings")
      
      if(state$showHCAplotPanel)
        ## update dendrogram plot
        drawDendrogramPlot(consoleInfo = "PCA loadings click output$plotDendrogram")
    })
    obsPCAloadingsHover <- observeEvent(input$plotPcaLoadings_hover, {
      hoverX <- input$plotPcaLoadings_hover$x
      hoverY <- input$plotPcaLoadings_hover$y
      plotWidth  <- session$clientData$output_plotPcaLoadings_width
      plotHeight <- session$clientData$output_plotPcaLoadings_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA Loadings hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaLoadingsX, poiCoordinatesY = pcaDataList$pcaLoadingsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaLoadingsPlotRange$xIntervalSize, plotRangeY = pcaLoadingsPlotRange$yIntervalSize
      )
      if(is.null(minimumIndex)){
        fragmentsXhovered <<- resultObj$fragmentMasses
        fragmentsYhovered <<- resultObj$fragmentAbundances
        #output$information <- renderText({
        #  print(paste("update output$information PCA Loadings hover ", precursorIndex, sep = ""))
        #  paste("", sep = "")
        #})
      } else {
        print(paste("Observe PCA Loadings hover", hoverX, hoverY, minimumIndex))
        
        #################################################
        ## fetch ms2 spectrum
        precursorIndex <- filterPca$filter[[minimumIndex]]
        resultObj <- getMS2spectrumOfPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        if(all(!is.null(selectionAnalysisPcaLoadingSet), minimumIndex == selectionAnalysisPcaLoadingSet)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
          fragmentsColorHovered <<- NULL
          fragmentsDiscriminativityHovered <<- NULL
        } else {
          fragmentsXhovered <<- resultObj$fragmentMasses
          fragmentsYhovered <<- resultObj$fragmentAbundances
          fragmentsColorHovered <<- resultObj$fragmentColor
          fragmentsDiscriminativityHovered <<- resultObj$fragmentDiscriminativity
        }
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings hover ", precursorIndex, sep = ""))
          paste(resultObj$infoText, sep = "")
        })
      }
      
      drawMS2Plot(consoleInfo = "loadings hover output$plotMS2")
      
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste(
          "Hover or click a loadings node to view information about the corresponding MS\u00B9 feature.", 
          "Brush and double-click to zoom in.", 
          "Double-click to zoom out.", 
          sep = "\n"
        )
      })
    })
    obsPCAloadingsDblClick <- observeEvent(input$plotPcaLoadings_dblclick, {
      brush <- input$plotPcaLoadings_brush
      
      print(paste("observe PCAloadings dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        pcaLoadingsPlotRange$xMin <<- brush$xmin
        pcaLoadingsPlotRange$xMax <<- brush$xmax
        pcaLoadingsPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        pcaLoadingsPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
        pcaLoadingsPlotRange$yMin <<- brush$ymin
        pcaLoadingsPlotRange$yMax <<- brush$ymax
        pcaLoadingsPlotRange$yInterval <<- c(brush$ymin, brush$ymax)
        pcaLoadingsPlotRange$yIntervalSize <<- brush$ymax - brush$ymin
      } else {
        ## reset range
        minX <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
        maxX <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
        minY <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
        maxY <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
        
        pcaLoadingsPlotRange$xMin <<- minX
        pcaLoadingsPlotRange$xMax <<- maxX
        pcaLoadingsPlotRange$xInterval <<- c(minX, maxX)
        pcaLoadingsPlotRange$xIntervalSize <<- maxX - minX
        pcaLoadingsPlotRange$yMin <<- minY
        pcaLoadingsPlotRange$yMax <<- maxY
        pcaLoadingsPlotRange$yInterval <<- c(minY, maxY)
        pcaLoadingsPlotRange$yIntervalSize <<- maxY - minY
      }
    })
    obsPCAloadingsRangeUpdate <- observe({
      xInterval <- pcaLoadingsPlotRange$xInterval
      yInterval <- pcaLoadingsPlotRange$yInterval
      
      if(!state$showPCAplotPanel)
        return()
      
      print(paste("observe pcaLoadingsPlotRange ", paste(xInterval, collapse = ", "), "; ", paste(yInterval, collapse = ", "), sep = ""))
      
      ## plot PCA
      #drawPcaLoadingsPlot(consoleInfo = "range update output$plotPcaLoadings")
    })
    ## fragment plot
    obsFragmentPlotdblClick <- observeEvent(input$fragmentPlot_dblclick, {
      brush <- input$fragmentPlot_brush
      
      print(paste("observe fragmentPlot dblclick", brush))
      
      if (!is.null(brush)) {
        ## set range
        min <- brush$xmin
        max <- brush$xmax
      } else {
        ## reset range
        min <- min(dataList$masses)
        max <- max(dataList$masses)
      }
      
      fragmentPlotRange$xMin <<- min
      fragmentPlotRange$xMax <<- max
      fragmentPlotRange$xInterval <<- c(min, max)
      fragmentPlotRange$xIntervalSize <<- max - min
    })
    obsFragmentPlotRangeUpdate <- observe({
      xInterval <- fragmentPlotRange$xInterval
      
      print(paste("observe fragmentPlotRange", xInterval))
      
      #drawFragmentPlot(consoleInfo = "update range output$fragmentPlotDendrogram")
    })
    ## listen to annotation events
    obsSetPresentAnnoPrimary <- observeEvent(input$setPresentAnnotationPrimary, {
      value     <- input$presentAnnotationValue
      
      drawPlots <- as.numeric(input$setPresentAnnotationPrimary)
      if(drawPlots == setPresentAnnotationPrimaryValue)
        return()
      setPresentAnnotationPrimaryValue <<- drawPlots
      print(paste("Observe setPresentAnnotationPrimary", drawPlots))
      
      setAnnotationPrimary(precursorSet = selectedPrecursorSet, annotationValue = value)
    })
    obsRemovePresentAnno <- observeEvent(input$removePresentAnnotation, {
      value     <- input$presentAnnotationValue
      
      drawPlots <- as.numeric(input$removePresentAnnotation)
      ## why?
      #if(drawPlots == removePresentAnnotationValue)
      #  return()
      #removePresentAnnotationValue <<- drawPlots
      print(paste("Observe removePresentAnnotation", drawPlots))
      
      removeAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value)
    })
    obsToggleAddNewAnnoButton <- observeEvent(input$newAnnotationValue, {
      value <- input$newAnnotationValue
      
      print(paste("Observe newAnnotationValue", nchar(value)))
      
      if(nchar(value) > 0)
        shinyjs::enable("submitNewAnnotation")
        #enableActionButton(session, "submitNewAnnotation")
      else
        shinyjs::disable("submitNewAnnotation")
        #disableActionButton(session, "submitNewAnnotation")
    })
    obsAddNewAnno <- observeEvent(input$submitNewAnnotation, {
      value <- input$newAnnotationValue
      color <- input$newAnnotationColor
      #color <- "grey"
      
      drawPlots <- as.numeric(input$submitNewAnnotation)
      print(paste("Observe submitNewAnnotation", drawPlots))
      ## XXX why?
      #if(drawPlots == submitNewAnnotationValue)
      #  return()
      #submitNewAnnotationValue <<- drawPlots
      
      addAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value, annotationColor = color)
    })
    obsAddPresentAnno <- observeEvent(input$submitPreviousAnnotation, {
      value <- input$previousAnnotationValue
      
      drawPlots <- as.numeric(input$submitPreviousAnnotation)
      if(drawPlots == submitPreviousAnnotationValue)
        return()
      submitPreviousAnnotationValue <<- drawPlots
      print(paste("Observe submitPreviousAnnotation", drawPlots))
      
      color <- dataList$annoPresentColorsList[[match(x = value, table = dataList$annoPresentAnnotationsList)]]
      addAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value, annotationColor = color)
    })
    obsIgnoreValueChanged <- observeEvent(input$updateArtifactsFromCheckboxes, {
      updateArtifactsFromCheckboxes <- as.numeric(input$updateArtifactsFromCheckboxes)
      
      #################################################
      ## check if button was hit
      if(updateArtifactsFromCheckboxes == updateArtifactsFromCheckboxesButtonValue)
        return()
      updateArtifactsFromCheckboxesButtonValue <<- updateArtifactsFromCheckboxes
      
      ## get and process values
      vals <- getInputValues(id = paste(state$selectedSelection, "Ignore", sep = "_"), counter = selectedTable_id, len = nrow(selectedTable))
      
      if(all(is.na(vals)))
        return()
      
      nas <- is.na(vals)
      vals[nas] <- dataList$annoArrayIsArtifact[selectedPrecursorSet][nas]
      vals <- as.logical(vals)
      
      setArtifactState(selectedPrecursorSet, vals)
    })
    ## export image buttons
    obsShowHCAplotPanel <- observe({
      print(paste("observe state$showHCAplotPanel", state$showHCAplotPanel))
      
      if(state$showHCAplotPanel){
        shinyjs::enable("downloadHcaImage")
      } else {
        shinyjs::disable("downloadHcaImage")
      }
    })
    obsShowPCAplotPanel <- observe({
      print(paste("observe state$showPCAplotPanel", state$showPCAplotPanel))
      
      if(state$showPCAplotPanel){
        shinyjs::enable("downloadPcaImage")
      } else {
        shinyjs::disable("downloadPcaImage")
      }
    })
    ## suspend observer
    session$onSessionEnded(function() {
      ## sidepanel
      obsTabs$suspend()
      obsChangePlot$suspend()
      obsChangeSelection$suspend()
      obsClearSelection$suspend()
      obsFile$suspend()
      obsLoadProjectData$suspend()
      #obsExampleDataSelection$suspend()
      obsLoadExampleData$suspend()
      obsImportMs1DataFile$suspend()
      obsImportMs2DataFile$suspend()
      obsImportMs1Ms2Data$suspend()
      obsApplyImportParameterFile$suspend()
      obsFileInputSelection$suspend()
      obsShowSideBar$suspend()
      obsShowClusterLabels$suspend()
      obsShowScoresLabels$suspend()
      #observePcaLoadingsProperties$suspend()
      obsShowLoadingsLabels$suspend()
      obsShowLoadingsAbundance$suspend()
      ## filter
      obsClearSearch$suspend()
      obsApplySearch$suspend()
      obsApplyGlobalMS2filters$suspend()
      obsApplyHcaFilters$suspend()
      obsApplyPcaFilters$suspend()
      observeGroupSet$suspend()
      ## draw
      obsDrawHCA$suspend()
      obsDrawPCA$suspend()
      ## dendrogram
      obsDendrogramHover$suspend()
      obsDendrogramClick$suspend()
      obsDendrogramdblClick$suspend()
      obsDendrogramRangeUpdate$suspend()
      ## heatmap
      obsHeatmaphover$suspend()
      ## MS2
      obsMS2hover$suspend()
      obsMS2click$suspend()
      obsMS2dblClick$suspend()
      obsMS2rangeUpdate$suspend()
      ## PCA
      obsPCAscoresHover$suspend()
      obsPCAscoresDblClick$suspend()
      obsPCAscoresRangeUpdate$suspend()
      obsPCAloadingsClick$suspend()
      obsPCAloadingsHover$suspend()
      obsPCAloadingsDblClick$suspend()
      obsPCAloadingsRangeUpdate$suspend()
      ## fragment plot
      obsFragmentPlotdblClick$suspend()
      obsFragmentPlotRangeUpdate$suspend()
      ## anno
      obsSetPresentAnnoPrimary$suspend()
      obsRemovePresentAnno$suspend()
      obsToggleAddNewAnnoButton$suspend()
      obsAddNewAnno$suspend()
      obsAddPresentAnno$suspend()
      obsIgnoreValueChanged$suspend()
      obsShowHCAplotPanel$suspend()
      obsShowPCAplotPanel$suspend()
    })
    
    #########################################################################################
    #########################################################################################
    ## direct output rendering
    output$fileInfo <- renderText({
      print(paste("init output$fileInfo"))
      paste("Please select a project file and press 'Load project data'")
    })
    output$information <- renderText({
      print(paste("init output$information"))
      ""
    })
    output$rInfo <- renderText({
      print(paste("init rInfo"))
      R.Version()$version.string
    })
    
    #########################################################################################
    #########################################################################################
    ## direct output values
    output$showGUI <- reactive({
      print("update output$showGUI")
      output$information <- renderText({
        print(paste("init information", sep = ""))
        paste("Please perform ploting.", sep = "")
      })
      return(!is.null(state$importedOrLoadedFile_s_))
    })
    output$analysisType <- reactive({
      print(paste("reactive update analysisType", state$analysisType))
      if(!is.null(state$analysisType)){
        if(state$analysisType == "HCA")
          plotToShow <<- "Display HCA"
        if(state$analysisType == "PCA")
          plotToShow <<- "Display PCA"
      }
      return(state$analysisType)
    })
    output$showSideBar <- reactive({
      print(paste("reactive update showSideBar", state$showSideBar))
      return(state$showSideBar)
    })
    output$showHCAplotPanel <- reactive({
      print(paste("reactive update showHCAplotPanel", state$showHCAplotPanel))
      updateChangePlotRadioButton()
      return(state$showHCAplotPanel)
    })
    output$showPCAplotPanel <- reactive({
      print(paste("reactive update showPCAplotPanel", state$showPCAplotPanel))
      updateChangePlotRadioButton()
      return(state$showPCAplotPanel)
    })
    output$globalMS2filterValid <- reactive({
      print(paste("reactive update globalMS2filterValid", state$globalMS2filterValid))
      return(state$globalMS2filterValid)
    })
    output$hcaFilterValid <- reactive({
      print(paste("reactive update hcaFilterValid", state$hcaFilterValid))
      return(state$hcaFilterValid)
    })
    output$pcaFilterValid <- reactive({
      print(paste("reactive update pcaFilterValid", state$pcaFilterValid))
      return(state$pcaFilterValid)
    })
    output$precursorSetSelected <- reactive({
      print(paste("reactive update precursorSetSelected", state$precursorSetSelected))
      return(state$precursorSetSelected)
    })
    output$plotHcaShown <- reactive({
      print(paste("reactive update plotHcaShown", state$plotHcaShown))
      return(state$plotHcaShown)
    })
    output$plotPcaShown <- reactive({
      print(paste("reactive update plotPcaShown", state$plotPcaShown))
      return(state$plotPcaShown)
    })
    output$searchfilterValid <- reactive({
      print(paste("reactive update searchfilterValid", state$searchfilterValid))
      return(state$searchfilterValid)
    })
    output$filterSearchActive <- reactive({
      print(paste("reactive update filterSearchActive", state$filterSearchActive))
      return(state$filterSearchActive)
    })
    output$selectedSelection <- reactive({
      print(paste("reactive update selectedSelection", state$selectedSelection))
      return(state$selectedSelection)
    })
    
    updateChangePlotRadioButton <- function(){
      if(state$showHCAplotPanel & state$showPCAplotPanel & !is.null(state$analysisType)){
        if(state$analysisType == "HCA")
          selectedItem <- "Display HCA"
        else
          selectedItem <- "Display PCA"
        updateRadioButtons(session = session, inputId = "changePlot", selected = selectedItem)
      }
    }
    
    #########################################################################################
    #########################################################################################
    ## properties
    options(shiny.maxRequestSize=500*1024^2) ## 500 mb file size
    outputOptions(output, 'showGUI', suspendWhenHidden=FALSE)
    outputOptions(output, 'showSideBar', suspendWhenHidden=FALSE)
    outputOptions(output, 'showHCAplotPanel', suspendWhenHidden=FALSE)
    outputOptions(output, 'showPCAplotPanel', suspendWhenHidden=FALSE)
    outputOptions(output, 'analysisType', suspendWhenHidden=FALSE)
    outputOptions(output, 'globalMS2filterValid', suspendWhenHidden=FALSE)
    outputOptions(output, 'hcaFilterValid', suspendWhenHidden=FALSE)
    outputOptions(output, 'pcaFilterValid', suspendWhenHidden=FALSE)
    outputOptions(output, 'precursorSetSelected', suspendWhenHidden=FALSE)
    
    outputOptions(output, 'searchfilterValid', suspendWhenHidden=FALSE)
    outputOptions(output, 'filterSearchActive', suspendWhenHidden=FALSE)
    outputOptions(output, 'plotHcaShown', suspendWhenHidden=FALSE)
    outputOptions(output, 'plotPcaShown', suspendWhenHidden=FALSE)
    outputOptions(output, 'selectedSelection', suspendWhenHidden=FALSE)
    
    #########################################################################################
    #########################################################################################
    ## ui generation
    output$runRightColumn <- renderUI({
      print(paste("### GUI ### Generate right GUI"))
      column(width = state$runRightColumnWidth,
         #########################################################################################
         ## show side bar and change plot
         conditionalPanel(
           condition = "(output.showHCAplotPanel && output.analysisType) == 'HCA' || (output.showPCAplotPanel && output.analysisType == 'PCA')",
           #condition = "output.showHCAplotPanel || output.showPCAplotPanel",
           fluidRow(
             column(width = 6,
                div(style="float:right",
                    bsTooltip(id = "showSideBar", title = "Display or hide the side bar", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "showSideBar", label = "Display side bar", value = showSideBar)
                )
             ),##column
             column(width = 6,
                conditionalPanel(
                  condition = "output.showHCAplotPanel && output.showPCAplotPanel",
                  div(style="float:left",
                      bsTooltip(id = "changePlot", title = "Switch between HCA plots and PCA plots", placement = "bottom", trigger = "hover"),
                      radioButtons(inputId = "changePlot", label = NULL, choices = c("Display HCA", "Display PCA"), inline = TRUE, selected = plotToShow)
                  )
                )##conditional
             )##column
           )##row
         ),##conditional
         ##############################################################################################
         ## HCA plots
         fluidRow(
           column(width = 12-state$legendColumnWidth,
             conditionalPanel(
               condition = 'output.analysisType == "HCA" && output.showHCAplotPanel',
               fluidRow(
                  plotOutput(height = 500, 
                             outputId = "plotDendrogram", 
                             #hover    = "plotDendrogram_hover", 
                             hover    = hoverOpts(
                               id = "plotDendrogram_hover",
                               delay = 50, 
                               delayType = "debounce"
                             ),
                             click    = "plotDendrogram_click",
                             dblclick = "plotDendrogram_dblclick",
                             #brush    = "plotDendrogram_brush"
                             brush    = brushOpts(
                               id = "plotDendrogram_brush",
                               resetOnNew = TRUE,
                               direction = "x",
                               delay = 00,
                               delayType = "debounce"
                             )
                  ),
                  plotOutput(height = 75, 
                             outputId = "plotHeatmap",
                             #hover    = "plotHeatmap_hover", 
                             hover    = hoverOpts(
                               id = "plotHeatmap_hover",
                               delay = 50, 
                               delayType = "debounce"
                             )
                             #click = "plotHeatmap_click"
                  )
               )## row
             ),## conditional
             ##############################################################################################
             ## PCA plots
             conditionalPanel(
               condition = 'output.analysisType == "PCA" && output.showPCAplotPanel',
               fluidRow(
                 column(width = 6,
                    plotOutput(height = 500, 
                               outputId = "plotPcaScores", 
                               #hover    = "plotPcaScores_hover",
                               hover    = hoverOpts(
                                 id = "plotPcaScores_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "plotPcaScores_click",
                               dblclick = "plotPcaScores_dblclick",
                               #brush    = "plotPcaScores_brush"
                               brush = brushOpts(
                                 id = "plotPcaScores_brush",
                                 resetOnNew = TRUE,
                                 direction = "xy",
                                 delay = 00,
                                 delayType = "debounce"
                               )
                    )
                 ),## column
                 column(width = 6,
                    plotOutput(height = 500, 
                               outputId = "plotPcaLoadings", 
                               #hover    = "plotPcaLoadings_hover",
                               hover    = hoverOpts(
                                 id = "plotPcaLoadings_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "plotPcaLoadings_click",
                               dblclick = "plotPcaLoadings_dblclick",
                               #brush    = "plotPcaScores_brush"
                               brush = brushOpts(
                                 id = "plotPcaLoadings_brush",
                                 resetOnNew = TRUE,
                                 direction = "xy",
                                 delay = 00,
                                 delayType = "debounce"
                               )
                    )
                 )## column
               )## row
             )## conditional
           ),##column
           column(width = state$legendColumnWidth,
              conditionalPanel(## scores / loadings properties
                condition = 'output.analysisType == "PCA" && output.showPCAplotPanel',
                #splitLayout(
                #  style = "border: 1px solid silver;padding: 0px 6px;",
                #  div(
                    h5("PCA scores"),
                    bsTooltip(id = "showScoresLabels", title = "Display scores labels", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "showScoresLabels", label = "Show labels", value = TRUE),
                #  )
                #),
                #splitLayout(
                #  style = "border: 1px solid silver;padding: 0px 6px;",
                #  div(
                    h5("PCA loadings"),
                    bsTooltip(id = "showLoadingsLabels", title = "Display loadings labels", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "showLoadingsLabels", label = "Show labels", value = FALSE),
                    bsTooltip(id = "showLoadingsAbundance", title = "Use abundance in MS\u00B9 to scale the size of loadings nodes", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "showLoadingsAbundance", label = "Show abundance", value = FALSE)
                    #bsTooltip(id = "pcaLoadingsProperties", title = "Check to display loadings labels or to use abundance in MS\u00B9 to scale the size of loadings nodes", placement = "bottom", trigger = "hover"),
                    #checkboxGroupInput(inputId = "pcaLoadingsProperties", label = NULL, choices = c("Show labels", "Show abundance"))
                #  )
                #)
              ),
              conditionalPanel(## dendrogram properties
                condition = 'output.analysisType == "HCA" && output.showHCAplotPanel',
                #splitLayout(
                #  style = "border: 1px solid silver;padding: 0px 6px;",
                #  div(
                    h5("HCA dendrogram"),
                    bsTooltip(id = "showClusterLabels", title = "Display the labels of cluster nodes and MS\u00B9 feature nodes representing the number of characteristic fragments", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "showClusterLabels", label = "Show labels", value = TRUE)
                #  )
                #)
              ),
              conditionalPanel(
                condition = 'output.analysisType == "HCA" && output.showHCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotAnnoLegendHCA", height = state$annotationLegendHeightHca)
                )
              ),## conditional
              conditionalPanel(
                condition = 'output.analysisType == "PCA" && output.showPCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotAnnoLegendPCA", height = state$annotationLegendHeightPca)
                )
              ),## conditional
              conditionalPanel(
                condition = '(output.analysisType == "HCA" && output.showHCAplotPanel) || (output.analysisType == "PCA" && output.showPCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "calcPlotDendrogramLegend", height = 80)
                )
              ),## conditional
              conditionalPanel(
                condition = 'output.analysisType == "HCA" && output.showHCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotHeatmapLegend", height = 150)
                )
              ),## conditional
              conditionalPanel(## loadings properties
                condition = 'output.analysisType == "PCA" && output.showPCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotScoresGroupsLegend", height = state$scoresGroupsLegendHeight)
                )
              ),
              conditionalPanel(
                condition = '(output.analysisType == "HCA" && output.showHCAplotPanel) || (output.analysisType == "PCA" && output.showPCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotMS2Legend", height = 80)
                )
              ),## conditional
              conditionalPanel(
                condition = '(output.analysisType == "HCA" && output.showHCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotFragmentDiscriminativityLegend", height = 100)
                )
              )## conditional
           )## column
         ),## row
         #########################################################################################
         ## MS2 plot and info
         conditionalPanel(
           condition = '(output.showHCAplotPanel && output.analysisType == "HCA") || (output.showPCAplotPanel && output.analysisType == "PCA")',
           fluidRow(
             plotOutput(height = 250, 
                        outputId = "plotMS2",
                        #hover    = "plotMS2_hover",
                        hover    = hoverOpts(
                          id = "plotMS2_hover",
                          delay = 50, 
                          delayType = "debounce"
                        ),
                        click    = "plotMS2_click",
                        dblclick = "plotMS2_dblclick",
                        #brush    = "plotMS2_brush",
                        brush = brushOpts(
                          id = "plotMS2_brush",
                          resetOnNew = TRUE,
                          direction = "x",
                          delay = 00,
                          delayType = "debounce"
                        )
             ),
             ##############################################################################################
             ## infos
             wellPanel(
               h4("Information"),
               bsTooltip(id = "information", title = "Information about items in the plot", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("information"),
               h4("Tip"),
               bsTooltip(id = "tip", title = "Information about operating options", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("tip")
             )## well
           )## row
         ),## conditional
         ##############################################################################################
         ## precursor set selection and annotation
         ## change selection
         conditionalPanel(
           condition = '(output.showHCAplotPanel && output.analysisType == "HCA") || (output.showPCAplotPanel && output.analysisType == "PCA")',
           fluidRow(
             wellPanel(
               h4("MS\u00B9 feature selections"),
               bsTooltip(id = "changeSelection", title = "Switch MS\u00B9 feature selection", placement = "bottom", trigger = "hover"),
               radioButtons(inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), selected = changeSelectionCurrentSelection, inline = TRUE),
               bsTooltip(id = "selectionInfo", title = "The number of MS\u00B9 features in the current selection", placement = "bottom", trigger = "hover"),
               hr(),
               verbatimTextOutput("selectionInfo"),
               conditionalPanel(
                 condition = 'output.precursorSetSelected',
                 tabsetPanel(id = "precursorSelectionTabs",
                   tabPanel(title = precursorSelectionTabSelection, 
                     wellPanel(
                       ## selection infos
                       bsTooltip(id = "metFragLink", title = "Press to send the current MS\u00B9 feature as well as the corresponding MS/MS spectrum to MetFrag", placement = "bottom", trigger = "hover"),
                       htmlOutput(outputId = "metFragLink"),
                       bsTooltip(id = "downloadSelectedPrecursors", title = "Download a project file which is reduced to the selected set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
                       downloadButton('downloadSelectedPrecursors', 'Download reduced project file'),
                       bsTooltip(id = "clearSelection", title = "Press to clear this selection", placement = "bottom", trigger = "hover"),
                       actionButton(inputId = "clearSelection", label = "Clear selection", class="btn-success")
                     )## well
                   ),## tab
                   tabPanel(title = precursorSelectionTabAnnotation, 
                      wellPanel(
                        h4("Present annotation(s)"),
                        fluidRow(
                          column(
                            width = 3,
                            bsTooltip(id = "presentAnnotationValue", title = "The set of present annotations for the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            selectInput(inputId = "presentAnnotationValue", label = NULL, choices = c("[init]"), selectize = FALSE)
                          ),## column
                          column(
                            width = 3,
                            bsTooltip(id = "setPresentAnnotationPrimary", title = "Sets the selected annotation primary for the set of selected MS\u00B9 features; i.e. this annotation will be used preferentially for coloring in HCA and PCA", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "setPresentAnnotationPrimary", label = "Set primary", class="btn-success")
                          ),## column
                          column(
                            width = 6,
                            bsTooltip(id = "removePresentAnnotation", title = "Removes the selected annotation from the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "removePresentAnnotation", label = "Remove annotation", class="btn-success")
                          )## column
                        ),##row
                        fluidRow(
                          column(
                            width = 6,
                            h4("Add new annotation"),
                            bsTooltip(id = "newAnnotationValue", title = "The name of this annotation", placement = "bottom", trigger = "hover"),
                            textInput(inputId = "newAnnotationValue", label = "Type new annotation"),
                            #colourInput(inputId = "newAnnotationColor", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = c("blue")),
                            bsTooltip(id = "newAnnotationColor", title = "The color of this annotation", placement = "bottom", trigger = "hover"),
                            colourInput(inputId = "newAnnotationColor", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = colorPalette()),
                            bsTooltip(id = "submitNewAnnotation", title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "submitNewAnnotation", label = "Add new annotation", class="btn-success")
                          ),
                          column(
                            width = 6,
                            h4("Add previous annotation"),
                            bsTooltip(id = "previousAnnotationValue", title = "The set of annotations which have been assigned before", placement = "bottom", trigger = "hover"),
                            selectInput(inputId = "previousAnnotationValue", label = "Select previous annotation", choices = c("Artifact"), selectize = FALSE),
                            bsTooltip(id = "submitPreviousAnnotation", title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "submitPreviousAnnotation", label = "Add previous annotation", class="btn-success")
                          )
                        )
                      )## well
                   ),## tab
                   tabPanel(title = precursorSelectionTabTable, 
                      wellPanel(
                        h4("Selected MS\u00B9 features"),
                        bsTooltip(id = "updateArtifactsFromCheckboxes", title = "Adds the annotation \\'ignore\\' to the set of checked MS\u00B9 features", placement = "bottom", trigger = "hover"),
                        actionButton(inputId = "updateArtifactsFromCheckboxes", label = "Apply annotation 'Ignore' to MS\u00B9 features", class="btn-success"),
                        DT::dataTableOutput("table")
                      )## well
                   )## tab
                 )## tab set
               )## conditional
             )## well
           )## row
         )##conditional
      )##column
    })
    
    #########################################################################################
    #########################################################################################
    ## download
    timeStampForFiles <- function(){
      timeStamp <- gsub(" ", "_", gsub(":", ".", Sys.time()))
      return(timeStamp)
    }
    createImportParameterSetExportFileName <- function(){
      fileProjectName <- dataList$importParameterSet$projectName
      fileProjectName <- gsub(" ", "_", gsub(":", ".", fileProjectName))
      fileName <- paste(timeStampForFiles(), "_", fileProjectName, "_import_parameters.txt", sep = "")
      return(fileName)
    }
    createExportMatrixName <- function(){
      fileName <- paste(timeStampForFiles(), "_selectedPrecursorMatrix.csv.gz", sep = "")
      return(fileName)
    }
    createExportImageName <- function(item, extension){
      fileName <- paste(timeStampForFiles(), "_", item, ".", extension, sep = "")
      return(fileName)
    }
    createExportDistanceMatrixName <- function(distanceMeasure){
      fileName <- paste(timeStampForFiles(), "_distanceMatrix_", distanceMeasure, ".csv", sep = "")
      return(fileName)
    }
    createExportMatrix <- function(precursorSet){
      
      numberOfRows    <- length(precursorSet)
      numberOfColumns <- ncol(dataList$featureMatrix)
      
      ###########################################################
      ## built reduced MS2 matrix
      fragmentMatrix      <- dataList$featureMatrix[precursorSet, ]
      fragmentCounts      <- apply(X = fragmentMatrix, MARGIN = 2, FUN = function(x){ sum(x != 0) })
      fragmentIntensities <- apply(X = fragmentMatrix, MARGIN = 2, FUN = function(x){ sum(x) }) / fragmentCounts
      fragmentMasses      <- dataList$fragmentMasses
      
      fragmentSelection   <- fragmentCounts != 0
      
      fragmentMatrix      <- fragmentMatrix[, fragmentSelection]
      fragmentCounts      <- fragmentCounts[fragmentSelection]
      fragmentIntensities <- fragmentIntensities[fragmentSelection]
      fragmentMasses      <- fragmentMasses[fragmentSelection]
      
      ## fragment matrix
      dgTMatrix <- as(fragmentMatrix, "dgTMatrix")
      matrixRows <- dgTMatrix@i + 1
      matrixCols <- dgTMatrix@j + 1
      matrixVals <- dgTMatrix@x
      
      numberOfColumns2 <- ncol(fragmentMatrix)
      
      fragmentMatrix <- matrix(data = rep(x = "", times = numberOfRows * numberOfColumns2), nrow = numberOfRows, ncol = numberOfColumns2)
      fragmentMatrix[cbind(matrixRows, matrixCols)] <- matrixVals
      
      ## box
      ms2Matrix     <- rbind(
        fragmentCounts,
        fragmentIntensities,
        fragmentMasses,
        fragmentMatrix
      )
      
      ###########################################################
      ## built MS1 matrix
      ms1Matrix     <- rbind(
        dataList$dataFrameMS1Header,
        dataList$dataFrameInfos[precursorSet, ]
      )
      ms1Matrix     <- as.matrix(ms1Matrix)
      
      ###########################################################
      ## export annotations
      
      ## process annotations
      annotations <- dataList$annoArrayOfLists
      for(i in 1:length(annotations))
        if(dataList$annoArrayIsArtifact[[i]])
          annotations[[i]] <- c(annotations[[i]], dataList$annotationValueIgnore)
      
      annotationStrings <- vector(mode = "character", length = length(annotations))
      for(i in 1:length(annotations)){
        if(length(annotations[[i]]) > 0)
          annotationStrings[[i]] <- paste(annotations[[i]], sep = ", ")
        else
          annotationStrings[[i]] <- ""
      }
      annotationStrings <- annotationStrings[precursorSet]
      
      ## process annotaiotn-color-map
      annoPresentAnnotations <- dataList$annoPresentAnnotationsList[-1]
      annoPresentColors      <- dataList$annoPresentColorsList[-1]
      
      if(length(annoPresentAnnotations) > 0)
        annotationColors <- paste(annoPresentAnnotations, annoPresentColors, sep = "=", collapse = ", ")
      else
        annotationColors <- ""
      annotationColors <- paste(dataList$annotationColorsName, "={", annotationColors, "}", sep = "")
      
      ## box
      annotationColumn <- c("", annotationColors, dataList$annotationColumnName, annotationStrings)
      ms1Matrix[, dataList$annotationColumnIndex] <- annotationColumn
      
      ###########################################################
      ## assemble
      dataFrame <- cbind(
        ms1Matrix,
        ms2Matrix
      )
      
      return(dataFrame)
    }
    writeTable <- function(precursorSet, file){
      #compressed
      #[1] "create 22
      #[1] "write 23
      
      #uncompressed
      #[1] "create 23
      #[1] "write 19
      
      ## 02 with zeros
      #[1] "create 4
      #[1] "write 38
      
      ## sparse matrix
      #[1] "create 5
      #[1] "write 19
      
      ## matrix
      #[1] "create 6
      #[1] "write 12
      
      #print(paste("create", Sys.time()))
      dataFrame <- createExportMatrix(precursorSet)
      #print(paste("create ready", Sys.time()))
      gz1 <- gzfile(description = file, open = "w")
      #print(paste("write", Sys.time()))
      write.table(x = dataFrame, file = gz1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      #write.table(x = dataFrame, file = file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      #print(paste("write ready", Sys.time()))
      close(gz1)
    }
    ## individual downloads
    output$downloadGlobalMS2filteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterGlobal$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadHcaFilteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterHca$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadPcaFilteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterPca$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadSearchPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterPca$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadSelectedPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- selectedPrecursorSet
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadAllPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- 1:dataList$numberOfPrecursors
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    ## download selected
    output$downloadHcaSelectedPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        ## get selected precursors
        if(is.null(selectionAnalysisTreeNodeSet)){
          ## all precursors
          precursorSet <- filterHca$filter
        } else {
          precursorSet <- getPrecursorSetFromTreeSelections(clusterDataList = clusterDataList, clusterLabels = selectionAnalysisTreeNodeSet)
        }
        
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadImportParameterSet <- downloadHandler(
      filename = function() {
        createImportParameterSetExportFileName()
      },
      content = function(file) {
        fileLines <- serializeParameterSetFile(dataList$importParameterSet, toolName, toolVersion)
        writeLines(text = fileLines, con = file)
      },
      contentType = 'text/csv'
    )
    ## download images
    output$downloadHcaImage <- downloadHandler(
      filename = function() {
        fileType <- input$downloadHcaImageType
        createExportImageName("HCA", fileType)
      },
      content = function(file) {
        ## 1 den    ## 2 hea
        ## 3 ms2    ## 4 l anno
        ## 5 l sel  ## 6 l hea
        ## 7 l ms2
        ## 
        ## 1 4
        ## 1 5
        ## 1 6
        ## 2 7
        ## 3 7
        ## 
        
        fileType <- input$downloadHcaImageType
        
        ## parameters
        widthInInch     <- 10
        heigthInInch    <- 7.5
        resolutionInDPI <- 600
        widthInPixel    <- widthInInch  * resolutionInDPI
        heightInPixel   <- heigthInInch * resolutionInDPI
        
        #tmpFile <- paste("/home/htreutle/Downloads/MetSWATH/", createExportImageName("PCA"), sep = "")
        #print(tmpFile)
        #png(filename = tmpFile, 3000, 2000, res = 300, bg = "white")
        switch(fileType,
               "png"={
                 png(filename = file, width = widthInPixel, height = heightInPixel, res = resolutionInDPI, bg = "white")
               },
               "svg"={
                 svg(filename = file)
               },
               "pdf"={
                 pdf(file = file, title = "PCA image export from MetFam")
               },
               stop(paste("Unknown file type (", fileType, ")!", sep = ""))
        )
        
        layout(
          mat = matrix(
            data = c(1, 1, 1, 1, 2, 3,
                     4, 5, 6, 7, 8, 8), 
            nrow = 6, ncol = 2), 
          widths = c(4, 1), 
          heights = c(0.6, 1.4, 0.6, 0.6, 0.5, 1.5)
        )
        
        drawDendrogramPlotImpl()
        drawHeatmapPlotImpl()
        drawMS2PlotImpl()
        
        drawDendrogramLegendImpl()
        drawHeatmapLegendImpl()
        drawMS2LegendImpl()
        drawFragmentDiscriminativityLegendImpl()
        calcPlotAnnoLegendForImage(state$annotationsHca$setOfAnnotations, state$annotationsHca$setOfColors, 15)
        #drawAnnotationLegendImpl()
        
        dev.off()
        ## cannot open file 'Rplots.pdf'
      }#,
      #contentType = 'image/png'
    )
    output$downloadPcaImage <- downloadHandler(
      filename = function() {
        fileType <- input$downloadPcaImageType
        createExportImageName("PCA", fileType)
      },
      content = function(file) {
        ## 1 score  ## 2 loadings
        ## 3 ms2    ## 4 l anno
        ## 5 l sel  ## 6 l hea
        ## 7 l ms2
        ## 
        ## 1 2 4
        ## 1 2 5
        ## 1 2 6
        ## 1 2 7
        ## 3 3 7
        ## 
        
        fileType <- input$downloadPcaImageType
        
        ## parameters
        widthInInch     <- 10
        heigthInInch    <- 7.5
        resolutionInDPI <- 600
        widthInPixel    <- widthInInch  * resolutionInDPI
        heightInPixel   <- heigthInInch * resolutionInDPI
        
        #tmpFile <- paste("/home/htreutle/Downloads/MetSWATH/", createExportImageName("PCA"), sep = "")
        #print(tmpFile)
        #png(filename = tmpFile, 3000, 2000, res = 300, bg = "white")
        
        switch(fileType,
               "png"={
                 png(filename = file, width = widthInPixel, height = heightInPixel, res = resolutionInDPI, bg = "white")
               },
               "svg"={
                 svg(filename = file)
               },
               "pdf"={
                 pdf(file = file, title = "PCA image export from MetFam")
               },
               stop(paste("Unknown file type (", fileType, ")!", sep = ""))
        )
        
        layout(
          mat = matrix(
            data = c(1, 1, 1, 1, 1, 3,
                     2, 2, 2, 2, 2, 3, 
                     4, 5, 6, 7, 8, 8), 
            nrow = 6, ncol = 3), 
          widths = c(2, 2, 1), 
          heights = c(0.7, 0.6, 0.6, 0.6, 1.2, 1.5)
        )
        
        drawPcaScoresPlotImpl()
        drawPcaLoadingsPlotImpl()
        drawMS2PlotImpl()
        
        calcPlotScoresGroupsLegendForImage(state$scoresGroups$groups, state$scoresGroups$colors, 5)
        drawDendrogramLegendImpl()
        drawMS2LegendImpl()
        drawFragmentDiscriminativityLegendImpl()
        calcPlotAnnoLegendForImage(state$annotationsPca$setOfAnnotations, state$annotationsPca$setOfColors, 20)
        #drawAnnotationLegendImpl()
        
        dev.off()
      }#,
      #contentType = 'image/png'
    )
    output$downloadDistanceMatrix <- downloadHandler(
      filename = function() {
        createExportDistanceMatrixName(currentDistanceMatrixObj$distanceMeasure)
      },
      content = function(file) {
        write.table(x = currentDistanceMatrixObj$distanceMatrix, file = file, sep = "\t", row.names = TRUE, quote = FALSE, col.names=NA)
      },
      contentType = 'text/csv'
    )
    ## download publication data
    output$downloadMsData <- downloadHandler(
      filename = function() {
        return("Metabolite_profile_showcase.txt")
      },
      content = function(file) {
        ## copy data for download
        file.copy(paste(shinyAppFolder, "files/Metabolite_profile_showcase.txt", sep = ""), file)
      },
      contentType = "application/zip"
    )
    output$downloadMsMsData <- downloadHandler(
      filename = function() {
        return("MSMS_library_showcase.msp")
      },
      content = function(file) {
        ## copy data for download
        file.copy(paste(shinyAppFolder, "files/MSMS_library_showcase.msp", sep = ""), file)
      },
      contentType = "application/zip"
    )
    output$downloadFragmentMatrix <- downloadHandler(
      filename = function() {
        return("Fragment_matrix_showcase.csv")
      },
      content = function(file) {
        ## copy data for download
        file.copy(paste(shinyAppFolder, "files/Fragment_matrix_showcase.csv", sep = ""), file)
      },
      contentType = "application/zip"
    )
    ### TODO
    #output$myImage <- renderImage({
    #  #list(src = "/home/htreutle/Downloads/MetSWATH/svg/svgTestScoresPlot01.svg",
    #  list(src = "/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg",
    #       contentType = 'image/svg',
    #       width = 800,
    #       height = 600
    #  )
    #}, deleteFile = FALSE)
    
    ## TODO
    #http://127.0.0.1:25805/library/utils/html/zip.html
    serialization <- function(){
      #######################################
      ## enlist
      paramsList <- list(
        ## global MS2 filter
        globalFilter_ms2_masses1          = input$globalFilter_ms2_masses1,
        globalFilter_ms2_masses2          = input$globalFilter_ms2_masses2,
        globalFilter_ms2_masses3          = input$globalFilter_ms2_masses3,
        globalFilter_ms2_ppm              = input$globalFilter_ms2_ppm,
        #input$applyGlobalMS2filters
        ## HCA
        hcaFilterGroupOne                 = input$hcaFilterGroupOne,
        hcaFilterGroupTwo                 = input$hcaFilterGroupTwo,
        hcaFilter_average                 = input$hcaFilter_average,
        hcaFilter_lfc                     = input$hcaFilter_lfc,
        hcaFilterIncludeIgnoredPrecursors = input$hcaFilterIncludeIgnoredPrecursors,
        #input$applyHcaFilters
        hcaDistanceFunction               = input$hcaDistanceFunction,
        hcaClusterMethod                  = input$hcaClusterMethod,
        #input$drawHCAplots
        ## PCA
        pcaGroups                         = input$pcaGroups,
        pcaFilter_average                 = input$pcaFilter_average,
        pcaFilter_lfc                     = input$pcaFilter_lfc,
        pcaFilterIncludeIgnoredPrecursors = input$pcaFilterIncludeIgnoredPrecursors,
        #input$applyPcaFilters
        pcaScaling                        = input$pcaScaling,
        pcaLogTransform                   = input$pcaLogTransform,
        pcaDimensionOne                   = input$pcaDimensionOne,
        pcaDimensionTwo                   = input$pcaDimensionTwo,
        #input$drawPCAplots
        ## plot properties
        showClusterLabels                 = input$showClusterLabels,
        showScoresLabels                  = input$showScoresLabels,
        showLoadingsLabels                = input$showLoadingsLabels,
        showLoadingsAbundance             = input$showLoadingsAbundance,
        #showLoadingsLabels                = "Show labels"    %in% input$pcaLoadingsProperties,
        #showLoadingsAbundance             = "Show abundance" %in% input$pcaLoadingsProperties,
        ## search
        searchMS1orMS2                    = input$searchMS1orMS2,
        searchMS1mass                     = input$searchMS1mass,
        searchMS1massPpm                  = input$searchMS1massPpm,
        search_ms2_masses1                = input$search_ms2_masses1,
        search_ms2_masses2                = input$search_ms2_masses2,
        search_ms2_masses3                = input$search_ms2_masses3,
        searchMS2massPpm                  = input$searchMS2massPpm,
        searchIncludeIgnoredPrecursors    = input$searchIncludeIgnoredPrecursors
      )
      
      #######################################
      ## serialize parameter list
      
      ## built matrix with param names and param values
      tempMatrix    <- matrix(data = c(names(unlist(paramsList)),paste("'", unlist(paramsList), "'", sep = "")), ncol = 2)
      ## quotation of strings
      #textEntries <- is.na(as.numeric(tempMatrix[, 2]))
      #tempMatrix[textEntries, 2] <- paste("'", tempMatrix[textEntries, 2], "'", sep = "")
      ## param name = param value
      paramStrings  <- apply(tempMatrix, MARGIN = 1, FUN = function(x) {paste(x, collapse = "=")})
      ## collapse
      serialization <- paste(paramStrings, collapse = ";")
      return(serialization)
    }
    deserialization <- function(serialization){
      #######################################
      ## deserialize parameters
      
      #paramStrings <- strsplit(x = strsplit(x = serialization, split = ";")[[1]], split = "=")
      #for(i in 1:length(paramStrings))
      #  paramStrings[[i]] <- paste(paramStrings[[i]][[1]], paste("'", paramStrings[[i]][[2]], "'", sep = ""), sep = "=")
      #paramString <- paste(paramStrings, collapse = ",")
      paramString <- paste(strsplit(x = serialization, split = ";")[[1]], collapse = ",")
      parseText <- paste("paramsList <- list(", paramString, ")", sep = "")
      
      #parseText <- paste("paramsList <- list(", paste(paramStrings, collapse = ","), ")", sep = "")
      eval(parse(text = parseText))
      
      #######################################
      ## update parameter fields
      
      #updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
      #updateRadioButtons(session = session, inputId = "hcaFilterGroupOne", choices = dataList$groups, selected = selectedOne)
      #updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
      #updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[init]"), selected = lalala)
      #updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = dataList$groups)
      
      ## global MS2 filter
      updateTextInput(         session = session, inputId = "globalFilter_ms2_masses1",          value = paramsList$globalFilter_ms2_masses1)
      updateTextInput(         session = session, inputId = "globalFilter_ms2_masses2",          value = paramsList$globalFilter_ms2_masses2)
      updateTextInput(         session = session, inputId = "globalFilter_ms2_masses3",          value = paramsList$globalFilter_ms2_masses3)
      updateTextInput(         session = session, inputId = "globalFilter_ms2_ppm",              value = paramsList$globalFilter_ms2_ppm)
      #input$applyGlobalMS2filters
      ## HCA
      updateRadioButtons(      session = session, inputId = "hcaFilterGroupOne",                 selected = paramsList$hcaFilterGroupOne)
      updateRadioButtons(      session = session, inputId = "hcaFilterGroupTwo",                 selected = paramsList$hcaFilterGroupTwo)
      updateTextInput(         session = session, inputId = "hcaFilter_average",                 value = paramsList$hcaFilter_average)
      updateTextInput(         session = session, inputId = "hcaFilter_lfc",                     value = paramsList$hcaFilter_lfc)
      updateCheckboxInput(     session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = as.logical(paramsList$hcaFilterIncludeIgnoredPrecursors))
      #input$applyHcaFilters
      updateSelectInput(       session = session, inputId = "hcaDistanceFunction",               selected = paramsList$hcaDistanceFunction)
      updateSelectInput(       session = session, inputId = "hcaClusterMethod",                  selected = paramsList$hcaClusterMethod)
      #input$drawHCAplots
      ## PCA
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",                         selected = paramsList$pcaGroups)
      updateTextInput(         session = session, inputId = "pcaFilter_average",                 value = paramsList$pcaFilter_average)
      updateTextInput(         session = session, inputId = "pcaFilter_lfc",                     value = paramsList$pcaFilter_lfc)
      updateSelectInput(       session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", selected = paramsList$pcaFilterIncludeIgnoredPrecursors)
      #input$applyPcaFilters
      updateSelectInput(       session = session, inputId = "pcaScaling",                        selected = paramsList$pcaScaling)
      updateCheckboxInput(     session = session, inputId = "pcaLogTransform",                   value = as.logical(paramsList$pcaLogTransform))
      updateSelectInput(       session = session, inputId = "pcaDimensionOne",                   selected = paramsList$pcaDimensionOne)
      updateSelectInput(       session = session, inputId = "pcaDimensionTwo",                   selected = paramsList$pcaDimensionTwo)
      #input$drawPCAplots
      ## plot properties
      updateCheckboxInput(     session = session, inputId = "showClusterLabels",                 value = as.logical(paramsList$showClusterLabels))
      updateCheckboxInput(     session = session, inputId = "showScoresLabels",                  value = as.logical(paramsList$showScoresLabels))
      #updateCheckboxGroupInput(session = session, inputId = "pcaLoadingsProperties",             selected = c(ifelse(as.logical(paramsList$showLoadingsLabels), "Show labels", NULL), ifelse(as.logical(paramsList$showLoadingsAbundance), "Show abundance", NULL)))
      updateCheckboxInput(     session = session, inputId = "showLoadingsLabels",                value = as.logical(paramsList$showLoadingsLabels))
      updateCheckboxInput(     session = session, inputId = "showLoadingsAbundance",             value = as.logical(paramsList$showLoadingsAbundance))
      ## search
      updateRadioButtons(      session = session, inputId = "searchMS1orMS2",                    selected = paramsList$searchMS1orMS2)
      updateTextInput(         session = session, inputId = "searchMS1mass",                     value = paramsList$searchMS1mass)
      updateTextInput(         session = session, inputId = "searchMS1massPpm",                  value = paramsList$searchMS1massPpm)
      updateTextInput(         session = session, inputId = "search_ms2_masses1",                value = paramsList$search_ms2_masses1)
      updateTextInput(         session = session, inputId = "search_ms2_masses2",                value = paramsList$search_ms2_masses2)
      updateTextInput(         session = session, inputId = "search_ms2_masses3",                value = paramsList$search_ms2_masses3)
      updateTextInput(         session = session, inputId = "searchMS2massPpm",                  value = paramsList$searchMS2massPpm)
      updateCheckboxInput(     session = session, inputId = "searchIncludeIgnoredPrecursors",    value = as.logical(paramsList$searchIncludeIgnoredPrecursors))
      #input$applySearch
      
      #######################################
      ## update GUI according to parameters
      
      ###################
      ## global MS2 filter
      filter_ms2_masses1  <- paramsList$globalFilter_ms2_masses1
      filter_ms2_masses2  <- paramsList$globalFilter_ms2_masses2
      filter_ms2_masses3  <- paramsList$globalFilter_ms2_masses3
      filter_ms2_ppm      <- paramsList$globalFilter_ms2_ppm
      
      groupSet        <- dataList$groups
      filter_average  <- NULL
      filter_lfc      <- NULL
      includeIgnoredPrecursors  <- TRUE
      filter_ms1_masses <- NULL
      filter_ms1_ppm <- NULL
      
      resultObj <- doPerformFiltering(
        groupSet, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, 
        includeIgnoredPrecursors
      )
      filterGlobal <<- resultObj$filter
      state$globalMS2filterValid <<- TRUE
      updateGlobalMS2filterInformation()
      
      ###################
      ## HCA filter
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      groupOne        <- paramsList$hcaFilterGroupOne
      groupTwo        <- paramsList$hcaFilterGroupTwo
      filter_average  <- paramsList$hcaFilter_average
      filter_lfc      <- paramsList$hcaFilter_lfc
      includeIgnoredPrecursors  <- paramsList$hcaFilterIncludeIgnoredPrecursors
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      groupSet        <- c(groupOne, groupTwo)
      
      resultObj <- doPerformFiltering(
        groupSet, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
      )
      filterHca <<- resultObj$filter
      updateHcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
      checkHcaFilterValidity(numberOfPrecursorsFiltered)
      
      ###################
      ## draw HCA
      # TODO
      
      ###################
      ## PCA filter
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      groupSet        <- paramsList$pcaGroups
      filter_average  <- paramsList$pcaFilter_average
      filter_lfc      <- paramsList$pcaFilter_lfc
      includeIgnoredPrecursors  <- paramsList$pcaFilterIncludeIgnoredPrecursors
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      
      if(length(groupSet) != 2)
        filter_lfc <- NULL
      
      resultObj <- doPerformFiltering(
        groupSet, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
      )
      filterPca <<- resultObj$filter
      updatePcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterPca$numberOfPrecursorsFiltered
      checkPcaFilterValidity(numberOfPrecursorsFiltered)
      
      ###################
      ## draw PCA
      # TODO
      
      ###################
      ## search
      searchMode <- paramsListsearchMS1orMS2
      if(searchMode == 'MS\u00B9 feature mass'){
        #################################################
        ## get inputs
        filter_ms1_masses <- paramsListsearchMS1mass
        filter_ms1_ppm  <- paramsListsearchMS1massPpm
        
        if(nchar(trimws(filter_ms1_masses)) == 0)
          return()
        
        filter_ms2_masses1  <- NULL
        filter_ms2_masses2  <- NULL
        filter_ms2_masses3  <- NULL
        filter_ms2_ppm      <- NULL
      }
      if(searchMode == 'Fragment mass'){
        #################################################
        ## get inputs
        filter_ms2_masses1 <- paramsListsearch_ms2_masses1
        filter_ms2_masses2 <- paramsListsearch_ms2_masses2
        filter_ms2_masses3 <- paramsListsearch_ms2_masses3
        filter_ms2_ppm     <- paramsListsearchMS2massPpm
        
        filter_ms1_masses <- NULL
        filter_ms1_ppm <- NULL
      }
      
      filter_lfc      <- NULL
      filter_average  <- NULL
      groupSet        <- dataList$groups
      includeIgnoredPrecursors  <- paramsListsearchIncludeIgnoredPrecursors
      
      #################################################
      ## do filtering
      resultObj <- doPerformFiltering(
        groupSet, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
      )
      processSearchFilterResult(resultObj)
    }
  }## function(input, output, session)
)## shinyServer
