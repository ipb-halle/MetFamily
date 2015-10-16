
#########################################################################################
#########################################################################################
## libraries
#devtools::install_github("rstudio/htmltools")
library(htmltools)
#install.packages("shiny")
library(shiny)
#install.packages("shinyjs")
library(shinyjs)

#library(devtools)
#install_github("shinyTable", "trestletech")
#library(shinyTable)
#install.packages("DT")
library(DT)
library("Matrix")

#source("/home/htreutle/Code/Java/MetSWATH/ClusteringMS2SpectraGUI.R")
source("ClusteringMS2SpectraGUI.R")
#source("/vol/R/shiny/srv/shiny-server/MetFam/ClusteringMS2SpectraGUI.R")

#########################################################################################
#########################################################################################
## global variables
shinyServer(
  func = function(input, output, session) {
    #########################################################################################
    #########################################################################################
    ## global variables per user
    
    ## constants
    artifact <- "Ignore"
    artifactColor <- "red"
    minimumNumberOfPrecursorsForHca <- 6
    maximumNumberOfPrecursorsForHca <- 1000
    ## data
    dataList <- NULL
    currentDistanceMeasure <- NULL
    filterGlobal <- NULL
    filterHca <- NULL
    filterPca <- NULL
    filterSearch <- NULL
    clusterDataList <- NULL
    pcaDataList <- NULL
    ## program state
    state <- reactiveValues(
      globalMS2filterValid = FALSE, 
      hcafilterValid = FALSE, 
      pcafilterValid = FALSE, 
      searchfilterValid = FALSE, 
      runRightColumnWidth = 8, 
      legendColumnWidth = 2,
      showSideBar = TRUE, 
      showHCAplotPanel = FALSE, 
      showPCAplotPanel = FALSE, 
      plotToShow = "Display HCA", 
      plotHcaShown = FALSE,
      plotPcaShown = FALSE,
      precursorSetSelected = FALSE,
      analysisType = NULL,
      showLoadingsLabels = FALSE,
      showLoadingsAbundance = FALSE
    )
    ## selection
    selectedFragmentIndex <- NULL
    treeNodeFromTreeSelection <- NULL
    treeNodeSetFromMS2Fragment <- NULL
    pcaLoadingSetFromMS2Fragment <- NULL
    pcaLoadingFromPcaLoadingSelection <- NULL
    pcaLoadingSetColors <- NULL
    ## buttons
    applyGlobalMS2filtersButtonValue <- 0
    applySearchMS1ButtonValue <- 0
    applySearchMS2ButtonValue <- 0
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
    ## MS2 peaks
    fragmentsX <- NULL
    fragmentsY <- NULL
    fragmentsColor <- NULL
    fragmentsXhovered <- NULL
    fragmentsYhovered <- NULL
    fragmentshoveredColor <- NULL
    ## plot ranges
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
    ## table data
    tableFeaturesIntersection <- NULL
    tablePrecursorLabels <- NULL
    tableMS1abundanceDataFrame <- NULL
    tableMS2fragmentDataFrame <- NULL
    tableAnnotationDataFrame <- NULL
    table <- reactiveValues(
      dataFrame = NULL
    )
    tableCheckboxIdCounter <- 0
    selectedTreeNodeLabel <- NULL
    selectedPrecursorSet <- NULL
    
    #########################################################################################
    #########################################################################################
    ## functions
    
    ## POI selection
    getSelectedPOI_X <- function(mouseX, poiCoordinatesX, plotWidth, plotRangeX){
      factorX <- plotWidth  / plotRangeX
      
      mouseX <- mouseX * factorX
      poiCoordinatesX <- poiCoordinatesX * factorX
      
      distances <- abs(poiCoordinatesX - mouseX)
      distanceThreshold <- factorX * plotRangeX / 35
      
      minimumIndex <- which.min(distances)
      minimumDistance <- distances[[minimumIndex]]
      
      print(paste("--##--", minimumDistance, "vs", distanceThreshold))
      
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    getSelectedPOI_XY <- function(mouseX, mouseY, poiCoordinatesX, poiCoordinatesY, plotWidth, plotHeight, plotRangeX, plotRangeY){
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
      minimumDistance <- distances[[minimumIndex]]
      
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    ## Parse the input file
    getClusterData <- function(file){
      print(paste("getClusterData", is.null(file)))
      
      #########################################################################################
      ## reset
      
      ## reset plots
      doClearPlots()
      
      ## reset variables
      clusterDataList <<- NULL
      pcaDataList <<- NULL
      
      pcaLoadingSetFromMS2Fragment <<- NULL
      pcaLoadingFromPcaLoadingSelection <<- NULL
      
      selectedFragmentIndex <<- NULL
      treeNodeFromTreeSelection <<- NULL
      treeNodeSetFromMS2Fragment <<- NULL
      
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      fragmentsColor <<- NULL
      fragmentsXhovered <<- NULL
      fragmentsYhovered <<- NULL
      fragmentshoveredColor <<- NULL
      ## precursor selections
      table$dataFrame <<- NULL
      tableFeaturesIntersection <<- NULL
      tablePrecursorLabels <<- NULL
      tableMS1abundanceDataFrame <<- NULL
      tableMS2fragmentDataFrame <<- NULL
      selectedTreeNodeLabel <<- NULL
      selectedPrecursorSet <<- NULL
      
      ## reset state
      state$analysisType <<- NULL
      state$showHCAplotPanel <<- FALSE
      state$showPCAplotPanel <<- FALSE
      state$precursorSetSelected <<- FALSE
      
      ## reset plot range
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
      ## read data
      withProgress(message = 'Reading file...', value = 0, {
        #dataList <<- calcClusterData(file = "/mnt/VOL1/ABT/Alle/Balcke/MetSWATH/data/MS-DIAL/UC Davis/Results/201558139_matrixPrecursorsVersusFragmentsDeisotoped_withoutZerosTest01.txt")
        dataList <<- readClusterData(file = file, progress = TRUE)
      })
      print(paste("getClusterData do data finished", dataList$minimumMass))
      
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
      
      ## input fields
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_ppm", value = "20")
      updateTextInput(session = session, inputId = "hcaFilter_average", value = "0")
      updateTextInput(session = session, inputId = "hcaFilter_lfc", value = "0")
      updateTextInput(session = session, inputId = "pcaFilter_average", value = "0")
      updateTextInput(session = session, inputId = "pcaFilter_lfc", value = "0")
      updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
      updateCheckboxInput(session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", value = FALSE)
      #updateColourInput(session = session, inputId = "newAnnotationColor", allowedCols = colorPalette())
      
      #########################################################################################
      ## update filter
      ## TODO search filter
      filter <- doPerformFiltering(dataList$groups, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      
      filterGlobal <<- filter
      filterHca    <<- filter
      filterPca    <<- filter
      
      updateGlobalMS2filterInformation()
      updateHcaFilterInformation()
      updatePcaFilterInformation()
      
      state$globalMS2filterValid <<- TRUE
      state$hcaFilterValid <<- TRUE
      state$pcaFilterValid <<- TRUE
      
      checkHcaFilterValidity(filter$numberOfPrecursorsFiltered)
      checkPcaFilterValidity(filter$numberOfPrecursorsFiltered)
      
      state$plotHcaShown <<- FALSE
      state$plotPcaShown <<- FALSE
      
      output$information <- renderText({
        print(paste("init output$information", sep = ""))
        paste("Number of filtered precursors: ", dataList$numberOfPrecursorsFiltered, sep = "")
      })
      
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
          paste("Number of filtered precursors: ", filterGlobal$numberOfPrecursorsFiltered, sep = "")
        })
      }
    }
    updateHcaFilterInformation <- function(){
      if(is.null(filterHca)){
        ## errors
        output$hcaFilteredPrecursors <- renderText({
          print(paste("update output$hcaFilteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
      } else {
        ## no errors
        numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
        
        if(numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
          ## filter valid
          output$hcaFilteredPrecursors <- renderText({
            print(paste("update output$hcaFilteredPrecursors ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
            paste("Number of filtered precursors: ", numberOfPrecursorsFiltered, sep = "")
          })
        } else {
          ## filter invalid
          
          ## update info
          if(numberOfPrecursorsFiltered == 0){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors # = 0", sep = ""))
              paste("There are no precursors which fulfill the given criteria.", sep = "")
            })
          }
          if(numberOfPrecursorsFiltered > 0 & numberOfPrecursorsFiltered < minimumNumberOfPrecursorsForHca){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors 0 < # < ", minimumNumberOfPrecursorsForHca, sep = ""))
              paste("There are only ", numberOfPrecursorsFiltered, " precursors which fulfill the given criteria. There must be at least more than five precursors to preceed.", sep = "")
            })
          }
          if(numberOfPrecursorsFiltered > maximumNumberOfPrecursorsForHca){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors # > ", maximumNumberOfPrecursorsForHca, sep = ""))
              paste("There are ", numberOfPrecursorsFiltered, " precursors which fulfill the given criteria. There must be at most ", maximumNumberOfPrecursorsForHca, " precursors to preceed.", sep = "")
            })
          }
        }
      }
    }
    updatePcaFilterInformation <- function(){
      if(is.null(filterPca)){
        output$pcaFilteredPrecursors <- renderText({
          print(paste("update output$pcaFilteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
      } else {
        if(filterPca$numberOfPrecursorsFiltered > 0){
          output$pcaFilteredPrecursors <- renderText({
            print(paste("update output$pcaFilteredPrecursors", sep = ""))
            paste("Number of filtered precursors: ", filterPca$numberOfPrecursorsFiltered, sep = "")
          })
        } else {
          output$pcaFilteredPrecursors <- renderText({
            print(paste("update output$pcaFilteredPrecursors", sep = ""))
            paste("There are no precursors which fulfill the given criteria", sep = "")
          })
        }
      }
    }
    updateSearchInformation <- function(){
      if(is.null(filterSearch)){
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo invalid search", sep = ""))
          paste("There are invalid search values", sep = "")
        })
      } else {
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo", sep = ""))
          paste("Number of precursors: ", filterSearch$numberOfPrecursorsFiltered, sep = "")
        })
      }
    }
    ## perform filtering
    doPerformFiltering <- function(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors, preFilter = NULL){
      print(paste("Observe applyFilters1", "gs", paste(groupSet, collapse = "-"), "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "ig", includeIgnoredPrecursors))
      
      groupSetOriginal                 <- groupSet
      filter_averageOriginal           <- filter_average
      filter_lfcOriginal               <- filter_lfc
      filter_ms2_masses1Original       <- filter_ms2_masses1
      filter_ms2_masses2Original       <- filter_ms2_masses2
      filter_ms2_masses3Original       <- filter_ms2_masses3
      filter_ms2_ppmOriginal           <- filter_ms2_ppm
      filter_ms1_massOriginal          <- filter_ms1_mass
      filter_ms1_ppmOriginal           <- filter_ms1_ppm
      includeIgnoredPrecursorsOriginal <- includeIgnoredPrecursors
      
      #################################################
      ## parse inputs
      if(any(is.null(filter_average), length(filter_average) == 0, nchar(filter_average) == 0))
        filter_average <- NULL
      else
        filter_average <- as.numeric(filter_average)
      
      if(any(is.null(filter_lfc), length(filter_lfc) == 0, nchar(filter_lfc) == 0))
        filter_lfc <- NULL
      else{
        if(nchar(filter_lfc) == 0)
          filter_lfc <- NULL
        else
          filter_lfc <- as.numeric(filter_lfc)
      }
      
      if(any(is.null(filter_ms2_masses1), length(filter_ms2_masses1) == 0, nchar(filter_ms2_masses1) == 0))
        filter_ms2_masses1 <- NULL
      if(any(is.null(filter_ms2_masses2), length(filter_ms2_masses2) == 0, nchar(filter_ms2_masses2) == 0))
        filter_ms2_masses2 <- NULL
      if(any(is.null(filter_ms2_masses3), length(filter_ms2_masses3) == 0, nchar(filter_ms2_masses3) == 0))
        filter_ms2_masses3 <- NULL
      
      if(any(is.null(filter_ms2_ppm), length(filter_ms2_ppm) == 0, nchar(filter_ms2_ppm) == 0))
        filter_ms2_ppm <- NULL
      else
        filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
      
      if(any(is.null(filter_ms1_mass), length(filter_ms1_mass) == 0, nchar(filter_ms1_mass) == 0))
        filter_ms1_mass <- NULL
      else
        filter_ms1_mass <- as.numeric(filter_ms1_mass)
      
      if(any(is.null(filter_ms1_ppm), length(filter_ms1_ppm) == 0, nchar(filter_ms1_ppm) == 0))
        filter_ms1_ppm <- NULL
      else
        filter_ms1_ppm <- as.numeric(filter_ms1_ppm)
      
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
      if(!is.null(filter_average))
        error <- error | is.na(filter_average)
      if(!is.null(filter_lfc))
        error <- error | is.na(filter_lfc)
      if(!is.null(filter_ms2_masses1)){
        ms2Masses <- strsplit(x = filter_ms2_masses1, split = "[,; ]+")[[1]]
        filter_ms2_masses1 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses1[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses1))
      }
      if(!is.null(filter_ms2_masses2)){
        ms2Masses <- strsplit(x = filter_ms2_masses2, split = "[,; ]+")[[1]]
        filter_ms2_masses2 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses2[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses2))
      }
      if(!is.null(filter_ms2_masses3)){
        ms2Masses <- strsplit(x = filter_ms2_masses3, split = "[,; ]+")[[1]]
        filter_ms2_masses3 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses3[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses3))
      }
      if(!is.null(filter_ms2_ppm))
        error <- error | is.na(filter_ms2_ppm)
      if(!is.null(filter_ms1_mass))
        error <- error | is.na(filter_ms1_mass)
      if(!is.null(filter_ms1_ppm))
        error <- error | is.na(filter_ms1_ppm)
      
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
          filter_ms1_mass = filter_ms1_mass, filter_ms1_ppm = filter_ms1_ppm,
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
        filterHere$filter_ms1_massOriginal          <- filter_ms1_massOriginal         
        filterHere$filter_ms1_ppmOriginal           <- filter_ms1_ppmOriginal          
        filterHere$includeIgnoredPrecursorsOriginal <- includeIgnoredPrecursorsOriginal
        
        resultObj$error  <- error
        resultObj$filter <- filterHere
        print(paste("Observe applyFilters5", "n", resultObj$filter$numberOfPrecursorsFiltered))
      }
      
      return(resultObj)
    }
    ## Plots
    doClearPlots <- function(){
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("reset output$clusterPlotDendrogram"))
        NULL
      })
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("reset output$clusterPlotHeatmap"))
        NULL
      })
      output$clusterPlotHeatmapLegend <- renderPlot({
        print(paste("reset output$clusterPlotHeatmapLegend"))
        NULL
      })
      output$clusterPlotMS2 <- renderPlot({
        print(paste("reset output$clusterPlotMS2"))
        NULL
      })
      output$pcaPlotScores <- renderPlot({
        print(paste("reset output$pcaPlotScores"))
        NULL
      })
      output$pcaPlotLoadings <- renderPlot({
        print(paste("reset output$pcaPlotLoadings"))
        NULL
       })
      output$calcClusterPlotAnnoLegend <- renderPlot({
        print(paste("reset output$calcClusterPlotAnnoLegend"))
        NULL
      })
      output$calcMS2PlotAnnoLegend <- renderPlot({
        print(paste("reset output$calcMS2PlotAnnoLegend"))
        NULL
      })
    }
    drawPcaPlots <- function(consoleInfo = NULL){
      drawPcaScoresPlot(consoleInfo = consoleInfo)
      drawPcaLoadingsPlot(consoleInfo = consoleInfo)
    }
    drawPcaScoresPlot <- function(consoleInfo = NULL){
      output$pcaPlotScores <- renderPlot({
        print(paste("### psc ###", consoleInfo))
        plotPCAscores(pcaObj = pcaDataList$pcaObj, dataList = dataList, filterObj = filterPca, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, xInterval = pcaScoresPlotRange$xInterval, yInterval = pcaScoresPlotRange$yInterval)
      })
    }
    drawPcaLoadingsPlot <- function(consoleInfo = NULL){
      output$pcaPlotLoadings <- renderPlot({
        print(paste("### psc ###", consoleInfo))
        plotPCAloadings(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filterPca$filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, pcaLoading = pcaLoadingFromPcaLoadingSelection, pcaLoadingSet = pcaLoadingSetFromMS2Fragment, showLoadingsLabels = state$showLoadingsLabels, showLoadingsAbundance = state$showLoadingsAbundance, xInterval = pcaLoadingsPlotRange$xInterval, yInterval = pcaLoadingsPlotRange$yInterval)
      })
    }
    drawHeatmapLegend <- function(consoleInfo = NULL){
      output$clusterPlotHeatmapLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcClusterPlotHeatmapLegend(dataList = dataList)
      })
    }
    drawAnnotationLegend <- function(consoleInfo = NULL){
      output$calcClusterPlotAnnoLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcClusterPlotAnnoLegend(dataList = dataList)
      })
    }
    drawMS2Legend <- function(consoleInfo = NULL){
      output$calcMS2PlotAnnoLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcMS2PlotLegend(dataList = dataList)
      })
    }
    drawFragmentPlot <- function(consoleInfo = NULL){
      output$fragmentPlot <- renderPlot({
        print(paste("### MS2 ### all", consoleInfo))
        plotFragments(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
      })
    }
    ## resets
    resetHcaStuff <- function(){
      treeNodeFromTreeSelection <<- NULL
      treeNodeSetFromMS2Fragment <<- NULL
      
      dendrogramPlotRange$xMin <<- 1
      dendrogramPlotRange$xMax <<- filterHca$numberOfPrecursorsFiltered
      dendrogramPlotRange$xInterval <<- c(1, filterHca$numberOfPrecursorsFiltered)
      dendrogramPlotRange$xIntervalSize <<- filterHca$numberOfPrecursorsFiltered - 1
      
      state$precursorSetSelected <<- FALSE
      tableFeaturesIntersection <<- NULL
      tablePrecursorLabels <<- NULL
      tableMS1abundanceDataFrame <<- NULL
      tableMS2fragmentDataFrame <<- NULL
      tableAnnotationDataFrame <<- NULL
      selectedTreeNodeLabel <<- NULL
      selectedPrecursorSet <<- NULL
      table$dataFrame <<- NULL
    }
    resetPcaStuff <- function(){
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
      
      pcaLoadingFromPcaLoadingSelection <<- NULL
    }
    resetMS2PlotRange <- function(){
      ms2PlotRange$xMin <<- dataList$minimumMass
      ms2PlotRange$xMax <<- dataList$maximumMass
      ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
    }
    ## draw plots
    drawMS2Plot <- function(consoleInfo = NULL){
      output$clusterPlotMS2 <- renderPlot({
        print(paste("### MS2 ###", consoleInfo))
        calcClusterPlotMS2(
          dataList = dataList, 
          fragmentsX = fragmentsX, 
          fragmentsY = fragmentsY, 
          fragmentsColor = fragmentsColor, 
          fragmentsX_02 = fragmentsXhovered, 
          fragmentsY_02 = fragmentsYhovered, 
          fragmentsColor_02 = fragmentshoveredColor, 
          xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex  
        )
      })
    }
    drawDendrogramPlot <- function(consoleInfo = NULL, withHeatmap = FALSE){
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("### den ###", consoleInfo))
        calcClusterPlotDendrogram(
          dataList = dataList, 
          filter = filterHca$filter, 
          clusterDataList = clusterDataList, 
          annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, 
          annoPresentColorsList = dataList$annoPresentColorsList, 
          distanceMeasure = currentDistanceMeasure, 
          nodeIndex = treeNodeFromTreeSelection, 
          nodeIndeces = treeNodeSetFromMS2Fragment, 
          xInterval = dendrogramPlotRange$xInterval
        )
      })
      if(!withHeatmap)
        return()
      
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("### hea ### update range output$clusterPlotHeatmap"))
        calcClusterPlotHeatmap(dataList = dataList, filterObj = filterHca, clusterDataList = clusterDataList, xInterval = dendrogramPlotRange$xInterval)
      })
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
      updatePlotsWithAnnotations()
    }
    setArtifactState <- function(precursorSet, isArtifact){
      ## add
      dataList$annoArrayIsArtifact[precursorSet] <<- isArtifact
      
      for(i in 1:length(precursorSet))
        print(paste(dataList$precursorLabels[[precursorSet[[i]]]], isArtifact[[i]]))
      
      ## update gui
      #print("setArtifact state updateAnnoGui")
      updateAnnoGui(precursorSet)
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
      updatePlotsWithAnnotations()
    }
    commonAnnotations <- function(precursorSet){
      if(all(dataList$annoArrayIsArtifact[precursorSet]))
        return(artifact)
      
      ## at least one non-artifact precursor present
      intersection <- unlist(dataList$annoPresentAnnotationsList)
      for(precursor in precursorSet)
        if(!dataList$annoArrayIsArtifact[[precursor]])
          intersection <- intersect(x = intersection, y = unlist(dataList$annoArrayOfLists[[precursor]]))
      return(intersection)
    }
    updateAnnoGui <- function(precursorSet){
      ## annotation
      commonAnnos <- commonAnnotations(precursorSet)
      
      if(length(commonAnnos) > 0)
        updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = commonAnnos)
      else
        updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[none]"))
      updateSelectInput(session = session, inputId = "previousAnnotationValue", choices = dataList$annoPresentAnnotationsList)
      
      if(length(commonAnnos) > 0){
        shinyjs::enable("removePresentAnnotation")
        shinyjs::enable("setPresentAnnotationPrimary")
      } else {
        shinyjs::disable("removePresentAnnotation")
        shinyjs::disable("setPresentAnnotationPrimary")
      }
      
      ## table
      if(!is.null(precursorSet)){
        resultObj2 <- getTableFromTreeSelection2(dataList = dataList, clusterDataList = clusterDataList, precursorSet = precursorSet, featuresIntersection = tableFeaturesIntersection)
        tableAnnotationDataFrame <<- resultObj2$annotationDataFrame
        table$dataFrame <<- createTable(precursorSet)
      } else
        table$dataFrame <<- NULL
    }
    updatePlotsWithAnnotations <- function(){
      ## plots
      if(state$showPCAplotPanel)
        drawPcaLoadingsPlot(consoleInfo = "updateAnnoGui output$pcaPlotLoadings")
      if(state$showHCAplotPanel)
        drawDendrogramPlot(consoleInfo = "updateAnnoGui output$clusterPlotDendrogram")
      
      drawAnnotationLegend(consoleInfo = "init output$calcClusterPlotAnnoLegend")
    }
    ## table functions
    createInputFields <- function(FUN, id, values) {
      ## running id
      tableCheckboxIdCounter <<- tableCheckboxIdCounter + 1
      id <- paste(id, "_", tableCheckboxIdCounter, sep = "")
      
      ## create a character vector of shiny inputs
      inputs <- character(length(values))
      for (i in 1:length(values)){
        itemId    <- paste(id, "_", i, sep = "")
        inputs[i] <- as.character(FUN(inputId = itemId, label = NULL, value = values[[i]]))
      }
      return(inputs)
    }
    getInputValues <- function(id, len) {
      id <- paste(id, "_", tableCheckboxIdCounter, sep = "")
      
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
    createTable <- function(precursorSet){
      if(is.null(precursorSet))
        return(NULL)
      
      #tablePrecursorLabels
      isArtifact <- dataList$annoArrayIsArtifact[precursorSet]
      dataFrameIgnore <- data.frame(Ignore = createInputFields(FUN = checkboxInput, id = 'Ignore', values = isArtifact))
      
      dataFrame <<- cbind(
        tableMS1abundanceDataFrame,
        dataFrameIgnore,
        tableAnnotationDataFrame,
        tableMS2fragmentDataFrame
      )
      
      return(dataFrame)
    }
    setTable <- function(){
      output$table <- DT::renderDataTable(
        expr = table$dataFrame,
        server = FALSE, escape = FALSE, selection = "none",
        options = list(
          preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
          drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }')
        )
      )
    }
    ## selection stuff
    selectionByMs2Reset <- function(){
      
    }
    selectionByMs2 <- function(){
      
    }
    selectionByHcaReset <- function(){
      
    }
    selectionByHca <- function(){
      
    }
    selectionByPcaReset <- function(){
      
    }
    selectionByPca <- function(){
      
    }
    selectionBySearchReset <- function(){
      
    }
    selectionBySearch <- function(){
      
    }
    
    #########################################################################################
    #########################################################################################
    ## observer
    
    ## side panel
    obsTabs <- observeEvent(input$runTabs, {
      tabId <- input$runTabs
      print(paste("Observe tab selection", tabId))
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
    })
    obsFile <- observeEvent(input$matrixFile$datapath, {
      file     <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      print(paste("Observe file for data", fileName))
      if(!is.null(file)){
        getClusterData(file)
        output$fileInfo <- renderText({fileName})
      }
    })
    obsShowSideBar <- observeEvent(input$showSideBar, {
      showSideBar <- input$showSideBar
      print(paste("Observe showSideBar", showSideBar))
      state$showSideBar <<- showSideBar
      if(showSideBar){
        state$runRightColumnWidth <<- 8
        state$legendColumnWidth <<- 2
      }
      else{
        state$runRightColumnWidth <<- 12
        state$legendColumnWidth <<- 1
      }
    })
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
    ## listen to submit filter button events
    obsApplySearchMS1 <- observeEvent(input$applySearchMS1, {
      applySearchMS1 <- as.numeric(input$applySearchMS1)
      
      print(paste("Observe applySearchMS1", applySearchMS1))
      
      #################################################
      ## check if button was hit
      if(applySearchMS1 == applySearchMS1ButtonValue)
        return()
      applySearchMS1ButtonValue <<- applySearchMS1
      
      #################################################
      ## get inputs
      filter_ms1_mass <- input$searchMS1mass
      filter_ms1_ppm  <- input$searchMS1massPpm
      includeIgnoredPrecursors  <- input$searchMS1includeIgnoredPrecursors
      
      
      filter_ms2_masses1  <- NULL
      filter_ms2_masses2  <- NULL
      filter_ms2_masses3  <- NULL
      filter_ms2_ppm      <- NULL
      
      groupSet        <- dataList$groups
      filter_average  <- NULL
      filter_lfc      <- NULL
      
      #################################################
      ## do filtering
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors)
      processSearchFilterResult(resultObj)
    })
    obsApplySearchMS2 <- observeEvent(input$applySearchMS2, {
      applySearchMS2 <- as.numeric(input$applySearchMS2)
      
      print(paste("Observe applySearchMS2", applySearchMS2))
      
      #################################################
      ## check if button was hit
      if(applySearchMS2 == applySearchMS2ButtonValue)
        return()
      applySearchMS2ButtonValue <<- applySearchMS2
      
      #################################################
      ## get inputs
      filter_ms2_masses1 <- input$search_ms2_masses1
      filter_ms2_masses2 <- input$search_ms2_masses2
      filter_ms2_masses3 <- input$search_ms2_masses3
      filter_ms2_ppm      <- input$searchMS2massPpm
      includeIgnoredPrecursors  <- input$searchMS2includeIgnoredPrecursors
      
      groupSet        <- NULL
      filter_average  <- NULL
      filter_lfc      <- NULL
      filter_ms1_mass <- NULL
      filter_ms1_ppm <- NULL
      
      #################################################
      ## do filtering
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors)
      processSearchFilterResult(resultObj)
    })
    processSearchFilterResult <- function(resultObj){
      #################################################
      ## info / error
      if(resultObj$error){
        filterSearch <<- NULL
        state$searchfilterValid <<- FALSE
        return()
      } else {
        filterSearch <<- resultObj$filter
        state$searchfilterValid <<- TRUE
      }
      
      updateSearchInformation()
      
      #################################################
      ## update plots
      
      ## TODO
      if(state$showHCAplotPanel)
        drawDendrogramPlot(consoleInfo = "update search output$clusterPlotDendrogram")
      
      if(state$showPCAplotPanel)
        ## update PCA plots
        drawPcaPlots(consoleInfo = "search output$pcaPlotScores")
    }
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
      filter_ms1_mass <- NULL
      filter_ms1_ppm <- NULL
      
      #################################################
      ## do filtering
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        filterGlobal <<- NULL
        state$globalMS2filterValid <<- FALSE
        return()
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
      filter_ms1_mass <- NULL
      filter_ms1_ppm  <- NULL
      
      groupSet        <- c(groupOne, groupTwo)
      
      #################################################
      ## do filtering and update
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors)
      
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
      filter_ms1_mass <- NULL
      filter_ms1_ppm  <- NULL
      
      if(length(groupSet) != 2)
        filter_lfc <- NULL
      
      #################################################
      ## do filtering and update
      resultObj <- doPerformFiltering(groupSet, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors)
      
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
    observeGroupSet <- observeEvent(input$pcaGroups, {
      print(paste("observe groups change", paste(input$pcaGroups, collapse = "-"), length(input$pcaGroups), length(input$pcaGroups) == 2))
      shinyjs::toggleState("pcaFilter_lfc", length(input$pcaGroups) == 2)
      #if(length(input$pcaGroups) == 0)
      #  shinyjs::disable("applyPcaFilters")
      #  #disableActionButton(session, "applyPcaFilters")
      #else
      #  shinyjs::enable("applyPcaFilters")
      #  #enableActionButton(session, "applyPcaFilters")
    })
    ## listen to draw button events
    obsDrawHCA <- observeEvent(input$drawHCAplots, {
      drawPlots <- as.numeric(input$drawHCAplots)
      
      print(paste("Observe draw HCA plots", drawPlots))
      
      #################################################
      ## check if button was hit
      if(drawPlots == drawHCAButtonValue)
        return()
      drawHCAButtonValue <<- drawPlots
      
      #################################################
      ## update plots
      
      distance <- input$hcaDistanceFunction
      clusterMethod <- input$hcaClusterMethod
      print(paste("Observe draw HCA plots", "D", distance, "M", clusterMethod))
      currentDistanceMeasure <<- distance
      
      resetHcaStuff()
      
      ##########################
      ## draw
      
      ## compute distance matrix
      withProgress(message = 'Calculating distances...', value = 0, {
        distanceMatrix <- calculateDistanceMatrix(dataList = dataList, filter = filterHca$filter, distance = currentDistanceMeasure, progress = TRUE)
      })
      ## compute cluster
      withProgress(message = 'Calculating cluster...', value = 0, {
        clusterDataList <<- calculateCluster(progress = TRUE, dataList = dataList, filter = filterHca$filter, distanceMatrix = distanceMatrix, method = clusterMethod)
      })
      
      ## clear info and tip
      output$information <- renderText({
        print(paste("update output$information clear"))
        paste("", sep = "")
      })
      output$tip <- renderText({
        print(paste("update output$tip clear"))
        paste("", sep = "")
      })
      
      ## tree seelction from fragment
      if(!is.null(selectedFragmentIndex))
        ## tree highlighting
        treeNodeSetFromMS2Fragment <<- getSetOfSubTreesFromRootForMass(dataList = dataList, fragmentMass = fragmentsX[[selectedFragmentIndex]], filter = filterHca$filter, clusterDataList = clusterDataList)
      
      drawDendrogramPlot(consoleInfo = "init output$clusterPlotDendrogram", withHeatmap = TRUE)
      drawHeatmapLegend(consoleInfo = "init output$clusterPlotHeatmapLegend")
      drawAnnotationLegend(consoleInfo = "init output$calcClusterPlotAnnoLegend")
      drawMS2Plot(consoleInfo = "init output$clusterPlotMS2")
      drawMS2Legend(consoleInfo = "init output$ms2LegendPlot")
      
      state$showHCAplotPanel <<- TRUE
    })
    obsDrawPCA <- observeEvent(input$drawPCAplots, {
      drawPlots <- as.numeric(input$drawPCAplots)
      
      print(paste("Observe draw PCA plots", drawPlots))
      
      #################################################
      ## check if button was hit
      if(drawPlots == drawPCAButtonValue)
        return()
      drawPCAButtonValue <<- drawPlots
      
      pcaScaling      <- input$pcaScaling
      pcaLogTransform <- input$pcaLogTransform
      pcaDimensionOne <<- as.numeric(input$pcaDimensionOne)
      pcaDimensionTwo <<- as.numeric(input$pcaDimensionTwo)
      
      ## calc PCA
      pca <- calcPCA(dataList = dataList, filterObj = filterPca, scaling = pcaScaling, logTransform = pcaLogTransform)
      
      pcaDataList <<- list()
      pcaDataList$pcaObj <<- pca
      pcaDataList$pcaScoresX <<- pca$scores[, pcaDimensionOne]
      pcaDataList$pcaScoresY <<- pca$scores[, pcaDimensionTwo]
      pcaDataList$pcaLoadingsX <<- pca$loadings[, pcaDimensionOne]
      pcaDataList$pcaLoadingsY <<- pca$loadings[, pcaDimensionTwo]
      pcaDataList$dimensionOne <<- pcaDimensionOne
      pcaDataList$dimensionTwo <<- pcaDimensionTwo
      
      resetPcaStuff()
      
      ## loadings highlighting from MS2
      if(!is.null(selectedFragmentIndex)){
        pcaLoadingSetFromMS2Fragment <<- getPrecursorsByFragment(dataList = dataList, filter = filterPca$filter, fragmentMass = fragmentsX[[selectedFragmentIndex]])
        
        #fragmentMass     <- fragmentsX[[selectedFragmentIndex]]
        #fragmentIndex    <- which(dataList$fragmentMasses == fragmentMass)
        #precursorMatches <- dataList$featureMatrix[filterPca$filter, fragmentIndex] != 0
        #pcaLoadingSetFromMS2Fragment <<- precursorMatches
      }
      ## loadings annotation
      #pcaLoadingSetColors <- lapply(X = dataList$)
      
      drawPcaPlots(consoleInfo = "drawPCA output$pcaPlotScores")
      drawMS2Plot(consoleInfo = "drawPCA output$clusterPlotMS2")
      
      state$showPCAplotPanel <<- TRUE
    })
    ## listen to dendrogram plot mouse events
    obsDendrogramHover <- observeEvent(input$clusterPlotDendrogram_hover, {
      hoverX <- input$clusterPlotDendrogram_hover$x
      hoverY <- input$clusterPlotDendrogram_hover$y
      
      plotWidth  <- session$clientData$output_clusterPlotDendrogram_width
      plotHeight <- session$clientData$output_clusterPlotDendrogram_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = 1
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
        
        if(all(!is.null(selectedTreeNodeLabel), minimumLabel == selectedTreeNodeLabel)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
        } else {
          #################################################
          ## fetch ms2 spectrum
          resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = minimumLabel)
          fragmentsXhovered <<- resultObj$fragmentMasses
          fragmentsYhovered <<- resultObj$fragmentAbundances
          fragmentshoveredColor <<- resultObj$fragmentColor
          
          #################################################
          ## output as message and plots
          output$information <- renderText({
            print(paste("update output$information", resultObj$infoText))
            paste(resultObj$infoText, sep = "; ")
          })
        }
      }
      
      ## MS2 plot
      drawMS2Plot(consoleInfo = "dendrogram hover output$clusterPlotMS2")
    })
    obsDendrogramClick <- observeEvent(input$clusterPlotDendrogram_click, {
      clickX <- input$clusterPlotDendrogram_click$x
      clickY <- input$clusterPlotDendrogram_click$y
      
      brush <- input$clusterPlotDendrogram_brush
      
      plotWidth  <- session$clientData$output_clusterPlotDendrogram_width
      plotHeight <- session$clientData$output_clusterPlotDendrogram_height
      
      if(is.null(clickX) | is.null(clickY))
        return()
      if(!is.null(brush))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = 1
      )
      print(paste("Observe dendrogram click", is.null(minimumIndex), minimumIndex))
      if(is.null(minimumIndex)){
        ## reset stuff
        pcaLoadingSetFromMS2Fragment <<- NULL
        pcaLoadingFromPcaLoadingSelection <<- NULL
        
        selectedFragmentIndex <<- NULL
        treeNodeFromTreeSelection <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentshoveredColor <<- NULL
        
        output$table <- DT::renderDataTable(NULL)
        ## precursor selections
        table$dataFrame <<- NULL
        tablePrecursorLabels <<- NULL
        tableMS1abundanceDataFrame <<- NULL
        tableMS2fragmentDataFrame <<- NULL
        tableAnnotationDataFrame <<- NULL
        selectedTreeNodeLabel <<- NULL
        selectedPrecursorSet <<- NULL
        state$precursorSetSelected <<- FALSE
        updateAnnoGui(NULL)
      } else {
        ## tree selection
        minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
        #minimumText <- clusterDataList$poiText[[minimumIndex]]
        
        #################################################
        ## fetch ms2 spectrum
        pcaLoadingSetFromMS2Fragment <<- NULL
        
        selectedFragmentIndex <<- NULL
        treeNodeFromTreeSelection <<- minimumLabel
        treeNodeSetFromMS2Fragment <<- NULL
        
        resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = minimumLabel)
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsColor <<- resultObj$fragmentColor
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentshoveredColor <<- NULL
        
        #################################################
        ## output as message
        output$information <- renderText({
          print(paste("update output$information", resultObj$numberOfPrecursors))
          paste("", resultObj$numberOfPrecursors, " precursors selected", sep = "")
        })
        
        ## MetFrag link
        if(!is.null(resultObj$landingPageUrl))
          output$metFragLink <- renderText({
            print(paste("update output$metFragLink I ", resultObj$landingPageUrl))
            paste("<a href=", gsub(pattern = " ", replacement = "%20", x = resultObj$landingPageUrl)[[1]], " target=\"_blank\">Send to MetFrag!</a>", sep = "")
          })
        else
          output$metFragLink <- renderText({
            print(paste("update output$metFragLink I empty"))
            paste("", sep = "")
          })
        
        ## table
        resultObj2 <- getTableFromTreeSelection(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = minimumLabel)
        setTable()
        
        #input$table_rows_selected
        state$precursorSetSelected <<- TRUE
        selectedTreeNodeLabel <<- minimumLabel
        selectedPrecursorSet <<- resultObj$precursorSet
        
        tableFeaturesIntersection  <<- resultObj2$featuresIntersection
        tablePrecursorLabels       <<- resultObj2$precursorLabels
        tableMS1abundanceDataFrame <<- resultObj2$ms1abundanceDataFrame
        tableMS2fragmentDataFrame  <<- resultObj2$ms2fragmentDataFrame
        tableAnnotationDataFrame   <<- resultObj2$annotationDataFrame
        
        updateAnnoGui(resultObj$precursorSet)
        
        ## pca selection
        if(minimumLabel < 0){
          ## leaf selected
          precursorIndex <- -minimumLabel
          precursorIndex <- filterHca$filter[[precursorIndex]]
          if(all(state$showPCAplotPanel, precursorIndex %in% filterPca$filter))
            pcaLoadingFromPcaLoadingSelection <<- match(x = precursorIndex, table = filterPca$filter)
          else
            pcaLoadingFromPcaLoadingSelection <<- NULL
        } else
          pcaLoadingFromPcaLoadingSelection <<- NULL
      }
      
      #################################################
      ## plots
      
      ## cluster dendrogram
      drawDendrogramPlot(consoleInfo = "dendrogram click output$clusterPlotDendrogram")
      
      ## MS2 plot
      resetMS2PlotRange()
      drawMS2Plot(consoleInfo = "dendrogram click output$clusterPlotMS2")
      
      if(state$showPCAplotPanel)
        ## update PCA plots
        drawPcaPlots(consoleInfo = "dendrogram click output$pcaPlotScores")
    })
    obsDendrogramdblClick <- observeEvent(input$clusterPlotDendrogram_dblclick, {
      brush <- input$clusterPlotDendrogram_brush
      
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
      
      drawDendrogramPlot(consoleInfo = "update range output$clusterPlotDendrogram", withHeatmap = TRUE)
    })
    ## listen to heatmap plot mouse events
    obsHeatmaphover <- observeEvent(input$clusterPlotHeatmap_hover, {
      hoverX <- input$clusterPlotHeatmap_hover$x
      hoverY <- input$clusterPlotHeatmap_hover$y
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(hoverX < 0.5 | hoverX > (clusterDataList$numberOfPrecursorsFiltered + 0.5))
        return()
      if(hoverY < 0 | hoverY > 3)
        return()
      
      print(paste("Observe heatmap hover", hoverX, hoverY))
      
      treeLeafIndex2 <- as.numeric(format(x = hoverX, digits = 0))
      treeLeafIndex  <- clusterDataList$cluster$order[[treeLeafIndex2]]
      precursorIndex <- filterHca$filter[[treeLeafIndex]]
      
      msg <- list()
      msg[[length(msg) + 1]] <- "Precursor: "
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
    obsMS2hover <- observeEvent(input$clusterPlotMS2_hover, {
      hoverX <- input$clusterPlotMS2_hover$x
      hoverY <- input$clusterPlotMS2_hover$y
      plotWidth  <- session$clientData$output_clusterPlotMS2_width
      plotHeight  <- session$clientData$output_clusterPlotMS2_height
      
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
        #  paste("")
        #})
      } else {
        ## point selected
        output$information <- renderText({
          print(paste("update output$information"))
          paste("Fragment: m/z = ", fragmentsX[[minimumIndex]], ", abundance = ", fragmentsY[[minimumIndex]], "", sep = "")
        })
      }
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Click a fragment peak to view information", "Brush and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    obsMS2click <- observeEvent(input$clusterPlotMS2_click, {
      clickX <- input$clusterPlotMS2_click$x
      clickY <- input$clusterPlotMS2_click$y
      
      brush  <- input$clusterPlotMS2_brush
      
      plotWidth  <- session$clientData$output_clusterPlotMS2_width
      plotHeight <- session$clientData$output_clusterPlotMS2_height
      
      if(is.null(clickX) | is.null(clickY))
        return()
      if(!is.null(brush))
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
        selectedFragmentIndex <<- NULL
        
        ## HCA
        if(!is.null(treeNodeSetFromMS2Fragment))
          treeNodeSetFromMS2Fragment <<- NULL
        ## PCA
        if(!is.null(pcaLoadingSetFromMS2Fragment))
          pcaLoadingSetFromMS2Fragment <<- NULL
      } else {
        ##########################################
        ## peak click
        selectedFragmentIndex <<- minimumIndex
        
        ## HCA
        ## fetch subroots of subtrees comprising the selected fragment
        if(state$showHCAplotPanel)
          ## tree highlighting
          treeNodeSetFromMS2Fragment <<- getSetOfSubTreesFromRootForMass(dataList = dataList, fragmentMass = fragmentsX[[minimumIndex]], filter = filterHca$filter, clusterDataList = clusterDataList)
        ## PCA
        if(state$showPCAplotPanel)
          pcaLoadingSetFromMS2Fragment <<- getPrecursorsByFragment(dataList = dataList, filter = filterPca$filter, fragmentMass = fragmentsX[[selectedFragmentIndex]])
      }
      
      ## HCA
      if(state$showHCAplotPanel)
        drawDendrogramPlot(consoleInfo = "MS2 click output$clusterPlotDendrogram")
      ## PCA
      if(state$showPCAplotPanel)
        drawPcaLoadingsPlot(consoleInfo = "MS2 click output$pcaPlotLoadings")
      
      ## update node selection
      drawMS2Plot(consoleInfo = "MS2 click output$clusterPlotMS2")
    })
    obsMS2dblClick <- observeEvent(input$clusterPlotMS2_dblclick, {
      brush <- input$clusterPlotMS2_brush
      
      if(any(is.null(fragmentsX), length(fragmentsX) == 0))
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
      
      drawMS2Plot(consoleInfo = "range update output$clusterPlotMS2")
    })
    ## listen to PCA events
    obsPCAscoresHover <- observeEvent(input$pcaPlotScores_hover, {
      hoverX <- input$pcaPlotScores_hover$x
      hoverY <- input$pcaPlotScores_hover$y
      plotWidth  <- session$clientData$output_pcaPlotScores_width
      plotHeight <- session$clientData$output_pcaPlotScores_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA scores hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaScoresX, poiCoordinatesY = pcaDataList$pcaScoresY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaScoresPlotRange$xIntervalSize, plotRangeY = pcaScoresPlotRange$yIntervalSize
      )
      if(is.null(minimumIndex))
        return()
      
      print(paste("Observe PCA scores hover", hoverX, hoverY, minimumIndex))
      
      dataColumnName <- dataList$dataColumnsNameFunctionFromNames(filterPca$groups)[[minimumIndex]]
      group <- dataList$groupNameFunctionFromDataColumnName(dataColumnName)
      output$information <- renderText({
        print(paste("update output$information PCA scores hover", sep = ""))
        paste(dataColumnName , " is a replicate of group ", group, sep = "")
      })
    })
    obsPCAscoresDblClick <- observeEvent(input$pcaPlotScores_dblclick, {
      brush <- input$pcaPlotScores_brush
      
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
      drawPcaScoresPlot(consoleInfo = "range update output$pcaPlotScores")
    })
    obsPCAloadingsClick <- observeEvent(input$pcaPlotLoadings_click, {
      clickX <- input$pcaPlotLoadings_click$x
      clickY <- input$pcaPlotLoadings_click$y
      plotWidth  <- session$clientData$output_pcaPlotLoadings_width
      plotHeight <- session$clientData$output_pcaPlotLoadings_height
      brush <- input$pcaPlotLoadings_brush
      
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
        pcaLoadingFromPcaLoadingSelection <<- NULL
        pcaLoadingSetFromMS2Fragment <<- NULL
        treeNodeFromTreeSelection <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentshoveredColor <<- NULL
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings click", sep = ""))
          paste("", sep = "")
        })
      } else {
        ## loadng selection
        precursorIndex <- filterPca$filter[[minimumIndex]]
        
        if(state$showHCAplotPanel & precursorIndex %in% filterHca$filter){
          #################################################
          ## fetch ms2 spectrum
          #minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
          #leafLabel <- -filterHca$filter[[precursorIndex]]
          leafLabel <- -match(x = precursorIndex, table = filterHca$filter)
          treeNodeFromTreeSelection <<- leafLabel
          
          ## MetFrag link
          resultObj3 <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = leafLabel)
          if(!is.null(resultObj3$landingPageUrl)){
            output$metFragLink <- renderText({
              print(paste("update output$metFragLink II", resultObj3$landingPageUrl))
              paste("<a href=", gsub(pattern = " ", replacement = "%20", x = resultObj3$landingPageUrl)[[1]], " target=\"_blank\">Send to MetFrag!</a>", sep = "")
            })
          } else {
            output$metFragLink <- renderText({
              print(paste("update output$metFragLink II empty"))
              paste("", sep = "")
            })
          }
          
          ## table
          resultObj2 <- getTableFromTreeSelection(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = leafLabel)
          setTable()
          
          #input$table_rows_selected
          state$precursorSetSelected <<- TRUE
          selectedTreeNodeLabel <<- leafLabel
          selectedPrecursorSet <<- resultObj$precursorSet
          
          tableFeaturesIntersection  <<- resultObj2$featuresIntersection
          tablePrecursorLabels       <<- resultObj2$precursorLabels
          tableMS1abundanceDataFrame <<- resultObj2$ms1abundanceDataFrame
          tableMS2fragmentDataFrame  <<- resultObj2$ms2fragmentDataFrame
          tableAnnotationDataFrame   <<- resultObj2$annotationDataFrame
          
          updateAnnoGui(resultObj$precursorSet)
        } else {
          treeNodeFromTreeSelection <<- NULL
        }
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings click", sep = ""))
          paste("Info for Precursor ", dataList$precursorLabels[[precursorIndex]], sep = "")
        })
        
        pcaLoadingFromPcaLoadingSelection <<- minimumIndex
        pcaLoadingSetFromMS2Fragment <<- NULL
        
        selectedFragmentIndex <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
        
        resultObj <- getMS2spectrumOfPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsColor <<- resultObj$fragmentColor
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentshoveredColor <<- NULL
      }
      
      ## MS2 plot
      resetMS2PlotRange()
      drawMS2Plot(consoleInfo = "PCA loadings click output$clusterPlotMS2")
      drawPcaLoadingsPlot(consoleInfo = "PCA loadings click output$pcaPlotLoadings")
      
      if(state$showHCAplotPanel)
        ## update dendrogram plot
        drawDendrogramPlot(consoleInfo = "PCA loadings click output$clusterPlotDendrogram")
    })
    obsPCAloadingsHover <- observeEvent(input$pcaPlotLoadings_hover, {
      hoverX <- input$pcaPlotLoadings_hover$x
      hoverY <- input$pcaPlotLoadings_hover$y
      plotWidth  <- session$clientData$output_pcaPlotLoadings_width
      plotHeight <- session$clientData$output_pcaPlotLoadings_height
      
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
      } else {
        print(paste("Observe PCA Loadings hover", hoverX, hoverY, minimumIndex))
        
        precursorIndex <- filterPca$filter[[minimumIndex]]
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings hover ", precursorIndex, sep = ""))
          paste("Info for Precursor ", dataList$precursorLabels[[precursorIndex]], sep = "")
        })
        
        resultObj <- getMS2spectrumOfPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        fragmentsXhovered <<- resultObj$fragmentMasses
        fragmentsYhovered <<- resultObj$fragmentAbundances
        fragmentshoveredColor <<- resultObj$fragmentColor
      }
      
      drawMS2Plot(consoleInfo = "loadings hover output$clusterPlotMS2")
    })
    obsPCAloadingsDblClick <- observeEvent(input$pcaPlotLoadings_dblclick, {
      brush <- input$pcaPlotLoadings_brush
      
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
      
      print(paste("observe pcaLoadingsPlotRange", paste(xInterval, collapse = ", "), paste(yInterval, collapse = ", ")))
      
      ## plot PCA
      drawPcaLoadingsPlot(consoleInfo = "range update output$pcaPlotLoadings")
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
      
      drawFragmentPlot(consoleInfo = "update range output$fragmentPlotDendrogram")
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
      if(drawPlots == removePresentAnnotationValue)
        return()
      removePresentAnnotationValue <<- drawPlots
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
      if(drawPlots == submitNewAnnotationValue)
        return()
      submitNewAnnotationValue <<- drawPlots
      print(paste("Observe submitNewAnnotation", drawPlots))
      
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
      
      vals <- getInputValues(id = "Ignore", len = nrow(table$dataFrame))
      
      if(all(is.na(vals)))
        return()
      
      nas <- is.na(vals)
      vals[nas] <- dataList$annoArrayIsArtifact[selectedPrecursorSet][nas]
      vals <- as.logical(vals)
      
      setArtifactState(selectedPrecursorSet, vals)
    })
    ## suspend observer
    session$onSessionEnded(function() {
      ## sidepanel
      obsTabs$suspend()
      obsFile$suspend()
      obsChangePlot$suspend()
      obsShowSideBar$suspend()
      obsShowLoadingsLabels$suspend()
      obsShowLoadingsAbundance$suspend()
      ## filter
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
      ## anno
      obsRemovePresentAnno$suspend()
      obsToggleAddNewAnnoButton$suspend()
      obsAddNewAnno$suspend()
      obsAddPresentAnno$suspend()
      obsIgnoreValueChanged$suspend()
    })
    #########################################################################################
    #########################################################################################
    ##direct output rendering
    output$fileInfo <- renderText({
      print(paste("init output$fileInfo"))
      paste("Please select a fragment matrix file in the right panel")
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
      return(!is.null(input$matrixFile))
    })
    output$analysisType <- reactive({
      print(paste("reactive update analysisType", state$analysisType))
      if(!is.null(state$analysisType)){
        if(state$analysisType == "HCA")
          state$plotToShow <<- "Display HCA"
        if(state$analysisType == "PCA")
          state$plotToShow <<- "Display PCA"
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
    output$plotToShow <- reactive({
      print(paste("reactive update plotToShow", state$plotToShow))
      return(state$plotToShow)
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
    ## download
    createExportMatrixName <- function(){
      fileName <- paste(gsub(" ", "_", gsub(":", ".", Sys.time())), "_selectedPrecursorMatrix.csv.gz", sep = "")
      return(fileName)
    }
    createExportMatrix <- function(precursorSet){
      #subMatrix <- matrix(data = dataList$featureMatrix[precursorSet, ], nrow = length(precursorSet))
      #
      #fragmentPresentinColumns <- apply(X = subMatrix, MARGIN = 2, FUN = sum) > 0
      #subMatrix <- subMatrix[, fragmentPresentinColumns]
      #
      #dataFrame <- as.data.frame(subMatrix)
      #
      #if(length(precursorSet) > 1){
      #  rownames(dataFrame) <- dataList$precursorLabels[precursorSet]
      #  colnames(dataFrame) <- dataList$fragmentMasses[fragmentPresentinColumns]
      #} else {
      #  names(dataFrame) <- dataList$precursorLabels[precursorSet]
      #}
      
      #dataFrame <- dataList$dataFrameOriginal[c(1:3, precursorSet + 3), ]
      
      fragmentMatrix      <- dataList$featureMatrix[precursorSet, ]
      fragmentCounts      <- apply(X = fragmentMatrix, MARGIN = 2, FUN = function(x){ sum(x != 0) })
      fragmentIntensities <- apply(X = fragmentMatrix, MARGIN = 2, FUN = function(x){ sum(x) }) / fragmentCounts
      fragmentMasses      <- dataList$fragmentMasses
      
      fragmentSelection   <- fragmentCounts != 0
      
      fragmentMatrix      <- fragmentMatrix[, fragmentSelection]
      fragmentCounts      <- fragmentCounts[fragmentSelection]
      fragmentIntensities <- fragmentIntensities[fragmentSelection]
      fragmentMasses      <- fragmentMasses[fragmentSelection]
      
      fragmentMatrix      <- as.matrix(x = fragmentMatrix)
      fragmentMatrix[fragmentMatrix == 0] <- ""
      ms2Matrix     <- rbind(
        fragmentCounts,
        fragmentIntensities,
        fragmentMasses,
        fragmentMatrix
      )
      ms1Matrix     <- rbind(
        dataList$dataFrameMS1Header,
        dataList$dataFrameInfos[precursorSet, ]
      )
      
      dataFrame <- cbind(
        ms1Matrix,
        ms2Matrix
      )
      
      #print(paste("v", ncol(dataList$featureMatrix), " --> n", ncol(ms2Matrix)))
      
      #as.data.frame(as.matrix(dataList$featureMatrix[precursorSet, ]))
      #dataFrame <- data.frame()
      
      return(dataFrame)
    }
    writeTable <- function(precursorSet, file){
      dataFrame <- createExportMatrix(file, precursorSet)
      gz1 <- gzfile(file, "w")
      write.table(x = dataFrame, file = gz1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      close(gz1)
    }
    ## download filtered
    output$downloadGlobalMS2filteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterGlobal$filter
        writeTable(precursorSet = precursorSet, file = file)
      }
    )
    output$downloadHcaFilteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterHca$filter
        writeTable(precursorSet = precursorSet, file = file)
      }
    )
    output$downloadPcaFilteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
        content = function(file) {
        precursorSet <- filterPca$filter
        writeTable(precursorSet = precursorSet, file = file)
      }
    )
    output$downloadSearchPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterPca$filter
        writeTable(precursorSet = precursorSet, file = file)
      }
    )
    ## download selected
    output$downloadHcaSelectedPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        ## get selected precursors
        if(is.null(treeNodeFromTreeSelection)){
          ## all precursors
          precursorSet <- filterHca$filter
        } else {
          if(treeNodeFromTreeSelection < 0){
            ## node
            precursorSet <- filterHca$filter[-treeNodeFromTreeSelection]
          } else {
            ## selected precursors
            precursorSet <- clusterDataList$innerNodeMembers[[treeNodeFromTreeSelection]]
          }
        }
        
        writeTable(precursorSet = precursorSet, file = file)
      }
    )
    
    #########################################################################################
    #########################################################################################
    ## properties
    options(shiny.maxRequestSize=500*1024^2) 
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
    outputOptions(output, 'plotHcaShown', suspendWhenHidden=FALSE)
    outputOptions(output, 'plotPcaShown', suspendWhenHidden=FALSE)
    
    #plotToShow = "Display HCA", 
    
    #########################################################################################
    #########################################################################################
    ## ui generation
    output$runRightColumn <- renderUI({
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
                    checkboxInput(inputId = "showSideBar", label = "Display side bar", value = state$showSideBar)
                )
             ),##column
             column(width = 6,
                conditionalPanel(
                  condition = "output.showHCAplotPanel && output.showPCAplotPanel",
                  div(style="float:left",
                      bsTooltip(id = "changePlot", title = "Switch between HCA plot and PCA plot", placement = "bottom", trigger = "hover"),
                      radioButtons(inputId = "changePlot", label = NULL, choices = c("Display HCA", "Display PCA"), inline = TRUE, selected = state$plotToShow)
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
                             outputId = "clusterPlotDendrogram", 
                             #hover    = "clusterPlotDendrogram_hover", 
                             hover    = hoverOpts(
                               id = "clusterPlotDendrogram_hover",
                               delay = 50, 
                               delayType = "debounce"
                             ),
                             click    = "clusterPlotDendrogram_click",
                             dblclick = "clusterPlotDendrogram_dblclick",
                             #brush    = "clusterPlotDendrogram_brush"
                             brush    = brushOpts(
                               id = "clusterPlotDendrogram_brush",
                               resetOnNew = TRUE,
                               direction = "x",
                               delay = 00,
                               delayType = "debounce"
                             )
                  ),
                  plotOutput(height = 75, 
                             outputId = "clusterPlotHeatmap",
                             #hover    = "clusterPlotHeatmap_hover", 
                             hover    = hoverOpts(
                               id = "clusterPlotHeatmap_hover",
                               delay = 50, 
                               delayType = "debounce"
                             )
                             #click = "clusterPlotHeatmap_click"
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
                               outputId = "pcaPlotScores", 
                               #hover    = "pcaPlotScores_hover",
                               hover    = hoverOpts(
                                 id = "pcaPlotScores_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "pcaPlotScores_click",
                               dblclick = "pcaPlotScores_dblclick",
                               #brush    = "pcaPlotScores_brush"
                               brush = brushOpts(
                                 id = "pcaPlotScores_brush",
                                 resetOnNew = TRUE,
                                 direction = "xy",
                                 delay = 00,
                                 delayType = "debounce"
                               )
                    )
                 ),## column
                 column(width = 6,
                    plotOutput(height = 500, 
                               outputId = "pcaPlotLoadings", 
                               #hover    = "pcaPlotLoadings_hover",
                               hover    = hoverOpts(
                                 id = "pcaPlotLoadings_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "pcaPlotLoadings_click",
                               dblclick = "pcaPlotLoadings_dblclick",
                               #brush    = "pcaPlotScores_brush"
                               brush = brushOpts(
                                 id = "pcaPlotLoadings_brush",
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
              conditionalPanel(
                condition = '(output.analysisType == "HCA" && output.showHCAplotPanel) || (output.analysisType == "PCA" && output.showPCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "calcClusterPlotAnnoLegend", height = 200)
                )
              ),## conditional
              conditionalPanel(
                condition = 'output.analysisType == "HCA" && output.showHCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "clusterPlotHeatmapLegend", height = 300)
                )
              ),## conditional
              conditionalPanel(## loadings properties
                condition = 'output.analysisType == "PCA" && output.showPCAplotPanel',
                bsTooltip(id = "showLoadingsLabels", title = "Display loadings labels", placement = "bottom", trigger = "hover"),
                checkboxInput(inputId = "showLoadingsLabels", label = "Show labels", value = FALSE),
                bsTooltip(id = "showLoadingsAbundance", title = "Size of loadings by abundance in MS1", placement = "bottom", trigger = "hover"),
                checkboxInput(inputId = "showLoadingsAbundance", label = "Show abundance", value = FALSE)
              ),
              conditionalPanel(
                condition = '(output.analysisType == "HCA" && output.showHCAplotPanel) || (output.analysisType == "PCA" && output.showPCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "calcMS2PlotAnnoLegend", height = 100)
                )
              )## conditional
           )## column
         ),## row
         #########################################################################################
         ## MS2 plot and info
         conditionalPanel(
           condition = '(output.analysisType == "HCA" && output.showHCAplotPanel) || (output.analysisType == "PCA" && output.showPCAplotPanel)',
           fluidRow(
             plotOutput(height = 250, 
                        outputId = "clusterPlotMS2",
                        #hover    = "clusterPlotMS2_hover",
                        hover    = hoverOpts(
                          id = "clusterPlotMS2_hover",
                          delay = 50, 
                          delayType = "debounce"
                        ),
                        click    = "clusterPlotMS2_click",
                        dblclick = "clusterPlotMS2_dblclick",
                        #brush    = "clusterPlotMS2_brush",
                        brush = brushOpts(
                          id = "clusterPlotMS2_brush",
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
               bsTooltip(id = "information", title = "Information about selected items in the plot", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("information"),
               htmlOutput(outputId = "metFragLink"),
               h4("Tip"),
               bsTooltip(id = "tip", title = "Information about operating options", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("tip"),
               bsTooltip(id = "downloadHcaSelectedPrecursors", title = "Download the set of selected precursors", placement = "bottom", trigger = "hover"),
               downloadButton('downloadHcaSelectedPrecursors', 'Download selected precursors')
             )## well
           )## row
         ),## conditional
         ##############################################################################################
         ## precursor set selection and annotation
         conditionalPanel(
           condition = 'output.precursorSetSelected && output.analysisType == "HCA"',
           fluidRow(
             wellPanel(
               ## annotation
               h4("Annotation"),
               h5("Present annotation(s)"),
               fluidRow(
                 column(
                   width = 3,
                   bsTooltip(id = "presentAnnotationValue", title = "The set of present annotations for the set of selected precursors", placement = "bottom", trigger = "hover"),
                   selectInput(inputId = "presentAnnotationValue", label = NULL, choices = c(""), selectize = FALSE)
                 ),## column
                 column(
                   width = 3,
                   bsTooltip(id = "setPresentAnnotationPrimary", title = "Sets the selected annotation primary for the set of selected precursors; i.e. this annotation will be used for coloring", placement = "bottom", trigger = "hover"),
                   actionButton(inputId = "setPresentAnnotationPrimary", label = "Set primary")
                 ),## column
                 column(
                   width = 6,
                   bsTooltip(id = "removePresentAnnotation", title = "Removes the selected annotation from the set of selected precursors", placement = "bottom", trigger = "hover"),
                   actionButton(inputId = "removePresentAnnotation", label = "Remove annotation")
                 )## column
               ),##row
               h4("Add annotation"),
               fluidRow(
                 column(
                   width = 4,
                   bsTooltip(id = "annotationSelection", title = "The user is able to add a new annotation or an annotation which has been used before", placement = "bottom", trigger = "hover"),
                   radioButtons(inputId = "annotationSelection", label = NULL, choices = c("New annotation", "Previous annotation"), selected = NULL, inline = FALSE)
                 ),## column
                 column(
                   width = 8,
                   #h4("Add new annotation"),
                   conditionalPanel(
                     condition = 'input.annotationSelection == "New annotation"',
                     bsTooltip(id = "newAnnotationValue", title = "The name of this annotation", placement = "bottom", trigger = "hover"),
                     textInput(inputId = "newAnnotationValue", label = "Type new annotation"),
                     #colourInput(inputId = "newAnnotationColor", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = c("blue")),
                     bsTooltip(id = "newAnnotationColor", title = "The color of this annotation", placement = "bottom", trigger = "hover"),
                     colourInput(inputId = "newAnnotationColor", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = colorPalette()),
                     bsTooltip(id = "submitNewAnnotation", title = "Adds this annotation to the set of selected precursors", placement = "bottom", trigger = "hover"),
                     actionButton(inputId = "submitNewAnnotation", label = "Add annotation")
                   ),## conditional
                   conditionalPanel(
                     condition = 'input.annotationSelection == "Previous annotation"',
                     bsTooltip(id = "previousAnnotationValue", title = "The set of annotations which have been added before", placement = "bottom", trigger = "hover"),
                     selectInput(inputId = "previousAnnotationValue", label = NULL, choices = c("Artifact"), selectize = FALSE),
                     bsTooltip(id = "submitPreviousAnnotation", title = "Adds this annotation to the set of selected precursors", placement = "bottom", trigger = "hover"),
                     actionButton(inputId = "submitPreviousAnnotation", label = "Add annotation")
                   )## conditional
                 )## column
               )## row
             )## well
           ),## row
           fluidRow(
             wellPanel(
               h4("Selected precursors"),
               DT::dataTableOutput("table"),
               bsTooltip(id = "updateArtifactsFromCheckboxes", title = "Adds the annotation \\'ignore\\' to the set of checked precursors", placement = "bottom", trigger = "hover"),
               actionButton(inputId = "updateArtifactsFromCheckboxes", label = "Update ignored precursors")
             )## well
           )## row
         )## conditional
      )##column
    })
  }## function(input, output, session)
)## shinyServer
