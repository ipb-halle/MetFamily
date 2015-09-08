
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#


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
    ## data
    dataList <- NULL
    currentDistanceMeasure <- NULL
    filter <- NULL
    clusterDataList <- NULL
    pcaDataList <- NULL
    ## program state
    state <- reactiveValues(
      filterValid = FALSE, 
      showHCAplotPanel = FALSE, 
      showPCAplotPanel = FALSE, 
      precursorSetSelected = FALSE,
      analysisType = NULL
    )
    ## selection
    selectedFragmentIndex <- NULL
    treeNodeFromTreeSelection <- NULL
    treeNodeSetFromMS2Fragment <- NULL
    pcaLoadingSetFromMS2Fragment <- NULL
    pcaLoadingFromPcaLoadingSelection <- NULL
    ## buttons
    submitButtonValue <- 0
    drawHCAButtonValue <- 0
    drawPCAButtonValue <- 0
    downloadMatrixButtonValue <- 0
    removePresentAnnotationValue <- 0
    submitNewAnnotationValue <- 0
    submitPreviousAnnotationValue <- 0
    updateArtifactsFromCheckboxesButtonValue <- 0
    
    ## MS2 peaks
    fragmentsX <- NULL
    fragmentsY <- NULL
    fragmentsXhovered <- NULL
    fragmentsYhovered <- NULL
    ## plot ranges
    dendrogramPlotRange  <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
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
    table <- reactiveValues(dataFrame = NULL)
    selectedPrecursorSet <- NULL
    
    #########################################################################################
    #########################################################################################
    ## functions
    disableActionButton <- function(session, id) {
      session$sendCustomMessage(type="jsCode", list(code = paste(
        "$('#", id, "').prop('disabled', true)", sep=""
      )))
    }
    enableActionButton <- function(session, id) {
      session$sendCustomMessage(type="jsCode", list(code = paste(
        "$('#", id, "').prop('disabled', false)", sep=""
      )))
    }
    getSelectedPOI_X <- function(mouseX, poiCoordinatesX, plotWidth, plotRangeX){
      factorX <- plotWidth  / plotRangeX
      
      mouseX <- mouseX * factorX
      poiCoordinatesX <- poiCoordinatesX * factorX
      
      distances <- abs(poiCoordinatesX - mouseX)
      distanceThreshold <- factorX * plotRangeX / 50
      
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
      distanceThreshold <- factorX * plotRangeX / 50
      
      minimumIndex <- which.min(distances)
      minimumDistance <- distances[[minimumIndex]]
      
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    #########################################################################################
    ##' Parse the input file
    getClusterData <- function(file){
      print(paste("getClusterData", is.null(file)))
      
      doClearPlots()
      
      clusterDataList <<- NULL
      pcaDataList <<- NULL
      
      state$analysisType <<- NULL
      state$filterValid <<- FALSE
      state$showHCAplotPanel <<- FALSE
      state$showPCAplotPanel <<- FALSE
      state$precursorSetSelected <<- FALSE
      
      selectedFragmentIndex <<- NULL
      treeNodeFromTreeSelection <<- NULL
      treeNodeSetFromMS2Fragment <<- NULL
      pcaLoadingSetFromMS2Fragment <<- NULL
      pcaLoadingFromPcaLoadingSelection <<- NULL
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      fragmentsXhovered <<- NULL
      fragmentsYhovered <<- NULL
      
      #ms2PlotRange$xMin <<- NULL
      #ms2PlotRange$xMax <<- NULL
      #ms2PlotRange$xInterval <<- NULL
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
      
      table$dataFrame <<- NULL
      selectedPrecursorSet <<- NULL
      
      withProgress(message = 'Reading file...', value = 0, {
        #dataList <<- calcClusterData(file = "/mnt/VOL1/ABT/Alle/Balcke/MetSWATH/data/MS-DIAL/UC Davis/Results/201558139_matrixPrecursorsVersusFragmentsDeisotoped_withoutZerosTest01.txt")
        dataList <<- readClusterData(file = file, progress = TRUE)
      })
      print(paste("getClusterData do data finished", dataList$minimumMass))
      
      ## update input values
      switch(as.character(length(dataList$groups)), 
        "0"={
          groupOne <- NA
          groupTwo <- NA
          selectedOne <- NULL
          selectedTwo <- NULL
        },
        "1"={
          groupOne <- dataList$groups[[1]]
          groupTwo <- dataList$groups[[1]]
          selectedOne <- dataList$groups[[1]]
          selectedTwo <- dataList$groups[[1]]
        },
        {
          groupOne <- dataList$groups[[1]]
          groupTwo <- dataList$groups[[2]]
          selectedOne <- dataList$groups[[1]]
          selectedTwo <- dataList$groups[[2]]
        }
      )
      
      updateSelectInput(session = session, inputId = "groupOne", choices = dataList$groups, selected = selectedOne)
      updateSelectInput(session = session, inputId = "groupTwo", choices = dataList$groups, selected = selectedTwo)
      updateSelectInput(session = session, inputId = "groups",   choices = dataList$groups, selected = dataList$groups)
      updateTextInput(session = session, inputId = "filter_average", value = "0")
      updateTextInput(session = session, inputId = "filter_lfc", value = "0")
      updateTextInput(session = session, inputId = "filter_ms2_masses1", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_ppm", value = "20")
      
      #updateColourInput(session = session, inputId = "newAnnotationColor", allowedCols = colorPalette())
      
      ## update filter
      filter <<- filterData(dataList = dataList, 
                            groups = dataList$groups, filter_average = NULL, filter_lfc = NULL, 
                            filterList_ms2_masses = NULL, filter_ms2_ppm = NULL, 
                            includeIgnoredPrecursors = TRUE, progress = FALSE)
      
      output$information <- renderText({
        print(paste("init output$information", sep = ""))
        paste("Number of filtered precursors: ", dataList$numberOfPrecursorsFiltered, sep = "")
      })
      
      ms2PlotRange$xMin <<- dataList$minimumMass
      ms2PlotRange$xMax <<- dataList$maximumMass
      ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
    }
    #########################################################################################
    ##' perform filtering
    doPerformFiltering <- function(){
      #################################################
      ## get inputs
      groupOne <- input$groupOne
      groupTwo <- input$groupTwo
      groups <- input$groups
      filter_average <- input$filter_average
      filter_lfc <- input$filter_lfc
      filter_ms2_masses1 <- input$filter_ms2_masses1
      filter_ms2_masses2 <- input$filter_ms2_masses2
      filter_ms2_masses3 <- input$filter_ms2_masses3
      filter_ms2_ppm <- input$filter_ms2_ppm
      includeIgnoredPrecursors <- input$filterIncludeIgnoredPrecursors
      #analysisType <- input$analysisType
      print(paste("Observe applyFilters1", "gs", paste(groups, collapse = "-"), "g1", groupOne, "g2", groupTwo, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "ig", includeIgnoredPrecursors))
      
      #if(analysisType == "PCA")
      #  filter_lfc <- NULL
      
      #################################################
      ## parse inputs
      if(is.null(filter_average) | nchar(filter_average) == 0)
        filter_average <- NULL
      else
        filter_average <- as.numeric(filter_average)
      
      if(is.null(filter_lfc))
        filter_lfc <- NULL
      else{
        if(nchar(filter_lfc) == 0)
          filter_lfc <- NULL
        else
          filter_lfc <- as.numeric(filter_lfc)
      }
      
      if(is.null(filter_ms2_masses1) | nchar(filter_ms2_masses1) == 0)
        filter_ms2_masses1 <- NULL
      if(is.null(filter_ms2_masses2) | nchar(filter_ms2_masses2) == 0)
        filter_ms2_masses2 <- NULL
      if(is.null(filter_ms2_masses3) | nchar(filter_ms2_masses3) == 0)
        filter_ms2_masses3 <- NULL
      
      if(is.null(filter_ms2_ppm) | nchar(filter_ms2_ppm) == 0)
        filter_ms2_ppm <- NULL
      else
        filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
      
      print(paste("Observe applyFilters2", "gs", paste(groups, collapse = "-"), "g1", groupOne, "g2", groupTwo, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm))
      print(paste("Observe applyFilters2", "gs", is.null(groups), "g1", is.null(groupOne), "g2", is.null(groupTwo), "a", is.null(filter_average), "lfc", is.null(filter_lfc), "ms2 1", is.null(filter_ms2_masses1), "ms2 2", is.null(filter_ms2_masses2), "ms2 3", is.null(filter_ms2_masses3), "ppm", is.null(filter_ms2_ppm)))
      
      #################################################
      ## check for errors in inputs amd process ms2
      error <- FALSE
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
      error <- error | ((!is.null(filter_ms2_masses1) | !is.null(filter_ms2_masses2) | !is.null(filter_ms2_masses3)) & is.null(filter_ms2_ppm))
      
      print(paste("Observe applyFilters3", "e", error, "gs", paste(groups, collapse = "-"), "g1", groupOne, "g2", groupTwo, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm))
      if(error){
        output$filteredPrecursors <- renderText({
          print(paste("update output$filteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
        return()
      }
      
      filterList_ms2_masses <- list()
      if(!is.null(filter_ms2_masses1))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses1
      if(!is.null(filter_ms2_masses2))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses2
      if(!is.null(filter_ms2_masses3))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses3
      
      #################################################
      ## do filtering
      #if(analysisType == "HCA"){
        groupSet <- c(groupOne, groupTwo)
      #} else {
      #  groupSet <- groups
      #}
      
      filterHere <- filterData(
        dataList = dataList, 
        #groupOne = groupOne, groupTwo = groupTwo, 
        groups = groupSet, filter_average = filter_average, filter_lfc = filter_lfc, 
        filterList_ms2_masses = filterList_ms2_masses, filter_ms2_ppm = filter_ms2_ppm, includeIgnoredPrecursors = includeIgnoredPrecursors,
        progress = FALSE
      )
      numberOfPrecursorsFiltered <- filterHere$numberOfPrecursorsFiltered
      print(paste("Observe applyFilters", numberOfPrecursorsFiltered))
      
      #################################################
      ## check filter validity
      minimumNumberOfPrecursors <- 6
      maximumNumberOfPrecursors <- 1000
      if(numberOfPrecursorsFiltered >= minimumNumberOfPrecursors & numberOfPrecursorsFiltered <= maximumNumberOfPrecursors){
        ## filter valid
        print(paste("Observe applyFilters ", minimumNumberOfPrecursors, " <= # <= ", maximumNumberOfPrecursors, sep = ""))
        output$filteredPrecursors <- renderText({
          print(paste("update output$filteredPrecursors ", minimumNumberOfPrecursors, " <= # <= ", maximumNumberOfPrecursors, sep = ""))
          paste("Number of filtered precursors: ", numberOfPrecursorsFiltered, sep = "")
        })
        
        enableActionButton(session, "drawHCAplots")
        enableActionButton(session, "drawPCAplots")
        filter <<- filterHere
        state$filterValid <<- TRUE
        state$showHCAplotPanel <<- FALSE
        state$showPCAplotPanel <<- FALSE
        state$precursorSetSelected <<- FALSE
        
        selectedFragmentIndex <<- NULL
        treeNodeFromTreeSelection <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
        pcaLoadingSetFromMS2Fragment <<- NULL
        pcaLoadingFromPcaLoadingSelection <<- NULL
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        doClearPlots()
      } else {
        ## filter invalid
        #filter <<- NULL
        disableActionButton(session, "drawHCAplots")
        disableActionButton(session, "drawPCAplots")
        #state$filterValid <<- FALSE
        
        ## update info
        if(numberOfPrecursorsFiltered == 0){
          print(paste("Observe applyFilters # = 0", sep = ""))
          output$filteredPrecursors <- renderText({
            print(paste("update output$filteredPrecursors # = 0", sep = ""))
            paste("There are no precursors which fulfill the given criteria.", sep = "")
          })
        }
        if(numberOfPrecursorsFiltered > 0 & numberOfPrecursorsFiltered < minimumNumberOfPrecursors){
          print(paste("Observe applyFilters 0 < # < ", minimumNumberOfPrecursors, sep = ""))
          output$filteredPrecursors <- renderText({
            print(paste("update output$filteredPrecursors 0 < # < ", minimumNumberOfPrecursors, sep = ""))
            paste("There are only ", numberOfPrecursorsFiltered, " precursors which fulfill the given criteria. There must be at least more than five precursors to preceed.", sep = "")
          })
        }
        if(numberOfPrecursorsFiltered > maximumNumberOfPrecursors){
          print(paste("Observe applyFilters # > ", maximumNumberOfPrecursors, sep = ""))
          output$filteredPrecursors <- renderText({
            print(paste("update output$filteredPrecursors # > ", maximumNumberOfPrecursors, sep = ""))
            paste("There are ", numberOfPrecursorsFiltered, " precursors which fulfill the given criteria. There must be at most ", maximumNumberOfPrecursors, " precursors to preceed.", sep = "")
          })
        }
      }
    }
    #########################################################################################
    ##' Clear plots
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
       })
    }
    #########################################################################################
    ##' annotation stuff
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
      updateAnnoGui(precursorSet)
      
      ## TODO
      
      #dataList$annoArrayOfLists
      #dataList$annoArrayIsArtifact
      #dataList$annoPresentAnnotationsList
      #dataList$annoPresentColorsList
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
      updateAnnoGui(precursorSet)
      
      ## TODO
    }
    setArtifactState <- function(precursorSet, isArtifact){
      ## add
      dataList$annoArrayIsArtifact[precursorSet] <<- isArtifact
      
      for(i in 1:length(precursorSet))
        print(paste(dataList$precursorLabels[[precursorSet[[i]]]], isArtifact[[i]]))
      
      ## update gui
      updateAnnoGui(precursorSet)
      
      ## TODO
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
      
      if(length(commonAnnos) > 0)
        enableActionButton(session, "removePresentAnnotation")
      else
        disableActionButton(session, "removePresentAnnotation")
      
      ## table
    }
    #########################################################################################
    ##' table functions
    createInputFields<-function(FUN, id, values) {
      ## create a character vector of shiny inputs
      inputs <- character(length(values))
      for (i in 1:length(values)){
        itemId    <- paste(id, i, sep = "")
        inputs[i] <- as.character(FUN(inputId = itemId, label = NULL, value = values[[i]]))
      }
      return(inputs)
    }
    getInputValues<-function(id, len) {
      ## obtain the values of inputs
      unlist(lapply(1:len, function(i) {
        itemId <- paste(id, i, sep = "")
        value  <- input[[itemId]]
        if (is.null(value))
          return(NA)
        else
          return(value)
      }))
    }
    #########################################################################################
    #########################################################################################
    ## observer
    #########################################################################################
    ##' listen to window resize events
    #obsClusterPlotResize <- observe({
    #  plotWidth  <- session$clientData$output_clusterPlotDendrogram_width
    #  plotHeight <- session$clientData$output_clusterPlotDendrogram_height
    #  print(paste("Observe dendrogram plot resize", plotWidth, plotHeight, is.null(dataList)))
    #  
    #  if(!is.null(dataList) & !is.null(clusterDataList))
    #    doCalculatePlots()
    #})
    #########################################################################################
    ##' listen to filtab selections
    obsTabs <- observeEvent(input$runTabs, {
      tabId <- input$runTabs
      print(paste("Observe tab selection", tabId))
      if(tabId == "HCA")
        state$analysisType <<- "HCA"
      if(tabId == "PCA")
        state$analysisType <<- "PCA"
    })
    #########################################################################################
    ##' listen to file input
    obsFile <- observeEvent(input$matrixFile$datapath, {
      file     <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      print(paste("Observe file for data", fileName))
      if(!is.null(file)){
        output$fileInfo <- renderText({paste("Processing file ", fileName, "...", sep = "")})
        getClusterData(file)
        output$fileInfo <- renderText({fileName})
      }
    })
    #########################################################################################
    ##' listen to submit button events
    obsApplyFilters <- observeEvent(input$applyFilters, {
      applyFilters <- as.numeric(input$applyFilters)
      
      print(paste("Observe applyFilters", applyFilters))
      
      #################################################
      ## check if button was hit
      if(applyFilters == submitButtonValue)
        return()
      submitButtonValue <<- applyFilters
      
      #################################################
      ## do filtering
      doPerformFiltering()
    })
    #########################################################################################
    ##' listen to draw HCA button events
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
      
      distance <- input$distanceFunction
      clusterMethod <- input$clusterMethod
      
      ## reset selections
      treeNodeFromTreeSelection <<- NULL
      treeNodeSetFromMS2Fragment <<- NULL
      
      dendrogramPlotRange$xMin <<- 1
      dendrogramPlotRange$xMax <<- filter$numberOfPrecursorsFiltered
      dendrogramPlotRange$xInterval <<- c(1, filter$numberOfPrecursorsFiltered)
      dendrogramPlotRange$xIntervalSize <<- filter$numberOfPrecursorsFiltered - 1
      
      ##########################
      ## draw
      
      ## update info
      print(paste("Observe do draw HCA plots", sep = ""))
      
      output$information <- renderText({
        print(paste("update output$information do draw HCA plots", sep = ""))
        paste("Number of filtered precursors: ", length(filter$filter), sep = "")
      })
      
      ## compute distance matrix
      print(paste("Calculating distances", sep = ""))
      currentDistanceMeasure <<- distance
      withProgress(message = 'Calculating distances...', value = 0, {
        distanceMatrix <- calculateDistanceMatrix(dataList = dataList, filter = filter, distance = currentDistanceMeasure, progress = TRUE)
      })
      ## compute cluster
      print(paste("Calculating cluster", sep = ""))
      withProgress(message = 'Calculating cluster...', value = 0, {
        clusterDataList <<- calculateCluster(progress = TRUE,
                                           dataList = dataList, filter = filter$filter, distanceMatrix = distanceMatrix, method = clusterMethod
        )
      })
      ## clear info
      output$information <- renderText({
        print(paste("update output$information clear"))
        paste("", sep = "")
      })
      output$tip <- renderText({
        print(paste("update output$tip clear"))
        paste("", sep = "")
      })
      
      ## compute plots
      print(paste("Calculating plots", sep = ""))
      
      if(!is.null(selectedFragmentIndex)){
        ## tree highlighting
        fragmentMass = fragmentsX[[selectedFragmentIndex]]
        yesNoFunction <- function(precursorIndex){
          features <- dataList$featureIndeces[[precursorIndex]]
          fragmentMasses <- dataList$fragmentMasses[features]
          comprisesFragmentMass <- any(fragmentMasses == fragmentMass)
          return(comprisesFragmentMass)
        }
        
        #treeNodeSetFromMS2Fragment <<- getSetOfSubTreesFromRootForMass(dataList = dataList, filter = filter, clusterDataList = clusterDataList, fragmentMass = fragmentsX[[selectedFragmentIndex]])
        treeNodeSetFromMS2Fragment <<- getSetOfSubTreesFromRoot(filter = filter, clusterDataList = clusterDataList, yesNoFunction = yesNoFunction)
      }
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(dataList = dataList, filter = filter$filter, clusterDataList = clusterDataList, annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, annoPresentColorsList = dataList$annoPresentColorsList, distanceMeasure = currentDistanceMeasure, nodeIndex = treeNodeFromTreeSelection, nodeIndeces = treeNodeSetFromMS2Fragment)
      })
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("output$clusterPlotHeatmap"))
        calcClusterPlotHeatmap(dataList = dataList, filter = filter, clusterDataList = clusterDataList)
      })
      output$clusterPlotHeatmapLegend <- renderPlot({
        print(paste("output$clusterPlotHeatmapLegend"))
        calcClusterPlotHeatmapLegend(dataList = dataList)
      })
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
      state$showHCAplotPanel <<- TRUE
    })
    #########################################################################################
    ##' listen to draw PCA button events
    obsDrawPCA <- observeEvent(input$drawPCAplots, {
      drawPlots <- as.numeric(input$drawPCAplots)
      
      print(paste("Observe draw PCA plots", drawPlots))
      
      #################################################
      ## check if button was hit
      if(drawPlots == drawPCAButtonValue)
        return()
      drawPCAButtonValue <<- drawPlots
      
      pcaScaling      <- input$pcaScaling
      #pcaMeanCentered <- FALSE
      #pcaMeanCentered <- input$pcaMeanCentered
      pcaLogTransform <- input$pcaLogTransform
      #pcaMethod       <- input$pcaMethod
      
      #scalingMethods  <- c("none", "pareto", "vector", "uv")
      #pcaScaling      <- scalingMethods[[match(x = pcaScaling, table = c("None", "Pareto", "Vector normalization", "Unit variance"))]]
      
      pcaDimensionOne <<- as.numeric(input$pcaDimensionOne)
      pcaDimensionTwo <<- as.numeric(input$pcaDimensionTwo)
      print(paste("Observe draw PCA plots", "SC", pcaScaling, "MC", "LT", pcaLogTransform, "D1", pcaDimensionOne, "D2", pcaDimensionTwo))
      
      ## calc PCA
      #pca <<- calcPCA(dataList = dataList, filter = filter, unitVariance = unitVariance, logTransform = logTransform)
      #pca <- calcPCA(dataList = dataList, filter = filter, scaling = pcaScaling, meanCentered = pcaMeanCentered, logTransform = pcaLogTransform, method = pcaMethod)
      pca <- calcPCA(dataList = dataList, filter = filter, scaling = pcaScaling, logTransform = pcaLogTransform)
      
      pcaDataList <<- list()
      pcaDataList$pcaObj <<- pca
      pcaDataList$pcaScoresX <<- pca$scores[, pcaDimensionOne]
      pcaDataList$pcaScoresY <<- pca$scores[, pcaDimensionTwo]
      pcaDataList$pcaLoadingsX <<- pca$loadings[, pcaDimensionOne]
      pcaDataList$pcaLoadingsY <<- pca$loadings[, pcaDimensionTwo]
      pcaDataList$dimensionOne <<- pcaDimensionOne
      pcaDataList$dimensionTwo <<- pcaDimensionTwo
      
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
      if(!is.null(selectedFragmentIndex)){
        ## loadings highlighting
        fragmentMass     <- fragmentsX[[selectedFragmentIndex]]
        fragmentIndex    <- which(dataList$fragmentMasses == fragmentMass)
        precursorMatches <- dataList$featureMatrix[filter$filter, fragmentIndex] != 0
        pcaLoadingSetFromMS2Fragment <<- precursorMatches
      }
      
      ## plot PCA
      output$pcaPlotScores <- renderPlot({
        print(paste("output$pcaPlotScores"))
        plotPCAscores(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo)
      })
      output$pcaPlotLoadings <- renderPlot({
        print(paste("output$pcaPlotLoadings"))
        plotPCAloadings(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, pcaLoading = pcaLoadingFromPcaLoadingSelection, pcaLoadingSet = pcaLoadingSetFromMS2Fragment)
      })
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
      state$showPCAplotPanel <<- TRUE
    })
    #########################################################################################
    ##' listen to table events
    obsIgnoreValueChanged <- observeEvent(input$updateArtifactsFromCheckboxes, {
      updateArtifactsFromCheckboxes <- as.numeric(input$updateArtifactsFromCheckboxes)
      
      print(paste("Observe update artifact", updateArtifactsFromCheckboxes))
      
      #################################################
      ## check if button was hit
      if(updateArtifactsFromCheckboxes == updateArtifactsFromCheckboxesButtonValue)
        return()
      updateArtifactsFromCheckboxesButtonValue <<- updateArtifactsFromCheckboxes
      
      vals <- getInputValues(id = "Ignore_", len = nrow(table$dataFrame))
      
      if(any(is.na(vals)))
        return()
      
      vals <- as.logical(vals)
      setArtifactState(selectedPrecursorSet, vals)
    })
    #########################################################################################
    ##' listen to dendrogram plot mouse events
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
        
        #################################################
        ## fetch ms2 spectrum
        resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = minimumLabel)
        fragmentsXhovered <<- resultObj$fragmentsX
        fragmentsYhovered <<- resultObj$fragmentsY
        
        #################################################
        ## output as message and plots
        output$information <- renderText({
          print(paste("update output$information", resultObj$infoText))
          paste(resultObj$infoText, sep = "; ")
        })
      }
      
      ## MS2 plot
      output$clusterPlotMS2 <- renderPlot({
        print(paste("update output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
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
        selectedFragmentIndex <<- NULL
        treeNodeFromTreeSelection <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        
        pcaLoadingSetFromMS2Fragment <<- NULL
        pcaLoadingFromPcaLoadingSelection <<- NULL
        output$table <- DT::renderDataTable(NULL)
        table$dataFrame <<- NULL
        selectedPrecursorSet <<- NULL
        state$precursorSetSelected <<- FALSE
        updateAnnoGui(NULL)
      } else {
        ## tree selection
        minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
        #minimumText <- clusterDataList$poiText[[minimumIndex]]
        
        #################################################
        ## fetch ms2 spectrum
        selectedFragmentIndex <<- NULL
        treeNodeFromTreeSelection <<- minimumLabel
        treeNodeSetFromMS2Fragment <<- NULL
        resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, clusterLabel = minimumLabel)
        fragmentsX <<- resultObj$fragmentsX
        fragmentsY <<- resultObj$fragmentsY
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        
        pcaLoadingSetFromMS2Fragment <<- NULL
        pcaLoadingFromPcaLoadingSelection <<- NULL
        
        #################################################
        ## output as message
        output$information <- renderText({
          print(paste("update output$information", resultObj$infoText))
          paste(resultObj$infoText, sep = "; ")
        })
        
        ## MetFrag link
        if(!is.null(resultObj$landingPageUrl))
          output$metFragLink <- renderText({
            print(paste("update output$metFragLink", resultObj$landingPageUrl))
            paste("<a href=", gsub(pattern = " ", replacement = "%20", x = resultObj$landingPageUrl)[[1]], " target=\"_blank\">Send to MetFrag!</a>", sep = "")
          })
        else
          output$metFragLink <- renderText({
            print(paste("update output$metFragLink", resultObj$landingPageUrl))
            paste("", sep = "")
          })
        
        ## table
        isArtifact <- dataList$annoArrayIsArtifact[resultObj$precursorSet]
        ## TODO
        output$table <- DT::renderDataTable(
          expr = cbind(
            table$dataFrame,
            data.frame(Ignore = createInputFields(FUN = checkboxInput, id = 'Ignore_', values = isArtifact))
          ),
          server=FALSE, escape=FALSE, selection = "none",
          options=list(
            preDrawCallback=JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
            drawCallback=JS('function() { Shiny.bindAll(this.api().table().node()); } ')
          )
        )
        #input$table_rows_selected
        state$precursorSetSelected <<- TRUE
        selectedPrecursorSet <<- resultObj$precursorSet
        updateAnnoGui(resultObj$precursorSet)
        table$dataFrame <<- resultObj$dataFrame
      }
      
      #################################################
      ## plots
      
      ## cluster dendrogram
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("update output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(dataList = dataList, filter = filter$filter, clusterDataList = clusterDataList, annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, annoPresentColorsList = dataList$annoPresentColorsList, distanceMeasure = currentDistanceMeasure, nodeIndex = treeNodeFromTreeSelection, nodeIndeces = treeNodeSetFromMS2Fragment, xInterval = dendrogramPlotRange$xInterval)
      })
      
      ## MS2 plot
      ms2PlotRange$xMin <<- dataList$minimumMass
      ms2PlotRange$xMax <<- dataList$maximumMass
      ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("update output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
      
      if(state$showPCAplotPanel){
        ## update PCA plots
        output$pcaPlotScores <- renderPlot({
          print(paste("output$pcaPlotScores"))
          plotPCAscores(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo)
        })
        output$pcaPlotLoadings <- renderPlot({
          print(paste("output$pcaPlotLoadings"))
          plotPCAloadings(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, pcaLoading = pcaLoadingFromPcaLoadingSelection, pcaLoadingSet = pcaLoadingSetFromMS2Fragment)
        })
      }
    })
    obsDendrogramdblClick <- observeEvent(input$clusterPlotDendrogram_dblclick, {
      print("huhu")
      brush <- input$clusterPlotDendrogram_brush
      
      print(paste("observe dendrogram dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        dendrogramPlotRange$xMin <<- brush$xmin
        dendrogramPlotRange$xMax <<- brush$xmax
        dendrogramPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        dendrogramPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
      } else {
        dendrogramPlotRange$xMin <<- 1
        dendrogramPlotRange$xMax <<- filter$numberOfPrecursorsFiltered
        dendrogramPlotRange$xInterval <<- c(1, filter$numberOfPrecursorsFiltered)
        dendrogramPlotRange$xIntervalSize <<- filter$numberOfPrecursorsFiltered - 1
      }
    })
    obsDendrogramRangeUpdate <- observe({
      xInterval <- dendrogramPlotRange$xInterval
      
      if(!state$showHCAplotPanel)
        return()
      
      print(paste("observe dendrogramPlotRange", xInterval))
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(dataList = dataList, filter = filter$filter, clusterDataList = clusterDataList, annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, annoPresentColorsList = dataList$annoPresentColorsList, distanceMeasure = currentDistanceMeasure, nodeIndex = treeNodeFromTreeSelection, nodeIndeces = treeNodeSetFromMS2Fragment, xInterval = dendrogramPlotRange$xInterval)
      })
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("output$clusterPlotHeatmap"))
        calcClusterPlotHeatmap(dataList = dataList, filter = filter, clusterDataList = clusterDataList, xInterval = dendrogramPlotRange$xInterval)
      })
    })
    #########################################################################################
    ##' listen to heatmap plot mouse events
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
      precursorIndex <- filter$filter[[treeLeafIndex]]
      
      msg <- list()
      msg[[length(msg) + 1]] <- "Precursor: "
      msg[[length(msg) + 1]] <- dataList$precursorLabels[[precursorIndex]]
      msg[[length(msg) + 1]] <- "\n"
      
      if(hoverY > 2){
        ## lcf
        msg[[length(msg) + 1]] <- paste("log fold change [ log_2( mean(group ", filter$groups[[2]], ") / mean(group ", filter$groups[[1]], ") ) ]: ", sep = "")
        lfc <- dataList$dataFrameMeasurements[precursorIndex, dataList$lfcColumnNameFunctionFromName(filter$groups[[1]], filter$groups[[2]])]
        msg[[length(msg) + 1]] <- as.numeric(format(x = lfc, digits = 2))
      } else {
        if(hoverY > 1){
          ## group 1
          groupHere <- filter$groups[[1]]
        } else {
          ## group 2
          groupHere <- filter$groups[[2]]
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
    #########################################################################################
    ##' listen to MS2 plot mouse events
    obsMS2hover <- observeEvent(input$clusterPlotMS2_hover, {
      hoverX <- input$clusterPlotMS2_hover$x
      hoverY <- input$clusterPlotMS2_hover$y
      plotWidth  <- session$clientData$output_clusterPlotMS2_width
      plotHeight  <- session$clientData$output_clusterPlotMS2_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(is.null(fragmentsX))
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
      if(is.null(fragmentsX))
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
        if(state$showHCAplotPanel){
          fragmentMass = fragmentsX[[minimumIndex]]
          yesNoFunction <- function(precursorIndex){
            features <- dataList$featureIndeces[[precursorIndex]]
            fragmentMasses <- dataList$fragmentMasses[features]
            comprisesFragmentMass <- any(fragmentMasses == fragmentMass)
            return(comprisesFragmentMass)
          }
          #treeNodeSetFromMS2Fragment <<- getSetOfSubTreesFromRootForMass(dataList = dataList, filter = filter, clusterDataList = clusterDataList, fragmentMass = fragmentsX[[minimumIndex]])
          treeNodeSetFromMS2Fragment <<- getSetOfSubTreesFromRoot(filter = filter, clusterDataList = clusterDataList, yesNoFunction = yesNoFunction)
        }
        
        ## PCA
        if(state$showPCAplotPanel){
          fragmentMass     <- fragmentsX[[minimumIndex]]
          fragmentIndex    <- which(dataList$fragmentMasses == fragmentMass)
          precursorMatches <- dataList$featureMatrix[filter$filter, fragmentIndex] != 0
          pcaLoadingSetFromMS2Fragment <<- precursorMatches
        }
      }
      
      ## HCA
      if(state$showHCAplotPanel){
        output$clusterPlotDendrogram <- renderPlot({
          print(paste("output$clusterPlotDendrogram"))
          calcClusterPlotDendrogram(dataList = dataList, filter = filter$filter, clusterDataList = clusterDataList, annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, annoPresentColorsList = dataList$annoPresentColorsList, distanceMeasure = currentDistanceMeasure, nodeIndex = treeNodeFromTreeSelection, nodeIndeces = treeNodeSetFromMS2Fragment, xInterval = dendrogramPlotRange$xInterval)
        })
      }
      
      ## PCA
      if(state$showPCAplotPanel){
        output$pcaPlotLoadings <- renderPlot({
          print(paste("output$pcaPlotLoadings"))
          plotPCAloadings(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, pcaLoading = pcaLoadingFromPcaLoadingSelection, pcaLoadingSet = pcaLoadingSetFromMS2Fragment, xInterval = pcaLoadingsPlotRange$xInterval, yInterval = pcaLoadingsPlotRange$yInterval)
        })
      }
      
      ## update node selection
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
    })
    obsMS2dblClick <- observeEvent(input$clusterPlotMS2_dblclick, {
      brush <- input$clusterPlotMS2_brush
      
      if(is.null(fragmentsX))
        return()
      
      print(paste("observe MS2 dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ms2PlotRange$xMin <<- brush$xmin
        ms2PlotRange$xMax <<- brush$xmax
        ms2PlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        ms2PlotRange$xIntervalSize <<- brush$xmax - brush$xmin
      } else {
        ms2PlotRange$xMin <<- dataList$minimumMass
        ms2PlotRange$xMax <<- dataList$maximumMass
        ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
        ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
      }
    })
    obsMS2rangeUpdate <- observe({
      xInterval <- ms2PlotRange$xInterval
      
      if(!(state$showHCAplotPanel | state$showPCAplotPanel))
        return()
      
      print(paste("observe ms2PlotRange", paste(xInterval, collapse = ", ")))
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
    })
    #########################################################################################
    ##' listen to PCA events
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
      
      dataColumnName <- dataList$dataColumnsNameFunctionFromNames(filter$groups)[[minimumIndex]]
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
        pcaScoresPlotRange$xMin <<- brush$xmin
        pcaScoresPlotRange$xMax <<- brush$xmax
        pcaScoresPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        pcaScoresPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
        pcaScoresPlotRange$yMin <<- brush$ymin
        pcaScoresPlotRange$yMax <<- brush$ymax
        pcaScoresPlotRange$yInterval <<- c(brush$ymin, brush$ymax)
        pcaScoresPlotRange$yIntervalSize <<- brush$ymax - brush$ymin
      } else {
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
      output$pcaPlotScores <- renderPlot({
        print(paste("observe pcaScoresPlotRange output$pcaPlotScores"))
        plotPCAscores(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, xInterval = xInterval, yInterval = yInterval)
      })
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
        selectedFragmentIndex <<- NULL
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        treeNodeFromTreeSelection <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings click", sep = ""))
          paste("", sep = "")
        })
      } else {
        ## loadng selection
        pcaLoadingFromPcaLoadingSelection <<- minimumIndex
        pcaLoadingSetFromMS2Fragment <<- NULL
        precursorIndex <- filter$filter[[minimumIndex]]
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings click", sep = ""))
          paste("Info for Precursor ", dataList$precursorLabels[[precursorIndex]], sep = "")
        })
        
        selectedFragmentIndex <<- NULL
        resultObj <- getMS2spectrumOfPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        treeNodeFromTreeSelection <<- NULL
        treeNodeSetFromMS2Fragment <<- NULL
      }
      
      ## MS2 plot
      ms2PlotRange$xMin <<- dataList$minimumMass
      ms2PlotRange$xMax <<- dataList$maximumMass
      ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
      
      output$pcaPlotLoadings <- renderPlot({
        print(paste("output$pcaPlotLoadings"))
        plotPCAloadings(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, pcaLoading = pcaLoadingFromPcaLoadingSelection, pcaLoadingSet = pcaLoadingSetFromMS2Fragment, xInterval = pcaLoadingsPlotRange$xInterval, yInterval = pcaLoadingsPlotRange$yInterval)
      })
      if(state$showHCAplotPanel){
        ## update dendrogram plot
        output$clusterPlotDendrogram <- renderPlot({
          print(paste("output$clusterPlotDendrogram"))
          calcClusterPlotDendrogram(dataList = dataList, filter = filter$filter, clusterDataList = clusterDataList, annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, annoPresentColorsList = dataList$annoPresentColorsList, distanceMeasure = currentDistanceMeasure, nodeIndex = treeNodeFromTreeSelection, nodeIndeces = treeNodeSetFromMS2Fragment)
        })
      }
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
        
        precursorIndex <- filter$filter[[minimumIndex]]
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings hover", sep = ""))
          paste("Info for Precursor ", dataList$precursorLabels[[precursorIndex]], sep = "")
        })
        
        resultObj <- getMS2spectrumOfPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        fragmentsXhovered <<- resultObj$fragmentMasses
        fragmentsYhovered <<- resultObj$fragmentAbundances
      }
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, fragmentsX_02 = fragmentsXhovered, fragmentsY_02 = fragmentsYhovered, xInterval = ms2PlotRange$xInterval, selectedFragmentIndex = selectedFragmentIndex)
      })
    })
    obsPCAloadingsDblClick <- observeEvent(input$pcaPlotLoadings_dblclick, {
      brush <- input$pcaPlotLoadings_brush
      
      print(paste("observe PCAloadings dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        pcaLoadingsPlotRange$xMin <<- brush$xmin
        pcaLoadingsPlotRange$xMax <<- brush$xmax
        pcaLoadingsPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        pcaLoadingsPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
        pcaLoadingsPlotRange$yMin <<- brush$ymin
        pcaLoadingsPlotRange$yMax <<- brush$ymax
        pcaLoadingsPlotRange$yInterval <<- c(brush$ymin, brush$ymax)
        pcaLoadingsPlotRange$yIntervalSize <<- brush$ymax - brush$ymin
      } else {
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
      output$pcaPlotLoadings <- renderPlot({
        print(paste("output$pcaPlotLoadings"))
        plotPCAloadings(pcaObj = pcaDataList$pcaObj, dataList = dataList, filter = filter, pcaDimensionOne = pcaDataList$dimensionOne, pcaDimensionTwo = pcaDataList$dimensionTwo, pcaLoading = pcaLoadingFromPcaLoadingSelection, pcaLoadingSet = pcaLoadingSetFromMS2Fragment, xInterval = xInterval, yInterval = yInterval)
      })
    })
    #########################################################################################
    ##' listen to annotation events
    ##' TODO
    obsRemovePresentAnno <- observeEvent(input$removePresentAnnotation, {
      value     <- input$presentAnnotationValue
      
      drawPlots <- as.numeric(input$removePresentAnnotation)
      if(drawPlots == removePresentAnnotationValue)
        return()
      removePresentAnnotationValue <<- drawPlots
      print(paste("Observe removePresentAnnotation", drawPlots))
      
      removeAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value)
    })
    #obs <- observeEvent(input$annotationSelection, {    })
    obsToggleAddNewAnnoButton <- observeEvent(input$newAnnotationValue, {
      value <- input$newAnnotationValue
      
      print(paste("Observe newAnnotationValue", nchar(value)))
      
      if(nchar(value) > 0)
        enableActionButton(session, "submitNewAnnotation")
      else
        disableActionButton(session, "submitNewAnnotation")
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
    #########################################################################################
    ##' suspend observer
    session$onSessionEnded(function() {
      ## sidepanel
      obsTabs$suspend()
      obsFile$suspend()
      obsApplyFilters$suspend()
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
      obsPCAloadingsHover$suspend()
      obsPCAloadingsDblClick$suspend()
      obsPCAloadingsRangeUpdate$suspend()
      ## anno
      obsRemovePresentAnno$suspend()
      obsAddNewAnno$suspend()
      obsAddPresentAnno$suspend()
      ## table
      obsIgnoreValueChanged$suspend()
    })
    #########################################################################################
    #########################################################################################
    ## direct output rendering
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
      disableActionButton(session, "drawHCAplots")
      disableActionButton(session, "drawPCAplots")
      #disableActionButton(session, "downloadMatrix")
      #shinyjs::toggleState("downloadMatrix", FALSE)
      output$filteredPrecursors <- renderText({
        print(paste("init filteredPrecursors", sep = ""))
        paste("Please perform filtering", sep = "")
      })
      output$information <- renderText({
        print(paste("init information", sep = ""))
        paste("Please perform ploting.", sep = "")
      })
      return(!is.null(input$matrixFile))
    })
    output$analysisType <- reactive({
      print(paste("reactive update analysisType", state$analysisType))
      return(state$analysisType)
    })
    output$showHCAplotPanel <- reactive({
      print(paste("reactive update showHCAplotPanel", state$showHCAplotPanel))
      return(state$showHCAplotPanel)
    })
    output$showPCAplotPanel <- reactive({
      print(paste("reactive update showPCAplotPanel", state$showPCAplotPanel))
      return(state$filterValid)
    })
    output$filterValid <- reactive({
      print(paste("reactive update filterValid", state$filterValid))
      return(state$filterValid)
    })
    output$precursorSetSelected <- reactive({
      print(paste("reactive update precursorSetSelected", state$precursorSetSelected))
      return(state$precursorSetSelected)
    })
    #########################################################################################
    #########################################################################################
    ## download
    createExportMatrixName <- function(){
      fileName <- paste(gsub(" ", "_", gsub(":", ".", Sys.time())), "_selectedPrecursorMatrix.txt", sep = "")
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
      
      dataFrame <- dataList$dataFrameOriginal[c(1:3, precursorSet + 3), ]
      
      return(dataFrame)
    }
    ## download filtered
    output$downloadFilteredPrecursors <- downloadHandler(
      # This function returns a string which tells the client browser what name to use when saving the file.
      filename = function() {
        createExportMatrixName()
      },
      # This function should write data to a file given to it by the argument 'file'.
      content = function(file) {
        ## all precursors
        precursorSet <- filter$filter
        
        ## create and write data
        dataFrame <- createExportMatrix(precursorSet)
        write.table(x = dataFrame, file = file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
      }
    )
    ## download selected
    output$downloadMatrix <- downloadHandler(
      # This function returns a string which tells the client browser what name to use when saving the file.
      filename = function() {
        createExportMatrixName()
      },
      # This function should write data to a file given to it by the argument 'file'.
      content = function(file) {
        ## get selected precursors
        if(is.null(treeNodeFromTreeSelection)){
          ## all precursors
          precursorSet <- filter$filter
        } else {
          if(treeNodeFromTreeSelection < 0){
            ## node
            precursorSet <- filter$filter[-treeNodeFromTreeSelection]
          } else {
            ## selected precursors
            precursorSet <- clusterDataList$innerNodeMembers[[treeNodeFromTreeSelection]]
          }
        }
        
        ## create and write data
        dataFrame <- createExportMatrix(precursorSet)
        write.table(x = dataFrame, file = file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
      }
    )
    #########################################################################################
    #########################################################################################
    ## properties
    options(shiny.maxRequestSize=500*1024^2) 
    outputOptions(output, 'showGUI', suspendWhenHidden=FALSE)
    outputOptions(output, 'showHCAplotPanel', suspendWhenHidden=FALSE)
    outputOptions(output, 'showPCAplotPanel', suspendWhenHidden=FALSE)
    outputOptions(output, 'analysisType', suspendWhenHidden=FALSE)
    outputOptions(output, 'filterValid', suspendWhenHidden=FALSE)
    outputOptions(output, 'precursorSetSelected', suspendWhenHidden=FALSE)
  }## function(input, output, session)
)## shinyServer
