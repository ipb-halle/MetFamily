
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

#install.packages("shiny")
library(shiny)
#install.packages("shinyjs")
library(shinyjs)
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
    dataList <- NULL
    currentDistanceMeasure <- NULL
    filter <- NULL
    clusterFilter <- NULL
    pca <- NULL
    
    state <- reactiveValues(showHCAplotPanel = FALSE, showPCAplotPanel = FALSE, analysisType = NULL, filterValid = FALSE)
    submitButtonValue <- 0
    drawHCAButtonValue <- 0
    drawPCAButtonValue <- 0
    downloadMatrixButtonValue <- 0
    
    pcaScoresPlotRange   <- reactiveValues(xMin = NULL, xMax = NULL, yMin = NULL, yMax = NULL)
    pcaLoadingsPlotRange <- reactiveValues(xMin = NULL, xMax = NULL, yMin = NULL, yMax = NULL)
    ms2PlotRange         <- reactiveValues(xMin = NULL, xMax = NULL)
    dendrogramPlotRange  <- reactiveValues(xMin = NULL, xMax = NULL)
    
    selectedTreeNode <- NULL
    selectedTreeNodes <- NULL
    fragmentSelectedInTree <- FALSE
    
    fragmentsX <- NULL
    fragmentsY <- NULL
    pcaScoresX <- NULL
    pcaScoresY <- NULL
    pcaLoadingsX <- NULL
    pcaLoadingsY <- NULL
    
    pcaDimensionOne <- NULL
    pcaDimensionTwo <- NULL
    
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
      
      filter <<- NULL
      clusterFilter <<- NULL
      pca <<- NULL
      pcaScoresX <<- NULL
      pcaScoresY <<- NULL
      pcaLoadingsX <<- NULL
      pcaLoadingsY <<- NULL
      
      state$analysisType <<- NULL
      state$filterValid <<- FALSE
      state$showHCAplotPanel <<- FALSE
      state$showPCAplotPanel <<- FALSE
      
      selectedTreeNode <<- NULL
      selectedTreeNodes <<- NULL
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      fragmentSelectedInTree <<- FALSE
      
      ms2PlotRange$xMin <<- NULL
      ms2PlotRange$xMax <<- NULL
      dendrogramPlotRange$xMin <<- NULL
      dendrogramPlotRange$xMax <<- NULL
      pcaScoresPlotRange$xMin <<- NULL
      pcaScoresPlotRange$xMax <<- NULL
      pcaScoresPlotRange$yMin <<- NULL
      pcaScoresPlotRange$yMax <<- NULL
      pcaLoadingsPlotRange$xMin <<- NULL
      pcaLoadingsPlotRange$xMax <<- NULL
      pcaLoadingsPlotRange$yMin <<- NULL
      pcaLoadingsPlotRange$yMax <<- NULL
      
      pcaDimensionOne <<- NULL
      pcaDimensionTwo <<- NULL
      
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
      updateTextInput(session = session, inputId = "filter_ms2_masses", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_ppm", value = "20")
      
      #doPerformFiltering()
      output$information <- renderText({
        print(paste("init output$information", sep = ""))
        paste("Number of filtered precursors: ", dataList$numberOfPrecursorsFiltered, sep = "")
      })
      
      doClearPlots()
    }
    
    #########################################################################################
    ##' perform filtering
    doPerformFiltering <- function(){
      #################################################
      ## process inputs
      
      ## get inputs
      groupOne <- input$groupOne
      groupTwo <- input$groupTwo
      groups <- input$groups
      filter_average <- input$filter_average
      filter_lfc <- input$filter_lfc
      filter_ms2_masses <- input$filter_ms2_masses
      filter_ms2_masses2 <- input$filter_ms2_masses2
      filter_ms2_masses3 <- input$filter_ms2_masses3
      filter_ms2_ppm <- input$filter_ms2_ppm
      #analysisType <- input$analysisType
      print(paste("Observe applyFilters1", "gs", paste(groups, collapse = "-"), "g1", groupOne, "g2", groupTwo, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm))
      
      #if(analysisType == "PCA")
      #  filter_lfc <- NULL
      
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
      
      if(is.null(filter_ms2_masses) | nchar(filter_ms2_masses) == 0)
        filter_ms2_masses <- NULL
      if(is.null(filter_ms2_masses2) | nchar(filter_ms2_masses2) == 0)
        filter_ms2_masses2 <- NULL
      if(is.null(filter_ms2_masses3) | nchar(filter_ms2_masses3) == 0)
        filter_ms2_masses3 <- NULL
      
      if(is.null(filter_ms2_ppm) | nchar(filter_ms2_ppm) == 0)
        filter_ms2_ppm <- NULL
      else
        filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
      
      print(paste("Observe applyFilters2", "gs", paste(groups, collapse = "-"), "g1", groupOne, "g2", groupTwo, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm))
      print(paste("Observe applyFilters2", "gs", is.null(groups), "g1", is.null(groupOne), "g2", is.null(groupTwo), "a", is.null(filter_average), "lfc", is.null(filter_lfc), "ms2 1", is.null(filter_ms2_masses), "ms2 2", is.null(filter_ms2_masses2), "ms2 3", is.null(filter_ms2_masses3), "ppm", is.null(filter_ms2_ppm)))
      
      ## check for errors in inputs amd process ms2
      error <- FALSE
      if(!is.null(filter_average))
        error <- error | is.na(filter_average)
      if(!is.null(filter_lfc))
        error <- error | is.na(filter_lfc)
      if(!is.null(filter_ms2_masses)){
        ms2Masses <- strsplit(x = filter_ms2_masses, split = "[,; ]+")[[1]]
        filter_ms2_masses <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses))
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
      error <- error | ((!is.null(filter_ms2_masses) | !is.null(filter_ms2_masses2) | !is.null(filter_ms2_masses3)) & is.null(filter_ms2_ppm))
      
      print(paste("Observe applyFilters3", "e", error, "gs", paste(groups, collapse = "-"), "g1", groupOne, "g2", groupTwo, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm))
      if(error){
        output$filteredPrecursors <- renderText({
          print(paste("update output$filteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
        return()
      }
      
      filterList_ms2_masses <- list()
      if(!is.null(filter_ms2_masses))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses
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
      
      filter <<- filterData(
        dataList = dataList, 
        #groupOne = groupOne, groupTwo = groupTwo, 
        groups = groupSet, filter_average = filter_average, filter_lfc = filter_lfc, 
        filterList_ms2_masses = filterList_ms2_masses, filter_ms2_ppm = filter_ms2_ppm, 
        progress = FALSE
      )
      numberOfPrecursorsFiltered <- length(filter$filter)
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
        #if(analysisType == "HCA")
          enableActionButton(session, "drawHCAplots")
        #if(analysisType == "PCA")
          enableActionButton(session, "drawPCAplots")
        
        print(paste("#+#1", state$filterValid))
        state$filterValid <<- TRUE
        print(paste("#+#2", state$filterValid))
      } else {
        ## filter invalid
        filter <<- NULL
        state$filterValid <<- FALSE
        
        #if(analysisType == "HCA")
          disableActionButton(session, "drawHCAplots")
        #if(analysisType == "PCA")
          disableActionButton(session, "drawPCAplots")
        
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
        NULL
      })
    }
    #########################################################################################
    ##' Calculate plots
    doCalculatePlots <- function(){
      print(paste("doCalcClusterPlot", is.null(dataList), is.null(clusterFilter)))
      #if(is.null(dataList) | is.null(clusterFilter)){
      #  return()
      #}
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(filterList = clusterFilter, distanceMeasure = currentDistanceMeasure, nodeIndex = selectedTreeNode, nodeIndeces = selectedTreeNodes)
      })
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("output$clusterPlotHeatmap"))
        calcClusterPlotHeatmap(dataList = dataList, filter = filter, filterList = clusterFilter)
      })
      output$clusterPlotHeatmapLegend <- renderPlot({
        print(paste("output$clusterPlotHeatmapLegend"))
        calcClusterPlotHeatmapLegend(dataList = dataList)
      })
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY)
      })
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
    #  if(!is.null(dataList) & !is.null(clusterFilter))
    #    doCalculatePlots()
    #})
    
    #########################################################################################
    ##' listen to filtab selections
    obsFile <- observe({
      tabId <- input$runTabs
      print(paste("Observe tab selection", tabId))
      if(tabId == "HCA")
        state$analysisType <<- "HCA"
      if(tabId == "PCA")
        state$analysisType <<- "PCA"
    })
    
    #########################################################################################
    ##' listen to file input
    obsFile <- observe({
      file <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      print(paste("Observe file for data", fileName))
      if(!is.null(file)){
        output$fileInfo <- renderText({paste("Processing file ", fileName, "...", sep = "")})
        getClusterData(file)
        output$fileInfo <- renderText({fileName})
      }
    })
    
    #########################################################################################
    ##' listen to heatmap plot mouse events
    obsHeatmaphover <- observe({
      hoverX <- input$clusterPlotHeatmap_hover$x
      hoverY <- input$clusterPlotHeatmap_hover$y
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(hoverX < 0.5 | hoverX > (clusterFilter$numberOfPrecursorsFiltered + 0.5))
        return()
      if(hoverY < 0 | hoverY > 3)
        return()
      
      print(paste("Observe heatmap hover", hoverX, hoverY))
      
      treeLeafIndex2 <- as.numeric(format(x = hoverX, digits = 0))
      treeLeafIndex  <- clusterFilter$cluster$order[[treeLeafIndex2]]
      precursorIndex <- filter$filter[[treeLeafIndex]]
      
      msg <- list()
      msg[[length(msg) + 1]] <- "Precursor: "
      msg[[length(msg) + 1]] <- dataList$precursorLabels[[precursorIndex]]
      msg[[length(msg) + 1]] <- "\n"
      
      if(hoverY > 2){
        ## lcf
        msg[[length(msg) + 1]] <- paste("log fold change [ log_2( mean(group ", filter$groups[[2]], ") / mean(group ", filter$groups[[1]], ") ) ]: ", sep = "")
        lfc <- dataList$dataFrameMeasurements[precursorIndex, dataList$lfcColumnNameFunction(filter$groups[[1]], filter$groups[[2]])]
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
        val <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunction(groupHere)]
        msg[[length(msg) + 1]] <- as.numeric(format(x = val, digits = 2))
        msg[[length(msg) + 1]] <- " = mean("
        vals <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataColumnsNameFunction(groupHere)]
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
    obsMS2click <- observe({
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
      if(is.null(ms2PlotRange$xMin)){
        massIntervalSize <- dataList$maximumMass - dataList$minimumMass
      } else {
        massIntervalSize <- ms2PlotRange$xMax - ms2PlotRange$xMin
      }
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = fragmentsX, poiCoordinatesY = fragmentsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = massIntervalSize, plotRangeY = 1
      )
      print(paste("Observe MS2 click", minimumIndex))
      if(is.null(minimumIndex))
        return()
      
      #analysisType <- input$analysisType
      
      if(state$showHCAplotPanel){
        ## HCA
        if(is.null(minimumIndex)){
          if(fragmentSelectedInTree){
            fragmentSelectedInTree <<- FALSE
            selectedTreeNodes <<- NULL
          }
        } else {
          ## fetch subroots of subtrees comprising the selected fragment
          selectedTreeNodes <<- getSetOfSubTreesFromRoot(dataList = dataList, clusterFilter = filter, clusterCluster = clusterFilter, fragmentMass = fragmentsX[[minimumIndex]])
          fragmentSelectedInTree <<- TRUE
        }
        
        if(is.null(dendrogramPlotRange$xMin)){
          xInterval <- NULL
        } else {
          xInterval = c(dendrogramPlotRange$xMin, dendrogramPlotRange$xMax)
        }
        output$clusterPlotDendrogram <- renderPlot({
          print(paste("output$clusterPlotDendrogram"))
          calcClusterPlotDendrogram(filterList = clusterFilter, distanceMeasure = currentDistanceMeasure, nodeIndex = selectedTreeNode, nodeIndeces = selectedTreeNodes, xInterval = xInterval)
        })
      }
      if(state$showPCAplotPanel){
        ## TODO PCA
        fragmentMass     <- fragmentsX[[minimumIndex]]
        fragmentIndex    <- which(dataList$fragmentMasses == fragmentMass)
        precursorMatches <- dataList$featureMatrix[filter$filter, fragmentIndex] != 0
        
        ## TODO
      }
    })
    obsMS2hover <- observe({
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
      if(is.null(ms2PlotRange$xMin)){
        massIntervalSize <- dataList$maximumMass - dataList$minimumMass
      } else {
        massIntervalSize <- ms2PlotRange$xMax - ms2PlotRange$xMin
      }
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = fragmentsX, poiCoordinatesY = fragmentsY, 
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = massIntervalSize, plotRangeY = 1
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
    obsMS2dblClick <- observeEvent(input$clusterPlotMS2_dblclick, {
      brush <- input$clusterPlotMS2_brush
      
      print(paste("observe MS2 dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ms2PlotRange$xMin <- brush$xmin
        ms2PlotRange$xMax <- brush$xmax
      } else {
        ms2PlotRange$xMin <- NULL
        ms2PlotRange$xMax <- NULL
      }
    })
    obsMS2rangeUpdate <- observe({
      xMin <- ms2PlotRange$xMin
      xMax <- ms2PlotRange$xMax
      
      if(is.null(dataList))
        return()
      
      if(is.null(xMin)){
        xInterval <- NULL
      } else {
        xInterval = c(xMin, xMax)
      }
      print(paste("observe ms2PlotRange", xMin, xMax))
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY, xInterval = xInterval)
      })
    })
    #########################################################################################
    ##' listen to dendrogram plot mouse events
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
      if(is.null(dendrogramPlotRange$xMin)){
        xIntervalSize <- clusterFilter$numberOfPrecursorsFiltered
      } else {
        xIntervalSize = dendrogramPlotRange$xMax - dendrogramPlotRange$xMin
      }
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = clusterFilter$poiCoordinatesX, poiCoordinatesY = clusterFilter$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = xIntervalSize, plotRangeY = 1
      )
      print(paste("Observe dendrogram click", minimumIndex))
      if(is.null(minimumIndex))
        return()
      
      minimumLabel <- clusterFilter$poiLabels[[minimumIndex]]
      minimumText <- clusterFilter$poiText[[minimumIndex]]
      
      #################################################
      ## fetch ms2 spectrum
      selectedTreeNodes <<- NULL
      selectedTreeNode <<- minimumLabel
      resultObj <- getMS2spectrum(dataList = dataList, filterList = clusterFilter, label = minimumLabel)
      fragmentsX <<- resultObj$fragmentsX
      fragmentsY <<- resultObj$fragmentsY
      
      #################################################
      ## output as message and plots
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
      
      ## cluster dendrogram
      if(is.null(dendrogramPlotRange$xMin)){
        xInterval <- NULL
      } else {
        xInterval = c(dendrogramPlotRange$xMin, dendrogramPlotRange$xMax)
      }
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("update output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(filterList = clusterFilter, distanceMeasure = currentDistanceMeasure, nodeIndex = minimumLabel, nodeIndeces = selectedTreeNodes, xInterval = xInterval)
      })
      
      ## MS2 plot
      ms2PlotRange$xMin <<- NULL
      ms2PlotRange$xMax <<- NULL
      output$clusterPlotMS2 <- renderPlot({
        print(paste("update output$clusterPlotMS2"))
        calcClusterPlotMS2(
          dataList = dataList, 
          fragmentsX = resultObj$fragmentsX, 
          fragmentsY = resultObj$fragmentsY
        )
      })
      
      ## enable button
      #enableActionButton(session, "downloadMatrix")
      #shinyjs::toggleState("downloadMatrix", TRUE)
    })
    obsDendrogramdblClick <- observeEvent(input$clusterPlotDendrogram_dblclick, {
      brush <- input$clusterPlotDendrogram_brush
      
      print(paste("observe dendrogram dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        dendrogramPlotRange$xMin <- brush$xmin
        dendrogramPlotRange$xMax <- brush$xmax
      } else {
        dendrogramPlotRange$xMin <- NULL
        dendrogramPlotRange$xMax <- NULL
      }
    })
    obsDendrogramRangeUpdate <- observe({
      xMin <- dendrogramPlotRange$xMin
      xMax <- dendrogramPlotRange$xMax
      
      if(is.null(dataList))
        return()
      
      if(is.null(xMin)){
        xInterval <- NULL
      } else {
        xInterval = c(xMin, xMax)
      }
      print(paste("observe dendrogramPlotRange", xMin, xMax))
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(filterList = clusterFilter, distanceMeasure = currentDistanceMeasure, nodeIndex = selectedTreeNode, nodeIndeces = selectedTreeNodes, xInterval = xInterval)
      })
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("output$clusterPlotHeatmap"))
        calcClusterPlotHeatmap(dataList = dataList, filter = filter, filterList = clusterFilter, xInterval = xInterval)
      })
    })
    #########################################################################################
    ##' listen to PCA events
    obsPCAscoresDblClick <- observeEvent(input$pcaPlotScores_dblclick, {
      brush <- input$pcaPlotScores_brush
      
      print(paste("observe PCAscores dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        pcaScoresPlotRange$xMin <- brush$xmin
        pcaScoresPlotRange$xMax <- brush$xmax
        pcaScoresPlotRange$yMin <- brush$ymin
        pcaScoresPlotRange$yMax <- brush$ymax
      } else {
        pcaScoresPlotRange$xMin <- NULL
        pcaScoresPlotRange$xMax <- NULL
        pcaScoresPlotRange$yMin <- NULL
        pcaScoresPlotRange$yMax <- NULL
      }
    })
    obsPCAscoresRangeUpdate <- observe({
      xMin <- pcaScoresPlotRange$xMin
      xMax <- pcaScoresPlotRange$xMax
      yMin <- pcaScoresPlotRange$yMin
      yMax <- pcaScoresPlotRange$yMax
      
      if(is.null(dataList))
        return()
      
      if(is.null(xMin)){
        xInterval <- NULL
        yInterval <- NULL
      } else {
        xInterval = c(xMin, xMax)
        yInterval = c(yMin, yMax)
      }
      print(paste("observe pcaScoresPlotRange", xMin, xMax, yMin, yMax))
      
      ## plot PCA
      output$pcaPlotScores <- renderPlot({
        print(paste("output$pcaPlotScores"))
        plotPCAscores(pca, dataList, filter, pcaDimensionOne, pcaDimensionTwo, xInterval, yInterval)
      })
    })
    obsPCAloadingsDblClick <- observeEvent(input$pcaPlotLoadings_dblclick, {
      brush <- input$pcaPlotLoadings_brush
      
      print(paste("observe PCAloadings dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        pcaLoadingsPlotRange$xMin <- brush$xmin
        pcaLoadingsPlotRange$xMax <- brush$xmax
        pcaLoadingsPlotRange$yMin <- brush$ymin
        pcaLoadingsPlotRange$yMax <- brush$ymax
      } else {
        pcaLoadingsPlotRange$xMin <- NULL
        pcaLoadingsPlotRange$xMax <- NULL
        pcaLoadingsPlotRange$yMin <- NULL
        pcaLoadingsPlotRange$yMax <- NULL
      }
    })
    obsPCAloadingsRangeUpdate <- observe({
      xMin <- pcaLoadingsPlotRange$xMin
      xMax <- pcaLoadingsPlotRange$xMax
      yMin <- pcaLoadingsPlotRange$yMin
      yMax <- pcaLoadingsPlotRange$yMax
      
      if(is.null(dataList))
        return()
      
      if(is.null(xMin)){
        xInterval <- NULL
        yInterval <- NULL
      } else {
        xInterval = c(xMin, xMax)
        yInterval = c(yMin, yMax)
      }
      print(paste("observe pcaLoadingsPlotRange", xMin, xMax, yMin, yMax))
      
      ## plot PCA
      output$pcaPlotLoadings <- renderPlot({
        print(paste("output$pcaPlotLoadings"))
        plotPCAloadings(pca, dataList, filter, pcaDimensionOne, pcaDimensionTwo, xInterval, yInterval)
      })
    })
    obsPCAplotScoresHover <- observe({
      hoverX <- input$pcaPlotScores_hover$x
      hoverY <- input$pcaPlotScores_hover$y
      plotWidth  <- session$clientData$output_pcaPlotScores_width
      plotHeight <- session$clientData$output_pcaPlotScores_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA scores hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      if(is.null(ms2PlotRange$xMin)){
        xIntervalSize <- max(pcaScoresX) - min(pcaScoresX)
        yIntervalSize <- max(pcaScoresY) - min(pcaScoresY)
      } else {
        xIntervalSize <- pcaScoresPlotRange$xMax - pcaScoresPlotRange$xMin
        yIntervalSize <- pcaScoresPlotRange$yMax - pcaScoresPlotRange$yMin
      }
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaScoresX, poiCoordinatesY = pcaScoresY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = xIntervalSize, plotRangeY = yIntervalSize
      )
      if(is.null(minimumIndex))
        return()
      
      print(paste("Observe PCA scores hover", hoverX, hoverY, minimumIndex))
      
      output$information <- renderText({
        print(paste("update output$information PCA scores hover", sep = ""))
        paste("Info for Replicate ", dataList$columnGroupOrgLabels[!is.na(match(x = dataList$columnGroupLabels, table = filter$groups))][[minimumIndex]], sep = "")
      })
    })
    obsPCAplotLoadingsHover <- observe({
      hoverX <- input$pcaPlotLoadings_hover$x
      hoverY <- input$pcaPlotLoadings_hover$y
      plotWidth  <- session$clientData$output_pcaPlotLoadings_width
      plotHeight <- session$clientData$output_pcaPlotLoadings_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA Loadings hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      if(is.null(ms2PlotRange$xMin)){
        xIntervalSize <- max(pcaLoadingsX) - min(pcaLoadingsX)
        yIntervalSize <- max(pcaLoadingsY) - min(pcaLoadingsY)
      } else {
        xIntervalSize <- pcaLoadingsPlotRange$xMax - pcaLoadingsPlotRange$xMin
        yIntervalSize <- pcaLoadingsPlotRange$yMax - pcaLoadingsPlotRange$yMin
      }
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaLoadingsX, poiCoordinatesY = pcaLoadingsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = xIntervalSize, plotRangeY = yIntervalSize
      )
      if(is.null(minimumIndex))
        return()
      
      print(paste("Observe PCA Loadings hover", hoverX, hoverY, minimumIndex))
      
      precursorIndex <- filter$filter[[minimumIndex]]
      
      output$information <- renderText({
        print(paste("update output$information PCA Loadings hover", sep = ""))
        paste("Info for Precursor ", dataList$precursorLabels[[precursorIndex]], sep = "")
      })
      
      resultObj <- getMS2spectrumOfPrecursor(dataList, precursorIndex)
      fragmentsX <<- resultObj$fragmentMasses
      fragmentsY <<- resultObj$fragmentAbundances
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY)
      })
    })
    
    #########################################################################################
    ##' listen to submit button events
    obsApplyFilters <- observe({
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
    obsDrawHCA <- observe({
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
      selectedTreeNode <<- NULL
      selectedTreeNodes <<- NULL
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      fragmentSelectedInTree <<- FALSE
      ms2PlotRange$xMin <<- NULL
      ms2PlotRange$xMax <<- NULL
      dendrogramPlotRange$xMin <<- NULL
      dendrogramPlotRange$xMax <<- NULL
      
      if(!is.null(filter)){
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
        withProgress(message = 'Calculating distances...', value = 0, {
          distanceMatrix <- calculateDistanceMatrix(dataList, filter, distance, progress = TRUE)
        })
        ## compute cluster
        print(paste("Calculating cluster", sep = ""))
        withProgress(message = 'Calculating cluster...', value = 0, {
          clusterFilter <<- calculateCluster(progress = TRUE,
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
        
        ## disable download button
        #disableActionButton(session, "downloadMatrix")
        #shinyjs::toggleState("downloadMatrix", FALSE)
        
        ## compute plots
        print(paste("Calculating plots", sep = ""))
        currentDistanceMeasure <<- distance
        state$showHCAplotPanel <<- TRUE
        #showPCAplotPanel$show <<- FALSE
        doCalculatePlots()
      } else {
        ##########################
        ## do not draw
        print(paste("Observe do not draw HCA plots", sep = ""))
        
        ## update info
        output$information <- renderText({
          print(paste("update do not draw HCA plots", sep = ""))
          paste("The given criteria do not result in a valid set of precursors. Please check the set of filters.", sep = "")
        })
        
        ## clear plots
        currentDistanceMeasure <<- NULL
        doClearPlots()
      }
    })
    
    #########################################################################################
    ##' listen to draw PCA button events
    obsDrawPCA <- observe({
      drawPlots <- as.numeric(input$drawPCAplots)
      
      print(paste("Observe draw PCA plots", drawPlots))
      
      #################################################
      ## check if button was hit
      if(drawPlots == drawPCAButtonValue)
        return()
      drawPCAButtonValue <<- drawPlots
      
      unitVariance <- input$pcaUnitVariance
      logTransform <- input$pcaLogTransform
      pcaDimensionOne <<- as.numeric(input$pcaDimensionOne)
      pcaDimensionTwo <<- as.numeric(input$pcaDimensionTwo)
      print(paste("Observe draw PCA plots", "UV", unitVariance, "LT", logTransform, "D1", pcaDimensionOne, "D2", pcaDimensionTwo))
      
      state$showPCAplotPanel <<- TRUE
      #showHCAplotPanel$show <<- FALSE
      
      ## calc PCA
      pcaScoresX <<- NULL
      pcaScoresY <<- NULL
      pcaLoadingsX <<- NULL
      pcaLoadingsY <<- NULL
      pca <<- calcPCA(dataList, filter, unitVariance, logTransform)
      
      pcaScoresX <<- pca$ind$coord[, pcaDimensionOne]
      pcaScoresY <<- pca$ind$coord[, pcaDimensionTwo]
      pcaLoadingsX <<- pca$var$coord[, pcaDimensionOne]
      pcaLoadingsY <<- pca$var$coord[, pcaDimensionTwo]
      
      ## plot PCA
      output$pcaPlotScores <- renderPlot({
        print(paste("output$pcaPlotScores"))
        plotPCAscores(pca, dataList, filter, pcaDimensionOne, pcaDimensionTwo)
      })
      output$pcaPlotLoadings <- renderPlot({
        print(paste("output$pcaPlotLoadings"))
        plotPCAloadings(pca, dataList, filter, pcaDimensionOne, pcaDimensionTwo)
      })
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = dataList, fragmentsX = fragmentsX, fragmentsY = fragmentsY)
      })
    })
    
    #########################################################################################
    ##' suspend observer
    session$onSessionEnded(function() {
      #obsClusterPlotResize$suspend()
      obsFile$suspend()
      obsMS2hover$suspend()
      obsHeatmaphover$suspend()
      obsMS2dblClick$suspend()
      obsMS2rangeUpdate$suspend()
      obsDendrogramClick$suspend()
      obsPCAplotScoresHover$suspend()
      obsPCAplotLoadingsHover$suspend()
      obsApplyFilters$suspend()
      obsDrawHCA$suspend()
      obsDrawPCA$suspend()
      obsPCAloadingsDblClick$suspend()
      obsPCAloadingsRangeUpdate$suspend()
      obsPCAscoresDblClick$suspend()
      obsPCAscoresRangeUpdate$suspend()
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
      return(state$showPCAplotPanel)
    })
    output$filterValid <- reactive({
      print(paste("reactive update filterValid", state$filterValid))
      return(state$filterValid)
    })
    
    #########################################################################################
    #########################################################################################
    ## download
    createExportMatrixName <- function(){
      fileName <- paste(gsub(" ", "_", gsub(":", ".", Sys.time())), "_selectedPrecursorMatrix.txt", sep = "")
      return(fileName)
    }
    createExportMatrix <- function(precursorSet){
      subMatrix <- matrix(data = dataList$featureMatrix[precursorSet, ], nrow = length(precursorSet))
      
      fragmentPresentinColumns <- apply(X = subMatrix, MARGIN = 2, FUN = sum) > 0
      subMatrix <- subMatrix[, fragmentPresentinColumns]
      
      dataFrame <- as.data.frame(subMatrix)
      
      if(length(precursorSet) > 1){
        rownames(dataFrame) <- dataList$precursorLabels[precursorSet]
        colnames(dataFrame) <- dataList$fragmentMasses[fragmentPresentinColumns]
        
        print(dataList$precursorLabels[precursorSet])
      } else {
        ## do not annotate or else error
      }
      
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
        if(is.null(selectedTreeNode)){
          ## all precursors
          precursorSet <- filter$filter
        } else {
          if(selectedTreeNode < 0){
            ## node
            precursorSet <- filter$filter[-selectedTreeNode]
          } else {
            ## selected precursors
            precursorSet <- clusterFilter$innerNodeMembers[[selectedTreeNode]]
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
    
    #observe({
    #  shinyjs::toggleState("filter_lfc", input$analysisType == "HCA")
    #})
  }## function(input, output, session)
)## shinyServer
