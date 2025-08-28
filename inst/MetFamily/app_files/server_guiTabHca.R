## HCA constants
minimumNumberOfPrecursorsForHca <- 6
maximumNumberOfPrecursorsForHca <- 5000
minimumNumberOfPrecursorsForDendrogramStatistics <- 2

minimumheatmapHeightPerRow <- 11
maximumheatmapHeightPerRow <- 25
minimumheatmapHeightRowCount <- 10
maximumheatmapHeightRowCount <- 3

currentDistanceMatrixObj <- NULL
clusterDataList <- NULL
columnsOfInterestForHeatmap <- NULL

dendrogramUserCoordinateRangeY <- NULL
dendrogramPlotRange  <- reactiveValues(
  xMin = NULL, 
  xMax = NULL, 
  xInterval = NULL, 
  xIntervalSize = NULL
)

## plotly dendrogram
curveNumberToCurveName <- NULL

## classifier table for dendrogram node
putativeAnnotationsTableFromAnalysisInputFieldIdCounter <- 0
putativeAnnotationsTableFromAnalysisCurrentlyShown <- NULL

## statistics for dendrogram node or pca loadings brush
dendrogramFragmentStatistics <- FALSE


# Extended task for HCA computation
hcaComputationTask <- ExtendedTask$new(function(dataList, filterHca, distanceMeasure,
                                                clusterMethod, minimumProportionOfLeafs,
                                                minimumProportionToShowFragment) {
  promise_result <- promises::future_promise({
    # global variables needed by HCA functions
    minimumProportionOfLeafs <<- minimumProportionOfLeafs
    minimumProportionToShowFragment <<- minimumProportionToShowFragment
    
    result <- tryCatch({
      ## compute distance matrix
      distanceMatrixObj <- calculateDistanceMatrix(
        dataList = dataList, 
        filter = filterHca$filter, 
        distanceMeasure = distanceMeasure, 
        progress = FALSE
      )
      
      ## compute cluster  
      clusterDataListObj <- calculateCluster(
        progress = FALSE, 
        dataList = dataList, 
        filterObj = filterHca, 
        distanceMatrix = distanceMatrixObj$distanceMatrix, 
        method = clusterMethod, 
        distanceMeasure = distanceMeasure
      )
      
      list(
        distanceMatrixObj = distanceMatrixObj,
        clusterDataList = clusterDataListObj
      )
      
    }, error = function(e) {
      stop(paste("HCA computation error:", e$message))
    })
    
    return(result)
    
  }, seed = TRUE)
  
  return(promise_result)
})

state_tabHca <- reactiveValues(
  ## plot dimensions
  #dendrogramHeatmapHeight = 1,#heatmapHeightPerRow * 3,
  heatmapHeight = 1,#heatmapHeightPerRow * 3,## plotly: remove
  ## plot controls
  showClusterLabels = TRUE,
  heatmapContent = "Log-fold-change",
  heatmapOrdering = "Specified order",
  hcaPrecursorLabels = "m/z / RT",
  ## plot annotations: $setOfAnnotations, $setOfColors
  annotationsHca = NULL,
  ## annotation legend height
  annotationLegendHeightHca = -1
)
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset tabHca state")
  state_tabHca$heatmapHeight <<- 1#heatmapHeightPerRow * 3,## plotly: remove
  ## plot controls
  state_tabHca$showClusterLabels <<- TRUE
  state_tabHca$heatmapContent <<- "Log-fold-change"
  state_tabHca$heatmapOrdering <<- "Specified order"
  state_tabHca$hcaPrecursorLabels <<- "m/z / RT"
  ## plot annotations: $setOfAnnotations, $setOfColors
  state_tabHca$annotationsHca <<- NULL
  ## annotation legend height
  state_tabHca$annotationLegendHeightHca <<- -1
  
  ## reset variables
  clusterDataList <<- NULL
  columnsOfInterestForHeatmap <<- NULL
  dendrogramFragmentStatistics <<- FALSE
  
  ## reset plot range
  # dendrogramPlotRangeY <<- NULL # never used
  dendrogramPlotRange$xMin <<- NULL
  dendrogramPlotRange$xMax <<- NULL
  dendrogramPlotRange$xInterval <<- NULL
  dendrogramPlotRange$xIntervalSize <<- NULL
  
  ## plotly dendrogram
  curveNumberToCurveName <<- NULL
  
  ## classifier table for dendrogram node
  putativeAnnotationsTableFromAnalysisCurrentlyShown <<- NULL
})


drawDendrogramPlot <- function(consoleInfo = NULL, withHeatmap = FALSE){
  output$plotDendrogram <- renderPlot({
    #output$plotDendrogram <- renderPlotly({
    print(paste("### den ###", consoleInfo))
    drawDendrogramPlotImpl()
  })
  
  if(!withHeatmap)
    return()
  
  ## plotly: remove
  output$plotHeatmap <- renderPlot({
    print(paste("### hea ### update range output$plotHeatmap"))
    drawHeatmapPlotImpl()
  })
}
drawDendrogramPlotImpl_forPlotly <- function(){
  heatmapContent <- state_tabHca$heatmapContent
  
  ## heatmap
  numberOfGroups <- -1
  switch(heatmapContent,
         "Log-fold-change"={## log-fold-change
           numberOfGroups <- 3
         },
         "Abundance by group"={## groups
           numberOfGroups <- length(dataList$sampleClasses)
         },
         "Abundance by sample"={## samples
           numberOfGroups <- length(dataList$dataColumnsNameFunctionFromGroupNames(sampleClasses = sampleClasses, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame)))
         },
         {## unknown state
           stop(paste("Unknown heatmapContent value", heatmapContent))
         }
  )## end switch
  
  ## heigth per row
  heatmapHeightPerRow <- -1
  if(numberOfGroups <= maximumheatmapHeightRowCount){
    heatmapHeightPerRow <- maximumheatmapHeightPerRow
  } else if(numberOfGroups >= minimumheatmapHeightRowCount){
    heatmapHeightPerRow <- minimumheatmapHeightPerRow
  } else {
    heatmapHeightPerRow <- 
      minimumheatmapHeightPerRow + 
      (numberOfGroups    - maximumheatmapHeightRowCount) / 
      (minimumheatmapHeightRowCount - maximumheatmapHeightRowCount) * 
      (maximumheatmapHeightPerRow - minimumheatmapHeightPerRow)
  }
  heatmapHeight <- heatmapHeightPerRow * numberOfGroups
  dendrogramHeatmapHeight <- 500 + heatmapHeight
  heatmapProportion <- heatmapHeight / dendrogramHeatmapHeight
  #state_tabHca$dendrogramHeatmapHeight <<- dendrogramHeatmapHeight
  
  resultObj <- calcPlotDendrogram(
    dataList = dataList, 
    filterObj = filterHca,#filter, 
    clusterDataList = clusterDataList, 
    #annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, 
    #annoPresentColorsList = dataList$annoPresentColorsList, 
    distanceMeasure = currentDistanceMatrixObj$distanceMeasure, 
    showClusterLabels = state_tabHca$showClusterLabels, 
    hcaPrecursorLabels = state_tabHca$hcaPrecursorLabels, 
    selectionFragmentTreeNodeSet = selectionFragmentTreeNodeSet,
    selectionAnalysisTreeNodeSet = selectionAnalysisTreeNodeSet,
    selectionSearchTreeNodeSet = selectionSearchTreeNodeSet,
    selectedSelection = state_selections$selectedSelection,
    heatmapContent = heatmapContent,
    heatmapOrdering = heatmapOrdering,
    heatmapProportion = heatmapProportion
    #xInterval = dendrogramPlotRange$xInterval
  )
  curveNumberToCurveName <<- resultObj$curveNumberToCurveName
  plotlyPlot             <- resultObj$plotlyPlot
  #columnsOfInterest      <- resultObj$columnsOfInterest
  
  dendrogramUserCoordinateRange <- par("usr")
  dendrogramUserCoordinateRangeY <<- dendrogramUserCoordinateRange[[4]] - dendrogramUserCoordinateRange[[3]]
  
  ## present annotations
  resultObjAnno <- getPrecursorColors(dataList = dataList, precursorSet = filterHca$filter)
  
  uniqueIndeces     <- which(!duplicated(resultObjAnno$setOfAnnotations))
  uniqueAnnotations <- resultObjAnno$setOfAnnotations[uniqueIndeces]
  uniqueColors      <- resultObjAnno$setOfColors[uniqueIndeces]
  
  state_tabHca$annotationsHca <<- list(
    setOfAnnotations = uniqueAnnotations,
    setOfColors      = uniqueColors
  )
  state_tabHca$annotationLegendHeightHca <<- annoLegendEntryHeight * (length(state_tabHca$annotationsHca$setOfAnnotations) + 1)
  
  return(plotlyPlot)
}
drawHeatmapLegend <- function(consoleInfo = NULL){
  output$plotHeatmapLegend <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawHeatmapLegendImpl()
  })
}
drawAnnotationLegendHCA <- function(consoleInfo = NULL){
  output$plotAnnoLegendHCA <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawAnnotationLegendHCAimpl()
    
  })
}
drawDendrogramLegend <- function(consoleInfo = NULL){
  output$calcPlotDendrogramLegend <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawDendrogramLegendImpl()
  })
}

resetHcaPlotRange <- function(){
  dendrogramPlotRange$xMin <<- 1
  dendrogramPlotRange$xMax <<- filterHca$numberOfPrecursorsFiltered
  dendrogramPlotRange$xInterval <<- c(1, filterHca$numberOfPrecursorsFiltered)
  dendrogramPlotRange$xIntervalSize <<- filterHca$numberOfPrecursorsFiltered - 1
}

obsDrawHCA <- observeEvent(input$drawHCAplots, {
  session$sendCustomMessage("disableButton", "drawHCAplots")
  #################################################
  ## get input
  drawPlots <- as.numeric(input$drawHCAplots)
  
  print(paste("Observe draw HCA plots", drawPlots))
  
  ## check if button was hit
  #if(drawPlots == drawHCAButtonValue)
  #  return()
  #drawHCAButtonValue <<- drawPlots
  
  distanceMeasure <- input$hcaDistanceFunction
  #clusterMethod <- input$hcaClusterMethod
  clusterMethod <- "ward.D"
  print(paste("Observe draw HCA plots", "D", distanceMeasure, "M", clusterMethod))
  
  ##########################
  ## Invoke async HCA computation (hcaComputationalTask)
  
  tryCatch({
    hcaComputationTask$invoke(
      dataList = dataList,
      filterHca = filterHca,
      distanceMeasure = distanceMeasure,
      clusterMethod = clusterMethod,
      minimumProportionOfLeafs = minimumProportionOfLeafs,
      minimumProportionToShowFragment = minimumProportionToShowFragment
    )
  }, error = function(e) {
    msg <- paste("Critical Error during HCA task invocation:", e$message)
    
    # Show error and cleanup
    output$information <- renderText(msg)
    session$sendCustomMessage("enableButton", "drawHCAplots")
    
    # Stop execution
    stop(e)
  })
  
  # processing after ExtendedTask is finished will be handled by the observer below
})

# Observer to handle ExtendedTask states
observe({
  status <- hcaComputationTask$status()
  
  if (status == "initial") {
    # waiting to be invoked
    
  } else if (status == "running") {
    # Task disable button and show modal spinner
    session$sendCustomMessage("disableButton", "drawHCAplots")
    
    # Show full-screen modal spinner
    shinybusy::show_modal_spinner(
      spin = "self-building-square",
      text = paste(
        "Computing HCA clusters...",
        paste("Distance measure:", input$hcaDistanceFunction),
        "This may take some time depending on data size.",
        "Please wait...",
        sep = "\n"
      ),
      color = "#5cb85c"
    )
    
    # Also update the info text (visible after spinner is removed)
    output$information <- renderText("⚙️ Computing HCA clusters... Please wait")
    output$tip <- renderText("HCA computation is running in the background. Please wait for completion.")
    
  } else if (status == "success") {
    # HCA computation was successfull
    # Remove modal spinner
    shinybusy::remove_modal_spinner()
    
    
    result <- tryCatch({
      hcaComputationTask$result()
    }, error = function(e) {
      stop(e)
    })
    
    # Extract results and assign to global variables
    currentDistanceMatrixObj <<- result$distanceMatrixObj
    clusterDataList <<- result$clusterDataList
    
    ##########################
    ## hca selections
    if(!is.null(listForTable_Fragment_PCA)){ ## selection from fragment
      fragmentIndex <- which(dataList$fragmentMasses == ms2PlotValues$fragmentListClicked$fragmentMasses[[selectionFragmentSelectedFragmentIndex]])
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
      ms2PlotValues$fragmentListClicked <<- NULL
      ms2PlotValues$fragmentListHovered <<- NULL
      dendrogramFragmentStatistics <<- FALSE
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
    updateChangePlotRadioButton()
    
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
    
    # Re-enable button
    session$sendCustomMessage("enableButton", "drawHCAplots")
    
  } else if (status == "error") {
    # Remove modal spinner
    shinybusy::remove_modal_spinner()
    
    # print error
    tryCatch({
      hcaComputationTask$result() 
    }, error = function(e) {
      msg <- e$message
      output$information <- renderText(msg)
      output$tip <- renderText("HCA computation failed. Please check the error message and try again.")
      
      # Re-enable button
      session$sendCustomMessage("enableButton", "drawHCAplots")
    })
  }
})

if(FALSE){
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
    plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramUserCoordinateRangeY
  )
  if(is.null(minimumIndex)){
    print(paste("Observe dendrogram hover", minimumIndex))
    if(!is.null(ms2PlotValues$fragmentListHovered)){
      ## reverse MS2 to clicked stuff
      ms2PlotValues$fragmentListHovered <<- NULL
    }
  } else {
    minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
    print(paste("Observe dendrogram hover i", minimumIndex, "l", minimumLabel))
    resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
    
    ## ## putative metabolite families statistics
    #putativeMetaboliteFamilies <- NULL
    #if(!is.null(classToSpectra_class)){
    #  putativeMetaboliteFamilies <- evaluatePutativeMetaboliteFamiliesOfDendrogramCluster(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel, classToSpectra_class = classToSpectra_class)
    #}
    
    #################################################
    ## output as message
    output$information <- renderText({
      print(paste("update output$information", resultObj$infoText))
      paste(
        resultObj$infoText, 
        #ifelse(test = is.null(putativeMetaboliteFamilies), yes = "", no = 
        #         paste("\n", paste(putativeMetaboliteFamilies, collapse = "\n"), sep = "")
        #), 
        sep = "")
    })
    
    if(all(!is.null(selectionAnalysisTreeNodeSet), minimumLabel == selectionAnalysisTreeNodeSet)){
      ## reverse MS2 to clicked stuff
      ms2PlotValues$fragmentListHovered <<- NULL
    } else {
      #################################################
      ## fetch ms2 spectrum
      ms2PlotValues$fragmentListHovered <<- resultObj
    }
  }
  
  ## MS2 plot
  #drawMS2Plot(consoleInfo = "dendrogram hover output$plotMS2")
  
  output$tip <- renderText({
    print(paste("update output$tip"))
    paste("Hover a cluster node or leaf node to view information about the corresponding MS\u00B9 feature cluster or MS\u00B9 feature respectively.", "Brush horizontally and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
  })
})
}
output$plotDendrogram_hover_info <- renderUI({
  hover <- input$plotDendrogram_hover
  hoverX <- hover$x
  hoverY <- hover$y
  plotWidth  <- session$clientData$output_plotDendrogram_width
  plotHeight <- session$clientData$output_plotDendrogram_height
  
  if(is.null(hoverX) | is.null(hoverY))
    return(NULL)
  
  output$tip <- renderText({
    print(paste("update output$tip"))
    paste("Hover a cluster node or leaf node to view information about the corresponding MS\u00B9 feature cluster or MS\u00B9 feature respectively.", "Brush horizontally and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
  })
  
  #################################################
  ## decide whether the hover is close enough to trigger event
  minimumIndex <- getSelectedPOI_XY(
    mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
    plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramUserCoordinateRangeY
  )
  print(paste("UI dendrogram hover", hoverX, hoverY, minimumIndex))
  
  if(is.null(minimumIndex)){
    print(paste("Observe dendrogram hover", minimumIndex))
    if(!is.null(ms2PlotValues$fragmentListHovered)){
      ## reverse MS2 to clicked stuff
      ms2PlotValues$fragmentListHovered <<- NULL
    }
    return(NULL)
  }
  
  minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
  print(paste("Observe dendrogram hover i", minimumIndex, "l", minimumLabel))
  resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
  
  ## ## putative metabolite families statistics
  #putativeMetaboliteFamilies <- NULL
  #if(!is.null(classToSpectra_class)){
  #  putativeMetaboliteFamilies <- evaluatePutativeMetaboliteFamiliesOfDendrogramCluster(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel, classToSpectra_class = classToSpectra_class)
  #}
  
  #################################################
  ## output as message
  #output$information <- renderText({
  #  print(paste("update output$information", resultObj$infoText))
  #  paste(
  #    resultObj$infoText, 
  #    #ifelse(test = is.null(putativeMetaboliteFamilies), yes = "", no = 
  #    #         paste("\n", paste(putativeMetaboliteFamilies, collapse = "\n"), sep = "")
  #    #), 
  #    sep = "")
  #})
  
  if(all(!is.null(selectionAnalysisTreeNodeSet), minimumLabel == selectionAnalysisTreeNodeSet)){
    ## reverse MS2 to clicked stuff
    ms2PlotValues$fragmentListHovered <<- NULL
  } else {
    #################################################
    ## fetch ms2 spectrum
    ms2PlotValues$fragmentListHovered <<- resultObj
  }
  
  ## compile information
  if(minimumLabel < 0){
    ## leaf node
    if(FALSE){
      resultObj$fragmentMasses <- fragmentsX
      resultObj$fragmentAbundances <- fragmentsY
      resultObj$fragmentColor <- fragmentsColor
      resultObj$fragmentDiscriminativity <- fragmentDiscriminativity
      resultObj$infoText <- infoText
      #resultObj$infoFeatureLabel <- featureID
      #resultObj$infoFragmentCount <- length(fragmentsX)
      resultObj$infoFamilies <- featureFamilies
      #resultObj$infoName <- featureName
      #resultObj$metFragLinkList <- metFragLinkList
      resultObj$precursorSet <- precursorSet
      resultObj$numberOfPrecursors <- numberOfPrecursors
    }
    info <- paste(
      "<b>MS\u00B9 feature: ", "</b>", resultObj$infoFeatureLabel, "<br>",
      "<b>Name: ", "</b>", resultObj$infoName, "<br>",
      "<b>Annotation: ", "</b>", paste(resultObj$infoFamilies, collapse = "; "), "<br>",
      "<b>Fragments: ", "</b>", resultObj$infoFragmentCount, "<br>",
      "<b>Max. cluster-disciminating power: ", "</b>", format(x = max(resultObj$fragmentDiscriminativity)*100, digits = 3, nsmall = 2), "%",
      sep = ""
    )
  } else {
    ## cluster node
    if(FALSE){
      resultObj$fragmentMasses <- fragmentsX
      resultObj$fragmentAbundances <- fragmentsY
      resultObj$fragmentColor <- fragmentsColor
      resultObj$fragmentDiscriminativity <- fragmentDiscriminativity
      #resultObj$clusterDiscriminativity <- clusterDiscriminativity
      resultObj$frequentFamilies <- frequentFamilies
      resultObj$infoText <- infoText
      #resultObj$metFragLinkList <- NULL
      resultObj$precursorSet <- precursorSet
      resultObj$numberOfPrecursors <- numberOfPrecursors
    }
    info <- paste(
      "<b>MS\u00B9 features: ", "</b>", resultObj$numberOfPrecursors, "<br>",
      "<b>Frequent annotations: ", "</b>", paste(resultObj$frequentFamilies, collapse = "; "), "<br>",
      "<b>Frequent fragments (", length(resultObj$fragmentMasses), "): ", "</b>", paste(format(x = resultObj$fragmentMasses, digits = 3, nsmall = 4), collapse = "; "), "<br>",
      "<b>Max. cluster-disciminating power: ", "</b>", format(x = resultObj$clusterDiscriminativity*100, digits = 3, nsmall = 2), "%",
      sep = ""
    )
  }
  
  panelWidth <- as.integer(plotWidth*0.6)
  showPlotTooltip(hover, info, panelWidth)
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
    plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramUserCoordinateRangeY
  )
  print(paste("Observe dendrogram click", is.null(minimumIndex), minimumIndex))
  
  if(is.null(minimumIndex)){
    ## reset stuff
    
    #ms2PlotValues$fragmentListClicked <<- NULL
    ms2PlotValues$fragmentListHovered <<- NULL
    
    selectionByAnalysisReset()
    selectionByFragmentReset()
    
    ##########################################################################################
    ## all fragments from the dendrogram
    
    ## get
    precursorSet <- filterHca$filter
    returnObj <- getSpectrumStatistics(dataList=dataList, precursorSet=precursorSet)
    fragmentMasses <- returnObj$fragmentMasses
    fragmentCounts <- returnObj$fragmentCounts
    
    ## filter
    theseFragments <- fragmentCounts >= minimumNumberOfPrecursorsForDendrogramStatistics
    fragmentMasses <- fragmentMasses[theseFragments]
    fragmentCounts <- fragmentCounts[theseFragments]
    
    ## set
    ms2PlotValues$fragmentListClicked <<- list(
      fragmentMasses = fragmentMasses,
      fragmentAbundances = fragmentCounts / length(precursorSet),
      fragmentColor = rep(x = "black", times = length(fragmentMasses)),
      fragmentDiscriminativity = rep(x = 0, times = length(fragmentMasses))
    )
    
    dendrogramFragmentStatistics <<- TRUE
    state$showPutativeAnnotationsTableFromAnalysis <<- FALSE
    putativeAnnotationsTableFromAnalysisCurrentlyShown <<- NULL
  } else {
    ## tree selection
    minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
    #minimumText <- clusterDataList$poiText[[minimumIndex]]
    
    ## fetch ms2 spectrum
    resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
    
    ## keep fragment selection
    selectionFragmentSelectedFragmentIndexNew <- NULL
    if(!is.null(selectionFragmentSelectedFragmentIndex)){
      fragmentMass <- ms2PlotValues$fragmentListClicked$fragmentMasses[[selectionFragmentSelectedFragmentIndex]]
      if(fragmentMass %in% resultObj$fragmentMasses)
        selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
    }
    
    ms2PlotValues$fragmentListClicked <<- resultObj
    ms2PlotValues$fragmentListHovered <<- NULL
    dendrogramFragmentStatistics <<- FALSE
    
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
    
    #################################################
    ## putative metabolite families statistics
    if(!is.null(classToSpectra_class)){
      precursorSet <- getPrecursorSetFromTreeSelection(clusterDataList = clusterDataList, clusterLabel = minimumLabel)
      setPutativeAnnotationsTableFromAnalysis(precursorSet)
      state$showPutativeAnnotationsTableFromAnalysis <<- TRUE
    }
  }## node selected
  
  #################################################
  ## plots
  
  ## cluster dendrogram
  ## TODO remove plot call
  drawDendrogramPlot(consoleInfo = "dendrogram click output$plotDendrogram", withHeatmap = TRUE)
  
  ## MS2 plot
  resetMS2PlotRange()
  #drawMS2Plot(consoleInfo = "dendrogram click output$plotMS2")
  
  #if(state$showPCAplotPanel)
  #  ## update PCA plots
  #  drawPcaPlots(consoleInfo = "dendrogram click output$plotPcaScores")
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

if(FALSE){
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
  treeLeafIndex2 <- as.numeric(format(x = hoverX, digits = 1))
  treeLeafIndex  <- clusterDataList$cluster$order[[treeLeafIndex2]]
  precursorIndex <- filterHca$filter[[treeLeafIndex]]
  
  msg <- list()
  msg[[length(msg) + 1]] <- "MS\u00B9 feature: "
  msg[[length(msg) + 1]] <- dataList$precursorLabels[[precursorIndex]]
  msg[[length(msg) + 1]] <- "\n"
  
  if(hoverY > 2){
    ## lcf
    groupOne <- filterHca$sampleClasses[[1]]
    groupTwo <- filterHca$sampleClasses[[2]]
    msg[[length(msg) + 1]] <- paste("log-fold-change = log_2( mean(group ", groupOne, ") / mean(group ", groupTwo, ") )", sep = "")
    
    valMeanOne <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupOne)]
    valMeanTwo <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupTwo)]
    msg[[length(msg) + 1]] <- " = log_2( "
    msg[[length(msg) + 1]] <- as.numeric(format(x = valMeanOne, digits = 2))
    msg[[length(msg) + 1]] <- " / "
    msg[[length(msg) + 1]] <- as.numeric(format(x = valMeanTwo, digits = 2))
    msg[[length(msg) + 1]] <- " ) = "
    
    lfc <- dataList$dataFrameMeasurements[precursorIndex, dataList$lfcColumnNameFunctionFromName(groupOne, groupTwo)]
    msg[[length(msg) + 1]] <- as.numeric(format(x = lfc, digits = 2))
  } else {
    if(hoverY > 1){
      ## group 1
      groupHere <- filterHca$sampleClasses[[1]]
    } else { ## hoverY <= 1
      ## group 2
      groupHere <- filterHca$sampleClasses[[2]]
    }
    msg[[length(msg) + 1]] <- paste("Mean abundance of group ", groupHere, ": ", sep = "")
    valMean <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupHere)]
    msg[[length(msg) + 1]] <- as.numeric(format(x = valMean, digits = 2))
    msg[[length(msg) + 1]] <- " = mean("
    #columnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = groupHere, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
    columnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = groupHere, sampleNamesToExclude = dataList$excludedSamples(groupSampleDataFrame = dataList$groupSampleDataFrame, sampleClasses = groupHere))
    vals <- dataList$dataFrameMeasurements[precursorIndex, columnNames]
    vals <- as.numeric(format(x = vals, digits = 2))
    msg[[length(msg) + 1]] <- paste(vals, collapse = ", ")
    msg[[length(msg) + 1]] <- ")"
  }
  
  output$information <- renderText({
    print(paste("update output$information heatmap hover", sep = ""))
    paste(msg, collapse = "")
  })
})
}  ## not finished

#print("entering the line ..line 680")
#print(columnsOfInterestForHeatmap)

output$plotHeatmap_hover_info <- renderUI({
  hover <- input$plotHeatmap_hover
  hoverX <- hover$x
  hoverY <- hover$y
  plotWidth  <- session$clientData$output_plotHeatmap_width
  plotHeight <- session$clientData$output_plotHeatmap_height
  
  columnsOfInterest <- columnsOfInterestForHeatmap
  
  if(is.null(hoverX) | is.null(hoverY))
    return(NULL)
  if(hoverX < 0.5 | hoverX > (clusterDataList$numberOfPrecursorsFiltered + 0.5))
    return(NULL)
  if(hoverY < 0 | hoverY > length(columnsOfInterest))
    return(NULL)
  
  #################################################
  treeLeafIndex2 <- as.numeric(format(x = hoverX, digits = 1))
  treeLeafIndex  <- clusterDataList$cluster$order[[treeLeafIndex2]]
  precursorIndex <- filterHca$filter[[treeLeafIndex]]
  
  ## differentiate heatmap content: LFC, samples, sampleClasses
  columnOfInterest <- columnsOfInterest[[ceiling(hoverY)]]
  #print("entering the line ...line 707")
  #print(columnOfInterest)
  #columnOfInterest <- rev(columnsOfInterest)[[round(hoverY)]]
  #####################
  #if(startsWith(x = columnOfInterest, prefix = "LFC_")){
  if(startsWith(columnOfInterest, "LFC_")){
    ####################################
    ## lcf
    #print("enter the line ...line 714")
    ########
    sampleClasses <- dataList$lfcColumnNameFunctionFromString(columnOfInterest)
    #######
    #print(columnOfInterest)
    #print(sampleClasses)
    #######
    groupOne <- sampleClasses[[1]]
    groupTwo <- sampleClasses[[2]]
    #msg[[length(msg) + 1]] <- paste("log-fold-change = log_2( mean(group ", groupOne, ") / mean(group ", groupTwo, ") )", sep = "")
    
    valMeanOne <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupOne)]
    valMeanTwo <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupTwo)]
    valMeanOne <- format(x = valMeanOne, digits = 2)
    valMeanTwo <- format(x = valMeanTwo, digits = 2)
    #msg[[length(msg) + 1]] <- " = log_2( "
    #msg[[length(msg) + 1]] <- valMeanOne
    #msg[[length(msg) + 1]] <- " / "
    #msg[[length(msg) + 1]] <- valMeanTwo
    #msg[[length(msg) + 1]] <- " ) = "
    
    lfc <- dataList$dataFrameMeasurements[precursorIndex, dataList$lfcColumnNameFunctionFromName(groupOne, groupTwo)]
    lfc <- format(x = lfc, digits = 2)
    #msg[[length(msg) + 1]] <- lfc
    
    #samplesOne <- dataList$dataColumnsNameFunctionFromGroupName(group = groupOne, sampleNamesToExclude = dataList$excludedSamples(groupSampleDataFrame = dataList$groupSampleDataFrame, sampleClasses = groupOne))
    #samplesTwo <- dataList$dataColumnsNameFunctionFromGroupName(group = groupTwo, sampleNamesToExclude = dataList$excludedSamples(groupSampleDataFrame = dataList$groupSampleDataFrame, sampleClasses = groupTwo))
    
    info <- paste(
      "<b>MS\u00B9 feature: ", "</b>", dataList$precursorLabels[[precursorIndex]], "<br>",
      "<b>Group ratio: ", "</b>", groupOne, " / ", groupTwo, "<br>",
      "<b>log2-fold-change: ", "</b>", lfc, "<br>",
      #"<b>Samples in groups: ", "</b>", paste(c(samplesOne, samplesTwo), collapse = ", "), "<br>",
      "<b>Mean abundances: ", "</b>", valMeanOne, " / ", valMeanTwo, #"<br>",
      #"<b>Sample abundances: ", "</b>", paste(vals, collapse = ", "),
      sep = ""
    )
  } #else if(startsWith(x = columnOfInterest, prefix = "HBR")){}
  else {
    #if(endsWith(x = columnOfInterest, suffix = "_mean")){
    if(endsWith(columnOfInterest, "_mean")){
      #####################################
      ## group
      #print("enter the line 750")
      #####
      groupHere <- dataList$dataMeanColumnNameFunctionFromString(columnOfInterest)
      
      #msg[[length(msg) + 1]] <- paste("Mean abundance of group ", groupHere, ": ", sep = "")
      valMean <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupHere)]
      valMean <- format(x = valMean, digits = 2)
      #msg[[length(msg) + 1]] <- valMean
      #msg[[length(msg) + 1]] <- " = mean("
      #columnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = groupHere, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
      columnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = groupHere, sampleNamesToExclude = dataList$excludedSamples(groupSampleDataFrame = dataList$groupSampleDataFrame, sampleClasses = groupHere))
      vals <- dataList$dataFrameMeasurements[precursorIndex, columnNames]
      vals <- format(x = vals, digits = 2)
      #msg[[length(msg) + 1]] <- paste(vals, collapse = ", ")
      #msg[[length(msg) + 1]] <- ")"
      
      info <- paste(
        "<b>MS\u00B9 feature: ", "</b>", dataList$precursorLabels[[precursorIndex]], "<br>",
        "<b>Sample group: ", "</b>", groupHere, "<br>",
        "<b>Mean abundance: ", "</b>", valMean, "<br>",
        "<b>Samples in group: ", "</b>", paste(columnNames, collapse = ", "), "<br>",
        "<b>Sample abundances: ", "</b>", paste(vals, collapse = ", "),
        sep = ""
      )
    } else {
      #####################################
      ## sample
      val <- dataList$dataFrameMeasurements[precursorIndex, columnOfInterest]
      #msg[[length(msg) + 1]] <- paste("Abundance in sample ", columnOfInterest, ": ", val, sep = "")
      
      info <- paste(
        "<b>MS\u00B9 feature: ", "</b>", dataList$precursorLabels[[precursorIndex]], "<br>",
        "<b>Sample: ", "</b>", columnOfInterest, "<br>",
        "<b>Abundance: ", "</b>", format(x = val, digits = 2),
        sep = ""
      )
    }
  }
 #################################### 
  #output$information <- renderText({
  #  print(paste("update output$information heatmap hover", sep = ""))
  #  paste(msg, collapse = "")
  #})
  
  panelWidth <- as.integer(plotWidth*0.4)
  showPlotTooltip(hover, info, panelWidth)
})

setPutativeAnnotationsTableFromAnalysis <- function(precursorSet){
  if(FALSE){
    dataList_ = dataList
    precursorSet_ = precursorSet
    classToSpectra_class_ = classToSpectra_class
  }
  if(FALSE){
    dataList = dataList_
    precursorSet = precursorSet_
    classToSpectra_class = classToSpectra_class_
  }
  
  ## data.frame("Class"=character(), "pValue"=numeric(), "ProportionInPercent" = numeric())
  returnObj  <- evaluatePutativeMetaboliteFamiliesOfPrecursorSet(dataList = dataList, precursorSet = precursorSet, classToSpectra_class = classToSpectra_class)
  overviewDf <- returnObj$overviewDf
  detailDf   <- returnObj$detailDf
  
  if(nrow(overviewDf) == 0){
    msg <- ifelse(test = length(precursorSet) == 1, yes = "There are no putative annotations.", no = "There are no frequent putative annotations.")
    showDf <- data.frame("Message" = msg, stringsAsFactors = FALSE)
    print(paste(length(precursorSet), "precursors without putative annotation"))
  } else {
    type <- "PutativeAnnotationsTableFromAnalysis"
    annotateColumn <- createActionButtonInputFields(
      FUN = actionButton,  id = paste(type, "Annotate", sep = "_"), itemCount=nrow(overviewDf), 
      label   = "Annotate", tableCounter = putativeAnnotationsTableFromAnalysisInputFieldIdCounter, 
      callFunction = putativeAnnotationsTableFromAnalysisAnnotateClicked
    )
    putativeAnnotationsTableFromAnalysisInputFieldIdCounter <<- putativeAnnotationsTableFromAnalysisInputFieldIdCounter + 1
    
    if(length(precursorSet) == 1)
      #overviewDf <- data.frame("Class" = df$Class, "pValue" = df$pValue, stringsAsFactors = FALSE)
      showDf <- cbind("Class" = overviewDf$Class, "pValue" = overviewDf$pValue, "Annotate" = annotateColumn)
    else
      showDf <- cbind("Class" = overviewDf$Class, "median pValue" = overviewDf$pValue, "Annotate" = annotateColumn)
    print(paste(length(precursorSet), "precursors with", nrow(detailDf), "putative annotations in", nrow(showDf), "classes"))
  }
  
  putativeAnnotationsTableFromAnalysisCurrentlyShown <<- showDf
  
  output$putativeAnnotationsTableFromAnalysis <- DT::renderDataTable(
    expr = showDf,
    server = FALSE, escape = FALSE, selection = ifelse(test = nrow(overviewDf)>0, yes = "single", no = "none"), #rownames = FALSE,
    options = list(
      #scrollY = "600px",
      scrollY = "40vh",
      preDrawCallback = DT::JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = DT::JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
      iDisplayLength=nrow(showDf),       # initial number of records
      #aLengthMenu = c(5,10),    # records/page options
      #bLengthChange =0,         # show/hide records per page dropdown
      #bFilter = 0,              # global search box on/off
      #bInfo = 0,                # information on/off (how many records filtered, etc)
      #bAutoWidth = 0,           # automatic column width calculation, disable if passing column width via aoColumnDefs
      #aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
      ordering = F,              # row ordering
      sDom  = 't'
      #sDom  = '<"top">rt<"bottom">ip'
    )
  )
}

#obsPutativeAnnotationsTableFromAnalysis_rows_selected <- observeEvent(eventExpr = input$putativeAnnotationsTableFromAnalysis_rows_selected, handlerExpr = {
#  print(paste("Observe putativeAnnotationsTableFromAnalysis_rows_selected", input$putativeAnnotationsTableFromAnalysis_rows_selected))
#  state$putativeAnnotationsTableFromAnalysisRowSelected <<- !is.null(input$putativeAnnotationsTableFromAnalysis_rows_selected)
#}, ignoreNULL = FALSE)
putativeAnnotationsTableFromAnalysisAnnotateClicked <- function(buttonId){
  ## PutativeAnnotationsTableFromAnalysis_Annotate_0_1
  rowIdx    <- as.integer(strsplit(x = buttonId, split = "_")[[1]][[4]])
  
  precursorSet <- getPrecursorSetFromTreeSelection(clusterDataList = clusterDataList, clusterLabel = selectionAnalysisTreeNodeSet)
  class <- putativeAnnotationsTableFromAnalysisCurrentlyShown[[rowIdx, "Class"]]
  subClass <- tail(x = strsplit(x = class, split = "; ")[[1]], n = 1)
  print(paste(length(precursorSet), "-->", class, "(", subClass, ")", "from row", rowIdx))
  
  callbackFunction <- function(value, color){
    print(paste(length(precursorSet), "-->", class, "from row", rowIdx, "-->", value, color))
    removeModal(session = session)
    
    ## add
    addAnnotation(precursorSet = precursorSet, annotationValue = value, annotationColor = color)
  }
  
  openAnnotaionNameColorDialog(predefinedClassName = subClass, callbackFunction = callbackFunction)
}

if(FALSE){
  ## plotly
  obsDendLabelsHeatmap <- observeEvent(event_data("plotly_click", source = "dendLabelsHeatmap"), {
    eventdata <- event_data("plotly_click", source = "dendLabelsHeatmap")
    #str(eventdata)
    ## 'data.frame':	1 obs. of  4 variables:
    ## $ curveNumber: int 7
    ## $ pointNumber: int 6
    ## $ x          : num 8.5
    ## $ y          : num 0.697
    
    if(is.null(eventdata))
      return()
    
    curveName <- curveNumberToCurveName[which(curveNumberToCurveName[, "curveNumber"] == eventdata$curveNumber), "name"]
    #print(paste("curveName:", curveName))
    if(curveName != "nodes")
      return()
    
    #nodeIndex <- eventdata$pointNumber + 1
    nodeIndex <- which((
      #eventdata$x == clusterDataList$poiCoordinatesX & 
      #eventdata$x %in% clusterDataList$poiCoordinatesX
      abs(eventdata$x - clusterDataList$poiCoordinatesX) <= 0.0001
    ) & (
      #eventdata$y == clusterDataList$poiCoordinatesY
      #eventdata$y %in% clusterDataList$poiCoordinatesY
      abs(eventdata$y - clusterDataList$poiCoordinatesY) <= 0.0001
    ))
    nodeLabel <- clusterDataList$poiLabels[[nodeIndex]]
    
    
    ## fetch ms2 spectrum
    resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = nodeLabel)
    
    ## keep fragment selection
    selectionFragmentSelectedFragmentIndexNew <- NULL
    if(!is.null(selectionFragmentSelectedFragmentIndex)){
      fragmentMass <- ms2PlotValues$fragmentListClicked$fragmentMasses[[selectionFragmentSelectedFragmentIndex]]
      if(fragmentMass %in% resultObj$fragmentMasses)
        selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
    }
    
    ms2PlotValues$fragmentListClicked <<- resultObj
    ms2PlotValues$fragmentListHovered <<- NULL
    dendrogramFragmentStatistics <<- FALSE
    
    print(paste("fragments:", paste(resultObj$fragmentMasses, collapse = "; ")))
    
    #################################################
    ## output as message
    selectionByHca(nodeLabel)
    
    ## update the selected fragment in case of overlapping spectra
    if(!is.null(selectionFragmentSelectedFragmentIndexNew))
      selectionByFragment(selectionFragmentSelectedFragmentIndexNew)
    else
      selectionByFragmentReset()
    
    output$information <- renderText({
      print(paste("update output$information", resultObj$infoText))
      paste(resultObj$infoText, sep = "")
    })
    
    ## TODO
    #drawMS2PlotImpl
    
    resetMS2PlotRange()
    #drawMS2Plot()
    
    #output$plotDendrogram <- renderPlotly({
  })
  
  #output$plotTmp <- renderPlotly({
  #  
  #  # Read in hover data
  #  eventdata <- event_data("plotly_click", source = "dendLabelsHeatmap")
  #  print("haha")
  #  str(eventdata)
  #  validate(need(!is.null(eventdata), "Hover over the time series chart to populate this heatmap"))
  #  
  #  curveName <- curveNumberToCurveName[which(curveNumberToCurveName[, "curveNumber"] == eventdata$curveNumber), "name"]
  #  print(curveName)
  #  
  #  plot_ly(x = 1:10, y = rnorm(10), type = "scatter", mode = "markers")
  #})
}

output$showPutativeAnnotationsTableFromAnalysis <- reactive({
  print(paste("reactive update showPutativeAnnotationsTableFromAnalysis", state$showPutativeAnnotationsTableFromAnalysis))
  return(state$showPutativeAnnotationsTableFromAnalysis)
})

outputOptions(output, 'showPutativeAnnotationsTableFromAnalysis',  suspendWhenHidden=FALSE)

output$ui_plotAnnoLegendHCA <- renderUI({
  print(paste("### GUI ### ui_plotAnnoLegendHCA"))
  if(state_tabHca$annotationLegendHeightHca != -1) {
    plotOutput(outputId = "plotAnnoLegendHCA", height = state_tabHca$annotationLegendHeightHca)
  }
})
output$ui_plotHeatmap <- renderUI({
  print(paste("### GUI ### ui_plotHeatmap"))
  plotOutput(height = state_tabHca$heatmapHeight, 
             outputId = "plotHeatmap",
             #hover    = "plotHeatmap_hover", 
             hover    = hoverOpts(
               id = "plotHeatmap_hover",
               delay = 50, 
               delayType = "debounce"
             )
             #click = "plotHeatmap_click"
  )
})


suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabHca observers")
  obsDrawHCA$suspend()
  obsDendrogramClick$suspend()
  obsDendrogramdblClick$suspend()
})
