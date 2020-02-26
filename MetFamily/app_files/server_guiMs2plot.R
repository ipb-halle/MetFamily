
##############################################
## MS2 peaks

## resultObj$fragmentMasses <- fragmentsX
## resultObj$fragmentAbundances <- fragmentsY
## resultObj$fragmentColor <- fragmentsColor
## resultObj$fragmentDiscriminativity <- fragmentDiscriminativity
## ?? resultObj$clusterDiscriminativity <- clusterDiscriminativity
## ? resultObj$infoText <- infoText
## ? resultObj$metFragLinkList <- NULL
## ? resultObj$precursorSet <- precursorSet
## ? resultObj$numberOfPrecursors <- numberOfPrecursors
ms2PlotValues <- reactiveValues(
  fragmentListClicked = NULL,
  fragmentListHovered = NULL
)
ms2PlotRange         <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)

resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset ms2plot state")
  ## MS2 plot range
  resetMS2PlotRange()
  ## fragments
  ms2PlotValues$fragmentListHovered <<- NULL
  ms2PlotValues$fragmentListClicked <<- NULL
})


drawMS2Plot <- function(consoleInfo = NULL){
  output$plotMS2 <- renderPlot({
    print(paste("### MS2 ###", consoleInfo))
    drawMS2PlotImpl()
  })
}
drawMS2Legend <- function(consoleInfo = NULL){
  output$plotMS2Legend <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawMS2LegendImpl()
  })
}
drawFragmentDiscriminativityLegend <- function(consoleInfo = NULL){
  output$plotFragmentDiscriminativityLegend <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawFragmentDiscriminativityLegendImpl()
  })
}

resetMS2PlotRange <- function(){
  ms2PlotRange$xMin <<- dataList$minimumMass
  ms2PlotRange$xMax <<- dataList$maximumMass
  ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
  ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
}

if(FALSE){
obsMS2hover <- observeEvent(input$plotMS2_hover, {
  hoverX <- input$plotMS2_hover$x
  hoverY <- input$plotMS2_hover$y
  plotWidth  <- session$clientData$output_plotMS2_width
  plotHeight  <- session$clientData$output_plotMS2_height
  
  if(is.null(hoverX) | is.null(hoverY))
    return()
  if(any(is.null(ms2PlotValues$fragmentListClicked), length(ms2PlotValues$fragmentListClicked$fragmentMasses) == 0))
    return()
  
  ################################################
  ## decide whether the click is close enough to trigger event
  
  minimumIndex <- getSelectedPOI_XY(
    mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = ms2PlotValues$fragmentListClicked$fragmentMasses, poiCoordinatesY = ms2PlotValues$fragmentListClicked$fragmentAbundances, 
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
    fragmentIndex <- which(dataList$fragmentMasses == ms2PlotValues$fragmentListClicked$fragmentMasses[[minimumIndex]])
    
    if(state$plotHcaShown)
      numberOfPrecursors <- sum(dataList$featureMatrix[filterHca$filter, fragmentIndex] != 0)
    if(state$plotPcaShown)
      numberOfPrecursors <- sum(dataList$featureMatrix[filterPca$filter, fragmentIndex] != 0)
    
    output$information <- renderText({
      print(paste("update output$information"))
      paste(
        "Fragment with m/z = ", ms2PlotValues$fragmentListClicked$fragmentMasses[[minimumIndex]], 
        " and (average) abundance = ", format(x = ms2PlotValues$fragmentListClicked$fragmentAbundances[[minimumIndex]], digits = 0, nsmall = 4), 
        " is present in ", numberOfPrecursors, " MS/MS spectra",
        "\nand has a cluster-discriminating power of ", format(x = ms2PlotValues$fragmentListClicked$fragmentDiscriminativity[[minimumIndex]]*100, digits = 3, nsmall = 2), "%.", 
        sep = ""
      )
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
}
output$plotMS2_hover_info <- renderUI({
  hover <- input$plotMS2_hover
  hoverX <- hover$x
  hoverY <- hover$y
  plotWidth  <- session$clientData$output_plotMS2_width
  plotHeight <- session$clientData$output_plotMS2_height
  
  if(is.null(hoverX) | is.null(hoverY))
    return()
  if(any(is.null(ms2PlotValues$fragmentListClicked), length(ms2PlotValues$fragmentListClicked$fragmentMasses) == 0))
    return()
  
  output$tip <- renderText({
    print(paste("update output$tip"))
    paste(
      "Hover or click a fragment node (only 'Fragments from selection') to view information about this fragment.", 
      "Brush horizontally and double-click to zoom in.", 
      "Double-click to zoom out.", 
      sep = "\n"
    )
  })
  ################################################
  ## decide whether the click is close enough to trigger event
  
  minimumIndex <- getSelectedPOI_XY(
    mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = ms2PlotValues$fragmentListClicked$fragmentMasses, poiCoordinatesY = ms2PlotValues$fragmentListClicked$fragmentAbundances, 
    plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = ms2PlotRange$xIntervalSize, plotRangeY = 1
  )
  print(paste("Observe MS2 hover", minimumIndex))
  if(is.null(minimumIndex)){
    ## nothing in selection range
    #output$information <- renderText({
    #  print(paste("update output$information"))
    #  paste("", sep = "")
    #})
    return()
  }
  
  ## point selected
  fragmentIndex <- which(dataList$fragmentMasses == ms2PlotValues$fragmentListClicked$fragmentMasses[[minimumIndex]])
  
  if(state$plotHcaShown){
    numberOfPrecursors <- sum(dataList$featureMatrix[filterHca$filter, fragmentIndex] != 0)
    currentlyShownAnalysisype <- "HCA"
  }
  if(state$plotPcaShown){
    numberOfPrecursors <- sum(dataList$featureMatrix[filterPca$filter, fragmentIndex] != 0)
    currentlyShownAnalysisype <- "PCA"
  }
  
  #output$information <- renderText({
  #  print(paste("update output$information"))
  #  paste(
  #    "Fragment with m/z = ", ms2PlotValues$fragmentListClicked$fragmentMasses[[minimumIndex]], 
  #    " and (average) abundance = ", format(x = ms2PlotValues$fragmentListClicked$fragmentAbundances[[minimumIndex]], digits = 0, nsmall = 4), 
  #    " is present in ", numberOfPrecursors, " MS/MS spectra",
  #    "\nand has a cluster-discriminating power of ", format(x = ms2PlotValues$fragmentListClicked$fragmentDiscriminativity[[minimumIndex]]*100, digits = 3, nsmall = 2), "%.", 
  #    sep = ""
  #  )
  #})
  
  fragmentMz <- ms2PlotValues$fragmentListClicked$fragmentMasses[[minimumIndex]]
  fragmentMzS <- format(x = fragmentMz, nsmall = 4)
  intensityS <- format(x = ms2PlotValues$fragmentListClicked$fragmentAbundances[[minimumIndex]], digits = 0, nsmall = 4)
  cdpS <- format(x = ms2PlotValues$fragmentListClicked$fragmentDiscriminativity[[minimumIndex]]*100, digits = 3, nsmall = 2)
  
  info <- paste(
    "<b>", ifelse(test = fragmentMz < 0, yes = "Neutral loss", no = "Fragment"), ": ", "</b>", fragmentMzS, "<br>",
    "<b>Average intensity: ", "</b>", intensityS, "<br>",
    "<b>MS/MS spectra: ", "</b>", numberOfPrecursors, " (in ", currentlyShownAnalysisype, ")", "<br>",
    "<b>Cluster-discriminating power: ", "</b>", cdpS, "%", 
    sep = ""
  )
  
  panelWidth <- as.integer(plotWidth*0.6)
  showPlotTooltip(hover, info, panelWidth)
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
  if(any(is.null(ms2PlotValues$fragmentListClicked), length(ms2PlotValues$fragmentListClicked$fragmentMasses) == 0))
    return()
  
  #################################################
  ## decide whether the click is close enough to trigger event
  minimumIndex <- getSelectedPOI_XY(
    mouseX = clickX, mouseY = clickY, poiCoordinatesX = ms2PlotValues$fragmentListClicked$fragmentMasses, poiCoordinatesY = ms2PlotValues$fragmentListClicked$fragmentAbundances,
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
    drawDendrogramPlot(consoleInfo = "MS2 click output$plotDendrogram", withHeatmap = TRUE)
  ## PCA
  if(state$showPCAplotPanel)
    drawPcaLoadingsPlot(consoleInfo = "MS2 click output$plotPcaLoadings")
  
  ## update node selection
  drawMS2Plot(consoleInfo = "MS2 click output$plotMS2")
})
obsMS2dblClick <- observeEvent(input$plotMS2_dblclick, {
  brush <- input$plotMS2_brush
  
  if(all(any(is.null(ms2PlotValues$fragmentListClicked), length(ms2PlotValues$fragmentListClicked$fragmentMasses) == 0), any(is.null(ms2PlotValues$fragmentListHovered), length(ms2PlotValues$fragmentListHovered$fragmentMasses) == 0)))
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


suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending ms2Plot observers")
  obsMS2click$suspend()
  obsMS2dblClick$suspend()
})
