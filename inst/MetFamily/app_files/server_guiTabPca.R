
scoresGroupsLegendEntryHeight <- 20

pcaDataList <- NULL

scoresGroups <- NULL

pcaScoresPlotRange   <- reactiveValues(
  xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL, 
  yMin = NULL, yMax = NULL, yInterval = NULL, yIntervalSize = NULL
)
pcaLoadingsPlotRange <- reactiveValues(
  xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL, 
  yMin = NULL, yMax = NULL, yInterval = NULL, yIntervalSize = NULL
)


state_tabPca <- reactiveValues(
  ## plot controls
  ## changing showScoresLabels = TRUE to FALSE
  showScoresLabels = FALSE,
  loadingsLabels = "None",
  showLoadingsFeaturesAnnotated = TRUE, 
  showLoadingsFeaturesUnannotated = TRUE, 
  showLoadingsFeaturesSelected = TRUE, 
  showLoadingsFeaturesUnselected = TRUE,
  showLoadingsAbundance = FALSE,
  ## annotation legend height
  annotationLegendHeightPca = -1,
  ## plot annotations: $setOfAnnotations, $setOfColors
  annotationsPca = NULL,
  scoresGroupsLegendHeight = -1
)
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset tabPca state")
  ## plot controls
  ###showScoresLabels from TRUE to FALSE
  state_tabPca$showScoresLabels <<- FALSE
  state_tabPca$loadingsLabels   <<- "None"
  state_tabPca$showLoadingsFeaturesAnnotated   <<- TRUE 
  state_tabPca$showLoadingsFeaturesUnannotated <<- TRUE 
  state_tabPca$showLoadingsFeaturesSelected    <<- TRUE 
  state_tabPca$showLoadingsFeaturesUnselected  <<- TRUE
  state_tabPca$showLoadingsAbundance           <<- FALSE
  ## annotation legend height
  state_tabPca$annotationLegendHeightPca <<- -1
  ## plot annotations: $setOfAnnotations, $setOfColors
  state_tabPca$annotationsPca <<- NULL
  state_tabPca$scoresGroupsLegendHeight <<- -1
  
  ## reset variables
  pcaDataList <<- NULL
  scoresGroups <<- NULL
  
  ## reset plot range
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
  
  ## number of components for PCA
  maximumNumberOfComponents <- length(dataList$includedSamples(dataList$groupSampleDataFrame)) - 1
  maximumNumberOfComponents <- min(maximumNumberOfComponents, 5)
  if(maximumNumberOfComponents < 2){
    updateSelectInput(session = session, inputId = "pcaDimensionOne", choices = "1")
    updateSelectInput(session = session, inputId = "pcaDimensionTwo", choices = "1")
    shinyjs::disable(id = "pcaDimensionOne")
    shinyjs::disable(id = "pcaDimensionTwo")
  } else {
    updateSelectInput(session = session, inputId = "pcaDimensionOne", choices = seq_len(maximumNumberOfComponents), selected = 1)
    updateSelectInput(session = session, inputId = "pcaDimensionTwo", choices = seq_len(maximumNumberOfComponents), selected = 2)
  }
})


obsDrawPCA <- observeEvent(input$drawPCAplots, {
  session$sendCustomMessage("disableButton", "drawPCAplots")
  #################################################
  ## get input
  drawPlots <- as.numeric(input$drawPCAplots)
  
  print(paste("Observe draw PCA plots", drawPlots))
  
  ## check if button was hit
  #if(drawPlots == drawPCAButtonValue)
  #  return()
  #drawPCAButtonValue <<- drawPlots
  
  ms1AnalysisMethod <- input$ms1AnalysisMethod
  pcaScaling      <- input$pcaScaling
  pcaLogTransform <- input$pcaLogTransform
  pcaDimensionOne <<- as.numeric(input$pcaDimensionOne)
  pcaDimensionTwo <<- as.numeric(input$pcaDimensionTwo)
  
  calculatePca(ms1AnalysisMethod, pcaScaling, pcaLogTransform, pcaDimensionOne, pcaDimensionTwo)
})

calculatePca <- function(ms1AnalysisMethod, pcaScaling, pcaLogTransform, pcaDimensionOne, pcaDimensionTwo){
  #################################################
  ## calc PCA
  pca <- calculatePCA(dataList = dataList, filterObj = filterPca, ms1AnalysisMethod = ms1AnalysisMethod, scaling = pcaScaling, logTransform = pcaLogTransform)
  print("enter the line ...line number 114")
  #print(pca)
  #print(colnames(as.data.frame(dataList$dataFrameInfos)))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"m/z"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"RT"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"Annotation"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"Average Rt(min)"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"Average Mz"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"Metabolite name"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"Adduct ion name"))
  #print(table(as.data.frame(dataList$dataFrameInfos)$"MySample"))
  ########
  #print(as.data.frame(table(unlist(dataList))))
  #########################################
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
    fragmentIndex <- which(dataList$fragmentMasses == ms2PlotValues$fragmentListClicked$fragmentMasses[[selectionFragmentSelectedFragmentIndex]])
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
    ms2PlotValues$fragmentListClicked <<- NULL
    ms2PlotValues$fragmentListHovered <<- NULL
    dendrogramFragmentStatistics <<- FALSE
  }
  
  ##########################
  ## draw
  resetPcaPlotRange()
  drawPcaPlots(consoleInfo = "drawPCA output$plotPcaScores")
  drawMS2Plot(consoleInfo = "drawPCA output$plotMS2")
  drawAnnotationLegendPCA(consoleInfo = "init output$plotAnnoLegend")
  
  scoresGroups <<- list(
    grouXXXps = filterPca$grouXXXps,
    colors = colorPaletteScores()[unlist(lapply(X = filterPca$grouXXXps, FUN = dataList$groupIdxFromGroupName))]
  )
  state_tabPca$scoresGroupsLegendHeight <<- scoresGroupsLegendEntryHeight * (length(scoresGroups$grouXXXps) + 1)
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
  updateChangePlotRadioButton()
  
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
  session$sendCustomMessage("enableButton", "drawPCAplots")
}
drawPcaPlots <- function(consoleInfo = NULL){
  drawPcaScoresPlot(consoleInfo = consoleInfo)
  drawPcaLoadingsPlot(consoleInfo = consoleInfo)
}
drawPcaScoresPlot <- function(consoleInfo = NULL){
  output$plotPcaScores <- renderPlot({
    print(paste("### psp ###", consoleInfo))
    drawPcaScoresPlotImpl()
  })
}
drawPcaLoadingsPlot <- function(consoleInfo = NULL){
  output$plotPcaLoadings <- renderPlot({
    print(paste("### plp ###", consoleInfo))
    drawPcaLoadingsPlotImpl()
  })
}

drawAnnotationLegendPCA <- function(consoleInfo = NULL){
  output$plotAnnoLegendPCA <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawAnnotationLegendPCAimpl()
  })
}
drawScoresGroupsLegend <- function(consoleInfo = NULL){
  output$plotScoresGroupsLegend <- renderPlot({
    print(paste("### leg ###", consoleInfo))
    drawScoresGroupsLegendImpl()
  })
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

if(FALSE){
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
    dataColumnName <- filterPca$sampleSet[[minimumIndex]]
    #dataColumnName <- dataList$dataColumnsNameFunctionFromGroupNames(grouXXXps = filterPca$grouXXXps, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))[[minimumIndex]]
    group <- dataList$groupNameFunctionFromDataColumnName(dataColumnName = dataColumnName, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
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
}
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
    
    ms2PlotValues$fragmentListClicked <<- NULL
    ms2PlotValues$fragmentListHovered <<- NULL
    dendrogramFragmentStatistics <<- FALSE
    
    selectionByAnalysisReset()
    selectionByFragmentReset()
  } else {
    ## loadng selection
    precursorIndex <- filterPca$filter[[minimumIndex]]
    ## fetch ms2 spectrum
    resultObj <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)
    
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
  
  ## TODO remove plot call
  drawPcaLoadingsPlot(consoleInfo = "PCA loadings click output$plotPcaLoadings")
  
  if(state$showHCAplotPanel)
    ## update dendrogram plot
    drawDendrogramPlot(consoleInfo = "PCA loadings click output$plotDendrogram", withHeatmap = TRUE)
})
if(FALSE){
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
    ## no loading hovered
    ms2PlotValues$fragmentListHovered <<- NULL
    #output$information <- renderText({
    #  print(paste("update output$information PCA Loadings hover ", minimumIndex, sep = ""))
    #  paste("", sep = "")
    #})
  } else {
    print(paste("Observe PCA Loadings hover", hoverX, hoverY, minimumIndex))
    
    #################################################
    ## fetch ms2 spectrum
    precursorIndex <- filterPca$filter[[minimumIndex]]
    resultObj <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)
    if(all(!is.null(selectionAnalysisPcaLoadingSet), minimumIndex == selectionAnalysisPcaLoadingSet)){
      ## reverse MS2 to clicked stuff
      ms2PlotValues$fragmentListHovered <<- NULL
    } else {
      ms2PlotValues$fragmentListHovered <<- resultObj
    }
    
    #output$information <- renderText({
    #  print(paste("update output$information PCA Loadings hover ", precursorIndex, sep = ""))
    #  paste(resultObj$infoText, sep = "")
    #})
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
}
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
obsPCAloadingsBrush <- observeEvent(input$plotPcaLoadings_brush, {
  if(is.null(classToSpectra_class))
    return()
  
  brush <- input$plotPcaLoadings_brush
  
  xInterval <- c(brush$xmin, brush$xmax)
  yInterval <- c(brush$ymin, brush$ymax)
  
  
  indeces <- which(
    pcaDataList$pcaLoadingsX >= brush$xmin & pcaDataList$pcaLoadingsX <= brush$xmax &
      pcaDataList$pcaLoadingsY >= brush$ymin & pcaDataList$pcaLoadingsY <= brush$ymax
  )
  precursorSet <- filterPca$filter[indeces]
  print(paste("observe PCAloadings brush", length(precursorSet)))
  if(length(precursorSet) == 0)
    return()
  
  ## putative metabolite families statistics
  #putativeMetaboliteFamilies <- evaluatePutativeMetaboliteFamiliesOfPrecursorSet(dataList = dataList, precursorSet = precursorSet, classToSpectra_class = classToSpectra_class)
  
  output$information <- renderText({
    print(paste("update output$information", resultObj$infoText))
    paste(
      resultObj$infoText, 
      #paste("\n", paste(putativeMetaboliteFamilies, collapse = "\n"), sep = ""),
      sep = "")
  })
})

observeGroupSet <- observeEvent(input$pcaGroups, {
  print(paste("observe groups change", paste(input$pcaGroups, collapse = "-"), length(input$pcaGroups), length(input$pcaGroups) == 2))
  shinyjs::toggleState("pcaFilter_lfc", length(input$pcaGroups) == 2)
  
  if(FALSE){
    ## update samples
    sampleNames <- dataList$dataColumnsNameFunctionFromGroupNames(grouXXXps = input$pcaGroups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
    updateCheckboxGroupInput(session = session, inputId = "pcaSamples",                         selected = sampleNames)
  }
})
observeSampleSet <- observeEvent(input$pcaSamples, {
  print(paste("observe samples change", paste(input$pcaSamples, collapse = "-"), length(input$pcaSamples)))
  
  groupsFromSamples <- unlist(lapply(X = dataList$grouXXXps, FUN = function(x){
    samplesOfGroups <- dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
    if(any(samplesOfGroups %in% input$pcaSamples))
      return(x)
    else
      return(NULL)
  }))
  
  shinyjs::toggleState(id = "pcaFilter_lfc", condition = length(groupsFromSamples) == 2)
  
  #########################################################################################
  ## update filter
  
  #if(!is.null(filterHca)){
  ##  applyHcaFilters(
  ##    filterHca$groupSetOriginal[[1]], 
  ##    filterHca$groupSetOriginal[[2]], 
  ##    filterHca$filter_averageOriginal, 
  ##    filterHca$filter_lfcOriginal, 
  ##    filterHca$includeIgnoredPrecursorsOriginal
  ##  )
  #  applyPcaFilters_default()
  ##}
  
  if(!is.null(filterPca))
    applyPcaFilters_default()
})

observeSelectAllPCAGroups <- observeEvent(input$selectAllPCAGroups, {
  print(paste("observe selectAllPCAGroups"))
  updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$grouXXXps, selected = dataList$grouXXXps)
})
observeSelectNoPCAGroups <- observeEvent(input$selectNoPCAGroups, {
  print(paste("observe selectNoPCAGroups"))
  updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$grouXXXps, selected = NULL)
})
observeSelectInvertedPCAGroups <- observeEvent(input$selectInvertedPCAGroups, {
  print(paste("observe selectInvertedPCAGroups"))
  updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$grouXXXps, selected = setdiff(dataList$grouXXXps, input$pcaGroups))
})

output$plotPcaScores_hover_info <- renderUI({
  hover <- input$plotPcaScores_hover
  hoverX <- hover$x
  hoverY <- hover$y
  plotWidth  <- session$clientData$output_plotPcaScores_width
  plotHeight <- session$clientData$output_plotPcaScores_height
  
  if(is.null(hoverX) | is.null(hoverY))
    return(NULL)
  #print(paste("UI PCA scores hover", hoverX, hoverY))
  
  output$tip <- renderText({
    print(paste("update output$tip"))
    paste("Hover a scores node to view information about the sample.", "Brush and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
  })
  
  #################################################
  ## decide whether the hover is close enough to trigger event
  minimumIndex <- getSelectedPOI_XY(
    mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaScoresX, poiCoordinatesY = pcaDataList$pcaScoresY,
    plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaScoresPlotRange$xIntervalSize, plotRangeY = pcaScoresPlotRange$yIntervalSize
  )
  print(paste("UI PCA scores hover", hoverX, hoverY, minimumIndex))
  
  if(is.null(minimumIndex))
    return(NULL)
  
  ## compile information
  dataColumnName <- filterPca$sampleSet[[minimumIndex]]
  #dataColumnName <- dataList$dataColumnsNameFunctionFromGroupNames(grouXXXps = filterPca$grouXXXps, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))[[minimumIndex]]
  group <- dataList$groupNameFunctionFromDataColumnName(dataColumnName = dataColumnName, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
  info <- paste(
    "<b>Sample: ", "</b>", dataColumnName, "<br>",
    "<b>Group: ", "</b>", group,
    sep = ""
  )
  #output$information <- renderText({
  #  print(paste("update output$information PCA scores hover", sep = ""))
  #  paste("Sample '", dataColumnName , "' is a replicate of group '", group, "'.", sep = "")
  #})
  
  panelWidth <- as.integer(plotWidth*0.6)
  showPlotTooltip(hover, info, panelWidth)
})
output$plotPcaLoadings_hover_info <- renderUI({
  hover <- input$plotPcaLoadings_hover
  hoverX <- hover$x
  hoverY <- hover$y
  plotWidth  <- session$clientData$output_plotPcaLoadings_width
  plotHeight <- session$clientData$output_plotPcaLoadings_height
  
  if(is.null(hoverX) | is.null(hoverY))
    return(NULL)
  #print(paste("UI PCA loadings hover", hoverX, hoverY))
  
  output$tip <- renderText({
    print(paste("update output$tip"))
    paste(
      "Hover or click a loadings node to view information about the corresponding MS\u00B9 feature.", 
      "Brush and double-click to zoom in.", 
      "Double-click to zoom out.", 
      sep = "\n"
    )
  })
  
  #################################################
  ## decide whether the hover is close enough to trigger event
  minimumIndex <- getSelectedPOI_XY(
    mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaLoadingsX, poiCoordinatesY = pcaDataList$pcaLoadingsY,
    plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaLoadingsPlotRange$xIntervalSize, plotRangeY = pcaLoadingsPlotRange$yIntervalSize
  )
  print(paste("UI PCA loadings hover", hoverX, hoverY, minimumIndex))
  
  if(is.null(minimumIndex)){
    ms2PlotValues$fragmentListHovered <<- NULL
    return(NULL)
  }
  
  ## compile information
  precursorIndex <- filterPca$filter[[minimumIndex]]
  resultObj <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)
  if(all(!is.null(selectionAnalysisPcaLoadingSet), minimumIndex == selectionAnalysisPcaLoadingSet)){
    ## reverse MS2 to clicked stuff
    ms2PlotValues$fragmentListHovered <<- NULL
  } else {
    ms2PlotValues$fragmentListHovered <<- resultObj
  }
  
  info <- paste(
    "<b>MS\u00B9 feature: ", "</b>", resultObj$infoFeatureLabel, "<br>",
    "<b>Name: ", "</b>", resultObj$infoName, "<br>",
    "<b>Annotation: ", "</b>", resultObj$infoFamilies, "<br>",
    "<b>Fragments: ", "</b>", resultObj$infoFragmentCount,
    sep = ""
  )
  
  panelWidth <- as.integer(plotWidth*0.6)
  showPlotTooltip(hover, info, panelWidth)
})

output$ui_plotAnnoLegendPCA <- renderUI({
  print(paste("### GUI ### ui_plotAnnoLegendPCA"))
  if(state_tabPca$annotationLegendHeightPca != -1) {
    plotOutput(outputId = "plotAnnoLegendPCA", height = state_tabPca$annotationLegendHeightPca)
  }
})
output$ui_plotScoresGroupsLegend <- renderUI({
  print(paste("### GUI ### ui_plotScoresGroupsLegend"))
  if(state_tabPca$scoresGroupsLegendHeight != -1) {
    plotOutput(outputId = "plotScoresGroupsLegend", height = state_tabPca$scoresGroupsLegendHeight)
  }
})

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabPca observers")
  obsDrawPCA$suspend()
  obsPCAscoresDblClick$suspend()
  obsPCAloadingsClick$suspend()
  obsPCAloadingsDblClick$suspend()
  obsPCAloadingsBrush$suspend()
  observeGroupSet$suspend()
  observeSampleSet$suspend()
  observeSelectAllPCAGroups$suspend()
  observeSelectNoPCAGroups$suspend()
  observeSelectInvertedPCAGroups$suspend()
})
