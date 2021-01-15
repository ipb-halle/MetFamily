
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset plot state")
  ## reset plots
  doClearPlots()
})

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
  ## see also nearPoints(...): http://shiny.rstudio.com/gallery/plot-interaction-selecting-points.html
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
  if(any(length(minimumIndex)==0, is.na(minimumIndex)))
    return(NULL)
  
  minimumDistance <- distances[[minimumIndex]]
  if(minimumDistance > distanceThreshold){
    return(NULL)
  } else {
    return(minimumIndex)
  }
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
    showClusterLabels = state_tabHca$showClusterLabels, 
    hcaPrecursorLabels = state_tabHca$hcaPrecursorLabels, 
    xInterval = dendrogramPlotRange$xInterval
  )
  
  dendrogramUserCoordinateRange <- par("usr")
  dendrogramUserCoordinateRangeY <<- dendrogramUserCoordinateRange[[4]] - dendrogramUserCoordinateRange[[3]]
  
  state_tabHca$annotationsHca <<- resultObj
  state_tabHca$annotationLegendHeightHca <<- annoLegendEntryHeight * (length(state_tabHca$annotationsHca$setOfAnnotations) + 1)
}
drawHeatmapPlotImpl <- function(consoleInfo = NULL){
  #calcPlotHeatmap(dataList = dataList, filterObj = filterHca, clusterDataList = clusterDataList, xInterval = dendrogramPlotRange$xInterval)
  if(!is.null(consoleInfo)) print(paste("### hea ###", consoleInfo))
  
  ## selections and selection colors
  selectedTreeNodeSet <- NULL
  frameColor <- NULL
  
  selectedSelection <- state_selections$selectedSelection
  if(!is.null(selectedSelection)){
    analysisSelection <- function(){
      selectedTreeNodeSet <<- selectionAnalysisTreeNodeSet
      frameColor <<- "blue"
    }
    fragmentSelection <- function(){
      selectedTreeNodeSet <<- selectionFragmentTreeNodeSet
      frameColor <<- "green"
    }
    searchSelection <- function(){
      selectedTreeNodeSet <<- selectionSearchTreeNodeSet
      frameColor <<- "red"
    }
    switch(selectedSelection,
           "Analysis_HCA"=analysisSelection(),
           "Analysis_PCA"=analysisSelection(),
           "Fragment_HCA"=fragmentSelection(),
           "Fragment_PCA"=fragmentSelection(),
           "Search_HCA"  =searchSelection(),
           "Search_PCA"  =searchSelection(),
           {## unknown state
             stop(paste("Unknown selectedSelection value", selectedSelection))
           }
    )## end switch
  }
  
  if(hcaHeatMapNew){
    #####################################################################################
    ## new heatmap functionality
    returnObj <- calcPlotHeatmap(
      dataList = dataList, 
      filterObj = filterHca, 
      clusterDataList = clusterDataList, 
      selectedTreeNodeSet = selectedTreeNodeSet, 
      frameColor = frameColor,
      heatmapContent = state_tabHca$heatmapContent,
      heatmapOrdering = state_tabHca$heatmapOrdering,
      xInterval = dendrogramPlotRange$xInterval
    )
    columnsOfInterest <- returnObj$columnsOfInterest
    columnsOfInterestForHeatmap <<- columnsOfInterest
    
    ## heigth per row
    heatmapHeightPerRow <- -1
    if(length(columnsOfInterest) <= maximumheatmapHeightRowCount){
      heatmapHeightPerRow <- maximumheatmapHeightPerRow
    } else if(length(columnsOfInterest) >= minimumheatmapHeightRowCount){
      heatmapHeightPerRow <- minimumheatmapHeightPerRow
    } else {
      heatmapHeightPerRow <- 
        minimumheatmapHeightPerRow + 
        (length(columnsOfInterest)    - maximumheatmapHeightRowCount) / 
        (minimumheatmapHeightRowCount - maximumheatmapHeightRowCount) * 
        (maximumheatmapHeightPerRow - minimumheatmapHeightPerRow)
    }
    
    state_tabHca$heatmapHeight <<- heatmapHeightPerRow * length(columnsOfInterest)
  } else {
    #####################################################################################
    ## old heatmap functionality
    columnsOfInterest <- calcPlotHeatmapOld(
      dataList = dataList, 
      filterObj = filterHca, 
      clusterDataList = clusterDataList, 
      xInterval = dendrogramPlotRange$xInterval
    )
    columnsOfInterestForHeatmap <<- columnsOfInterest
    state_tabHca$heatmapHeight <<- heatmapHeightPerRow * 3
  }
}
drawHeatmapLegendImpl <- function(){
  calcPlotHeatmapLegend(dataList = dataList)
}
drawDendrogramLegendImpl <- function(){
  calcPlotDendrogramLegend()
}
drawAnnotationLegendHCAimpl <- function(){
  calcPlotAnnoLegend(state_tabHca$annotationsHca$setOfAnnotations, state_tabHca$annotationsHca$setOfColors)
}
drawAnnotationLegendForImageHCAimpl <- function(){
  calcPlotAnnoLegendForImage(state_tabHca$annotationsHca$setOfAnnotations, state_tabHca$annotationsHca$setOfColors)
}

drawPcaScoresPlotImpl <- function(){
  #resultObj <- calcPlotPCAscores(
  resultObj <- calcPlotPCAscores(
    pcaObj = pcaDataList$pcaObj, 
    dataList = dataList, 
    filterObj = filterPca, 
    pcaDimensionOne = pcaDataList$dimensionOne, 
    pcaDimensionTwo = pcaDataList$dimensionTwo, 
    showScoresLabels = state_tabPca$showScoresLabels, 
    xInterval = pcaScoresPlotRange$xInterval, 
    yInterval = pcaScoresPlotRange$yInterval
  )
}
#### I am adding this new 
drawPcaScoresPlotImpl1 <- function(){
  #resultObj <- calcPlotPCAscores(
  resultObj <- calcPlotPCAscores1(
    pcaObj = pcaDataList$pcaObj, 
    dataList = dataList, 
    filterObj = filterPca, 
    pcaDimensionOne = pcaDataList$dimensionOne, 
    pcaDimensionTwo = pcaDataList$dimensionTwo, 
    showScoresLabels = state_tabPca$showScoresLabels, 
    xInterval = pcaScoresPlotRange$xInterval, 
    yInterval = pcaScoresPlotRange$yInterval
  )
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
    selectionSearchPcaLoadingSet   = selectionSearchPcaLoadingSet,
    xInterval = pcaLoadingsPlotRange$xInterval, 
    yInterval = pcaLoadingsPlotRange$yInterval,
    loadingsLabels = state_tabPca$loadingsLabels, 
    showLoadingsAbundance = state_tabPca$showLoadingsAbundance, 
    showLoadingsFeaturesAnnotated   = state_tabPca$showLoadingsFeaturesAnnotated,
    showLoadingsFeaturesUnannotated = state_tabPca$showLoadingsFeaturesUnannotated,
    showLoadingsFeaturesSelected    = state_tabPca$showLoadingsFeaturesSelected,
    showLoadingsFeaturesUnselected  = state_tabPca$showLoadingsFeaturesUnselected
  )
  
  state_tabPca$annotationsPca <<- resultObj
  state_tabPca$annotationLegendHeightPca <<- annoLegendEntryHeight * (length(state_tabPca$annotationsPca$setOfAnnotations) + 1)
}
drawAnnotationLegendPCAimpl <- function(){
  calcPlotAnnoLegend(state_tabPca$annotationsPca$setOfAnnotations, state_tabPca$annotationsPca$setOfColors)
}
drawAnnotationLegendForImagePCAimpl <- function(){
  #### testing this with 1
  calcPlotAnnoLegendForImage1(state_tabPca$annotationsPca$setOfAnnotations, state_tabPca$annotationsPca$setOfColors)
}
drawScoresGroupsLegendImpl <- function(){
  calcPlotScoresGroupsLegend(scoresGroups$groups, scoresGroups$colors)
}

drawMS2PlotImpl <- function(){
  calcPlotMS2(
    dataList = dataList, 
    fragmentsX = ms2PlotValues$fragmentListClicked$fragmentMasses, 
    fragmentsY = ms2PlotValues$fragmentListClicked$fragmentAbundances, 
    fragmentsColor = ms2PlotValues$fragmentListClicked$fragmentColor, 
    fragmentsDiscriminativity = ms2PlotValues$fragmentListClicked$fragmentDiscriminativity, 
    fragmentsX_02 = ms2PlotValues$fragmentListHovered$fragmentMasses, 
    fragmentsY_02 = ms2PlotValues$fragmentListHovered$fragmentAbundances, 
    fragmentsColor_02 = ms2PlotValues$fragmentListHovered$fragmentColor, 
    xInterval = ms2PlotRange$xInterval, 
    selectedFragmentIndex = selectionFragmentSelectedFragmentIndex,
    dendrogramFragmentStatistics = dendrogramFragmentStatistics
  )
}
drawMS2LegendImpl <- function(){
  calcPlotMS2Legend(dataList = dataList)
}
drawFragmentDiscriminativityLegendImpl <- function(){
  calcPlotDiscriminativityLegend()
}

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

## taken from https://gitlab.com/snippets/16220
## or https://github.com/OHDSI/StudyProtocols/blob/master/LargeScalePopEst/extras/ShinyApp/server.R
## or https://stackoverflow.com/questions/48044543/custom-tooltip-in-ggplot-using-geom-polygoncoord-map
showPlotTooltip <- function(hover, info, panelWidth){
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right  - hover$range$left)
  top_px  <- hover$range$top  + top_pct  * (hover$range$bottom - hover$range$top)
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  #style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
  style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                  "left:", left_px + 7, "px; top:", top_px + 7, "px; width:", panelWidth, "px;")
  
  # actual tooltip created as wellPanel
  wellPanel(
    style = style,
    p(HTML(info))
  )
}
