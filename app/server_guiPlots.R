
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
    hcaPrecursorLabels = state$hcaPrecursorLabels, 
    xInterval = dendrogramPlotRange$xInterval
  )
  
  dendrogramUserCoordinateRange <- par("usr")
  dendrogramUserCoordinateRangeY <<- dendrogramUserCoordinateRange[[4]] - dendrogramUserCoordinateRange[[3]]
  
  state$annotationsHca <<- resultObj
  state$annotationLegendHeightHca <<- annoLegendEntryHeight * (length(state$annotationsHca$setOfAnnotations) + 1)
}
drawHeatmapPlotImpl <- function(consoleInfo = NULL){
  #calcPlotHeatmap(dataList = dataList, filterObj = filterHca, clusterDataList = clusterDataList, xInterval = dendrogramPlotRange$xInterval)
  if(!is.null(consoleInfo)) print(paste("### hea ###", consoleInfo))
  
  ## selections and selection colors
  selectedTreeNodeSet <- NULL
  frameColor <- NULL
  
  selectedSelection <- state$selectedSelection
  if(!is.null(selectedSelection)){
    analysisSelection <- function(){
      selectedTreeNodeSet <- selectionAnalysisTreeNodeSet
      frameColor <- "blue"
    }
    fragmentSelection <- function(){
      selectedTreeNodeSet <- selectionFragmentTreeNodeSet
      frameColor <- "green"
    }
    searchSelection <- function(){
      selectedTreeNodeSet <- selectionSearchTreeNodeSet
      frameColor <- "red"
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
    heatmapContent  <- state$heatmapContent
    heatmapOrdering <- state$heatmapOrdering
    
    returnObj <- calcPlotHeatmap(
      dataList = dataList, 
      filterObj = filterHca, 
      clusterDataList = clusterDataList, 
      selectedTreeNodeSet = selectedTreeNodeSet, 
      frameColor = frameColor,
      heatmapContent = heatmapContent,
      heatmapOrdering = heatmapOrdering,
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
    
    state$heatmapHeight <<- heatmapHeightPerRow * length(columnsOfInterest)
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
    state$heatmapHeight <<- heatmapHeightPerRow * 3
  }
}
drawHeatmapLegendImpl <- function(){
  calcPlotHeatmapLegend(dataList = dataList)
}
drawDendrogramLegendImpl <- function(){
  calcPlotDendrogramLegend()
}
drawAnnotationLegendHCAimpl <- function(){
  calcPlotAnnoLegend(state$annotationsHca$setOfAnnotations, state$annotationsHca$setOfColors)
}

drawPcaScoresPlotImpl <- function(){
  #resultObj <- calcPlotPCAscores(
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
    loadingsLabels = state$loadingsLabels, 
    showLoadingsAbundance = state$showLoadingsAbundance, 
    showLoadingsFeaturesAnnotated   = state$showLoadingsFeaturesAnnotated,
    showLoadingsFeaturesUnannotated = state$showLoadingsFeaturesUnannotated,
    showLoadingsFeaturesSelected    = state$showLoadingsFeaturesSelected,
    showLoadingsFeaturesUnselected  = state$showLoadingsFeaturesUnselected
  )
  
  state$annotationsPca <<- resultObj
  state$annotationLegendHeightPca <<- annoLegendEntryHeight * (length(state$annotationsPca$setOfAnnotations) + 1)
}
drawAnnotationLegendPCAimpl <- function(){
  calcPlotAnnoLegend(state$annotationsPca$setOfAnnotations, state$annotationsPca$setOfColors)
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
                  "left:", left_px + 2, "px; top:", top_px + 2, "px; width:", panelWidth, "px;")
  
  # actual tooltip created as wellPanel
  wellPanel(
    style = style,
    p(HTML(info))
  )
}
