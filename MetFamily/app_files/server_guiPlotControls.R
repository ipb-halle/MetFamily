
#TODO need this?
obsShowPlotControls <- observeEvent(input$showPlotControls, {
  showPlotControls <- input$showPlotControls
  print(paste("Observe showPlotControls", showPlotControls))
  state$showPlotControls <<- showPlotControls
})
obsShowClusterLabels <- observeEvent(input$showClusterLabels, {
  showClusterLabels <- input$showClusterLabels
  print(paste("Observe showClusterLabels", showClusterLabels))
  state_tabHca$showClusterLabels <<- showClusterLabels
  #drawDendrogramPlot(consoleInfo = "showClusterLabels")
})
obsHeatmapContent <- observeEvent(input$heatmapContent, {
  heatmapContent <- input$heatmapContent
  
  if(is.null(dataList))
    return()
  
  print(paste("Observe heatmapContent", heatmapContent))
  state_tabHca$heatmapContent <<- heatmapContent
})
obsHeatmapOrdering <- observeEvent(input$heatmapOrdering, {
  if(is.null(dataList))
    return()
  
  heatmapOrdering <- input$heatmapOrdering
  print(paste("Observe heatmapOrdering", heatmapOrdering))
  state_tabHca$heatmapOrdering <<- heatmapOrdering
})
obsHcaPrecursorLabels <- observeEvent(input$hcaPrecursorLabels, {
  hcaPrecursorLabels <- input$hcaPrecursorLabels
  print(paste("Observe hcaPrecursorLabels", hcaPrecursorLabels))
  state_tabHca$hcaPrecursorLabels <<- hcaPrecursorLabels
})
obsShowScoresLabels <- observeEvent(input$showScoresLabels, {
  showScoresLabels <- input$showScoresLabels
  print(paste("Observe showScoresLabels", showScoresLabels))
  state_tabPca$showScoresLabels <<- showScoresLabels
})
obsLoadingsLabels <- observeEvent(input$loadingsLabels, {
  loadingsLabels <- input$loadingsLabels
  print(paste("Observe loadingsLabels", loadingsLabels))
  state_tabPca$loadingsLabels <<- loadingsLabels
})
obsShowLoadingsFeatures <- observeEvent(input$showLoadingsFeatures, {
  showLoadingsFeatures <- input$showLoadingsFeatures
  print(paste("Observe showLoadingsFeatures", paste(showLoadingsFeatures, collapse = ";")))
  
  {
    state_tabPca$showLoadingsFeaturesAnnotated   <<- "Annotated"     %in% showLoadingsFeatures
    state_tabPca$showLoadingsFeaturesUnannotated <<- "Not Annotated" %in% showLoadingsFeatures
    state_tabPca$showLoadingsFeaturesSelected    <<- "Selected"      %in% showLoadingsFeatures
    state_tabPca$showLoadingsFeaturesUnselected  <<- "Not Selected"  %in% showLoadingsFeatures
  }
})
obsShowLoadingsAbundance <- observeEvent(input$showLoadingsAbundance, {
  showLoadingsAbundance <- input$showLoadingsAbundance
  print(paste("Observe showLoadingsAbundance", showLoadingsAbundance))
  state_tabPca$showLoadingsAbundance <<- showLoadingsAbundance
})

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending plotControls observers")
  obsShowPlotControls$suspend()
  obsShowClusterLabels$suspend()
  obsHeatmapContent$suspend()
  obsHeatmapOrdering$suspend()
  obsHcaPrecursorLabels$suspend()
  observeGroupSet$suspend()
  obsShowScoresLabels$suspend()
  obsLoadingsLabels$suspend()
  obsShowLoadingsFeatures$suspend()
  obsShowLoadingsAbundance$suspend()
})
