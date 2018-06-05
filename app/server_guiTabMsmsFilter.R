
fragmentPlotRange    <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)

obsFragmentPlotdblClick <- observeEvent(input$fragmentPlot_dblclick, {
  brush <- input$fragmentPlot_brush
  print(paste("observe fragmentPlot dblclick", paste(brush, collapse = "; ")))
  
  if (!is.null(brush)) {
    ## set range
    min <- brush$xmin
    max <- brush$xmax
  } else {
    ## reset range
    min <- min(dataList$ms2_masses)
    max <- max(dataList$ms2_masses)
  }
  
  fragmentPlotRange$xMin <<- min
  fragmentPlotRange$xMax <<- max
  fragmentPlotRange$xInterval <<- c(min, max)
  fragmentPlotRange$xIntervalSize <<- max - min
})
obsApplyGlobalMS2filters <- observeEvent(input$applyGlobalMS2filters, {
  session$sendCustomMessage("disableButton", "applyGlobalMS2filters")
  applyGlobalMS2filters <- as.numeric(input$applyGlobalMS2filters)
  
  print(paste("Observe applyGlobalMS2filters", applyGlobalMS2filters))
  
  #################################################
  ## check if button was hit
  #if(applyGlobalMS2filters == applyGlobalMS2filtersButtonValue)
  #  return()
  #applyGlobalMS2filtersButtonValue <<- applyGlobalMS2filters
  
  #################################################
  ## get inputs
  filter_ms2_masses1  <- input$globalFilter_ms2_masses1
  filter_ms2_masses2  <- input$globalFilter_ms2_masses2
  filter_ms2_masses3  <- input$globalFilter_ms2_masses3
  filter_ms2_ppm      <- input$globalFilter_ms2_ppm
  
  applyGlobalMS2filters(filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm)
  session$sendCustomMessage("enableButton", "applyGlobalMS2filters")
})
obsClearGlobalMS2filters <- observeEvent(input$clearGlobalMS2filters, {
  session$sendCustomMessage("disableButton", "clearGlobalMS2filters")
  clearGlobalMS2filters <- as.numeric(input$clearGlobalMS2filters)
  
  print(paste("Observe clearGlobalMS2filters", clearGlobalMS2filters))
  
  #################################################
  ## check if button was hit
  #if(clearGlobalMS2filters == clearGlobalMS2filtersButtonValue)
  #  return()
  #clearGlobalMS2filtersButtonValue <<- clearGlobalMS2filters
  
  #################################################
  ## get inputs
  filter_ms2_masses1  <- ""
  filter_ms2_masses2  <- ""
  filter_ms2_masses3  <- ""
  filter_ms2_ppm      <- 20
  
  applyGlobalMS2filters(filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm)
  session$sendCustomMessage("enableButton", "clearGlobalMS2filters")
})

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabMsmsFilter observers")
  obsFragmentPlotdblClick$suspend()
  obsApplyGlobalMS2filters$suspend()
  obsClearGlobalMS2filters$suspend()
})
