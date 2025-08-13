
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset tabSearch state")
  ## input fields: MS1 search
  updateTextInput(session = session, inputId = "searchMS1mass", value = "")
  updateTextInput(session = session, inputId = "searchMS1massPpm", value = 20)
  ## input fields: MS2 search
  updateTextInput(session = session, inputId = "search_ms2_masses1", value = "")
  updateTextInput(session = session, inputId = "search_ms2_masses2", value = "")
  updateTextInput(session = session, inputId = "search_ms2_masses3", value = "")
  updateTextInput(session = session, inputId = "searchMS2massPpm", value = 20)
  updateCheckboxInput(session = session, inputId = "searchIncludeIgnoredPrecursors", value = FALSE)
})

obsClearSearch <- observeEvent(input$clearSearch, {
  session$sendCustomMessage("disableButton", "clearSearch")
  clearSearch <- as.numeric(input$clearSearch)
  
  print(paste("Observe clearSearch", clearSearch))
  
  #################################################
  ## check if button was hit
  #if(clearSearch == clearSearchButtonValue)
  #  return()
  #clearSearchButtonValue <<- clearSearch
  
  filterSearch <<- NULL
  state_filters$filterSearchActive <<- FALSE
  state_filters$searchFilterValid <<- TRUE
  
  selectionBySearchReset()
  updateSearchInformation()
  
  #################################################
  ## update plots
  if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "clear search", withHeatmap = TRUE)
  if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "clear search")
  session$sendCustomMessage("enableButton", "clearSearch")
})
obsApplySearch <- observeEvent(input$applySearch, {
  session$sendCustomMessage("disableButton", "applySearch")
  applySearch <- as.numeric(input$applySearch)
  
  print(paste("Observe applySearch", applySearch))
  
  #################################################
  ## check if button was hit
  #if(applySearch == applySearchButtonValue)
  #  return()
  #applySearchButtonValue <<- applySearch
  
  #################################################
  ## MS1 or MS2?
  includeIgnoredPrecursors  <- input$searchIncludeIgnoredPrecursors
  
  searchMode <- input$searchMS1orMS2
  switch (searchMode,
          'MS1 feature m/z' = {
            filter_ms1_masses <- input$searchMS1mass
            filter_ms1_ppm    <- input$searchMS1massPpm
            
            if(nchar(trimws(filter_ms1_masses)) == 0)
              return()
            
            resultObj <- doApplySearch(filter_ms1_masses = filter_ms1_masses, filter_ms1_ppm = filter_ms1_ppm, includeIgnoredPrecursors = includeIgnoredPrecursors)
          },
          'Fragment m/z' = {
            filter_ms2_masses1 <- input$search_ms2_masses1
            filter_ms2_masses2 <- input$search_ms2_masses2
            filter_ms2_masses3 <- input$search_ms2_masses3
            filter_ms2_ppm     <- input$searchMS2massPpm
            
            resultObj <- doApplySearch(filter_ms2_masses1 = filter_ms2_masses1, filter_ms2_masses2 = filter_ms2_masses2, filter_ms2_masses3 = filter_ms2_masses3, filter_ms2_ppm = filter_ms2_ppm, includeIgnoredPrecursors = includeIgnoredPrecursors)
          }
  )
  
  processSearchFilterResult(resultObj)
  
  session$sendCustomMessage("enableButton", "applySearch")
})
doApplySearch <- function(filter_ms2_masses1 = NULL, filter_ms2_masses2 = NULL, filter_ms2_masses3 = NULL, filter_ms2_ppm = NULL, filter_ms1_masses = NULL, filter_ms1_ppm = NULL, includeIgnoredPrecursors){
  filter_lfc      <- NULL
  filter_average  <- NULL
  groupSet        <- dataList$sampleClasses
  
  #################################################
  ## do filtering
  sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterBySamples <- TRUE
  
  print(paste("Observe applySearch", "1m", filter_ms1_masses, "1p", filter_ms1_ppm, "i", includeIgnoredPrecursors, "gs", paste(groupSet, collapse = "-")))
  resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
  return(resultObj)
}
processSearchFilterResult <- function(resultObj){
  state_filters$filterSearchActive <<- TRUE
  #################################################
  ## info / error
  if(resultObj$error){
    filterSearch <<- NULL
    state_filters$searchFilterValid <<- FALSE
    selectionBySearch(NULL)
    updateSearchInformation()
  } else {
    filterSearch <<- resultObj$filter
    state_filters$searchFilterValid <<- TRUE
    selectionBySearch(filterSearch$filter)
    updateSearchInformation()
    
    #################################################
    ## update plots
    if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "update search", withHeatmap = TRUE)
    if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "update search")
  }
}
updateSearchInformation <- function(){
  if(!state_filters$filterSearchActive)
    output$searchInfo <- renderText({
      print(paste("update output$searchInfo inactive search", sep = ""))
      paste("Please search for MS\u00B9 features", sep = "")
    })
  if(state_filters$filterSearchActive & is.null(filterSearch))
    output$searchInfo <- renderText({
      print(paste("update output$searchInfo invalid search", sep = ""))
      paste("There are invalid or missing search values", sep = "")
    })
  if(state_filters$filterSearchActive & !is.null(filterSearch)){
    str1 <- ""
    str2 <- ""
    if(state$showHCAplotPanel & !is.null(listForTable_Search_HCA))
      str1 <- paste(length(listForTable_Search_HCA$precursorSet), " in HCA", sep = "")
    if(state$showPCAplotPanel & !is.null(listForTable_Search_PCA))
      str2 <- paste(length(listForTable_Search_PCA$precursorSet), " in PCA", sep = "")
    
    if(nchar(str1) > 0 & nchar(str2) > 0)
      val <- paste(str1, str2, sep = ", ")
    if(nchar(str1) > 0 & !(nchar(str2) > 0))
      val <- str1
    if(!(nchar(str1) > 0) & nchar(str2) > 0)
      val <- str2
    if(!(nchar(str1) > 0) & !(nchar(str2) > 0))
      val <- "None"
    
    output$searchInfo <- renderText({
      print(paste("update output$searchInfo", sep = ""))
      paste("Number of hits among MS\u00B9 features: ", val, sep = "")
    })
  }
}

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabSearch observers")
  obsClearSearch$suspend()
  obsApplySearch$suspend()
})
