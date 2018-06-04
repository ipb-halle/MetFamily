
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
  state$filterSearchActive <<- FALSE
  state$searchfilterValid <<- TRUE
  
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
  searchMode <- input$searchMS1orMS2
  if(searchMode == 'MS1 feature m/z'){
    #################################################
    ## get inputs
    filter_ms1_masses <- input$searchMS1mass
    filter_ms1_ppm    <- input$searchMS1massPpm
    
    if(nchar(trimws(filter_ms1_masses)) == 0)
      return()
    
    filter_ms2_masses1  <- NULL
    filter_ms2_masses2  <- NULL
    filter_ms2_masses3  <- NULL
    filter_ms2_ppm      <- NULL
  }
  if(searchMode == 'Fragment m/z'){
    #################################################
    ## get inputs
    filter_ms2_masses1 <- input$search_ms2_masses1
    filter_ms2_masses2 <- input$search_ms2_masses2
    filter_ms2_masses3 <- input$search_ms2_masses3
    filter_ms2_ppm     <- input$searchMS2massPpm
    
    filter_ms1_masses <- NULL
    filter_ms1_ppm    <- NULL
  }
  
  filter_lfc      <- NULL
  filter_average  <- NULL
  groupSet        <- dataList$groups
  includeIgnoredPrecursors  <- input$searchIncludeIgnoredPrecursors
  
  #################################################
  ## do filtering
  sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterBySamples <- TRUE
  
  print(paste("Observe applySearch", "1m", filter_ms1_masses, "1p", filter_ms1_ppm, "i", includeIgnoredPrecursors, "gs", paste(groupSet, collapse = "-")))
  resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
  processSearchFilterResult(resultObj)
  session$sendCustomMessage("enableButton", "applySearch")
})
processSearchFilterResult <- function(resultObj){
  state$filterSearchActive <<- TRUE
  #################################################
  ## info / error
  if(resultObj$error){
    filterSearch <<- NULL
    state$searchfilterValid <<- FALSE
    selectionBySearch(NULL)
    updateSearchInformation()
  } else {
    filterSearch <<- resultObj$filter
    state$searchfilterValid <<- TRUE
    selectionBySearch(filterSearch$filter)
    updateSearchInformation()
    
    #################################################
    ## update plots
    if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "update search", withHeatmap = TRUE)
    if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "update search")
  }
}
updateSearchInformation <- function(){
  if(!state$filterSearchActive)
    output$searchInfo <- renderText({
      print(paste("update output$searchInfo inactive search", sep = ""))
      paste("Please search for MS\u00B9 features", sep = "")
    })
  if(state$filterSearchActive & is.null(filterSearch))
    output$searchInfo <- renderText({
      print(paste("update output$searchInfo invalid search", sep = ""))
      paste("There are invalid or missing search values", sep = "")
    })
  if(state$filterSearchActive & !is.null(filterSearch)){
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
selectionBySearchReset <- function(){
  if(!is.null(selectionSearchTreeNodeSet)){ ## HCA
    selectionSearchTreeNodeSet <<- NULL
    listForTable_Search_HCA <<- NULL
    table_Search_HCA_id <<- NULL
    table$df_Search_HCA <<- NULL
    if(state$selectedSelection == selectionSearchHcaName)
      updateSelectedPrecursorSet()
  }
  if(!is.null(selectionSearchPcaLoadingSet)){ ## HCA
    selectionSearchPcaLoadingSet <<- NULL
    listForTable_Search_PCA <<- NULL
    table_Search_PCA_id <<- NULL
    table$df_Search_PCA <<- NULL
    if(state$selectedSelection == selectionSearchHcaName)
      updateSelectedPrecursorSet()
  }
}
selectionBySearch <- function(precursorSet){
  if(state$showHCAplotPanel)
    selectionBySearchInitHca(precursorSet)
  if(state$showPCAplotPanel)
    selectionBySearchInitPca(precursorSet)
  
  if(input$changeSelection != selectionSearchName){
    updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionSearchName)
    updateSelectedSelection()
  }
}
