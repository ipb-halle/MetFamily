
#https://stat.ethz.ch/R-manual/R-devel/library/utils/html/zip.html
serialization <- function(){
  #######################################
  ## enlist
  paramsList <- list(
    ## global MS2 filter
    globalFilter_ms2_masses1          = input$globalFilter_ms2_masses1,
    globalFilter_ms2_masses2          = input$globalFilter_ms2_masses2,
    globalFilter_ms2_masses3          = input$globalFilter_ms2_masses3,
    globalFilter_ms2_ppm              = input$globalFilter_ms2_ppm,
    #input$applyGlobalMS2filters
    ## HCA
    hcaFilterGroupOne                 = input$hcaFilterGroupOne,
    hcaFilterGroupTwo                 = input$hcaFilterGroupTwo,
    hcaFilter_average                 = input$hcaFilter_average,
    hcaFilter_lfc                     = input$hcaFilter_lfc,
    hcaFilterIncludeIgnoredPrecursors = input$hcaFilterIncludeIgnoredPrecursors,
    #input$applyHcaFilters
    hcaDistanceFunction               = input$hcaDistanceFunction,
    #hcaClusterMethod                  = input$hcaClusterMethod,
    hcaClusterMethod                  = "ward.D",
    #input$drawHCAplots
    ## PCA
    pcaGroups                         = input$pcaGroups,
    pcaSamples                        = input$pcaSamples,
    pcaFilter_average                 = input$pcaFilter_average,
    pcaFilter_lfc                     = input$pcaFilter_lfc,
    pcaFilterIncludeIgnoredPrecursors = input$pcaFilterIncludeIgnoredPrecursors,
    #input$applyPcaFilters
    pcaScaling                        = input$pcaScaling,
    pcaLogTransform                   = input$pcaLogTransform,
    pcaDimensionOne                   = input$pcaDimensionOne,
    pcaDimensionTwo                   = input$pcaDimensionTwo,
    #input$drawPCAplots
    ## plot properties
    showPlotControls                  = input$showPlotControls,
    showClusterLabels                 = input$showClusterLabels,
    heatmapContent                    = input$heatmapContent,
    heatmapOrdering                   = input$heatmapOrdering,
    hcaPrecursorLabels                = input$hcaPrecursorLabels,
    showScoresLabels                  = input$showScoresLabels,
    loadingsLabels                    = input$loadingsLabels,
    showLoadingsFeatures              = input$showLoadingsFeatures,
    #showLoadingsFeaturesAnnotated     = input$showLoadingsFeaturesAnnotated,
    #showLoadingsFeaturesUnannotated   = input$showLoadingsFeaturesUnannotated,
    #showLoadingsFeaturesSelected      = input$showLoadingsFeaturesSelected,
    #showLoadingsFeaturesUnselected    = input$showLoadingsFeaturesUnselected,
    showLoadingsAbundance             = input$showLoadingsAbundance,
    #showLoadingsLabels                = "Show labels"    %in% input$pcaLoadingsProperties,
    #showLoadingsAbundance             = "Show abundance" %in% input$pcaLoadingsProperties,
    ## search
    searchMS1orMS2                    = input$searchMS1orMS2,
    searchMS1mass                     = input$searchMS1mass,
    searchMS1massPpm                  = input$searchMS1massPpm,
    search_ms2_masses1                = input$search_ms2_masses1,
    search_ms2_masses2                = input$search_ms2_masses2,
    search_ms2_masses3                = input$search_ms2_masses3,
    searchMS2massPpm                  = input$searchMS2massPpm,
    searchIncludeIgnoredPrecursors    = input$searchIncludeIgnoredPrecursors
  )
  
  #######################################
  ## serialize parameter list
  
  ## built matrix with param names and param values
  tempMatrix    <- matrix(data = c(names(unlist(paramsList)),paste("'", unlist(paramsList), "'", sep = "")), ncol = 2)
  ## quotation of strings
  #textEntries <- is.na(as.numeric(tempMatrix[, 2]))
  #tempMatrix[textEntries, 2] <- paste("'", tempMatrix[textEntries, 2], "'", sep = "")
  ## param name = param value
  paramStrings  <- apply(tempMatrix, MARGIN = 1, FUN = function(x) {paste(x, collapse = "=")})
  ## collapse
  serialization <- paste(paramStrings, collapse = ";")
  return(serialization)
}
## https://github.com/daattali/advanced-shiny/blob/master/update-input/app.R
deserialization <- function(serialization){
  #######################################
  ## deserialize parameters
  
  #paramStrings <- strsplit(x = strsplit(x = serialization, split = ";")[[1]], split = "=")
  #for(i in 1:length(paramStrings))
  #  paramStrings[[i]] <- paste(paramStrings[[i]][[1]], paste("'", paramStrings[[i]][[2]], "'", sep = ""), sep = "=")
  #paramString <- paste(paramStrings, collapse = ",")
  paramString <- paste(strsplit(x = serialization, split = ";")[[1]], collapse = ",")
  parseText <- paste("paramsList <- list(", paramString, ")", sep = "")
  
  #parseText <- paste("paramsList <- list(", paste(paramStrings, collapse = ","), ")", sep = "")
  eval(parse(text = parseText))
  
  #######################################
  ## update parameter fields
  
  #updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
  #updateRadioButtons(session = session, inputId = "hcaFilterGroupOne", choices = dataList$grouXXXps, selected = selectedOne)
  #updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
  #updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[init]"), selected = lalala)
  #updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$grouXXXps, selected = dataList$grouXXXps)
  
  ## global MS2 filter
  updateTextInput(         session = session, inputId = "globalFilter_ms2_masses1",          value = paramsList$globalFilter_ms2_masses1)
  updateTextInput(         session = session, inputId = "globalFilter_ms2_masses2",          value = paramsList$globalFilter_ms2_masses2)
  updateTextInput(         session = session, inputId = "globalFilter_ms2_masses3",          value = paramsList$globalFilter_ms2_masses3)
  updateTextInput(         session = session, inputId = "globalFilter_ms2_ppm",              value = paramsList$globalFilter_ms2_ppm)
  #input$applyGlobalMS2filters
  ## HCA
  updateRadioButtons(      session = session, inputId = "hcaFilterGroupOne",                 selected = paramsList$hcaFilterGroupOne)
  updateRadioButtons(      session = session, inputId = "hcaFilterGroupTwo",                 selected = paramsList$hcaFilterGroupTwo)
  updateTextInput(         session = session, inputId = "hcaFilter_average",                 value = paramsList$hcaFilter_average)
  updateTextInput(         session = session, inputId = "hcaFilter_lfc",                     value = paramsList$hcaFilter_lfc)
  updateCheckboxInput(     session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = as.logical(paramsList$hcaFilterIncludeIgnoredPrecursors))
  #input$applyHcaFilters
  updateSelectInput(       session = session, inputId = "hcaDistanceFunction",               selected = paramsList$hcaDistanceFunction)
  #updateSelectInput(       session = session, inputId = "hcaClusterMethod",                  selected = paramsList$hcaClusterMethod)
  #input$drawHCAplots
  ## PCA
  updateCheckboxGroupInput(session = session, inputId = "pcaGroups",                         selected = paramsList$pcaGroups)
  updateCheckboxGroupInput(session = session, inputId = "pcaSamples",                        selected = paramsList$pcaSamples)
  updateTextInput(         session = session, inputId = "pcaFilter_average",                 value = paramsList$pcaFilter_average)
  updateTextInput(         session = session, inputId = "pcaFilter_lfc",                     value = paramsList$pcaFilter_lfc)
  updateSelectInput(       session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", selected = paramsList$pcaFilterIncludeIgnoredPrecursors)
  #input$applyPcaFilters
  updateSelectInput(       session = session, inputId = "pcaScaling",                        selected = paramsList$pcaScaling)
  updateCheckboxInput(     session = session, inputId = "pcaLogTransform",                   value = as.logical(paramsList$pcaLogTransform))
  updateSelectInput(       session = session, inputId = "pcaDimensionOne",                   selected = paramsList$pcaDimensionOne)
  updateSelectInput(       session = session, inputId = "pcaDimensionTwo",                   selected = paramsList$pcaDimensionTwo)
  #input$drawPCAplots
  ## plot properties
  updateCheckboxInput(     session = session, inputId = "showPlotControls",                  value = as.logical(paramsList$showPlotControls))
  updateCheckboxInput(     session = session, inputId = "showClusterLabels",                 value = as.logical(paramsList$showClusterLabels))
  updateRadioButtons(      session = session, inputId = "heatmapContent",                    selected = paramsList$heatmapContent)
  updateRadioButtons(      session = session, inputId = "heatmapOrdering",                   selected = paramsList$heatmapOrdering)
  updateRadioButtons(      session = session, inputId = "hcaPrecursorLabels",                selected = paramsList$hcaPrecursorLabels)
  updateCheckboxInput(     session = session, inputId = "showScoresLabels",                  value = as.logical(paramsList$showScoresLabels))
  #updateCheckboxGroupInput(session = session, inputId = "pcaLoadingsProperties",             selected = c(ifelse(as.logical(paramsList$showLoadingsLabels), "Show labels", NULL), ifelse(as.logical(paramsList$showLoadingsAbundance), "Show abundance", NULL)))
  updateRadioButtons(      session = session, inputId = "loadingsLabels",                    selected = paramsList$loadingsLabels)
  updateCheckboxGroupInput(session = session, inputId = "showLoadingsFeatures",              selected = paramsList$showLoadingsFeatures)
  #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesAnnotated",     value = as.logical(paramsList$showLoadingsFeaturesAnnotated))
  #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesUnannotated",   value = as.logical(paramsList$showLoadingsFeaturesUnannotated))
  #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesSelected",      value = as.logical(paramsList$showLoadingsFeaturesSelected))
  #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesUnselected",    value = as.logical(paramsList$showLoadingsFeaturesUnselected))
  updateCheckboxInput(     session = session, inputId = "showLoadingsAbundance",             value = as.logical(paramsList$showLoadingsAbundance))
  ## search
  updateRadioButtons(      session = session, inputId = "searchMS1orMS2",                    selected = paramsList$searchMS1orMS2)
  updateTextInput(         session = session, inputId = "searchMS1mass",                     value = paramsList$searchMS1mass)
  updateTextInput(         session = session, inputId = "searchMS1massPpm",                  value = paramsList$searchMS1massPpm)
  updateTextInput(         session = session, inputId = "search_ms2_masses1",                value = paramsList$search_ms2_masses1)
  updateTextInput(         session = session, inputId = "search_ms2_masses2",                value = paramsList$search_ms2_masses2)
  updateTextInput(         session = session, inputId = "search_ms2_masses3",                value = paramsList$search_ms2_masses3)
  updateTextInput(         session = session, inputId = "searchMS2massPpm",                  value = paramsList$searchMS2massPpm)
  updateCheckboxInput(     session = session, inputId = "searchIncludeIgnoredPrecursors",    value = as.logical(paramsList$searchIncludeIgnoredPrecursors))
  #input$applySearch
  
  #######################################
  ## update GUI according to parameters
  
  ###################
  ## global MS2 filter
  filter_ms2_masses1  <- paramsList$globalFilter_ms2_masses1
  filter_ms2_masses2  <- paramsList$globalFilter_ms2_masses2
  filter_ms2_masses3  <- paramsList$globalFilter_ms2_masses3
  filter_ms2_ppm      <- paramsList$globalFilter_ms2_ppm
  
  groupSet        <- dataList$grouXXXps
  sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterBySamples <- TRUE
  filter_average  <- NULL
  filter_lfc      <- NULL
  includeIgnoredPrecursors  <- TRUE
  filter_ms1_masses <- NULL
  filter_ms1_ppm <- NULL
  
  resultObj <- doPerformFiltering(
    groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
    filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
    filter_ms1_masses, filter_ms1_ppm, 
    includeIgnoredPrecursors
  )
  filterGlobal <<- resultObj$filter
  state_filters$globalMS2filterValid <<- TRUE
  updateGlobalMS2filterInformation()
  
  ###################
  ## HCA filter
  filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
  filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
  filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
  filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
  
  groupOne        <- paramsList$hcaFilterGroupOne
  groupTwo        <- paramsList$hcaFilterGroupTwo
  filter_average  <- paramsList$hcaFilter_average
  filter_lfc      <- paramsList$hcaFilter_lfc
  includeIgnoredPrecursors  <- paramsList$hcaFilterIncludeIgnoredPrecursors
  filter_ms1_masses <- NULL
  filter_ms1_ppm  <- NULL
  groupSet        <- c(groupOne, groupTwo)
  sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterBySamples <- TRUE
  resultObj <- doPerformFiltering(
    groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
    filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
    filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
  )
  filterHca <<- resultObj$filter
  updateHcaFilterInformation()
  
  numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
  checkHcaFilterValidity(numberOfPrecursorsFiltered)
  
  ###################
  ## draw HCA
  # TODO?
  
  ###################
  ## PCA filter
  filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
  filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
  filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
  filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
  
  groupSet        <- paramsList$pcaGroups
  sampleSet       <- paramsList$pcaSamples
  filterBySamples <- paramsList$filterByPCAgroupSamples
  filter_average  <- paramsList$pcaFilter_average
  filter_lfc      <- paramsList$pcaFilter_lfc
  includeIgnoredPrecursors  <- paramsList$pcaFilterIncludeIgnoredPrecursors
  filter_ms1_masses <- NULL
  filter_ms1_ppm  <- NULL
  
  if(length(groupSet) != 2) {
    filter_lfc <- NULL
  }
  
  resultObj <- doPerformFiltering(
    groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
    filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
    filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
  )
  filterPca <<- resultObj$filter
  updatePcaFilterInformation()
  
  numberOfPrecursorsFiltered <- filterPca$numberOfPrecursorsFiltered
  checkPcaFilterValidity(numberOfPrecursorsFiltered)
  
  ###################
  ## draw PCA
  # TODO
  
  ###################
  ## search
  searchMode <- paramsListsearchMS1orMS2
  if(searchMode == 'MS1 feature m/z'){
    #################################################
    ## get inputs
    filter_ms1_masses <- paramsListsearchMS1mass
    filter_ms1_ppm  <- paramsListsearchMS1massPpm
    
    if(nchar(trimws(filter_ms1_masses)) == 0) {
      return()
    }
    
    filter_ms2_masses1  <- NULL
    filter_ms2_masses2  <- NULL
    filter_ms2_masses3  <- NULL
    filter_ms2_ppm      <- NULL
  }
  
  if(searchMode == 'Fragment m/z'){
    #################################################
    ## get inputs
    filter_ms2_masses1 <- paramsListsearch_ms2_masses1
    filter_ms2_masses2 <- paramsListsearch_ms2_masses2
    filter_ms2_masses3 <- paramsListsearch_ms2_masses3
    filter_ms2_ppm     <- paramsListsearchMS2massPpm
    
    filter_ms1_masses <- NULL
    filter_ms1_ppm <- NULL
  }
  
  filter_lfc      <- NULL
  filter_average  <- NULL
  groupSet        <- dataList$grouXXXps
  sampleSet       <- dataList$includedSamples(dataList$groupSampleDataFrame)
  #sampleSet       <- dataList$dataColumnsNameFunctionFromGroupNames(grouXXXps = groupSet, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
  filterBySamples <- TRUE
  includeIgnoredPrecursors  <- paramsListsearchIncludeIgnoredPrecursors
  
  #################################################
  ## do filtering
  resultObj <- doPerformFiltering(
    groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
    filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
    filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
  )
  processSearchFilterResult(resultObj)
}
