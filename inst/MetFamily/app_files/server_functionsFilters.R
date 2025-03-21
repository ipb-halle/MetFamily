
## filter
filterGlobal <- NULL
filterHca <- NULL
filterPca <- NULL
filterSearch <- NULL


state_filters <- reactiveValues(
  ## filter stuff
  globalMS2filterValid = FALSE, 
  hcaFilterValid = FALSE, 
  pcaFilterValid = FALSE,
  searchFilterValid = TRUE,
  filterSearchActive = FALSE
)
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset filters state")
  
  #########################################################################################
  ## update filter
  sampleSet <- dataList$grouXXXpsampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filter <- doPerformFiltering(dataList$grouXXXps, sampleSet, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
  if(length(dataList$grouXXXps) == 1) {
    filter2 <- doPerformFiltering(c(dataList$grouXXXps[[1]], dataList$grouXXXps[[1]]), NULL, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
  } else{
    filter2 <- filter
  }
  
  filterGlobal <<- filter
  filterHca    <<- filter2
  filterPca    <<- filter
  state_filters$filterSearchActive <<- FALSE
  state_filters$searchFilterValid <<- TRUE
  filterSearch    <<- NULL
  
  updateGlobalMS2filterInformation()
  updateHcaFilterInformation()
  updatePcaFilterInformation()
  updateSearchInformation()
  
  state_filters$globalMS2filterValid <<- TRUE
  state_filters$hcaFilterValid <<- TRUE
  state_filters$pcaFilterValid <<- TRUE
  
  checkHcaFilterValidity(filter2$numberOfPrecursorsFiltered)
  checkPcaFilterValidity(filter$numberOfPrecursorsFiltered)
  
  #########################################################################################
  ## update filter input values
  
  ## grouXXXps
  switch(as.character(length(dataList$grouXXXps)), 
         "0"={
           stop("No grouXXXps available")
         },
         "1"={
           selectedOne <- dataList$grouXXXps[[1]]
           selectedTwo <- dataList$grouXXXps[[1]]
         },
         {
           selectedOne <- dataList$grouXXXps[[1]]
           selectedTwo <- dataList$grouXXXps[[2]]
         }
  )
  
  sampleNames <- dataList$groupSampleDataFrame[, "Sample"]
  
  ## input fields: global MS/MS filter
  updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
  updateTextInput(session = session, inputId = "globalFilter_ms2_masses2", value = "")
  updateTextInput(session = session, inputId = "globalFilter_ms2_masses3", value = "")
  updateTextInput(session = session, inputId = "globalFilter_ms2_ppm", value = "20")
  
  ## input fields: HCA filter
  updateRadioButtons(session = session, inputId = "hcaFilterGroupOne", choices = dataList$grouXXXps, selected = selectedOne)
  updateRadioButtons(session = session, inputId = "hcaFilterGroupTwo", choices = dataList$grouXXXps, selected = selectedTwo)
  updateTextInput(session = session, inputId = "hcaFilter_average", value = "0")
  updateTextInput(session = session, inputId = "hcaFilter_lfc", value = "0")
  updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
  
  ## input fields: PCA filter
  updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$grouXXXps, selected = dataList$grouXXXps)
  updateCheckboxGroupInput(session = session, inputId = "pcaSamples",  choices = sampleNames,     selected = sampleNames)
  updateTextInput(session = session, inputId = "pcaFilter_average", value = "0")
  updateTextInput(session = session, inputId = "pcaFilter_lfc", value = "0")
  updateCheckboxInput(session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", value = FALSE)
})


## filter info
updateGlobalMS2filterInformation <- function(){
  if(is.null(filterGlobal)){
    output$globalMS2filteredPrecursors <- renderText({
      print(paste("update output$globalMS2filteredPrecursors invalid filters", sep = ""))
      paste("There are invalid filter values", sep = "")
    })
  } else {
    output$globalMS2filteredPrecursors <- renderText({
      print(paste("update output$globalMS2filteredPrecursors", sep = ""))
      paste("Number of filtered MS1 features: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, sep = "")
    })
  }
}
updateHcaFilterInformation <- function(){
  if(is.null(filterHca)){
    ## errors
    output$hcaFilteredPrecursors <- renderText({
      print(paste("update output$hcaFilteredPrecursors invalid filters", sep = ""))
      paste("There are invalid filter values.", sep = "")
    })
  } else {
    ## no errors
    globalMs2Filter <- ifelse(filterGlobal$numberOfPrecursorsFiltered != dataList$numberOfPrecursors, paste("\n(Global MS/MS filter: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, " precursors)", sep = ""), "")
    
    if(filterHca$numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & filterHca$numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
      ## filter valid
      output$hcaFilteredPrecursors <- renderText({
        print(paste("update output$hcaFilteredPrecursors ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
        paste("Number of filtered MS1 features: ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, globalMs2Filter, sep = "")
      })
    } else {
      ## filter invalid
      
      ## update info
      if(filterHca$numberOfPrecursorsFiltered == 0){
        output$hcaFilteredPrecursors <- renderText({
          print(paste("update output$hcaFilteredPrecursors # = 0", sep = ""))
          paste("There are no MS\u00B9 features which fulfill the given criteria.", globalMs2Filter, sep = "")
        })
      }
      if(filterHca$numberOfPrecursorsFiltered > 0 & filterHca$numberOfPrecursorsFiltered < minimumNumberOfPrecursorsForHca){
        output$hcaFilteredPrecursors <- renderText({
          print(paste("update output$hcaFilteredPrecursors 0 < # < ", minimumNumberOfPrecursorsForHca, sep = ""))
          paste("There are only ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, " MS\u00B9 features which fulfill the given criteria. There must be at least more than five MS\u00B9 features to proceed.", globalMs2Filter, sep = "")
        })
      }
      if(filterHca$numberOfPrecursorsFiltered > maximumNumberOfPrecursorsForHca){
        output$hcaFilteredPrecursors <- renderText({
          print(paste("update output$hcaFilteredPrecursors # > ", maximumNumberOfPrecursorsForHca, sep = ""))
          paste("There are ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, " MS\u00B9 features which fulfill the given criteria. There must be at most ", maximumNumberOfPrecursorsForHca, " MS\u00B9 features to proceed.", globalMs2Filter, sep = "")
        })
      }
    }
  }
}
updatePcaFilterInformation <- function(){
  if(is.null(filterPca)){
    output$pcaFilteredPrecursors <- renderText({
      print(paste("update output$pcaFilteredPrecursors invalid filters", sep = ""))
      paste("There are invalid filter values.", sep = "")
    })
  } else {
    globalMs2Filter <- ifelse(filterGlobal$numberOfPrecursorsFiltered != dataList$numberOfPrecursors, paste("\n(Global MS/MS filter: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, " precursors)", sep = ""), "")
    if(filterPca$numberOfPrecursorsFiltered > 0){
      output$pcaFilteredPrecursors <- renderText({
        print(paste("update output$pcaFilteredPrecursors", sep = ""))
        paste("Number of filtered MS1 features: ", filterPca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, globalMs2Filter, sep = "")
      })
    } else {
      output$pcaFilteredPrecursors <- renderText({
        print(paste("update output$pcaFilteredPrecursors", sep = ""))
        paste("There are no MS\u00B9 features which fulfill the given criteria.", globalMs2Filter, sep = "")
      })
    }
  }
}

## perform filtering
doPerformFiltering <- function(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, preFilter = NULL){
  suppressWarnings(
    doPerformFiltering_impl(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, preFilter)
  )
}

#' Main internal filtering function
#' 
#' Called from `doPerformFiltering()`. Takes dataList from global environment.
#'
#' @returns filterObject
#' @noRd
doPerformFiltering_impl <- function(
    groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, 
    filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, 
    filter_ms1_ppm, includeIgnoredPrecursors, preFilter = NULL){
  print(paste("Observe applyFilters1", "gs", paste(groupSet, collapse = "-"), "a", filter_average, "lfc", filter_lfc, "ms2_1", filter_ms2_masses1, "ms2_2", filter_ms2_masses2, "ms2_3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "ig", includeIgnoredPrecursors))
  
  groupSetOriginal                 <- groupSet
  sampleSetOriginal                <- sampleSet
  filterBySamplesOriginal          <- filterBySamples
  filter_averageOriginal           <- filter_average
  filter_lfcOriginal               <- filter_lfc
  filter_ms2_masses1Original       <- filter_ms2_masses1
  filter_ms2_masses2Original       <- filter_ms2_masses2
  filter_ms2_masses3Original       <- filter_ms2_masses3
  filter_ms2_ppmOriginal           <- filter_ms2_ppm
  filter_ms1_massesOriginal        <- filter_ms1_masses
  filter_ms1_ppmOriginal           <- filter_ms1_ppm
  includeIgnoredPrecursorsOriginal <- includeIgnoredPrecursors
  
  #################################################
  ## parse inputs
  print(paste("Observe applyFilters2", "gs", paste(groupSet, collapse = "-"), "ss", paste(sampleSet, collapse = "-"), "fbss", filterBySamples, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "i", includeIgnoredPrecursors))
  print(paste("Observe applyFilters3", "gs", is.null(groupSet), "ss", is.null(sampleSet), "gs", is.null(filterBySamples), "a", is.null(filter_average), "lfc", is.null(filter_lfc), "ms2 1", is.null(filter_ms2_masses1), "ms2 2", is.null(filter_ms2_masses2), "ms2 3", is.null(filter_ms2_masses3), "ppm", is.null(filter_ms2_ppm), "i", includeIgnoredPrecursors))
  
  #################################################
  ## sanity checks
  if(all(!is.null(filter_lfc), !is.na(filter_lfc), filter_lfc != 0) & length(groupSet) != 2) {
    stop("lfc filter for not exactly two groups")
  }
  
  #################################################
  ## check for errors in inputs amd process ms2
  error <- FALSE
  if(any(is.null(groupSet), is.na(groupSet), length(groupSet) == 0, any(nchar(groupSet) == 0))) {
    error <- TRUE
  }
  
  if(any(is.null(filter_average), is.na(filter_average), length(filter_average) == 0, nchar(filter_average) == 0)) {
    filter_average <- NULL
  } else {
    filter_average <- as.numeric(filter_average)
    error <- error | is.na(filter_average)
  }
  
  if(any(is.null(filter_lfc), is.na(filter_lfc), length(filter_lfc) == 0, nchar(filter_lfc) == 0)) {
    filter_lfc <- NULL
  } else {
    filter_lfc <- as.numeric(filter_lfc)
    error <- error | is.na(filter_lfc)
  }
  
  if(any(is.null(filter_ms2_masses1), is.na(filter_ms2_masses1), length(filter_ms2_masses1) == 0, nchar(filter_ms2_masses1) == 0)) {
    filter_ms2_masses1 <- NULL
  } else {
    ms2Masses <- strsplit(x = filter_ms2_masses1, split = "[,; ]+")[[1]]
    filter_ms2_masses1 <- vector(mode = "numeric", length = length(ms2Masses))
    for(idx in 1:length(ms2Masses)) {
      filter_ms2_masses1[[idx]] <- as.numeric(ms2Masses[[idx]])
    }
    error <- error | any(is.na(filter_ms2_masses1))
  }
  if(any(is.null(filter_ms2_masses2), is.na(filter_ms2_masses2), length(filter_ms2_masses2) == 0, nchar(filter_ms2_masses2) == 0)){
    filter_ms2_masses2 <- NULL
  } else {
    ms2Masses <- strsplit(x = filter_ms2_masses2, split = "[,; ]+")[[1]]
    filter_ms2_masses2 <- vector(mode = "numeric", length = length(ms2Masses))
    for(idx in 1:length(ms2Masses))
      filter_ms2_masses2[[idx]] <- as.numeric(ms2Masses[[idx]])
    error <- error | any(is.na(filter_ms2_masses2))
  }
  if(any(is.null(filter_ms2_masses3), is.na(filter_ms2_masses3), length(filter_ms2_masses3) == 0, nchar(filter_ms2_masses3) == 0)) {
    filter_ms2_masses3 <- NULL
  } else {
    ms2Masses <- strsplit(x = filter_ms2_masses3, split = "[,; ]+")[[1]]
    filter_ms2_masses3 <- vector(mode = "numeric", length = length(ms2Masses))
    for(idx in 1:length(ms2Masses))
      filter_ms2_masses3[[idx]] <- as.numeric(ms2Masses[[idx]])
    error <- error | any(is.na(filter_ms2_masses3))
  }
  
  if(any(is.null(filter_ms2_ppm), is.na(filter_ms2_ppm), length(filter_ms2_ppm) == 0, nchar(filter_ms2_ppm) == 0)) {
    filter_ms2_ppm <- NULL
  } else {
    filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
    error <- error | is.na(filter_ms2_ppm)
  }
  
  if(any(is.null(filter_ms1_masses), is.na(filter_ms1_masses), length(filter_ms1_masses) == 0, nchar(filter_ms1_masses) == 0)) {
    filter_ms1_masses <- NULL
  } else {
    ms1Masses <- strsplit(x = filter_ms1_masses, split = "[,; ]+")[[1]]
    filter_ms1_masses <- vector(mode = "numeric", length = length(ms1Masses))
    for(idx in 1:length(ms1Masses))
      filter_ms1_masses[[idx]] <- as.numeric(ms1Masses[[idx]])
    error <- error | any(is.na(filter_ms1_masses))
  }
  
  if(any(is.null(filter_ms1_ppm), is.na(filter_ms1_ppm), length(filter_ms1_ppm) == 0, nchar(filter_ms1_ppm) == 0)) {
    filter_ms1_ppm <- NULL
  } else {
    filter_ms1_ppm <- as.numeric(filter_ms1_ppm)
    error <- error | is.na(filter_ms1_ppm)
  }
  
  ## sanity check
  error <- error | (!is.null(filter_ms1_masses) & any(is.null(filter_ms1_ppm), is.na(filter_ms1_ppm)))
  error <- error | (any(!is.null(filter_ms2_masses1), !is.null(filter_ms2_masses2), !is.null(filter_ms2_masses3)) & any(is.null(filter_ms2_ppm), is.na(filter_ms2_ppm)))
  
  print(paste("Observe applyFilters4", "e", error, "gs", paste(groupSet, collapse = "-"), "ss", paste(sampleSet, collapse = "-"), "fbss", filterBySamples, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "i", includeIgnoredPrecursors))
  
  ## collect ms2 masses
  filterList_ms2_masses <- list()
  if(!is.null(filter_ms2_masses1))
    filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses1
  if(!is.null(filter_ms2_masses2))
    filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses2
  if(!is.null(filter_ms2_masses3))
    filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses3
  
  #################################################
  ## do filtering and wrap results
  resultObj <- list()
  resultObj$error  <- error
  
  if(error){
    resultObj$filter <- NULL
  } else {
    filterHere <- filterData(
      dataList = dataList, 
      #groupOne = groupOne, groupTwo = groupTwo, 
      grouXXXps = groupSet, sampleSet, filterBySamples, filter_average = filter_average, filter_lfc = filter_lfc, 
      filterList_ms2_masses = filterList_ms2_masses, filter_ms2_ppm = filter_ms2_ppm, 
      filter_ms1_masses = filter_ms1_masses, filter_ms1_ppm = filter_ms1_ppm,
      includeIgnoredPrecursors = includeIgnoredPrecursors,
      progress = FALSE
    )
    
    ## set original values
    if(is.null(groupSetOriginal)){
      filterHere$groupSetOriginal <- list()
    } else {
      filterHere$groupSetOriginal  <- groupSetOriginal
      filterHere$sampleSetOriginal <- sampleSetOriginal
      filterHere$filterBySamplesOriginal <- filterBySamplesOriginal
    }
    #filterHere$groupSetOriginal                 <- ifelse(test = is.null(groupSetOriginal),                 yes = NA, no = groupSetOriginal)
    filterHere$filter_averageOriginal           <- ifelse(test = is.null(filter_averageOriginal),           yes = 0,  no = filter_averageOriginal)
    filterHere$filter_lfcOriginal               <- ifelse(test = is.null(filter_lfcOriginal),               yes = 0,  no = filter_lfcOriginal)
    if(is.null(filter_ms2_masses1Original)){
      filterHere$filter_ms2_masses1Original <- list()
    } else {
      filterHere$filter_ms2_masses1Original <- filter_ms2_masses1Original
    }
    #filterHere$filter_ms2_masses1Original       <- ifelse(test = is.null(filter_ms2_masses1Original),       yes = "", no = filter_ms2_masses1Original)
    if(is.null(filter_ms2_masses2Original)){
      filterHere$filter_ms2_masses2Original <- list()
    } else {
      filterHere$filter_ms2_masses2Original <- filter_ms2_masses2Original
    }
    #filterHere$filter_ms2_masses2Original       <- ifelse(test = is.null(filter_ms2_masses2Original),       yes = "", no = filter_ms2_masses2Original)
    if(is.null(filter_ms2_masses3Original)){
      filterHere$filter_ms2_masses3Original <- list()
    } else {
      filterHere$filter_ms2_masses3Original <- filter_ms2_masses3Original
    }
    #filterHere$filter_ms2_masses3Original       <- ifelse(test = is.null(filter_ms2_masses3Original),       yes = "", no = filter_ms2_masses3Original)
    filterHere$filter_ms2_ppmOriginal           <- ifelse(test = is.null(filter_ms2_ppmOriginal),           yes = "", no = filter_ms2_ppmOriginal)
    if(is.null(filter_ms1_massesOriginal)){
      filterHere$filter_ms1_massesOriginal <- list()
    } else {
      filterHere$filter_ms1_massesOriginal <- filter_ms1_massesOriginal
    }
    #filterHere$filter_ms1_massesOriginal        <- ifelse(test = is.null(filter_ms1_massesOriginal),        yes = "", no = filter_ms1_massesOriginal)
    filterHere$filter_ms1_ppmOriginal           <- ifelse(test = is.null(filter_ms1_ppmOriginal),           yes = "", no = filter_ms1_ppmOriginal)
    filterHere$includeIgnoredPrecursorsOriginal <- ifelse(test = is.null(includeIgnoredPrecursorsOriginal), yes = NA, no = includeIgnoredPrecursorsOriginal)
    
    resultObj$error  <- error
    resultObj$filter <- filterHere
    print(paste("Observe applyFilters5", "n", resultObj$filter$numberOfPrecursorsFiltered))
  }
  
  return(resultObj)
}

checkHcaFilterValidity <- function(numberOfPrecursorsFiltered){
  if(numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
    ## filter valid
    print(paste("Observe applyFilters ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
    
    shinyjs::enable("drawHCAplots")
    #enableActionButton(session, "drawHCAplots")
    state_filters$hcaFilterValid <<- TRUE
  } else {
    ## filter invalid
    
    shinyjs::disable("drawHCAplots")
    #disableActionButton(session, "drawHCAplots")
    state_filters$hcaFilterValid <<- FALSE
  }
}
checkPcaFilterValidity <- function(numberOfPrecursorsFiltered){
  if(numberOfPrecursorsFiltered > 0){
    ## filter valid
    print(paste("Observe applyFilters # > 0", sep = ""))
    
    shinyjs::enable("drawPCAplots")
    state_filters$pcaFilterValid <<- TRUE
  } else {
    ## filter invalid
    print(paste("Observe applyFilters # = 0", sep = ""))
    
    shinyjs::disable("drawPCAplots")
    state_filters$pcaFilterValid <<- FALSE
  }
}

applyGlobalMS2filters <- function(filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm){
  groupSet        <- dataList$grouXXXps
  filter_average  <- NULL
  filter_lfc      <- NULL
  includeIgnoredPrecursors  <- TRUE
  filter_ms1_masses <- NULL
  filter_ms1_ppm <- NULL
  
  #################################################
  ## do filtering
  sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterBySamples <- TRUE
  resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
  
  if(resultObj$error){
    filterGlobal <<- NULL
    state_filters$globalMS2filterValid <<- FALSE
  } else {
    filterGlobal <<- resultObj$filter
    state_filters$globalMS2filterValid <<- TRUE
  }
  updateGlobalMS2filterInformation()
}
obsApplyHcaFilters <- observeEvent(input$applyHcaFilters, {
  session$sendCustomMessage("disableButton", "applyHcaFilters")
  applyHcaFilters <- as.numeric(input$applyHcaFilters)
  
  print(paste("Observe applyHcaFilters", applyHcaFilters))
  
  #################################################
  ## check if button was hit
  #if(applyHcaFilters == applyHcaFiltersButtonValue)
  #  return()
  #applyHcaFiltersButtonValue <<- applyHcaFilters
  
  #################################################
  ## get inputs
  groupOne        <- input$hcaFilterGroupOne
  groupTwo        <- input$hcaFilterGroupTwo
  filter_average  <- input$hcaFilter_average
  filter_lfc      <- input$hcaFilter_lfc
  includeIgnoredPrecursors  <- input$hcaFilterIncludeIgnoredPrecursors
  
  applyHcaFilters(groupOne, groupTwo, filter_average, filter_lfc, includeIgnoredPrecursors)
  session$sendCustomMessage("enableButton", "applyHcaFilters")
})
obsClearHcaFilters <- observeEvent(input$clearHcaFilters, {
  session$sendCustomMessage("disableButton", "clearHcaFilters")
  clearHcaFilters <- as.numeric(input$clearHcaFilters)
  
  print(paste("Observe clearHcaFilters", clearHcaFilters))
  
  #################################################
  ## check if button was hit
  #if(clearHcaFilters == clearHcaFiltersButtonValue)
  #  return()
  #clearHcaFiltersButtonValue <<- clearHcaFilters
  
  #################################################
  ## get inputs
  groupOne        <- input$hcaFilterGroupOne
  groupTwo        <- input$hcaFilterGroupTwo
  filter_average  <- ""
  filter_lfc      <- ""
  includeIgnoredPrecursors  <- TRUE
  
  applyHcaFilters(groupOne, groupTwo, filter_average, filter_lfc, includeIgnoredPrecursors)
  session$sendCustomMessage("enableButton", "clearHcaFilters")
})
applyHcaFilters <- function(groupOne, groupTwo, filter_average, filter_lfc, includeIgnoredPrecursors){
  filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
  filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
  filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
  filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
  
  filter_ms1_masses <- NULL
  filter_ms1_ppm  <- NULL
  groupSet        <- c(groupOne, groupTwo)
  sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterBySamples <- TRUE
  #################################################
  ## do filtering and update
  resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
  
  if(resultObj$error){
    shinyjs::disable("drawHCAplots")
    #disableActionButton(session, "drawHCAplots")
    filterHca <<- NULL
    updateHcaFilterInformation()
    state_filters$hcaFilterValid <<- FALSE
    return()
  }
  
  #################################################
  ## check filter validity
  filterHca <<- resultObj$filter
  updateHcaFilterInformation()
  
  numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
  checkHcaFilterValidity(numberOfPrecursorsFiltered)
}
obsApplyPcaFilters <- observeEvent(input$applyPcaFilters, {
  session$sendCustomMessage("disableButton", "applyPcaFilters")
  applyPcaFilters <- as.numeric(input$applyPcaFilters)
  
  print(paste("Observe applyPcaFilters", applyPcaFilters))
  
  #################################################
  ## check if button was hit
  #if(applyPcaFilters == applyPcaFiltersButtonValue)
  #  return()
  #applyPcaFiltersButtonValue <<- applyPcaFilters
  
  #################################################
  ## get inputs
  groupSet        <- input$pcaGroups
  sampleSet       <- input$pcaSamples
  filterBySamples <- input$filterByPCAgroupSamples
  filter_average  <- input$pcaFilter_average
  filter_lfc      <- input$pcaFilter_lfc
  includeIgnoredPrecursors  <- input$pcaFilterIncludeIgnoredPrecursors
  
  if(filterBySamples){
    ## update grouXXXps and samples mutually
    
    ## grouXXXps which are covered by at least one sample
    groupsFromSamples <- unlist(lapply(X = dataList$grouXXXps, FUN = function(x){
      samplesOfGroups <- dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
      if(any(samplesOfGroups %in% sampleSet))
        return(x)
      else
        return(NULL)
    }))
    
    ## samples which ae covered by a group
    samplesFromGroups <- dataList$dataColumnsNameFunctionFromGroupNames(grouXXXps = groupSet, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
    groupSet  <- intersect(groupSet, groupsFromSamples)
    sampleSet <- intersect(sampleSet, samplesFromGroups)
  } else {
    sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  }
  
  if(length(groupSet) != 2)
    filter_lfc <- NULL
  
  applyPcaFilters(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, includeIgnoredPrecursors)
  session$sendCustomMessage("enableButton", "applyPcaFilters")
})
obsClearPcaFilters <- observeEvent(input$clearPcaFilters, {
  session$sendCustomMessage("disableButton", "clearPcaFilters")
  clearPcaFilters <- as.numeric(input$clearPcaFilters)
  
  print(paste("Observe clearPcaFilters", clearPcaFilters))
  
  #################################################
  ## check if button was hit
  #if(clearPcaFilters == clearPcaFiltersButtonValue)
  #  return()
  #clearPcaFiltersButtonValue <<- clearPcaFilters
  
  #################################################
  ## get inputs
  groupSet        <- dataList$grouXXXps
  sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
  filterByPCAgroupSamples <- TRUE
  filter_average  <- ""
  filter_lfc      <- ""
  includeIgnoredPrecursors  <- TRUE
  
  if(length(groupSet) != 2)
    filter_lfc <- NULL
  
  applyPcaFilters(groupSet, sampleSet, filterByPCAgroupSamples, filter_average, filter_lfc, includeIgnoredPrecursors)
  session$sendCustomMessage("enableButton", "clearPcaFilters")
})

applyPcaFilters_default <- function(){
  if(!is.null(filterPca)){
    applyPcaFilters(
      groupSet = filterPca$groupSetOriginal, 
      sampleSet = filterPca$sampleSetOriginal, 
      filterBySamples = filterPca$filterBySamplesOriginal, 
      filter_average = filterPca$filter_averageOriginal, 
      filter_lfc = filterPca$filter_lfcOriginal, 
      includeIgnoredPrecursors = filterPca$includeIgnoredPrecursorsOriginal
    )
  } else {
    stop("Tried to apply pca filtering without filterPca obj")
  }
}
applyPcaFilters <- function(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, includeIgnoredPrecursors){
  filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
  filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
  filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
  filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
  
  filter_ms1_masses <- NULL
  filter_ms1_ppm  <- NULL
  
  #################################################
  ## do filtering and update
  resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
  
  if(resultObj$error){
    shinyjs::disable("drawPCAplots")
    #disableActionButton(session, "drawPCAplots")
    filterPca <<- NULL
    updatePcaFilterInformation()
    state_filters$pcaFilterValid <<- FALSE
    return()
  }
  
  filterPca <<- resultObj$filter
  updatePcaFilterInformation()
  
  numberOfPrecursorsFiltered <- filterPca$numberOfPrecursorsFiltered
  checkPcaFilterValidity(numberOfPrecursorsFiltered)
}

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending filters observers")
  obsApplyHcaFilters$suspend()
  obsClearHcaFilters$suspend()
  obsApplyPcaFilters$suspend()
  obsClearPcaFilters$suspend()
})

output$globalMS2filterValid <- reactive({
  print(paste("reactive update globalMS2filterValid", state_filters$globalMS2filterValid))
  return(state_filters$globalMS2filterValid)
})
output$hcaFilterValid <- reactive({
  print(paste("reactive update hcaFilterValid", state_filters$hcaFilterValid))
  return(state_filters$hcaFilterValid)
})
output$pcaFilterValid <- reactive({
  print(paste("reactive update pcaFilterValid", state_filters$pcaFilterValid))
  return(state_filters$pcaFilterValid)
})
output$searchFilterValid <- reactive({
  print(paste("reactive update searchFilterValid", state_filters$searchFilterValid))
  return(state_filters$searchFilterValid)
})
output$filterSearchActive <- reactive({
  print(paste("reactive update filterSearchActive", state_filters$filterSearchActive))
  return(state_filters$filterSearchActive)
})

outputOptions(output, 'globalMS2filterValid',    suspendWhenHidden=FALSE)
outputOptions(output, 'hcaFilterValid',          suspendWhenHidden=FALSE)
outputOptions(output, 'pcaFilterValid',          suspendWhenHidden=FALSE)
outputOptions(output, 'searchFilterValid',       suspendWhenHidden=FALSE)
outputOptions(output, 'filterSearchActive',      suspendWhenHidden=FALSE)
