
## create and set table for all six types
createCheckboxInputFields <- function(FUN, id, values, tableCounter) {
  ## running id
  id <- paste(id, tableCounter, sep = "_")
  
  ## create a character vector of shiny inputs
  inputs <- character(length(values))
  for (i in 1:length(values)){
    itemId    <- paste(id, "_", i, sep = "")
    inputs[[i]] <- as.character(FUN(inputId = itemId, label = NULL, value = values[[i]]))
  }
  return(inputs)
}
createCheckboxInputFields2 <- function(FUN, id, values, tableCounter, triggerSampleExclusionClick) {
  ## running id
  id <- paste(id, tableCounter, sep = "_")
  
  ## create a character vector of shiny inputs
  inputs <- character(length(values))
  for (i in 1:length(values)){
    itemId    <- paste(id, "_", i, sep = "")
    inputs[[i]] <- as.character(FUN(inputId = itemId, label = NULL, value = values[[i]]))
    
    ## trigger event on button-click
    if(triggerSampleExclusionClick)
      lapply(c(itemId), function(i){
        observeEvent(input[[i]], {
          sampleExcludeClicked()
        })
      })
  }
  return(inputs)
}
createActionButtonInputFields <- function(FUN, id, itemCount, label = NULL, icon = NULL, tableCounter, callFunction) {
  ## running id
  id <- paste(id, tableCounter, sep = "_")
  
  ## create a character vector of shiny inputs
  inputs <- character(length = itemCount)
  for (i in seq_along(inputs)){
    itemId    <- paste(id, "_", i, sep = "")
    inputs[[i]] <- as.character(FUN(inputId = itemId, label = label, icon = icon))
    
    ## trigger event on button-click
    #lapply(itemId, function(itemId){
    #  observeEvent(input[[itemId]], {
    #    print("huhu")
    #    callFunction(itemId)
    #  })
    #})
    
    #observeEvent(input[[itemId]], {
    #  print("huhu")
    #  print(itemId)
    #  callFunction(itemId)
    #})
  }
  
  itemIds <- paste(id, "_", seq_along(inputs), sep = "")
  lapply(itemIds, function(itemId){
    observeEvent(input[[itemId]], {
      print("hoho")
      callFunction(itemId)
    })
  })
  
  return(inputs)
}
createActionButtonInputFields2 <- function(FUN, id, itemCount, iconUp, iconDown, tableCounter) {
  ## running id
  id1 <- paste(id, "Up",   tableCounter, sep = "_")
  id2 <- paste(id, "Down", tableCounter, sep = "_")
  
  ## create a character vector of shiny inputs
  inputs <- character(length = itemCount)
  for (i in seq_along(inputs)){
    itemId1    <- paste(id1, "_", i, sep = "")
    itemId2    <- paste(id2, "_", i, sep = "")
    
    inputs[[i]] <- as.character(
      fluidRow(
        column(width = 6,
               div(style="float:right",
                   bsTooltip(id = itemId1, title = "Move sample up", placement = "bottom", trigger = "hover"),
                   FUN(inputId = itemId1, label = NULL, icon = iconUp)
               )
        ),##column
        column(width = 6,
               div(style="float:left",
                   bsTooltip(id = itemId2, title = "Move sample down", placement = "bottom", trigger = "hover"),
                   FUN(inputId = itemId2, label = NULL, icon = iconDown)
               )
        )##column
      )##row
    )
    
    ## trigger event on button-click
    lapply(c(itemId1, itemId2), function(i){
      observeEvent(input[[i]], {
        sampleMoveClicked(i)
      })
    })
  }
  return(inputs)
}

## for classifier
createPlotOutput <- function(dataList, id, frequentFragments, characteristicFragments, precursorIndeces, tableCounter, mappingSpectraToClassDf) {
  
  if(FALSE){
    dataList_ <<- dataList
    id_ <<- id
    frequentFragments_ <<- frequentFragments
    characteristicFragments_ <<- characteristicFragments
    precursorIndeces_ <<- precursorIndeces
    tableCounter_ <<- tableCounter
  }
  if(FALSE){
    dataList <<- dataList_
    id <<- id_
    frequentFragments <<- frequentFragments_
    characteristicFragments <<- characteristicFragments_
    precursorIndeces <<- precursorIndeces_
    tableCounter <<- tableCounter_
  }
  
  ## running id
  id <- paste(id, tableCounter, sep = "_")
  
  xInterval <- c(dataList$minimumMass, dataList$maximumMass)
  
  ## class statistics for class plot
  returnObj <- preprocessClassPlot(frequentFragments, characteristicFragments)
  masses_class    <- returnObj$masses_class
  frequency_class <- returnObj$frequency_class
  colors_class    <- returnObj$colors_class
  
  ## create a character vector of shiny inputs
  inputs_i <- character(length(precursorIndeces))
  numberOfMatchingMasses_i <- integer(length(precursorIndeces))
  for (i in seq_along(precursorIndeces)){
    itemId    <- paste(id, "_", i, sep = "")
    assign(x = eval(itemId), value = itemId)
    
    precursorIndex <- precursorIndeces[[i]]
    
    ## match spectrum masses for spectrum plot
    returnObj <- preprocessSpectrumVsClassPlot(dataList, precursorIndex, masses_class, mappingSpectraToClassDf, "Intensity")
    masses_spec <- returnObj$masses_spec
    intensity_spec <- returnObj$intensity_spec
    colors_spec <- returnObj$colors_spec
    numberOfMatchingMasses_i[[i]] <- returnObj$numberOfMatchingMasses
    
    plotAsString <- paste(
      "output$", itemId, " <- renderPlot({",
      #"  print(paste('### SvC ###', ", i, "));",
      #"  plot_ly(x=1:10,y=1:10);",
      #"  plot(1:10);",
      "  calcPlotSpectrumVsClass_small(",
      "masses_spec",     "=", numericVectorToStringForEval(masses_spec    ), ", ",
      "intensity_spec",  "=", numericVectorToStringForEval(intensity_spec ), ", ",
      "colors_spec",     "=", colorVectorToStringForEval(  colors_spec    ), ", ",
      "masses_class",    "=", numericVectorToStringForEval(masses_class   ), ", ",
      "frequency_class", "=", numericVectorToStringForEval(frequency_class), ", ",
      "colors_class",    "=", colorVectorToStringForEval(  colors_class   ), ", ",
      "xInterval",       "=", numericVectorToStringForEval(xInterval      ),
      ");",
      "})",
      sep = ""
    )
    
    ## draw plot
    eval(parse(text=plotAsString), envir = environment())
    
    inputs_i[[i]] <- as.character(
      plotOutput(height = 50, 
                 outputId = itemId
      )
      #plotlyOutput(height = 400, 
      #             outputId = itemId
      #)
    )
  }
  
  returnObj <- list(
    numberOfMatchingMasses_i = numberOfMatchingMasses_i,
    inputs_i = inputs_i
  )
  
  return(returnObj)
}
getNumberOfHits <- function(dataList, frequentFragments, characteristicFragments, precursorIndeces, mappingSpectraToClassDf, properties_class) {
  if(TRUE){
    dataList_ <<- dataList
    frequentFragments_ <<- frequentFragments
    characteristicFragments_ <<- characteristicFragments
    precursorIndeces_ <<- precursorIndeces
    mappingSpectraToClassDf_ <<- mappingSpectraToClassDf
    properties_class_ <<- properties_class
  }
  if(FALSE){
    dataList <- dataList_
    frequentFragments <- frequentFragments_
    characteristicFragments <- characteristicFragments_
    precursorIndeces <- precursorIndeces_
    mappingSpectraToClassDf <- mappingSpectraToClassDf_
    properties_class <- properties_class_
  }
  
  ## class statistics for class plot
  frequentMasses       <- as.numeric(names(frequentFragments)) 
  characteristicMasses <- as.numeric(names(characteristicFragments))
  masses_class <- unique(c(frequentMasses, characteristicMasses))
  
  if(FALSE){
    is_equal_vec <- function(x, y_array, tol = .Machine$double.eps ^ 0.5) {
      abs(x - y_array) < tol
    }
    is_equal_mat <- function(x_array, y_array, tol = .Machine$double.eps ^ 0.5) {
      sapply(X = x_array, FUN = is_equal_vec, y_array)
    }
  }
  
  ## create a character vector of shiny inputs
  tolerance <- .Machine$double.eps ^ 0.5 ## default in function all.equal
  
  
  #mappingSpectraToClassDf <- data.frame(
  #  "SpectraMasses"      = fragmentMasses           [mappedFragmentIndeces_source],
  #  "ClassMasses"        = fragmentMasses_classifier[mappedFragmentIndeces_target],
  #  "SpectraMassIndeces" = mappedFragmentIndeces_source,
  #  "ClassMassIndeces"   = mappedFragmentIndeces_target#,
  ##  "fragmentMasses"            = fragmentMasses,
  ##  "fragmentMasses_classifier" = fragmentMasses_classifier
  #)
  
  ## sum(match(x = masses_spec, table = dataList$fragmentMasses) %in% mappingSpectraToClassDf$SpectraMassIndeces)
  
  if(FALSE){
  numberOfMatchingMasses_i <- sapply(X = precursorIndeces, FUN = function(precursorIndex){
    masses_spec <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)$fragmentMasses
    matchingMassRowIndeces <- which(
      apply(X = outer(X = mappingSpectraToClassDf$SpectraMasses, Y = masses_spec , FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)}) &
      apply(X = outer(X = mappingSpectraToClassDf$ClassMasses  , Y = masses_class, FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)})
    )
    return(length(matchingMassRowIndeces))
  })
  }
  
  spectraMasses <- as.character(mappingSpectraToClassDf$SpectraMasses)
  tmpResult <- as.character(mappingSpectraToClassDf$ClassMasses  ) %in% as.character(masses_class)
  
  numberOfMatchingMasses_i <- sapply(X = precursorIndeces, FUN = function(precursorIndex){
    masses_spec <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)$fragmentMasses
    
    matchingMassRowIndeces <- which(
      spectraMasses %in% as.character(masses_spec) &
      tmpResult
      
      
      #match(x = masses_spec,  table = dataList$fragmentMasses) &
      #match(x = masses_class, table = properties_class[[1]]$fragmentMasses)
      
      #apply(X = outer(X = mappingSpectraToClassDf$SpectraMasses, Y = masses_spec , FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)}) &
      #apply(X = outer(X = mappingSpectraToClassDf$ClassMasses  , Y = masses_class, FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)})
    )
    return(length(matchingMassRowIndeces))
  })
  
  if(FALSE){
  numberOfMatchingMasses_i <- integer(length(precursorIndeces))
  for (i in seq_along(precursorIndeces)){
    #precursorIndex <- precursorIndeces[[i]]
    masses_spec <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndeces[[i]])
    
    #if(length(precursorIndeces) == 1){
    #  resultObj   <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndeces)
    #  masses_spec <- resultObj$fragmentMasses
    #} else {
    #  returnObj   <- getSpectrumStatistics(dataList = dataList, precursorSet = precursorIndeces)
    #  masses_spec <- returnObj$fragmentMasses
    #}
    
    tolerance <- .Machine$double.eps ^ 0.5 ## default in function all.equal
    
    #matchingMassRowIndeces <- 1
    matchingMassRowIndeces <- which(
      apply(X = outer(X = mappingSpectraToClassDf$SpectraMasses, Y = masses_spec , FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)}) &
      apply(X = outer(X = mappingSpectraToClassDf$ClassMasses  , Y = masses_class, FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)})
    )
    numberOfMatchingMasses_i[[i]] <- length(matchingMassRowIndeces)
  }
  }
  
  return(numberOfMatchingMasses_i)
}
getInputValues <- function(id, counter, len) {
  id <- paste(id, "_", counter, sep = "")
  
  ## obtain the values of inputs
  res <- unlist(lapply(1:len, function(i) {
    itemId <- paste(id, "_", i, sep = "")
    value  <- input[[itemId]]
    if (is.null(value))
      return(NA)
    else
      return(value)
  }))
  
  return(res)
}
