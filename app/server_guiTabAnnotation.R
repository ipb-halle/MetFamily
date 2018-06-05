
## metabolite families
allAnnotationNames <- NULL

fragmentMassesForFamily <- NULL
fragmentPlot2Range    <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)

updateAnnotationOverview <- function(){
  allAnnotationNames  <- unlist(dataList$annoPresentAnnotationsList)
  allAnnotationColors <- unlist(dataList$annoPresentColorsList)
  annoCounts          <- sapply(X = allAnnotationNames, FUN = function(x){sum(unlist(lapply(X = dataList$annoArrayOfLists, FUN = function(y){any(y==x)})))})
  
  these <- annoCounts > 0
  allAnnotationNames  <- allAnnotationNames[these]
  allAnnotationColors <- allAnnotationColors[these]
  annoCounts          <- annoCounts[these]
  
  allAnnotationNames <<- allAnnotationNames
  
  metaboliteFamiliesDf <- data.frame(
    "Family" = allAnnotationNames,
    "Color"  = allAnnotationColors,
    "Count"  = annoCounts
  )
  dataTable <- datatable(
    data = metaboliteFamiliesDf, 
    escape = FALSE, rownames = FALSE,
    selection = list(mode = "single"),#, selected = selected), 
    options = list(
      #scrollY = "600px",
      scrollY = "20vh",
      preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
      #rowCallback = JS(
      #  "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
      #  "var full_text = classes[[nRow]]",
      #  "$('td:eq(1)', nRow).attr('title', full_text);",
      #  "}"),
      iDisplayLength=nrow(metaboliteFamiliesDf),       # initial number of records
      #aLengthMenu = c(5,10),    # records/page options
      #bLengthChange =0,        # show/hide records per page dropdown
      #bFilter = 0,              # global search box on/off
      #bInfo = 0,                # information on/off (how many records filtered, etc)
      #bAutoWidth = 0,           # automatic column width calculation, disable if passing column width via aoColumnDefs
      #aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
      ordering = F,              # row ordering
      sDom  = 't'
      #sDom  = '<"top">rt<"bottom">ip'
    ))
  
  output$familySelectionTable <- DT::renderDataTable(
    expr = formatStyle(table = dataTable, columns = "Color", target = "cell", backgroundColor = styleEqual(levels = allAnnotationColors, values = allAnnotationColors)),
    server = FALSE
  )
  
  output$familyCount2 <- renderText({
    print(paste("init output$familyCount2"))
    paste("Available metabolite families:", nrow(metaboliteFamiliesDf))
  })
}
obsFamilySelectionTable_rows_selected <- observeEvent(eventExpr = input$familySelectionTable_rows_selected, handlerExpr = {
  print(paste("Observe familySelectionTable_rows_selected", input$familySelectionTable_rows_selected))
  
  state$metaboliteFamilySelected <<- !is.null(input$familySelectionTable_rows_selected)
  
  if(is.null(input$familySelectionTable_rows_selected))
    return()
  
  annotation <- allAnnotationNames[[input$familySelectionTable_rows_selected]]
  precursorSet <- which(unlist(lapply(X = dataList$annoArrayOfLists, FUN = function(y){any(y==annotation)})))
  numberOfPrecursors <- length(precursorSet)
  
  ###################################################
  ## family info
  output$featureCountForFamily <- renderText({
    print(paste("init output$featureCountForFamily"))
    paste(numberOfPrecursors, " MS\u00B9 features for metabolite family '", annotation, "'", sep = "")
  })
  
  ###################################################
  ## consensus spectrum (vs class spectrum)
  plotMetaboliteFamilyVersusClass(precursorSet)
  
  ###################################################
  ## ms1 feature table
  ms1FeatureTable <- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
  ms1FeatureDf <- createMS1FeatureTable2(ms1FeatureTable)
  
  output$ms1FeatureTableForAnnotation <- DT::renderDataTable(
    expr = ms1FeatureDf,
    server = FALSE, escape = FALSE, selection = "none",
    options = list(
      preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }')
    )
  )
}, ignoreNULL = FALSE)
plotMetaboliteFamilyVersusClass <- function(precursorSet){
  selectedClassifierClass <- isolate(input$metaboliteFamilyComparisonClass)
  
  addClassifierConsensusSpectrum <- state$classifierLoaded & selectedClassifierClass != selectionNone
  returnObj <- metaboliteFamilyVersusClass(dataList = dataList, precursorSet = precursorSet, classToSpectra_class = classToSpectra_class, classifierClass = selectedClassifierClass, addClassifierConsensusSpectrum = addClassifierConsensusSpectrum, mappingSpectraToClassDf = mappingSpectraToClassDf)
  masses_spec     <- returnObj$masses_spec
  frequency_spec  <- returnObj$intensity_spec
  colors_spec     <- returnObj$colors_spec
  masses_class    <- returnObj$masses_class
  frequency_class <- returnObj$frequency_class
  colors_class    <- returnObj$colors_class
  
  fragmentMassesForFamily <<- unique(c(masses_spec, masses_class))
  
  min <- min(fragmentMassesForFamily)
  max <- max(fragmentMassesForFamily)
  fragmentPlot2Range$xMin <<- min
  fragmentPlot2Range$xMax <<- max
  fragmentPlot2Range$xInterval <<- c(min, max)
  fragmentPlot2Range$xIntervalSize <<- max - min
  
  output$fragmentPlot2 <- renderPlot({
    print(paste("### MS2 ### family init"))
    #plotFragments2(masses = masses_spec, numberOfFragments = frequency_spec, numberOfPrecursors = numberOfPrecursors, xInterval = fragmentPlot2Range$xInterval)
    #plotFragments2(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
    #}, bg = "transparent")
    calcPlotSpectrumVsClass_big(
      masses_spec     = masses_spec, 
      intensity_spec  = frequency_spec, 
      colors_spec     = colors_spec, 
      masses_class    = masses_class, 
      frequency_class = frequency_class, 
      colors_class    = colors_class, 
      singleSpec      = FALSE,
      xInterval       = fragmentPlot2Range$xInterval
    )
  })
}

createMS1FeatureTable2 <- function(list){
  if(is.null(list$precursorSet))
    return(NULL)
  
  ## get (possibly truncated) data
  if(length(list$precursorSet) > maximumNumberOfTableEntries){
    precursorSet          <- list$precursorSet         [1:maximumNumberOfTableEntries]
    ms1abundanceDataFrame <- list$ms1abundanceDataFrame[1:maximumNumberOfTableEntries, , drop=FALSE]
    annotationDataFrame   <- list$annotationDataFrame  [1:maximumNumberOfTableEntries, , drop=FALSE]
    ms2fragmentDataFrame  <- list$ms2fragmentDataFrame [1:maximumNumberOfTableEntries, , drop=FALSE]
  } else {
    precursorSet          <- list$precursorSet
    ms1abundanceDataFrame <- list$ms1abundanceDataFrame
    annotationDataFrame   <- list$annotationDataFrame
    ms2fragmentDataFrame  <- list$ms2fragmentDataFrame
  }
  
  ## assemble
  dataFrame <<- cbind(
    #ms1abundanceDataFrame,
    annotationDataFrame
    #ms2fragmentDataFrame
  )
  
  return(dataFrame)
}

obsFragmentPlot2dblClick <- observeEvent(input$fragmentPlot2_dblclick, {
  brush <- input$fragmentPlot2_brush
  print(paste("observe fragmentPlot2 dblclick", paste(brush, collapse = "; ")))
  
  if (!is.null(brush)) {
    ## set range
    min <- brush$xmin
    max <- brush$xmax
  } else {
    ## reset range
    min <- min(fragmentMassesForFamily)
    max <- max(fragmentMassesForFamily)
  }
  
  fragmentPlot2Range$xMin <<- min
  fragmentPlot2Range$xMax <<- max
  fragmentPlot2Range$xInterval <<- c(min, max)
  fragmentPlot2Range$xIntervalSize <<- max - min
})
obsClassSelection <- observeEvent(input$metaboliteFamilyComparisonClass, {
  print(paste("Observe metaboliteFamilyComparisonClass", input$metaboliteFamilyComparisonClass))
  
  annotationIdx <- isolate(input$familySelectionTable_rows_selected)
  if(is.null(annotationIdx))
    return()
  
  annotation <- allAnnotationNames[[annotationIdx]]
  precursorSet <- which(unlist(lapply(X = dataList$annoArrayOfLists, FUN = function(y){any(y==annotation)})))
  plotMetaboliteFamilyVersusClass(precursorSet)
})

obsSelectMetaboliteFamily <- observeEvent(input$selectMetaboliteFamily, {
  selectMetaboliteFamily <- as.numeric(input$selectMetaboliteFamily)
  print(paste("Observe selectMetaboliteFamily", selectMetaboliteFamily))
  
  ## TODO
  selectedAnnotationTableRowIdx <- isolate(input$familySelectionTable_rows_selected)
  annotationName  <- allAnnotationNames[[selectedAnnotationTableRowIdx]]
  
  annoThere <- lapply(X = dataList$annoArrayOfLists, FUN = function(x){ match(x = annotationName, table = x) })
  precursorSet <- which(!is.na(annoThere))
  
  #selectionBySearch(filterSearch$filter)
  
  #################################################
  ## update plots
  if(state$plotHcaShown){
    precursorSet <- intersect(precursorSet, filterHca$filter)
    selectionAnalysisTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
    
    selectionByHca2(precursorSet)
    
    ## TODO remove plot call
    drawDendrogramPlot(consoleInfo = "annotation selection", withHeatmap = TRUE)
  }
  if(state$plotPcaShown){
    precursorSet <- intersect(precursorSet, filterPca$filter)
    selectionAnalysisPcaLoadingSet <<- which(filterPca$filter %in% precursorSet)
    selectionByPca2(precursorSet)
    
    ## TODO remove plot call
    drawPcaLoadingsPlot(consoleInfo = "PCA loadings click output$plotPcaLoadings")
  }
})
obsRenameMetaboliteFamily <- observeEvent(input$renameMetaboliteFamily, {
  renameMetaboliteFamily <- as.numeric(input$renameMetaboliteFamily)
  print(paste("Observe renameMetaboliteFamily", renameMetaboliteFamily))
  
  renameMetaboliteFamily()
})
obsRemoveMetaboliteFamily <- observeEvent(input$removeMetaboliteFamily, {
  removeMetaboliteFamily <- as.numeric(input$removeMetaboliteFamily)
  print(paste("Observe removeMetaboliteFamily", removeMetaboliteFamily))
  
  ## TODO
  selectedAnnotationTableRowIdx <- isolate(input$familySelectionTable_rows_selected)
  annotationName  <- allAnnotationNames[[selectedAnnotationTableRowIdx]]
  
  annoThere <- lapply(X = dataList$annoArrayOfLists, FUN = function(x){ match(x = annotationName, table = x) })
  precursorSet <- which(!is.na(annoThere))
  
  removeAnnotation(precursorSet = precursorSet, annotationValue = annotationName)
})

annotationDialogCounter2 <- 1
renameMetaboliteFamily <- function(){
  selectedAnnotationTableRowIdx <- isolate(input$familySelectionTable_rows_selected)
  annotationName  <- allAnnotationNames[[selectedAnnotationTableRowIdx]]
  annotationIdx   <- which(dataList$annoPresentAnnotationsList == annotationName)
  annotationColor <- dataList$annoPresentColorsList[[annotationIdx]]
  
  callbackFunction <- function(value, color){
    print(paste("Selected annotation properties:", value, color))
    removeModal(session = session)
    
    ## add
    renameAnnotation(annotationValueOld = annotationName, annotationColorOld = annotationColor, annotationValueNew = value, annotationColorNew = color)
    #selectRows(proxy = dataTableProxy(outputId = "familySelectionTable", session = session), selected = NULL)
  }
  
  openAnnotaionNameColorDialog(predefinedClassName = annotationName, predefinedClassColor = annotationColor, callbackFunction = callbackFunction, confirmButtonText = "Rename annotation")
}

annotationDialogCounter <- 1
openAnnotaionNameColorDialog <- function(predefinedClassName, predefinedClassColor = colorPalette()[[1]], callbackFunction, confirmButtonText = "Add new annotation"){
  
  itemIdValue  <- paste("newAnnotationValue3",  annotationDialogCounter, sep = "_")
  itemIdColor  <- paste("newAnnotationColor3",  annotationDialogCounter, sep = "_")
  itemIdButton <- paste("submitNewAnnotation3", annotationDialogCounter, sep = "_")
  
  modalDialog <- modalDialog(title = "Add new annotation", 
                             #HTML("huhu"),
                             h4("Please specifiy name and color of the new annotation"),
                             bsTooltip(id = itemIdValue, title = "The name of this annotation", placement = "bottom", trigger = "hover"),
                             textInput(inputId = itemIdValue, placeholder = 'Metabolite family name here', label = "Type new annotation", value = predefinedClassName),
                             bsTooltip(id = itemIdColor, title = "The color of this annotation", placement = "bottom", trigger = "hover"),
                             colourInput(inputId = itemIdColor, label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = colorPalette(), value = predefinedClassColor),
                             bsTooltip(id = itemIdButton, title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                             actionButton(inputId = itemIdButton, label = confirmButtonText, class="btn-success")
  )
  
  lapply(c(itemIdButton), function(buttonId){
    observeEvent(input[[buttonId]], {
      value <- input[[itemIdValue]]
      color <- input[[itemIdColor]]
      
      callbackFunction(value, color)
    })
  })
  
  showUiDialog(modalDialog)
  
  annotationDialogCounter <<- annotationDialogCounter + 1
}

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabAnnotation observers")
  obsFamilySelectionTable_rows_selected$suspend()
  obsFragmentPlot2dblClick$suspend()
  obsClassSelection$suspend()
  obsSelectMetaboliteFamily$suspend()
  obsRenameMetaboliteFamily$suspend()
  obsRemoveMetaboliteFamily$suspend()
})
