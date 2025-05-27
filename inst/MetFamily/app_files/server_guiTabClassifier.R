
## classifier annotation
availableClassifiersDf <- NULL
availableClassifiersDfProperties <- NULL
classToSpectra_class <- NULL
mappingSpectraToClassDf <- NULL
properties_class <- NULL
selectedClassPrecursorIndeces <- NULL
selectedClassRowIdx <- -1
selectedClassFeatureRowIdx <- -1

putativeAnnotationsResultTableClassInputFieldIdCounter <- 1

specVsClassRange     <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
ms1FeatureVsClassTableCounter <- 0

state_tabClassifier <- reactiveValues(
  classifierLoaded = FALSE,
  classifierClassSelected = FALSE,
  classifierClassMS1featureSelected = FALSE
)
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset tabClassifier state")
  state_tabClassifier$classifierLoaded <<- FALSE
  state_tabClassifier$classifierClassSelected <<- FALSE
  state_tabClassifier$classifierClassMS1featureSelected <<- FALSE
  
  ## reset plot range
  specVsClassRange$xMin <<- NULL
  specVsClassRange$xMax <<- NULL
  
  ## classifier annotation
  classToSpectra_class <<- NULL
  mappingSpectraToClassDf <<- NULL
  properties_class <<- NULL
  selectedClassPrecursorIndeces <<- NULL
  selectedClassRowIdx <<- -1
  selectedClassFeatureRowIdx <<- -1
})

obsClassifierSelectionTable_rows_selected <- observeEvent(eventExpr = input$classifierSelectionTable_rows_selected, handlerExpr = {
  print(paste("Observe classifierSelectionTable_rows_selected", input$classifierSelectionTable_rows_selected))
  if(is.null(input$classifierSelectionTable_rows_selected)){
    session$sendCustomMessage("disableButton", "doAnnotation")
  } else {
    session$sendCustomMessage("enableButton", "doAnnotation")
  }
}, ignoreNULL = FALSE)
loadClassifier <- function(){
  #################################################
  ## do 999
  selectedRow <- isolate(expr = input$classifierSelectionTable_rows_selected)
  #if(is.null(selectedRow))
  #  selectedRow <- 1
  
  #selectedFile <- rownames(availableClassifiersDf)[[selectedRow]]
  #classifierProperties <- availableClassifiersDfProperties[selectedRow]
  filePath <- availableClassifiersDfProperties[selectedRow, "FilePath"]
  #method   <- availableClassifiersDfProperties[selectedRow, "algoName"]
  propertiesList <- as.list(availableClassifiersDfProperties[selectedRow,])
  
  withProgress(message = 'Calculating cluster...', value = 0, {
    resultObj <- doAnnotation(filePath, propertiesList, dataList$featureMatrix, dataList$importParameterSet, progress = TRUE)
  })
  
  classToSpectra_class    <<- resultObj$classToSpectra_class
  properties_class        <<- resultObj$properties_class
  mappingSpectraToClassDf <<- resultObj$mappingSpectraToClassDf
  state_tabClassifier$classifierLoaded <<- TRUE
  updateSelectInput(session = session, inputId = "metaboliteFamilyComparisonClass", choices = c(selectionNone, names(classToSpectra_class)), selected = selectionNone)
}
obsDoAnnotation <- observeEvent(input$doAnnotation, {
  session$sendCustomMessage("disableButton", "doAnnotation")
  doAnnotation <- as.numeric(input$doAnnotation)
  print(paste("Observe doAnnotation", doAnnotation))
  
  #################################################
  ## check if button was hit
  #if(doAnnotation == doAnnotationButtonValue)
  #  return()
  #doAnnotationButtonValue <<- doAnnotation
  
  loadClassifier()
  
  ## order by pValue
  bestPValues <- sapply(X = classToSpectra_class, FUN = head, n=1)
  
  if(FALSE){
    ## order by pValue
    order <- order(bestPValues, decreasing = F)
    classToSpectra_class <- classToSpectra_class[order]
    bestPValues          <- bestPValues[order]
    
    classToSpectra_class <<- classToSpectra_class
  }
  
  if(length(classToSpectra_class) == 0){
    output$noAnnotationsPopupDialog <- renderUI({
      shinyBS::bsModal(id = "noAnnotationsPopupDialog", title = "No annotations detected", trigger = "", size = "large",
              HTML("No MS\u00B9 feature could be annotated with one or more metabolite families.")
      )
    })
    shinyBS::toggleModal(session = session, modalId = "noAnnotationsPopupDialog", toggle = "open")
    
    return()
  }
  
  bestPValues <- format(bestPValues, digits=4)
  numPValues  <- sapply(X = classToSpectra_class, FUN = length)
  classes     <- names(classToSpectra_class)
  
  classesShort  <- classes
  
  chemOntClassNodeRegEx <- "[- ,/'0-9a-zA-Z]+"
  chemOntClassRegEx <- paste("^(", chemOntClassNodeRegEx, "; )*(", chemOntClassNodeRegEx, ")$", sep = "")
  
  ## shrink class names of the form "Organic compounds; Benzenoids; Benzene and substituted derivatives; Benzoic acids and derivatives; Halobenzoic acids and derivatives; 2-halobenzoic acids and derivatives"
  if(all(grepl(x = classes, pattern = chemOntClassRegEx))){
    maxLength <- 70
    tokenSeparator <- "; "
    abbreviator <- "."
    abbrevAbbrevs <- FALSE
    abbrevEndClass <- FALSE
    classesShort  <- sapply(X = seq_along(classes), FUN = function(x){
      len <- nchar(classes[[x]])
      if(len > maxLength){
        tokens   <- str_split(string = classes[[x]], pattern = tokenSeparator)[[1]]
        #abbrev   <- seq(from = 3+3+2, by = 3+3+2, length.out = length(tokens))
        abbrev   <- rep(3+3+2, times = length(tokens))
        noAbbrev <- nchar(tokens) + 2
        lengths  <- sapply(X = seq_len(length(tokens) - 1), FUN = function(y){
          sum(c(abbrev[1:y], noAbbrev[(y + 1):length(tokens)]))
        })
        if((lengths[[length(tokens) - 1]] > maxLength) & abbrevAbbrevs){
          if((3+3+3+1+2+nchar(tokens[[length(tokens)]]) > maxLength) & abbrevEndClass){
            classShort <- paste(c(
              paste(substr(x = c(tokens[[1]], tokens[[length(tokens) - 1]]), start = 1, stop = 3), abbreviator, sep = "", collapse = ".."), 
              paste("...", substr(x = tokens[[length(tokens)]], start = nchar(tokens[[length(tokens)]]) - (maxLength - (3+3+3+1+2+3) - 1), stop = nchar(tokens[[length(tokens)]])), sep = "")), 
              collapse = tokenSeparator)
          } else {
            classShort <- paste(c(
              paste(substr(x = c(tokens[[1]], tokens[[length(tokens) - 1]]), start = 1, stop = 3), abbreviator, sep = "", collapse = ".."), 
              tokens[[length(tokens)]]), 
              collapse = tokenSeparator)
          }
        } else {
          maxIdxToAbbrev <- ifelse(test = lengths[[length(tokens) - 1]] > maxLength, yes = length(tokens) - 1, no = min(which(lengths <= maxLength)))
          classShort <- paste(c(
            paste(substr(x = tokens[1:maxIdxToAbbrev], start = 1, stop = 3), abbreviator, sep = "", collapse = ""), 
            tokens[(maxIdxToAbbrev + 1):length(tokens)]), 
            collapse = tokenSeparator)
        }
      } else {
        classShort <- classes[[x]]
      }
      return(classShort)
    })
  }
  
  ## annotation column
  type <- "PutativeAnnotationsResultTableClass"
  annotateColumn <- createActionButtonInputFields(
    FUN = actionButton,  id = paste(type, "Annotate", sep = "_"), itemCount=length(classesShort), 
    label   = "Annotate", tableCounter = putativeAnnotationsResultTableClassInputFieldIdCounter, 
    callFunction = putativeAnnotationsResultTableClassAnnotateClicked
  )
  putativeAnnotationsResultTableClassInputFieldIdCounter <<- putativeAnnotationsResultTableClassInputFieldIdCounter + 1
  
  
  resultTable_classes <- data.frame(
    "Class"       = classesShort, 
    "Hits"        = numPValues,
    "Best pValue" = bestPValues,
    "Annotate"    = annotateColumn
    #row.names = rep("")
  )
  
  output$familyCount <- renderText({
    print(paste("init output$familyCount"))
    paste("Number of putative metabolite families: ", nrow(resultTable_classes),
          "\nNumber of putatively annotated features: ", length(unique(unlist(lapply(X = classToSpectra_class, FUN = names)))), " / ", dataList$numberOfPrecursors, 
          sep = "")
  })
  
  #selected <- ifelse(test = nrow(resultTable_classes) > 0, yes = 1, no = NULL)
  output$annotationResultTableClass <- DT::renderDataTable(
    expr = resultTable_classes,
    server = FALSE, escape = FALSE, rownames = FALSE,
    selection = list(mode = "single"),#, selected = selected), 
    options = list(
      #scrollY = "600px",
      scrollY = "30vh",
      preDrawCallback = DT::JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = DT::JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
      #rowCallback = JS(
      #  "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
      #  "var full_text = classes[[nRow]]",
      #  "$('td:eq(1)', nRow).attr('title', full_text);",
      #  "}"),
      iDisplayLength=nrow(resultTable_classes),       # initial number of records
      #aLengthMenu = c(5,10),    # records/page options
      #bLengthChange =0,        # show/hide records per page dropdown
      #bFilter = 0,              # global search box on/off
      #bInfo = 0,                # information on/off (how many records filtered, etc)
      #bAutoWidth = 0,           # automatic column width calculation, disable if passing column width via aoColumnDefs
      #aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
      ordering = F,              # row ordering
      sDom  = 'ft'
      #sDom  = '<"top">rt<"bottom">ip'
    )
  )
  
  ## state
  state$showAnnotationplotPanel <<- TRUE
  state$plotAnnotationShown <<- TRUE
  updateChangePlotRadioButton()
  
  #if(nrow(resultTable_classes) > 0)
  #  classSelected(1)
  
  session$sendCustomMessage("enableButton", "doAnnotation")
})
putativeAnnotationsResultTableClassAnnotateClicked <- function(buttonId){
  ## PutativeAnnotationsResultTableClass_Annotate_0_1
  rowIdx    <- as.integer(strsplit(x = buttonId, split = "_")[[1]][[4]])
  
  
  #selectedRow <- isolate(input$annotationResultTableClass_rows_selected)
  #print("#####")
  #print(selectedRow)
  #print(rowIdx)
  #if(length(selectedRow) == 0)
  #  selectRows(proxy = dataTableProxy("annotationResultTableClass"), selected = rowIdx)
  #else
  #  selectRows(proxy = dataTableProxy("annotationResultTableClass"), selected = NULL)
  
  
  precursorSet <- as.integer(names(classToSpectra_class[[rowIdx]]))
  class <- names(classToSpectra_class)[[rowIdx]]
  
  subClass <- tail(x = strsplit(x = class, split = "; ")[[1]], n = 1)
  print(paste(length(precursorSet), "-->", class, "(", subClass, ")", "from row", rowIdx))
  
  callbackFunction <- function(value, color){
    print(paste(length(precursorSet), "-->", class, "from row", rowIdx, "-->", value, color))
    removeModal(session = session)
    
    ## add
    addAnnotation(precursorSet = precursorSet, annotationValue = value, annotationColor = color)
  }
  
  openAnnotaionNameColorDialog(predefinedClassName = subClass, callbackFunction = callbackFunction)
}

obsAnnotationResultTableClass_selection <- observeEvent(ignoreNULL = FALSE, eventExpr = input$annotationResultTableClass_rows_selected, handlerExpr = {
  selectedRowIdx <- input$annotationResultTableClass_rows_selected
  print(paste("Selected class row:", selectedRowIdx))
  if(is.null(selectedRowIdx) && !initialGuiUpdatePerformed) return()
  
  if(is.null(selectedRowIdx)){
    state_tabClassifier$classifierClassSelected <<- FALSE
    return()
  }
  
  selectedClassRowIdx <<- selectedRowIdx
  state_tabClassifier$classifierClassSelected <<- TRUE
  
  classSelected(selectedRowIdx)
})
classSelected <- function(selectedRowIdx){
  ms1FeatureVsClassTableCounter <<- ms1FeatureVsClassTableCounter + 1
  
  #"classifierName",          #: chr "library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE"
  #"numberOfSpectra",         #: int 1355
  #"numberOfPositiveSpectra", #: int 20
  #"numberOfNegativeSpectra", #: int 1335
  #"class",                   #: chr "Organic compounds; Alkaloids and derivatives"
  #"fragmentMasses",          #: num [1:12207] -970 -969 -965 -965 -963 ...
  #"classOfClass",            #: chr "ChemOnt_SubstanceClass"
  #"frequentFragments",       #: Named num [1:29] 0.571 0.286 0.214 0.214 0.214 ...
  #"characteristicFragments", #: Named num [1:28] 0.539 0.275 0.203 0.203 0.143 ...
  #"AUC"                      #: num 0.912
  #"algoName"
  #"methodName"
  #"paramsString"
  
  classProperties         <- properties_class[[selectedRowIdx]]
  frequentFragments       <- classProperties$frequentFragments
  characteristicFragments <- classProperties$characteristicFragments
  
  class     <- names(classToSpectra_class)[[selectedRowIdx]]
  classToSpectra <- classToSpectra_class[[selectedRowIdx]]
  
  precursorIndecesAll <- as.integer(names(classToSpectra))
  pValuesAll <- format(unname(classToSpectra), digits=4)
  precursorIndeces <- precursorIndecesAll
  pValues <- pValuesAll
  
  maximumNumberOfDisplayedSpectrumHits <- 100
  if(length(precursorIndecesAll) > maximumNumberOfDisplayedSpectrumHits){
    precursorIndeces <- precursorIndeces[seq_len(maximumNumberOfDisplayedSpectrumHits)]
    pValues          <- pValues         [seq_len(maximumNumberOfDisplayedSpectrumHits)]
  }
  
  precursorLabels <- dataList$precursorLabels[precursorIndeces]
  selectedClassPrecursorIndeces <<- precursorIndeces
  
  
  output$classToSpectraCount <- renderText({
    print(paste("init output$classToSpectraCount"))
    paste("Number of putatively annotated features: ", ifelse(test = length(precursorIndecesAll) > maximumNumberOfDisplayedSpectrumHits, yes = paste(length(precursorIndeces), "/", length(precursorIndecesAll), sep = ""), no = length(precursorIndecesAll)),
          "\nEstimated true positive rate: ", format(classProperties$TPR_for_FPR_of_5Percent, digits=4), sep = "")
  })
  
  numberOfMatchingMasses_i <- getNumberOfHits(dataList, frequentFragments, characteristicFragments, precursorIndeces, mappingSpectraToClassDf, properties_class)
  
  returnObj <- createPlotOutput(
    dataList = dataList,
    id = "plotSpectrumVsClass", 
    frequentFragments = frequentFragments, 
    characteristicFragments = characteristicFragments, 
    precursorIndeces = precursorIndeces, 
    tableCounter = ms1FeatureVsClassTableCounter,
    mappingSpectraToClassDf = mappingSpectraToClassDf
  )
  
  numberOfMatchingMasses_i <- returnObj$numberOfMatchingMasses_i
  #inputs_i <- returnObj$inputs_i
  
  checkboxes   <- createCheckboxInputFields2(
    FUN = checkboxInput, 
    id = "MS1_feature_confirm", 
    values = rep(x = FALSE, times = length(precursorIndeces)), 
    tableCounter = ms1FeatureVsClassTableCounter,
    triggerSampleExclusionClick = FALSE
  )
  
  resultTable_features <- data.frame(check.names = FALSE,
                                     "Mz/RT"    = precursorLabels,
                                     "pValue"   = pValues,
                                     "Hits"     = numberOfMatchingMasses_i,
                                     "Confirm"  = checkboxes#,
                                     #"Spectrum" = inputs_i
  )
  selected <- ifelse(test = nrow(resultTable_features) > 0, yes = 1, no = NULL)
  output$annotationResultTableFeature <- DT::renderDataTable(
    expr = resultTable_features,
    server = FALSE, escape = FALSE, rownames = FALSE,
    selection = list(mode = "single"),#, selected = selected), 
    options = list(
      #scrollY = "600px",
      scrollY = "30vh",
      preDrawCallback = DT::JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = DT::JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
      #rowCallback = JS(
      #  "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
      #  "var full_text = classes[[nRow]]",
      #  "$('td:eq(1)', nRow).attr('title', full_text);",
      #  "}"),
      #selection = 'single',      # single row selection
      iDisplayLength=nrow(resultTable_features),       # initial number of records
      #aLengthMenu = c(5,10),    # records/page options
      #bLengthChange =0,        # show/hide records per page dropdown
      #bFilter = 0,              # global search box on/off
      #bInfo = 0,                # information on/off (how many records filtered, etc)
      #bAutoWidth = 0,           # automatic column width calculation, disable if passing column width via aoColumnDefs
      #aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
      autoWidth = TRUE,
      #columnDefs = list(list(width = '45px', targets = c(0,1,2))),#, list(width = '20px', targets = c(1))),
      columnDefs = list(
        list(width = '45px', targets = c(0)), 
        list(width = '20px', targets = c(1,2)),
        list(width = '15px', targets = c(3))
      ),
      ordering = F,              # row ordering
      sDom  = 't'
      #sDom  = '<"top">rt<"bottom">ip'
    )
  )
  
  ## update default class name for annotation
  #class <- names(classToSpectra_class)[[selectedClassRowIdx]]
  subClass <- tail(x = strsplit(x = class, split = "; ")[[1]], n = 1)
  updateTextInput(session = session, inputId = "newAnnotationValue2", value = subClass)
}
obsAnnotationResultTableFeature_selection <- observeEvent(ignoreNULL = FALSE, eventExpr = input$annotationResultTableFeature_rows_selected, handlerExpr = {
  selectedRowIdx <- input$annotationResultTableFeature_rows_selected
  print(paste("Selected feature row:", selectedRowIdx))
  
  if(is.null(selectedRowIdx)){
    state_tabClassifier$classifierClassMS1featureSelected <<- FALSE
    return()
  }
  
  selectedClassFeatureRowIdx <<- selectedRowIdx
  state_tabClassifier$classifierClassMS1featureSelected <<- TRUE
  
  classMS1FeatureSelected(selectedRowIdx)
})
classMS1FeatureSelected <- function(selectedRowIdx){
  classProperties         <- properties_class[[selectedClassRowIdx]]
  frequentFragments       <- classProperties$frequentFragments
  characteristicFragments <- classProperties$characteristicFragments
  
  precursorIndex = selectedClassPrecursorIndeces[[selectedRowIdx]]
  
  resetMS2VsClassPlotRange()
  
  drawMS2vsClassPlot(
    consoleInfo = "Feature selected", 
    frequentFragments=frequentFragments, 
    characteristicFragments=characteristicFragments, 
    precursorIndex=precursorIndex,
    mappingSpectraToClassDf=mappingSpectraToClassDf
  )
}
drawMS2vsClassPlot <- function(consoleInfo = NULL, frequentFragments, characteristicFragments, precursorIndex, mappingSpectraToClassDf){
  output$plotMS2vsClass <- renderPlot({
    print(paste("### SvC ###", consoleInfo))
    drawMS2vsClassPlotImpl(
      dataList=dataList,
      frequentFragments=frequentFragments, 
      characteristicFragments=characteristicFragments, 
      precursorIndex=precursorIndex,
      mappingSpectraToClassDf=mappingSpectraToClassDf
    )
  })
}
drawMS2vsClassPlotImpl <- function(dataList, frequentFragments, characteristicFragments, precursorIndex, mappingSpectraToClassDf){
  ## class statistics for class plot
  returnObj <- preprocessClassPlot(frequentFragments, characteristicFragments)
  masses_class    <- returnObj$masses_class
  frequency_class <- returnObj$frequency_class
  colors_class    <- returnObj$colors_class
  
  ## match spectrum masses for spectrum plot
  returnObj <- preprocessSpectrumVsClassPlot(dataList, precursorIndex, masses_class, mappingSpectraToClassDf, "Intensity")
  masses_spec     <- returnObj$masses_spec
  intensity_spec  <- returnObj$intensity_spec
  colors_spec     <- returnObj$colors_spec
  #numberOfMatchingMasses <- returnObj$numberOfMatchingMasses
  
  #xInterval <- c(dataList$minimumMass, dataList$maximumMass)
  
  calcPlotSpectrumVsClass_big(
    masses_spec    = masses_spec, 
    intensity_spec = intensity_spec, 
    colors_spec    = colors_spec, 
    masses_class   = masses_class, 
    frequency_class= frequency_class, 
    colors_class   = colors_class, 
    singleSpec     = TRUE,
    xInterval      = specVsClassRange$xInterval
  )
}
## plot range resets
resetMS2VsClassPlotRange <- function(){
  specVsClassRange$xMin <<- dataList$minimumMass
  specVsClassRange$xMax <<- dataList$maximumMass
  specVsClassRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
  specVsClassRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
}

obsApplyConfirmedAnnotations <- observeEvent(input$confirmAnnotation, {
  confirmAnnotation <- as.numeric(input$confirmAnnotation)
  
  #################################################
  ## check if button was hit
  #if(confirmAnnotation == confirmAnnotationButtonValue)
  #  return()
  #confirmAnnotationButtonValue <<- confirmAnnotation
  
  value <- input$newAnnotationValue2
  color <- input$newAnnotationColor2
  
  ## fetch data
  confirm <- getInputValues(id = paste("MS1_feature_confirm", sep = "_"), counter = ms1FeatureVsClassTableCounter, len = length(selectedClassPrecursorIndeces))
  precursorIndeces <- selectedClassPrecursorIndeces[confirm]
  print(paste("Annotate: Class ", value, " to ", length(precursorIndeces), " precursors (", paste(precursorIndeces, collapse = ", "), ")", sep = ""))
  
  ## apply annotation and update
  addAnnotation(precursorSet = precursorIndeces, annotationValue = value, annotationColor = color)
  ## updates e.g. plots automatically
})
obsToggleConfirmAnnoButton <- observeEvent(input$newAnnotationValue2, {
  value <- input$newAnnotationValue2
  
  print(paste("Observe newAnnotationValue2", nchar(value)))
  
  if(nchar(value) > 0)
    shinyjs::enable("confirmAnnotation")
  else
    shinyjs::disable("confirmAnnotation")
})

obsMS2VsClassdblClick <- observeEvent(input$plotMS2vsClass_dblclick, {
  brush <- input$plotMS2vsClass_brush
  print(paste("observe MS2vsClass dblclick", is.null(brush)))
  
  if (!is.null(brush)) {
    ## set range
    specVsClassRange$xMin <<- brush$xmin
    specVsClassRange$xMax <<- brush$xmax
    specVsClassRange$xInterval <<- c(brush$xmin, brush$xmax)
    specVsClassRange$xIntervalSize <<- brush$xmax - brush$xmin
  } else {
    ## reset range
    resetMS2VsClassPlotRange()
  }
})

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabClassifier observers")
  obsClassifierSelectionTable_rows_selected$suspend()
  obsDoAnnotation$suspend()
  obsAnnotationResultTableClass_selection$suspend()
  obsAnnotationResultTableFeature_selection$suspend()
  obsApplyConfirmedAnnotations$suspend()
  obsToggleConfirmAnnoButton$suspend()
  obsMS2VsClassdblClick$suspend()
})

output$classifierLoaded <- reactive({
  print(paste("reactive update classifierLoaded", state_tabClassifier$classifierLoaded))
  return(state_tabClassifier$classifierLoaded)
})
output$classifierClassSelected <- reactive({
  print(paste("reactive update classifierClassSelected", state_tabClassifier$classifierClassSelected))
  return(state_tabClassifier$classifierClassSelected)
})
output$classifierClassMS1featureSelected <- reactive({
  print(paste("reactive update classifierClassMS1featureSelected", state_tabClassifier$classifierClassMS1featureSelected))
  return(state_tabClassifier$classifierClassMS1featureSelected)
})
outputOptions(output, 'classifierLoaded',        suspendWhenHidden=FALSE)
outputOptions(output, 'classifierClassSelected', suspendWhenHidden=FALSE)
outputOptions(output, 'classifierClassMS1featureSelected',  suspendWhenHidden=FALSE)
