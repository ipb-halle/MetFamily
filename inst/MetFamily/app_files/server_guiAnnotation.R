
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset annotation state")
  ## anno
  #updateColourInput(session = session, inputId = "newAnnotationColor", allowedCols = colorPalette())
  updateTextInput(session = session, inputId = "newAnnotationValue", value = "")
  updateTextInput(session = session, inputId = "newAnnotationValue2", value = "")
  updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[init]"))
  updateSelectInput(session = session, inputId = "previousAnnotationValue", choices = c("Artifact"))
})

checkAnnotationValue <- function(annotationValue){
  ## no ', '
  ## no '='
  ## not empty
  ## not only whitespaces (space, tab, newline, carriage return, vertical tab, form feed)
  return(!grepl(x = annotationValue, pattern = "(, )|(=)|(^$)|(^[ \t\n\r\v\f]+$)"))
}
checkAnnotationValueWithModalErrorDialog <- function(annotationValue){
  if(!checkAnnotationValue(annotationValue = annotationValue)){
    msg <- paste0(
      "The annotation name is not valid.", "<br>",
      "<br>",
      "The annotation name is not allowed to: <br> - contain comma ','",
      "<br> - contain '=' <br> - be empty <br> - only contain whitespaces",
      "<br><br>",
      "Please make sure that the above rules apply and try again."
    )
    showErrorDialog(msg)
    return(FALSE)
  }
  return(TRUE)
}
addAnnotation <- function(precursorSet, annotationValue, annotationColor){
  ## validate
  if(!checkAnnotationValueWithModalErrorDialog(annotationValue = annotationValue))
    return()
  
  if(is.na(match(x = annotationValue, table = dataList$annoPresentAnnotationsList))){
    ## new annotation
    annoIdx <- length(dataList$annoPresentAnnotationsList) + 1
    dataList$annoPresentAnnotationsList[[annoIdx]] <<- annotationValue
    dataList$annoPresentColorsList[[annoIdx]] <<- annotationColor
  }
  ## add
  for(precursor in precursorSet){
    if(annotationValue == artifactName)
      ## is artifact!
      dataList$annoArrayIsArtifact[[precursor]] <<- TRUE
    else{
      if(is.na(match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])))
        ## new anno
        dataList$annoArrayOfLists[[precursor]] <<- c(dataList$annoArrayOfLists[[precursor]], annotationValue)
    }
  }
  
  ## update gui
  #print("addAnnotation updateAnnoGui")
  checkForDoublyNamedAnnotations(annotationValue)
  updateAnnotationDependentGui(precursorSet)
  
  ##################################################################
  ## feedback
  msg <- paste(
    length(precursorSet), " MS\u00B9 features have been annotated as '", annotationValue, "'.",
    sep = ""
  )
  showInfoDialog(msg)
}
removeAnnotation <- function(precursorSet, annotationValue){
  ## remove
  for(precursor in precursorSet){
    if(annotationValue == artifactName)
      dataList$annoArrayIsArtifact[[precursor]] <<- FALSE
    else{
      idx <- match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])
      dataList$annoArrayOfLists[[precursor]] <<- dataList$annoArrayOfLists[[precursor]][-idx]
    }
  }
  
  ## remove anno completely?
  annoThere <- lapply(X = dataList$annoArrayOfLists, FUN = function(x){ match(x = annotationValue, table = x) })
  if(annotationValue != artifactName & all(is.na(annoThere))){
    idx <- match(x = annotationValue, table = dataList$annoPresentAnnotationsList)
    dataList$annoPresentAnnotationsList <<- dataList$annoPresentAnnotationsList[-idx]
    dataList$annoPresentColorsList      <<- dataList$annoPresentColorsList[-idx]
  }
  
  ## update gui
  #print("removeAnnotation updateAnnoGui")
  updateAnnotationDependentGui(precursorSet)
}
renameAnnotation <- function(annotationValueOld, annotationColorOld, annotationValueNew, annotationColorNew){
  
  # check name is valid
  if(!checkAnnotationValueWithModalErrorDialog(annotationValue = annotationValueNew)) return()
  
  ## update summary
  annotationIdx   <- which(dataList$annoPresentAnnotationsList == annotationValueOld)
  dataList$annoPresentAnnotationsList[[annotationIdx]] <<- annotationValueNew
  dataList$annoPresentColorsList     [[annotationIdx]] <<- annotationColorNew
  
  ## update feature ano
  precursorList <- list()
  for(featureIdx in seq_len(dataList$numberOfPrecursors)){
    if(annotationValueOld %in% dataList$annoArrayOfLists[[featureIdx]]){
      precursorList[[length(precursorList) + 1]] <- featureIdx
      dataList$annoArrayOfLists[[featureIdx]][[which(dataList$annoArrayOfLists[[featureIdx]] == annotationValueOld)]] <<- annotationValueNew
    }
  }
  
  checkForDoublyNamedAnnotations(annotationValueNew)
  updateAnnotationDependentGui(unlist(precursorList))
}
checkForDoublyNamedAnnotations <- function(annotationName){
  annotationIndeces   <- which(dataList$annoPresentAnnotationsList == annotationName)
  if(length(annotationIndeces) > 1){
    dataList$annoPresentAnnotationsList <<- dataList$annoPresentAnnotationsList[-annotationIndeces[-1]]
    dataList$annoPresentColorsList      <<- dataList$annoPresentColorsList[-annotationIndeces[-1]]
  }
}
setArtifactState <- function(precursorSet, isArtifact){
  ## add
  dataList$annoArrayIsArtifact[precursorSet] <<- isArtifact
  
  ## update gui
  updateAnnotationDependentGui(precursorSet)
}
setAnnotationPrimary <- function(precursorSet, annotationValue){
  ## remove
  for(precursor in precursorSet){
    idx <- match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])
    dataList$annoArrayOfLists[[precursor]] <<- c(dataList$annoArrayOfLists[[precursor]][[idx]], dataList$annoArrayOfLists[[precursor]][-idx])
  }
  
  ## update gui
  #print("setAnnotation primary updateAnnoGui")
  updateAnnotationDependentGui(precursorSet)
}
updateAnnotationDependentGui <- function(precursorSet){
  updateAnnoGui(precursorSet)
  updateMS1FeatureTableGui(precursorSet)
  updatePlotsWithAnnotations()
  updateAnnotationOverview()
}
commonAnnotations <- function(precursorSet){
  if(is.null(precursorSet))
    return(NULL)
  if(all(unlist(dataList$annoArrayIsArtifact[precursorSet])))
    return(artifactName)
  
  ## at least one non-artifact precursor present
  intersection <- unlist(dataList$annoPresentAnnotationsList)
  for(precursor in precursorSet)
    if(!dataList$annoArrayIsArtifact[[precursor]])
      intersection <- intersect(x = intersection, y = unlist(dataList$annoArrayOfLists[[precursor]]))
  return(intersection)
}
## table update
commonAnnotation <- character()
updateAnnoGui <- function(precursorSet){
  if(length(precursorSet) == 0)
    return()
  
  ## annotation
  commonAnnos <- commonAnnotations(precursorSet)
  
  if(length(commonAnnos) > 0){
    commonAnnotation <<- commonAnnos
    updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = commonAnnos)
    shinyjs::enable(id="presentAnnotationValue")
  } else {
    commonAnnotation <<- character()
    updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[none]"))
    shinyjs::disable(id="presentAnnotationValue")
  }
  
  updateSelectInput(session = session, inputId = "previousAnnotationValue", choices = dataList$annoPresentAnnotationsList)
  updateTabsetPanel(session = session, inputId = "precursorSelectionTabs", selected = precursorSelectionTabCurrentTab)
  
  if(length(commonAnnos) > 0){
    shinyjs::enable("removePresentAnnotation")
    shinyjs::enable("setPresentAnnotationPrimary")
  } else {
    shinyjs::disable("removePresentAnnotation")
    shinyjs::disable("setPresentAnnotationPrimary")
  }
}

## listen to annotation events
obsSetPresentAnnoPrimary <- observeEvent(input$setPresentAnnotationPrimary, {
  session$sendCustomMessage("disableButton", "setPresentAnnotationPrimary")
  value     <- input$presentAnnotationValue
  
  drawPlots <- as.numeric(input$setPresentAnnotationPrimary)
  #if(drawPlots == setPresentAnnotationPrimaryValue)
  #  return()
  #setPresentAnnotationPrimaryValue <<- drawPlots
  print(paste("Observe setPresentAnnotationPrimary", drawPlots))
  
  setAnnotationPrimary(precursorSet = selectedPrecursorSet, annotationValue = value)
  
  session$sendCustomMessage("enableButton", "setPresentAnnotationPrimary")
})
obsRemovePresentAnno <- observeEvent(input$removePresentAnnotation, {
  #session$sendCustomMessage("disableButton", "removePresentAnnotation")
  #session$sendCustomMessage(type="jsCode", list(code= paste("$('#","removePresentAnnotation","').prop('disabled',true)",sep="")))
  shinyjs::disable(id="removePresentAnnotation")
  
  value     <- input$presentAnnotationValue
  
  removePresentAnnotationValue <- as.numeric(input$removePresentAnnotation)
  ## XXX why?
  #if(drawPlots == removePresentAnnotationValue)
  #  return()
  #removePresentAnnotationValue <<- drawPlots
  print(paste("Observe removePresentAnnotation", removePresentAnnotationValue))
  
  removeAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value)
  #session$sendCustomMessage("enableButton", "removePresentAnnotation")
  #session$sendCustomMessage(type="jsCode", list(code= paste("$('#","removePresentAnnotation","').prop('enabled',true)",sep="")))
  
  #if(length(commonAnnotation) > 0)
  #  shinyjs::enable(id="removePresentAnnotation")
})#, priority = -1)
obsToggleAddNewAnnoButton <- observeEvent(input$newAnnotationValue, {
  value <- input$newAnnotationValue
  
  print(paste("Observe newAnnotationValue", nchar(value)))
  
  if(nchar(value) > 0)
    shinyjs::enable("submitNewAnnotation")
  #enableActionButton(session, "submitNewAnnotation")
  else
    shinyjs::disable("submitNewAnnotation")
  #disableActionButton(session, "submitNewAnnotation")
})
obsAddNewAnno <- observeEvent(input$submitNewAnnotation, {
  session$sendCustomMessage("disableButton", "submitNewAnnotation")
  value <- input$newAnnotationValue
  color <- input$newAnnotationColor
  
  submitNewAnnotation <- as.numeric(input$submitNewAnnotation)
  print(paste("Observe submitNewAnnotation", submitNewAnnotation))
  
  #if(submitNewAnnotation == submitNewAnnotationValue)
  #  return()
  #submitNewAnnotationValue <<- submitNewAnnotation
  
  addAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value, annotationColor = color)
  session$sendCustomMessage("enableButton", "submitNewAnnotation")
})
obsAddPresentAnno <- observeEvent(input$submitPreviousAnnotation, {
  session$sendCustomMessage("disableButton", "submitPreviousAnnotation")
  value <- input$previousAnnotationValue
  
  submitPreviousAnnotation <- as.numeric(input$submitPreviousAnnotation)
  #if(submitPreviousAnnotation == submitPreviousAnnotationValue)
  #  return()
  #submitPreviousAnnotationValue <<- submitPreviousAnnotation
  print(paste("Observe submitPreviousAnnotation", submitPreviousAnnotation))
  
  color <- dataList$annoPresentColorsList[[match(x = value, table = dataList$annoPresentAnnotationsList)]]
  addAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value, annotationColor = color)
  session$sendCustomMessage("enableButton", "submitPreviousAnnotation")
})
obsIgnoreValueChanged <- observeEvent(input$updateArtifactsFromCheckboxes, {
  session$sendCustomMessage("disableButton", "updateArtifactsFromCheckboxes")
  updateArtifactsFromCheckboxes <- as.numeric(input$updateArtifactsFromCheckboxes)
  
  #################################################
  ## check if button was hit
  #if(updateArtifactsFromCheckboxes == updateArtifactsFromCheckboxesButtonValue)
  #  return()
  #updateArtifactsFromCheckboxesButtonValue <<- updateArtifactsFromCheckboxes
  
  ## get and process values
  vals <- getInputValues(id = paste(state_selections$selectedSelection, artifactName, sep = "_"), counter = selectedTable_id, len = nrow(selectedTable))
  
  if(all(is.na(vals)))
    return()
  
  nas <- is.na(vals)
  vals[nas] <- dataList$annoArrayIsArtifact[selectedPrecursorSet][nas]
  vals <- as.logical(vals)
  
  setArtifactState(selectedPrecursorSet, vals)
  session$sendCustomMessage("enableButton", "updateArtifactsFromCheckboxes")
})

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabPca observers")
  obsSetPresentAnnoPrimary$suspend()
  obsRemovePresentAnno$suspend()
  obsToggleAddNewAnnoButton$suspend()
  obsAddNewAnno$suspend()
  obsAddPresentAnno$suspend()
  obsIgnoreValueChanged$suspend()
})
