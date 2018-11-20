
## selections
selectionAnalysisName <- "Selection by HCA/PCA"
selectionFragmentName <- "Selection by fragment"
selectionSearchName   <- "Selection by search"
precursorSelectionTabSelection  <- "Selection"
precursorSelectionTabAnnotation <- "Annotation"
precursorSelectionTabTable      <- "Table"
selectionAnalysisHcaName <- "Analysis_HCA"
selectionAnalysisPcaName <- "Analysis_PCA"
selectionFragmentHcaName <- "Fragment_HCA"
selectionFragmentPcaName <- "Fragment_PCA"
selectionSearchHcaName   <- "Search_HCA"
selectionSearchPcaName   <- "Search_PCA"

changeSelectionCurrentSelection <- selectionAnalysisName
precursorSelectionTabCurrentTab <- precursorSelectionTabSelection

## selection MS2
selectionFragmentSelectedFragmentIndex <- NULL
selectionFragmentTreeNodeSet <- NULL
selectionFragmentPcaLoadingSet <- NULL
## selection analysis
selectionAnalysisTreeNodeSet <- NULL
selectionAnalysisPcaLoadingSet <- NULL
## selection search
selectionSearchTreeNodeSet <- NULL
selectionSearchPcaLoadingSet <- NULL

## table data
ms1FeatureTableInputFieldIdCounter <- 0
selectedPrecursorSet <- NULL

selectedTable <- NULL
selectedTable_id <- NULL
table <- reactiveValues(
  df_Fragment_HCA = NULL,
  df_Fragment_PCA = NULL,
  df_Search_HCA = NULL,
  df_Search_PCA = NULL,
  df_Analysis_HCA = NULL,
  df_Analysis_PCA = NULL
)
listForTable_Fragment_HCA = NULL
listForTable_Fragment_PCA = NULL
listForTable_Analysis_HCA = NULL
listForTable_Analysis_PCA = NULL
listForTable_Search_HCA = NULL
listForTable_Search_PCA = NULL

table_Fragment_HCA_id = NULL
table_Fragment_PCA_id = NULL
table_Analysis_HCA_id = NULL
table_Analysis_PCA_id = NULL
table_Search_HCA_id = NULL
table_Search_PCA_id = NULL


state_selections <- reactiveValues(
  ## precursor set selections
  precursorSetSelected = FALSE,
  selectedSelection = NULL
)
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset selection state")
  state_selections$precursorSetSelected <<- FALSE
  state_selections$selectedSelection <<- NULL
  
  ## selection
  selectedPrecursorSet <<- NULL
  selectedTable <<- NULL
  
  selectionFragmentSelectedFragmentIndex <- NULL
  
  selectionFragmentTreeNodeSet <- NULL
  selectionAnalysisTreeNodeSet <- NULL
  selectionSearchTreeNodeSet <- NULL
  
  selectionFragmentPcaLoadingSet <- NULL
  selectionAnalysisPcaLoadingSet <- NULL
  selectionSearchPcaLoadingSet <- NULL
  
  selectionByFragmentReset()
  selectionByAnalysisReset()
  selectionBySearchReset()
})


selectionByFragmentReset <- function(){
  selectionFragmentSelectedFragmentIndex <<- NULL
  
  if(!is.null(selectionFragmentTreeNodeSet)){ ## HCA
    selectionFragmentSelectedFragmentIndex <<- NULL
    selectionFragmentTreeNodeSet <<- NULL
    listForTable_Fragment_HCA <<- NULL
    table_Fragment_HCA_id <<- NULL
    table$df_Fragment_HCA <<- NULL
    #output$dt_Fragment_HCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionFragmentHcaName)
      updateSelectedPrecursorSet()
  }
  if(!is.null(selectionFragmentPcaLoadingSet)){ ## PCA
    selectionFragmentPcaLoadingSet <<- NULL
    listForTable_Fragment_PCA <<- NULL
    table_Fragment_PCA_id <<- NULL
    table$df_Fragment_PCA <<- NULL
    #output$dt_Fragment_PCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionFragmentPcaName)
      updateSelectedPrecursorSet()
  }
}
selectionByFragment <- function(minimumIndex){
  selectionFragmentSelectedFragmentIndex <<- minimumIndex
  
  fragmentMass  <- ms2PlotValues$fragmentListClicked$fragmentMasses[[minimumIndex]]
  fragmentIndex <- which(dataList$fragmentMasses == fragmentMass)
  precursorSet  <- which(dataList$featureMatrix[, fragmentIndex] != 0)
  
  if(state$showHCAplotPanel)
    selectionByFragmentInitHca(precursorSet)
  if(state$showPCAplotPanel)
    selectionByFragmentInitPca(precursorSet)
  
  if(input$changeSelection != selectionFragmentName){
    updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionFragmentName)
    updateSelectedSelection()
  }
}
selectionByFragmentInitHca <- function(precursorSet){
  ## HCA - fetch subroots of subtrees comprising the selected fragment
  selectionFragmentTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
  precursorSetHca <- intersect(x = precursorSet, y = filterHca$filter)
  
  if(length(precursorSetHca) > 0){
    listForTable_Fragment_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetHca)
    table$df_Fragment_HCA <<- createMS1FeatureTable(listForTable_Fragment_HCA, selectionFragmentHcaName)
    table_Fragment_HCA_id <<- ms1FeatureTableInputFieldIdCounter
    #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
    if(state_selections$selectedSelection == selectionFragmentHcaName)
      updateSelectedPrecursorSet()
  } else {
    listForTable_Fragment_HCA <<- NULL
    table_Fragment_HCA_id <<- NULL
    table$df_Fragment_HCA <<- NULL
    #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
    if(state_selections$selectedSelection == selectionFragmentHcaName)
      updateSelectedPrecursorSet()
  }
}
selectionByFragmentInitPca <- function(precursorSet){
  ## PCA
  #selectionFragmentPcaLoadingSet <<- which(dataList$featureMatrix[filterPca$filter, fragmentIndex] != 0)
  selectionFragmentPcaLoadingSet <<- which(filterPca$filter %in% precursorSet)
  precursorSetPca <- intersect(x = precursorSet, y = filterPca$filter)
  
  if(length(precursorSetPca) > 0){
    listForTable_Fragment_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetPca)
    table$df_Fragment_PCA <<- createMS1FeatureTable(listForTable_Fragment_PCA, selectionFragmentPcaName)
    table_Fragment_PCA_id <<- ms1FeatureTableInputFieldIdCounter
    #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
    if(state_selections$selectedSelection == selectionFragmentPcaName)
      updateSelectedPrecursorSet()
  } else {
    listForTable_Fragment_PCA <<- NULL
    table_Fragment_PCA_id <<- NULL
    table$df_Fragment_PCA <<- NULL
    #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
    if(state_selections$selectedSelection == selectionFragmentPcaName)
      updateSelectedPrecursorSet()
  }
}

selectionByAnalysisReset <- function(){
  if(!is.null(selectionAnalysisTreeNodeSet)){ ## HCA
    selectionAnalysisTreeNodeSet <<- NULL
    listForTable_Analysis_HCA <<- NULL
    table_Analysis_HCA_id <<- NULL
    table$df_Analysis_HCA <<- NULL
    #output$dt_Analysis_HCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionAnalysisHcaName)
      updateSelectedPrecursorSet()
  }
  if(!is.null(selectionAnalysisPcaLoadingSet)){ ## PCA
    selectionAnalysisPcaLoadingSet <<- NULL
    listForTable_Analysis_PCA <<- NULL
    table_Analysis_PCA_id <<- NULL
    table$df_Analysis_PCA <<- NULL
    #output$dt_Analysis_PCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionAnalysisPcaName)
      updateSelectedPrecursorSet()
  }
}
selectionByAnalysisInitHca <- function(precursorSet){
  if(any(precursorSet %in% filterHca$filter)){
    #selectionAnalysisTreeNode <<- -match(x = precursorIndex, table = filterHca$filter)
    selectionAnalysisTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
    listForTable_Analysis_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
    table$df_Analysis_HCA <<- createMS1FeatureTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
    table_Analysis_HCA_id <<- ms1FeatureTableInputFieldIdCounter
    #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
    if(state_selections$selectedSelection == selectionAnalysisHcaName)
      updateSelectedPrecursorSet()
  } else {
    selectionAnalysisTreeNodeSet <<- NULL
    listForTable_Analysis_HCA <<- NULL
    table_Analysis_HCA_id <<- NULL
    table$df_Analysis_HCA <<- NULL
    #output$dt_Analysis_HCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionAnalysisHcaName)
      updateSelectedPrecursorSet()
  }
}
selectionByHca <- function(minimumLabel){
  selectionAnalysisTreeNodeSet <<- minimumLabel
  precursorSet <- getPrecursorSetFromTreeSelection(clusterDataList = clusterDataList, clusterLabel = minimumLabel)
  selectionByHca2(precursorSet)
}
selectionByHca2 <- function(precursorSet){
  listForTable_Analysis_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
  table$df_Analysis_HCA <<- createMS1FeatureTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
  #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
  if(state_selections$selectedSelection == selectionAnalysisHcaName)
    updateSelectedPrecursorSet()
  
  ## pca selection
  if(state$showPCAplotPanel)
    selectionByAnalysisInitPca(precursorSet)
  
  if(input$changeSelection != selectionAnalysisName){
    updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionAnalysisName)
    updateSelectedSelection()
  }
}
selectionByAnalysisInitPca <- function(precursorSet){
  precursorSetPca <- intersect(x = precursorSet, y = filterPca$filter)
  
  if(length(precursorSetPca) > 0){
    selectionAnalysisPcaLoadingSet <<- which(filterPca$filter %in% precursorSetPca)
    
    listForTable_Analysis_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetPca)
    table$df_Analysis_PCA <<- createMS1FeatureTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
    table_Analysis_PCA_id <<- ms1FeatureTableInputFieldIdCounter
    #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
    if(state_selections$selectedSelection == selectionAnalysisPcaName)
      updateSelectedPrecursorSet()
  } else {
    selectionAnalysisPcaLoadingSet <<- NULL
    listForTable_Analysis_PCA <<- NULL
    table_Analysis_PCA_id <<- NULL
    table$df_Analysis_PCA <<- NULL
    #output$dt_Analysis_PCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionAnalysisPcaName)
      updateSelectedPrecursorSet()
  }
}
selectionByPca <- function(minimumIndex){
  selectionAnalysisPcaLoadingSet <<- minimumIndex
  precursorIndex <- filterPca$filter[[minimumIndex]]
  selectionByPca2(precursorIndex)
}
selectionByPca2 <- function(precursorSet){
  listForTable_Analysis_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
  table$df_Analysis_PCA <<- createMS1FeatureTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
  #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
  if(state_selections$selectedSelection == selectionAnalysisPcaName)
    updateSelectedPrecursorSet()
  
  if(state$showHCAplotPanel)
    selectionByAnalysisInitHca(precursorSet)
  
  if(input$changeSelection != selectionAnalysisName){
    updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionAnalysisName)
    updateSelectedSelection()
  }
}

selectionBySearchReset <- function(){
  if(!is.null(selectionSearchTreeNodeSet)){ ## HCA
    selectionSearchTreeNodeSet <<- NULL
    listForTable_Search_HCA <<- NULL
    table_Search_HCA_id <<- NULL
    table$df_Search_HCA <<- NULL
    if(state_selections$selectedSelection == selectionSearchHcaName)
      updateSelectedPrecursorSet()
  }
  if(!is.null(selectionSearchPcaLoadingSet)){ ## HCA
    selectionSearchPcaLoadingSet <<- NULL
    listForTable_Search_PCA <<- NULL
    table_Search_PCA_id <<- NULL
    table$df_Search_PCA <<- NULL
    if(state_selections$selectedSelection == selectionSearchHcaName)
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
selectionBySearchInitHca <- function(precursorSet){
  precursorSetHca <- intersect(x = precursorSet, y = filterHca$filter)
  
  if(length(precursorSetHca) > 0){
    selectionSearchTreeNodeSet <<- getSetOfSubTreesFromRootForPrecursorSet(dataList = dataList, precursorSet = precursorSet, filter = filterHca$filter, clusterDataList = clusterDataList)
    listForTable_Search_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetHca)
    table$df_Search_HCA <<- createMS1FeatureTable(listForTable_Search_HCA, selectionSearchHcaName)
    table_Search_HCA_id <<- ms1FeatureTableInputFieldIdCounter
    #output$dt_Search_HCA <- DT::renderDataTable(table$df_Search_HCA)
    if(state_selections$selectedSelection == selectionSearchHcaName)
      updateSelectedPrecursorSet()
  } else {
    selectionSearchTreeNodeSet <<- NULL
    listForTable_Search_HCA <<- NULL
    table_Search_HCA_id <<- NULL
    table$df_Search_HCA <<- NULL
    #output$dt_Search_HCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionSearchHcaName)
      updateSelectedPrecursorSet()
  }
}
selectionBySearchInitPca <- function(precursorSet){
  precursorSetPca <- intersect(x = precursorSet, y = filterPca$filter)
  if(length(precursorSetPca) > 0){
    selectionSearchPcaLoadingSet <<- which(filterPca$filter %in% precursorSetPca)
    listForTable_Search_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSetPca)
    table$df_Search_PCA <<- createMS1FeatureTable(listForTable_Search_PCA, selectionSearchPcaName)
    table_Search_PCA_id <<- ms1FeatureTableInputFieldIdCounter
    #output$dt_Search_PCA <- DT::renderDataTable(table$df_Search_PCA)
    if(state_selections$selectedSelection == selectionSearchPcaName)
      updateSelectedPrecursorSet()
  } else {
    selectionSearchPcaLoadingSet <<- NULL
    listForTable_Search_PCA <<- NULL
    table_Search_PCA_id <<- NULL
    table$df_Search_PCA <<- NULL
    #output$dt_Search_PCA <- DT::renderDataTable(NULL)
    if(state_selections$selectedSelection == selectionSearchPcaName)
      updateSelectedPrecursorSet()
  }
}

createMS1FeatureTable <- function(list, type){
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
  
  ## checkboxes
  ms1FeatureTableInputFieldIdCounter <<- ms1FeatureTableInputFieldIdCounter + 1
  isArtifact <- dataList$annoArrayIsArtifact[precursorSet]
  #iconUp   <- icon(name = "chevron-up",   lib = "font-awesome")
  #iconDown <- icon(name = "chevron-down", lib = "font-awesome")
  checkboxes <- createCheckboxInputFields(       FUN = checkboxInput, id = paste(type, artifactName,   sep = "_"), values = isArtifact, tableCounter = ms1FeatureTableInputFieldIdCounter)
  #buttonUp   <- createActionButtonInputFields(   FUN = actionButton,  id = paste(type, "MoveUp",   sep = "_"), itemCount=length(isArtifact), icon   = iconUp, tableCounter = ms1FeatureTableInputFieldIdCounter)
  #buttonDown <- createActionButtonInputFields(   FUN = actionButton,  id = paste(type, "MoveDown", sep = "_"), itemCount=length(isArtifact), icon   = iconDown, tableCounter = ms1FeatureTableInputFieldIdCounter)
  #buttonUpDown <- createActionButtonInputFields2(FUN = actionButton,  id = paste(type, "Move", sep = "_"), itemCount=length(isArtifact), iconUp = iconUp, iconDown = iconDown, tableCounter = ms1FeatureTableInputFieldIdCounter)
  dataFrameIgnore <- data.frame(check.names = F,
                                Ignore = checkboxes
                                #"Move \u2191\u2193" = buttonUpDown
                                #"Change order" = buttonUpDown
  )
  
  ## assemble
  dataFrame <<- cbind(
    ms1abundanceDataFrame,
    dataFrameIgnore,
    annotationDataFrame,
    ms2fragmentDataFrame
  )
  
  return(dataFrame)
}
setMS1FeatureTable <- function(){
  output$ms1FeatureTable <- DT::renderDataTable(
    expr = selectedTable,
    server = FALSE, escape = FALSE, selection = "none",
    options = list(
      preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }')
    )
  )
}

updateMS1FeatureTableGui <- function(precursorSet){
  print("updateMS1FeatureTableGui")
  ## table update with new annotations
  if(all(length(precursorSet) > 0, !is.null(selectionFragmentTreeNodeSet), precursorSet %in% listForTable_Fragment_HCA$precursorSet)){ ## HCA
    table$df_Fragment_HCA <<- createMS1FeatureTable(listForTable_Fragment_HCA, selectionFragmentHcaName)
    #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
    table_Fragment_HCA_id <<- ms1FeatureTableInputFieldIdCounter
  }
  if(all(length(precursorSet) > 0, !is.null(selectionFragmentPcaLoadingSet), precursorSet %in% listForTable_Fragment_PCA$precursorSet)){ ## PCA
    table$df_Fragment_PCA <<- createMS1FeatureTable(listForTable_Fragment_PCA, selectionFragmentPcaName)
    #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
    table_Fragment_PCA_id <<- ms1FeatureTableInputFieldIdCounter
  }
  
  if(all(length(precursorSet) > 0, !is.null(selectionAnalysisTreeNodeSet), precursorSet %in% listForTable_Analysis_HCA$precursorSet)){ ## HCA
    table$df_Analysis_HCA <<- createMS1FeatureTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
    #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
    table_Analysis_HCA_id <<- ms1FeatureTableInputFieldIdCounter
  }
  if(all(length(precursorSet) > 0, !is.null(selectionAnalysisPcaLoadingSet), precursorSet %in% listForTable_Analysis_PCA$precursorSet)){ ## PCA
    table$df_Analysis_PCA <<- createMS1FeatureTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
    #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
    table_Analysis_PCA_id <<- ms1FeatureTableInputFieldIdCounter
  }
  
  if(all(length(precursorSet) > 0, !is.null(selectionSearchTreeNodeSet), precursorSet %in% listForTable_Search_HCA$precursorSet)){ ## HCA
    table$df_Search_HCA <<- createMS1FeatureTable(listForTable_Search_HCA, selectionSearchHcaName)
    #output$dt_Search_HCA <- DT::renderDataTable(table$df_Search_HCA)
    table_Search_HCA_id <<- ms1FeatureTableInputFieldIdCounter
  }
  if(all(length(precursorSet) > 0, !is.null(selectionSearchPcaLoadingSet), precursorSet %in% listForTable_Search_PCA$precursorSet)){ ## PCA
    table$df_Search_PCA <<- createMS1FeatureTable(listForTable_Search_PCA, selectionSearchPcaName)
    #output$dt_Search_PCA <- DT::renderDataTable(table$df_Search_PCA)
    table_Search_PCA_id <<- ms1FeatureTableInputFieldIdCounter
  }
  ## update
  updateTableAssignment()
}
updateTableAssignment <- function(){
  print(paste("updateTableAssignment '", state_selections$selectedSelection, "'", sep = ""))
  switch(state_selections$selectedSelection, 
         "Analysis_HCA"={ 
           #selectionAnalysisHcaName={ 
           selectedTable_id <<- table_Analysis_HCA_id
           selectedTable <<- table$df_Analysis_HCA
         },"Analysis_PCA"={  
           #},selectionAnalysisPcaName={  
           selectedTable_id <<- table_Analysis_PCA_id
           selectedTable <<- table$df_Analysis_PCA
         },"Fragment_HCA"={  
           #},selectionFragmentHcaName={  
           selectedTable_id <<- table_Fragment_HCA_id
           selectedTable <<- table$df_Fragment_HCA
           state_selections$precursorSetSelected <<- !is.null(listForTable_Fragment_HCA)
         },"Fragment_PCA"={  
           #},selectionFragmentPcaName={  
           selectedTable_id <<- table_Fragment_PCA_id
           selectedTable <<- table$df_Fragment_PCA
         },"Search_HCA"  ={  
           #},selectionSearchHcaName  ={  
           selectedTable_id <<- table_Search_HCA_id
           selectedTable <<- table$df_Search_HCA
         },"Search_PCA"  ={  
           #},selectionSearchPcaName  ={  
           selectedTable_id <<- table_Search_PCA_id
           selectedTable <<- table$df_Search_PCA
         },{
           print(paste("### unknown state_selections$selectedSelection: '", state_selections$selectedSelection, "'", sep = ""))
         }
  )
  setMS1FeatureTable()
}
updateSelectedSelection <- function(){
  selection <- input$changeSelection
  
  if(state$analysisType == "Annotation")
    return()
  
  if(selection == selectionAnalysisName & state$analysisType == "HCA")
    selectedSelection <- selectionAnalysisHcaName
  if(selection == selectionAnalysisName & state$analysisType == "PCA")
    selectedSelection <- selectionAnalysisPcaName
  if(selection == selectionFragmentName & state$analysisType == "HCA")
    selectedSelection <- selectionFragmentHcaName
  if(selection == selectionFragmentName & state$analysisType == "PCA")
    selectedSelection <- selectionFragmentPcaName
  if(selection == selectionSearchName & state$analysisType == "HCA")
    selectedSelection <- selectionSearchHcaName
  if(selection == selectionSearchName & state$analysisType == "PCA")
    selectedSelection <- selectionSearchPcaName
  
  state_selections$selectedSelection <<- selectedSelection
  updateSelectedPrecursorSet()
}
updateSelectedPrecursorSet <- function(){
  print(paste("updateSelectionAssignment '", state_selections$selectedSelection, "'", sep = ""))
  switch(state_selections$selectedSelection, 
         "Analysis_HCA"={ 
           state_selections$precursorSetSelected <<- !is.null(listForTable_Analysis_HCA)
           if(!is.null(listForTable_Analysis_HCA)) selectedPrecursorSet <<- listForTable_Analysis_HCA$precursorSet
           else                                    selectedPrecursorSet <<- NULL
         },"Analysis_PCA"={  
           state_selections$precursorSetSelected <<- !is.null(listForTable_Analysis_PCA)
           if(!is.null(listForTable_Analysis_PCA)) selectedPrecursorSet <<- listForTable_Analysis_PCA$precursorSet
           else                                    selectedPrecursorSet <<- NULL
         },"Fragment_HCA"={  
           state_selections$precursorSetSelected <<- !is.null(listForTable_Fragment_HCA)
           if(!is.null(listForTable_Fragment_HCA)) selectedPrecursorSet <<- listForTable_Fragment_HCA$precursorSet
           else                                    selectedPrecursorSet <<- NULL
         },"Fragment_PCA"={  
           state_selections$precursorSetSelected <<- !is.null(listForTable_Fragment_PCA)
           if(!is.null(listForTable_Fragment_PCA)) selectedPrecursorSet <<- listForTable_Fragment_PCA$precursorSet
           else                                    selectedPrecursorSet <<- NULL
         },"Search_HCA"  ={  
           state_selections$precursorSetSelected <<- !is.null(listForTable_Search_HCA)
           if(!is.null(listForTable_Search_HCA)) selectedPrecursorSet <<- listForTable_Search_HCA$precursorSet
           else                                    selectedPrecursorSet <<- NULL
         },"Search_PCA"  ={  
           state_selections$precursorSetSelected <<- !is.null(listForTable_Search_PCA)
           if(!is.null(listForTable_Search_PCA)) selectedPrecursorSet <<- listForTable_Search_PCA$precursorSet
           else                                    selectedPrecursorSet <<- NULL
         },{
           print(paste("### unknown state_selections$selectedSelection: '", state_selections$selectedSelection, "'", sep = ""))
         }
  )
  precursorSelectionChanged()
}
precursorSelectionChanged <- function(){
  selectionPresent <- !is.null(selectedPrecursorSet)
  
  ####################
  ## anno, table
  updateAnnoGui(selectedPrecursorSet)
  updateMS1FeatureTableGui(selectedPrecursorSet)
  
  ####################
  ## MetFrag link
  if(length(selectedPrecursorSet) == 1){
    metFragLinkList <- getMetFragLink(dataList, selectedPrecursorSet)
    
    if(is.null(metFragLinkList$error)){
      output$metFragLink <- renderText({
        print(paste("update output$metFragLink II", metFragLinkList$landingPageUrl))
        paste("<a href=", gsub(pattern = " ", replacement = "%20", x = metFragLinkList$landingPageUrl)[[1]], " target=\"_blank\">Send to MetFrag!</a>", sep = "")
      })
    } else {
      output$metFragLink <- renderText({
        print(paste("update output$metFragLink II", metFragLinkList$error))
        paste(metFragLinkList$error, sep = "")
      })
    }
  }
  
  ####################
  ## selection info
  selection <- state_selections$selectedSelection
  selectionInfo <- ""
  if(selectionPresent){
    switch(as.character(length(selectedPrecursorSet)), 
           "0"={ selectionInfo <- paste("The set of selected MS\u00B9 features is empty", sep = "")         },
           "1"={ selectionInfo <- paste(length(selectedPrecursorSet), " MS\u00B9 feature selected", sep = "")  },
           {     selectionInfo <- paste(length(selectedPrecursorSet), " MS\u00B9 features selected", sep = "") }
    )
  } else {
    if(selection == selectionAnalysisHcaName)
      selectionInfo <- paste("Please select a cluster or MS\u00B9 feature in the HCA plot", sep = "")
    if(selection == selectionAnalysisPcaName)
      selectionInfo <- paste("Please select a loading in the PCA plot", sep = "")
    if(selection == selectionFragmentHcaName | selection == selectionFragmentPcaName)
      selectionInfo <- paste("Please select a fragment in the Fragment plot above", sep = "")
    if(selection == selectionSearchHcaName | selection == selectionSearchPcaName)
      selectionInfo <- paste("Please select a set of MS\u00B9 features in the 'Search' tab of the sidebar panel", sep = "")
  }
  
  output$selectionInfo <- renderText({
    print(paste("update output$selectionInfo '", selectionInfo, "'", sep = ""))
    selectionInfo
  })
}
updatePlotsWithAnnotations <- function(){
  ## plots
  if(state$showHCAplotPanel){
    drawDendrogramPlot(consoleInfo = "updatePlotsWithAnnotations")
    #drawAnnotationLegendHCA(consoleInfo = "updatePlotsWithAnnotations")
  }
  if(state$showPCAplotPanel){
    drawPcaLoadingsPlot(consoleInfo = "updatePlotsWithAnnotations")
    #drawAnnotationLegendPCA(consoleInfo = "updatePlotsWithAnnotations")
  }
}

obsChangeSelection <- observeEvent(input$changeSelection, {
  selection <- input$changeSelection
  print(paste("Observe changeSelection", selection, "for", state$analysisType, ""))
  
  changeSelectionCurrentSelection <<- selection
  
  updateSelectedSelection()
})
obsPrecursorSelectionTabs <- observeEvent(input$precursorSelectionTabs, {
  selectedTab <- input$precursorSelectionTabs
  print(paste("Observe selectedTab", selectedTab))
  
  precursorSelectionTabCurrentTab <<- selectedTab
})
obsClearSelection <- observeEvent(input$clearSelection, {
  clearSelection  <- as.numeric(input$clearSelection)
  selection       <- input$changeSelection
  
  print(paste("Observe clearSelection", clearSelection))
  
  #################################################
  ## check if button was hit
  #if(clearSelection == clearSelectionButtonValue)
  #  return()
  #clearSelectionButtonValue <<- clearSelection
  
  switch(selection, 
         "Selection by HCA/PCA"={
           #selectionAnalysisName={
           selectionByAnalysisReset()
         },"Selection by fragment"={
           #},selectionFragmentName={
           selectionByFragmentReset()
         },"Selection by search"={
           #},selectionSearchName={
           selectionBySearchReset()
         },{
           print(paste("### unknown selection '", selection, "'", sep = ""))
         }
  )
  
  #################################################
  ## update plots
  if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "clear selection", withHeatmap = TRUE)
  if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "clear selection")
})

output$precursorSetSelected <- reactive({
  print(paste("reactive update precursorSetSelected", state_selections$precursorSetSelected))
  return(state_selections$precursorSetSelected)
})
output$selectedSelection <- reactive({
  print(paste("reactive update selectedSelection", state_selections$selectedSelection))
  return(state_selections$selectedSelection)
})

outputOptions(output, 'precursorSetSelected',    suspendWhenHidden=FALSE)
outputOptions(output, 'selectedSelection',       suspendWhenHidden=FALSE)

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending selections observers")
  obsChangeSelection$suspend()
  obsPrecursorSelectionTabs$suspend()
  obsClearSelection$suspend()
})
