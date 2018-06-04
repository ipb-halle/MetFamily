
## for saving the table state while rearranging
sampleOrder_tmp <- NULL
sampleExclusion_tmp <- NULL

sampleTableInputFieldIdCounter <- 0
sampleTable <- NULL

createSampleTable <- function(){
  sampleTableInputFieldIdCounter <<- sampleTableInputFieldIdCounter + 1
  
  groupSampleDataFrame <- dataList$groupSampleDataFrame[sampleOrder_tmp, c("Group", "Sample")]
  isExcluded <- sampleExclusion_tmp[sampleOrder_tmp]
  
  iconUp   <- icon(name = "chevron-up",   lib = "font-awesome")
  iconDown <- icon(name = "chevron-down", lib = "font-awesome")
  
  checkboxes   <- createCheckboxInputFields2(    FUN = checkboxInput, id = "Sample_Exclude", values = isExcluded, tableCounter = sampleTableInputFieldIdCounter, triggerSampleExclusionClick = TRUE)
  buttonUpDown <- createActionButtonInputFields2(FUN = actionButton,  id = "Sample_Move", itemCount = length(isExcluded), iconUp = iconUp, iconDown = iconDown, tableCounter = sampleTableInputFieldIdCounter)
  dataFrameExcludeMove <- data.frame(check.names = F,
                                     "Exclude" = checkboxes,
                                     #"Move \u2191\u2193" = buttonUpDown
                                     "Change order" = buttonUpDown
  )
  
  ## assemble
  dataFrame <<- cbind(
    groupSampleDataFrame,
    dataFrameExcludeMove
  )
  
  return(dataFrame)
}
obsSampleValuesChanged <- observeEvent(input$updateSampleTable, {
  session$sendCustomMessage("disableButton", "updateSampleTable")
  updateSampleTable <- as.numeric(input$updateSampleTable)
  
  #################################################
  ## check if button was hit
  #if(updateSampleTable == updateSamplesButtonValue)
  #  return()
  #updateSamplesButtonValue <<- updateSampleTable
  
  updateSampleOrderAndExclusion()
  session$sendCustomMessage("enableButton", "updateSampleTable")
})
updateSampleOrderAndExclusion <- function(){
  ## check list validity
  if(length(unique(dataList$groupSampleDataFrame[!sampleExclusion_tmp, "Group"])) < length(unique(dataList$groupSampleDataFrame[, "Group"]))){
    showModal(modalDialog(
      title = "Invalid sample exclusion",
      #footer = tagList(
      #  actionButton(inputId = "ok", label = "OK")
      #),
      HTML(
        paste(
          "It is not supported to exclude all samples of one sample group.",
          "<br>",
          "Please retain at least one sample of each sample group."
        )
      )
    ))
    
    return()
  }
  
  ## update dataList
  dataList$groupSampleDataFrame[, "Order"]   <<- sampleOrder_tmp
  dataList$groupSampleDataFrame[, "Exclude"] <<- sampleExclusion_tmp
  
  sampleNamesToExclude <- dataList$groupSampleDataFrame[sampleExclusion_tmp, "Sample"]
  
  returnObj <- processMS1data(sampleNamesToExclude=sampleNamesToExclude, numberOfMS1features=dataList$numberOfPrecursors, precursorLabels=dataList$precursorLabels, 
                              groups=dataList$groups, metaboliteProfileColumnNames=dataList$metaboliteProfileColumnNames, tagsSector = dataList$tagsSector, 
                              dataColumnIndecesFunctionFromGroupIndex=dataList$dataColumnIndecesFunctionFromGroupIndex, dataColumnsNameFunctionFromGroupIndex=dataList$dataColumnsNameFunctionFromGroupIndex, dataColumnsNameFunctionFromGroupName=dataList$dataColumnsNameFunctionFromGroupName, dataColumnsNameFunctionFromGroupNames=dataList$dataColumnsNameFunctionFromGroupNames, groupNameFunctionFromDataColumnName=dataList$groupNameFunctionFromDataColumnName,
                              metaboliteProfile=dataList$dataFrameInfos, progress=FALSE)
  
  ## name functions
  dataList$dataMeanColumnNameFunctionFromIndex <<- returnObj$dataMeanColumnNameFunctionFromIndex
  dataList$dataMeanColumnNameFunctionFromName  <<- returnObj$dataMeanColumnNameFunctionFromName
  dataList$lfcColumnNameFunctionFromIndex      <<- returnObj$lfcColumnNameFunctionFromIndex
  dataList$lfcColumnNameFunctionFromName       <<- returnObj$lfcColumnNameFunctionFromName
  dataList$groupNameFromGroupIndex             <<- returnObj$groupNameFromGroupIndex
  dataList$groupIdxFromGroupName               <<- returnObj$groupIdxFromGroupName
  ## data and names
  dataList$dataFrameMeasurements               <<- returnObj$dataFrameMeasurements
  #dataMeanColumnNames                 <- returnObj$dataMeanColumnNames
  #lfcColumnNames                      <- returnObj$lfcColumnNames
  ## colors
  #matrixDataFrame                     <- returnObj$matrixDataFrame
  dataList$colorMatrixDataFrame                <<- returnObj$colorMatrixDataFrame
  dataList$colorMapAbsoluteData                <<- returnObj$colorMapAbsoluteData
  dataList$colorMapLogFoldChange               <<- returnObj$colorMapLogFoldChange
  dataList$columnGroupLabels                   <<- returnObj$columnGroupLabels
  ## constants
  dataList$meanAllMax       <<- returnObj$meanAllMax
  dataList$logFoldChangeMax <<- returnObj$logFoldChangeMax
  dataList$logAbsMax        <<- returnObj$logAbsMax
  
  ## update filter
  sampleNames <- dataList$groupSampleDataFrame[sampleOrder_tmp, "Sample"][!sampleExclusion_tmp[sampleOrder_tmp]]
  #selectedSampleNames <- intersect(input$pcaSamples, sampleNames)
  #updateCheckboxGroupInput(session = session, inputId = "pcaSamples",  choices = sampleNames,     selected = selectedSampleNames)
  updateCheckboxGroupInput(session = session, inputId = "pcaSamples",  choices = sampleNames,     selected = sampleNames)
}
sampleExcludeClicked <- function(){
  ## get exclusion
  sampleExclude <- getInputValues(id = paste("Sample_Exclude", sep = "_"), counter = sampleTableInputFieldIdCounter, len = nrow(sampleTable))
  sampleExclude <- sampleExclude[order(sampleOrder_tmp)]
  
  ## update exclusion
  sampleExclusion_tmp <<- sampleExclude
}
sampleMoveClicked <- function(buttonId){
  ## get source
  buttonUp <-            strsplit(x = buttonId, split = "_")[[1]][[3]] == "Up"
  rowId    <- as.numeric(strsplit(x = buttonId, split = "_")[[1]][[5]])
  
  ## ignore two cases
  if(rowId == 1 & buttonUp)
    return()
  if(rowId == nrow(dataList$groupSampleDataFrame) & !buttonUp)
    return()
  
  ## update sample order
  rowId2 <- ifelse(buttonUp, rowId - 1, rowId + 1)
  sampleOrder <- sampleOrder_tmp
  
  ## switch
  tmp <- sampleOrder[[rowId]]
  sampleOrder[[rowId]] <- sampleOrder[[rowId2]]
  sampleOrder[[rowId2]] <- tmp
  
  sampleOrder_tmp <<- sampleOrder
  
  ## update table
  sampleTable <<- createSampleTable()
  setSampleTable()
}
setSampleTable <- function(){
  output$sampleTable <- DT::renderDataTable(
    expr = sampleTable,
    server = FALSE, escape = FALSE, selection = "none", #rownames = FALSE,
    options = list(
      #scrollY = "600px",
      scrollY = "65vh",
      preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
      drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
      iDisplayLength=nrow(sampleTable),       # initial number of records
      #aLengthMenu = c(5,10),    # records/page options
      #bLengthChange =0,        # show/hide records per page dropdown
      #bFilter = 0,              # global search box on/off
      #bInfo = 0,                # information on/off (how many records filtered, etc)
      #bAutoWidth = 0,           # automatic column width calculation, disable if passing column width via aoColumnDefs
      #aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
      ordering = F,              # row ordering
      sDom  = 't'
      #sDom  = '<"top">rt<"bottom">ip'
    )
  )
}
