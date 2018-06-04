
## data import: fixed parameters
proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- 0.9
mzDeviationAbsolute_mapping <- 0.01
#minimumNumberOfMS2PeaksPerGroup <- 1

## data
dataList <- NULL

enableLoadButtons <- function(){
  session$sendCustomMessage("enableButton", "loadProjectData")
  session$sendCustomMessage("enableButton", "loadExampleData")
  session$sendCustomMessage("enableButton", "importMs1Ms2Data")
  session$sendCustomMessage("enableButton", "importMs2Data")
}
disableLoadButtons <- function(){
  session$sendCustomMessage("disableButton", "loadProjectData")
  session$sendCustomMessage("disableButton", "loadExampleData")
  session$sendCustomMessage("disableButton", "importMs1Ms2Data")
  session$sendCustomMessage("disableButton", "importMs2Data")
}
obsFile <- observeEvent(input$matrixFile$datapath, {
  filePath <- input$matrixFile$datapath
  fileName <- input$matrixFile$name
  print(paste("Observe file for data", fileName))
  if(!is.null(filePath))
    shinyjs::enable("loadProjectData")
  
  updateFileInputInfo()
})
obsLoadProjectData <- observeEvent(input$loadProjectData, {
  disableLoadButtons()
  loadProjectData <- as.numeric(input$loadProjectData)
  print(paste("Observe loadProjectData", loadProjectData))
  
  #################################################
  ## check if button was hit
  #if(loadProjectData == loadProjectDataButtonValue)
  #  return()
  #loadProjectDataButtonValue <<- loadProjectData
  
  #################################################
  ## files
  filePath <- input$matrixFile$datapath
  loadProjectFile(filePath = filePath)
  enableLoadButtons()
})
obsLoadExampleData <- observeEvent(input$loadExampleData, {
  disableLoadButtons()
  loadExampleData <- as.numeric(input$loadExampleData)
  print(paste("Observe loadExampleData", loadExampleData))
  
  #################################################
  ## check if button was hit
  #if(loadExampleData == loadExampleDataButtonValue)
  #  return()
  #loadExampleDataButtonValue <<- loadExampleData
  
  #################################################
  ## files
  filePath <- getFile("Project_file_showcase_annotated.csv.gz")
  loadProjectFile(filePath = filePath)
  enableLoadButtons()
})
loadProjectFile <- function(filePath){
  fileName <- basename(filePath)
  #########################################################################################
  ## read data
  
  error <<- NULL
  withProgress(message = 'Reading file...', value = 0, {
    dataList <<- tryCatch(
      {
        readClusterDataFromProjectFile(file = filePath, progress = TRUE)
      }, 
      error = function(e) {
        print(e)
        error <<- e
      }
    )
  })
  
  if(!is.null(error)){
    print(paste("readClusterDataFromProjectFile resulted in error:", error))
    msg <- paste("An error occurred while reading the input files. Please check the file format and content and try again. The error was", error)
    output$fileInfo <- renderText({msg})
    #session$sendCustomMessage("enableButton", buttonId)
    #shinyBS::addPopover(session = session, id = "fileInputSelection", title = "Error", content = "huhu")
    
    msg <- paste(
      "An error occurred while reading the input files.",
      "Please check the file format and content and try again.",
      "The error was:",
      "<br>", 
      error
    )
    showErrorDialog(msg)
    
    return()
  }
  print(paste("readClusterDataFromProjectFile finished", dataList$minimumMass))
  
  resetWorkspace()
  
  state$importedOrLoadedFile_s_ <<- fileName
  updateFileInputInfo()
}
obsImportMs1DataFile <- observeEvent(input$ms1DataFile$datapath, {
  fileMs1Path <- input$ms1DataFile$datapath
  fileMs1Name <- input$ms1DataFile$name
  fileMs2Path <- input$ms2DataFile$datapath
  fileMs2Name <- input$ms2DataFile$name
  print(paste("Observe import MS1 file", fileMs1Name))
  
  if(all(!is.null(fileMs1Path), !is.null(fileMs2Path)))
    shinyjs::enable("importMs1Ms2Data")
  else
    shinyjs::disable("importMs1Ms2Data")
  
  updateFileInputInfo()
})
obsImportMs2DataFile <- observeEvent(input$ms2DataFile$datapath, {
  setImportState()
})
setImportState <- function(){
  fileMs1Path <- input$ms1DataFile$datapath
  fileMs1Name <- input$ms1DataFile$name
  fileMs2Path <- input$ms2DataFile$datapath
  fileMs2Name <- input$ms2DataFile$name
  print(paste("Observe import MS2 file", fileMs2Name))
  
  if(all(!is.null(fileMs1Path), !is.null(fileMs2Path)))
    shinyjs::enable("importMs1Ms2Data")
  else
    shinyjs::disable("importMs1Ms2Data")
  
  if(!is.null(fileMs2Path))
    shinyjs::enable("importMs2Data")
  else
    shinyjs::disable("importMs2Data")
  
  updateFileInputInfo()
}
obsImportMs1Ms2Data <- observeEvent(input$importMs1Ms2Data, {
  enableLoadButtons()
  importMs1Ms2Data <- as.numeric(input$importMs1Ms2Data)
  
  print(paste("Observe importMs1Ms2Data", importMs1Ms2Data))
  
  #################################################
  ## check if button was hit
  #if(importMs1Ms2Data == importMs1Ms2DataButtonValue)
  #  return()
  #importMs1Ms2DataButtonValue <<- importMs1Ms2Data
  
  importData(TRUE)
  disableLoadButtons()
})
obsImportMs2Data <- observeEvent(input$importMs2Data, {
  enableLoadButtons()
  importMs2Data <- as.numeric(input$importMs2Data)
  
  print(paste("Observe importMs2Data", importMs2Data))
  
  #################################################
  ## check if button was hit
  #if(importMs2Data == importMs2DataButtonValue)
  #  return()
  #importMs2DataButtonValue <<- importMs2Data
  
  importData(FALSE)
  disableLoadButtons()
})
importData <- function(importMS1andMS2data){
  #################################################
  ## files
  if(importMS1andMS2data){
    fileMs1Path <- input$ms1DataFile$datapath
    fileMs1Name <- input$ms1DataFile$name
  } else {
    fileMs1Path <- NULL
    fileMs1Name <- NULL
  }
  fileMs2Path <- input$ms2DataFile$datapath
  fileMs2Name <- input$ms2DataFile$name
  
  #################################################
  ## params
  
  ## project name
  projectName <- input$projectName
  projectName <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectName)))
  projectDescription <- input$projectDescription
  projectDescription <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectDescription)))
  
  ## minimum MS2 peak intensity
  minimumIntensityOfMaximalMS2peak <- input$minimumIntensityOfMaximalMS2peak
  minimumProportionOfMS2peaks <- input$minimumProportionOfMS2peaks
  ## grouping of MS2 peaks
  mzDeviationAbsolute_grouping <- input$mzDeviationAbsolute_grouping
  mzDeviationInPPM_grouping <- input$mzDeviationInPPM_grouping
  ## precursor deisotoping
  doPrecursorDeisotoping <- input$doPrecursorDeisotoping
  mzDeviationAbsolute_precursorDeisotoping <- input$mzDeviationAbsolute_precursorDeisotoping
  mzDeviationInPPM_precursorDeisotoping <- input$mzDeviationInPPM_precursorDeisotoping
  maximumRtDifference <- input$maximumRtDifference
  ## fragment deisotoping
  doMs2PeakGroupDeisotoping <- input$doMs2PeakGroupDeisotoping
  mzDeviationAbsolute_ms2PeakGroupDeisotoping <- input$mzDeviationAbsolute_ms2PeakGroupDeisotoping
  mzDeviationInPPM_ms2PeakGroupDeisotoping <- input$mzDeviationInPPM_ms2PeakGroupDeisotoping
  ## neutral losses
  neutralLossesPrecursorToFragments <- input$neutralLossesPrecursorToFragments
  neutralLossesFragmentsToFragments <- input$neutralLossesFragmentsToFragments
  #neutralLossesPrecursorToFragments <- TRUE
  #neutralLossesFragmentsToFragments <- FALSE
  
  ## fixed
  proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere <- proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
  mzDeviationAbsolute_mappingHere <- mzDeviationAbsolute_mapping
  #minimumNumberOfMS2PeaksPerGroupHere <- minimumNumberOfMS2PeaksPerGroup
  
  #################################################
  ## check params
  error <- FALSE
  if(any(is.null(minimumIntensityOfMaximalMS2peak), length(minimumIntensityOfMaximalMS2peak) == 0, nchar(minimumIntensityOfMaximalMS2peak) == 0))
    error <- TRUE
  else{
    minimumIntensityOfMaximalMS2peak <- as.numeric(minimumIntensityOfMaximalMS2peak)
    error <- error | is.na(minimumIntensityOfMaximalMS2peak)
  }
  if(any(is.null(minimumProportionOfMS2peaks), length(minimumProportionOfMS2peaks) == 0, nchar(minimumProportionOfMS2peaks) == 0))
    error <- TRUE
  else{
    minimumProportionOfMS2peaks <- as.numeric(minimumProportionOfMS2peaks)
    error <- error | is.na(minimumProportionOfMS2peaks)
  }
  if(any(is.null(mzDeviationAbsolute_grouping), length(mzDeviationAbsolute_grouping) == 0, nchar(mzDeviationAbsolute_grouping) == 0))
    error <- TRUE
  else{
    mzDeviationAbsolute_grouping <- as.numeric(mzDeviationAbsolute_grouping)
    error <- error | is.na(mzDeviationAbsolute_grouping)
  }
  if(any(is.null(mzDeviationInPPM_grouping), length(mzDeviationInPPM_grouping) == 0, nchar(mzDeviationInPPM_grouping) == 0))
    error <- TRUE
  else{
    mzDeviationInPPM_grouping <- as.numeric(mzDeviationInPPM_grouping)
    error <- error | is.na(mzDeviationInPPM_grouping)
  }
  if(doPrecursorDeisotoping){
    if(any(is.null(mzDeviationAbsolute_precursorDeisotoping), length(mzDeviationAbsolute_precursorDeisotoping) == 0, nchar(mzDeviationAbsolute_precursorDeisotoping) == 0))
      error <- TRUE
    else{
      mzDeviationAbsolute_precursorDeisotoping <- as.numeric(mzDeviationAbsolute_precursorDeisotoping)
      error <- error | is.na(mzDeviationAbsolute_precursorDeisotoping)
    }
    if(any(is.null(mzDeviationInPPM_precursorDeisotoping), length(mzDeviationInPPM_precursorDeisotoping) == 0, nchar(mzDeviationInPPM_precursorDeisotoping) == 0))
      error <- TRUE
    else{
      mzDeviationInPPM_precursorDeisotoping <- as.numeric(mzDeviationInPPM_precursorDeisotoping)
      error <- error | is.na(mzDeviationInPPM_precursorDeisotoping)
    }
  }
  if(any(is.null(maximumRtDifference), length(maximumRtDifference) == 0, nchar(maximumRtDifference) == 0))
    error <- TRUE
  else{
    maximumRtDifference <- as.numeric(maximumRtDifference)
    error <- error | is.na(maximumRtDifference)
  }
  if(doMs2PeakGroupDeisotoping){
    if(any(is.null(mzDeviationAbsolute_ms2PeakGroupDeisotoping), length(mzDeviationAbsolute_ms2PeakGroupDeisotoping) == 0, nchar(mzDeviationAbsolute_ms2PeakGroupDeisotoping) == 0))
      error <- TRUE
    else{
      mzDeviationAbsolute_ms2PeakGroupDeisotoping <- as.numeric(mzDeviationAbsolute_ms2PeakGroupDeisotoping)
      error <- error | is.na(mzDeviationAbsolute_ms2PeakGroupDeisotoping)
    }
    if(any(is.null(mzDeviationInPPM_ms2PeakGroupDeisotoping), length(mzDeviationInPPM_ms2PeakGroupDeisotoping) == 0, nchar(mzDeviationInPPM_ms2PeakGroupDeisotoping) == 0))
      error <- TRUE
    else{
      mzDeviationInPPM_ms2PeakGroupDeisotoping <- as.numeric(mzDeviationInPPM_ms2PeakGroupDeisotoping)
      error <- error | is.na(mzDeviationInPPM_ms2PeakGroupDeisotoping)
    }
  }
  
  if(any(is.null(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere), length(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere) == 0, nchar(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere) == 0))
    error <- TRUE
  else{
    proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere <- as.numeric(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere)
    error <- error | is.na(proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere)
  }
  if(any(is.null(mzDeviationAbsolute_mappingHere), length(mzDeviationAbsolute_mappingHere) == 0, nchar(mzDeviationAbsolute_mappingHere) == 0))
    error <- TRUE
  else{
    mzDeviationAbsolute_mappingHere <- as.numeric(mzDeviationAbsolute_mappingHere)
    error <- error | is.na(mzDeviationAbsolute_mappingHere)
  }
  # if(any(is.null(minimumNumberOfMS2PeaksPerGroupHere), length(minimumNumberOfMS2PeaksPerGroupHere) == 0, nchar(minimumNumberOfMS2PeaksPerGroupHere) == 0))
  #   error <- TRUE
  # else{
  #   minimumNumberOfMS2PeaksPerGroupHere <- as.numeric(minimumNumberOfMS2PeaksPerGroupHere)
  #   error <- error | is.na(minimumNumberOfMS2PeaksPerGroupHere)
  # }
  
  if(error){
    setImportState()
    output$fileInfo <- renderText({paste("There are invalid parameter values. Please check the parameters and press 'Import MS\u00B9 and MS/MS data' again.")})
    return()
  }
  
  ## box parameters
  print(paste("Observe importMs1Ms2Data", "e", error, "mi", minimumIntensityOfMaximalMS2peak, "mp", minimumProportionOfMS2peaks, "ga", mzDeviationAbsolute_grouping, "gr", mzDeviationInPPM_grouping, "pd", doPrecursorDeisotoping, "pa", mzDeviationAbsolute_precursorDeisotoping, "pr", mzDeviationInPPM_precursorDeisotoping, "mr", maximumRtDifference, "fd", doMs2PeakGroupDeisotoping, "fa", mzDeviationAbsolute_ms2PeakGroupDeisotoping, "fr", mzDeviationInPPM_ms2PeakGroupDeisotoping, "pm", proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere, "ma", mzDeviationAbsolute_mappingHere))
  parameterSet <- list()
  parameterSet$projectName                                       <- projectName
  parameterSet$projectDescription                                <- projectDescription
  parameterSet$toolVersion                                       <- paste(toolName, toolVersion, sep = " ")
  parameterSet$minimumIntensityOfMaximalMS2peak                  <- minimumIntensityOfMaximalMS2peak
  parameterSet$minimumProportionOfMS2peaks                       <- minimumProportionOfMS2peaks
  parameterSet$mzDeviationAbsolute_grouping                      <- mzDeviationAbsolute_grouping
  parameterSet$mzDeviationInPPM_grouping                         <- mzDeviationInPPM_grouping
  parameterSet$doPrecursorDeisotoping                            <- doPrecursorDeisotoping
  parameterSet$mzDeviationAbsolute_precursorDeisotoping          <- mzDeviationAbsolute_precursorDeisotoping
  parameterSet$mzDeviationInPPM_precursorDeisotoping             <- mzDeviationInPPM_precursorDeisotoping
  parameterSet$maximumRtDifference                               <- maximumRtDifference
  parameterSet$doMs2PeakGroupDeisotoping                         <- doMs2PeakGroupDeisotoping
  parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping       <- mzDeviationAbsolute_ms2PeakGroupDeisotoping
  parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping          <- mzDeviationInPPM_ms2PeakGroupDeisotoping
  parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- proportionOfMatchingPeaks_ms2PeakGroupDeisotopingHere
  parameterSet$mzDeviationAbsolute_mapping                       <- mzDeviationAbsolute_mappingHere
  #parameterSet$minimumNumberOfMS2PeaksPerGroup                   <- minimumNumberOfMS2PeaksPerGroupHere
  parameterSet$neutralLossesPrecursorToFragments                 <- neutralLossesPrecursorToFragments
  parameterSet$neutralLossesFragmentsToFragments                 <- neutralLossesFragmentsToFragments
  
  #################################################
  ## convert to project file
  
  ## built matrix
  error <- NULL
  withProgress(message = 'Generating matrix...', value = 0, {
    resultObj <- tryCatch(
      {
        convertToProjectFile(
          filePeakMatrix = fileMs1Path, 
          fileSpectra = fileMs2Path, 
          parameterSet = parameterSet, 
          progress = TRUE
        )
      }, error = function(e) {
        error <<- e
      }
    )
  })
  
  if(!is.null(error)){
    msg <- paste(
      "There occurred an error while processing the input file. Please check the file format and content and try again.", "\n",
      "Occurred error: ", error, sep = ""
    )
    output$fileInfo <- renderText({msg})
    showErrorDialog(msg)
    setImportState()
    return()
  }
  if(length(resultObj) == 1){
    if(resultObj == "Number of spectra is zero"){
      msg <- paste("There are no MS/MS spectra which fulfill the given criteria. Please refine parameter 'Spectrum intensity' and try 'Import MS\u00B9 and MS/MS data' again.")
      output$fileInfo <- renderText({msg})
      showErrorDialog(msg)
      setImportState()
      return()
    }
  }
  
  error <- NULL
  withProgress(message = 'Processing matrix...', value = 0, {
    lines <- sparseMatrixToString(matrixRows = resultObj$matrixRows, matrixCols = resultObj$matrixCols, matrixVals = resultObj$matrixVals, parameterSet = parameterSet)
    
    #################################################
    ## process project file
    
    dataList <<- tryCatch({
        readProjectData(fileLines = lines, progress = TRUE)
      }, error = function(e) {
        error <<- e
      }
    )
  })
  
  if(!is.null(error)){
    msg <- paste(
      "There occurred an error while processing the input file. Please check the file format and content and try again.", "\n",
      "Occurred error: ", error, sep = ""
    )
    output$fileInfo <- renderText({msg})
    showErrorDialog(msg)
    setImportState()
    return()
  }
  
  print(paste("readProjectData do data finished", dataList$minimumMass))
  
  spectraImport  <- paste(
    ## spectra
    resultObj$numberOfParsedSpectra, " / ", resultObj$numberOfSpectraOriginal, " spectra were imported successfully.",
    ifelse(test = resultObj$numberOfParsedSpectra < resultObj$numberOfSpectraOriginal, yes = paste(" (",paste( Filter(nchar, c(
      ifelse(test = resultObj$numberOfSpectraDiscardedDueToNoPeaks      > 0, yes = paste(resultObj$numberOfSpectraDiscardedDueToNoPeaks,      " empty", sep = ""), no = ""), 
      ifelse(test = resultObj$numberOfSpectraDiscardedDueToMaxIntensity > 0, yes = paste(resultObj$numberOfSpectraDiscardedDueToMaxIntensity, " low intensity", sep = ""), no = ""), 
      ifelse(test = resultObj$numberOfSpectraDiscardedDueToTooHeavy     > 0, yes = paste(resultObj$numberOfSpectraDiscardedDueToTooHeavy,     " too heavy", sep = ""), no = "")
    )), collapse = ", "), ")", sep = ""), no = ""),
    sep = ""
  )
  spectraMapping <- paste(
    ## mapping
    resultObj$numberOfPrecursors, " / ", resultObj$numberOfParsedSpectra, " spectra were successfully mapped to MS\u00B9 features.", 
    ifelse(test = resultObj$numberOfPrecursors < resultObj$numberOfParsedSpectra, yes = paste(" (",paste( Filter(nchar, c(
      #ifelse(test = resultObj$numberOfUnmappedPrecursorsMz > 0, yes = paste(resultObj$numberOfUnmappedPrecursorsMz, " with m/z deviation", sep = ""), no = ""), 
      #ifelse(test = resultObj$numberOfUnmappedPrecursorsRt > 0, yes = paste(resultObj$numberOfUnmappedPrecursorsRt, " with RT deviation",  sep = ""), no = "")
      ifelse(test = resultObj$numberOfUnmappedSpectra > 0, yes = paste(resultObj$numberOfUnmappedSpectra, " unmapped",  sep = ""), no = "")
    )), collapse = ", "), ")", sep = ""), no = ""),
    sep = ""
  )
  fragmentImport <- paste(
    ## fragments
    resultObj$numberOfMS2PeaksAboveThreshold, " / ", resultObj$numberOfMS2PeaksOriginal, " fragments were successfully imported.", 
    ifelse(test = resultObj$numberOfMS2PeaksAboveThreshold < resultObj$numberOfMS2PeaksOriginal, yes = paste(" (",paste( Filter(nchar, c(
      ifelse(test = resultObj$numberOfTooHeavyFragments      > 0, yes = paste(resultObj$numberOfTooHeavyFragments,      " too heavy",      sep = ""), no = ""), 
      ifelse(test = resultObj$numberOfMS2PeaksBelowThreshold > 0, yes = paste(resultObj$numberOfMS2PeaksBelowThreshold, " low intensity",  sep = ""), no = "")
    )), collapse = ", "), ")", sep = ""), no = ""),
    sep = ""
  )
  featureImport  <- paste(
    ## MS1 features
    resultObj$numberOfPrecursors, " / ", resultObj$numberOfParsedMs1Features, " MS\u00B9 features were successfully imported.",
    ifelse(test = resultObj$numberOfPrecursors < resultObj$numberOfParsedMs1Features, yes = paste(" (",paste( Filter(nchar, c(
      ifelse(test = resultObj$numberOfRemovedPrecursorIsotopePeaks > 0, yes = paste(resultObj$numberOfRemovedPrecursorIsotopePeaks, " were isotopes",   sep = ""), no = ""), 
      ifelse(test = resultObj$numberOfUnmappedPrecursors           > 0, yes = paste(resultObj$numberOfUnmappedPrecursors,           " without spectra", sep = ""), no = ""), 
      ifelse(test = resultObj$numberOfDuplicatedPrecursors         > 0, yes = paste(resultObj$numberOfDuplicatedPrecursors,         " duplicated",      sep = ""), no = "")
    )), collapse = ", "), ")", sep = ""), no = ""),
    sep = ""
  )
  
  msg <- paste(
    "The data import was successful.", "<br>",
    "<br>",
    spectraImport,  "<br>",
    ifelse(test = resultObj$numberOfParsedMs1Features!=-1, yes = paste(spectraMapping, "<br>", sep = ""), no = ""), 
    fragmentImport, "<br>",
    ifelse(test = resultObj$numberOfParsedMs1Features!=-1, yes = featureImport,  no = ""),
    sep = ""
  )
  showInfoDialog(msg)
  
  ## MS2
  # + returnObj$numberOfSpectraOriginal
  # + returnObj$numberOfMS2PeaksOriginal
  # - returnObj$numberOfMS2PeaksWithNeutralLosses
  # + returnObj$numberOfMS2PeaksAboveThreshold
  # + returnObj$numberOfMS2PeaksBelowThreshold
  # + returnObj$numberOfTooHeavyFragments
  # + returnObj$numberOfSpectraDiscardedDueToNoPeaks
  # + returnObj$numberOfSpectraDiscardedDueToMaxIntensity
  # + returnObj$numberOfSpectraDiscardedDueToTooHeavy
  #
  ## MS1
  # + returnObj$numberOfPrecursors
  # 
  # + returnObj$numberOfDuplicatedPrecursors
  # + returnObj$numberOfUnmappedPrecursors
  # + returnObj$numberOfUnmappedPrecursorsMz
  # + returnObj$numberOfUnmappedPrecursorsRt
  # + returnObj$numberOfParsedSpectra
  # + returnObj$numberOfParsedMs1Features
  # + returnObj$numberOfRemovedPrecursorIsotopePeaks
  
  
  resetWorkspace()
  
  if(importMS1andMS2data)
    state$importedOrLoadedFile_s_ <<- c(fileMs1Name, fileMs2Name)
  else
    state$importedOrLoadedFile_s_ <<- c(fileMs2Name)
  updateFileInputInfo()
  
  setImportState()
}
obsFileInputSelection <- observeEvent(input$fileInputSelection, {
  updateFileInputInfo()
})
obsApplyImportParameterFile <- observeEvent(input$importParameterFileInput$datapath, {
  filePath <- input$importParameterFileInput$datapath
  fileName <- input$importParameterFileInput$name
  print(paste("Observe importParameterFile", fileName))
  if(is.null(filePath))
    return()
  
  ## read and parse
  fileContent <- readLines(con = filePath)
  parameterSet <- deserializeParameterSetFile(fileContent)
  
  ## apply
  #projectName2 <- parameterSet$projectName
  #projectName2 <- paste(projectName2, " adopted", sep = "")
  #updateTextInput    (session = session, inputId = "projectName",                                 value = projectName2)
  #parameterSet$toolVersion
  updateTextInput    (session = session, inputId = "minimumIntensityOfMaximalMS2peak",            value = parameterSet$minimumIntensityOfMaximalMS2peak)
  updateTextInput    (session = session, inputId = "minimumProportionOfMS2peaks",                 value = parameterSet$minimumProportionOfMS2peaks)
  updateTextInput    (session = session, inputId = "mzDeviationAbsolute_grouping",                value = parameterSet$mzDeviationAbsolute_grouping)
  updateTextInput    (session = session, inputId = "mzDeviationInPPM_grouping",                   value = parameterSet$mzDeviationInPPM_grouping)
  updateCheckboxInput(session = session, inputId = "doPrecursorDeisotoping",                      value = parameterSet$doPrecursorDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationAbsolute_precursorDeisotoping",    value = parameterSet$mzDeviationAbsolute_precursorDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationInPPM_precursorDeisotoping",       value = parameterSet$mzDeviationInPPM_precursorDeisotoping)
  updateTextInput    (session = session, inputId = "maximumRtDifference",                         value = parameterSet$maximumRtDifference)
  updateCheckboxInput(session = session, inputId = "doMs2PeakGroupDeisotoping",                   value = parameterSet$doMs2PeakGroupDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping", value = parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping",    value = parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping)
  #parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
  #parameterSet$mzDeviationAbsolute_mapping
  #parameterSet$minimumNumberOfMS2PeaksPerGroup
  updateCheckboxInput(session = session, inputId = "neutralLossesPrecursorToFragments",           value = parameterSet$neutralLossesPrecursorToFragments)
  updateCheckboxInput(session = session, inputId = "neutralLossesFragmentsToFragments",           value = parameterSet$neutralLossesFragmentsToFragments)
})
updateFileInputInfo <- function(){
  fileInputSelection <- input$fileInputSelection
  filePath <- input$matrixFile$datapath
  fileName <- input$matrixFile$name
  fileMs1Path <- input$ms1DataFile$datapath
  fileMs1Name <- input$ms1DataFile$name
  fileMs2Path <- input$ms2DataFile$datapath
  fileMs2Name <- input$ms2DataFile$name
  #exampleDataSelection <- input$exampleDataSelection
  
  if(all(fileInputSelection == "Example data"))
    output$fileInfo <- renderText({paste("Please press 'Load example data' to load the full example data set")})
  if(all(fileInputSelection == "Load project", is.null(filePath)))
    output$fileInfo <- renderText({paste("Please select a project file and press 'Load project data'")})
  if(all(fileInputSelection == "Load project", !is.null(filePath), any(is.null(state$importedOrLoadedFile_s_), fileName != state$importedOrLoadedFile_s_)))
    output$fileInfo <- renderText({paste("Please press 'Load project data'")})
  if(all(fileInputSelection == "Load project", !is.null(filePath), !is.null(state$importedOrLoadedFile_s_), fileName == state$importedOrLoadedFile_s_))
    output$fileInfo <- renderText({paste(fileName)})
  if(all(fileInputSelection == "Import data", is.null(fileMs1Path), is.null(fileMs2Path)))
    output$fileInfo <- renderText({paste("Please select a metabolite profile and a MS/MS library and press 'Import MS1 and MS/MS data' or select a MS/MS library and press 'Import MS/MS data'")})
  if(all(fileInputSelection == "Import data", is.null(fileMs1Path), !is.null(fileMs2Path)))
    output$fileInfo <- renderText({paste("Please press 'Import MS/MS data' or select a metabolite profile and press 'Import MS1 and MS/MS data'")})
  if(all(fileInputSelection == "Import data", !is.null(fileMs1Path), is.null(fileMs2Path)))
    output$fileInfo <- renderText({paste("Please select a MS/MS library and press 'Import MS1 and MS/MS data'")})
  if(all(fileInputSelection == "Import data", !is.null(fileMs1Path), !is.null(fileMs2Path), any(is.null(state$importedOrLoadedFile_s_), c(fileMs1Name,fileMs2Name) != state$importedOrLoadedFile_s_)))
    output$fileInfo <- renderText({paste("Please press 'Import MS1 and MS/MS data'")})
  if(all(fileInputSelection == "Import data", !is.null(fileMs1Path), !is.null(fileMs2Path), !is.null(state$importedOrLoadedFile_s_), c(fileMs1Name,fileMs2Name) == state$importedOrLoadedFile_s_))
    output$fileInfo <- renderText({paste(fileMs1Name, "\n", fileMs2Name, sep = "")})
}

output$fileInfo <- renderText({
  print(paste("init output$fileInfo"))
  paste("Please select a project file and press 'Load project data'")
})