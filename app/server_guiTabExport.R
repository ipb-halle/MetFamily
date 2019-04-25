
resetWorkspaceFunctions <- c(resetWorkspaceFunctions, function(){
  print("Reset tabProject state")
  ## gui input fields
  shinyjs::toggleState("minimumIntensityOfMaximalMS2peak2", FALSE)
  shinyjs::toggleState("minimumProportionOfMS2peaks2", FALSE)
  shinyjs::toggleState("mzDeviationAbsolute_grouping2", FALSE)
  shinyjs::toggleState("mzDeviationInPPM_grouping2", FALSE)
  shinyjs::toggleState("doPrecursorDeisotoping2", FALSE)
  shinyjs::toggleState("mzDeviationAbsolute_precursorDeisotoping2", FALSE)
  shinyjs::toggleState("mzDeviationInPPM_precursorDeisotoping2", FALSE)
  shinyjs::toggleState("maximumRtDifference2", FALSE)
  shinyjs::toggleState("doMs2PeakGroupDeisotoping2", FALSE)
  shinyjs::toggleState("mzDeviationAbsolute_ms2PeakGroupDeisotoping2", FALSE)
  shinyjs::toggleState("mzDeviationInPPM_ms2PeakGroupDeisotoping2", FALSE)
  shinyjs::toggleState("neutralLossesPrecursorToFragments2", FALSE)
  shinyjs::toggleState("neutralLossesFragmentsToFragments2", FALSE)
  
  ## project infos
  updateTextInput    (session = session, inputId = "projectName2",           value = dataList$importParameterSet$projectName)
  shinyjs::toggleState("projectName2", FALSE)
  updateTextInput    (session = session, inputId = "projectDescription2",    value = dataList$importParameterSet$projectDescription)
  
  updateTextInput    (session = session, inputId = "minimumIntensityOfMaximalMS2peak2",            value = dataList$importParameterSet$minimumIntensityOfMaximalMS2peak)
  updateTextInput    (session = session, inputId = "minimumProportionOfMS2peaks2",                 value = dataList$importParameterSet$minimumProportionOfMS2peaks)
  updateTextInput    (session = session, inputId = "mzDeviationAbsolute_grouping2",                value = dataList$importParameterSet$mzDeviationAbsolute_grouping)
  updateTextInput    (session = session, inputId = "mzDeviationInPPM_grouping2",                   value = dataList$importParameterSet$mzDeviationInPPM_grouping)
  updateCheckboxInput(session = session, inputId = "doPrecursorDeisotoping2",                      value = dataList$importParameterSet$doPrecursorDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationAbsolute_precursorDeisotoping2",    value = dataList$importParameterSet$mzDeviationAbsolute_precursorDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationInPPM_precursorDeisotoping2",       value = dataList$importParameterSet$mzDeviationInPPM_precursorDeisotoping)
  updateTextInput    (session = session, inputId = "maximumRtDifference2",                         value = dataList$importParameterSet$maximumRtDifference)
  updateCheckboxInput(session = session, inputId = "doMs2PeakGroupDeisotoping2",                   value = dataList$importParameterSet$doMs2PeakGroupDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping2", value = dataList$importParameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping)
  updateTextInput    (session = session, inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping2",    value = dataList$importParameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping)
  #dataList$importParameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
  #dataList$importParameterSet$mzDeviationAbsolute_mapping
  #dataList$importParameterSet$minimumNumberOfMS2PeaksPerGroup
  updateCheckboxInput(session = session, inputId = "neutralLossesPrecursorToFragments2",           value = dataList$importParameterSet$neutralLossesPrecursorToFragments)
  updateCheckboxInput(session = session, inputId = "neutralLossesFragmentsToFragments2",           value = dataList$importParameterSet$neutralLossesFragmentsToFragments)
})

obsUpdateProjectDescription <- observeEvent(input$updateProjectDescription, {
  session$sendCustomMessage("disableButton", "updateProjectDescription")
  updateProjectDescription <- as.numeric(input$updateProjectDescription)
  
  print(paste("Observe updateProjectDescription", updateProjectDescription))
  
  #################################################
  ## check if button was hit
  #if(updateProjectDescription == updateProjectDescriptionButtonValue)
  #  return()
  #updateProjectDescriptionButtonValue <<- updateProjectDescription
  
  projectDescription <- input$projectDescription2
  projectDescription <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectDescription)))
  dataList$importParameterSet$projectDescription <<- projectDescription
  session$sendCustomMessage("enableButton", "updateProjectDescription")
})

## export image buttons
obsShowHCAplotPanel <- observe({
  print(paste("observe state$showHCAplotPanel", state$showHCAplotPanel))
  
  if(state$showHCAplotPanel){
    shinyjs::enable("downloadHcaImage")
  } else {
    shinyjs::disable("downloadHcaImage")
  }
})
obsShowPCAplotPanel <- observe({
  print(paste("observe state$showPCAplotPanel", state$showPCAplotPanel))
  
  if(state$showPCAplotPanel){
    shinyjs::enable("downloadPcaImage")
  } else {
    shinyjs::disable("downloadPcaImage")
  }
})

suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
  print("Suspending tabExport observers")
  obsShowHCAplotPanel$suspend()
  obsShowPCAplotPanel$suspend()
})
