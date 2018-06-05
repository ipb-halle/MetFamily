
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
