
options(shiny.sanitize.errors = FALSE)

#########################################################################################
#########################################################################################
## libraries and functions

sourceFolder <- getwd()
#isDevelopment <- TRUE
isDevelopment <- grepl(pattern = "htreutle", x = sourceFolder) | grepl(pattern = "Treutler", x = sourceFolder)
errorHunting <- FALSE
hcaHeatMapNew <- TRUE

#####################################################################################################
## handling of errors and warnings
if(errorHunting){
  #options(warn = 2)
  options(warn = 2, shiny.error = recover)
  options(shiny.trace=TRUE)
  options(shiny.fullstacktrace=TRUE)
  options(shiny.testmode=TRUE)
} else {
  options(warn = 1, shiny.error = NULL)
  options(shiny.trace=FALSE)
  options(shiny.fullstacktrace=FALSE)
  options(shiny.testmode=FALSE)
}

#####################################################################################################
## file paths
getFile <- function(files){
  isPackage <- "MetFamily" %in% rownames(installed.packages())
  #isPackage <- FALSE
  
  resultFiles <- vector(mode = "character", length = length(files))
  for(idx in seq_along(files)){
    switch(files[[idx]], 
           "logo_ipb_en.png"={                                file <- ifelse(test = isPackage, yes = system.file("www/logo_ipb_en.png",                                           package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "www/logo_ipb_en.png",                                          sep = "/"))},
           
           "MetFamily_Input_Specification.pdf"={              file <- ifelse(test = isPackage, yes = system.file("doc/MetFamily_Input_Specification.pdf",                         package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/doc/MetFamily_Input_Specification.pdf",                        sep = "/"))},
           "MetFamily_user_guide.pdf"={                       file <- ifelse(test = isPackage, yes = system.file("doc/MetFamily_user_guide.pdf",                                  package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/doc/MetFamily_user_guide.pdf",                                 sep = "/"))},
           "MetFamily_Showcase_protocol.pdf"={                file <- ifelse(test = isPackage, yes = system.file("doc/MetFamily_Showcase_protocol.pdf",                           package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/doc/MetFamily_Showcase_protocol.pdf",                          sep = "/"))},
           "Fragment_matrix_showcase.csv"={                   file <- ifelse(test = isPackage, yes = system.file("data/showcase/Fragment_matrix_showcase.csv",                    package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/showcase/Fragment_matrix_showcase.csv",                   sep = "/"))},
           "Metabolite_profile_showcase.txt"={                file <- ifelse(test = isPackage, yes = system.file("data/showcase/Metabolite_profile_showcase.txt",                 package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/showcase/Metabolite_profile_showcase.txt",                sep = "/"))},
           "MSMS_library_showcase.msp"={                      file <- ifelse(test = isPackage, yes = system.file("data/showcase/MSMS_library_showcase.msp",                       package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/showcase/MSMS_library_showcase.msp",                      sep = "/"))},
           "Project_file_showcase_annotated.csv.gz"={         file <- ifelse(test = isPackage, yes = system.file("data/showcase/Project_file_showcase_annotated.csv.gz",          package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/showcase/Project_file_showcase_annotated.csv.gz",         sep = "/"))},
           "Project_file_showcase_annotated_reduced.csv.gz"={ file <- ifelse(test = isPackage, yes = system.file("data/showcase/Project_file_showcase_annotated_reduced.csv.gz",  package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/showcase/Project_file_showcase_annotated_reduced.csv.gz", sep = "/"))},
           "Project_file_showcase_reduced.csv.gz"={           file <- ifelse(test = isPackage, yes = system.file("data/showcase/Project_file_showcase_reduced.csv.gz",            package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/showcase/Project_file_showcase_reduced.csv.gz",           sep = "/"))},
           
           "Classifiers"={                                    file <- ifelse(test = isPackage, yes = system.file("data/classifier/",                                              package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "inst/data/classifier/",                                             sep = "/"))},
           
           "Analysis.R"={                                     file <- ifelse(test = isPackage, yes = "", no = paste("R/Analysis.R",                sep = "/"))},
           "DataProcessing.R"={                               file <- ifelse(test = isPackage, yes = "", no = paste("R/DataProcessing.R",          sep = "/"))},
           "FragmentMatrixFunctions.R"={                      file <- ifelse(test = isPackage, yes = "", no = paste("R/FragmentMatrixFunctions.R", sep = "/"))},
           "Plots.R"={                                        file <- ifelse(test = isPackage, yes = "", no = paste("R/Plots.R",                   sep = "/"))},
           "R_packages.R"={                                   file <- ifelse(test = isPackage, yes = "", no = paste("R/R_packages.R",              sep = "/"))},
           "StartApp.R"={                                     file <- ifelse(test = isPackage, yes = "", no = paste("R/StartApp.R",                sep = "/"))},
           "Annotation.R"={                                   file <- ifelse(test = isPackage, yes = "", no = paste("R/Annotation.R",              sep = "/"))},
           "Classifiers.R"={                                  file <- ifelse(test = isPackage, yes = "", no = paste("R/Classifiers.R",             sep = "/"))},
           "TreeAlgorithms.R"={                               file <- ifelse(test = isPackage, yes = "", no = paste("R/TreeAlgorithms.R",          sep = "/"))},
           {## unknown state
             stop(paste("Unknown file", file))
           }
    )
    resultFiles[[idx]] <- file
  }
  return(resultFiles)
}

#####################################################################################################
## source code
getSourceFileNames <- function(){
  return(c(
    "R_packages.R",
    "FragmentMatrixFunctions.R",
    "DataProcessing.R",
    "TreeAlgorithms.R",
    "Analysis.R",
    "Annotation.R",
    "Classifiers.R",
    "Plots.R"
  ))
}
getSourceFiles <- function(){
  return(getFile(getSourceFileNames()))
}
sourceTheCode <- function(){
  print("Testing whether files must be sourced")
  isPackage <- "MetFamily" %in% rownames(installed.packages())
  
  if(!isPackage){
    sourceFiles <- getSourceFiles()
    print(paste("Sourcing", length(sourceFiles), "files"))
    
    for(sourceFile in sourceFiles)
      source(sourceFile)
  } else {
    library("MetFamily")
  }
}
sourceTheCode()

#########################################################################################
#########################################################################################
## global variables
# none

#########################################################################################
#########################################################################################
## server-side logic of the Shiny app
shinyServer(
  func = function(input, output, session) {
    #########################################################################################
    #########################################################################################
    ## global variables per user
    
    ##############################################
    ## constants
    source("version.R")
      
    ## annotation constants
    artifactName   <- "Ignore"
    artifactColor  <- "red"
    selectionNone  <- "None"
    
    ## Download names
    ExportMatrixName <- NULL
    
    ## GUI constants
    runRightColumnWidthFull <- 12
    legendColumnWidthFull <- 2
    runRightColumnWidthPart <- 8
    legendColumnWidthPart <- 2
    
    annoLegendEntryHeight <- 20
    maximumNumberOfTableEntries <- 50
    
    ##############################################
    ## program state
    initialGuiUpdatePerformed <- FALSE
    state <- reactiveValues(
      ## side bar handling
      runRightColumnWidth = runRightColumnWidthPart, 
      legendColumnWidth = legendColumnWidthPart,
      showSideBar = TRUE, 
      ## HCA vs PCA plots handling
      analysisType = "HCA",
      anyPlotDrawn = FALSE,
      ## HCA / PCA / classifier handling
      showHCAplotPanel = FALSE, 
      showPCAplotPanel = FALSE, 
      showAnnotationplotPanel = FALSE, 
      plotHcaShown = FALSE,
      plotPcaShown = FALSE,
      plotAnnotationShown = FALSE,
      ## plot controls
      showPlotControls = FALSE
    )
    plotToShow  <- "Display HCA"
    plotsToShow <- "Display HCA"
    showSideBar <- TRUE
    
    #########################################################################################
    #########################################################################################
    ## functions
    
    #########################################################################################
    ## source all server stuff
    resetWorkspaceFunctions <- list()
    suspendOnExitFunctions <- list()
    
    source(file = "app_files/server_functionsFilters.R", local = TRUE)$value
    source(file = "app_files/server_functionsSelections.R", local = TRUE)$value
    source(file = "app_files/server_functionsTableGui.R", local = TRUE)$value
    source(file = "app_files/server_functionsDownloads.R", local = TRUE)$value
    source(file = "app_files/server_functionsSerialization.R", local = TRUE)$value
    source(file = "app_files/server_guiDialogs.R", local = TRUE)$value
    source(file = "app_files/server_guiPlots.R", local = TRUE)$value
    source(file = "app_files/server_guiAnnotation.R", local = TRUE)$value
    source(file = "app_files/server_guiTabInput.R", local = TRUE)$value
    source(file = "app_files/server_guiTabAnnotation.R", local = TRUE)$value
    source(file = "app_files/server_guiTabClassifier.R", local = TRUE)$value
    source(file = "app_files/server_guiTabSampleFilter.R", local = TRUE)$value
    source(file = "app_files/server_guiTabMsmsFilter.R", local = TRUE)$value
    source(file = "app_files/server_guiTabPca.R", local = TRUE)$value
    source(file = "app_files/server_guiTabHca.R", local = TRUE)$value
    source(file = "app_files/server_guiTabSearch.R", local = TRUE)$value
    source(file = "app_files/server_guiTabExport.R", local = TRUE)$value
    source(file = "app_files/server_guiPlotControls.R", local = TRUE)$value
    source(file = "app_files/server_guiMs2plot.R", local = TRUE)$value
    ## ui generation
    #source(file = "app_files/ui_rightColumn.R", local = TRUE)
    
    ## Parse the input file
    resetWorkspace <- function(){
      print(paste("resetWorkspace"))
      
      for(resetWorkspaceFunction in resetWorkspaceFunctions)
        resetWorkspaceFunction()
      
      #########################################################################################
      ## panels
      state$showHCAplotPanel <<- FALSE
      state$showPCAplotPanel <<- FALSE
      state$showAnnotationplotPanel <<- FALSE
      
      state$plotHcaShown <<- FALSE
      state$plotPcaShown <<- FALSE
      state$plotAnnotationShown <<- FALSE
      
      state$analysisType <<- "HCA"
      state$anyPlotDrawn <<- FALSE
      ## plot controls
      showPlotControls <<- FALSE
    }
    
    #########################################################################################
    #########################################################################################
    ## observer
    
    ## controls
    obsTabs <- observeEvent(input$runTabs, {
      tabId <- input$runTabs
      print(paste("observe tabs", tabId))
      if(tabId == "HCA"){
        state$analysisType <<- "HCA"
        if(state$showHCAplotPanel){
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
          state$plotHcaShown <<- TRUE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      if(tabId == "PCA"){
        state$analysisType <<- "PCA"
        if(state$showPCAplotPanel){
          state$plotHcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
          state$plotPcaShown <<- TRUE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      if(tabId == "Classifiers"){
        state$analysisType <<- "Annotation"
        if(state$showAnnotationplotPanel){
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- TRUE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      
      #########################################################
      ## initial gui update
      if(tabId == "Input" & !initialGuiUpdatePerformed){
        print(paste("update GUI initially", tabId))
        
        #shinyjs::disable("importMs1Ms2Data")
        #shinyjs::disable("importMs2Data")
        #shinyjs::disable("loadProjectData")
        
        ## annotation classifier selection
        #session$sendCustomMessage("disableButton", "doAnnotation")
        
        filePath <- getFile("Classifiers")
        resultObj <- getAvailableClassifiers(filePath)
        availableClassifiersDf           <<- resultObj$availableClassifiersDf
        availableClassifiersDfProperties <<- resultObj$availableClassifiersDfProperties
        
        output$classifierCount <- renderText({
          print(paste("init output$classifierCount"))
          paste("Available classifiers:", nrow(availableClassifiersDf))
        })
        
        if(nrow(availableClassifiersDf) > 0){
          output$classifierSelectionTable <- DT::renderDataTable(
            expr = availableClassifiersDf,
            server = FALSE, escape = FALSE, rownames = FALSE,
            selection = list(mode = "single"),#, selected = 1), 
            options = list(
              #scrollY = "600px",
              scrollX = "40vh",
              scrollY = "60vh",
              preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
              drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
              iDisplayLength=nrow(availableClassifiersDf),       # initial number of records
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
        
        if(nrow(availableClassifiersDf) > 0)
          shinyjs::toggleState("doAnnotation", TRUE)
        else
          shinyjs::toggleState("doAnnotation", FALSE)
        
        
        initialGuiUpdatePerformed <<- TRUE
      }
      
      if(tabId == "HCA" | tabId == "PCA")
        updateSelectedSelection()
    })
    obsChangePlot <- observeEvent(input$changePlot, {
      plot <- input$changePlot
      print(paste("Observe changePlot", plot))
      if(plot == "Display HCA"){
        analysisType <- "HCA"
        if(state$showHCAplotPanel){
          state$plotHcaShown <<- TRUE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      if(plot == "Display PCA") {
        analysisType <- "PCA"
        if(state$showPCAplotPanel){
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- TRUE
          state$plotAnnotationShown <<- FALSE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      if(plot == "Display Annotation") {
        analysisType <- "Annotation"
        if(state$showAnnotationplotPanel){
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- TRUE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      state$analysisType <<- analysisType
      
      updateSelectedSelection()
    })
    obsShowSideBar <- observeEvent(input$showSideBar, {
      showSideBar <- input$showSideBar
      print(paste("Observe showSideBar", showSideBar))
      
      state$showSideBar <<- showSideBar
      showSideBar <<- showSideBar
      
      if(showSideBar){
        state$runRightColumnWidth <<- runRightColumnWidthPart
        state$legendColumnWidth   <<- legendColumnWidthPart
      }
      else{
        state$runRightColumnWidth <<- runRightColumnWidthFull
        state$legendColumnWidth   <<- legendColumnWidthFull
      }
      
      ## restore gui state
      ## TODO check if necessary
      updateCheckboxInput(session = session, inputId = "showPlotControls",      value    = state$showPlotControls)
      updateCheckboxInput(session = session, inputId = "showClusterLabels",     value    = state_tabHca$showClusterLabels)
      updateRadioButtons( session = session, inputId = "heatmapContent",        selected = state_tabHca$heatmapContent)
      updateRadioButtons( session = session, inputId = "heatmapOrdering",       selected = state_tabHca$heatmapOrdering)
      updateRadioButtons( session = session, inputId = "hcaPrecursorLabels",    selected = state_tabHca$hcaPrecursorLabels)
      updateCheckboxInput(session = session, inputId = "showScoresLabels",      value    = state_tabPca$showScoresLabels)
      updateRadioButtons( session = session, inputId = "loadingsLabels",        selected = state_tabPca$loadingsLabels)
      updateCheckboxGroupInput(session = session, inputId = "showLoadingsFeatures", selected = c(
        ifelse(test = state_tabPca$showLoadingsFeaturesAnnotated,   yes = "Annotated",     no = NULL),
        ifelse(test = state_tabPca$showLoadingsFeaturesUnannotated, yes = "Not Annotated", no = NULL),
        ifelse(test = state_tabPca$showLoadingsFeaturesSelected,    yes = "Selected",      no = NULL),
        ifelse(test = state_tabPca$showLoadingsFeaturesUnselected,  yes = "Not Selected",  no = NULL)
      ))
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesAnnotated",    selected = state_tabPca$showLoadingsFeaturesAnnotated)
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesUnannotated",  selected = state_tabPca$showLoadingsFeaturesUnannotated)
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesSelected",     selected = state_tabPca$showLoadingsFeaturesSelected)
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesUnselected",   selected = state_tabPca$showLoadingsFeaturesUnselected)
      updateCheckboxInput(session = session, inputId = "showLoadingsAbundance", value    = state_tabPca$showLoadingsAbundance)
    })
    
    ## display of tabs
    observe({
      #toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Filter']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='MS/MS filter']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Sample filter']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='PCA']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='HCA']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Search']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Classifiers']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Annotations']")
      toggle(condition = !is.null(state_tabInput$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Project']")
    })
    
    
    suspendOnExitFunctions <- c(suspendOnExitFunctions, function(){
      print("Suspending server observers")
      obsTabs$suspend()
      obsChangePlot$suspend()
      obsShowSideBar$suspend()
    })
    
    ## suspend observer
    session$onSessionEnded(function() {
      print("Suspending observers")
      for(suspendOnExitFunction in suspendOnExitFunctions)
        suspendOnExitFunction()
    })
    
    #########################################################################################
    #########################################################################################
    ## direct output rendering
    output$information <- renderText({
      print(paste("init output$information"))
      ""
    })
    ## about page
    output$rInfo <- renderText({
      print(paste("init rInfo"))
      #paste(R.Version()$version.string, "\nwd: ", getwd(), sep = "")
      paste(
        R.Version()$version.string, 
        "\nMetFamily build: ", metFamilyBuilt, "-", system(command = "hostname", intern = TRUE),
        "\nMetFamily package: ", packageVersion,
        #"\n", getwd(),
        #"\n", paste(list.files(path = getwd(), all.files = F, full.names = F, recursive = T, include.dirs = F), collapse = "\n"),
        sep = ""
      )
      #R.Version()$version.string
    })
    output$ipbImage <- renderImage({
      file <- getFile("logo_ipb_en.png")
      
      list(src = file,
           alt = "IPB Halle"
      )
    }, deleteFile = FALSE)
    
    #########################################################################################
    #########################################################################################
    ## reactive output values
    output$showGUI <- reactive({
      print("update output$showGUI")
      output$information <- renderText({
        print(paste("init information", sep = ""))
        paste("Please perform ploting.", sep = "")
      })
      return(!is.null(state_tabInput$importedOrLoadedFile_s_))
    })
    output$analysisType <- reactive({
      print(paste("reactive update analysisType", state$analysisType))
      if(!is.null(state$analysisType)){
        if(state$analysisType == "HCA")
          plotToShow <<- "Display HCA"
        if(state$analysisType == "PCA")
          plotToShow <<- "Display PCA"
        if(state$analysisType == "Annotation")
          plotToShow <<- "Display Annotation"
      }
      return(state$analysisType)
    })
    output$showSideBar <- reactive({
      print(paste("reactive update showSideBar", state$showSideBar))
      return(state$showSideBar)
    })
    output$showHCAplotPanel <- reactive({
      print(paste("reactive update showHCAplotPanel", state$showHCAplotPanel))
      updateChangePlotRadioButton()
      return(state$showHCAplotPanel)
    })
    output$showPCAplotPanel <- reactive({
      print(paste("reactive update showPCAplotPanel", state$showPCAplotPanel))
      updateChangePlotRadioButton()
      return(state$showPCAplotPanel)
    })
    output$showAnnotationplotPanel <- reactive({
      print(paste("reactive update showAnnotationplotPanel", state$showAnnotationplotPanel))
      updateChangePlotRadioButton()
      return(state$showAnnotationplotPanel)
    })
    output$plotHcaShown <- reactive({
      print(paste("reactive update plotHcaShown", state$plotHcaShown))
      return(state$plotHcaShown)
    })
    output$plotPcaShown <- reactive({
      print(paste("reactive update plotPcaShown", state$plotPcaShown))
      return(state$plotPcaShown)
    })
    output$plotAnnotationShown <- reactive({
      print(paste("reactive update plotAnnotationShown", state$plotAnnotationShown))
      return(state$plotAnnotationShown)
    })
    #output$putativeAnnotationsTableFromAnalysisRowSelected <- reactive({
    #  print(paste("reactive update putativeAnnotationsTableFromAnalysisRowSelected", state$putativeAnnotationsTableFromAnalysisRowSelected))
    #  return(state$putativeAnnotationsTableFromAnalysisRowSelected)
    #})
    
    updateChangePlotRadioButton <- function(){
      if((sum(c(state$showHCAplotPanel, state$showPCAplotPanel, state$showAnnotationplotPanel)) > 1) & !is.null(state$analysisType)){
        if(state$analysisType == "HCA")
          selectedItem <- "Display HCA"
        if(state$analysisType == "PCA")
          selectedItem <- "Display PCA"
        if(state$analysisType == "Annotation")
          selectedItem <- "Display Annotation"
        
        shownPlots <- NULL
        if(state$showHCAplotPanel)
          shownPlots <- c(shownPlots, "Display HCA")
        if(state$showPCAplotPanel)
          shownPlots <- c(shownPlots, "Display PCA")
        if(state$showAnnotationplotPanel)
          shownPlots <- c(shownPlots, "Display Annotation")
        
        plotsToShow <<- shownPlots
        updateRadioButtons(session = session, inputId = "changePlot", selected = selectedItem, choices = shownPlots, inline = TRUE)
      }
    }
    
    #########################################################################################
    #########################################################################################
    ## properties
    options(shiny.maxRequestSize=1024*1024^2) ## 500 mb file size
    outputOptions(output, 'showGUI',                 suspendWhenHidden=FALSE)
    outputOptions(output, 'showSideBar',             suspendWhenHidden=FALSE)
    outputOptions(output, 'showHCAplotPanel',        suspendWhenHidden=FALSE)
    outputOptions(output, 'showPCAplotPanel',        suspendWhenHidden=FALSE)
    outputOptions(output, 'showAnnotationplotPanel', suspendWhenHidden=FALSE)
    outputOptions(output, 'analysisType',            suspendWhenHidden=FALSE)
    outputOptions(output, 'plotHcaShown',            suspendWhenHidden=FALSE)
    outputOptions(output, 'plotPcaShown',            suspendWhenHidden=FALSE)
    outputOptions(output, 'plotAnnotationShown',     suspendWhenHidden=FALSE)
    #outputOptions(output, 'putativeAnnotationsTableFromAnalysisRowSelected',  suspendWhenHidden=FALSE)
  }## function(input, output, session)
)## shinyServer
