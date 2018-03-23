
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
## redirect console output to file


if(FALSE & !isDevelopment){
  logPath <- "var/log/shiny-server/"
  timeStamp <- gsub(" ", "_", gsub(":", ".", Sys.time()))
  logFileName <- paste(timeStamp, "MetFamily", "Console", "Output", sep = "_")
  logFilePath <- paste(logPath, logFileName, sep = "/")
  
  writeToFile <- TRUE
  if(!file.exists(logPath))
    if(!dir.create(path = logPath, recursive = TRUE))
      writeToFile <- FALSE
  if(writeToFile)
    sink(file = logFilePath, type=c("output", "message"))
}

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
  #isPackage <- "MetFamily" %in% rownames(installed.packages())
  isPackage <- FALSE
  
  resultFiles <- vector(mode = "character", length = length(files))
  for(idx in seq_along(files)){
    switch(files[[idx]], 
           "logo_ipb_en.png"={                                file <- ifelse(test = isPackage, yes = system.file("www/logo_ipb_en.png",                                           package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/www/logo_ipb_en.png",                                          sep = "/"))},
           "MetFamily_Input_Specification.pdf"={              file <- ifelse(test = isPackage, yes = system.file("doc/MetFamily_Input_Specification.pdf",                         package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/doc/MetFamily_Input_Specification.pdf",                        sep = "/"))},
           "MetFamily_user_guide.pdf"={                       file <- ifelse(test = isPackage, yes = system.file("doc/MetFamily_user_guide.pdf",                                  package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/doc/MetFamily_user_guide.pdf",                                 sep = "/"))},
           "MetFamily_Showcase_protocol.pdf"={                file <- ifelse(test = isPackage, yes = system.file("doc/MetFamily_Showcase_protocol.pdf",                           package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/doc/MetFamily_Showcase_protocol.pdf",                          sep = "/"))},
           "Fragment_matrix_showcase.csv"={                   file <- ifelse(test = isPackage, yes = system.file("data/showcase/Fragment_matrix_showcase.csv",                    package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/showcase/Fragment_matrix_showcase.csv",                   sep = "/"))},
           "Metabolite_profile_showcase.txt"={                file <- ifelse(test = isPackage, yes = system.file("data/showcase/Metabolite_profile_showcase.txt",                 package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/showcase/Metabolite_profile_showcase.txt",                sep = "/"))},
           "MSMS_library_showcase.msp"={                      file <- ifelse(test = isPackage, yes = system.file("data/showcase/MSMS_library_showcase.msp",                       package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/showcase/MSMS_library_showcase.msp",                      sep = "/"))},
           "Project_file_showcase_annotated.csv.gz"={         file <- ifelse(test = isPackage, yes = system.file("data/showcase/Project_file_showcase_annotated.csv.gz",          package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/showcase/Project_file_showcase_annotated.csv.gz",         sep = "/"))},
           "Project_file_showcase_annotated_reduced.csv.gz"={ file <- ifelse(test = isPackage, yes = system.file("data/showcase/Project_file_showcase_annotated_reduced.csv.gz",  package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/showcase/Project_file_showcase_annotated_reduced.csv.gz", sep = "/"))},
           "Project_file_showcase_reduced.csv.gz"={           file <- ifelse(test = isPackage, yes = system.file("data/showcase/Project_file_showcase_reduced.csv.gz",            package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/showcase/Project_file_showcase_reduced.csv.gz",           sep = "/"))},
           
           #"Classifiers"={                                    file <- ifelse(test = isPackage, yes = system.file("data/classifiers/",                                             package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/classifiers/",                                            sep = "/"))},
           "Classifiers"={                                    file <- ifelse(test = isPackage, yes = system.file("data/classifier/",                                              package = "MetFamily", lib.loc=.libPaths()), no = paste(sourceFolder, "../inst/data/classifier/",                                             sep = "/"))},
           #resultFolderForClassifiers <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Classifiers"
           #resultFolderForClassifiers <- "/mnt/Treutler/Data/SubstanceClasses/Classifier_ROC_Analysis/Classifiers"
           
           "Analysis.R"={                                     file <- ifelse(test = isPackage, yes = "", no = paste("Analysis.R",                sep = "/"))},
           "DataProcessing.R"={                               file <- ifelse(test = isPackage, yes = "", no = paste("DataProcessing.R",          sep = "/"))},
           "FragmentMatrixFunctions.R"={                      file <- ifelse(test = isPackage, yes = "", no = paste("FragmentMatrixFunctions.R", sep = "/"))},
           "Plots.R"={                                        file <- ifelse(test = isPackage, yes = "", no = paste("Plots.R",                   sep = "/"))},
           "R_packages.R"={                                   file <- ifelse(test = isPackage, yes = "", no = paste("R_packages.R",              sep = "/"))},
           "StartApp.R"={                                     file <- ifelse(test = isPackage, yes = "", no = paste("StartApp.R",                sep = "/"))},
           "Annotation.R"={                                   file <- ifelse(test = isPackage, yes = "", no = paste("Annotation.R",              sep = "/"))},
           "Classifiers.R"={                                  file <- ifelse(test = isPackage, yes = "", no = paste("Classifiers.R",             sep = "/"))},
           "TreeAlgorithms.R"={                               file <- ifelse(test = isPackage, yes = "", no = paste("TreeAlgorithms.R",          sep = "/"))},
           {## unknown state
             stop(paste("Unknown file", file))
           }
    )
    #print(paste(isPackage,files[[idx]],file))
    if(!file.exists(file))
      file <- gsub(x = file, pattern = "\\.\\./inst/", replacement = "")
    
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
    print(paste("Sourcing", length(sourceFiles), "files"))
    sourceFiles <- getSourceFiles()
    
    for(sourceFile in sourceFiles)
      source(sourceFile)
  } else {
    library("MetFamily")
  }
}

if(isDevelopment){
  sourceFiles <- paste(sourceFolder, "../R", getSourceFileNames(), sep = "/")
  for(sourceFile in sourceFiles)
    source(sourceFile)
} else {
  for(sourceFile in getSourceFileNames())
    source(sourceFile)
  #sourceTheCode()
}

if(!isDevelopment)  setwd("/var/log/shiny-server")

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
    
    ## MetFamily properties
    toolName    <- "MetFamily"
    toolVersion <- "1.0"
    metFamilyBuilt <- "1.1.0"
    
    ## data import
    proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- 0.9
    mzDeviationAbsolute_mapping <- 0.01
    #minimumNumberOfMS2PeaksPerGroup <- 1
    ## annotation
    artifact <- "Ignore"
    artifactColor <- "red"
    ## HCA
    minimumNumberOfPrecursorsForHca <- 6
    maximumNumberOfPrecursorsForHca <- 1000000
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
    ## GUI
    runRightColumnWidthFull <- 12
    legendColumnWidthFull <- 2
    runRightColumnWidthPart <- 8
    legendColumnWidthPart <- 2
    minimumheatmapHeightPerRow <- 11
    maximumheatmapHeightPerRow <- 25
    minimumheatmapHeightRowCount <- 10
    maximumheatmapHeightRowCount <- 3
    annoLegendEntryHeight <- 20
    scoresGroupsLegendEntryHeight <- 20
    maximumNumberOfTableEntries <- 50
    
    ##############################################
    ## data
    dataList <- NULL
    currentDistanceMatrixObj <- NULL
    clusterDataList <- NULL
    pcaDataList <- NULL
    ## filter
    filterGlobal <- NULL
    filterHca <- NULL
    filterPca <- NULL
    filterSearch <- NULL
    
    ##############################################
    ## program state
    initialGuiUpdatePerformed <- FALSE
    state <- reactiveValues(
      importedOrLoadedFile_s_ = NULL, 
      globalMS2filterValid = FALSE, 
      hcafilterValid = FALSE, 
      pcafilterValid = FALSE, 
      searchfilterValid = FALSE, 
      filterSearchActive = FALSE, 
      runRightColumnWidth = runRightColumnWidthPart, 
      legendColumnWidth = legendColumnWidthPart,
      dendrogramHeatmapHeight = 1,#heatmapHeightPerRow * 3,
      heatmapHeight = 1,#heatmapHeightPerRow * 3,## plotly: remove
      showSideBar = TRUE, 
      analysisType = "HCA",
      anyPlotDrawn = FALSE,
      showHCAplotPanel = FALSE, 
      showPCAplotPanel = FALSE, 
      showAnnotationplotPanel = FALSE, 
      plotHcaShown = FALSE,
      plotPcaShown = FALSE,
      plotAnnotationShown = FALSE,
      precursorSetSelected = FALSE,
      selectedSelection = NULL,
      showPlotControls = FALSE,
      showClusterLabels = TRUE,
      heatmapContent = "Log-fold-change",
      heatmapOrdering = "Specified order",
      hcaPrecursorLabels = TRUE,
      showScoresLabels = TRUE,
      loadingsLabels = "None",
      showLoadingsFeaturesAnnotated = TRUE, 
      showLoadingsFeaturesUnannotated = TRUE, 
      showLoadingsFeaturesSelected = TRUE, 
      showLoadingsFeaturesUnselected = TRUE,
      showLoadingsAbundance = FALSE,
      annotationLegendHeightHca = -1,
      annotationLegendHeightPca = -1,
      ## plot annotations
      annotationsHca = NULL,
      annotationsPca = NULL,
      scoresGroupsLegendHeight = -1
    )
    scoresGroups <- NULL
    changeSelectionCurrentSelection <- selectionAnalysisName
    precursorSelectionTabCurrentTab <- precursorSelectionTabSelection
    plotToShow  <- "Display HCA"
    plotsToShow <- "Display HCA"
    showSideBar <- TRUE
    
    sampleOrder_tmp <- NULL
    sampleExclusion_tmp <- NULL
    curveNumberToCurveName <- NULL
    
    ## classifier annotation
    availableClassifiersDf <- NULL
    availableClassifiersDfProperties <- NULL
    classToSpectra_class <- NULL
    mappingSpectraToClassDf <- NULL
    properties_class <- NULL
    selectedClassPrecursorIndeces <- NULL
    selectedClassRowIdx <- -1
    selectedClassFeatureRowIdx <- -1
    
    ##############################################
    ## buttons
    updateProjectDescriptionButtonValue <- 0
    applyGlobalMS2filtersButtonValue <- 0
    clearGlobalMS2filtersButtonValue <- 0
    applySearchButtonValue <- 0
    clearSearchButtonValue <- 0
    applyHcaFiltersButtonValue <- 0
    clearHcaFiltersButtonValue <- 0
    applyPcaFiltersButtonValue <- 0
    clearPcaFiltersButtonValue <- 0
    drawHCAButtonValue <- 0
    drawPCAButtonValue <- 0
    downloadMatrixButtonValue <- 0
    removePresentAnnotationValue <- 0
    setPresentAnnotationPrimaryValue <- 0
    submitNewAnnotationValue <- 0
    submitPreviousAnnotationValue <- 0
    updateArtifactsFromCheckboxesButtonValue <- 0
    updateSamplesButtonValue <- 0
    clearSelectionButtonValue <- 0
    importMs1Ms2DataButtonValue <- 0
    importMs2DataButtonValue <- 0
    loadProjectDataButtonValue <- 0
    loadExampleDataButtonValue <- 0
    confirmAnnotationButtonValue <- 0
    
    ##############################################
    ## MS2 peaks
    fragmentsX <- NULL
    fragmentsY <- NULL
    fragmentsColor <- NULL
    fragmentsDiscriminativity <- NULL
    fragmentsXhovered <- NULL
    fragmentsYhovered <- NULL
    fragmentsColorHovered <- NULL
    fragmentsDiscriminativityHovered <- NULL
    
    ##############################################
    ## plot ranges
    dendrogramPlotRangeY <- NULL
    dendrogramPlotRange  <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
    fragmentPlotRange    <- reactiveValues(xMin = NULL, xMax = NULL)
    ms2PlotRange         <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
    pcaScoresPlotRange   <- reactiveValues(
      xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL, 
      yMin = NULL, yMax = NULL, yInterval = NULL, yIntervalSize = NULL
    )
    pcaLoadingsPlotRange <- reactiveValues(
      xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL, 
      yMin = NULL, yMax = NULL, yInterval = NULL, yIntervalSize = NULL
    )
    specVsClassRange     <- reactiveValues(xMin = NULL, xMax = NULL, xInterval = NULL, xIntervalSize = NULL)
    
    ##############################################
    ## selections
    
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
    sampleTableInputFieldIdCounter <- 0
    ms1FeatureVsClassTableCounter <- 0
    sampleTable <- NULL
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
    
    #########################################################################################
    #########################################################################################
    ## functions
    
    ## POI selection
    getSelectedPOI_X <- function(mouseX, poiCoordinatesX, plotWidth, plotRangeX){
      if(any(is.na(c(poiCoordinatesX))))
        return(NULL)
      
      factorX <- plotWidth  / plotRangeX
      
      mouseX <- mouseX * factorX
      poiCoordinatesX <- poiCoordinatesX * factorX
      
      distances <- abs(poiCoordinatesX - mouseX)
      distanceThreshold <- factorX * plotRangeX / 35
      
      minimumIndex <- which.min(distances)
      minimumDistance <- distances[[minimumIndex]]
      
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    getSelectedPOI_XY <- function(mouseX, mouseY, poiCoordinatesX, poiCoordinatesY, plotWidth, plotHeight, plotRangeX, plotRangeY){
      ## see also nearPoints(...): http://shiny.rstudio.com/gallery/plot-interaction-selecting-points.html
      if(any(is.na(c(poiCoordinatesX, poiCoordinatesY))))
        return(NULL)
      
      factorX <- plotWidth  / plotRangeX
      factorY <- plotHeight / plotRangeY
      
      mouseX <- mouseX * factorX
      mouseY <- mouseY * factorY
      poiCoordinatesX <- poiCoordinatesX * factorX
      poiCoordinatesY <- poiCoordinatesY * factorY
      
      distancesX <- poiCoordinatesX - mouseX
      distancesY <- poiCoordinatesY - mouseY
      distances <- sqrt(distancesX * distancesX + distancesY * distancesY)
      distanceThreshold <- factorX * plotRangeX / 35
      
      minimumIndex <- which.min(distances)
      if(length(minimumIndex)==0)
        return(NULL)
      if(is.na(minimumIndex))
        return(NULL)
      minimumDistance <- distances[[minimumIndex]]
      if(minimumDistance > distanceThreshold){
        return(NULL)
      } else {
        return(minimumIndex)
      }
    }
    ## Parse the input file
    resetWorkspace <- function(){
      print(paste("resetWorkspace"))
      
      #########################################################################################
      ## reset
      
      ## reset plots
      doClearPlots()
      
      ## reset variables
      clusterDataList <<- NULL
      pcaDataList <<- NULL
      
      ## selection
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
      
      
      ## fragments
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      fragmentsColor <<- NULL
      fragmentsDiscriminativity <<- NULL
      fragmentsXhovered <<- NULL
      fragmentsYhovered <<- NULL
      fragmentsColorHovered <<- NULL
      fragmentsDiscriminativityHovered <<- NULL
      
      ## reset state
      state$importedOrLoadedFile_s_ <<- NULL
      state$analysisType <<- "HCA"
      state$showHCAplotPanel <<- FALSE
      state$showPCAplotPanel <<- FALSE
      state$showAnnotationplotPanel <<- FALSE
      state$plotHcaShown <<- FALSE
      state$plotPcaShown <<- FALSE
      state$plotAnnotationShown <<- FALSE
      state$precursorSetSelected <<- FALSE
      state$anyPlotDrawn <<- FALSE
      selectedPrecursorSet <<- NULL
      selectedTable <<- NULL
      state$annotationsHca <<- NULL
      state$annotationsPca <<- NULL
      scoresGroups <<- NULL
      
      ## reset plot range
      dendrogramPlotRangeY <<- NULL
      dendrogramPlotRange$xMin <<- NULL
      dendrogramPlotRange$xMax <<- NULL
      dendrogramPlotRange$xInterval <<- NULL
      dendrogramPlotRange$xIntervalSize <<- NULL
      pcaScoresPlotRange$xMin <<- NULL
      pcaScoresPlotRange$xMax <<- NULL
      pcaScoresPlotRange$xInterval <<- NULL
      pcaScoresPlotRange$xIntervalSize <<- NULL
      pcaScoresPlotRange$yMin <<- NULL
      pcaScoresPlotRange$yMax <<- NULL
      pcaScoresPlotRange$yInterval <<- NULL
      pcaScoresPlotRange$yIntervalSize <<- NULL
      pcaLoadingsPlotRange$xMin <<- NULL
      pcaLoadingsPlotRange$xMax <<- NULL
      pcaLoadingsPlotRange$xInterval <<- NULL
      pcaLoadingsPlotRange$xIntervalSize <<- NULL
      pcaLoadingsPlotRange$yMin <<- NULL
      pcaLoadingsPlotRange$yMax <<- NULL
      pcaLoadingsPlotRange$yInterval <<- NULL
      pcaLoadingsPlotRange$yIntervalSize <<- NULL
      specVsClassRange$xMin <<- NULL
      specVsClassRange$xMax <<- NULL
      
      #########################################################################################
      ## update fragment plot
      
      min <- min(dataList$masses)
      max <- max(dataList$masses)
      
      fragmentPlotRange$xMin <<- min
      fragmentPlotRange$xMax <<- max
      fragmentPlotRange$xInterval <<- c(min, max)
      fragmentPlotRange$xIntervalSize <<- max - min
      
      output$fragmentPlot <- renderPlot({
        print(paste("### MS2 ### all init"))
        plotFragments(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
      #}, bg = "transparent")
      })
      
      #########################################################################################
      ## update filter input values
      
      ## groups
      switch(as.character(length(dataList$groups)), 
        "0"={
          stop("No groups available")
        },
        "1"={
          selectedOne <- dataList$groups[[1]]
          selectedTwo <- dataList$groups[[1]]
        },
        {
          selectedOne <- dataList$groups[[1]]
          selectedTwo <- dataList$groups[[2]]
        }
      )
      
      sampleNames <- dataList$groupSampleDataFrame[, "Sample"]
      
      updateRadioButtons(session = session, inputId = "hcaFilterGroupOne", choices = dataList$groups, selected = selectedOne)
      updateRadioButtons(session = session, inputId = "hcaFilterGroupTwo", choices = dataList$groups, selected = selectedTwo)
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = dataList$groups)
      updateCheckboxGroupInput(session = session, inputId = "pcaSamples",  choices = sampleNames,     selected = sampleNames)
      
      ## input fields: HCA
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "globalFilter_ms2_ppm", value = "20")
      
      ## input fields: HCA
      updateTextInput(session = session, inputId = "hcaFilter_average", value = "0")
      updateTextInput(session = session, inputId = "hcaFilter_lfc", value = "0")
      updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
      ## input fields: PCA
      updateTextInput(session = session, inputId = "pcaFilter_average", value = "0")
      updateTextInput(session = session, inputId = "pcaFilter_lfc", value = "0")
      updateCheckboxInput(session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", value = FALSE)
      ## input fields: search MS1
      updateTextInput(session = session, inputId = "searchMS1mass", value = "")
      updateTextInput(session = session, inputId = "searchMS1massPpm", value = 20)
      ## input fields: search MS2
      updateTextInput(session = session, inputId = "search_ms2_masses1", value = "")
      updateTextInput(session = session, inputId = "search_ms2_masses2", value = "")
      updateTextInput(session = session, inputId = "search_ms2_masses3", value = "")
      updateTextInput(session = session, inputId = "searchMS2massPpm", value = 20)
      updateCheckboxInput(session = session, inputId = "searchIncludeIgnoredPrecursors", value = FALSE)
      #updateColourInput(session = session, inputId = "newAnnotationColor", allowedCols = colorPalette())
      ## anno
      updateTextInput(session = session, inputId = "newAnnotationValue", value = "")
      updateTextInput(session = session, inputId = "newAnnotationValue2", value = "")
      updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[init]"))
      updateSelectInput(session = session, inputId = "previousAnnotationValue", choices = c("Artifact"))
      
      #########################################################################################
      ## update filter
      sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filter <- doPerformFiltering(dataList$groups, sampleSet, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      if(length(dataList$groups) == 1)
        filter2 <- doPerformFiltering(c(dataList$groups[[1]], dataList$groups[[1]]), NULL, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      else
        filter2 <- filter
      
      filterGlobal <<- filter
      filterHca    <<- filter2
      filterPca    <<- filter
      state$filterSearchActive <<- FALSE
      state$searchfilterValid <<- TRUE
      filterSearch    <<- NULL
      
      updateGlobalMS2filterInformation()
      updateHcaFilterInformation()
      updatePcaFilterInformation()
      updateSearchInformation()
      
      state$globalMS2filterValid <<- TRUE
      state$hcaFilterValid <<- TRUE
      state$pcaFilterValid <<- TRUE
      
      checkHcaFilterValidity(filter2$numberOfPrecursorsFiltered)
      checkPcaFilterValidity(filter$numberOfPrecursorsFiltered)
      
      state$plotHcaShown <<- FALSE
      state$plotPcaShown <<- FALSE
      state$plotAnnotationShown <<- FALSE
      
      ## project infos
      updateTextInput    (session = session, inputId = "projectName2",           value = dataList$importParameterSet$projectName)
      shinyjs::toggleState("projectName2", FALSE)
      updateTextInput    (session = session, inputId = "projectDescription2",    value = dataList$importParameterSet$projectDescription)
      #shinyjs::toggleState("projectName2", FALSE)
      
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
      
      ## project tab
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
      
      ## MS2 plot range
      resetMS2PlotRange()
      
      ## sample table
      sampleOrder_tmp     <<- dataList$groupSampleDataFrame[, "Order"]
      sampleExclusion_tmp <<- dataList$groupSampleDataFrame[, "Exclude"]
      curveNumberToCurveName <<- NULL
      
      ## classifier annotation
      classToSpectra_class <<- NULL
      mappingSpectraToClassDf <<- NULL
      properties_class <<- NULL
      selectedClassPrecursorIndeces <<- NULL
      selectedClassRowIdx <<- -1
      selectedClassFeatureRowIdx <<- -1
      
      sampleTable <<- createSampleTable()
      setSampleTable()
      #updateAnnotationOverview()
    }
    ## filter info
    updateGlobalMS2filterInformation <- function(){
      if(is.null(filterGlobal)){
        output$globalMS2filteredPrecursors <- renderText({
          print(paste("update output$globalMS2filteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
      } else {
        output$globalMS2filteredPrecursors <- renderText({
          print(paste("update output$globalMS2filteredPrecursors", sep = ""))
          paste("Number of filtered MS1 features: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, sep = "")
        })
      }
    }
    updateHcaFilterInformation <- function(){
      if(is.null(filterHca)){
        ## errors
        output$hcaFilteredPrecursors <- renderText({
          print(paste("update output$hcaFilteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values.", sep = "")
        })
      } else {
        ## no errors
        globalMs2Filter <- ifelse(filterGlobal$numberOfPrecursorsFiltered != dataList$numberOfPrecursors, paste("\n(Global MS/MS filter: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, " precursors)", sep = ""), "")
        
        if(filterHca$numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & filterHca$numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
          ## filter valid
          output$hcaFilteredPrecursors <- renderText({
            print(paste("update output$hcaFilteredPrecursors ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
            paste("Number of filtered MS1 features: ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, globalMs2Filter, sep = "")
          })
        } else {
          ## filter invalid
          
          ## update info
          if(filterHca$numberOfPrecursorsFiltered == 0){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors # = 0", sep = ""))
              paste("There are no MS\u00B9 features which fulfill the given criteria.", globalMs2Filter, sep = "")
            })
          }
          if(filterHca$numberOfPrecursorsFiltered > 0 & filterHca$numberOfPrecursorsFiltered < minimumNumberOfPrecursorsForHca){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors 0 < # < ", minimumNumberOfPrecursorsForHca, sep = ""))
              paste("There are only ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, " MS\u00B9 features which fulfill the given criteria. There must be at least more than five MS\u00B9 features to proceed.", globalMs2Filter, sep = "")
            })
          }
          if(filterHca$numberOfPrecursorsFiltered > maximumNumberOfPrecursorsForHca){
            output$hcaFilteredPrecursors <- renderText({
              print(paste("update output$hcaFilteredPrecursors # > ", maximumNumberOfPrecursorsForHca, sep = ""))
              paste("There are ", filterHca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, " MS\u00B9 features which fulfill the given criteria. There must be at most ", maximumNumberOfPrecursorsForHca, " MS\u00B9 features to proceed.", globalMs2Filter, sep = "")
            })
          }
        }
      }
    }
    updatePcaFilterInformation <- function(){
      if(is.null(filterPca)){
        output$pcaFilteredPrecursors <- renderText({
          print(paste("update output$pcaFilteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values.", sep = "")
        })
      } else {
        globalMs2Filter <- ifelse(filterGlobal$numberOfPrecursorsFiltered != dataList$numberOfPrecursors, paste("\n(Global MS/MS filter: ", filterGlobal$numberOfPrecursorsFiltered, " / ", dataList$numberOfPrecursors, " precursors)", sep = ""), "")
        if(filterPca$numberOfPrecursorsFiltered > 0){
          output$pcaFilteredPrecursors <- renderText({
            print(paste("update output$pcaFilteredPrecursors", sep = ""))
            paste("Number of filtered MS1 features: ", filterPca$numberOfPrecursorsFiltered, " / ", filterGlobal$numberOfPrecursorsFiltered, globalMs2Filter, sep = "")
          })
        } else {
          output$pcaFilteredPrecursors <- renderText({
            print(paste("update output$pcaFilteredPrecursors", sep = ""))
            paste("There are no MS\u00B9 features which fulfill the given criteria.", globalMs2Filter, sep = "")
          })
        }
      }
    }
    updateSearchInformation <- function(){
      if(!state$filterSearchActive)
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo inactive search", sep = ""))
          paste("Please search for MS\u00B9 features", sep = "")
        })
      if(state$filterSearchActive & is.null(filterSearch))
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo invalid search", sep = ""))
          paste("There are invalid or missing search values", sep = "")
        })
      if(state$filterSearchActive & !is.null(filterSearch)){
        str1 <- ""
        str2 <- ""
        if(state$showHCAplotPanel & !is.null(listForTable_Search_HCA))
          str1 <- paste(length(listForTable_Search_HCA$precursorSet), " in HCA", sep = "")
        if(state$showPCAplotPanel & !is.null(listForTable_Search_PCA))
          str2 <- paste(length(listForTable_Search_PCA$precursorSet), " in PCA", sep = "")
        
        if(nchar(str1) > 0 & nchar(str2) > 0)
          val <- paste(str1, str2, sep = ", ")
        if(nchar(str1) > 0 & !(nchar(str2) > 0))
          val <- str1
        if(!(nchar(str1) > 0) & nchar(str2) > 0)
          val <- str2
        if(!(nchar(str1) > 0) & !(nchar(str2) > 0))
          val <- "None"
        
        output$searchInfo <- renderText({
          print(paste("update output$searchInfo", sep = ""))
          paste("Number of hits among MS\u00B9 features: ", val, sep = "")
        })
      }
    }
    ## perform filtering
    doPerformFiltering <- function(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, preFilter = NULL){
      suppressWarnings(
        doPerformFiltering_impl(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, preFilter)
      )
    }
    doPerformFiltering_impl <- function(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, preFilter = NULL){
      print(paste("Observe applyFilters1", "gs", paste(groupSet, collapse = "-"), "a", filter_average, "lfc", filter_lfc, "ms2_1", filter_ms2_masses1, "ms2_2", filter_ms2_masses2, "ms2_3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "ig", includeIgnoredPrecursors))
      
      groupSetOriginal                 <- groupSet
      sampleSetOriginal                <- sampleSet
      filterBySamplesOriginal          <- filterBySamples
      filter_averageOriginal           <- filter_average
      filter_lfcOriginal               <- filter_lfc
      filter_ms2_masses1Original       <- filter_ms2_masses1
      filter_ms2_masses2Original       <- filter_ms2_masses2
      filter_ms2_masses3Original       <- filter_ms2_masses3
      filter_ms2_ppmOriginal           <- filter_ms2_ppm
      filter_ms1_massesOriginal        <- filter_ms1_masses
      filter_ms1_ppmOriginal           <- filter_ms1_ppm
      includeIgnoredPrecursorsOriginal <- includeIgnoredPrecursors
      
      #################################################
      ## parse inputs
      print(paste("Observe applyFilters2", "gs", paste(groupSet, collapse = "-"), "ss", paste(sampleSet, collapse = "-"), "fbss", filterBySamples, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "i", includeIgnoredPrecursors))
      print(paste("Observe applyFilters3", "gs", is.null(groupSet), "ss", is.null(sampleSet), "gs", is.null(filterBySamples), "a", is.null(filter_average), "lfc", is.null(filter_lfc), "ms2 1", is.null(filter_ms2_masses1), "ms2 2", is.null(filter_ms2_masses2), "ms2 3", is.null(filter_ms2_masses3), "ppm", is.null(filter_ms2_ppm), "i", includeIgnoredPrecursors))
      
      #################################################
      ## sanity checks
      if(all(!is.null(filter_lfc), !is.na(filter_lfc), filter_lfc != 0) & length(groupSet) != 2)
        stop("lfc filter for not exactly two groups")
      
      #################################################
      ## check for errors in inputs amd process ms2
      error <- FALSE
      if(any(is.null(groupSet), is.na(groupSet), length(groupSet) == 0, any(nchar(groupSet) == 0)))
        error <- TRUE
      
      if(any(is.null(filter_average), is.na(filter_average), length(filter_average) == 0, nchar(filter_average) == 0))
        filter_average <- NULL
      else{
        filter_average <- as.numeric(filter_average)
        error <- error | is.na(filter_average)
      }
      
      if(any(is.null(filter_lfc), is.na(filter_lfc), length(filter_lfc) == 0, nchar(filter_lfc) == 0))
        filter_lfc <- NULL
      else{
        filter_lfc <- as.numeric(filter_lfc)
        error <- error | is.na(filter_lfc)
      }
      
      if(any(is.null(filter_ms2_masses1), is.na(filter_ms2_masses1), length(filter_ms2_masses1) == 0, nchar(filter_ms2_masses1) == 0))
        filter_ms2_masses1 <- NULL
      else{
        ms2Masses <- strsplit(x = filter_ms2_masses1, split = "[,; ]+")[[1]]
        filter_ms2_masses1 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses1[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses1))
      }
      if(any(is.null(filter_ms2_masses2), is.na(filter_ms2_masses2), length(filter_ms2_masses2) == 0, nchar(filter_ms2_masses2) == 0))
        filter_ms2_masses2 <- NULL
      else{
        ms2Masses <- strsplit(x = filter_ms2_masses2, split = "[,; ]+")[[1]]
        filter_ms2_masses2 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses2[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses2))
      }
      if(any(is.null(filter_ms2_masses3), is.na(filter_ms2_masses3), length(filter_ms2_masses3) == 0, nchar(filter_ms2_masses3) == 0))
        filter_ms2_masses3 <- NULL
      else{
        ms2Masses <- strsplit(x = filter_ms2_masses3, split = "[,; ]+")[[1]]
        filter_ms2_masses3 <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses3[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses3))
      }
      
      if(any(is.null(filter_ms2_ppm), is.na(filter_ms2_ppm), length(filter_ms2_ppm) == 0, nchar(filter_ms2_ppm) == 0))
        filter_ms2_ppm <- NULL
      else{
        filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
        error <- error | is.na(filter_ms2_ppm)
      }
      
      if(any(is.null(filter_ms1_masses), is.na(filter_ms1_masses), length(filter_ms1_masses) == 0, nchar(filter_ms1_masses) == 0))
        filter_ms1_masses <- NULL
      else{
        ms1Masses <- strsplit(x = filter_ms1_masses, split = "[,; ]+")[[1]]
        filter_ms1_masses <- vector(mode = "numeric", length = length(ms1Masses))
        for(idx in 1:length(ms1Masses))
          filter_ms1_masses[[idx]] <- as.numeric(ms1Masses[[idx]])
        error <- error | any(is.na(filter_ms1_masses))
      }
      
      if(any(is.null(filter_ms1_ppm), is.na(filter_ms1_ppm), length(filter_ms1_ppm) == 0, nchar(filter_ms1_ppm) == 0))
        filter_ms1_ppm <- NULL
      else{
        filter_ms1_ppm <- as.numeric(filter_ms1_ppm)
        error <- error | is.na(filter_ms1_ppm)
      }
      
      ## sanity check
      error <- error | (!is.null(filter_ms1_masses) & any(is.null(filter_ms1_ppm), is.na(filter_ms1_ppm)))
      error <- error | (any(!is.null(filter_ms2_masses1), !is.null(filter_ms2_masses2), !is.null(filter_ms2_masses3)) & any(is.null(filter_ms2_ppm), is.na(filter_ms2_ppm)))
      
      print(paste("Observe applyFilters4", "e", error, "gs", paste(groupSet, collapse = "-"), "ss", paste(sampleSet, collapse = "-"), "fbss", filterBySamples, "a", filter_average, "lfc", filter_lfc, "ms2 1", filter_ms2_masses1, "ms2 2", filter_ms2_masses2, "ms2 3", filter_ms2_masses3, "ppm", filter_ms2_ppm, "i", includeIgnoredPrecursors))
      
      ## collect ms2 masses
      filterList_ms2_masses <- list()
      if(!is.null(filter_ms2_masses1))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses1
      if(!is.null(filter_ms2_masses2))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses2
      if(!is.null(filter_ms2_masses3))
        filterList_ms2_masses[[length(filterList_ms2_masses) + 1]] <- filter_ms2_masses3
      
      #################################################
      ## do filtering and wrap results
      resultObj <- list()
      resultObj$error  <- error
      
      if(error){
        resultObj$filter <- NULL
      } else {
        filterHere <- filterData(
          dataList = dataList, 
          #groupOne = groupOne, groupTwo = groupTwo, 
          groups = groupSet, sampleSet, filterBySamples, filter_average = filter_average, filter_lfc = filter_lfc, 
          filterList_ms2_masses = filterList_ms2_masses, filter_ms2_ppm = filter_ms2_ppm, 
          filter_ms1_masses = filter_ms1_masses, filter_ms1_ppm = filter_ms1_ppm,
          includeIgnoredPrecursors = includeIgnoredPrecursors,
          progress = FALSE
        )
        
        ## set original values
        if(is.null(groupSetOriginal)){
          filterHere$groupSetOriginal <- list()
        } else {
          filterHere$groupSetOriginal  <- groupSetOriginal
          filterHere$sampleSetOriginal <- sampleSetOriginal
          filterHere$filterBySamplesOriginal <- filterBySamplesOriginal
        }
        #filterHere$groupSetOriginal                 <- ifelse(test = is.null(groupSetOriginal),                 yes = NA, no = groupSetOriginal)
        filterHere$filter_averageOriginal           <- ifelse(test = is.null(filter_averageOriginal),           yes = 0,  no = filter_averageOriginal)
        filterHere$filter_lfcOriginal               <- ifelse(test = is.null(filter_lfcOriginal),               yes = 0,  no = filter_lfcOriginal)
        if(is.null(filter_ms2_masses1Original)){
          filterHere$filter_ms2_masses1Original <- list()
        } else {
          filterHere$filter_ms2_masses1Original <- filter_ms2_masses1Original
        }
        #filterHere$filter_ms2_masses1Original       <- ifelse(test = is.null(filter_ms2_masses1Original),       yes = "", no = filter_ms2_masses1Original)
        if(is.null(filter_ms2_masses2Original)){
          filterHere$filter_ms2_masses2Original <- list()
        } else {
          filterHere$filter_ms2_masses2Original <- filter_ms2_masses2Original
        }
        #filterHere$filter_ms2_masses2Original       <- ifelse(test = is.null(filter_ms2_masses2Original),       yes = "", no = filter_ms2_masses2Original)
        if(is.null(filter_ms2_masses3Original)){
          filterHere$filter_ms2_masses3Original <- list()
        } else {
          filterHere$filter_ms2_masses3Original <- filter_ms2_masses3Original
        }
        #filterHere$filter_ms2_masses3Original       <- ifelse(test = is.null(filter_ms2_masses3Original),       yes = "", no = filter_ms2_masses3Original)
        filterHere$filter_ms2_ppmOriginal           <- ifelse(test = is.null(filter_ms2_ppmOriginal),           yes = "", no = filter_ms2_ppmOriginal)
        if(is.null(filter_ms1_massesOriginal)){
          filterHere$filter_ms1_massesOriginal <- list()
        } else {
          filterHere$filter_ms1_massesOriginal <- filter_ms1_massesOriginal
        }
        #filterHere$filter_ms1_massesOriginal        <- ifelse(test = is.null(filter_ms1_massesOriginal),        yes = "", no = filter_ms1_massesOriginal)
        filterHere$filter_ms1_ppmOriginal           <- ifelse(test = is.null(filter_ms1_ppmOriginal),           yes = "", no = filter_ms1_ppmOriginal)
        filterHere$includeIgnoredPrecursorsOriginal <- ifelse(test = is.null(includeIgnoredPrecursorsOriginal), yes = NA, no = includeIgnoredPrecursorsOriginal)
        
        resultObj$error  <- error
        resultObj$filter <- filterHere
        print(paste("Observe applyFilters5", "n", resultObj$filter$numberOfPrecursorsFiltered))
      }
      
      return(resultObj)
    }
    processSearchFilterResult <- function(resultObj){
      state$filterSearchActive <<- TRUE
      #################################################
      ## info / error
      if(resultObj$error){
        filterSearch <<- NULL
        state$searchfilterValid <<- FALSE
        selectionBySearch(NULL)
        updateSearchInformation()
      } else {
        filterSearch <<- resultObj$filter
        state$searchfilterValid <<- TRUE
        selectionBySearch(filterSearch$filter)
        updateSearchInformation()
        
        #################################################
        ## update plots
        if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "update search", withHeatmap = TRUE)
        if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "update search")
      }
    }
    checkHcaFilterValidity <- function(numberOfPrecursorsFiltered){
      if(numberOfPrecursorsFiltered >= minimumNumberOfPrecursorsForHca & numberOfPrecursorsFiltered <= maximumNumberOfPrecursorsForHca){
        ## filter valid
        print(paste("Observe applyFilters ", minimumNumberOfPrecursorsForHca, " <= # <= ", maximumNumberOfPrecursorsForHca, sep = ""))
        
        shinyjs::enable("drawHCAplots")
        #enableActionButton(session, "drawHCAplots")
        state$hcaFilterValid <<- TRUE
      } else {
        ## filter invalid
        
        shinyjs::disable("drawHCAplots")
        #disableActionButton(session, "drawHCAplots")
        state$hcaFilterValid <<- FALSE
      }
    }
    checkPcaFilterValidity <- function(numberOfPrecursorsFiltered){
      if(numberOfPrecursorsFiltered > 0){
        ## filter valid
        print(paste("Observe applyFilters # > 0", sep = ""))
        
        shinyjs::enable("drawPCAplots")
        state$pcaFilterValid <<- TRUE
      } else {
        ## filter invalid
        print(paste("Observe applyFilters # = 0", sep = ""))
        
        shinyjs::disable("drawPCAplots")
        state$pcaFilterValid <<- FALSE
      }
    }
    ## input
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
    ## plots
    doClearPlots <- function(){
      output$plotDendrogram <- renderPlot({
        print(paste("reset output$plotDendrogram"))
        NULL
      })
      output$plotHeatmap <- renderPlot({
        print(paste("reset output$plotHeatmap"))
        NULL
      })
      output$plotHeatmapLegend <- renderPlot({
        print(paste("reset output$plotHeatmapLegend"))
        NULL
      })
      output$plotMS2 <- renderPlot({
        print(paste("reset output$plotMS2"))
        NULL
      })
      output$plotPcaScores <- renderPlot({
        print(paste("reset output$plotPcaScores"))
        NULL
      })
      output$plotPcaLoadings <- renderPlot({
        print(paste("reset output$plotPcaLoadings"))
        NULL
       })
      output$plotAnnoLegendHCA <- renderPlot({
        print(paste("reset output$plotAnnoLegendHCA"))
        NULL
      })
      output$plotAnnoLegendPCA <- renderPlot({
        print(paste("reset output$plotAnnoLegendPCA"))
        NULL
      })
      output$plotMS2Legend <- renderPlot({
        print(paste("reset output$plotMS2Legend"))
        NULL
      })
      output$plotFragmentDiscriminativityLegend <- renderPlot({
        print(paste("reset output$plotFragmentDiscriminativityLegend"))
        NULL
      })
    }
    drawMS2Plot <- function(consoleInfo = NULL){
      output$plotMS2 <- renderPlot({
        print(paste("### MS2 ###", consoleInfo))
        drawMS2PlotImpl()
      })
    }
    drawMS2PlotImpl <- function(){
      calcPlotMS2(
        dataList = dataList, 
        fragmentsX = fragmentsX, 
        fragmentsY = fragmentsY, 
        fragmentsColor = fragmentsColor, 
        fragmentsDiscriminativity = fragmentsDiscriminativity, 
        fragmentsX_02 = fragmentsXhovered, 
        fragmentsY_02 = fragmentsYhovered, 
        fragmentsColor_02 = fragmentsColorHovered, 
        fragmentsDiscriminativity_02 = fragmentsDiscriminativityHovered, 
        xInterval = ms2PlotRange$xInterval, 
        selectedFragmentIndex = selectionFragmentSelectedFragmentIndex
      )
    }
    drawDendrogramPlot <- function(consoleInfo = NULL, withHeatmap = FALSE){
      output$plotDendrogram <- renderPlot({
      #output$plotDendrogram <- renderPlotly({
        print(paste("### den ###", consoleInfo))
        drawDendrogramPlotImpl()
      })
      
      if(!withHeatmap)
        return()
      
      ## plotly: remove
      output$plotHeatmap <- renderPlot({
        print(paste("### hea ### update range output$plotHeatmap"))
        drawHeatmapPlotImpl()
      })
    }
    drawDendrogramPlotImpl <- function(){
      resultObj <- calcPlotDendrogram(
        dataList = dataList, 
        filter = filterHca$filter, 
        clusterDataList = clusterDataList, 
        annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, 
        annoPresentColorsList = dataList$annoPresentColorsList, 
        distanceMeasure = currentDistanceMatrixObj$distanceMeasure, 
        selectionFragmentTreeNodeSet = selectionFragmentTreeNodeSet,
        selectionAnalysisTreeNodeSet = selectionAnalysisTreeNodeSet,
        selectionSearchTreeNodeSet = selectionSearchTreeNodeSet,
        showClusterLabels = state$showClusterLabels, 
        hcaPrecursorLabels = state$hcaPrecursorLabels, 
        xInterval = dendrogramPlotRange$xInterval
      )
      
      dendrogramPlotRange <- par("usr")
      dendrogramPlotRangeY <<- dendrogramPlotRange[[4]] - dendrogramPlotRange[[3]]
      
      state$annotationsHca <<- resultObj
      state$annotationLegendHeightHca <<- annoLegendEntryHeight * (length(state$annotationsHca$setOfAnnotations) + 1)
    }
    drawDendrogramPlotImpl_forPlotly <- function(){
      heatmapContent <- state$heatmapContent
      
      ## heatmap
      numberOfGroups <- -1
      switch(heatmapContent,
             "Log-fold-change"={## log-fold-change
               numberOfGroups <- 3
             },
             "Abundance by group"={## groups
               numberOfGroups <- length(dataList$groups)
             },
             "Abundance by sample"={## samples
               numberOfGroups <- length(dataList$dataColumnsNameFunctionFromGroupNames(groups = groups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame)))
             },
             {## unknown state
               stop(paste("Unknown heatmapContent value", heatmapContent))
             }
      )## end switch
      
      ## heigth per row
      heatmapHeightPerRow <- -1
      if(numberOfGroups <= maximumheatmapHeightRowCount){
        heatmapHeightPerRow <- maximumheatmapHeightPerRow
      } else if(numberOfGroups >= minimumheatmapHeightRowCount){
        heatmapHeightPerRow <- minimumheatmapHeightPerRow
      } else {
        heatmapHeightPerRow <- 
          minimumheatmapHeightPerRow + 
          (numberOfGroups    - maximumheatmapHeightRowCount) / 
          (minimumheatmapHeightRowCount - maximumheatmapHeightRowCount) * 
          (maximumheatmapHeightPerRow - minimumheatmapHeightPerRow)
      }
      heatmapHeight <- heatmapHeightPerRow * numberOfGroups
      dendrogramHeatmapHeight <- 500 + heatmapHeight
      heatmapProportion <- heatmapHeight / dendrogramHeatmapHeight
      state$dendrogramHeatmapHeight <<- dendrogramHeatmapHeight
      
      ## heatmap selections
      selectedSelection <- state$selectedSelection
      
      resultObj <- calcPlotDendrogram(
        dataList = dataList, 
        filterObj = filterHca,#filter, 
        clusterDataList = clusterDataList, 
        #annoPresentAnnotationsList = dataList$annoPresentAnnotationsList, 
        #annoPresentColorsList = dataList$annoPresentColorsList, 
        distanceMeasure = currentDistanceMatrixObj$distanceMeasure, 
        showClusterLabels = state$showClusterLabels, 
        hcaPrecursorLabels = state$hcaPrecursorLabels, 
        selectionFragmentTreeNodeSet = selectionFragmentTreeNodeSet,
        selectionAnalysisTreeNodeSet = selectionAnalysisTreeNodeSet,
        selectionSearchTreeNodeSet = selectionSearchTreeNodeSet,
        selectedSelection = selectedSelection,
        heatmapContent = heatmapContent,
        heatmapOrdering = heatmapOrdering,
        heatmapProportion = heatmapProportion
        #xInterval = dendrogramPlotRange$xInterval
      )
      curveNumberToCurveName <<- resultObj$curveNumberToCurveName
      plotlyPlot             <- resultObj$plotlyPlot
      #columnsOfInterest      <- resultObj$columnsOfInterest
      
      dendrogramPlotRange <- par("usr")
      dendrogramPlotRangeY <<- dendrogramPlotRange[[4]] - dendrogramPlotRange[[3]]
      
      ## present annotations
      resultObjAnno <- getPrecursorColors(dataList = dataList, precursorSet = filterHca$filter)
      
      uniqueIndeces     <- which(!duplicated(resultObjAnno$setOfAnnotations))
      uniqueAnnotations <- resultObjAnno$setOfAnnotations[uniqueIndeces]
      uniqueColors      <- resultObjAnno$setOfColors[uniqueIndeces]
      
      state$annotationsHca <<- list(
        setOfAnnotations = uniqueAnnotations,
        setOfColors      = uniqueColors
      )
      state$annotationLegendHeightHca <<- annoLegendEntryHeight * (length(state$annotationsHca$setOfAnnotations) + 1)
      
      return(plotlyPlot)
    }
    drawHeatmapPlotImpl <- function(consoleInfo = NULL){
      #calcPlotHeatmap(dataList = dataList, filterObj = filterHca, clusterDataList = clusterDataList, xInterval = dendrogramPlotRange$xInterval)
      if(!is.null(consoleInfo)) print(paste("### hea ###", consoleInfo))
      
      ## selections and selection colors
      selectedTreeNodeSet <- NULL
      frameColor <- NULL
      
      selectedSelection <- state$selectedSelection
      if(!is.null(selectedSelection)){
        analysisSelection <- function(){
          selectedTreeNodeSet <- selectionAnalysisTreeNodeSet
          frameColor <- "blue"
        }
        fragmentSelection <- function(){
          selectedTreeNodeSet <- selectionFragmentTreeNodeSet
          frameColor <- "green"
        }
        searchSelection <- function(){
          selectedTreeNodeSet <- selectionSearchTreeNodeSet
          frameColor <- "red"
        }
        switch(selectedSelection,
               "Analysis_HCA"=analysisSelection(),
               "Analysis_PCA"=analysisSelection(),
               "Fragment_HCA"=fragmentSelection(),
               "Fragment_PCA"=fragmentSelection(),
               "Search_HCA"  =searchSelection(),
               "Search_PCA"  =searchSelection(),
               {## unknown state
                 stop(paste("Unknown selectedSelection value", selectedSelection))
               }
        )## end switch
      }
      
      if(hcaHeatMapNew){
        heatmapContent  <- state$heatmapContent
        heatmapOrdering <- state$heatmapOrdering
        columnsOfInterest <- calcPlotHeatmap(
          dataList = dataList, 
          filterObj = filterHca, 
          clusterDataList = clusterDataList, 
          selectedTreeNodeSet = selectedTreeNodeSet, 
          frameColor = frameColor,
          heatmapContent = heatmapContent,
          heatmapOrdering = heatmapOrdering,
          xInterval = dendrogramPlotRange$xInterval
        )
        
        ## heigth per row
        heatmapHeightPerRow <- -1
        if(length(columnsOfInterest) <= maximumheatmapHeightRowCount){
          heatmapHeightPerRow <- maximumheatmapHeightPerRow
        } else if(length(columnsOfInterest) >= minimumheatmapHeightRowCount){
          heatmapHeightPerRow <- minimumheatmapHeightPerRow
        } else {
          heatmapHeightPerRow <- 
             minimumheatmapHeightPerRow + 
            (length(columnsOfInterest)    - maximumheatmapHeightRowCount) / 
            (minimumheatmapHeightRowCount - maximumheatmapHeightRowCount) * 
            (maximumheatmapHeightPerRow - minimumheatmapHeightPerRow)
        }
        
        ## TODO 12345
        #print(columnsOfInterest)
        #print(paste(heatmapHeightPerRow, length(columnsOfInterest), heatmapHeightPerRow * length(columnsOfInterest)))
        state$heatmapHeight <<- heatmapHeightPerRow * length(columnsOfInterest)
      } else {
        calcPlotHeatmapOld(
          dataList = dataList, 
          filterObj = filterHca, 
          clusterDataList = clusterDataList, 
          xInterval = dendrogramPlotRange$xInterval
        )
        state$heatmapHeight <<- heatmapHeightPerRow * 3
      }
    }
    drawPcaPlots <- function(consoleInfo = NULL){
      drawPcaScoresPlot(consoleInfo = consoleInfo)
      drawPcaLoadingsPlot(consoleInfo = consoleInfo)
    }
    drawPcaScoresPlot <- function(consoleInfo = NULL){
      output$plotPcaScores <- renderPlot({
        print(paste("### psp ###", consoleInfo))
        drawPcaScoresPlotImpl()
      })
    }
    drawPcaScoresPlotImpl <- function(){
      #resultObj <- calcPlotPCAscores(
      resultObj <- calcPlotPCAscores(
        pcaObj = pcaDataList$pcaObj, 
        dataList = dataList, 
        filterObj = filterPca, 
        pcaDimensionOne = pcaDataList$dimensionOne, 
        pcaDimensionTwo = pcaDataList$dimensionTwo, 
        showScoresLabels = state$showScoresLabels, 
        xInterval = pcaScoresPlotRange$xInterval, 
        yInterval = pcaScoresPlotRange$yInterval
      )
    }
    drawPcaLoadingsPlot <- function(consoleInfo = NULL){
      output$plotPcaLoadings <- renderPlot({
        print(paste("### plp ###", consoleInfo))
        drawPcaLoadingsPlotImpl()
      })
    }
    drawPcaLoadingsPlotImpl <- function(){
      resultObj <- calcPlotPCAloadings(
        pcaObj = pcaDataList$pcaObj, 
        dataList = dataList, 
        filter = filterPca$filter, 
        pcaDimensionOne = pcaDataList$dimensionOne, 
        pcaDimensionTwo = pcaDataList$dimensionTwo, 
        selectionFragmentPcaLoadingSet = selectionFragmentPcaLoadingSet,
        selectionAnalysisPcaLoadingSet = selectionAnalysisPcaLoadingSet,
        selectionSearchPcaLoadingSet   = selectionSearchPcaLoadingSet,
        xInterval = pcaLoadingsPlotRange$xInterval, 
        yInterval = pcaLoadingsPlotRange$yInterval,
        loadingsLabels = state$loadingsLabels, 
        showLoadingsAbundance = state$showLoadingsAbundance, 
        showLoadingsFeaturesAnnotated   = state$showLoadingsFeaturesAnnotated,
        showLoadingsFeaturesUnannotated = state$showLoadingsFeaturesUnannotated,
        showLoadingsFeaturesSelected    = state$showLoadingsFeaturesSelected,
        showLoadingsFeaturesUnselected  = state$showLoadingsFeaturesUnselected
      )
      
      state$annotationsPca <<- resultObj
      state$annotationLegendHeightPca <<- annoLegendEntryHeight * (length(state$annotationsPca$setOfAnnotations) + 1)
    }
    drawHeatmapLegend <- function(consoleInfo = NULL){
      output$plotHeatmapLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawHeatmapLegendImpl()
      })
    }
    drawHeatmapLegendImpl <- function(){
      calcPlotHeatmapLegend(dataList = dataList)
    }
    drawAnnotationLegendHCA <- function(consoleInfo = NULL){
      output$plotAnnoLegendHCA <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcPlotAnnoLegend(state$annotationsHca$setOfAnnotations, state$annotationsHca$setOfColors)
      })
    }
    drawAnnotationLegendPCA <- function(consoleInfo = NULL){
      output$plotAnnoLegendPCA <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcPlotAnnoLegend(state$annotationsPca$setOfAnnotations, state$annotationsPca$setOfColors)
      })
    }
    drawScoresGroupsLegend <- function(consoleInfo = NULL){
      output$plotScoresGroupsLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        calcPlotScoresGroupsLegend(scoresGroups$groups, scoresGroups$colors)
      })
    }
    drawMS2Legend <- function(consoleInfo = NULL){
      output$plotMS2Legend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawMS2LegendImpl()
      })
    }
    drawMS2LegendImpl <- function(){
      calcPlotMS2Legend(dataList = dataList)
    }
    drawDendrogramLegend <- function(consoleInfo = NULL){
      output$calcPlotDendrogramLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawDendrogramLegendImpl()
      })
    }
    drawDendrogramLegendImpl <- function(){
      calcPlotDendrogramLegend()
    }
    drawFragmentDiscriminativityLegend <- function(consoleInfo = NULL){
      output$plotFragmentDiscriminativityLegend <- renderPlot({
        print(paste("### leg ###", consoleInfo))
        drawFragmentDiscriminativityLegendImpl()
      })
    }
    drawFragmentDiscriminativityLegendImpl <- function(){
      calcPlotDiscriminativityLegend()
    }
    drawFragmentPlot <- function(consoleInfo = NULL){
      output$fragmentPlot <- renderPlot({
        print(paste("### MS2 ### all", consoleInfo))
        drawFragmentPlotImpl()
      })
    }
    drawFragmentPlotImpl <- function(){
      plotFragments(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
    }
    ## plot range resets
    resetHcaPlotRange <- function(){
      dendrogramPlotRange$xMin <<- 1
      dendrogramPlotRange$xMax <<- filterHca$numberOfPrecursorsFiltered
      dendrogramPlotRange$xInterval <<- c(1, filterHca$numberOfPrecursorsFiltered)
      dendrogramPlotRange$xIntervalSize <<- filterHca$numberOfPrecursorsFiltered - 1
    }
    resetPcaPlotRange <- function(){
      minX <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
      maxX <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
      minY <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
      maxY <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
      
      if(any(is.na(c(minX, maxX, minY, maxY)))){
        minX <- -1
        maxX <- 1
        minY <- -1
        maxY <- 1
      }
      
      pcaScoresPlotRange$xMin <<- minX
      pcaScoresPlotRange$xMax <<- maxX
      pcaScoresPlotRange$xInterval <<- c(minX, maxX)
      pcaScoresPlotRange$xIntervalSize <<- maxX - minX
      pcaScoresPlotRange$yMin <<- minY
      pcaScoresPlotRange$yMax <<- maxY
      pcaScoresPlotRange$yInterval <<- c(minY, maxY)
      pcaScoresPlotRange$yIntervalSize <<- maxY - minY
      
      minX <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
      maxX <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
      minY <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
      maxY <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
      
      if(any(is.na(c(minX, maxX, minY, maxY)))){
        minX <- -1
        maxX <- 1
        minY <- -1
        maxY <- 1
      }
      
      pcaLoadingsPlotRange$xMin <<- minX
      pcaLoadingsPlotRange$xMax <<- maxX
      pcaLoadingsPlotRange$xInterval <<- c(minX, maxX)
      pcaLoadingsPlotRange$xIntervalSize <<- maxX - minX
      pcaLoadingsPlotRange$yMin <<- minY
      pcaLoadingsPlotRange$yMax <<- maxY
      pcaLoadingsPlotRange$yInterval <<- c(minY, maxY)
      pcaLoadingsPlotRange$yIntervalSize <<- maxY - minY
    }
    resetMS2PlotRange <- function(){
      ms2PlotRange$xMin <<- dataList$minimumMass
      ms2PlotRange$xMax <<- dataList$maximumMass
      ms2PlotRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      ms2PlotRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
    }
    resetMS2VsClassPlotRange <- function(){
      specVsClassRange$xMin <<- dataList$minimumMass
      specVsClassRange$xMax <<- dataList$maximumMass
      specVsClassRange$xInterval <<- c(dataList$minimumMass, dataList$maximumMass)
      specVsClassRange$xIntervalSize <<- dataList$maximumMass - dataList$minimumMass
    }
    ## annotation stuff
    addAnnotation <- function(precursorSet, annotationValue, annotationColor){
      if(is.na(match(x = annotationValue, table = dataList$annoPresentAnnotationsList))){
        ## new annotation
        annoIdx <- length(dataList$annoPresentAnnotationsList) + 1
        dataList$annoPresentAnnotationsList[[annoIdx]] <<- annotationValue
        dataList$annoPresentColorsList[[annoIdx]] <<- annotationColor
      }
      ## add
      for(precursor in precursorSet){
        if(annotationValue == artifact)
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
      updateAnnoGui(precursorSet)
      updateMS1FeatureTableGui(precursorSet)
      updatePlotsWithAnnotations()
      #updateAnnotationOverview()
    }
    removeAnnotation <- function(precursorSet, annotationValue){
      ## remove
      for(precursor in precursorSet){
        if(annotationValue == artifact)
          dataList$annoArrayIsArtifact[[precursor]] <<- FALSE
        else{
          idx <- match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])
          dataList$annoArrayOfLists[[precursor]] <<- dataList$annoArrayOfLists[[precursor]][-idx]
        }
      }
      
      ## remove anno completely?
      annoThere <- lapply(X = dataList$annoArrayOfLists, FUN = function(x){ match(x = annotationValue, table = x) })
      if(annotationValue != artifact & all(is.na(annoThere))){
        idx <- match(x = annotationValue, table = dataList$annoPresentAnnotationsList)
        dataList$annoPresentAnnotationsList <<- dataList$annoPresentAnnotationsList[-idx]
        dataList$annoPresentColorsList      <<- dataList$annoPresentColorsList[-idx]
      }
      
      ## update gui
      #print("removeAnnotation updateAnnoGui")
      updateAnnoGui(precursorSet)
      updateMS1FeatureTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    setArtifactState <- function(precursorSet, isArtifact){
      ## add
      dataList$annoArrayIsArtifact[precursorSet] <<- isArtifact
      
      ## update gui
      updateAnnoGui(precursorSet)
      updateMS1FeatureTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    setAnnotationPrimary <- function(precursorSet, annotationValue){
      ## remove
      for(precursor in precursorSet){
        idx <- match(x = annotationValue, table = dataList$annoArrayOfLists[[precursor]])
        dataList$annoArrayOfLists[[precursor]] <<- c(dataList$annoArrayOfLists[[precursor]][[idx]], dataList$annoArrayOfLists[[precursor]][-idx])
      }
      
      ## update gui
      #print("setAnnotation primary updateAnnoGui")
      updateAnnoGui(precursorSet)
      updateMS1FeatureTableGui(precursorSet)
      updatePlotsWithAnnotations()
    }
    commonAnnotations <- function(precursorSet){
      if(is.null(precursorSet))
        return(NULL)
      if(all(dataList$annoArrayIsArtifact[precursorSet]))
        return(artifact)
      
      ## at least one non-artifact precursor present
      intersection <- unlist(dataList$annoPresentAnnotationsList)
      for(precursor in precursorSet)
        if(!dataList$annoArrayIsArtifact[[precursor]])
          intersection <- intersect(x = intersection, y = unlist(dataList$annoArrayOfLists[[precursor]]))
      return(intersection)
    }
    ## table update
    updateAnnoGui <- function(precursorSet){
      ## annotation
      commonAnnos <- commonAnnotations(precursorSet)
      
      if(!is.null(commonAnnos))
        updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = commonAnnos)
      else
        updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[none]"))
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
    updateMS1FeatureTableGui <- function(precursorSet){
      print("updateMS1FeatureTableGui")
      ## table update with new annotations
      if(all(!is.null(selectionFragmentTreeNodeSet), precursorSet %in% listForTable_Fragment_HCA$precursorSet)){ ## HCA
        table$df_Fragment_HCA <<- createMS1FeatureTable(listForTable_Fragment_HCA, selectionFragmentHcaName)
        #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
        table_Fragment_HCA_id <<- ms1FeatureTableInputFieldIdCounter
      }
      if(all(!is.null(selectionFragmentPcaLoadingSet), precursorSet %in% listForTable_Fragment_PCA$precursorSet)){ ## PCA
        table$df_Fragment_PCA <<- createMS1FeatureTable(listForTable_Fragment_PCA, selectionFragmentPcaName)
        #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
        table_Fragment_PCA_id <<- ms1FeatureTableInputFieldIdCounter
      }
      
      if(all(!is.null(selectionAnalysisTreeNodeSet), precursorSet %in% listForTable_Analysis_HCA$precursorSet)){ ## HCA
        table$df_Analysis_HCA <<- createMS1FeatureTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
        #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
        table_Analysis_HCA_id <<- ms1FeatureTableInputFieldIdCounter
      }
      if(all(!is.null(selectionAnalysisPcaLoadingSet), precursorSet %in% listForTable_Analysis_PCA$precursorSet)){ ## PCA
        table$df_Analysis_PCA <<- createMS1FeatureTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
        #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
        table_Analysis_PCA_id <<- ms1FeatureTableInputFieldIdCounter
      }
      
      if(all(!is.null(selectionSearchTreeNodeSet), precursorSet %in% listForTable_Search_HCA$precursorSet)){ ## HCA
        table$df_Search_HCA <<- createMS1FeatureTable(listForTable_Search_HCA, selectionSearchHcaName)
        #output$dt_Search_HCA <- DT::renderDataTable(table$df_Search_HCA)
        table_Search_HCA_id <<- ms1FeatureTableInputFieldIdCounter
      }
      if(all(!is.null(selectionSearchPcaLoadingSet), precursorSet %in% listForTable_Search_PCA$precursorSet)){ ## PCA
        table$df_Search_PCA <<- createMS1FeatureTable(listForTable_Search_PCA, selectionSearchPcaName)
        #output$dt_Search_PCA <- DT::renderDataTable(table$df_Search_PCA)
        table_Search_PCA_id <<- ms1FeatureTableInputFieldIdCounter
      }
      ## update
      updateTableAssignment()
    }
    updateTableAssignment <- function(){
      print(paste("updateTableAssignment '", state$selectedSelection, "'", sep = ""))
      switch(state$selectedSelection, 
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
               state$precursorSetSelected <<- !is.null(listForTable_Fragment_HCA)
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
               print(paste("### unknown state$selectedSelection: '", state$selectedSelection, "'", sep = ""))
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
      
      state$selectedSelection <<- selectedSelection
      updateSelectedPrecursorSet()
    }
    updateSelectedPrecursorSet <- function(){
      print(paste("updateSelectionAssignment '", state$selectedSelection, "'", sep = ""))
      switch(state$selectedSelection, 
             "Analysis_HCA"={ 
               state$precursorSetSelected <<- !is.null(listForTable_Analysis_HCA)
               if(!is.null(listForTable_Analysis_HCA)) selectedPrecursorSet <<- listForTable_Analysis_HCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Analysis_PCA"={  
               state$precursorSetSelected <<- !is.null(listForTable_Analysis_PCA)
               if(!is.null(listForTable_Analysis_PCA)) selectedPrecursorSet <<- listForTable_Analysis_PCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Fragment_HCA"={  
               state$precursorSetSelected <<- !is.null(listForTable_Fragment_HCA)
               if(!is.null(listForTable_Fragment_HCA)) selectedPrecursorSet <<- listForTable_Fragment_HCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Fragment_PCA"={  
               state$precursorSetSelected <<- !is.null(listForTable_Fragment_PCA)
               if(!is.null(listForTable_Fragment_PCA)) selectedPrecursorSet <<- listForTable_Fragment_PCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Search_HCA"  ={  
               state$precursorSetSelected <<- !is.null(listForTable_Search_HCA)
               if(!is.null(listForTable_Search_HCA)) selectedPrecursorSet <<- listForTable_Search_HCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },"Search_PCA"  ={  
               state$precursorSetSelected <<- !is.null(listForTable_Search_PCA)
               if(!is.null(listForTable_Search_PCA)) selectedPrecursorSet <<- listForTable_Search_PCA$precursorSet
               else                                    selectedPrecursorSet <<- NULL
             },{
               print(paste("### unknown state$selectedSelection: '", state$selectedSelection, "'", sep = ""))
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
      selection <- state$selectedSelection
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
    updateAnnotationOverview <- function(){
      
      allAnnotations      <- unlist(dataList$annoPresentAnnotationsList)
      allAnnotationColors <- unlist(dataList$annoPresentColorsList)
      annoCounts <- sapply(X = allAnnotations, FUN = function(x){sum(unlist(lapply(X = dataList$annoArrayOfLists, FUN = function(y){any(y==x)})))})
      df <- data.frame(
        "Family" = allAnnotations,
        "Color"  = allAnnotationColors,
        "Count"  = annoCounts
      )
      output$familySelectionTable <- DT::renderDataTable(
        expr = datatable(df) %>% formatStyle(columns = "Color", target = "cell", backgroundColor = styleEqual(levels = allAnnotations, values = allAnnotationColors)),
        server = FALSE, escape = FALSE, rownames = FALSE,
        selection = list(mode = "single"),#, selected = selected), 
        options = list(
          #scrollY = "600px",
          scrollY = "30vh",
          preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
          drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
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
      #dataList$annoArrayOfLists[[precursor]] <<- c(dataList$annoArrayOfLists[[precursor]], annotationValue)
    }
    obsFamilySelectionTable_rows_selected <- observeEvent(input$familySelectionTable_rows_selected, {
      print(paste("Observe familySelectionTable_rows_selected", input$familySelectionTable_rows_selected))
      
      ## TODO
      
    })
    
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
    createCheckboxInputFields2 <- function(FUN, id, values, tableCounter) {
      ## running id
      id <- paste(id, tableCounter, sep = "_")
      
      ## create a character vector of shiny inputs
      inputs <- character(length(values))
      for (i in 1:length(values)){
        itemId    <- paste(id, "_", i, sep = "")
        inputs[[i]] <- as.character(FUN(inputId = itemId, label = NULL, value = values[[i]]))
        
        ## trigger event on button-click
        lapply(c(itemId), function(i){
          observeEvent(input[[i]], {
            sampleExcludeClicked(i)
          })
        })
      }
      return(inputs)
    }
    createActionButtonInputFields <- function(FUN, id, itemCount, icon, tableCounter) {
      ## running id
      id <- paste(id, tableCounter, sep = "_")
      
      ## create a character vector of shiny inputs
      inputs <- character(length = itemCount)
      for (i in seq_along(inputs)){
        itemId    <- paste(id, "_", i, sep = "")
        inputs[[i]] <- as.character(FUN(inputId = itemId, label = NULL, icon = icon))
      }
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
    
    preprocessClassPlot <- function(frequentFragments, characteristicFragments){
      frequentMasses       <- as.numeric(names(frequentFragments)) 
      characteristicMasses <- as.numeric(names(characteristicFragments))
      
      masses_class <- unique(c(frequentMasses, characteristicMasses))
      
      frequency_class <- rep(x = 0.01, times = length(masses_class))
      frequency_class[match(x = frequentMasses, table = masses_class)] <- unname(frequentFragments)
      
      characteristics_class <- rep(x = 0., times = length(masses_class))
      characteristics_class[match(x = characteristicMasses, table = masses_class)] <- unname(characteristicFragments)
      
      classDataColorMapFragmentData  <- makecmap(
        x = c(0, 1), n = 100, 
        colFn = colorRampPalette(c('grey', 'black'))
      )
      colors_class <- cmap(x = characteristics_class, map = classDataColorMapFragmentData)
      
      returnObj <- list(
        masses_class = masses_class,
        frequency_class = frequency_class,
        colors_class = colors_class
      )
      return(returnObj)
    }
    preprocessSpectrumVsClassPlot <- function(dataList, precursorIndex, masses_class){
      resultObj <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)
      masses_spec <- resultObj$fragmentMasses
      intensity_spec <- resultObj$fragmentAbundances
      
      tolerance <- .Machine$double.eps ^ 0.5 ## default in function all.equal
      
      matchingMassRowIndeces <- which(
        apply(X = outer(X = mappingSpectraToClassDf$SpectraMasses, Y = masses_spec , FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)}) &
        apply(X = outer(X = mappingSpectraToClassDf$ClassMasses  , Y = masses_class, FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)})
      )
      #matchingMassRowIndeces <- which(
      #  mappingSpectraToClassDf$SpectraMasses %in% masses_spec & 
      #  mappingSpectraToClassDf$ClassMasses   %in% masses_class
      
      specIndeces  <- match(x = mappingSpectraToClassDf$SpectraMasses[matchingMassRowIndeces], table = masses_spec )
      classIndeces <- match(x = mappingSpectraToClassDf$ClassMasses  [matchingMassRowIndeces], table = masses_class)
      
      colors_spec <- rep(x = "grey", times = length(masses_spec))
      #colors_spec[specIndeces] <- colors_class[classIndeces]
      colors_spec[specIndeces] <- "black"
      
      numberOfMatchingMasses <- length(matchingMassRowIndeces)
      
      returnObj <- list(
        masses_spec = masses_spec,
        intensity_spec = intensity_spec,
        colors_spec = colors_spec,
        numberOfMatchingMasses = numberOfMatchingMasses
      )
      return(returnObj)
    }
    
    createPlotOutput <- function(dataList, id, frequentFragments, characteristicFragments, precursorIndeces, tableCounter) {
      
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
        returnObj <- preprocessSpectrumVsClassPlot(dataList, precursorIndex, masses_class)
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
        
        #print("############################")
        #print(paste(i,itemId))
        #print(plotAsString)
        
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
        
        #print(inputs_i[i])
      }
      
      returnObj <- list(
        numberOfMatchingMasses_i = numberOfMatchingMasses_i,
        inputs_i = inputs_i
      )
      
      return(returnObj)
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
    
    ## shinyBS uses jQuery to select page elements based on their id. When shiny passes the id attribute to the output page, it prepends a "#" to indicate an id. 
    ## There is nothing stopping you from expanding that id attribute to contain additional css selectors. for example, lets say your side table has an id of sideTable. 
    ## If you wanted to add a tooltip to the first column header you could do something like this:
    ## addTooltip(session, "sideTable th:nth-child(1)", "This is a tooltip")
    ## The problem is the elements of the table aren't there when the page first loads so there isn't anything for the selector to select... As it stands now, you would
    ## need to write some sort of observer that fires after the table is rendered in order for jQuery to find the table cells.
    createSampleTable <- function(){
      sampleTableInputFieldIdCounter <<- sampleTableInputFieldIdCounter + 1
      
      groupSampleDataFrame <- dataList$groupSampleDataFrame[sampleOrder_tmp, c("Group", "Sample")]
      isExcluded <- sampleExclusion_tmp[sampleOrder_tmp]
      
      iconUp   <- icon(name = "chevron-up",   lib = "font-awesome")
      iconDown <- icon(name = "chevron-down", lib = "font-awesome")
      
      checkboxes   <- createCheckboxInputFields2(    FUN = checkboxInput, id = "Sample_Exclude", values = isExcluded, tableCounter = sampleTableInputFieldIdCounter)
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
      iconUp   <- icon(name = "chevron-up",   lib = "font-awesome")
      iconDown <- icon(name = "chevron-down", lib = "font-awesome")
      checkboxes <- createCheckboxInputFields(       FUN = checkboxInput, id = paste(type, "Ignore",   sep = "_"), values = isArtifact, tableCounter = ms1FeatureTableInputFieldIdCounter)
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
    ## selection stuff
    selectionByFragmentReset <- function(){
      selectionFragmentSelectedFragmentIndex <<- NULL
      
      if(!is.null(selectionFragmentTreeNodeSet)){ ## HCA
        selectionFragmentSelectedFragmentIndex <<- NULL
        selectionFragmentTreeNodeSet <<- NULL
        listForTable_Fragment_HCA <<- NULL
        table_Fragment_HCA_id <<- NULL
        table$df_Fragment_HCA <<- NULL
        #output$dt_Fragment_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionFragmentHcaName)
          updateSelectedPrecursorSet()
      }
      if(!is.null(selectionFragmentPcaLoadingSet)){ ## PCA
        selectionFragmentPcaLoadingSet <<- NULL
        listForTable_Fragment_PCA <<- NULL
        table_Fragment_PCA_id <<- NULL
        table$df_Fragment_PCA <<- NULL
        #output$dt_Fragment_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionFragmentPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByFragment <- function(minimumIndex){
      selectionFragmentSelectedFragmentIndex <<- minimumIndex
      
      fragmentMass  <- fragmentsX[[minimumIndex]]
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
        if(state$selectedSelection == selectionFragmentHcaName)
          updateSelectedPrecursorSet()
      } else {
        listForTable_Fragment_HCA <<- NULL
        table_Fragment_HCA_id <<- NULL
        table$df_Fragment_HCA <<- NULL
        #output$dt_Fragment_HCA <- DT::renderDataTable(table$df_Fragment_HCA)
        if(state$selectedSelection == selectionFragmentHcaName)
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
        if(state$selectedSelection == selectionFragmentPcaName)
          updateSelectedPrecursorSet()
      } else {
        listForTable_Fragment_PCA <<- NULL
        table_Fragment_PCA_id <<- NULL
        table$df_Fragment_PCA <<- NULL
        #output$dt_Fragment_PCA <- DT::renderDataTable(table$df_Fragment_PCA)
        if(state$selectedSelection == selectionFragmentPcaName)
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
        if(state$selectedSelection == selectionAnalysisHcaName)
          updateSelectedPrecursorSet()
      }
      if(!is.null(selectionAnalysisPcaLoadingSet)){ ## PCA
        selectionAnalysisPcaLoadingSet <<- NULL
        listForTable_Analysis_PCA <<- NULL
        table_Analysis_PCA_id <<- NULL
        table$df_Analysis_PCA <<- NULL
        #output$dt_Analysis_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByHca <- function(minimumLabel){
      selectionAnalysisTreeNodeSet <<- minimumLabel
      precursorSet <- getPrecursorSetFromTreeSelection(clusterDataList = clusterDataList, clusterLabel = minimumLabel)
      listForTable_Analysis_HCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorSet)
      table$df_Analysis_HCA <<- createMS1FeatureTable(listForTable_Analysis_HCA, selectionAnalysisHcaName)
      #output$dt_Analysis_HCA <- DT::renderDataTable(table$df_Analysis_HCA)
      if(state$selectedSelection == selectionAnalysisHcaName)
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
        if(state$selectedSelection == selectionAnalysisPcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionAnalysisPcaLoadingSet <<- NULL
        listForTable_Analysis_PCA <<- NULL
        table_Analysis_PCA_id <<- NULL
        table$df_Analysis_PCA <<- NULL
        #output$dt_Analysis_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisPcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionByPca <- function(minimumIndex){
      precursorIndex <- filterPca$filter[[minimumIndex]]
      
      selectionAnalysisPcaLoadingSet <<- minimumIndex
      listForTable_Analysis_PCA <<- getTableFromPrecursorSet(dataList = dataList, precursorSet = precursorIndex)
      table$df_Analysis_PCA <<- createMS1FeatureTable(listForTable_Analysis_PCA, selectionAnalysisPcaName)
      #output$dt_Analysis_PCA <- DT::renderDataTable(table$df_Analysis_PCA)
      if(state$selectedSelection == selectionAnalysisPcaName)
        updateSelectedPrecursorSet()
      
      if(state$showHCAplotPanel)
        selectionByAnalysisInitHca(precursorIndex)
      
      if(input$changeSelection != selectionAnalysisName){
        updateRadioButtons(session = session, inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), inline = TRUE, selected = selectionAnalysisName)
        updateSelectedSelection()
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
        if(state$selectedSelection == selectionAnalysisHcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionAnalysisTreeNodeSet <<- NULL
        listForTable_Analysis_HCA <<- NULL
        table_Analysis_HCA_id <<- NULL
        table$df_Analysis_HCA <<- NULL
        #output$dt_Analysis_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionAnalysisHcaName)
          updateSelectedPrecursorSet()
      }
    }
    selectionBySearchReset <- function(){
      if(!is.null(selectionSearchTreeNodeSet)){ ## HCA
        selectionSearchTreeNodeSet <<- NULL
        listForTable_Search_HCA <<- NULL
        table_Search_HCA_id <<- NULL
        table$df_Search_HCA <<- NULL
        if(state$selectedSelection == selectionSearchHcaName)
          updateSelectedPrecursorSet()
      }
      if(!is.null(selectionSearchPcaLoadingSet)){ ## HCA
        selectionSearchPcaLoadingSet <<- NULL
        listForTable_Search_PCA <<- NULL
        table_Search_PCA_id <<- NULL
        table$df_Search_PCA <<- NULL
        if(state$selectedSelection == selectionSearchHcaName)
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
        if(state$selectedSelection == selectionSearchHcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionSearchTreeNodeSet <<- NULL
        listForTable_Search_HCA <<- NULL
        table_Search_HCA_id <<- NULL
        table$df_Search_HCA <<- NULL
        #output$dt_Search_HCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionSearchHcaName)
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
        if(state$selectedSelection == selectionSearchPcaName)
          updateSelectedPrecursorSet()
      } else {
        selectionSearchPcaLoadingSet <<- NULL
        listForTable_Search_PCA <<- NULL
        table_Search_PCA_id <<- NULL
        table$df_Search_PCA <<- NULL
        #output$dt_Search_PCA <- DT::renderDataTable(NULL)
        if(state$selectedSelection == selectionSearchPcaName)
          updateSelectedPrecursorSet()
      }
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
          state$plotHcaShown <<- TRUE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
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
          state$plotPcaShown <<- TRUE
          state$plotAnnotationShown <<- FALSE
        } else {
          state$plotHcaShown <<- FALSE
          state$plotPcaShown <<- FALSE
          state$plotAnnotationShown <<- FALSE
        }
      }
      if(tabId == "Annotation"){
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
      if(tabId == "Input" & !initialGuiUpdatePerformed){
        ## initial gui update
        print(paste("update GUI initially", tabId))
        
        shinyjs::disable("importMs1Ms2Data")
        shinyjs::disable("importMs2Data")
        shinyjs::disable("loadProjectData")
        
        ## annotation classifier selection
        session$sendCustomMessage("disableButton", "doAnnotation")
        
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
      if(clearSelection == clearSelectionButtonValue)
        return()
      clearSelectionButtonValue <<- clearSelection
      
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
    obsFile <- observeEvent(input$matrixFile$datapath, {
      filePath <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      print(paste("Observe file for data", fileName))
      if(!is.null(filePath))
        shinyjs::enable("loadProjectData")
      
      updateFileInputInfo()
    })
    obsLoadProjectData <- observeEvent(input$loadProjectData, {
      loadProjectData <- as.numeric(input$loadProjectData)
      
      print(paste("Observe loadProjectData", loadProjectData))
      
      #################################################
      ## check if button was hit
      if(loadProjectData == loadProjectDataButtonValue)
        return()
      loadProjectDataButtonValue <<- loadProjectData
      
      #################################################
      ## files
      filePath <- input$matrixFile$datapath
      loadProjectFile(filePath = filePath, buttonId = "loadProjectData")
    })
    obsLoadExampleData <- observeEvent(input$loadExampleData, {
      loadExampleData <- as.numeric(input$loadExampleData)
      
      print(paste("Observe loadExampleData", loadExampleData))
      
      #################################################
      ## check if button was hit
      if(loadExampleData == loadExampleDataButtonValue)
        return()
      loadExampleDataButtonValue <<- loadExampleData
      
      #################################################
      ## files
      filePath <- getFile("Project_file_showcase_annotated.csv.gz")
      loadProjectFile(filePath = filePath, buttonId = "loadExampleData")
    })
    loadProjectFile <- function(filePath, buttonId){
      fileName <- basename(filePath)
      #########################################################################################
      ## read data
      session$sendCustomMessage("disableButton", buttonId)
      
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
        session$sendCustomMessage("enableButton", buttonId)
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
      
      session$sendCustomMessage("enableButton", buttonId)
    }
    showErrorDialog <- function(msg){
      showDialog("An error occurred", msg)
    }
    showInfoDialog <- function(msg){
      showDialog("Information", msg)
    }
    showDialog <- function(title, msg){
      print("Show dialog")
      #output$infoPopupDialog <- renderUI({
      #  bsModal(id = "modalInfoPopupDialog", title = "Information", trigger = "", size = "large", HTML(msg))
      #})
      #toggleModal(session = session, modalId = "modalInfoPopupDialog", toggle = "open")
      showModal(modalDialog(title = title, HTML(msg)))
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
      importMs1Ms2Data <- as.numeric(input$importMs1Ms2Data)
      
      print(paste("Observe importMs1Ms2Data", importMs1Ms2Data))
      
      #################################################
      ## check if button was hit
      if(importMs1Ms2Data == importMs1Ms2DataButtonValue)
        return()
      importMs1Ms2DataButtonValue <<- importMs1Ms2Data
      
      importData(TRUE)
    })
    obsImportMs2Data <- observeEvent(input$importMs2Data, {
      importMs2Data <- as.numeric(input$importMs2Data)
      
      print(paste("Observe importMs2Data", importMs2Data))
      
      #################################################
      ## check if button was hit
      if(importMs2Data == importMs2DataButtonValue)
        return()
      importMs2DataButtonValue <<- importMs2Data
      
      importData(FALSE)
    })
    importData <- function(importMS1andMS2data){
      session$sendCustomMessage("disableButton", "importMs1Ms2Data")
      session$sendCustomMessage("disableButton", "importMs2Data")
      
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
      
      print(parameterSet)
      
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
        matrixRows <- resultObj$matrixRows
        matrixCols <- resultObj$matrixCols
        matrixVals <- resultObj$matrixVals
        
        matrixRows <- c(matrixRows, 1)
        matrixCols <- c(matrixCols, 1)
        matrixVals <- c(matrixVals, serializeParameterSet(parameterSet))
        
        ## TODO performance 25s
        ## convert matrix to dataframe
        numberOfRows    <- max(matrixRows)
        numberOfColumns <- max(matrixCols)
        
        lines <- vector(mode = "character", length = numberOfRows)
        for(rowIdx in seq_len(numberOfRows)){
          indeces <- matrixRows == rowIdx
          tokens  <- vector(mode = "character", length = numberOfColumns)
          tokens[matrixCols[indeces]] <- matrixVals[indeces]
          lines[[rowIdx]] <- paste(tokens, collapse = "\t")
        }
        
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
      
      
      msg <- paste(
        "The data import was successful.",
        "<br>",
        "<br>",
        ## spectra
        resultObj$numberOfParsedSpectra, " / ", resultObj$numberOfSpectraOriginal, " spectra were imported successfully.",
        ifelse(test = resultObj$numberOfParsedSpectra < resultObj$numberOfSpectraOriginal, yes = paste(" (",paste( Filter(nchar, c(
          ifelse(test = resultObj$numberOfSpectraDiscardedDueToNoPeaks      > 0, yes = paste(resultObj$numberOfSpectraDiscardedDueToNoPeaks,      " empty", sep = ""), no = ""), 
          ifelse(test = resultObj$numberOfSpectraDiscardedDueToMaxIntensity > 0, yes = paste(resultObj$numberOfSpectraDiscardedDueToMaxIntensity, " low intensity", sep = ""), no = ""), 
          ifelse(test = resultObj$numberOfSpectraDiscardedDueToTooHeavy     > 0, yes = paste(resultObj$numberOfSpectraDiscardedDueToTooHeavy,     " too heavy", sep = ""), no = "")
        )), collapse = ", "), ")", sep = ""), no = ""),
        "<br>",
        ## mapping
        resultObj$numberOfPrecursors, " / ", resultObj$numberOfParsedSpectra, " spectra were successfully mapped to MS\u00B9 features.", 
        ifelse(test = resultObj$numberOfPrecursors < resultObj$numberOfParsedSpectra, yes = paste(" (",paste( Filter(nchar, c(
          #ifelse(test = resultObj$numberOfUnmappedPrecursorsMz > 0, yes = paste(resultObj$numberOfUnmappedPrecursorsMz, " with m/z deviation", sep = ""), no = ""), 
          #ifelse(test = resultObj$numberOfUnmappedPrecursorsRt > 0, yes = paste(resultObj$numberOfUnmappedPrecursorsRt, " with RT deviation",  sep = ""), no = "")
          ifelse(test = resultObj$numberOfUnmappedSpectra > 0, yes = paste(resultObj$numberOfUnmappedSpectra, " unmapped",  sep = ""), no = "")
        )), collapse = ", "), ")", sep = ""), no = ""),
        "<br>",
        ## fragments
        resultObj$numberOfMS2PeaksAboveThreshold, " / ", resultObj$numberOfMS2PeaksOriginal, " fragments were successfully imported.", 
        ifelse(test = resultObj$numberOfMS2PeaksAboveThreshold < resultObj$numberOfMS2PeaksOriginal, yes = paste(" (",paste( Filter(nchar, c(
          ifelse(test = resultObj$numberOfTooHeavyFragments      > 0, yes = paste(resultObj$numberOfTooHeavyFragments,      " too heavy",      sep = ""), no = ""), 
          ifelse(test = resultObj$numberOfMS2PeaksBelowThreshold > 0, yes = paste(resultObj$numberOfMS2PeaksBelowThreshold, " low intensity",  sep = ""), no = "")
        )), collapse = ", "), ")", sep = ""), no = ""),
        "<br>",
        ## MS1 features
        resultObj$numberOfPrecursors, " / ", resultObj$numberOfParsedMs1Features, " MS\u00B9 features were successfully imported.",
        ifelse(test = resultObj$numberOfPrecursors < resultObj$numberOfParsedMs1Features, yes = paste(" (",paste( Filter(nchar, c(
          ifelse(test = resultObj$numberOfRemovedPrecursorIsotopePeaks > 0, yes = paste(resultObj$numberOfRemovedPrecursorIsotopePeaks, " were isotopes",   sep = ""), no = ""), 
          ifelse(test = resultObj$numberOfUnmappedPrecursors           > 0, yes = paste(resultObj$numberOfUnmappedPrecursors,           " without spectra", sep = ""), no = ""), 
          ifelse(test = resultObj$numberOfDuplicatedPrecursors         > 0, yes = paste(resultObj$numberOfDuplicatedPrecursors,         " duplicated",      sep = ""), no = "")
        )), collapse = ", "), ")", sep = ""), no = ""),
        sep = ""
      )
      showInfoDialog(msg)
      
      "
      125 / 509 spectra were imported successfully. (380 low intensity, 3 too heavy)
      99 / 125 spectra were successfully mapped to MS¹ features. (26 unmapped)
      774 / 18218 fragments were successfully imported. (2609 too heavy, 14835 low intensity)
      99 / 502 MS¹ features were successfully imported. (1 were isotopes, 402 without spectra)
      "
      
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
    obsUpdateProjectDescription <- observeEvent(input$updateProjectDescription, {
      updateProjectDescription <- as.numeric(input$updateProjectDescription)
      
      print(paste("Observe updateProjectDescription", updateProjectDescription))
      
      #################################################
      ## check if button was hit
      if(updateProjectDescription == updateProjectDescriptionButtonValue)
        return()
      updateProjectDescriptionButtonValue <<- updateProjectDescription
      
      projectDescription <- input$projectDescription2
      projectDescription <- gsub(";", "_", gsub(",", "_", gsub("\t", "_", projectDescription)))
      dataList$importParameterSet$projectDescription <<- projectDescription
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
      updateCheckboxInput(session = session, inputId = "showPlotControls",      value    = state$showPlotControls)
      updateCheckboxInput(session = session, inputId = "showClusterLabels",     value    = state$showClusterLabels)
      updateRadioButtons( session = session, inputId = "heatmapContent",        selected = state$heatmapContent)
      updateRadioButtons( session = session, inputId = "heatmapOrdering",       selected = state$heatmapOrdering)
      updateRadioButtons( session = session, inputId = "hcaPrecursorLabels",    selected = state$hcaPrecursorLabels)
      updateCheckboxInput(session = session, inputId = "showScoresLabels",      value    = state$showScoresLabels)
      updateRadioButtons( session = session, inputId = "loadingsLabels",        selected = state$loadingsLabels)
      updateCheckboxGroupInput(session = session, inputId = "showLoadingsFeatures", selected = c(
        ifelse(test = state$showLoadingsFeaturesAnnotated,   yes = "Annotated",     no = NULL),
        ifelse(test = state$showLoadingsFeaturesUnannotated, yes = "Not Annotated", no = NULL),
        ifelse(test = state$showLoadingsFeaturesSelected,    yes = "Selected",      no = NULL),
        ifelse(test = state$showLoadingsFeaturesUnselected,  yes = "Not Selected",  no = NULL)
      ))
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesAnnotated",    selected = state$showLoadingsFeaturesAnnotated)
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesUnannotated",  selected = state$showLoadingsFeaturesUnannotated)
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesSelected",     selected = state$showLoadingsFeaturesSelected)
      #updateRadioButtons( session = session, inputId = "showLoadingsFeaturesUnselected",   selected = state$showLoadingsFeaturesUnselected)
      updateCheckboxInput(session = session, inputId = "showLoadingsAbundance", value    = state$showLoadingsAbundance)
    })
    obsShowPlotControls <- observeEvent(input$showPlotControls, {
      showPlotControls <- input$showPlotControls
      print(paste("Observe showPlotControls", showPlotControls))
      state$showPlotControls <<- showPlotControls
    })
    obsShowClusterLabels <- observeEvent(input$showClusterLabels, {
      showClusterLabels <- input$showClusterLabels
      print(paste("Observe showClusterLabels", showClusterLabels))
      state$showClusterLabels <<- showClusterLabels
      drawDendrogramPlot(consoleInfo = "showClusterLabels")
    })
    obsHeatmapContent <- observeEvent(input$heatmapContent, {
      heatmapContent <- input$heatmapContent
      
      if(is.null(dataList))
        return()
      
      print(paste("Observe heatmapContent", heatmapContent))
      state$heatmapContent <<- heatmapContent
      
      drawDendrogramPlotImpl()
      #drawHeatmapPlotImpl(consoleInfo = "heatmapContent")
      
      ## TODO does not work because "### GUI ### Generate right GUI" resets inputs
      #if(heatmapContent == "Log-fold-change") {
      #  shinyjs::disable("heatmapOrdering")
      #} else {
      #  shinyjs::enable("heatmapOrdering")
      #}
    })
    obsHeatmapOrdering <- observeEvent(input$heatmapOrdering, {
      if(is.null(dataList))
        return()
      
      heatmapOrdering <- input$heatmapOrdering
      print(paste("Observe heatmapOrdering", heatmapOrdering))
      state$heatmapOrdering <<- heatmapOrdering
      
      #drawHeatmapPlotImpl(consoleInfo = "heatmapOrdering")
    })
    obsHcaPrecursorLabels <- observeEvent(input$hcaPrecursorLabels, {
      hcaPrecursorLabels <- input$hcaPrecursorLabels
      print(paste("Observe hcaPrecursorLabels", hcaPrecursorLabels))
      state$hcaPrecursorLabels <<- hcaPrecursorLabels
      drawDendrogramPlot(consoleInfo = "hcaPrecursorLabels")
    })
    obsShowScoresLabels <- observeEvent(input$showScoresLabels, {
      showScoresLabels <- input$showScoresLabels
      print(paste("Observe showScoresLabels", showScoresLabels))
      state$showScoresLabels <<- showScoresLabels
      drawPcaScoresPlot(consoleInfo = "showScoresLabels")
    })
    obsLoadingsLabels <- observeEvent(input$loadingsLabels, {
      loadingsLabels <- input$loadingsLabels
      print(paste("Observe loadingsLabels", loadingsLabels))
      state$loadingsLabels <<- loadingsLabels
      drawPcaLoadingsPlot(consoleInfo = "loadingsLabels")
    })
    
    obsShowLoadingsFeatures <- observeEvent(input$showLoadingsFeatures, {
      showLoadingsFeatures <- input$showLoadingsFeatures
      print(paste("Observe showLoadingsFeatures", paste(showLoadingsFeatures, collapse = ";")))
      
      {
      state$showLoadingsFeaturesAnnotated   <<- "Annotated"     %in% showLoadingsFeatures
      state$showLoadingsFeaturesUnannotated <<- "Not Annotated" %in% showLoadingsFeatures
      state$showLoadingsFeaturesSelected    <<- "Selected"      %in% showLoadingsFeatures
      state$showLoadingsFeaturesUnselected  <<- "Not Selected"  %in% showLoadingsFeatures
      }
      
      #drawPcaLoadingsPlot(consoleInfo = "showLoadingsFeatures")
    })
    #obsShowLoadingsFeaturesAnnotated <- observeEvent(input$showLoadingsFeaturesAnnotated, {
    #  showLoadingsFeaturesAnnotated <- input$showLoadingsFeaturesAnnotated
    #  print(paste("Observe showLoadingsFeaturesAnnotated", showLoadingsFeaturesAnnotated))
    #  state$showLoadingsFeaturesAnnotated <<- showLoadingsFeaturesAnnotated
    #  drawPcaLoadingsPlot(consoleInfo = "showLoadingsFeaturesAnnotated")
    #})
    #obsShowLoadingsFeaturesUnannotated <- observeEvent(input$showLoadingsFeaturesUnannotated, {
    #  showLoadingsFeaturesUnannotated <- input$showLoadingsFeaturesUnannotated
    #  print(paste("Observe showLoadingsFeaturesUnannotated", showLoadingsFeaturesUnannotated))
    #  state$showLoadingsFeaturesUnannotated <<- showLoadingsFeaturesUnannotated
    #  drawPcaLoadingsPlot(consoleInfo = "showLoadingsFeaturesUnannotated")
    #})
    #obsShowLoadingsFeaturesSelected <- observeEvent(input$showLoadingsFeaturesSelected, {
    #  showLoadingsFeaturesSelected <- input$showLoadingsFeaturesSelected
    #  print(paste("Observe showLoadingsFeaturesSelected", showLoadingsFeaturesSelected))
    #  state$showLoadingsFeaturesSelected <<- showLoadingsFeaturesSelected
    #  drawPcaLoadingsPlot(consoleInfo = "showLoadingsFeaturesSelected")
    #})
    #obsShowLoadingsFeaturesUnselected <- observeEvent(input$showLoadingsFeaturesUnselected, {
    #  showLoadingsFeaturesUnselected <- input$showLoadingsFeaturesUnselected
    #  print(paste("Observe showLoadingsFeaturesUnselected", showLoadingsFeaturesUnselected))
    #  state$showLoadingsFeaturesUnselected <<- showLoadingsFeaturesUnselected
    #  drawPcaLoadingsPlot(consoleInfo = "showLoadingsFeaturesUnselected")
    #})
    
    obsShowLoadingsAbundance <- observeEvent(input$showLoadingsAbundance, {
      showLoadingsAbundance <- input$showLoadingsAbundance
      print(paste("Observe showLoadingsAbundance", showLoadingsAbundance))
      state$showLoadingsAbundance <<- showLoadingsAbundance
      drawPcaLoadingsPlot(consoleInfo = "showLoadingsAbundance")
    })
    ## listen to filter events
    obsClearSearch <- observeEvent(input$clearSearch, {
      clearSearch <- as.numeric(input$clearSearch)
      
      print(paste("Observe clearSearch", clearSearch))
      
      #################################################
      ## check if button was hit
      if(clearSearch == clearSearchButtonValue)
        return()
      clearSearchButtonValue <<- clearSearch
      
      filterSearch <<- NULL
      state$filterSearchActive <<- FALSE
      state$searchfilterValid <<- TRUE
      
      selectionBySearchReset()
      updateSearchInformation()
      
      #################################################
      ## update plots
      if(state$showHCAplotPanel)  drawDendrogramPlot( consoleInfo = "clear search", withHeatmap = TRUE)
      if(state$showPCAplotPanel)  drawPcaPlots(       consoleInfo = "clear search")
    })
    obsApplySearch <- observeEvent(input$applySearch, {
      applySearch <- as.numeric(input$applySearch)
      
      print(paste("Observe applySearch", applySearch))
      
      #################################################
      ## check if button was hit
      if(applySearch == applySearchButtonValue)
        return()
      applySearchButtonValue <<- applySearch
      
      #################################################
      ## MS1 or MS2?
      searchMode <- input$searchMS1orMS2
      if(searchMode == 'MS1 feature m/z'){
        #################################################
        ## get inputs
        filter_ms1_masses <- input$searchMS1mass
        filter_ms1_ppm    <- input$searchMS1massPpm
        
        if(nchar(trimws(filter_ms1_masses)) == 0)
          return()
        
        filter_ms2_masses1  <- NULL
        filter_ms2_masses2  <- NULL
        filter_ms2_masses3  <- NULL
        filter_ms2_ppm      <- NULL
      }
      if(searchMode == 'Fragment m/z'){
        #################################################
        ## get inputs
        filter_ms2_masses1 <- input$search_ms2_masses1
        filter_ms2_masses2 <- input$search_ms2_masses2
        filter_ms2_masses3 <- input$search_ms2_masses3
        filter_ms2_ppm     <- input$searchMS2massPpm
        
        filter_ms1_masses <- NULL
        filter_ms1_ppm    <- NULL
      }
      
      filter_lfc      <- NULL
      filter_average  <- NULL
      groupSet        <- dataList$groups
      includeIgnoredPrecursors  <- input$searchIncludeIgnoredPrecursors
      
      #################################################
      ## do filtering
      sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filterBySamples <- TRUE
      
      print(paste("Observe applySearch", "1m", filter_ms1_masses, "1p", filter_ms1_ppm, "i", includeIgnoredPrecursors, "gs", paste(groupSet, collapse = "-")))
      resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      processSearchFilterResult(resultObj)
    })
    obsApplyGlobalMS2filters <- observeEvent(input$applyGlobalMS2filters, {
      applyGlobalMS2filters <- as.numeric(input$applyGlobalMS2filters)
      
      print(paste("Observe applyGlobalMS2filters", applyGlobalMS2filters))
      
      #################################################
      ## check if button was hit
      if(applyGlobalMS2filters == applyGlobalMS2filtersButtonValue)
        return()
      applyGlobalMS2filtersButtonValue <<- applyGlobalMS2filters
      
      #################################################
      ## get inputs
      filter_ms2_masses1  <- input$globalFilter_ms2_masses1
      filter_ms2_masses2  <- input$globalFilter_ms2_masses2
      filter_ms2_masses3  <- input$globalFilter_ms2_masses3
      filter_ms2_ppm      <- input$globalFilter_ms2_ppm
      
      applyGlobalMS2filters(filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm)
    })
    obsClearGlobalMS2filters <- observeEvent(input$clearGlobalMS2filters, {
      clearGlobalMS2filters <- as.numeric(input$clearGlobalMS2filters)
      
      print(paste("Observe clearGlobalMS2filters", clearGlobalMS2filters))
      
      #################################################
      ## check if button was hit
      if(clearGlobalMS2filters == clearGlobalMS2filtersButtonValue)
        return()
      clearGlobalMS2filtersButtonValue <<- clearGlobalMS2filters
      
      #################################################
      ## get inputs
      filter_ms2_masses1  <- ""
      filter_ms2_masses2  <- ""
      filter_ms2_masses3  <- ""
      filter_ms2_ppm      <- 20
      
      applyGlobalMS2filters(filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm)
    })
    applyGlobalMS2filters <- function(filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm){
      groupSet        <- dataList$groups
      filter_average  <- NULL
      filter_lfc      <- NULL
      includeIgnoredPrecursors  <- TRUE
      filter_ms1_masses <- NULL
      filter_ms1_ppm <- NULL
      
      #################################################
      ## do filtering
      sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filterBySamples <- TRUE
      resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        filterGlobal <<- NULL
        state$globalMS2filterValid <<- FALSE
      } else {
        filterGlobal <<- resultObj$filter
        state$globalMS2filterValid <<- TRUE
      }
      updateGlobalMS2filterInformation()
    }
    obsApplyHcaFilters <- observeEvent(input$applyHcaFilters, {
      applyHcaFilters <- as.numeric(input$applyHcaFilters)
      
      print(paste("Observe applyHcaFilters", applyHcaFilters))
      
      #################################################
      ## check if button was hit
      if(applyHcaFilters == applyHcaFiltersButtonValue)
        return()
      applyHcaFiltersButtonValue <<- applyHcaFilters
      
      #################################################
      ## get inputs
      groupOne        <- input$hcaFilterGroupOne
      groupTwo        <- input$hcaFilterGroupTwo
      filter_average  <- input$hcaFilter_average
      filter_lfc      <- input$hcaFilter_lfc
      includeIgnoredPrecursors  <- input$hcaFilterIncludeIgnoredPrecursors
      
      applyHcaFilters(groupOne, groupTwo, filter_average, filter_lfc, includeIgnoredPrecursors)
    })
    obsClearHcaFilters <- observeEvent(input$clearHcaFilters, {
      clearHcaFilters <- as.numeric(input$clearHcaFilters)
      
      print(paste("Observe clearHcaFilters", clearHcaFilters))
      
      #################################################
      ## check if button was hit
      if(clearHcaFilters == clearHcaFiltersButtonValue)
        return()
      clearHcaFiltersButtonValue <<- clearHcaFilters
      
      #################################################
      ## get inputs
      groupOne        <- input$hcaFilterGroupOne
      groupTwo        <- input$hcaFilterGroupTwo
      filter_average  <- ""
      filter_lfc      <- ""
      includeIgnoredPrecursors  <- TRUE
      
      applyHcaFilters(groupOne, groupTwo, filter_average, filter_lfc, includeIgnoredPrecursors)
    })
    applyHcaFilters <- function(groupOne, groupTwo, filter_average, filter_lfc, includeIgnoredPrecursors){
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      groupSet        <- c(groupOne, groupTwo)
      sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filterBySamples <- TRUE
      #################################################
      ## do filtering and update
      resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        shinyjs::disable("drawHCAplots")
        #disableActionButton(session, "drawHCAplots")
        filterHca <<- NULL
        updateHcaFilterInformation()
        state$hcaFilterValid <<- FALSE
        return()
      }
      
      #################################################
      ## check filter validity
      filterHca <<- resultObj$filter
      updateHcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
      checkHcaFilterValidity(numberOfPrecursorsFiltered)
    }
    obsApplyPcaFilters <- observeEvent(input$applyPcaFilters, {
      applyPcaFilters <- as.numeric(input$applyPcaFilters)
      
      print(paste("Observe applyPcaFilters", applyPcaFilters))
      
      #################################################
      ## check if button was hit
      if(applyPcaFilters == applyPcaFiltersButtonValue)
        return()
      applyPcaFiltersButtonValue <<- applyPcaFilters
      
      #################################################
      ## get inputs
      groupSet        <- input$pcaGroups
      sampleSet       <- input$pcaSamples
      filterBySamples <- input$filterByPCAgroupSamples
      filter_average  <- input$pcaFilter_average
      filter_lfc      <- input$pcaFilter_lfc
      includeIgnoredPrecursors  <- input$pcaFilterIncludeIgnoredPrecursors
      
      if(filterBySamples){
        ## update groups and samples mutually
        
        ## groups which are covered by at least one sample
        groupsFromSamples <- unlist(lapply(X = dataList$groups, FUN = function(x){
          samplesOfGroups <- dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
          if(any(samplesOfGroups %in% sampleSet))
            return(x)
          else
            return(NULL)
        }))
        
        ## samples which ae covered by a group
        samplesFromGroups <- dataList$dataColumnsNameFunctionFromGroupNames(groups = groupSet, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
        groupSet  <- intersect(groupSet, groupsFromSamples)
        sampleSet <- intersect(sampleSet, samplesFromGroups)
      } else {
        sampleSet <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      }
      
      if(length(groupSet) != 2)
        filter_lfc <- NULL
      
      applyPcaFilters(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, includeIgnoredPrecursors)
    })
    obsClearPcaFilters <- observeEvent(input$clearPcaFilters, {
      clearPcaFilters <- as.numeric(input$clearPcaFilters)
      
      print(paste("Observe clearPcaFilters", clearPcaFilters))
      
      #################################################
      ## check if button was hit
      if(clearPcaFilters == clearPcaFiltersButtonValue)
        return()
      clearPcaFiltersButtonValue <<- clearPcaFilters
      
      #################################################
      ## get inputs
      groupSet        <- dataList$groups
      sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filterByPCAgroupSamples <- TRUE
      filter_average  <- ""
      filter_lfc      <- ""
      includeIgnoredPrecursors  <- TRUE
      
      if(length(groupSet) != 2)
        filter_lfc <- NULL
      
      applyPcaFilters(groupSet, sampleSet, filterByPCAgroupSamples, filter_average, filter_lfc, includeIgnoredPrecursors)
    })
    applyPcaFilters <- function(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, includeIgnoredPrecursors){
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      
      #################################################
      ## do filtering and update
      resultObj <- doPerformFiltering(groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors)
      
      if(resultObj$error){
        shinyjs::disable("drawPCAplots")
        #disableActionButton(session, "drawPCAplots")
        filterPca <<- NULL
        updatePcaFilterInformation()
        state$pcaFilterValid <<- FALSE
        return()
      }
      
      filterPca <<- resultObj$filter
      updatePcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterPca$numberOfPrecursorsFiltered
      checkPcaFilterValidity(numberOfPrecursorsFiltered)
    }
    observeGroupSet <- observeEvent(input$pcaGroups, {
      print(paste("observe groups change", paste(input$pcaGroups, collapse = "-"), length(input$pcaGroups), length(input$pcaGroups) == 2))
      shinyjs::toggleState("pcaFilter_lfc", length(input$pcaGroups) == 2)
      
      if(FALSE){
      ## update samples
      sampleNames <- dataList$dataColumnsNameFunctionFromGroupNames(groups = input$pcaGroups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
      updateCheckboxGroupInput(session = session, inputId = "pcaSamples",                         selected = sampleNames)
      }
    })
    observeSampleSet <- observeEvent(input$pcaSamples, {
      print(paste("observe samples change", paste(input$pcaSamples, collapse = "-"), length(input$pcaSamples)))
      
      groupsFromSamples <- unlist(lapply(X = dataList$groups, FUN = function(x){
        samplesOfGroups <- dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
        if(any(samplesOfGroups %in% input$pcaSamples))
          return(x)
        else
          return(NULL)
      }))
      
      shinyjs::toggleState("pcaFilter_lfc", length(groupsFromSamples) == 2)
      
      #########################################################################################
      ## update filter
      
      ## TODO 888
      
      if(FALSE){
      filter <- doPerformFiltering(dataList$groups, NULL, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      if(length(dataList$groups) == 1)
        filter2 <- doPerformFiltering(c(dataList$groups[[1]], dataList$groups[[1]]), NULL, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, TRUE)$filter
      else
        filter2 <- filter
      
      filterHca    <<- filter2
      filterPca    <<- filter
      
      updateHcaFilterInformation()
      updatePcaFilterInformation()
      
      state$hcaFilterValid <<- TRUE
      state$pcaFilterValid <<- TRUE
      
      checkHcaFilterValidity(filter2$numberOfPrecursorsFiltered)
      checkPcaFilterValidity(filter$numberOfPrecursorsFiltered)
      }
      
      if(!is.null(filterHca)){
        applyHcaFilters(
          filterHca$groupSetOriginal[[1]], 
          filterHca$groupSetOriginal[[2]], 
          filterHca$filter_averageOriginal, 
          filterHca$filter_lfcOriginal, 
          filterHca$includeIgnoredPrecursorsOriginal
        )
      }
      if(!is.null(filterPca)){
        applyPcaFilters(
          filterPca$groupSetOriginal, 
          filterPca$sampleSetOriginal, 
          filterPca$filterBySamplesOriginal, 
          filterPca$filter_averageOriginal, 
          filterPca$filter_lfcOriginal, 
          filterPca$includeIgnoredPrecursorsOriginal
        )
      }
    })
    
    observeSelectAllPCAGroups <- observeEvent(input$selectAllPCAGroups, {
      print(paste("observe selectAllPCAGroups"))
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = dataList$groups)
    })
    observeSelectNoPCAGroups <- observeEvent(input$selectNoPCAGroups, {
      print(paste("observe selectNoPCAGroups"))
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = NULL)
    })
    observeSelectInvertedPCAGroups <- observeEvent(input$selectInvertedPCAGroups, {
      print(paste("observe selectInvertedPCAGroups"))
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = setdiff(dataList$groups, input$pcaGroups))
    })
    
    ## listen to draw button events
    obsDrawHCA <- observeEvent(input$drawHCAplots, {
      #################################################
      ## get input
      drawPlots <- as.numeric(input$drawHCAplots)
      
      print(paste("Observe draw HCA plots", drawPlots))
      
      ## check if button was hit
      if(drawPlots == drawHCAButtonValue)
        return()
      drawHCAButtonValue <<- drawPlots
      
      distanceMeasure <- input$hcaDistanceFunction
      #clusterMethod <- input$hcaClusterMethod
      clusterMethod <- "ward.D"
      print(paste("Observe draw HCA plots", "D", distanceMeasure, "M", clusterMethod))
      
      ##########################
      ## calc
      
      ## compute distance matrix
      withProgress(message = 'Calculating distances...', value = 0, {
        currentDistanceMatrixObj <<- calculateDistanceMatrix(dataList = dataList, filter = filterHca$filter, distanceMeasure = distanceMeasure, progress = TRUE)
      })
      ## compute cluster
      withProgress(message = 'Calculating cluster...', value = 0, {
        clusterDataList <<- calculateCluster(progress = TRUE, dataList = dataList, filterObj = filterHca, distanceMatrix = currentDistanceMatrixObj$distanceMatrix, method = clusterMethod, distanceMeasure = distanceMeasure)
      })
      
      ##########################
      ## hca selections
      if(!is.null(listForTable_Fragment_PCA)){ ## selection from fragment
        fragmentIndex <- which(dataList$fragmentMasses == fragmentsX[[selectionFragmentSelectedFragmentIndex]])
        precursorSet  <- which(dataList$featureMatrix[, fragmentIndex] != 0)
        selectionByFragmentInitHca(precursorSet)
      } else {
        selectionByFragmentReset()
      }
      if(!is.null(listForTable_Analysis_PCA)){ ## selection from PCA
        precursorSet <- filterPca$filter[selectionAnalysisPcaLoadingSet]
        selectionByAnalysisInitHca(precursorSet)
      } else {
        selectionByAnalysisReset()
      }
      if(!is.null(filterSearch)){ ## selection from search
        selectionBySearchInitHca(filterSearch$filter)
      } else {
        selectionBySearchReset()
      }
      
      ##########################
      ## reset MS2 stuff
      if(!state$showPCAplotPanel){
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
      }
      
      ##########################
      ## draw
      resetHcaPlotRange()
      drawDendrogramPlot(consoleInfo = "init output$plotDendrogram", withHeatmap = TRUE)
      drawMS2Plot(consoleInfo = "init output$plotMS2")
      drawAnnotationLegendHCA(consoleInfo = "init output$plotAnnoLegend")
      
      if(!state$anyPlotDrawn){
        drawMS2Legend(consoleInfo = "init output$ms2LegendPlot")
        drawFragmentDiscriminativityLegend(consoleInfo = "init output$plotFragmentDiscriminativityLegend")
        drawHeatmapLegend(consoleInfo = "init output$plotHeatmapLegend")
        drawDendrogramLegend(consoleInfo = "init output$calcPlotDendrogramLegend")
        state$anyPlotDrawn <<- TRUE
      }
      
      ## state
      state$showHCAplotPanel <<- TRUE
      state$plotHcaShown <<- TRUE
      updateChangePlotRadioButton()
      
      ##########################
      ## update info and tip
      output$information <- renderText({
        print(paste("update output$information clear"))
        paste("", sep = "")
      })
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste(
          "Hover or select a cluster node or leaf node to view information about the corresponding MS\u00B9 feature cluster or MS\u00B9 feature respectively.", 
          "Brush horizontally and double-click to zoom in.", 
          "Double-click to zoom out.", 
          sep = "\n"
        )
      })
    })
    obsDrawPCA <- observeEvent(input$drawPCAplots, {
      #################################################
      ## get input
      drawPlots <- as.numeric(input$drawPCAplots)
      
      print(paste("Observe draw PCA plots", drawPlots))
      
      ## check if button was hit
      if(drawPlots == drawPCAButtonValue)
        return()
      drawPCAButtonValue <<- drawPlots
      
      ms1AnalysisMethod <- input$ms1AnalysisMethod
      pcaScaling      <- input$pcaScaling
      pcaLogTransform <- input$pcaLogTransform
      pcaDimensionOne <<- as.numeric(input$pcaDimensionOne)
      pcaDimensionTwo <<- as.numeric(input$pcaDimensionTwo)
      
      #################################################
      ## calc PCA
      pca <- calculatePCA(dataList = dataList, filterObj = filterPca, ms1AnalysisMethod = ms1AnalysisMethod, scaling = pcaScaling, logTransform = pcaLogTransform)
      
      pcaDataList <<- list()
      pcaDataList$pcaObj <<- pca
      pcaDataList$pcaScoresX <<- pca$scores[, pcaDimensionOne]
      pcaDataList$pcaScoresY <<- pca$scores[, pcaDimensionTwo]
      pcaDataList$pcaLoadingsX <<- pca$loadings[, pcaDimensionOne]
      pcaDataList$pcaLoadingsY <<- pca$loadings[, pcaDimensionTwo]
      pcaDataList$dimensionOne <<- pcaDimensionOne
      pcaDataList$dimensionTwo <<- pcaDimensionTwo
      
      ##########################
      ## hca selections
      if(!is.null(listForTable_Fragment_HCA)){ ## selection from fragment
        fragmentIndex <- which(dataList$fragmentMasses == fragmentsX[[selectionFragmentSelectedFragmentIndex]])
        precursorSet  <- which(dataList$featureMatrix[, fragmentIndex] != 0)
        selectionByFragmentInitPca(precursorSet)
      } else {
        selectionByFragmentReset()
      }
      if(!is.null(listForTable_Analysis_HCA)){ ## selection from HCA
        precursorSet <- getPrecursorSetFromTreeSelections(clusterDataList = clusterDataList, clusterLabels = selectionAnalysisTreeNodeSet)
        selectionByAnalysisInitPca(precursorSet)
      } else {
        selectionByAnalysisReset()
      }
      if(!is.null(filterSearch)){ ## selection from search
        selectionBySearchInitPca(filterSearch$filter)
      } else {
        selectionBySearchReset()
      }
      
      ##########################
      ## reset MS2 stuff
      if(!state$showHCAplotPanel){
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
      }
      
      ##########################
      ## draw
      resetPcaPlotRange()
      drawPcaPlots(consoleInfo = "drawPCA output$plotPcaScores")
      drawMS2Plot(consoleInfo = "drawPCA output$plotMS2")
      drawAnnotationLegendPCA(consoleInfo = "init output$plotAnnoLegend")
      
      scoresGroups <<- list(
        groups = filterPca$groups,
        colors = colorPaletteScores()[unlist(lapply(X = filterPca$groups, FUN = dataList$groupIdxFromGroupName))]
      )
      state$scoresGroupsLegendHeight <<- scoresGroupsLegendEntryHeight * (length(scoresGroups$groups) + 1)
      drawScoresGroupsLegend(consoleInfo = "init output$plotScoresGroupsLegend")
      
      if(!state$anyPlotDrawn){
        drawMS2Legend(consoleInfo = "init output$ms2LegendPlot")
        drawFragmentDiscriminativityLegend(consoleInfo = "init output$plotFragmentDiscriminativityLegend")
        drawHeatmapLegend(consoleInfo = "init output$plotHeatmapLegend")
        drawDendrogramLegend(consoleInfo = "init output$calcPlotDendrogramLegend")
        state$anyPlotDrawn <<- TRUE
      }
      
      ## state
      state$showPCAplotPanel <<- TRUE
      state$plotPcaShown <<- TRUE
      updateChangePlotRadioButton()
      
      ##########################
      ## update info and tip
      output$information <- renderText({
        print(paste("update output$information clear"))
        paste("", sep = "")
      })
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Hover a score node in the scores plot or a loadings node in the loadings plot to view information about the corresponding sample or MS\u00B9 feature respectively.", "Brush and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    ## listen to dendrogram/heatmap plot mouse events
    
    obsClassifierSelectionTable_rows_selected <- observeEvent(input$classifierSelectionTable_rows_selected, {
      print(paste("Observe classifierSelectionTable_rows_selected", input$classifierSelectionTable_rows_selected))
      if(is.null(input$classifierSelectionTable_rows_selected)){
        session$sendCustomMessage("disableButton", "doAnnotation")
      } else {
        session$sendCustomMessage("enableButton", "doAnnotation")
      }
    })
    obsDoAnnotation <- observeEvent(input$doAnnotation, {
      doAnnotation <- as.numeric(input$doAnnotation)
      print(paste("Observe doAnnotation", doAnnotation))
      
      #################################################
      ## check if button was hit
      #if(doAnnotation == doAnnotationButtonValue)
      #  return()
      #doAnnotationButtonValue <<- doAnnotation
      
      session$sendCustomMessage("disableButton", "doAnnotation")
      
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
          bsModal(id = "noAnnotationsPopupDialog", title = "No annotations detected", trigger = "", size = "large",
                  HTML("No MS\u00B9 feature could be annotated with one or more metabolite families.")
          )
        })
        toggleModal(session = session, modalId = "noAnnotationsPopupDialog", toggle = "open")
        
        return()
      }
      
      bestPValues <- format(bestPValues, digits=4)
      numPValues  <- sapply(X = classToSpectra_class, FUN = length)
      classes     <- names(classToSpectra_class)
      
      ## shrink class names
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
      
      resultTable_classes <- data.frame(
        "Class"       = classesShort, 
        "Hits"        = numPValues,
        "Best pValue" = bestPValues#,
        #row.names = rep("")
      )
      
      output$familyCount <- renderText({
        print(paste("init output$familyCount"))
        paste("Number of putative metabolite families:", nrow(resultTable_classes))
      })
      
      #selected <- ifelse(test = nrow(resultTable_classes) > 0, yes = 1, no = NULL)
      output$annotationResultTableClass <- DT::renderDataTable(
        expr = resultTable_classes,
        server = FALSE, escape = FALSE, rownames = FALSE,
        selection = list(mode = "single"),#, selected = selected), 
        options = list(
          #scrollY = "600px",
          scrollY = "30vh",
          preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
          drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
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
    obsAnnotationResultTableClass_selection <- observeEvent(input$annotationResultTableClass_rows_selected, {
      selectedRowIdx <- input$annotationResultTableClass_rows_selected
      print(paste("Selected class row:", selectedRowIdx))
      selectedClassRowIdx <<- selectedRowIdx
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
        paste("Number of spectrum hits:", ifelse(test = length(precursorIndecesAll) > maximumNumberOfDisplayedSpectrumHits, yes = paste(length(precursorIndeces), "/", length(precursorIndecesAll), sep = ""), no = length(precursorIndecesAll)))
      })
      
      
      returnObj <- createPlotOutput(
        dataList = dataList,
        id = "plotSpectrumVsClass", 
        frequentFragments = frequentFragments, 
        characteristicFragments = characteristicFragments, 
        precursorIndeces = precursorIndeces, 
        tableCounter = ms1FeatureVsClassTableCounter
      )
      numberOfMatchingMasses_i <- returnObj$numberOfMatchingMasses_i
      inputs_i <- returnObj$inputs_i
      
      checkboxes   <- createCheckboxInputFields2(
        FUN = checkboxInput, 
        id = "MS1_feature_confirm", 
        values = rep(x = FALSE, times = length(precursorIndeces)), 
        tableCounter = ms1FeatureVsClassTableCounter
      )
      
      resultTable_features <- data.frame(check.names = FALSE,
        "Mz/RT"    = precursorLabels,
        "pValue"   = pValues,
        "Hits"     = numberOfMatchingMasses_i,
        "Confirm"  = checkboxes,
        "Spectrum" = inputs_i
      )
      selected <- ifelse(test = nrow(resultTable_features) > 0, yes = 1, no = NULL)
      output$annotationResultTableFeature <- DT::renderDataTable(
        expr = resultTable_features,
        server = FALSE, escape = FALSE, rownames = FALSE,
        selection = list(mode = "single"),#, selected = selected), 
        options = list(
          #scrollY = "600px",
          scrollY = "30vh",
          preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
          drawCallback    = JS('function() { Shiny.bindAll(  this.api().table().node()); }'),
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
    obsAnnotationResultTableFeature_selection <- observeEvent(input$annotationResultTableFeature_rows_selected, {
      selectedRowIdx <- input$annotationResultTableFeature_rows_selected
      print(paste("Selected feature row:", selectedRowIdx))
      selectedClassFeatureRowIdx <<- selectedRowIdx
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
        precursorIndex = precursorIndex
      )
    }
    drawMS2vsClassPlot <- function(consoleInfo = NULL, frequentFragments, characteristicFragments, precursorIndex){
      output$plotMS2vsClass <- renderPlot({
        print(paste("### SvC ###", consoleInfo))
        drawMS2vsClassPlotImpl(
          dataList=dataList,
          frequentFragments=frequentFragments, 
          characteristicFragments=characteristicFragments, 
          precursorIndex=precursorIndex
        )
      })
    }
    drawMS2vsClassPlotImpl <- function(dataList, frequentFragments, characteristicFragments, precursorIndex){
      ## class statistics for class plot
      returnObj <- preprocessClassPlot(frequentFragments, characteristicFragments)
      masses_class    <- returnObj$masses_class
      frequency_class <- returnObj$frequency_class
      colors_class    <- returnObj$colors_class
      
      ## match spectrum masses for spectrum plot
      returnObj <- preprocessSpectrumVsClassPlot(dataList, precursorIndex, masses_class)
      masses_spec     <- returnObj$masses_spec
      intensity_spec  <- returnObj$intensity_spec
      colors_spec     <- returnObj$colors_spec
      #numberOfMatchingMasses <- returnObj$numberOfMatchingMasses
      
      #xInterval <- c(dataList$minimumMass, dataList$maximumMass)
      
      calcPlotSpectrumVsClass_big(
        masses_spec=masses_spec, 
        intensity_spec=intensity_spec, 
        colors_spec=colors_spec, 
        masses_class=masses_class, 
        frequency_class=frequency_class, 
        colors_class=colors_class, 
        xInterval = specVsClassRange$xInterval
      )
    }
    
    ## TODO 999
    obsApplyConfirmedAnnotations <- observeEvent(input$confirmAnnotation, {
      confirmAnnotation <- as.numeric(input$confirmAnnotation)
      
      #################################################
      ## check if button was hit
      if(confirmAnnotation == confirmAnnotationButtonValue)
        return()
      confirmAnnotationButtonValue <<- confirmAnnotation
      
      value <- input$newAnnotationValue2
      color <- input$newAnnotationColor2
      
      ## fetch data
      confirm <- getInputValues(id = paste("MS1_feature_confirm", sep = "_"), counter = ms1FeatureVsClassTableCounter, len = length(selectedClassPrecursorIndeces))
      precursorIndeces <- selectedClassPrecursorIndeces[confirm]
      
      print(paste("Annotate: Class ", value, " to ", paste(precursorIndeces, collapse = ", "), sep = ""))
      
      ## apply annotation and update
      addAnnotation(precursorSet = precursorIndeces, annotationValue = value, annotationColor = color)
      
      
      ## update e.g. plots
    })
    obsToggleConfirmAnnoButton <- observeEvent(input$newAnnotationValue2, {
      value <- input$newAnnotationValue2
      
      print(paste("Observe newAnnotationValue2", nchar(value)))
      
      if(nchar(value) > 0)
        shinyjs::enable("confirmAnnotation")
      else
        shinyjs::disable("confirmAnnotation")
    })
    
    #obsDendLabelsHeatmap1 <- observe({
    #  eventdata <- event_data("plotly_click", source = "dendLabelsHeatmap")
    #  print("hihi")
    #  str(eventdata)
    #})
    if(FALSE){
     ## plotly
    obsDendLabelsHeatmap <- observeEvent(event_data("plotly_click", source = "dendLabelsHeatmap"), {
      eventdata <- event_data("plotly_click", source = "dendLabelsHeatmap")
      #str(eventdata)
      ## 'data.frame':	1 obs. of  4 variables:
      ## $ curveNumber: int 7
      ## $ pointNumber: int 6
      ## $ x          : num 8.5
      ## $ y          : num 0.697
      
      if(is.null(eventdata))
        return()
      
      curveName <- curveNumberToCurveName[which(curveNumberToCurveName[, "curveNumber"] == eventdata$curveNumber), "name"]
      #print(paste("curveName:", curveName))
      if(curveName != "nodes")
        return()
      
      #nodeIndex <- eventdata$pointNumber + 1
      nodeIndex <- which((
        #eventdata$x == clusterDataList$poiCoordinatesX & 
        #eventdata$x %in% clusterDataList$poiCoordinatesX
        abs(eventdata$x - clusterDataList$poiCoordinatesX) <= 0.0001
      ) & (
        #eventdata$y == clusterDataList$poiCoordinatesY
        #eventdata$y %in% clusterDataList$poiCoordinatesY
        abs(eventdata$y - clusterDataList$poiCoordinatesY) <= 0.0001
      ))
      nodeLabel <- clusterDataList$poiLabels[[nodeIndex]]
      
      
      ## fetch ms2 spectrum
      resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = nodeLabel)
      
      ## keep fragment selection
      selectionFragmentSelectedFragmentIndexNew <- NULL
      if(!is.null(selectionFragmentSelectedFragmentIndex)){
        fragmentMass <- fragmentsX[[selectionFragmentSelectedFragmentIndex]]
        if(fragmentMass %in% resultObj$fragmentMasses)
          selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
      }
      
      fragmentsX <<- resultObj$fragmentMasses
      fragmentsY <<- resultObj$fragmentAbundances
      fragmentsColor <<- resultObj$fragmentColor
      fragmentsDiscriminativity <<- resultObj$fragmentDiscriminativity
      fragmentsXhovered <<- NULL
      fragmentsYhovered <<- NULL
      fragmentsColorHovered <<- NULL
      fragmentsDiscriminativityHovered <<- NULL
      
      print(paste("fragments:", paste(resultObj$fragmentMasses, collapse = "; ")))
      
      #################################################
      ## output as message
      selectionByHca(nodeLabel)
      
      ## update the selected fragment in case of overlapping spectra
      if(!is.null(selectionFragmentSelectedFragmentIndexNew))
        selectionByFragment(selectionFragmentSelectedFragmentIndexNew)
      else
        selectionByFragmentReset()
      
      output$information <- renderText({
        print(paste("update output$information", resultObj$infoText))
        paste(resultObj$infoText, sep = "")
      })
      
      ## TODO
      #drawMS2PlotImpl
      
      resetMS2PlotRange()
      drawMS2Plot()
      
      #output$plotDendrogram <- renderPlotly({
    })
    }
    
    #output$plotTmp <- renderPlotly({
    #  
    #  # Read in hover data
    #  eventdata <- event_data("plotly_click", source = "dendLabelsHeatmap")
    #  print("haha")
    #  str(eventdata)
    #  validate(need(!is.null(eventdata), "Hover over the time series chart to populate this heatmap"))
    #  
    #  curveName <- curveNumberToCurveName[which(curveNumberToCurveName[, "curveNumber"] == eventdata$curveNumber), "name"]
    #  print(curveName)
    #  
    #  plot_ly(x = 1:10, y = rnorm(10), type = "scatter", mode = "markers")
    #})
    
    obsDendrogramHover <- observeEvent(input$plotDendrogram_hover, {
      hoverX <- input$plotDendrogram_hover$x
      hoverY <- input$plotDendrogram_hover$y
      
      plotWidth  <- session$clientData$output_plotDendrogram_width
      plotHeight <- session$clientData$output_plotDendrogram_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramPlotRangeY
      )
      if(is.null(minimumIndex)){
        print(paste("Observe dendrogram hover", minimumIndex))
        if(!is.null(fragmentsXhovered)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
        }
      } else {
        minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
        print(paste("Observe dendrogram hover i", minimumIndex, "l", minimumLabel))
        resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
        
        ## TODO 8765
        #putativeMetaboliteFamilies <- NULL
        #if(!is.null(classToSpectra_class)){
        #  putativeMetaboliteFamilies <- evaluateDendrogramCluster(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel, classToSpectra_class = classToSpectra_class)
        #}
        
        #################################################
        ## output as message
        output$information <- renderText({
          print(paste("update output$information", resultObj$infoText))
          paste(
            resultObj$infoText, 
            #ifelse(test = is.null(putativeMetaboliteFamilies), yes = "", no = 
            #         paste("\n", paste(putativeMetaboliteFamilies, collapse = "\n"), sep = "")
            #), 
          sep = "")
        })
        
        if(all(!is.null(selectionAnalysisTreeNodeSet), minimumLabel == selectionAnalysisTreeNodeSet)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
          fragmentsColorHovered <<- NULL
          fragmentsDiscriminativityHovered <<- NULL
        } else {
          #################################################
          ## fetch ms2 spectrum
          fragmentsXhovered <<- resultObj$fragmentMasses
          fragmentsYhovered <<- resultObj$fragmentAbundances
          fragmentsColorHovered <<- resultObj$fragmentColor
          fragmentsDiscriminativityHovered <<- resultObj$fragmentDiscriminativity
        }
      }
      
      ## MS2 plot
      drawMS2Plot(consoleInfo = "dendrogram hover output$plotMS2")
      
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Hover a cluster node or leaf node to view information about the corresponding MS\u00B9 feature cluster or MS\u00B9 feature respectively.", "Brush horizontally and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    obsDendrogramClick <- observeEvent(input$plotDendrogram_click, {
      clickX <- input$plotDendrogram_click$x
      clickY <- input$plotDendrogram_click$y
      
      brush <- input$plotDendrogram_brush
      
      plotWidth  <- session$clientData$output_plotDendrogram_width
      plotHeight <- session$clientData$output_plotDendrogram_height
      
      if(is.null(clickX) | is.null(clickY))
        return()
      if(!is.null(brush))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = clusterDataList$poiCoordinatesX, poiCoordinatesY = clusterDataList$poiCoordinatesY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = dendrogramPlotRange$xIntervalSize, plotRangeY = dendrogramPlotRangeY
      )
      print(paste("Observe dendrogram click", is.null(minimumIndex), minimumIndex))
      
      if(is.null(minimumIndex)){
        ## reset stuff
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        selectionByAnalysisReset()
        selectionByFragmentReset()
      } else {
        ## tree selection
        minimumLabel <- clusterDataList$poiLabels[[minimumIndex]]
        #minimumText <- clusterDataList$poiText[[minimumIndex]]
        
        ## fetch ms2 spectrum
        resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel)
        
        ## keep fragment selection
        selectionFragmentSelectedFragmentIndexNew <- NULL
        if(!is.null(selectionFragmentSelectedFragmentIndex)){
          fragmentMass <- fragmentsX[[selectionFragmentSelectedFragmentIndex]]
          if(fragmentMass %in% resultObj$fragmentMasses)
            selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
        }
        
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsColor <<- resultObj$fragmentColor
        fragmentsDiscriminativity <<- resultObj$fragmentDiscriminativity
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        #################################################
        ## output as message
        selectionByHca(minimumLabel)
        
        ## update the selected fragment in case of overlapping spectra
        if(!is.null(selectionFragmentSelectedFragmentIndexNew))
          selectionByFragment(selectionFragmentSelectedFragmentIndexNew)
        else
          selectionByFragmentReset()
        
        ## TODO 8765
        putativeMetaboliteFamilies <- NULL
        if(!is.null(classToSpectra_class)){
          putativeMetaboliteFamilies <- evaluateDendrogramCluster(dataList = dataList, clusterDataList = clusterDataList, treeLabel = minimumLabel, classToSpectra_class = classToSpectra_class)
        }
        
        output$information <- renderText({
          print(paste("update output$information", resultObj$infoText))
          paste(
            resultObj$infoText, 
            ifelse(test = is.null(putativeMetaboliteFamilies), yes = "", no = 
                     paste("\n", paste(putativeMetaboliteFamilies, collapse = "\n"), sep = "")
            ), 
            sep = "")
        })
      }
      
      #################################################
      ## plots
      
      ## cluster dendrogram
      drawDendrogramPlot(consoleInfo = "dendrogram click output$plotDendrogram", withHeatmap = TRUE)
      
      ## MS2 plot
      resetMS2PlotRange()
      drawMS2Plot(consoleInfo = "dendrogram click output$plotMS2")
      
      if(state$showPCAplotPanel)
        ## update PCA plots
        drawPcaPlots(consoleInfo = "dendrogram click output$plotPcaScores")
    })
    obsDendrogramdblClick <- observeEvent(input$plotDendrogram_dblclick, {
      brush <- input$plotDendrogram_brush
      
      print(paste("observe dendrogram dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        dendrogramPlotRange$xMin <<- brush$xmin
        dendrogramPlotRange$xMax <<- brush$xmax
        dendrogramPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        dendrogramPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
      } else {
        ## reset range
        dendrogramPlotRange$xMin <<- 1
        dendrogramPlotRange$xMax <<- filterHca$numberOfPrecursorsFiltered
        dendrogramPlotRange$xInterval <<- c(1, filterHca$numberOfPrecursorsFiltered)
        dendrogramPlotRange$xIntervalSize <<- filterHca$numberOfPrecursorsFiltered - 1
      }
    })
    #obsDendrogramRangeUpdate <- observe({
    #  xInterval <- dendrogramPlotRange$xInterval
    #  
    #  if(!state$showHCAplotPanel)
    #    return()
    #  
    #  print(paste("observe dendrogramPlotRange", xInterval))
    #  
    #  #drawDendrogramPlot(consoleInfo = "update range output$plotDendrogram", withHeatmap = TRUE)
    #})
    obsHeatmaphover <- observeEvent(input$plotHeatmap_hover, {
      hoverX <- input$plotHeatmap_hover$x
      hoverY <- input$plotHeatmap_hover$y
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(hoverX < 0.5 | hoverX > (clusterDataList$numberOfPrecursorsFiltered + 0.5))
        return()
      if(hoverY < 0 | hoverY > 3)
        return()
      
      print(paste("Observe heatmap hover", hoverX, hoverY))
      
      #################################################
      ## info
      treeLeafIndex2 <- as.numeric(format(x = hoverX, digits = 0))
      treeLeafIndex  <- clusterDataList$cluster$order[[treeLeafIndex2]]
      precursorIndex <- filterHca$filter[[treeLeafIndex]]
      
      msg <- list()
      msg[[length(msg) + 1]] <- "MS\u00B9 feature: "
      msg[[length(msg) + 1]] <- dataList$precursorLabels[[precursorIndex]]
      msg[[length(msg) + 1]] <- "\n"
      
      if(hoverY > 2){
        ## lcf
        groupOne <- filterHca$groups[[1]]
        groupTwo <- filterHca$groups[[2]]
        msg[[length(msg) + 1]] <- paste("log-fold-change = log_2( mean(group ", groupOne, ") / mean(group ", groupTwo, ") )", sep = "")
        
        valMeanOne <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupOne)]
        valMeanTwo <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupTwo)]
        msg[[length(msg) + 1]] <- " = log_2( "
        msg[[length(msg) + 1]] <- as.numeric(format(x = valMeanOne, digits = 2))
        msg[[length(msg) + 1]] <- " / "
        msg[[length(msg) + 1]] <- as.numeric(format(x = valMeanTwo, digits = 2))
        msg[[length(msg) + 1]] <- " ) = "
        
        lfc <- dataList$dataFrameMeasurements[precursorIndex, dataList$lfcColumnNameFunctionFromName(groupOne, groupTwo)]
        msg[[length(msg) + 1]] <- as.numeric(format(x = lfc, digits = 2))
      } else {
        if(hoverY > 1){
          ## group 1
          groupHere <- filterHca$groups[[1]]
        } else { ## hoverY <= 1
          ## group 2
          groupHere <- filterHca$groups[[2]]
        }
        msg[[length(msg) + 1]] <- paste("Mean abundance of group ", groupHere, ": ", sep = "")
        valMean <- dataList$dataFrameMeasurements[precursorIndex, dataList$dataMeanColumnNameFunctionFromName(groupHere)]
        msg[[length(msg) + 1]] <- as.numeric(format(x = valMean, digits = 2))
        msg[[length(msg) + 1]] <- " = mean("
        #columnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = groupHere, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
        columnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = groupHere, sampleNamesToExclude = dataList$excludedSamples(groupSampleDataFrame = dataList$groupSampleDataFrame, groups = groupHere))
        vals <- dataList$dataFrameMeasurements[precursorIndex, columnNames]
        vals <- as.numeric(format(x = vals, digits = 2))
        msg[[length(msg) + 1]] <- paste(vals, collapse = ", ")
        msg[[length(msg) + 1]] <- ")"
      }
      
      output$information <- renderText({
        print(paste("update output$information heatmap hover", sep = ""))
        paste(msg, collapse = "")
      })
    })
    ## listen to MS2 plot mouse events
    obsMS2hover <- observeEvent(input$plotMS2_hover, {
      hoverX <- input$plotMS2_hover$x
      hoverY <- input$plotMS2_hover$y
      plotWidth  <- session$clientData$output_plotMS2_width
      plotHeight  <- session$clientData$output_plotMS2_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      if(any(is.null(fragmentsX), length(fragmentsX) == 0))
        return()
      
      ################################################
      ## decide whether the click is close enough to trigger event
      
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = fragmentsX, poiCoordinatesY = fragmentsY, 
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = ms2PlotRange$xIntervalSize, plotRangeY = 1
      )
      print(paste("Observe MS2 hover", minimumIndex))
      if(is.null(minimumIndex)){
        ## nothing in selection range
        #output$information <- renderText({
        #  print(paste("update output$information"))
        #  paste("", sep = "")
        #})
      } else {
        ## point selected
        fragmentIndex <- which(dataList$fragmentMasses == fragmentsX[[minimumIndex]])
        
        if(state$plotHcaShown)
          numberOfPrecursors <- sum(dataList$featureMatrix[filterHca$filter, fragmentIndex] != 0)
        if(state$plotPcaShown)
          numberOfPrecursors <- sum(dataList$featureMatrix[filterPca$filter, fragmentIndex] != 0)
        
        output$information <- renderText({
          print(paste("update output$information"))
          paste(
            "Fragment with m/z = ", fragmentsX[[minimumIndex]], 
            " and (average) abundance = ", format(x = fragmentsY[[minimumIndex]], digits = 0, nsmall = 4), 
            " is present in ", numberOfPrecursors, " MS/MS spectra",
            "\nand has a cluster-discriminating power of ", format(x = fragmentsDiscriminativity[[minimumIndex]]*100, digits = 3, nsmall = 2), "%.", 
            sep = ""
          )
        })
      }
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste(
          "Hover or click a fragment node (only 'Fragments from selection') to view information about this fragment.", 
          "Brush horizontally and double-click to zoom in.", 
          "Double-click to zoom out.", 
          sep = "\n"
        )
      })
    })
    obsMS2click <- observeEvent(input$plotMS2_click, {
      clickX <- input$plotMS2_click$x
      clickY <- input$plotMS2_click$y
      
      brush  <- input$plotMS2_brush
      
      plotWidth  <- session$clientData$output_plotMS2_width
      plotHeight <- session$clientData$output_plotMS2_height
      
      if(is.null(clickX) | is.null(clickY))
        return()
      if(!is.null(brush))
        ## ongoing brushing
        return()
      if(any(is.null(fragmentsX), length(fragmentsX) == 0))
        return()
      
      #################################################
      ## decide whether the click is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = fragmentsX, poiCoordinatesY = fragmentsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = ms2PlotRange$xIntervalSize, plotRangeY = 1
      )
      print(paste("Observe MS2 click", minimumIndex))
      
      if(is.null(minimumIndex)){
        ##########################################
        ## reset click
        selectionByFragmentReset()
      } else {
        ##########################################
        ## peak click
        selectionByFragment(minimumIndex)
      }
      
      ## HCA
      if(state$showHCAplotPanel)
        drawDendrogramPlot(consoleInfo = "MS2 click output$plotDendrogram", withHeatmap = TRUE)
      ## PCA
      if(state$showPCAplotPanel)
        drawPcaLoadingsPlot(consoleInfo = "MS2 click output$plotPcaLoadings")
      
      ## update node selection
      drawMS2Plot(consoleInfo = "MS2 click output$plotMS2")
    })
    obsMS2dblClick <- observeEvent(input$plotMS2_dblclick, {
      brush <- input$plotMS2_brush
      
      if(all(any(is.null(fragmentsX), length(fragmentsX) == 0), any(is.null(fragmentsXhovered), length(fragmentsXhovered) == 0)))
        return()
      
      print(paste("observe MS2 dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        ms2PlotRange$xMin <<- brush$xmin
        ms2PlotRange$xMax <<- brush$xmax
        ms2PlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        ms2PlotRange$xIntervalSize <<- brush$xmax - brush$xmin
      } else {
        ## reset range
        resetMS2PlotRange()
      }
    })
    #obsMS2rangeUpdate <- observe({
    #  xInterval <- ms2PlotRange$xInterval
    #  
    #  if(!(state$showHCAplotPanel | state$showPCAplotPanel))
    #    return()
    #  
    #  print(paste("observe ms2PlotRange", paste(xInterval, collapse = ", ")))
    #  
    #  #drawMS2Plot(consoleInfo = "range update output$plotMS2")
    #})
    ## listen to PCA events
    obsPCAscoresHover <- observeEvent(input$plotPcaScores_hover, {
      hoverX <- input$plotPcaScores_hover$x
      hoverY <- input$plotPcaScores_hover$y
      plotWidth  <- session$clientData$output_plotPcaScores_width
      plotHeight <- session$clientData$output_plotPcaScores_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA scores hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaScoresX, poiCoordinatesY = pcaDataList$pcaScoresY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaScoresPlotRange$xIntervalSize, plotRangeY = pcaScoresPlotRange$yIntervalSize
      )
      print(paste("Observe PCA scores hover", hoverX, hoverY, minimumIndex))
      
      if(is.null(minimumIndex)){
        #output$information <- renderText({
        #  print(paste("update output$information PCA scores hover", sep = ""))
        #  paste("", group, sep = "")
        #})
      }
      else{
        dataColumnName <- filterPca$sampleSet[[minimumIndex]]
        #dataColumnName <- dataList$dataColumnsNameFunctionFromGroupNames(groups = filterPca$groups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))[[minimumIndex]]
        group <- dataList$groupNameFunctionFromDataColumnName(dataColumnName = dataColumnName, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
        output$information <- renderText({
          print(paste("update output$information PCA scores hover", sep = ""))
          paste("Sample '", dataColumnName , "' is a replicate of group '", group, "'.", sep = "")
        })
      }
      
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste("Hover a scores node to view information about the sample.", "Brush and double-click to zoom in.", "Double-click to zoom out.", sep = "\n")
      })
    })
    obsPCAscoresDblClick <- observeEvent(input$plotPcaScores_dblclick, {
      brush <- input$plotPcaScores_brush
      
      print(paste("observe PCAscores dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        pcaScoresPlotRange$xMin <<- brush$xmin
        pcaScoresPlotRange$xMax <<- brush$xmax
        pcaScoresPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        pcaScoresPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
        pcaScoresPlotRange$yMin <<- brush$ymin
        pcaScoresPlotRange$yMax <<- brush$ymax
        pcaScoresPlotRange$yInterval <<- c(brush$ymin, brush$ymax)
        pcaScoresPlotRange$yIntervalSize <<- brush$ymax - brush$ymin
      } else {
        ## reset range
        minX <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
        maxX <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionOne])
        minY <- min(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
        maxY <- max(pcaDataList$pcaObj$scores[, pcaDataList$dimensionTwo])
        
        pcaScoresPlotRange$xMin <<- minX
        pcaScoresPlotRange$xMax <<- maxX
        pcaScoresPlotRange$xInterval <<- c(minX, maxX)
        pcaScoresPlotRange$xIntervalSize <<- maxX - minX
        pcaScoresPlotRange$yMin <<- minY
        pcaScoresPlotRange$yMax <<- maxY
        pcaScoresPlotRange$yInterval <<- c(minY, maxY)
        pcaScoresPlotRange$yIntervalSize <<- maxY - minY
      }
    })
    #obsPCAscoresRangeUpdate <- observe({
    #  xInterval <- pcaScoresPlotRange$xInterval
    #  yInterval <- pcaScoresPlotRange$yInterval
    #  
    #  if(!state$showPCAplotPanel)
    #    return()
    #  
    #  print(paste("observe pcaScoresPlotRange", paste(xInterval, collapse = ", "), paste(yInterval, collapse = ", ")))
    #  
    #  ## plot PCA
    #  #drawPcaScoresPlot(consoleInfo = "range update output$plotPcaScores")
    #})
    obsPCAloadingsClick <- observeEvent(input$plotPcaLoadings_click, {
      clickX <- input$plotPcaLoadings_click$x
      clickY <- input$plotPcaLoadings_click$y
      plotWidth  <- session$clientData$output_plotPcaLoadings_width
      plotHeight <- session$clientData$output_plotPcaLoadings_height
      brush <- input$plotPcaLoadings_brush
      
      if(!is.null(brush))
        return()
      if(is.null(clickX) | is.null(clickY))
        return()
      print(paste("Observe PCA Loadings click", clickX, clickY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = clickX, mouseY = clickY, poiCoordinatesX = pcaDataList$pcaLoadingsX, poiCoordinatesY = pcaDataList$pcaLoadingsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaLoadingsPlotRange$xIntervalSize, plotRangeY = pcaLoadingsPlotRange$yIntervalSize
      )
      print(paste("Observe PCA Loadings hover", clickX, clickY, minimumIndex))
      
      if(is.null(minimumIndex)){
        ## reset stuff
        
        fragmentsX <<- NULL
        fragmentsY <<- NULL
        fragmentsColor <<- NULL
        fragmentsDiscriminativity <<- NULL
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        selectionByAnalysisReset()
        selectionByFragmentReset()
      } else {
        ## loadng selection
        precursorIndex <- filterPca$filter[[minimumIndex]]
        ## fetch ms2 spectrum
        resultObj <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        
        ## keep fragment selection
        selectionFragmentSelectedFragmentIndexNew <- NULL
        if(!is.null(selectionFragmentSelectedFragmentIndex)){
          fragmentMass <- fragmentsX[[selectionFragmentSelectedFragmentIndex]]
          if(fragmentMass %in% resultObj$fragmentMasses)
            selectionFragmentSelectedFragmentIndexNew <- which(resultObj$fragmentMasses %in% fragmentMass)
        }
        
        fragmentsX <<- resultObj$fragmentMasses
        fragmentsY <<- resultObj$fragmentAbundances
        fragmentsColor <<- resultObj$fragmentColor
        fragmentsDiscriminativity <<- resultObj$fragmentDiscriminativity
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        fragmentsColorHovered <<- NULL
        fragmentsDiscriminativityHovered <<- NULL
        
        selectionByPca(minimumIndex)
        
        ## update the selected fragment in case of overlapping spectra
        if(!is.null(selectionFragmentSelectedFragmentIndexNew))
          selectionByFragment(selectionFragmentSelectedFragmentIndexNew)
        else
          selectionByFragmentReset()
      }
      
      ## MS2 plot
      resetMS2PlotRange()
      drawMS2Plot(consoleInfo = "PCA loadings click output$plotMS2")
      drawPcaLoadingsPlot(consoleInfo = "PCA loadings click output$plotPcaLoadings")
      
      if(state$showHCAplotPanel)
        ## update dendrogram plot
        drawDendrogramPlot(consoleInfo = "PCA loadings click output$plotDendrogram", withHeatmap = TRUE)
    })
    obsPCAloadingsHover <- observeEvent(input$plotPcaLoadings_hover, {
      hoverX <- input$plotPcaLoadings_hover$x
      hoverY <- input$plotPcaLoadings_hover$y
      plotWidth  <- session$clientData$output_plotPcaLoadings_width
      plotHeight <- session$clientData$output_plotPcaLoadings_height
      
      if(is.null(hoverX) | is.null(hoverY))
        return()
      print(paste("Observe PCA Loadings hover", hoverX, hoverY))
      
      #################################################
      ## decide whether the hover is close enough to trigger event
      minimumIndex <- getSelectedPOI_XY(
        mouseX = hoverX, mouseY = hoverY, poiCoordinatesX = pcaDataList$pcaLoadingsX, poiCoordinatesY = pcaDataList$pcaLoadingsY,
        plotWidth = plotWidth, plotHeight = plotHeight, plotRangeX = pcaLoadingsPlotRange$xIntervalSize, plotRangeY = pcaLoadingsPlotRange$yIntervalSize
      )
      if(is.null(minimumIndex)){
        ## no loading hovered
        fragmentsXhovered <<- NULL
        fragmentsYhovered <<- NULL
        output$information <- renderText({
          print(paste("update output$information PCA Loadings hover ", minimumIndex, sep = ""))
          paste("", sep = "")
        })
      } else {
        print(paste("Observe PCA Loadings hover", hoverX, hoverY, minimumIndex))
        
        #################################################
        ## fetch ms2 spectrum
        precursorIndex <- filterPca$filter[[minimumIndex]]
        resultObj <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndex)
        if(all(!is.null(selectionAnalysisPcaLoadingSet), minimumIndex == selectionAnalysisPcaLoadingSet)){
          ## reverse MS2 to clicked stuff
          fragmentsXhovered <<- NULL
          fragmentsYhovered <<- NULL
          fragmentsColorHovered <<- NULL
          fragmentsDiscriminativityHovered <<- NULL
        } else {
          fragmentsXhovered <<- resultObj$fragmentMasses
          fragmentsYhovered <<- resultObj$fragmentAbundances
          fragmentsColorHovered <<- resultObj$fragmentColor
          fragmentsDiscriminativityHovered <<- resultObj$fragmentDiscriminativity
        }
        
        output$information <- renderText({
          print(paste("update output$information PCA Loadings hover ", precursorIndex, sep = ""))
          paste(resultObj$infoText, sep = "")
        })
      }
      
      drawMS2Plot(consoleInfo = "loadings hover output$plotMS2")
      
      output$tip <- renderText({
        print(paste("update output$tip"))
        paste(
          "Hover or click a loadings node to view information about the corresponding MS\u00B9 feature.", 
          "Brush and double-click to zoom in.", 
          "Double-click to zoom out.", 
          sep = "\n"
        )
      })
    })
    obsPCAloadingsDblClick <- observeEvent(input$plotPcaLoadings_dblclick, {
      brush <- input$plotPcaLoadings_brush
      
      print(paste("observe PCAloadings dblclick", is.null(brush)))
      
      if (!is.null(brush)) {
        ## set range
        pcaLoadingsPlotRange$xMin <<- brush$xmin
        pcaLoadingsPlotRange$xMax <<- brush$xmax
        pcaLoadingsPlotRange$xInterval <<- c(brush$xmin, brush$xmax)
        pcaLoadingsPlotRange$xIntervalSize <<- brush$xmax - brush$xmin
        pcaLoadingsPlotRange$yMin <<- brush$ymin
        pcaLoadingsPlotRange$yMax <<- brush$ymax
        pcaLoadingsPlotRange$yInterval <<- c(brush$ymin, brush$ymax)
        pcaLoadingsPlotRange$yIntervalSize <<- brush$ymax - brush$ymin
      } else {
        ## reset range
        minX <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
        maxX <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionOne])
        minY <- min(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
        maxY <- max(pcaDataList$pcaObj$loadings[, pcaDataList$dimensionTwo])
        
        pcaLoadingsPlotRange$xMin <<- minX
        pcaLoadingsPlotRange$xMax <<- maxX
        pcaLoadingsPlotRange$xInterval <<- c(minX, maxX)
        pcaLoadingsPlotRange$xIntervalSize <<- maxX - minX
        pcaLoadingsPlotRange$yMin <<- minY
        pcaLoadingsPlotRange$yMax <<- maxY
        pcaLoadingsPlotRange$yInterval <<- c(minY, maxY)
        pcaLoadingsPlotRange$yIntervalSize <<- maxY - minY
      }
    })
    #obsPCAloadingsRangeUpdate <- observe({
    #  xInterval <- pcaLoadingsPlotRange$xInterval
    #  yInterval <- pcaLoadingsPlotRange$yInterval
    #  
    #  if(!state$showPCAplotPanel)
    #    return()
    #  
    #  print(paste("observe pcaLoadingsPlotRange ", paste(xInterval, collapse = ", "), "; ", paste(yInterval, collapse = ", "), sep = ""))
    #  
    #  ## plot PCA
    #  #drawPcaLoadingsPlot(consoleInfo = "range update output$plotPcaLoadings")
    #})
    ## fragment plot
    obsFragmentPlotdblClick <- observeEvent(input$fragmentPlot_dblclick, {
      brush <- input$fragmentPlot_brush
      print(paste("observe fragmentPlot dblclick", paste(brush, collapse = "; ")))
      
      if (!is.null(brush)) {
        ## set range
        min <- brush$xmin
        max <- brush$xmax
      } else {
        ## reset range
        min <- min(dataList$masses)
        max <- max(dataList$masses)
      }
      
      fragmentPlotRange$xMin <<- min
      fragmentPlotRange$xMax <<- max
      fragmentPlotRange$xInterval <<- c(min, max)
      fragmentPlotRange$xIntervalSize <<- max - min
    })
    #obsFragmentPlotRangeUpdate <- observe({
    #  xInterval <- fragmentPlotRange$xInterval
    #  
    #  print(paste("observe fragmentPlotRange", xInterval))
    #  
    #  #drawFragmentPlot(consoleInfo = "update range output$fragmentPlotDendrogram")
    #})
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
    ## listen to annotation events
    obsSetPresentAnnoPrimary <- observeEvent(input$setPresentAnnotationPrimary, {
      value     <- input$presentAnnotationValue
      
      drawPlots <- as.numeric(input$setPresentAnnotationPrimary)
      if(drawPlots == setPresentAnnotationPrimaryValue)
        return()
      setPresentAnnotationPrimaryValue <<- drawPlots
      print(paste("Observe setPresentAnnotationPrimary", drawPlots))
      
      setAnnotationPrimary(precursorSet = selectedPrecursorSet, annotationValue = value)
    })
    obsRemovePresentAnno <- observeEvent(input$removePresentAnnotation, {
      value     <- input$presentAnnotationValue
      
      drawPlots <- as.numeric(input$removePresentAnnotation)
      ## XXX why?
      #if(drawPlots == removePresentAnnotationValue)
      #  return()
      #removePresentAnnotationValue <<- drawPlots
      print(paste("Observe removePresentAnnotation", drawPlots))
      
      removeAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value)
    })
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
      value <- input$newAnnotationValue
      color <- input$newAnnotationColor
      
      submitNewAnnotation <- as.numeric(input$submitNewAnnotation)
      print(paste("Observe submitNewAnnotation", submitNewAnnotation))
      
      if(submitNewAnnotation == submitNewAnnotationValue)
        return()
      submitNewAnnotationValue <<- submitNewAnnotation
      
      addAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value, annotationColor = color)
    })
    obsAddPresentAnno <- observeEvent(input$submitPreviousAnnotation, {
      value <- input$previousAnnotationValue
      
      submitPreviousAnnotation <- as.numeric(input$submitPreviousAnnotation)
      if(submitPreviousAnnotation == submitPreviousAnnotationValue)
        return()
      submitPreviousAnnotationValue <<- submitPreviousAnnotation
      print(paste("Observe submitPreviousAnnotation", submitPreviousAnnotation))
      
      color <- dataList$annoPresentColorsList[[match(x = value, table = dataList$annoPresentAnnotationsList)]]
      addAnnotation(precursorSet = selectedPrecursorSet, annotationValue = value, annotationColor = color)
    })
    obsIgnoreValueChanged <- observeEvent(input$updateArtifactsFromCheckboxes, {
      updateArtifactsFromCheckboxes <- as.numeric(input$updateArtifactsFromCheckboxes)
      
      #################################################
      ## check if button was hit
      if(updateArtifactsFromCheckboxes == updateArtifactsFromCheckboxesButtonValue)
        return()
      updateArtifactsFromCheckboxesButtonValue <<- updateArtifactsFromCheckboxes
      
      ## get and process values
      vals <- getInputValues(id = paste(state$selectedSelection, "Ignore", sep = "_"), counter = selectedTable_id, len = nrow(selectedTable))
      
      if(all(is.na(vals)))
        return()
      
      nas <- is.na(vals)
      vals[nas] <- dataList$annoArrayIsArtifact[selectedPrecursorSet][nas]
      vals <- as.logical(vals)
      
      setArtifactState(selectedPrecursorSet, vals)
    })
    obsSampleValuesChanged <- observeEvent(input$updateSampleTable, {
      updateSampleTable <- as.numeric(input$updateSampleTable)
      
      #################################################
      ## check if button was hit
      if(updateSampleTable == updateSamplesButtonValue)
        return()
      updateSamplesButtonValue <<- updateSampleTable
      
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
      ## TODO filter again...?
    })
    sampleExcludeClicked <- function(sampleId){
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
    ## display of tabs
    observe({
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Filter']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='MS/MS filter']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Sample filter']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='PCA']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='HCA']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Search']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Annotation']")
      #toggle(condition = FALSE, selector = "#runTabs li a[data-value='Annotation']")
      toggle(condition = !is.null(state$importedOrLoadedFile_s_), selector = "#runTabs li a[data-value='Project']")
    })
    
    ## suspend observer
    session$onSessionEnded(function() {
      ## sidepanel
      obsTabs$suspend()
      obsChangePlot$suspend()
      obsChangeSelection$suspend()
      obsClearSelection$suspend()
      obsFile$suspend()
      obsLoadProjectData$suspend()
      #obsExampleDataSelection$suspend()
      obsLoadExampleData$suspend()
      obsImportMs1DataFile$suspend()
      obsImportMs2DataFile$suspend()
      obsImportMs1Ms2Data$suspend()
      obsImportMs2Data$suspend()
      obsApplyImportParameterFile$suspend()
      obsFileInputSelection$suspend()
      obsShowSideBar$suspend()
      obsShowPlotControls$suspend()
      obsShowClusterLabels$suspend()
      obsHcaPrecursorLabels$suspend()
      obsShowScoresLabels$suspend()
      #observePcaLoadingsProperties$suspend()
      obsLoadingsLabels$suspend()
      obsShowLoadingsFeatures$suspend()
      #obsShowLoadingsFeaturesAnnotated$suspend()
      #obsShowLoadingsFeaturesUnannotated$suspend()
      #obsShowLoadingsFeaturesSelected$suspend()
      #obsShowLoadingsFeaturesUnselected$suspend()
      obsShowLoadingsAbundance$suspend()
      ## filter
      obsClearSearch$suspend()
      obsApplySearch$suspend()
      obsApplyGlobalMS2filters$suspend()
      obsClearGlobalMS2filters$suspend()
      obsApplyHcaFilters$suspend()
      obsClearHcaFilters$suspend()
      obsApplyPcaFilters$suspend()
      obsClearPcaFilters$suspend()
      observeGroupSet$suspend()
      ## draw
      obsDrawHCA$suspend()
      obsDrawPCA$suspend()
      ## dendrogram
      if(FALSE) obsDendLabelsHeatmap$suspend() ## plotly
      obsDendrogramHover$suspend()
      obsDendrogramClick$suspend()
      obsDendrogramdblClick$suspend()
      #obsDendrogramRangeUpdate$suspend()
      ## heatmap
      obsHeatmaphover$suspend()
      ## MS2
      obsMS2hover$suspend()
      obsMS2click$suspend()
      obsMS2dblClick$suspend()
      #obsMS2rangeUpdate$suspend()
      ## PCA
      obsPCAscoresHover$suspend()
      obsPCAscoresDblClick$suspend()
      #obsPCAscoresRangeUpdate$suspend()
      obsPCAloadingsClick$suspend()
      obsPCAloadingsHover$suspend()
      obsPCAloadingsDblClick$suspend()
      #obsPCAloadingsRangeUpdate$suspend()
      ## fragment plot
      obsFragmentPlotdblClick$suspend()
      #obsFragmentPlotRangeUpdate$suspend()
      ## anno
      obsSetPresentAnnoPrimary$suspend()
      obsRemovePresentAnno$suspend()
      obsToggleAddNewAnnoButton$suspend()
      obsAddNewAnno$suspend()
      obsAddPresentAnno$suspend()
      obsIgnoreValueChanged$suspend()
      obsSampleValuesChanged$suspend()
      obsShowHCAplotPanel$suspend()
      obsShowPCAplotPanel$suspend()
      obsApplyConfirmedAnnotations$suspend()
      obsToggleConfirmAnnoButton$suspend()
      stopApp()
    })
    
    #########################################################################################
    #########################################################################################
    ## direct output rendering
    output$fileInfo <- renderText({
      print(paste("init output$fileInfo"))
      paste("Please select a project file and press 'Load project data'")
    })
    output$information <- renderText({
      print(paste("init output$information"))
      ""
    })
    output$rInfo <- renderText({
      print(paste("init rInfo"))
      #paste(R.Version()$version.string, "\nwd: ", getwd(), sep = "")
      paste(
        R.Version()$version.string, 
        "\nMetFamily build: ", metFamilyBuilt, "-", system(command = "hostname", intern = TRUE),
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
    ## direct output values
    output$showGUI <- reactive({
      print("update output$showGUI")
      output$information <- renderText({
        print(paste("init information", sep = ""))
        paste("Please perform ploting.", sep = "")
      })
      return(!is.null(state$importedOrLoadedFile_s_))
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
    output$globalMS2filterValid <- reactive({
      print(paste("reactive update globalMS2filterValid", state$globalMS2filterValid))
      return(state$globalMS2filterValid)
    })
    output$hcaFilterValid <- reactive({
      print(paste("reactive update hcaFilterValid", state$hcaFilterValid))
      return(state$hcaFilterValid)
    })
    output$pcaFilterValid <- reactive({
      print(paste("reactive update pcaFilterValid", state$pcaFilterValid))
      return(state$pcaFilterValid)
    })
    output$precursorSetSelected <- reactive({
      print(paste("reactive update precursorSetSelected", state$precursorSetSelected))
      return(state$precursorSetSelected)
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
    output$searchfilterValid <- reactive({
      print(paste("reactive update searchfilterValid", state$searchfilterValid))
      return(state$searchfilterValid)
    })
    output$filterSearchActive <- reactive({
      print(paste("reactive update filterSearchActive", state$filterSearchActive))
      return(state$filterSearchActive)
    })
    output$selectedSelection <- reactive({
      print(paste("reactive update selectedSelection", state$selectedSelection))
      return(state$selectedSelection)
    })
    
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
    outputOptions(output, 'globalMS2filterValid',    suspendWhenHidden=FALSE)
    outputOptions(output, 'hcaFilterValid',          suspendWhenHidden=FALSE)
    outputOptions(output, 'pcaFilterValid',          suspendWhenHidden=FALSE)
    outputOptions(output, 'precursorSetSelected',    suspendWhenHidden=FALSE)
    outputOptions(output, 'searchfilterValid',       suspendWhenHidden=FALSE)
    outputOptions(output, 'filterSearchActive',      suspendWhenHidden=FALSE)
    outputOptions(output, 'plotHcaShown',            suspendWhenHidden=FALSE)
    outputOptions(output, 'plotPcaShown',            suspendWhenHidden=FALSE)
    outputOptions(output, 'plotAnnotationShown',     suspendWhenHidden=FALSE)
    outputOptions(output, 'selectedSelection',       suspendWhenHidden=FALSE)
    
    #########################################################################################
    #########################################################################################
    ## ui generation
    output$runRightColumn <- renderUI({
      print(paste("### GUI ### Generate right GUI"))
      
      column(width = state$runRightColumnWidth,
         #########################################################################################
         ## show side bar and change plot
         conditionalPanel(
           condition = "(output.showHCAplotPanel & output.analysisType == 'HCA') | (output.showPCAplotPanel & output.analysisType == 'PCA') | (output.showAnnotationplotPanel & output.analysisType == 'Annotation')",
           #condition = "output.showHCAplotPanel | output.showPCAplotPanel",
           fluidRow(
             column(width = 6,
                div(style="float:right",
                    bsTooltip(id = "showSideBar", title = "Display or hide the side bar", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "showSideBar", label = "Display side bar", value = showSideBar)
                )
             ),##column
             column(width = 6,
                conditionalPanel(
                  condition = "(output.showHCAplotPanel & output.showPCAplotPanel) | (output.showHCAplotPanel & output.showAnnotationplotPanel) | (output.showAnnotationplotPanel & output.showPCAplotPanel)",
                  div(style="float:left",
                      bsTooltip(id = "changePlot", title = "Switch between HCA plots and PCA plots", placement = "bottom", trigger = "hover"),
                      #c("Display HCA", "Display PCA", "Display Annotation")
                      radioButtons(inputId = "changePlot", label = NULL, choices = plotsToShow, inline = TRUE, selected = plotToShow)
                  )
                )##conditional
             )##column
           )##row
         ),##conditional
         ##############################################################################################
         ##############################################################################################
         ## plots
         
         ##############################################################################################
         ## plots controls
         conditionalPanel(## plot properties
           condition = '(output.analysisType == "HCA" & output.showHCAplotPanel) | (output.analysisType == "PCA" & output.showPCAplotPanel)',
           wellPanel(
             fluidRow(
               column(width = 6,
                      div(style="float:left",
                          h4("Plot controls")
                      )
               ),##column
               column(width = 6,
                      div(style="float:right",
                          bsTooltip(id = "showPlotControls", title = "Display control panels for the plots below", placement = "bottom", trigger = "hover"),
                          checkboxInput(inputId = "showPlotControls", label = "Show plot controls", value = input$showPlotControls)
                      )
               )##column
             ),##row
             fluidRow(
               ################################################
               ## HCA plot controls
               conditionalPanel(## dendrogram properties
                 condition = 'output.analysisType == "HCA" & output.showHCAplotPanel & input.showPlotControls',
                 column(width = 3,
                        h5("Heatmap content"),
                        tags$div(
                          title="Please select the abundance information you would like to display in the heatmap below the dendrogram",
                          radioButtons(inputId = "heatmapContent", label = NULL, choices = c("Log-fold-change", "Abundance by group", "Abundance by sample"), selected = input$heatmapContent)
                        )
                 ),
                 column(width = 3,
                        h5("Heatmap ordering"),
                        tags$div(
                          title="Please select the mode of ordering the heatmap rows below the dendrogram",
                          radioButtons(inputId = "heatmapOrdering", label = NULL, choices = c("Specified order", "MS1 clustering"), selected = input$heatmapOrdering)
                        )
                 ),
                 column(width = 3,
                        h5("HCA dendrogram"),
                        bsTooltip(id = "showClusterLabels", title = "Display the labels of cluster nodes and MS\u00B9 feature nodes representing the number of characteristic fragments", placement = "bottom", trigger = "hover"),
                        checkboxInput(inputId = "showClusterLabels", label = "Show node labels", value = input$showClusterLabels)
                 ),
                 column(width = 3,
                        h5("Precursor label"),
                        tags$div(
                          title="Please select the information you would like to display in the labels below the precursors",
                          radioButtons(inputId = "hcaPrecursorLabels", label = NULL, choices = c("m/z / RT", "Metabolite name", "Metabolite family"), selected = input$hcaPrecursorLabels)
                        )
                 )
               ),## conditional
               ################################################
               ## PCA plot controls
               conditionalPanel(## scores / loadings properties
                 condition = 'output.analysisType == "PCA" & output.showPCAplotPanel & input.showPlotControls',
                 column(width = 3,
                        h5("PCA scores"),
                        bsTooltip(id = "showScoresLabels", title = "Display scores labels", placement = "bottom", trigger = "hover"),
                        checkboxInput(inputId = "showScoresLabels", label = "Show labels", value = input$showScoresLabels)
                 ),
                 column(width = 3,
                        h5("PCA loadings labels"),
                        tags$div(
                          title="Please select the information you would like to display in the precursor labels of the loadings",
                          radioButtons(inputId = "loadingsLabels", label = NULL, choices = c("None", "m/z / RT", "Metabolite name", "Metabolite family"), selected = input$loadingsLabels)
                        )
                 ),
                 column(width = 3,
                        h5("PCA loadings shown"),
                        tags$div(
                          title="Please select the MS\u00B9 features you would like to display in the loadings plot",
                          bsTooltip(id = "showLoadingsAbundance", title = "Use abundance in MS\u00B9 to scale the size of loadings nodes", placement = "bottom", trigger = "hover"),
                          checkboxGroupInput(inputId = "showLoadingsFeatures", label = NULL, 
                                             choices = c("Annotated", "Not Annotated", "Selected", "Not Selected"), 
                                             selected = c(
                                               ifelse(test = state$showLoadingsFeaturesAnnotated,   yes = "Annotated", no = ""),
                                               ifelse(test = state$showLoadingsFeaturesUnannotated, yes = "Not Annotated", no = ""),
                                               ifelse(test = state$showLoadingsFeaturesSelected,    yes = "Selected", no = ""),
                                               ifelse(test = state$showLoadingsFeaturesUnselected,  yes = "Not Selected", no = "")
                                             )
                                            )
                          #checkboxInput(inputId = "showLoadingsFeaturesAnnotated",   label = "Annotated",     value = input$showLoadingsFeaturesAnnotated),
                          #checkboxInput(inputId = "showLoadingsFeaturesUnannotated", label = "Not Annotated", value = input$showLoadingsFeaturesUnannotated),
                          #checkboxInput(inputId = "showLoadingsFeaturesSelected",    label = "Selected",      value = input$showLoadingsFeaturesSelected),
                          #checkboxInput(inputId = "showLoadingsFeaturesUnselected",  label = "Not Selected",  value = input$showLoadingsFeaturesUnselected)
                        )
                 ),
                 column(width = 3,
                        h5("PCA loadings size"),
                        bsTooltip(id = "showLoadingsAbundance", title = "Use abundance in MS\u00B9 to scale the size of loadings nodes", placement = "bottom", trigger = "hover"),
                        checkboxInput(inputId = "showLoadingsAbundance", label = "Scale by abundance", value = input$showLoadingsAbundance)
                 )
               )## conditional
             )## row
           )## well
         ),## conditional
         ##############################################################################################
         ## plots
         fluidRow(
           column(width = 12-state$legendColumnWidth,
             ##############################################################################################
             ## HCA plots
             conditionalPanel(
               condition = 'output.analysisType == "HCA" & output.showHCAplotPanel',
               fluidRow(
                 #plotlyOutput(
                #   height = state$dendrogramHeatmapHeight, 
                #   #height = 500, 
                #   outputId = "plotDendrogram"
                # )
                 plotOutput(height = 500, 
                            outputId = "plotDendrogram", 
                            #hover    = "plotDendrogram_hover", 
                            hover    = hoverOpts(
                              id = "plotDendrogram_hover",
                              delay = 50, 
                              delayType = "debounce"
                            ),
                            click    = "plotDendrogram_click",
                            dblclick = "plotDendrogram_dblclick",
                            #brush    = "plotDendrogram_brush"
                            brush    = brushOpts(
                              id = "plotDendrogram_brush",
                              resetOnNew = TRUE,
                              direction = "x",
                              delay = 00,
                              delayType = "debounce"
                            )
                 ),
                 plotOutput(height = state$heatmapHeight, 
                            outputId = "plotHeatmap",
                            #hover    = "plotHeatmap_hover", 
                            hover    = hoverOpts(
                              id = "plotHeatmap_hover",
                              delay = 50, 
                              delayType = "debounce"
                            )
                            #click = "plotHeatmap_click"
                 )
               )## row
             ),## conditional
             ##############################################################################################
             ## PCA plots
             conditionalPanel(
               condition = 'output.analysisType == "PCA" & output.showPCAplotPanel',
               fluidRow(
                 column(width = 6,
                    plotOutput(height = 500, 
                               outputId = "plotPcaScores", 
                               #hover    = "plotPcaScores_hover",
                               hover    = hoverOpts(
                                 id = "plotPcaScores_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "plotPcaScores_click",
                               dblclick = "plotPcaScores_dblclick",
                               #brush    = "plotPcaScores_brush"
                               brush = brushOpts(
                                 id = "plotPcaScores_brush",
                                 resetOnNew = TRUE,
                                 direction = "xy",
                                 delay = 00,
                                 delayType = "debounce"
                               )
                    )
                 ),## column
                 column(width = 6,
                    plotOutput(height = 500, 
                               outputId = "plotPcaLoadings", 
                               #hover    = "plotPcaLoadings_hover",
                               hover    = hoverOpts(
                                 id = "plotPcaLoadings_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "plotPcaLoadings_click",
                               dblclick = "plotPcaLoadings_dblclick",
                               #brush    = "plotPcaScores_brush"
                               brush = brushOpts(
                                 id = "plotPcaLoadings_brush",
                                 resetOnNew = TRUE,
                                 direction = "xy",
                                 delay = 00,
                                 delayType = "debounce"
                               )
                    )
                 )## column
               )## row
             )## conditional
           ),##column for plot controls and plots
           column(width = state$legendColumnWidth,
              conditionalPanel(
                condition = 'output.analysisType == "HCA" & output.showHCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotAnnoLegendHCA", height = state$annotationLegendHeightHca)
                )
              ),## conditional
              conditionalPanel(
                condition = 'output.analysisType == "PCA" & output.showPCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotAnnoLegendPCA", height = state$annotationLegendHeightPca)
                )
              ),## conditional
              conditionalPanel(
                condition = '(output.analysisType == "HCA" & output.showHCAplotPanel) | (output.analysisType == "PCA" & output.showPCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "calcPlotDendrogramLegend", height = 80)
                )
              ),## conditional
              conditionalPanel(
                condition = 'output.analysisType == "HCA" & output.showHCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotHeatmapLegend", height = 150)
                )
              ),## conditional
              conditionalPanel(## loadings properties
                condition = 'output.analysisType == "PCA" & output.showPCAplotPanel',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotScoresGroupsLegend", height = state$scoresGroupsLegendHeight)
                )
              ),
              conditionalPanel(
                condition = '(output.analysisType == "HCA" & output.showHCAplotPanel) | (output.analysisType == "PCA" & output.showPCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotMS2Legend", height = 80)
                )
              ),## conditional
              conditionalPanel(
                condition = '(output.analysisType == "HCA" & output.showHCAplotPanel)',
                splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput(outputId = "plotFragmentDiscriminativityLegend", height = 100)
                )
              )## conditional
           )## column
         ),## row
         #########################################################################################
         ## MS2 plot and info
         conditionalPanel(
           condition = '(output.showHCAplotPanel & output.analysisType == "HCA") | (output.showPCAplotPanel & output.analysisType == "PCA")',
           fluidRow(
             plotOutput(height = 250, 
                        outputId = "plotMS2",
                        #hover    = "plotMS2_hover",
                        hover    = hoverOpts(
                          id = "plotMS2_hover",
                          delay = 50, 
                          delayType = "debounce"
                        ),
                        click    = "plotMS2_click",
                        dblclick = "plotMS2_dblclick",
                        #brush    = "plotMS2_brush",
                        brush = brushOpts(
                          id = "plotMS2_brush",
                          resetOnNew = TRUE,
                          direction = "x",
                          delay = 00,
                          delayType = "debounce"
                        )
             ),
             ##############################################################################################
             ## infos
             wellPanel(
               h4("Information"),
               bsTooltip(id = "information", title = "Information about items in the plot", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("information"),
               h4("Tip"),
               bsTooltip(id = "tip", title = "Information about operating options", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("tip")
             )## well
           )## row
         ),## conditional
         ##############################################################################################
         ## precursor set selection and annotation
         ## change selection
         conditionalPanel(
           condition = '(output.showHCAplotPanel & output.analysisType == "HCA") | (output.showPCAplotPanel & output.analysisType == "PCA")',
           fluidRow(
             wellPanel(
               h4("MS\u00B9 feature selections"),
               bsTooltip(id = "changeSelection", title = "Switch MS\u00B9 feature selection", placement = "bottom", trigger = "hover"),
               radioButtons(inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), selected = changeSelectionCurrentSelection, inline = TRUE),
               bsTooltip(id = "selectionInfo", title = "The number of MS\u00B9 features in the current selection", placement = "bottom", trigger = "hover"),
               hr(),
               verbatimTextOutput("selectionInfo"),
               conditionalPanel(
                 condition = 'output.precursorSetSelected',
                 tabsetPanel(id = "precursorSelectionTabs",
                   tabPanel(title = precursorSelectionTabSelection, 
                     wellPanel(
                       ## selection infos
                       tags$div(
                         style="margin-bottom:5px;",
                         bsTooltip(id = "metFragLink", title = "Press to send the current MS\u00B9 feature as well as the corresponding MS/MS spectrum to MetFrag", placement = "bottom", trigger = "hover"),
                         htmlOutput(outputId = "metFragLink")
                       ),
                       bsTooltip(id = "downloadSelectedPrecursors", title = "Download a project file which is reduced to the selected set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
                       downloadButton('downloadSelectedPrecursors', 'Download reduced project file'),
                       bsTooltip(id = "clearSelection", title = "Press to clear this selection", placement = "bottom", trigger = "hover"),
                       actionButton(inputId = "clearSelection", label = "Clear selection", class="btn-success")
                     )## well
                   ),## tab
                   tabPanel(title = precursorSelectionTabAnnotation, 
                      wellPanel(
                        h4("Present annotation(s)"),
                        fluidRow(
                          column(
                            width = 3,
                            bsTooltip(id = "presentAnnotationValue", title = "The set of present annotations for the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            selectInput(inputId = "presentAnnotationValue", label = NULL, choices = c("[init]"), selectize = FALSE)
                          ),## column
                          column(
                            width = 3,
                            bsTooltip(id = "setPresentAnnotationPrimary", title = "Sets the selected annotation primary for the set of selected MS\u00B9 features; i.e. this annotation will be used preferentially for coloring in HCA and PCA", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "setPresentAnnotationPrimary", label = "Set primary", class="btn-success")
                          ),## column
                          column(
                            width = 6,
                            bsTooltip(id = "removePresentAnnotation", title = "Removes the selected annotation from the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "removePresentAnnotation", label = "Remove annotation", class="btn-success")
                          )## column
                        ),##row
                        fluidRow(
                          column(
                            width = 6,
                            h4("Add new annotation"),
                            bsTooltip(id = "newAnnotationValue", title = "The name of this annotation", placement = "bottom", trigger = "hover"),
                            textInput(inputId = "newAnnotationValue", label = "Type new annotation"),
                            bsTooltip(id = "newAnnotationColor", title = "The color of this annotation", placement = "bottom", trigger = "hover"),
                            colourInput(inputId = "newAnnotationColor", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = colorPalette()),
                            bsTooltip(id = "submitNewAnnotation", title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "submitNewAnnotation", label = "Add new annotation", class="btn-success")
                          ),
                          column(
                            width = 6,
                            h4("Add previous annotation"),
                            bsTooltip(id = "previousAnnotationValue", title = "The set of annotations which have been assigned before", placement = "bottom", trigger = "hover"),
                            selectInput(inputId = "previousAnnotationValue", label = "Select previous annotation", choices = c("Artifact"), selectize = FALSE),
                            bsTooltip(id = "submitPreviousAnnotation", title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            actionButton(inputId = "submitPreviousAnnotation", label = "Add previous annotation", class="btn-success")
                          )
                        )
                      )## well
                   ),## tab
                   tabPanel(title = precursorSelectionTabTable, 
                      wellPanel(
                        h4("Selected MS\u00B9 features"),
                        bsTooltip(id = "updateArtifactsFromCheckboxes", title = "Adds the annotation \\'ignore\\' to the set of checked MS\u00B9 features", placement = "bottom", trigger = "hover"),
                        actionButton(inputId = "updateArtifactsFromCheckboxes", label = "Apply annotation 'Ignore' to MS\u00B9 features", class="btn-success"),
                        DT::dataTableOutput("ms1FeatureTable")
                      )## well
                   )## tab
                 )## tab set
               )## conditional
             )## well
           )## row
         ),##conditional
         ##############################################################################################
         ## precursor set selection and annotation
         ## change selection
         conditionalPanel(
           condition = "(output.showAnnotationplotPanel & output.analysisType == 'Annotation')",
           fluidRow(
             wellPanel(
               ## TODO 999
               h4("Metabolite family selection"),
               bsTooltip(id = "familyCount", title = "The number of metabolite families with one or more potential MS\u00B9 features", placement = "bottom", trigger = "hover"),
               verbatimTextOutput(outputId = "familyCount"),
               DT::dataTableOutput("annotationResultTableClass"),
               h4("MS\u00B9 feature annotation"),
               bsTooltip(id = "classToSpectraCount", title = "The number of potential MS\u00B9 feature hits", placement = "bottom", trigger = "hover"),
               verbatimTextOutput(outputId = "classToSpectraCount"),
               DT::dataTableOutput("annotationResultTableFeature"),
               h4("MS\u00B9 feature spectrum versus Metabolite family"),
               plotOutput(height = 250, 
                          outputId = "plotMS2vsClass",
                          #hover    = "plotMS2_hover",
                          #hover    = hoverOpts(
                          #  id = "plotMS2_hover",
                          #  delay = 50, 
                          #  delayType = "debounce"
                          #),
                          #click    = "plotMS2_click",
                          dblclick = "plotMS2vsClass_dblclick",
                          ##brush    = "plotMS2_brush",
                          brush = brushOpts(
                            id = "plotMS2vsClass_brush",
                            resetOnNew = TRUE,
                            direction = "x",
                            delay = 00,
                            delayType = "debounce"
                          )
               ),
               bsTooltip(id = "newAnnotationValue2", title = "The name of this annotation", placement = "bottom", trigger = "hover"),
               textInput(inputId = "newAnnotationValue2", label = "Type new annotation"),
               bsTooltip(id = "newAnnotationColor2", title = "The color of this annotation", placement = "bottom", trigger = "hover"),
               colourInput(inputId = "newAnnotationColor2", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = colorPalette()),
               bsTooltip(id = "confirmAnnotation", title = "Applies the metabolite family annotation to all confirmed MS\u00B9 features", placement = "bottom", trigger = "hover"),
               actionButton(inputId = "confirmAnnotation", label = "Apply confirmed annotations", class="btn-success")
             )## well
           )## row
         )##conditional
      )##column
    })
    
    #########################################################################################
    #########################################################################################
    ## download
    timeStampForFiles <- function(){
      timeStamp <- gsub(" ", "_", gsub(":", ".", Sys.time()))
      return(timeStamp)
    }
    createImportParameterSetExportFileName <- function(){
      fileProjectName <- dataList$importParameterSet$projectName
      fileProjectName <- gsub(" ", "_", gsub(":", ".", fileProjectName))
      fileName <- paste(timeStampForFiles(), "_", fileProjectName, "_import_parameters.txt", sep = "")
      return(fileName)
    }
    createExportMatrixName <- function(){
      fileName <- paste(timeStampForFiles(), "_selectedPrecursorMatrix.csv.gz", sep = "")
      return(fileName)
    }
    createExportImageName <- function(item, extension){
      fileName <- paste(timeStampForFiles(), "_", item, ".", extension, sep = "")
      return(fileName)
    }
    createExportDistanceMatrixName <- function(distanceMeasure){
      fileName <- paste(timeStampForFiles(), "_distanceMatrix_", distanceMeasure, ".csv", sep = "")
      return(fileName)
    }
    createExportMatrix <- function(precursorSet){
      ################################################################################
      ## fragment matrix
      fragmentMatrix      <- dataList$featureMatrix[precursorSet, ]
      dgTMatrix <- as(fragmentMatrix, "dgTMatrix")
      matrixRows <- dgTMatrix@i + 1
      matrixCols <- dgTMatrix@j + 1
      matrixVals <- dgTMatrix@x
      
      numberOfColumns <- ncol(fragmentMatrix)
      numberOfRows <- nrow(fragmentMatrix)
      chunkSize <- 1000
      numberOfChunks <- ceiling(numberOfColumns / chunkSize)
      
      fragmentCounts      <- vector(mode = "integer", length = numberOfColumns)
      fragmentIntensities <- vector(mode = "numeric", length = numberOfColumns)
      fragmentMasses      <- dataList$fragmentMasses
      linesMatrix <- matrix(nrow = numberOfRows, ncol = numberOfChunks)
      
      for(chunkIdx in seq_len(numberOfChunks)){
        colStart <- 1 + (chunkIdx - 1) * chunkSize
        colEnd <- colStart + chunkSize - 1
        if(chunkIdx == numberOfChunks)
          colEnd <- numberOfColumns
        
        numberOfColumnsHere <- colEnd - colStart + 1
        numberOfRowsHere <- max(matrixRows)
        indeces <- matrixCols >= colStart & matrixCols <= colEnd
        
        fragmentMatrixPart <- matrix(data = rep(x = "", times = numberOfRowsHere * numberOfColumnsHere), nrow = numberOfRowsHere, ncol = numberOfColumnsHere)
        fragmentMatrixPart[cbind(matrixRows[indeces], matrixCols[indeces] - colStart + 1)] <- matrixVals[indeces]
        
        fragmentCountsPart      <- apply(X = fragmentMatrixPart, MARGIN = 2, FUN = function(x){ sum(x != "") })
        fragmentIntensitiesPart <- apply(X = fragmentMatrixPart, MARGIN = 2, FUN = function(x){ sum(as.numeric(x), na.rm = TRUE) }) / fragmentCountsPart
        
        linesPart <- apply(X = fragmentMatrixPart, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
        
        fragmentCounts[colStart:colEnd] <- fragmentCountsPart
        fragmentIntensities[colStart:colEnd] <- fragmentIntensitiesPart
        linesMatrix[, chunkIdx] <- linesPart
      }
      
      ## assemble
      linesFragmentMatrixWithHeader <- c(
        paste(fragmentCounts, collapse = "\t"),
        paste(fragmentIntensities, collapse = "\t"),
        paste(fragmentMasses, collapse = "\t"),
        apply(X = linesMatrix, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
      )
      
      ################################################################################
      ## MS1 matrix
      ms1Matrix     <- rbind(
        dataList$dataFrameMS1Header,
        dataList$dataFrameInfos[precursorSet, ]
      )
      ms1Matrix     <- as.matrix(ms1Matrix)
      
      ###########################################################
      ## export annotations
      
      ## process annotations
      annotations <- dataList$annoArrayOfLists
      for(i in 1:length(annotations))
        if(dataList$annoArrayIsArtifact[[i]])
          annotations[[i]] <- c(annotations[[i]], dataList$annotationValueIgnore)
      
      annotationStrings <- vector(mode = "character", length = length(annotations))
      for(i in 1:length(annotations)){
        if(length(annotations[[i]]) > 0)
          annotationStrings[[i]] <- paste(annotations[[i]], collapse = ", ")
        else
          annotationStrings[[i]] <- ""
      }
      annotationStrings <- annotationStrings[precursorSet]
      
      ## process annotaiotn-color-map
      annoPresentAnnotations <- dataList$annoPresentAnnotationsList[-1]
      annoPresentColors      <- dataList$annoPresentColorsList[-1]
      
      if(length(annoPresentAnnotations) > 0){
        annotationColors <- paste(annoPresentAnnotations, annoPresentColors, sep = "=", collapse = ", ")
      } else {
        annotationColors <- ""
      }
      annotationColors <- paste(dataList$annotationColorsName, "={", annotationColors, "}", sep = "")
      
      ## box
      annotationColumn <- c("", annotationColors, dataList$annotationColumnName, annotationStrings)
      
      ms1Matrix[, dataList$annotationColumnIndex] <- annotationColumn
      
      ################################################################################
      ## assemble
      #dataFrame <- cbind(
      #  ms1Matrix,
      #  ms2Matrix
      #)
      linesMS1MatrixWithHeader <- apply(X = ms1Matrix, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
      lines <- paste(linesMS1MatrixWithHeader, linesFragmentMatrixWithHeader, sep = "\t")
      
      return(lines)
    }
    createExportMatrixOld <- function(precursorSet){
      numberOfRows    <- length(precursorSet)
      numberOfColumns <- ncol(dataList$featureMatrix)
      
      ###########################################################
      ## built reduced MS2 matrix
      fragmentMatrix      <- dataList$featureMatrix[precursorSet, ]
      fragmentCounts      <- apply(X = fragmentMatrix, MARGIN = 2, FUN = function(x){ sum(x != 0) })
      fragmentIntensities <- apply(X = fragmentMatrix, MARGIN = 2, FUN = function(x){ sum(x) }) / fragmentCounts
      fragmentMasses      <- dataList$fragmentMasses
      
      fragmentSelection   <- fragmentCounts != 0
      
      fragmentMatrix      <- fragmentMatrix[, fragmentSelection]
      fragmentCounts      <- fragmentCounts[fragmentSelection]
      fragmentIntensities <- fragmentIntensities[fragmentSelection]
      fragmentMasses      <- fragmentMasses[fragmentSelection]
      
      ## fragment matrix
      dgTMatrix <- as(fragmentMatrix, "dgTMatrix")
      matrixRows <- dgTMatrix@i + 1
      matrixCols <- dgTMatrix@j + 1
      matrixVals <- dgTMatrix@x
      
      numberOfColumns2 <- ncol(fragmentMatrix)
      
      fragmentMatrix <- matrix(data = rep(x = "", times = numberOfRows * numberOfColumns2), nrow = numberOfRows, ncol = numberOfColumns2)
      fragmentMatrix[cbind(matrixRows, matrixCols)] <- matrixVals
      
      ## box
      ms2Matrix     <- rbind(
        fragmentCounts,
        fragmentIntensities,
        fragmentMasses,
        fragmentMatrix
      )
      
      ###########################################################
      ## built MS1 matrix
      ms1Matrix     <- rbind(
        dataList$dataFrameMS1Header,
        dataList$dataFrameInfos[precursorSet, ]
      )
      ms1Matrix     <- as.matrix(ms1Matrix)
      
      ###########################################################
      ## export annotations
      
      ## process annotations
      annotations <- dataList$annoArrayOfLists
      for(i in 1:length(annotations))
        if(dataList$annoArrayIsArtifact[[i]])
          annotations[[i]] <- c(annotations[[i]], dataList$annotationValueIgnore)
      
      annotationStrings <- vector(mode = "character", length = length(annotations))
      for(i in 1:length(annotations)){
        if(length(annotations[[i]]) > 0)
          annotationStrings[[i]] <- paste(annotations[[i]], sep = ", ")
        else
          annotationStrings[[i]] <- ""
      }
      annotationStrings <- annotationStrings[precursorSet]
      
      ## process annotaiotn-color-map
      annoPresentAnnotations <- dataList$annoPresentAnnotationsList[-1]
      annoPresentColors      <- dataList$annoPresentColorsList[-1]
      
      if(length(annoPresentAnnotations) > 0){
        annotationColors <- paste(annoPresentAnnotations, annoPresentColors, sep = "=", collapse = ", ")
      } else {
        annotationColors <- ""
      }
      annotationColors <- paste(dataList$annotationColorsName, "={", annotationColors, "}", sep = "")
      
      ## box
      annotationColumn <- c("", annotationColors, dataList$annotationColumnName, annotationStrings)
      
      ms1Matrix[, dataList$annotationColumnIndex] <- annotationColumn
      
      ###########################################################
      ## assemble
      dataFrame <- cbind(
        ms1Matrix,
        ms2Matrix
      )
      
      return(dataFrame)
    }
    writeTable <- function(precursorSet, file){
      #dataFrame <- createExportMatrix(precursorSet)
      #gz1 <- gzfile(description = file, open = "w")
      #write.table(x = dataFrame, file = gz1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      #close(gz1)
      lines <- createExportMatrix(precursorSet)
      gz1 <- gzfile(description = file, open = "w")
      writeLines(text = lines, con = gz1)
      close(gz1)
    }
    ## individual downloads
    output$downloadGlobalMS2filteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterGlobal$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadHcaFilteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterHca$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadPcaFilteredPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterPca$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadSearchPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- filterPca$filter
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadSelectedPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- selectedPrecursorSet
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadAllPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        precursorSet <- 1:dataList$numberOfPrecursors
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    ## download selected
    output$downloadHcaSelectedPrecursors <- downloadHandler(
      filename = function() {
        createExportMatrixName()
      },
      content = function(file) {
        ## get selected precursors
        if(is.null(selectionAnalysisTreeNodeSet)){
          ## all precursors
          precursorSet <- filterHca$filter
        } else {
          precursorSet <- getPrecursorSetFromTreeSelections(clusterDataList = clusterDataList, clusterLabels = selectionAnalysisTreeNodeSet)
        }
        
        writeTable(precursorSet = precursorSet, file = file)
      },
      contentType = 'text/csv'
    )
    output$downloadImportParameterSet <- downloadHandler(
      filename = function() {
        createImportParameterSetExportFileName()
      },
      content = function(file) {
        fileLines <- serializeParameterSetFile(dataList$importParameterSet, toolName, toolVersion)
        writeLines(text = fileLines, con = file)
      },
      contentType = 'text/csv'
    )
    ## download images
    output$downloadHcaImage <- downloadHandler(
      filename = function() {
        fileType <- input$downloadHcaImageType
        createExportImageName("HCA", fileType)
      },
      content = function(file) {
        fileType <- input$downloadHcaImageType
        plotHCA(file, fileType)
      }#,
      #contentType = 'image/png'
    )
    plotHCA <- function(file = NULL, fileType = NULL){
      ## 1 den    ## 2 hea
      ## 3 ms2    ## 4 l anno
      ## 5 l sel  ## 6 l hea
      ## 7 l ms2
      ## 
      ## 1 4
      ## 1 5
      ## 1 6
      ## 2 7
      ## 3 7
      ## 
      
      ## parameters
      widthInInch     <- 10
      heigthInInch    <- 7.5
      resolutionInDPI <- 600
      widthInPixel    <- widthInInch  * resolutionInDPI
      heightInPixel   <- heigthInInch * resolutionInDPI
      
      if(!is.null(file)){
        switch(fileType,
               "png"={
                 png(filename = file, width = widthInPixel, height = heightInPixel, res = resolutionInDPI, bg = "white")
               },
               "svg"={
                 svg(filename = file)
               },
               "pdf"={
                 pdf(file = file, title = "PCA image export from MetFam")
               },
               stop(paste("Unknown file type (", fileType, ")!", sep = ""))
        )
      }
      
      
      graphics::layout(
        mat = matrix(
          data = c(1, 1, 1, 1, 2, 3,
                   4, 5, 6, 7, 8, 8), 
          nrow = 6, ncol = 2), 
        widths = c(4, 1), 
        heights = c(0.6, 1.4, 0.6, 0.6, 0.5, 1.5)
      )
      
      #cex <- par("cex")
      #par(cex = 0.42)
      drawDendrogramPlotImpl()
      #par(cex = cex)
      drawHeatmapPlotImpl() ## out for plotly and adapt layout
      drawMS2PlotImpl()
      
      drawDendrogramLegendImpl()
      drawHeatmapLegendImpl()
      drawMS2LegendImpl()
      drawFragmentDiscriminativityLegendImpl()
      calcPlotAnnoLegendForImage(state$annotationsHca$setOfAnnotations, state$annotationsHca$setOfColors, 15)
      #drawAnnotationLegendImpl()
      
      if(!is.null(file)){
        dev.off()
      }
    }
    output$downloadPcaImage <- downloadHandler(
      filename = function() {
        fileType <- input$downloadPcaImageType
        createExportImageName("PCA", fileType)
      },
      content = function(file) {
        fileType <- input$downloadPcaImageType
        plotPCA(file, fileType)
      }#,
      #contentType = 'image/png'
    )
    plotPCA <- function(file = NULL, fileType = NULL){
      ## 1 score  ## 2 loadings
      ## 3 ms2    ## 4 l anno
      ## 5 l sel  ## 6 l hea
      ## 7 l ms2
      ## 
      ## 1 2 4
      ## 1 2 5
      ## 1 2 6
      ## 1 2 7
      ## 3 3 7
      ## 
      
      ## parameters
      widthInInch     <- 10
      #widthInInch     <- 10 * 4 / 5
      heigthInInch    <- 6
      #heigthInInch    <- 6 * 4 / 5
      resolutionInDPI <- 600
      widthInPixel    <- widthInInch  * resolutionInDPI
      heightInPixel   <- heigthInInch * resolutionInDPI
      
      if(!is.null(file)){
        switch(fileType,
               "png"={
                 png(filename = file, width = widthInPixel, height = heightInPixel, res = resolutionInDPI, bg = "white")
               },
               "svg"={
                 svg(filename = file)
               },
               "pdf"={
                 pdf(file = file, title = "PCA image export from MetFam")
               },
               stop(paste("Unknown file type (", fileType, ")!", sep = ""))
        )
      }
      
      graphics::layout(
        mat = matrix(
          data = c(1, 1, 1, 1, 1, 3,
                   2, 2, 2, 2, 2, 3, 
                   4, 5, 6, 7, 8, 8), 
          nrow = 6, ncol = 3), 
        widths = c(2, 2, 1), 
        heights = c(0.7, 0.6, 0.6, 0.6, 1.2, 1.5)
      )
      
      drawPcaScoresPlotImpl()
      drawPcaLoadingsPlotImpl()
      drawMS2PlotImpl()
      
      calcPlotScoresGroupsLegendForImage(scoresGroups$groups, scoresGroups$colors, 5)
      #calcPlotScoresGroupsLegendForImage(c("Glandular trichomes", "Trichome-free leaves"), scoresGroups$colors, 5)
      drawDendrogramLegendImpl()
      drawMS2LegendImpl()
      drawFragmentDiscriminativityLegendImpl()
      calcPlotAnnoLegendForImage(state$annotationsPca$setOfAnnotations, state$annotationsPca$setOfColors, 20)
      #drawAnnotationLegendImpl()
      
      if(!is.null(file)){
        dev.off()
      }
    }
    output$downloadDistanceMatrix <- downloadHandler(
      filename = function() {
        createExportDistanceMatrixName(currentDistanceMatrixObj$distanceMeasure)
      },
      content = function(file) {
        write.table(x = currentDistanceMatrixObj$distanceMatrix, file = file, sep = "\t", row.names = TRUE, quote = FALSE, col.names=NA)
      },
      contentType = 'text/csv'
    )
    ## download publication data
    output$downloadMsData <- downloadHandler(
      filename = function() {
        return("Metabolite_profile_showcase.txt")
      },
      content = function(file) {
        ## copy data for download
        file.copy(getFile("Metabolite_profile_showcase.txt"), file)
      },
      contentType = "application/zip"
    )
    output$downloadMsMsData <- downloadHandler(
      filename = function() {
        return("MSMS_library_showcase.msp")
      },
      content = function(file) {
        ## copy data for download
        file.copy(getFile("MSMS_library_showcase.msp"), file)
      },
      contentType = "application/zip"
    )
    output$downloadFragmentMatrix <- downloadHandler(
      filename = function() {
        return("Fragment_matrix_showcase.csv")
      },
      content = function(file) {
        ## copy data for download
        file.copy(getFile("Fragment_matrix_showcase.csv"), file)
      },
      contentType = "application/zip"
    )
    output$downloadDocShowcaseProtocol <- downloadHandler(
      filename = function() {
        return("MetFamily_Showcase_protocol.pdf")
      },
      content = function(file) {
        ## copy data for download
        file.copy(getFile("MetFamily_Showcase_protocol.pdf"), file)
      },
      contentType = "application/pdf"
    )
    output$downloadDocUserGuide <- downloadHandler(
      filename = function() {
        return("MetFamily_user_guide.pdf")
      },
      content = function(file) {
        ## copy data for download
        file.copy(getFile("MetFamily_user_guide.pdf"), file)
      },
      contentType = "application/pdf"
    )
    output$downloadDocInputSpecification <- downloadHandler(
      filename = function() {
        return("MetFamily_Input_Specification.pdf")
      },
      content = function(file) {
        ## copy data for download
        file.copy(getFile("MetFamily_Input_Specification.pdf"), file)
      },
      contentType = "application/pdf"
    )
    
    #########################################################################################
    ## report
    output$downloadReport2 <- downloadHandler(
      filename = function() {
        return("MetFamilyReport.pdf")
        #return("MetFamilyReport.html")
      },
      content = function(file) {
        ## TODO 999
        
        ##########################################################################
        ## files
        
        ## source file and tmp file
        tempReportFile <- file.path(tempdir(), "MetFamilyReport.Rmd")
        reportSourceFile <- "/home/htreutle/Code/Java/MetFam/inst/report/Report.Rmd"
        
        ## copy template to tmp dir for reasons of file permissions
        file.copy(from = reportSourceFile, to = tempReportFile, overwrite = TRUE)
        
        ##########################################################################
        ## HCA analyses
        drawHCA = state$showHCAplotPanel
        drawPCA = state$showPCAplotPanel
        
        ##########################################################################
        ## HCA analysis
        if(drawHCA){
          imageFileHCA <- file.path(tempdir(), "HcaTmpFile.png")
          plotHCA(file = imageFileHCA, fileType = "png")
        } else {
          imageFileHCA <- ""
        }
        
        ##########################################################################
        ## PCA analysis
        if(drawPCA){
          imageFilePCA <- file.path(tempdir(), "PcaTmpFile.png")
          plotPCA(file = imageFilePCA, fileType = "png")
        } else {
          imageFilePCA <- ""
        }
        
        # Set up parameters to pass to Rmd document
        params <- list(
          creationTime = date(),
          importParameterSet = dataList$importParameterSet,
          drawHCA = drawHCA,
          drawPCA = drawPCA,
          imageFileHCA = imageFileHCA,
          imageFilePCA = imageFilePCA,
          clusterDataList = clusterDataList,
          pcaDataList = pcaDataList
        )
        
        # Knit the document and eval it in a child of the global environment (this isolates the code in the document from the code in this app)
        rmarkdown::render(input = tempReportFile, output_file = file,
                          output_format = "pdf_document",
                          params = params,
                          envir = new.env(parent = globalenv()),
                          quiet = FALSE
        )
      },
      contentType = "application/pdf"
    )
    
    ## TODO
    #http://127.0.0.1:25805/library/utils/html/zip.html
    serialization <- function(){
      #######################################
      ## enlist
      paramsList <- list(
        ## global MS2 filter
        globalFilter_ms2_masses1          = input$globalFilter_ms2_masses1,
        globalFilter_ms2_masses2          = input$globalFilter_ms2_masses2,
        globalFilter_ms2_masses3          = input$globalFilter_ms2_masses3,
        globalFilter_ms2_ppm              = input$globalFilter_ms2_ppm,
        #input$applyGlobalMS2filters
        ## HCA
        hcaFilterGroupOne                 = input$hcaFilterGroupOne,
        hcaFilterGroupTwo                 = input$hcaFilterGroupTwo,
        hcaFilter_average                 = input$hcaFilter_average,
        hcaFilter_lfc                     = input$hcaFilter_lfc,
        hcaFilterIncludeIgnoredPrecursors = input$hcaFilterIncludeIgnoredPrecursors,
        #input$applyHcaFilters
        hcaDistanceFunction               = input$hcaDistanceFunction,
        #hcaClusterMethod                  = input$hcaClusterMethod,
        hcaClusterMethod                  = "ward.D",
        #input$drawHCAplots
        ## PCA
        pcaGroups                         = input$pcaGroups,
        pcaSamples                        = input$pcaSamples,
        pcaFilter_average                 = input$pcaFilter_average,
        pcaFilter_lfc                     = input$pcaFilter_lfc,
        pcaFilterIncludeIgnoredPrecursors = input$pcaFilterIncludeIgnoredPrecursors,
        #input$applyPcaFilters
        pcaScaling                        = input$pcaScaling,
        pcaLogTransform                   = input$pcaLogTransform,
        pcaDimensionOne                   = input$pcaDimensionOne,
        pcaDimensionTwo                   = input$pcaDimensionTwo,
        #input$drawPCAplots
        ## plot properties
        showPlotControls                  = input$showPlotControls,
        showClusterLabels                 = input$showClusterLabels,
        heatmapContent                    = input$heatmapContent,
        heatmapOrdering                   = input$heatmapOrdering,
        hcaPrecursorLabels                = input$hcaPrecursorLabels,
        showScoresLabels                  = input$showScoresLabels,
        loadingsLabels                    = input$loadingsLabels,
        showLoadingsFeatures              = input$showLoadingsFeatures,
        #showLoadingsFeaturesAnnotated     = input$showLoadingsFeaturesAnnotated,
        #showLoadingsFeaturesUnannotated   = input$showLoadingsFeaturesUnannotated,
        #showLoadingsFeaturesSelected      = input$showLoadingsFeaturesSelected,
        #showLoadingsFeaturesUnselected    = input$showLoadingsFeaturesUnselected,
        showLoadingsAbundance             = input$showLoadingsAbundance,
        #showLoadingsLabels                = "Show labels"    %in% input$pcaLoadingsProperties,
        #showLoadingsAbundance             = "Show abundance" %in% input$pcaLoadingsProperties,
        ## search
        searchMS1orMS2                    = input$searchMS1orMS2,
        searchMS1mass                     = input$searchMS1mass,
        searchMS1massPpm                  = input$searchMS1massPpm,
        search_ms2_masses1                = input$search_ms2_masses1,
        search_ms2_masses2                = input$search_ms2_masses2,
        search_ms2_masses3                = input$search_ms2_masses3,
        searchMS2massPpm                  = input$searchMS2massPpm,
        searchIncludeIgnoredPrecursors    = input$searchIncludeIgnoredPrecursors
      )
      
      #######################################
      ## serialize parameter list
      
      ## built matrix with param names and param values
      tempMatrix    <- matrix(data = c(names(unlist(paramsList)),paste("'", unlist(paramsList), "'", sep = "")), ncol = 2)
      ## quotation of strings
      #textEntries <- is.na(as.numeric(tempMatrix[, 2]))
      #tempMatrix[textEntries, 2] <- paste("'", tempMatrix[textEntries, 2], "'", sep = "")
      ## param name = param value
      paramStrings  <- apply(tempMatrix, MARGIN = 1, FUN = function(x) {paste(x, collapse = "=")})
      ## collapse
      serialization <- paste(paramStrings, collapse = ";")
      return(serialization)
    }
    ## https://github.com/daattali/advanced-shiny/blob/master/update-input/app.R
    deserialization <- function(serialization){
      #######################################
      ## deserialize parameters
      
      #paramStrings <- strsplit(x = strsplit(x = serialization, split = ";")[[1]], split = "=")
      #for(i in 1:length(paramStrings))
      #  paramStrings[[i]] <- paste(paramStrings[[i]][[1]], paste("'", paramStrings[[i]][[2]], "'", sep = ""), sep = "=")
      #paramString <- paste(paramStrings, collapse = ",")
      paramString <- paste(strsplit(x = serialization, split = ";")[[1]], collapse = ",")
      parseText <- paste("paramsList <- list(", paramString, ")", sep = "")
      
      #parseText <- paste("paramsList <- list(", paste(paramStrings, collapse = ","), ")", sep = "")
      eval(parse(text = parseText))
      
      #######################################
      ## update parameter fields
      
      #updateTextInput(session = session, inputId = "globalFilter_ms2_masses1", value = "")
      #updateRadioButtons(session = session, inputId = "hcaFilterGroupOne", choices = dataList$groups, selected = selectedOne)
      #updateCheckboxInput(session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = FALSE)
      #updateSelectInput(session = session, inputId = "presentAnnotationValue", choices = c("[init]"), selected = lalala)
      #updateCheckboxGroupInput(session = session, inputId = "pcaGroups",   choices = dataList$groups, selected = dataList$groups)
      
      ## global MS2 filter
      updateTextInput(         session = session, inputId = "globalFilter_ms2_masses1",          value = paramsList$globalFilter_ms2_masses1)
      updateTextInput(         session = session, inputId = "globalFilter_ms2_masses2",          value = paramsList$globalFilter_ms2_masses2)
      updateTextInput(         session = session, inputId = "globalFilter_ms2_masses3",          value = paramsList$globalFilter_ms2_masses3)
      updateTextInput(         session = session, inputId = "globalFilter_ms2_ppm",              value = paramsList$globalFilter_ms2_ppm)
      #input$applyGlobalMS2filters
      ## HCA
      updateRadioButtons(      session = session, inputId = "hcaFilterGroupOne",                 selected = paramsList$hcaFilterGroupOne)
      updateRadioButtons(      session = session, inputId = "hcaFilterGroupTwo",                 selected = paramsList$hcaFilterGroupTwo)
      updateTextInput(         session = session, inputId = "hcaFilter_average",                 value = paramsList$hcaFilter_average)
      updateTextInput(         session = session, inputId = "hcaFilter_lfc",                     value = paramsList$hcaFilter_lfc)
      updateCheckboxInput(     session = session, inputId = "hcaFilterIncludeIgnoredPrecursors", value = as.logical(paramsList$hcaFilterIncludeIgnoredPrecursors))
      #input$applyHcaFilters
      updateSelectInput(       session = session, inputId = "hcaDistanceFunction",               selected = paramsList$hcaDistanceFunction)
      #updateSelectInput(       session = session, inputId = "hcaClusterMethod",                  selected = paramsList$hcaClusterMethod)
      #input$drawHCAplots
      ## PCA
      updateCheckboxGroupInput(session = session, inputId = "pcaGroups",                         selected = paramsList$pcaGroups)
      updateCheckboxGroupInput(session = session, inputId = "pcaSamples",                        selected = paramsList$pcaSamples)
      updateTextInput(         session = session, inputId = "pcaFilter_average",                 value = paramsList$pcaFilter_average)
      updateTextInput(         session = session, inputId = "pcaFilter_lfc",                     value = paramsList$pcaFilter_lfc)
      updateSelectInput(       session = session, inputId = "pcaFilterIncludeIgnoredPrecursors", selected = paramsList$pcaFilterIncludeIgnoredPrecursors)
      #input$applyPcaFilters
      updateSelectInput(       session = session, inputId = "pcaScaling",                        selected = paramsList$pcaScaling)
      updateCheckboxInput(     session = session, inputId = "pcaLogTransform",                   value = as.logical(paramsList$pcaLogTransform))
      updateSelectInput(       session = session, inputId = "pcaDimensionOne",                   selected = paramsList$pcaDimensionOne)
      updateSelectInput(       session = session, inputId = "pcaDimensionTwo",                   selected = paramsList$pcaDimensionTwo)
      #input$drawPCAplots
      ## plot properties
      updateCheckboxInput(     session = session, inputId = "showPlotControls",                  value = as.logical(paramsList$showPlotControls))
      updateCheckboxInput(     session = session, inputId = "showClusterLabels",                 value = as.logical(paramsList$showClusterLabels))
      updateRadioButtons(      session = session, inputId = "heatmapContent",                    selected = paramsList$heatmapContent)
      updateRadioButtons(      session = session, inputId = "heatmapOrdering",                   selected = paramsList$heatmapOrdering)
      updateRadioButtons(      session = session, inputId = "hcaPrecursorLabels",                selected = paramsList$hcaPrecursorLabels)
      updateCheckboxInput(     session = session, inputId = "showScoresLabels",                  value = as.logical(paramsList$showScoresLabels))
      #updateCheckboxGroupInput(session = session, inputId = "pcaLoadingsProperties",             selected = c(ifelse(as.logical(paramsList$showLoadingsLabels), "Show labels", NULL), ifelse(as.logical(paramsList$showLoadingsAbundance), "Show abundance", NULL)))
      updateRadioButtons(      session = session, inputId = "loadingsLabels",                    selected = paramsList$loadingsLabels)
      updateCheckboxGroupInput(session = session, inputId = "showLoadingsFeatures",              selected = paramsList$showLoadingsFeatures)
      #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesAnnotated",     value = as.logical(paramsList$showLoadingsFeaturesAnnotated))
      #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesUnannotated",   value = as.logical(paramsList$showLoadingsFeaturesUnannotated))
      #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesSelected",      value = as.logical(paramsList$showLoadingsFeaturesSelected))
      #updateCheckboxInput(     session = session, inputId = "showLoadingsFeaturesUnselected",    value = as.logical(paramsList$showLoadingsFeaturesUnselected))
      updateCheckboxInput(     session = session, inputId = "showLoadingsAbundance",             value = as.logical(paramsList$showLoadingsAbundance))
      ## search
      updateRadioButtons(      session = session, inputId = "searchMS1orMS2",                    selected = paramsList$searchMS1orMS2)
      updateTextInput(         session = session, inputId = "searchMS1mass",                     value = paramsList$searchMS1mass)
      updateTextInput(         session = session, inputId = "searchMS1massPpm",                  value = paramsList$searchMS1massPpm)
      updateTextInput(         session = session, inputId = "search_ms2_masses1",                value = paramsList$search_ms2_masses1)
      updateTextInput(         session = session, inputId = "search_ms2_masses2",                value = paramsList$search_ms2_masses2)
      updateTextInput(         session = session, inputId = "search_ms2_masses3",                value = paramsList$search_ms2_masses3)
      updateTextInput(         session = session, inputId = "searchMS2massPpm",                  value = paramsList$searchMS2massPpm)
      updateCheckboxInput(     session = session, inputId = "searchIncludeIgnoredPrecursors",    value = as.logical(paramsList$searchIncludeIgnoredPrecursors))
      #input$applySearch
      
      #######################################
      ## update GUI according to parameters
      
      ###################
      ## global MS2 filter
      filter_ms2_masses1  <- paramsList$globalFilter_ms2_masses1
      filter_ms2_masses2  <- paramsList$globalFilter_ms2_masses2
      filter_ms2_masses3  <- paramsList$globalFilter_ms2_masses3
      filter_ms2_ppm      <- paramsList$globalFilter_ms2_ppm
      
      groupSet        <- dataList$groups
      sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filterBySamples <- TRUE
      filter_average  <- NULL
      filter_lfc      <- NULL
      includeIgnoredPrecursors  <- TRUE
      filter_ms1_masses <- NULL
      filter_ms1_ppm <- NULL
      
      resultObj <- doPerformFiltering(
        groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, 
        includeIgnoredPrecursors
      )
      filterGlobal <<- resultObj$filter
      state$globalMS2filterValid <<- TRUE
      updateGlobalMS2filterInformation()
      
      ###################
      ## HCA filter
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      groupOne        <- paramsList$hcaFilterGroupOne
      groupTwo        <- paramsList$hcaFilterGroupTwo
      filter_average  <- paramsList$hcaFilter_average
      filter_lfc      <- paramsList$hcaFilter_lfc
      includeIgnoredPrecursors  <- paramsList$hcaFilterIncludeIgnoredPrecursors
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      groupSet        <- c(groupOne, groupTwo)
      sampleSet       <- dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
      filterBySamples <- TRUE
      resultObj <- doPerformFiltering(
        groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
      )
      filterHca <<- resultObj$filter
      updateHcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterHca$numberOfPrecursorsFiltered
      checkHcaFilterValidity(numberOfPrecursorsFiltered)
      
      ###################
      ## draw HCA
      # TODO
      
      ###################
      ## PCA filter
      filter_ms2_masses1  <- filterGlobal$filter_ms2_masses1Original   
      filter_ms2_masses2  <- filterGlobal$filter_ms2_masses2Original   
      filter_ms2_masses3  <- filterGlobal$filter_ms2_masses3Original   
      filter_ms2_ppm      <- filterGlobal$filter_ms2_ppmOriginal
      
      groupSet        <- paramsList$pcaGroups
      sampleSet       <- paramsList$pcaSamples
      filterBySamples <- paramsList$filterByPCAgroupSamples
      filter_average  <- paramsList$pcaFilter_average
      filter_lfc      <- paramsList$pcaFilter_lfc
      includeIgnoredPrecursors  <- paramsList$pcaFilterIncludeIgnoredPrecursors
      filter_ms1_masses <- NULL
      filter_ms1_ppm  <- NULL
      
      if(length(groupSet) != 2)
        filter_lfc <- NULL
      
      resultObj <- doPerformFiltering(
        groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
      )
      filterPca <<- resultObj$filter
      updatePcaFilterInformation()
      
      numberOfPrecursorsFiltered <- filterPca$numberOfPrecursorsFiltered
      checkPcaFilterValidity(numberOfPrecursorsFiltered)
      
      ###################
      ## draw PCA
      # TODO
      
      ###################
      ## search
      searchMode <- paramsListsearchMS1orMS2
      if(searchMode == 'MS1 feature m/z'){
        #################################################
        ## get inputs
        filter_ms1_masses <- paramsListsearchMS1mass
        filter_ms1_ppm  <- paramsListsearchMS1massPpm
        
        if(nchar(trimws(filter_ms1_masses)) == 0)
          return()
        
        filter_ms2_masses1  <- NULL
        filter_ms2_masses2  <- NULL
        filter_ms2_masses3  <- NULL
        filter_ms2_ppm      <- NULL
      }
      if(searchMode == 'Fragment m/z'){
        #################################################
        ## get inputs
        filter_ms2_masses1 <- paramsListsearch_ms2_masses1
        filter_ms2_masses2 <- paramsListsearch_ms2_masses2
        filter_ms2_masses3 <- paramsListsearch_ms2_masses3
        filter_ms2_ppm     <- paramsListsearchMS2massPpm
        
        filter_ms1_masses <- NULL
        filter_ms1_ppm <- NULL
      }
      
      filter_lfc      <- NULL
      filter_average  <- NULL
      groupSet        <- dataList$groups
      sampleSet       <- dataList$includedSamples(dataList$groupSampleDataFrame)
      #sampleSet       <- dataList$dataColumnsNameFunctionFromGroupNames(groups = groupSet, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
      filterBySamples <- TRUE
      includeIgnoredPrecursors  <- paramsListsearchIncludeIgnoredPrecursors
      
      #################################################
      ## do filtering
      resultObj <- doPerformFiltering(
        groupSet, sampleSet, filterBySamples, filter_average, filter_lfc, 
        filter_ms2_masses1, filter_ms2_masses2, filter_ms2_masses3, filter_ms2_ppm, 
        filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors
      )
      processSearchFilterResult(resultObj)
    }
  }## function(input, output, session)
)## shinyServer
