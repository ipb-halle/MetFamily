
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

## set log path: broken
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

#if(!isDevelopment)  setwd("/var/log/shiny-server")

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
    metFamilyBuilt <- "1.2.0"
    
    ## annotation
    artifactName <- "Ignore"
    artifactColor <- "red"
    selectionNone <- "None"
    
    ## GUI
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
      #importedOrLoadedFile_s_ = NULL, 
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
      
      ## panels
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
      scoresGroupsLegendHeight = -1,
      ## metabolite families
      #metaboliteFamilySelected = FALSE,
      #classifierLoaded = FALSE,
      #classifierClassSelected = FALSE,
      #classifierClassMS1featureSelected = FALSE,
      putativeAnnotationsTableFromAnalysisRowSelected = FALSE
    )
    plotToShow  <- "Display HCA"
    plotsToShow <- "Display HCA"
    showSideBar <- TRUE
    
    ##############################################
    ## MS2 peaks
    
    ## resultObj$fragmentMasses <- fragmentsX
    ## resultObj$fragmentAbundances <- fragmentsY
    ## resultObj$fragmentColor <- fragmentsColor
    ## resultObj$fragmentDiscriminativity <- fragmentDiscriminativity
    ## ?? resultObj$clusterDiscriminativity <- clusterDiscriminativity
    ## ? resultObj$infoText <- infoText
    ## ? resultObj$metFragLinkList <- NULL
    ## ? resultObj$precursorSet <- precursorSet
    ## ? resultObj$numberOfPrecursors <- numberOfPrecursors
    ms2PlotValues <- reactiveValues(
      fragmentListClicked = NULL,
      fragmentListHovered = NULL
    )
    
    ## statistics for dendrogram node or pca loadings brush
    dendrogramFragmentStatistics <- FALSE
    
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
    ## TODO decompose
    resetWorkspaceFunctions <- list()
    resetWorkspace <- function(){
      print(paste("resetWorkspace"))
      
      for(resetWorkspaceFunction in resetWorkspaceFunctions)
        resetWorkspaceFunction()
      
      
      #########################################################################################
      ## reset
      
      ## reset plots
      doClearPlots()
      
      ## reset variables
      clusterDataList <<- NULL
      pcaDataList <<- NULL
      columnsOfInterestForHeatmap <<- NULL
      
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
      ms2PlotValues$fragmentListHovered <<- NULL
      ms2PlotValues$fragmentListClicked <<- NULL
      dendrogramFragmentStatistics <<- FALSE
      
      ## reset state
      #state$importedOrLoadedFile_s_ <<- NULL
      state$analysisType <<- "HCA"
      
      ## panels
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
      #state$metaboliteFamilySelected <<- FALSE
      #state$classifierLoaded <<- FALSE
      #state$classifierClassSelected <<- FALSE
      #state$classifierClassMS1featureSelected <<- FALSE
      state$putativeAnnotationsTableFromAnalysisRowSelected <<- FALSE
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
      
      min <- min(dataList$ms2_masses)
      max <- max(dataList$ms2_masses)
      
      fragmentPlotRange$xMin <<- min
      fragmentPlotRange$xMax <<- max
      fragmentPlotRange$xInterval <<- c(min, max)
      fragmentPlotRange$xIntervalSize <<- max - min
      
      output$fragmentPlot <- renderPlot({
        print(paste("### MS2 ### all init"))
        plotFragmentsFromDataList(dataList = dataList, xInterval = fragmentPlotRange$xInterval)
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
      putativeAnnotationsTableFromAnalysisCurrentlyShown <<- NULL
      
      ## annotation
      allAnnotationNames <<- NULL
      
      ## set state
      sampleTable <<- createSampleTable()
      setSampleTable()
      updateAnnotationOverview()
    }
    
    #########################################################################################
    ## source all server stuff
    suspendOnExitFunctions <- list()
    
    source(file = "server_functionsFilters.R", local = TRUE)$value
    source(file = "server_functionsSelections.R", local = TRUE)$value
    source(file = "server_functionsTableGui.R", local = TRUE)$value
    source(file = "server_functionsDownloads.R", local = TRUE)$value
    source(file = "server_functionsSerialization.R", local = TRUE)$value
    source(file = "server_guiPlots.R", local = TRUE)$value
    source(file = "server_guiAnnotation.R", local = TRUE)$value
    source(file = "server_guiTabInput.R", local = TRUE)$value
    source(file = "server_guiTabAnnotation.R", local = TRUE)$value
    source(file = "server_guiTabClassifier.R", local = TRUE)$value
    source(file = "server_guiTabSampleFilter.R", local = TRUE)$value
    source(file = "server_guiTabMsmsFilter.R", local = TRUE)$value
    source(file = "server_guiTabPca.R", local = TRUE)$value
    source(file = "server_guiTabHca.R", local = TRUE)$value
    source(file = "server_guiTabSearch.R", local = TRUE)$value
    source(file = "server_guiTabExport.R", local = TRUE)$value
    source(file = "server_guiPlotControls.R", local = TRUE)$value
    source(file = "server_guiMs2plot.R", local = TRUE)$value
    ## ui generation
    source(file = "ui_rightColumn.R", local = TRUE)
    
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
      showModal(session = session, ui = modalDialog(title = title, HTML(msg)))
    }
    showUiDialog <- function(modalDialog){
      print("Show ui dialog")
      showModal(session = session, ui = modalDialog)
    }
    
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
      
      print("Exiting")
      stopApp()
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
    output$metaboliteFamilySelected <- reactive({
      print(paste("reactive update metaboliteFamilySelected", state_tabAnnotation$metaboliteFamilySelected))
      return(state_tabAnnotation$metaboliteFamilySelected)
    })
    output$showPutativeAnnotationsTableFromAnalysis <- reactive({
      print(paste("reactive update showPutativeAnnotationsTableFromAnalysis", state$showPutativeAnnotationsTableFromAnalysis))
      return(state$showPutativeAnnotationsTableFromAnalysis)
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
    outputOptions(output, 'metaboliteFamilySelected',suspendWhenHidden=FALSE)
    outputOptions(output, 'showPutativeAnnotationsTableFromAnalysis',  suspendWhenHidden=FALSE)
    #outputOptions(output, 'putativeAnnotationsTableFromAnalysisRowSelected',  suspendWhenHidden=FALSE)
  }## function(input, output, session)
)## shinyServer
