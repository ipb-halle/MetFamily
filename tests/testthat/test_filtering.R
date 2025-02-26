library(testthat)
library(MetFamily)

test_that("Filtering functions work with annotation data", {
  # Set up test data with annotation
  filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
  fileSpectra <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")
  fileAnnotation <- system.file("extdata/testdata/canopus/canopusShort.txt", package = "MetFamily")
  
  # Load parameter set and filter object
  parameterSetPath <- system.file("extdata/testdata/parameterSet.RData", package = "MetFamily")
  filterObjPath <- system.file("extdata/testdata/filterObj.Rdata", package = "MetFamily")
  load(parameterSetPath)
  load(filterObjPath)

  # Convert files to project with annotation
  resultObj <- convertToProjectFile(
    filePeakMatrix = filePeakMatrixPath,
    fileSpectra = fileSpectra, 
    fileAnnotation = fileAnnotation,
    parameterSet = parameterSet,
    progress = FALSE
  )
  
  # Convert matrix to string and read project data
  lines <- sparseMatrixToString(
    matrixRows = resultObj$matrixRows,
    matrixCols = resultObj$matrixCols,
    matrixVals = resultObj$matrixVals,
    parameterSet = parameterSet
  )
  dataList <- readProjectData(
    fileLines = lines,
    progress = FALSE,
    qfeatures = resultObj$qfeatures
  )

  # Test filtering using filterObj
  filteredResult <- filterDataX(
    dataList = dataList,
    grouXXXps = filterObj$groupSet,
    sampleSet = filterObj$sampleSet,
    filterBySamples = filterObj$filterBySamples,
    filter_average = filterObj$filter_average,
    filter_lfc = filterObj$filter_lfc,
    filterList_ms2_masses = filterObj$filter_ms2_masses,
    filter_ms2_ppm = filterObj$filter_ms2_ppm,
    filter_ms1_masses = filterObj$filter_ms1_masses,
    filter_ms1_ppm = filterObj$filter_ms1_ppm,
    includeIgnoredPrecursors = filterObj$includeIgnoredPrecursors,
    progress = FALSE
  )

  # Test that filtering worked
  expect_equal(filteredResult$numberOfPrecursorsFiltered, 209)
})

test_that("Filtering functions work without annotation data", {
  # Set up test data without annotation
  filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
  fileSpectra <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")
  
  # Load parameter set and filter object
  parameterSetPath <- system.file("extdata/testdata/parameterSet.RData", package = "MetFamily")
  filterObjPath <- system.file("extdata/testdata/filterObj.Rdata", package = "MetFamily")
  load(parameterSetPath)
  load(filterObjPath)

  # Convert files to project without annotation
  resultObj <- convertToProjectFile(
    filePeakMatrix = filePeakMatrixPath,
    fileSpectra = fileSpectra,
    fileAnnotation = NULL,  # No annotation
    parameterSet = parameterSet,
    progress = FALSE
  )
  
  # Convert matrix to string and read project data
  lines <- sparseMatrixToString(
    matrixRows = resultObj$matrixRows,
    matrixCols = resultObj$matrixCols,
    matrixVals = resultObj$matrixVals,
    parameterSet = parameterSet
  )
  dataList <- readProjectData(
    fileLines = lines,
    progress = FALSE,
    qfeatures = resultObj$qfeatures
  )

  # Test filtering using filterObj
  filteredResult <- filterDataX(
    dataList = dataList,
    grouXXXps = filterObj$groupSet,
    sampleSet = filterObj$sampleSet,
    filterBySamples = filterObj$filterBySamples,
    filter_average = filterObj$filter_average,
    filter_lfc = filterObj$filter_lfc,
    filterList_ms2_masses = filterObj$filter_ms2_masses,
    filter_ms2_ppm = filterObj$filter_ms2_ppm,
    filter_ms1_masses = filterObj$filter_ms1_masses,
    filter_ms1_ppm = filterObj$filter_ms1_ppm,
    includeIgnoredPrecursors = filterObj$includeIgnoredPrecursors,
    progress = FALSE
  )

  # Test that filtering worked as expected
  expect_equal(filteredResult$numberOfPrecursorsFiltered, 209)
})



#####debug custom function 

filterDataX <- function(dataList, grouXXXps, sampleSet, filterBySamples, filter_average, filter_lfc, filterList_ms2_masses, filter_ms2_ppm, filter_ms1_masses, filter_ms1_ppm, includeIgnoredPrecursors, progress = FALSE){
  
  ##########################################
  ## filter
  filter <- rep(x = TRUE, times = dataList$numberOfPrecursors)
   print(paste("Initial features:", sum(filter)))
  
  ## filter_average
  if(!is.null(filter_average)){
    filter <- filter & apply(
      X = as.data.frame(dataList$dataFrameMeasurements[, sapply(as.vector(grouXXXps), FUN = dataList$dataMeanColumnNameFunctionFromName)]), 
      MARGIN = 1, 
      FUN = mean
    ) >= filter_average
     print(paste("After average filter:", sum(filter)))
  }
  
  ## filter_lfc
  if(!is.null(filter_lfc)){
    if(filter_lfc != 0){
      if(length(grouXXXps) != 2){  stop("The number of grouXXXps for LFC is not equal to two!") }
      if(filter_lfc > 0)
        filter <- filter & dataList$dataFrameMeasurements[, dataList$lfcColumnNameFunctionFromName(grouXXXps[[1]], grouXXXps[[2]])] >= filter_lfc
      else
        filter <- filter & dataList$dataFrameMeasurements[, dataList$lfcColumnNameFunctionFromName(grouXXXps[[1]], grouXXXps[[2]])] <= filter_lfc
       print(paste("After LFC filter:", sum(filter)))
    }
  }
  
  ## filter_ms2_masses, filter_ms2_ppm
  if(!is.null(filterList_ms2_masses) & !is.null(filter_ms2_ppm) & length(filterList_ms2_masses) > 0){
    error <- abs(dataList$fragmentMasses) * filter_ms2_ppm / 1E6
    filterFragmentLists <- rep(x = FALSE, times = dataList$numberOfPrecursors)
    for(filter_ms2_masses in filterList_ms2_masses){
      filterTempFragmentList <- rep(x = TRUE, times = dataList$numberOfPrecursors)
      for(fragmentIndex in seq_len(length(filter_ms2_masses))){
        distances <- abs(dataList$fragmentMasses - filter_ms2_masses[[fragmentIndex]])
        filterTempFragment <- rep(x = FALSE, times = dataList$numberOfPrecursors)
        columns <- which(distances <= error)
        for(column in columns)
          filterTempFragment <- filterTempFragment | dataList$featureMatrix[, column] != 0
        filterTempFragmentList <- filterTempFragmentList & filterTempFragment
      }
      filterFragmentLists <- filterFragmentLists | filterTempFragmentList
    }
    filter <- filter & filterFragmentLists
    print(paste("After MS2 filter:", sum(filter)))
  }
  
  ## filter_ms1_masses, filter_ms1_ppm
  if(!is.null(filter_ms1_masses) & !is.null(filter_ms1_ppm) & length(filter_ms1_masses) > 0){
    precursorMasses <- as.numeric(dataList$dataFrameInfos[["m/z"]])
    error <- precursorMasses * filter_ms1_ppm / 1E6
    
    filterMS1masses <- rep(x = FALSE, times = dataList$numberOfPrecursors)
    for(precursorMassIndex in seq_len(length(filter_ms1_masses))){
      distances <- abs(precursorMasses - filter_ms1_masses[[precursorMassIndex]])
      filterPart <- distances <= error
      filterPart[is.na(filterPart)] <- FALSE
      filterMS1masses <- filterMS1masses | filterPart
    }
    filter <- filter & filterMS1masses
     print(paste("After MS1 filter:", sum(filter)))
  }
  
  ## include ignored precursors
  if(!includeIgnoredPrecursors)
    filter <- filter & !dataList$annoArrayIsArtifact
   print(paste("After removing ignored precursors:", sum(filter)))
  
  filter <- which(filter)
  
  resultObj <- list()
  resultObj$filter <- filter
  resultObj$numberOfPrecursors <- dataList$numberOfPrecursors
  resultObj$numberOfPrecursorsFiltered <- length(filter)
  if(is.null(grouXXXps)){
    resultObj$grouXXXps    <- list()
    resultObj$sampleSet <- list()
    resultObj$filterBySamples <- NA
  } else {
    resultObj$grouXXXps    <- grouXXXps
    resultObj$sampleSet <- sampleSet
    resultObj$filterBySamples <- filterBySamples
  }
  resultObj$filter_average           <- ifelse(is.null(filter_average), 0, filter_average)
  resultObj$filter_lfc               <- ifelse(is.null(filter_lfc), 0, filter_lfc)
  resultObj$filterList_ms2_masses    <- if(is.null(filterList_ms2_masses)) list() else filterList_ms2_masses
  resultObj$filter_ms2_ppm           <- ifelse(is.null(filter_ms2_ppm), "", filter_ms2_ppm)
  resultObj$filter_ms1_masses        <- if(is.null(filter_ms1_masses)) list() else filter_ms1_masses
  resultObj$filter_ms1_ppm           <- ifelse(is.null(filter_ms1_ppm), "", filter_ms1_ppm)
  resultObj$includeIgnoredPrecursors <- ifelse(is.null(includeIgnoredPrecursors), NA, includeIgnoredPrecursors)
  
  return (resultObj)
}
