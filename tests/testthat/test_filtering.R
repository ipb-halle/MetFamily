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
  filteredResult <- filterData(
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
  filteredResult <- filterData(
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
