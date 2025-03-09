
filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")

fileSpectra <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")

fileAnnotation <- system.file("extdata/testdata/canopus/canopusShort.txt", package = "MetFamily")

parameterSetPath <- system.file("extdata/testdata/parameterSet.RData",  package = "MetFamily" )
load(parameterSetPath)

resultObj <- convertToProjectFile(filePeakMatrixPath, 
                                  fileSpectra,
                                  parameterSet = parameterSet, 
                                  progress = FALSE)

test_that("Data import produces the expected numbers", {
  dataList <<- list()
  
  # Test the message
  #2460 / 5824 spectra were imported successfully 
  expect_equal(resultObj$numberOfParsedSpectra, 2640)        
  expect_equal(resultObj$numberOfSpectraOriginal, 5824)       
  
  #(15 empty, 3163 low intensity, 5 too heavy)
  expect_equal(resultObj$numberOfSpectraDiscardedDueToNoPeaks, 15) 
  expect_equal(resultObj$numberOfSpectraDiscardedDueToMaxIntensity, 3163)
  expect_equal(resultObj$numberOfSpectraDiscardedDueToTooHeavy, 5)
  
  #2414 / 2640 spectra were successfully mapped to MS¹ features
  expect_equal(resultObj$numberOfPrecursors,2414)  
  
  #34369 / 145973 fragments were successfully imported.
  expect_equal(resultObj$numberOfMS2PeaksAboveThreshold, 34369)
  expect_equal(resultObj$numberOfMS2PeaksOriginal, 145973)
  
  #(1807 too heavy, 109797 low intensity)
  expect_equal(resultObj$numberOfMS2PeaksBelowThreshold,109797 )
  expect_equal(resultObj$numberOfTooHeavyFragments, 1807) 
  
  #2414 / 5823 MS¹ features were successfully imported
  expect_equal(resultObj$numberOfParsedMs1Features, 5823) 
  
  #(420 were isotopes, 2989 without spectra)
  expect_equal(resultObj$numberOfRemovedPrecursorIsotopePeaks, 420)
  expect_equal(resultObj$numberOfUnmappedPrecursors, 2989) 
  
})

lines <- sparseMatrixToString(matrixRows = resultObj$matrixRows, matrixCols = resultObj$matrixCols,
                              matrixVals = resultObj$matrixVals, parameterSet = parameterSet)

dataList0 <- readProjectData(fileLines = lines, progress = FALSE)

dataList <- add_qfeatures(dataList0, qfeatures = resultObj$qfeatures, fileAnnotation)

filterObjPath <- system.file("extdata/testdata/canopus/filterPca_canopus.rds", package = "MetFamily")
filterObj <- readRDS(filterObjPath)

test_that("Filtering functions work with annotation data", {

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
  expect_equal(rownames(dataList$dataFrameInfos)[32], "  215.127 /   8.30")
  expect_equal(filteredResult$numberOfPrecursorsFiltered, 209)
})

test_that("Filtering functions work without annotation data", {

  # Test filtering using filterObj
  filteredResult <- filterData(
    dataList = dataList0,
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
  expect_equal(rownames(dataList$dataFrameInfos)[32], "  215.127 /   8.30")

  
  
})


test_that("Example data properly integrates with filtering functions", {
  
  example_file_path <- system.file("extdata/showcase/Project_file_showcase_annotated.csv.gz", 
                                   package = "MetFamily")

  dataList <- readClusterDataFromProjectFile(
    file = example_file_path,
    progress = FALSE
  )
  
  load(filterObjPath)
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
  
  expect_true(filteredResult$numberOfPrecursorsFiltered > 0)
  expect_true(filteredResult$numberOfPrecursorsFiltered, 209)
})



