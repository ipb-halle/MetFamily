
if (interactive()) {
  library(MetFamily)
  library(testthat)
}


filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")

fileSpectra <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")

fileAnnotation <- system.file("extdata/testdata/canopus/canopus1680.txt", package = "MetFamily")

# system.file("extdata/testdata/parameterSet.RData",  package = "MetFamily" )
parameterSet <- parameterSetDefault()
parameterSet$minimumIntensityOfMaximalMS2peak <- 2000
parameterSet$minimumProportionOfMS2peaks <- 0.05
parameterSet$minimumNumbersOfFragments <- 2
parameterSet$minimumAbsoluteMS2peaks <- 500

resultObj <- convertToProjectFile(filePeakMatrixPath, 
                                  fileSpectra,
                                  parameterSet = parameterSet, 
                                  progress = FALSE)

# Note: Care, changing default import parameters will affect these values
test_that("Data import produces the expected numbers", {
  
  # Test the message
  #2688 / 5823 spectra were imported successfully 
  expect_equal(resultObj$numberOfParsedSpectra, 2688)        
  expect_equal(resultObj$numberOfSpectraOriginal, 5823)       

  #(63 too few peaks, 3072 low intensity, 0 too heavy)
  expect_equal(resultObj$numberOfSpectraDiscardedDueToNoPeaks, 63) 
  expect_equal(resultObj$numberOfSpectraDiscardedDueToMaxIntensity, 3072)
  expect_equal(resultObj$numberOfSpectraDiscardedDueToTooHeavy, 0)
  
  #2463 / 2640 spectra were successfully mapped to MS¹ features
  expect_equal(resultObj$numberOfPrecursors, 2463)  
  
  #34369 / 145973 fragments were successfully imported.
  expect_equal(resultObj$numberOfMS2PeaksAboveThreshold, 48307)
  expect_equal(resultObj$numberOfMS2PeaksOriginal, 145973)
  
  #(44 too heavy, 97622 low intensity)
  expect_equal(resultObj$numberOfMS2PeaksBelowThreshold, 97622)
  expect_equal(resultObj$numberOfTooHeavyFragments, 44)
  
  #2414 / 5823 MS¹ features were successfully imported
  expect_equal(resultObj$numberOfParsedMs1Features, 5823) 
  
  #(420 were isotopes, 2940 without spectra)
  expect_equal(resultObj$numberOfRemovedPrecursorIsotopePeaks, 420)
  expect_equal(resultObj$numberOfUnmappedPrecursors, 2940) 
  
})


lines <- sparseMatrixToString(matrixRows = resultObj$matrixRows, matrixCols = resultObj$matrixCols,
                              matrixVals = resultObj$matrixVals, parameterSet = parameterSet)

dataList0 <- readProjectData(fileLines = lines, progress = FALSE)

dataList <- add_qfeatures(dataList = dataList0, qfeatures = resultObj$qfeatures, fileAnnotation)

# old object doesn't work well anymore
# filterObjPath <- system.file("extdata/testdata/canopus/filterPca_canopus.rds", package = "MetFamily")
# filterObj <- readRDS(filterObjPath)

filterObj <- makeFilterObj(dataList, filter_average = 10000, filter_lfc = 2)

test_that("Filtering functions work with annotation data", {

  # Test filtering using filterObj
  filteredResult <- filterData(
    dataList = dataList,
    sampleClasses = filterObj$groupSet,
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
  expect_equal(rownames(dataList$dataFrameInfos)[32], "  231.137 / 14.34")
  expect_equal(filteredResult$numberOfPrecursorsFiltered, 214)
})

test_that("Filtering functions work without annotation data", {

  # Test filtering using filterObj
  filteredResult <- filterData(
    dataList = dataList0,
    sampleClasses = filterObj$groupSet,
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
  expect_equal(filteredResult$numberOfPrecursorsFiltered, 214)
  expect_equal(rownames(dataList$dataFrameInfos)[32], "  231.137 / 14.34")

  
  
})


test_that("Example data properly integrates with filtering functions", {
  
  example_file_path <- system.file("extdata/showcase/Project_file_showcase_annotated.csv.gz", 
                                   package = "MetFamily")

  dataList <- readClusterDataFromProjectFile(
    file = example_file_path,
    progress = FALSE
  )
  
  filteredResult <- filterData(
    dataList = dataList,
    sampleClasses = filterObj$groupSet,
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
  # there is a known difference between loading the project file (219) and the
  # example data from separate files (209)
  expect_equal(filteredResult$numberOfPrecursorsFiltered, 219)
})
