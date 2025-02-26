library(testthat)



test_that("Data import with annotations produces the expected numbers", {
  dataList <<- list()
  
  # Set the file paths to test data 
  filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
  fileSpectra        <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")
  fileAnnotation     <- system.file("extdata/testdata/canopus/canopusShort.txt", package = "MetFamily")
  
  # Load the saved parameter set.
  parameterSetPath <- system.file("extdata/testdata/parameterSet.RData", package = "MetFamily")
  load(parameterSetPath)  # Loads 'parameterSet'
  
  # Process data as in the shiny app
  resultObj <- convertToProjectFile(
    filePeakMatrix = filePeakMatrixPath, 
    fileSpectra      = fileSpectra,
    fileAnnotation   = fileAnnotation,
    parameterSet     = parameterSet, 
    progress         = FALSE
  )
  
  lines <- sparseMatrixToString(
    matrixRows = resultObj$matrixRows, 
    matrixCols = resultObj$matrixCols, 
    matrixVals = resultObj$matrixVals, 
    parameterSet = parameterSet
  )
  
  qfeatures <- resultObj$qfeatures
  
  dataList <<- readProjectData(
    fileLines = lines, 
    progress  = FALSE, 
    qfeatures = qfeatures
  )
  
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

test_that("Data import without annotations produces the expected numbers", {
  dataList <<- list()
  
  # Set the file paths to test data 
  filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
  fileSpectra        <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")
  
  # Load the saved parameter set.
  parameterSetPath <- system.file("extdata/testdata/parameterSet.RData", package = "MetFamily")
  load(parameterSetPath)  # Loads 'parameterSet'
  
  # Process data as in the shiny app
  resultObj <- convertToProjectFile(
    filePeakMatrix = filePeakMatrixPath, 
    fileSpectra      = fileSpectra,
    fileAnnotation   =  NULL,  # No annotation
    parameterSet     = parameterSet, 
    progress         = FALSE
  )
  
  lines <- sparseMatrixToString(
    matrixRows = resultObj$matrixRows, 
    matrixCols = resultObj$matrixCols, 
    matrixVals = resultObj$matrixVals, 
    parameterSet = parameterSet
  )
  
  qfeatures <- resultObj$qfeatures
  
  dataList <<- readProjectData(
    fileLines = lines, 
    progress  = FALSE, 
    qfeatures = qfeatures
  )
  
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

##TODO: Expot Import test 