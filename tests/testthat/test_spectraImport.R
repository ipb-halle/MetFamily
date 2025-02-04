library(testthat)

dataList <<- list()

# Set the file paths to test data 
filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
fileSpectra        <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")
fileAnnotation     <- system.file("extdata/testdata/canopus/canopusShort.txt", package = "MetFamily")

# Load a saved parameter set.
parameterSetPath <- system.file("extdata/testdata/parameterSet.RData", package = "MetFamily")
load(parameterSetPath)  # Loads 'parameterSet'

# Process data as in the app
# Convert the raw files to a project file 
resultObj <- convertToProjectFile(
  filePeakMatrix = filePeakMatrixPath, 
  fileSpectra      = fileSpectra,
  fileAnnotation   = fileAnnotation,
  parameterSet     = parameterSet, 
  progress         = FALSE
)

# Convert the sparse matrix to a string
lines <- sparseMatrixToString(
  matrixRows = resultObj$matrixRows, 
  matrixCols = resultObj$matrixCols, 
  matrixVals = resultObj$matrixVals, 
  parameterSet = parameterSet
)
qfeatures <- resultObj$qfeatures

# Read the project data
dataList <<- readProjectData(
  fileLines = lines, 
  progress  = FALSE, 
  qfeatures = qfeatures
)

test_that("Data import produces the expected numbers", {
  expect_equal(resultObj$numberOfParsedSpectra, 2640)        
  expect_equal(resultObj$numberOfSpectraOriginal, 5824) #2460 / 5824 spectra were imported successfully        
  expect_equal(resultObj$numberOfSpectraDiscardedDueToNoPeaks, 15) 
  #expect_equal(resultObj$numberOfSpectraDiscardedDueToMaxIntensity, 5)  
  expect_equal(resultObj$numberOfSpectraDiscardedDueToTooHeavy, 5)     


  #Check that the final dataList has been populated correctly.
  expect_true(is.list(dataList))
  expect_true("minimumMass" %in% names(dataList))
})
