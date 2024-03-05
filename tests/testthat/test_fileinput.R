test_that("trivial", {
    expect_equal(1, 1)
})

test_that("exampledata", {
  filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
  data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping=TRUE,
                                   mzDeviationInPPM_precursorDeisotoping=10, mzDeviationAbsolute_precursorDeisotoping=0.01,
                                   maximumRtDifference=0.05,
                                   progress=FALSE)

  ## Test dimensions  
  expect_equal(nrow(data$dataFrame1), 5403)
  expect_equal(ncol(data$dataFrame1),   20)

  ## Test some values
  expect_true(all(summary(t(data$dataFrame1[1,c("TRI03", "TRI02", "TRI01", "LVS03", "LVS02", "LVS01")])) == c("Min.   : 236.0  ", "1st Qu.: 306.2  ", "Median : 357.5  ", "Mean   : 501.0  ", "3rd Qu.: 554.2  ", "Max.   :1146.0  ")))
  expect_true(all(round(summary(data$dataFrame1[, "TRI01"])) == c(0,     809,    1537,    9818,    3207, 4407926)))
})

test_that("MS-Dial 4.X", {
#   filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
#     data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping, 
# 	mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, 
# 	maximumRtDifference, 
# 	progress=FALSE)
# 
#     expect_equal(ncol(data$dataFrame), 41)
})

test_that("MSP reading works", {
  
  if(FALSE) {
  fileSpectra <- system.file("extdata/showcase/Metabolite_profile_showcase.msp", package = "MetFamily")
  #load(system.file("extdata/testdata/readMSPreturnObj.Rdata", package = "MetFamily"))
  
  returnObj <- parseMSP(fileSpectra = fileSpectra, 
                        minimumIntensityOfMaximalMS2peak = 2000, 
                        minimumProportionOfMS2peaks = 0.05, 
                        neutralLossesPrecursorToFragments = TRUE,
                        neutralLossesFragmentsToFragments = FALSE,
                        progress = FALSE)

  
  expect_equal(returnObj$numberOfSpectra, 2640)
  expect_equal(returnObj$numberOfSpectraOriginal, 5824)
  expect_equal(returnObj$numberOfMS2PeaksOriginal, 145973)
  expect_equal(returnObj$numberOfMS2PeaksWithNeutralLosses, 68738)
  expect_equal(returnObj$numberOfMS2PeaksAboveThreshold, 34369)
  expect_equal(returnObj$numberOfMS2PeaksBelowThreshold, 109797)
  expect_equal(returnObj$numberOfTooHeavyFragments, 1807)
  expect_equal(returnObj$numberOfSpectraDiscardedDueToNoPeaks, 15)
  expect_equal(returnObj$numberOfSpectraDiscardedDueToMaxIntensity, 3163)
  expect_equal(returnObj$numberOfSpectraDiscardedDueToTooHeavy, 5)
  expect_equal(length(returnObj$precursorMz), 2640)
  expect_equal(returnObj$precursorMz[1], 85)
  expect_equal(length(returnObj$precursorRt), 2640)
  expect_equal(returnObj$precursorRt[1], 10.98)
  
  }
  
  # > str(returnObj, max.level = 1)
  # List of 14
  # $ fileSpectra                              : chr "/tmp/RtmpL82I4E/ff16bb92ac5e87cb202eec2e/0.msp"
  # $ spectraList                              :List of 2640
  # $ numberOfSpectra                          : int 2640
  # $ numberOfSpectraOriginal                  : int 5824
  # $ numberOfMS2PeaksOriginal                 : num 145973
  # $ numberOfMS2PeaksWithNeutralLosses        : num 68738
  # $ numberOfMS2PeaksAboveThreshold           : num 34369
  # $ numberOfMS2PeaksBelowThreshold           : num 109797
  # $ numberOfTooHeavyFragments                : num 1807
  # $ numberOfSpectraDiscardedDueToNoPeaks     : num 15
  # $ numberOfSpectraDiscardedDueToMaxIntensity: num 3163
  # $ numberOfSpectraDiscardedDueToTooHeavy    : num 5
  # $ precursorMz                              : num [1:2640] 85 85 85 85 85 ...
  # $ precursorRt                              : num [1:2640] 10.98 7.42 12.4 16.26 13 ...
  # 
  
  
  
  ##
  ## returnObj contains the following elements:
  ##
  
  # returnObj <- list()
  # returnObj$fileSpectra <- fileSpectra
  # returnObj$spectraList <- list()
  # returnObj$numberOfSpectra <- 0
  # returnObj$numberOfMS2PeaksOriginal <- 0
  # returnObj$numberOfMS2PeaksWithNeutralLosses <- 0
  # returnObj$numberOfMS2PeaksAboveThreshold <- 0
  # returnObj$numberOfMS2PeaksBelowThreshold <- 0
  # returnObj$precursorMz <- vector(mode = "numeric")
  # returnObj$precursorRt <- vector(mode = "numeric")
  
})
