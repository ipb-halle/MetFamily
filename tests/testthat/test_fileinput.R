test_that("trivial", {
    expect_equal(1, 1)
})

test_that("exampledata", {
#   filePeakMatrix <- system.file("data/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
#     data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping, 
# 	mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, 
# 	maximumRtDifference, 
# 	progress=FALSE)
# 
#     expect_equal(nrow(data), 17)
#     expect_equal(ncol(data), 34)
# 
})

test_that("MS-Dial 4.X", {
#   filePeakMatrix <- system.file("data/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
#     data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping, 
# 	mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, 
# 	maximumRtDifference, 
# 	progress=FALSE)
# 
#     expect_equal(ncol(data$dataFrame), 41)
})

