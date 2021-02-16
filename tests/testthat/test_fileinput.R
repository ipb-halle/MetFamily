test_that("trivial", {
    expect_equal(1, 1)
})

test_test("exampledata", {

    data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping, 
	mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, 
	maximumRtDifference, 
	progress=FALSE)

    expect_equal(nrow(data), 17)
    expect_equal(ncol(data), 34)

})

test_test("MS-Dial 4.X", {

    data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping, 
	mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, 
	maximumRtDifference, 
	progress=FALSE)

    expect_equal(nrow(data), 17)
})


