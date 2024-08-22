filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
qf <- readMSDial(filePeakMatrix)

test_that("exampledata", {
  data <- parsePeakAbundanceMatrixQF(qf, doPrecursorDeisotoping=TRUE,
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


test_that("Sum of Intensities is correct", {
  
  sumQf <- sum(colSums(assay(qf)))
  expect_equal(sumQf, 232301678)

})

test_that("Number of Rows and Columns are correct", {
  
  # check if number of rows is identical
  nrowQf <- nrow(assay(qf))
  expect_equal(nrowQf, 5823L)
  
  # check if number of cols is identical
  ncolQf <- ncol(assay(qf))+ ncols(rowData(qf))
  expect_equal(as.integer(ncolQf), 20)
  
})


