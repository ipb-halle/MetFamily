filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
qf <- readMSDial(filePeakMatrix)

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


