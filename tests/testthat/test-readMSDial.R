library(purrr)

filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
qf <- readMSDial(filePeakMatrix)

test_that("exampledata", {
  data <- parsePeakAbundanceMatrixQF(qf, doPrecursorDeisotoping=TRUE,
                                   mzDeviationInPPM_precursorDeisotoping=10, mzDeviationAbsolute_precursorDeisotoping=0.01,
                                   maximumRtDifference=0.05,
                                   progress=FALSE)
  
  
  ## Test dimensions  
  expect_equal(nrow(data$dataFrame), 5403)
  expect_equal(ncol(data$dataFrame),   20)
  
  ## Test some values
  expect_true(all(summary(t(data$dataFrame[1,c("TRI03", "TRI02", "TRI01", "LVS03", "LVS02", "LVS01")])) == c("Min.   : 236.0  ", "1st Qu.: 306.2  ", "Median : 357.5  ", "Mean   : 501.0  ", "3rd Qu.: 554.2  ", "Max.   :1146.0  ")))
  expect_true(all(round(summary(data$dataFrame[, "TRI01"])) == c(0,     809,    1537,    9818,    3207, 4407926)))
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



# gp: testing different MS-Dial data formats

files <- system.file(package = "MetFamily",
                     c("extdata/showcase/Metabolite_profile_showcase.txt", 
                       "extdata/testdata/Height_MSDIAL.ver5.2.240424-short.txt",
                       "extdata/testdata/Height_202502_wide-msms.txt"))

rowDataDefaultNames2020 <- c("Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", 
                             "Adduct ion name", "Fill %", "MS/MS included", "INCHIKEY", "SMILES", 
                             "LINK", "Dot product", "Reverse dot product", "Fragment presence %", 
                             "Spectrum reference file name")

colDataDefaultNames2020 <- c("Class", "Type", "Injection order")


check_qf <- function(q1) {
  # q1 <- ms_reads[[1]]
  
  rD <- rowData(q1)[[1]]
  cD <- colData(q1)
  
  c(identical(names(rD), rowDataDefaultNames2020),
    identical(names(cD), colDataDefaultNames2020),
    identical(dim(assay(q1)), c(nrow(rD), nrow(cD))))
}


test_that("Different MS-Dial formats can be read in properly", {
  
  ms_reads <- map(files, readMSDial)
  
  t3x3 <- rep(list(rep(TRUE, 3)), 3)
  
  expect_equal(map(ms_reads, check_qf), t3x3)
      
})
  