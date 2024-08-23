test_that("metaboliteProfileParsing works", {
  load(system.file("extdata/testdata/processMS1data.Rdata", package = "MetFamily"))

  d <- dataColumnsNameFunctionFromGroupIndex(groupIdx = 2, sampleNamesToExclude = NA)
  
  ## Test dimensions  
  expect_equal(length(d), 3)
  
  result <- processMS1data(sampleNamesToExclude=sampleNamesToExclude,
                           numberOfMS1features=numberOfMS1features,
                           precursorLabels=precursorLabels,
                           grouXXXps=c("TRI", "LVS"),
                           metaboliteProfileColumnNames=metaboliteProfileColumnNames,
                           dataColumnIndecesFunctionFromGroupIndex=dataColumnIndecesFunctionFromGroupIndex,
                           dataColumnsNameFunctionFromGroupIndex=dataColumnsNameFunctionFromGroupIndex,
                           dataColumnsNameFunctionFromGroupName=dataColumnsNameFunctionFromGroupName,
                           dataColumnsNameFunctionFromGroupNames=dataColumnsNameFunctionFromGroupNames,
                           groupNameFunctionFromDataColumnName=groupNameFunctionFromDataColumnName,
                           tagsSector=tagsSector,
                           metaboliteProfile=metaboliteProfile,
                           progress=FALSE)

expect_equal(min(result$dataFrameMeasurements[,1]), 0)

})

test_that("HCA works", {
  load(system.file("extdata/testdata/calculateDistanceMatrix.Rdata", package = "MetFamily"))
  
  result <- calculateDistanceMatrix(dataList=dataList, 
                                    filter=filter, 
                                    distanceMeasure = "Jaccard", progress = FALSE)
  expect_equal(sum(result$distanceMatrix), 46522.7, tolerance = 0.005)
  expect_equal(sum(result$filter), 294508)
  expect_true(result$distanceMeasure == "Jaccard")
  
})
