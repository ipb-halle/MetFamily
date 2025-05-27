test_that("metaboliteProfileParsing works", {
  load(system.file("extdata/testdata/processMS1data.Rdata", package = "MetFamily"))

  d <- dataColumnsNameFunctionFromGroupIndex(groupIdx = 2, sampleNamesToExclude = NA)
  
  ## Test dimensions  
  expect_equal(length(d), 3)
  
  result <- processMS1data(sampleNamesToExclude=sampleNamesToExclude,
                           numberOfMS1features=numberOfMS1features,
                           precursorLabels=precursorLabels,
                           sampleClasses=c("TRI", "LVS"),
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
