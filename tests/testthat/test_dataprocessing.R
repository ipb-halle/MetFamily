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

test_that("calculateDistanceMatrix works", {
  load(system.file("extdata/testdata/calculateDistanceMatrix.Rdata", package = "MetFamily"))
  
  cluster <- parallel::makeForkCluster()
  doParallel::registerDoParallel(cl=cluster)
  result <- calculateDistanceMatrix(dataList=dataList, 
                                    filter=filter, 
                                    distanceMeasure = "Jaccard", progress = FALSE)
  parallel::stopCluster(cluster)
  expect_equal(sum(result$distanceMatrix), 46522.7, tolerance = 0.005)
  expect_equal(sum(result$filter), 294508)
  expect_true(result$distanceMeasure == "Jaccard")
  expect_true(all(diag(result$distanceMatrix)==0))
  
  i<-c(1, 215, 39, 107, 48, 49, 219)
  j<-c(219, 139, 130, 147, 13, 90, 1)
  a <- c(1, 0.96875, 0.875, 1, 1, 0.888888888888889, 0, 1, 0.918032786885246,  
         1, 1, 1, 1, 0.978260869565217, 1, 0.925925925925926, 1, 1, 1,  1, 
         0.973684210526316, 1, 0.968253968253968, 1, 1, 1, 1, 1, 1,  1, 1, 1, 
         1, 1, 1, 1, 0.976744186046512, 0.947368421052632, 1,  0.95, 0.95, 
         0.91304347826087, 0, 1, 1, 1, 1, 1, 1)
  
  expect_equal(as.numeric(result$distanceMatrix[i,j]), a)
  
})
