# Load example data
example_file_path <- system.file("extdata/showcase/Project_file_showcase_annotated.csv.gz", 
                                 package = "MetFamily")
dataList <- readClusterDataFromProjectFile(file = example_file_path, progress = FALSE)

# Load existing filter object
filterObjPath <- system.file("extdata/testdata/canopus/filterPca_canopus.rds", package = "MetFamily")
filterObj <- readRDS(filterObjPath)

# Apply filtering
filteredResult <- filterData(
  dataList = dataList,
  sampleClasses = filterObj$groupSet,
  sampleSet = filterObj$sampleSet,
  filterBySamples = filterObj$filterBySamples,
  filter_average = filterObj$filter_average,
  filter_lfc = filterObj$filter_lfc,
  filterList_ms2_masses = filterObj$filter_ms2_masses,
  filter_ms2_ppm = filterObj$filter_ms2_ppm,
  filter_ms1_masses = filterObj$filter_ms1_masses,
  filter_ms1_ppm = filterObj$filter_ms1_ppm,
  includeIgnoredPrecursors = filterObj$includeIgnoredPrecursors,
  progress = FALSE
)


test_that("HCA analysis produces expected values", {
  # Use original subset size for testing
  test_filter <- filteredResult$filter[1:100]
  
  # Calculate distance matrix
  distanceResult <- calculateDistanceMatrix(
    dataList = dataList, 
    filter = test_filter, 
    distanceMeasure = "Jaccard", 
    progress = FALSE
  )
  
  # Calculate cluster
  filterHca <- list(
    filter = test_filter,
    numberOfPrecursorsFiltered = length(test_filter),
    sampleClasses = filteredResult$sampleClasses
  )
  
  clusterResult <- calculateCluster(
    dataList = dataList,
    filterObj = filterHca,
    distanceMatrix = distanceResult$distanceMatrix,
    method = "ward.D",
    distanceMeasure = "Jaccard",
    progress = FALSE
  )
  
  # Test basic structure
  expect_equal(nrow(distanceResult$distanceMatrix), length(test_filter))
  expect_equal(clusterResult$numberOfPrecursorsFiltered, length(test_filter))
  expect_true(isSymmetric(distanceResult$distanceMatrix))
  
  # Test specific distance pairs 
  expect_equal(round(distanceResult$distanceMatrix[25, 72], 6), 0.978723)
  expect_equal(round(distanceResult$distanceMatrix[42, 43], 6), 0.866667)  
  expect_equal(round(distanceResult$distanceMatrix[9, 16], 6), 0.954545)
  expect_equal(round(distanceResult$distanceMatrix[11, 99], 6), 0.982456)
  
  # Test some distance values 
  expect_equal(round(min(distanceResult$distanceMatrix[upper.tri(distanceResult$distanceMatrix)]), 6), 0.454545)
  expect_equal(round(median(distanceResult$distanceMatrix[upper.tri(distanceResult$distanceMatrix)]), 6), 1.0)
   
  # Test cluster results 
  expect_equal(clusterResult$numberOfPois, 199)
  expect_equal(round(clusterResult$cluster$height[1], 6), 0.454545)
  expect_equal(round(clusterResult$cluster$height[10], 6), 0.619048)
  expect_equal(clusterResult$cluster$order[1], 83)
  expect_equal(clusterResult$cluster$order[5], 95)
  expect_equal(round(max(clusterResult$cluster$height), 6), 2.773333)
})
test_that("Different distance measures work", {
  # Test with small subset
  test_filter_small <- filteredResult$filter[1:6]
  
  # Test Jaccard vs NDP (which should give different results)
  jaccard_result <- calculateDistanceMatrix(
    dataList = dataList, 
    filter = test_filter_small, 
    distanceMeasure = "Jaccard", 
    progress = FALSE
  )
  
  ndp_result <- calculateDistanceMatrix(
    dataList = dataList, 
    filter = test_filter_small, 
    distanceMeasure = "NDP (Normalized dot product)", 
    progress = FALSE
  )
  
  # Test basic structure
  expect_equal(nrow(jaccard_result$distanceMatrix), length(test_filter_small))
  expect_true(isSymmetric(jaccard_result$distanceMatrix))
  expect_true(isSymmetric(ndp_result$distanceMatrix))
  expect_true(all(abs(diag(jaccard_result$distanceMatrix)) < 1e-10))
  
  # Test that different measures give different results
  expect_false(identical(jaccard_result$distanceMatrix, ndp_result$distanceMatrix))
})

test_that("Tree analysis works", {
  # Test with minimal dataset
  test_filter_tiny <- filteredResult$filter[1:6]
  
  distanceResult <- calculateDistanceMatrix(
    dataList = dataList, 
    filter = test_filter_tiny, 
    distanceMeasure = "Jaccard", 
    progress = FALSE
  )
  
  # Create cluster for tree analysis
  dist <- stats::as.dist(m = distanceResult$distanceMatrix)
  cluster <- hclust(d = dist, method = "ward.D")
  cluster$labels <- dataList$precursorLabels[test_filter_tiny]
  
  # Test tree analysis
  resultObj <- analyzeTreeFromRoot(dataList, cluster = cluster, test_filter_tiny)
  
  # Basic structure tests
  expect_equal(length(resultObj$innerNodeMembersTreeLeaves), length(test_filter_tiny) - 1)
  expect_equal(length(resultObj$leafHeights), length(test_filter_tiny))
  expect_true(is.list(resultObj$innerNodeFeaturesIntersection))
  expect_true(is.list(resultObj$innerNodeFeaturesUnion))
}) 