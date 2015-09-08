
#install.packages("squash")
library("squash")
#install.packages("FactoMineR")
#library("FactoMineR")
#require(devtools)
#install_github('rCharts', 'ramnathv')
#library("rCharts")
#source("http://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")
#library("pcaMethods")
#install.packages("cba")
#library("cba")
#install.packages("matrixStats")
library("matrixStats")

#########################################################################################
## helper
analyzeTreeFromRoot <- function(dataList, cluster, filter){
  numberOfPrecursorsFiltered <- length(filter)
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  rootIndex <- length(cluster$height)
  
  ## create fields and compute stuff
  innerNodeMembersTreeHere <<- list()
  innerNodeMembersHere <<- list()
  innerNodeFeaturesIntersectionHere <<- list()
  innerNodeFeaturesAnnotations <<- list()
  innerNodeFeaturesUnionHere <<- list()
  innerNodeFeaturesCountsHere <<- list()
  innerNodeFeaturesMeanAbundanceHere <<- list()
  innerNodeFeaturesCounterIntersectionHere <<- list()
  innerNodeFeaturesCounterUnionHere <<- list()
  innerNodePositionHere <<- vector(mode = "numeric", length = numberOfInnerNodes)
  leafHeightsHere <<- vector(mode = "numeric", length = numberOfPrecursorsFiltered)
  for(i in 1:numberOfInnerNodes)
    innerNodeMembersHere[[i]] <<- NA
  
  analyzeTree(dataList, cluster, filter, rootIndex)
  
  ## box
  resultObj <- list()
  resultObj$innerNodeMembersTree <- innerNodeMembersTreeHere
  resultObj$innerNodeMembers <- innerNodeMembersHere
  resultObj$innerNodeFeaturesIntersection <- innerNodeFeaturesIntersectionHere
  resultObj$innerNodeFeaturesAnnotations <- innerNodeFeaturesAnnotations
  resultObj$innerNodeFeaturesUnion <- innerNodeFeaturesUnionHere
  resultObj$innerNodeFeaturesCounts <- innerNodeFeaturesCountsHere
  resultObj$innerNodeFeaturesMeanAbundance <- innerNodeFeaturesMeanAbundanceHere
  resultObj$innerNodeFeaturesIntersectionCounter <- innerNodeFeaturesCounterIntersectionHere
  resultObj$innerNodeFeaturesUnionCounter <- innerNodeFeaturesCounterUnionHere
  resultObj$innerNodePosition <- innerNodePositionHere
  resultObj$leafHeights <- leafHeightsHere
  
  return(resultObj)
}
analyzeTree <- function(dataList, cluster, filter, nodeIdx){
  if(nodeIdx < 0){
    ###################################
    ## leaf
    leafIdx <- -nodeIdx
    precursorIndex <- filter[leafIdx]
    featureValues <- dataList$featureMatrix[precursorIndex, ]
    featuresBinary <- featureValues != 0
    featuresBinaryIntersection <- featuresBinary
    featuresBinaryUnion        <- featuresBinary
    featureCounts <- vector(length = length(featuresBinary))
    featureCounts[featuresBinary] <- 1
    featuresCounter <- sum(featuresBinary)
    featuresCounterIntersection <- featuresCounter
    featuresCounterUnion        <- featuresCounter
    position <- match(x = leafIdx, table = cluster$order)
    featuresAnnotations <- dataList$annoArrayOfLists[[precursorIndex]]
    
    ## box
    resultObj <- list()
    resultObj$membersTree <- leafIdx
    resultObj$members <- precursorIndex
    resultObj$featuresBinaryIntersection <- featuresBinaryIntersection
    resultObj$featuresBinaryUnion <- featuresBinaryUnion
    #resultObj$features <- which(featuresBinary)
    #resultObj$featuresCounts <- rep(x = 1, times = resultObj$featuresCounter)
    resultObj$featuresCounts <- featureCounts
    #resultObj$featuresMeanAbundance <- featureValues[resultObj$features]
    resultObj$featuresMeanAbundance <- featureValues
    resultObj$featuresCounterIntersection <- featuresCounterIntersection
    resultObj$featuresCounterUnion        <- featuresCounterUnion
    #resultObj$featuresCounter <- length(resultObj$features)
    resultObj$position <- position
    resultObj$featuresAnnotations <- featuresAnnotations
    
    return(resultObj)
  } else {
    ###################################
    ## inner node
    resultObj.l <- analyzeTree(dataList, cluster, filter, cluster$merge[nodeIdx, 1])
    resultObj.r <- analyzeTree(dataList, cluster, filter, cluster$merge[nodeIdx, 2])
    
    membersTree <- c(resultObj.l$membersTree, resultObj.r$membersTree)
    members <- c(resultObj.l$members, resultObj.r$members)
    featuresBinaryIntersection <- resultObj.l$featuresBinaryIntersection & resultObj.r$featuresBinaryIntersection
    featuresBinaryUnion        <- resultObj.l$featuresBinaryUnion        | resultObj.r$featuresBinaryUnion
    #features <- which(featuresBinary)
    featuresCounts <- resultObj.l$featuresCounts + resultObj.r$featuresCounts
    featuresMeanAbundance <- (resultObj.l$featuresMeanAbundance * resultObj.l$featuresCounts + resultObj.r$featuresMeanAbundance * resultObj.r$featuresCounts) / featuresCounts
    featuresCounterIntersection <- sum(featuresBinaryIntersection)
    featuresCounterUnion        <- sum(featuresBinaryUnion)
    position <- (resultObj.l$position + resultObj.r$position) / 2
    featuresAnnotations <- intersect(x = unlist(resultObj.l$featuresAnnotations), y = unlist(resultObj.r$featuresAnnotations))
    
    ## box
    resultObj <- list()
    resultObj$membersTree <- membersTree
    resultObj$members <- members
    resultObj$featuresBinaryIntersection <- featuresBinaryIntersection
    resultObj$featuresBinaryUnion <- featuresBinaryUnion
    #resultObj$features <- union(x = resultObj.l$features, y = resultObj.r$features)
    #resultObj$features <- features
    resultObj$featuresCounts <- featuresCounts
    resultObj$featuresMeanAbundance <- featuresMeanAbundance
    resultObj$featuresCounterIntersection <- featuresCounterIntersection
    resultObj$featuresCounterUnion <- featuresCounterUnion
    resultObj$position <- position
    resultObj$featuresAnnotations <- featuresAnnotations
    
    ## set values
    innerNodeMembersTreeHere[[nodeIdx]] <<- resultObj$membersTree
    innerNodeMembersHere[[nodeIdx]] <<- resultObj$members
    innerNodePositionHere[[nodeIdx]] <<- resultObj$position
    innerNodeFeaturesIntersectionHere[[nodeIdx]] <<- which(resultObj$featuresBinaryIntersection)
    innerNodeFeaturesAnnotations[[nodeIdx]] <<- resultObj$featuresAnnotations
    innerNodeFeaturesUnionHere[[nodeIdx]] <<- which(resultObj$featuresBinaryUnion)
    #innerNodeFeaturesBinaryHere[[nodeIdx]] <<- resultObj$featuresBinary
    innerNodeFeaturesCountsHere[[nodeIdx]] <<- resultObj$featuresCounts
    innerNodeFeaturesMeanAbundanceHere[[nodeIdx]] <<- resultObj$featuresMeanAbundance
    innerNodeFeaturesCounterIntersectionHere[[nodeIdx]] <<- resultObj$featuresCounterIntersection
    innerNodeFeaturesCounterUnionHere[[nodeIdx]] <<- resultObj$featuresCounterUnion
    
    ## leaf heights in case of leafs
    if(cluster$merge[nodeIdx, 1] < 0)
      leafHeightsHere[[-cluster$merge[nodeIdx, 1]]] <<- cluster$height[[nodeIdx]]
    if(cluster$merge[nodeIdx, 2] < 0)
      leafHeightsHere[[-cluster$merge[nodeIdx, 2]]] <<- cluster$height[[nodeIdx]]
    
    return(resultObj)
  }
}
getSetOfSubTreesFromRoot <- function(filter, clusterDataList, yesNoFunction){
  rootIndex <- length(clusterDataList$cluster$height)
  result <- getSetOfSubTrees(filter, clusterDataList, yesNoFunction, rootIndex)
  if(result$criterionFulfilled)
    result$results <- c(result$results, rootIndex)
  return(result$results)
}
getSetOfSubTrees <- function(filter, clusterDataList, yesNoFunction, index){
  if(index<0){ # it is a leaf
    leafIndex <- -index
    precursorIndex <- filter$filter[[leafIndex]]
    criterionFulfilled <- yesNoFunction(precursorIndex)
    
    result <- list()
    result$results <- c()
    result$criterionFulfilled <- criterionFulfilled
    
    return(result)
  }
  
  result.l  <- getSetOfSubTrees(filter, clusterDataList, yesNoFunction, clusterDataList$cluster$merge[index, 1])
  result.r  <- getSetOfSubTrees(filter, clusterDataList, yesNoFunction, clusterDataList$cluster$merge[index, 2])
  
  result <- list()
  result$results <- c(result.l$results, result.r$results)
  if(result.l$criterionFulfilled & result.r$criterionFulfilled){
    ## both children do comprise the fragment
    result$criterionFulfilled <- TRUE
  } else {
    result$criterionFulfilled <- FALSE
    if((!result.l$criterionFulfilled & result.r$criterionFulfilled) | (result.l$criterionFulfilled & !result.r$criterionFulfilled)){
      ## only one child does comprise the fragment
      if(result.l$criterionFulfilled)
        result$results <- c(result$results, clusterDataList$cluster$merge[index, 1])
      else
        result$results <- c(result$results, clusterDataList$cluster$merge[index, 2])
    } else {
      ## no child does comprise the fragment
    }
  }
  
  return(result)
}
colorSubTree <- function(cluster, index, lwd = 1, lty = 1, col = "black"){
  if(index<0){ # it is a leaf
    a2r_counter <<- a2r_counter + 1
    return(list(
      x = a2r_counter
    ))       
  }
  
  h.m   <- cluster$height[index]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- cluster$merge[index,1]
  
  h.l <- if(index.l<0) 0 else cluster$height[index.l]
  
  out.l   <- colorSubTree(cluster = cluster, index = index.l, col=col, lty=lty, lwd=lwd)
  x.l     <- out.l$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- cluster$merge[index,2]
  h.r <- if(index.r<0) 0 else cluster$height[index.r]
  out.r   <- colorSubTree(cluster = cluster, index = index.r, col=col, lty=lty, lwd=lwd)
  x.r     <- out.r$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  x.m  <- (x.r + x.l) / 2  
  
  segments(
    x0  = c(x.l, x.l, x.r),
    x1  = c(x.l, x.r, x.r),
    y0  = c(h.l, h.m, h.r),
    y1  = c(h.m, h.m, h.m),
    col = col,
    lty = lty,
    lwd = lwd
  )
  
  list(x=x.m)
}
colorSubTreeForAnnotations <- function(cluster, index, innerNodeAnnotations, innerNodeColors, parentIndex, lwd = 1, lty = 1){
  if(index<0){ # it is a leaf
    a2r_counter <<- a2r_counter + 1
    return(list(
      x = a2r_counter
    ))       
  }
  
  if(is.null(parentIndex))
    parentAnnotations <- NULL
  else{
    if(length(innerNodeAnnotations) < parentIndex)
      parentAnnotations <- NULL
    else
      parentAnnotations <- innerNodeAnnotations[[parentIndex]]
  }
  
  if(length(innerNodeAnnotations) < index)
    currentAnnotations <- NULL
  else
    currentAnnotations <- innerNodeAnnotations[[index]]
  
  newAnnotations <- setdiff(x = currentAnnotations, y = parentAnnotations)
  if(length(newAnnotations) == 0){
    if(length(parentAnnotations) == 0)
      color <- "black"
    else
      color <- innerNodeColors[[parentIndex]][[1]]
  } else
    color <- innerNodeColors[[index]][[1]]
  col <- color
  
  print(paste("i", index, "pi", parentIndex, "a", currentAnnotations, "na", newAnnotations, "pa", parentAnnotations, "c", color))
  
  h.m   <- cluster$height[index]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- cluster$merge[index,1]
  
  h.l <- if(index.l<0) 0 else cluster$height[index.l]
  
  out.l   <- colorSubTreeForAnnotations(cluster = cluster, index = index.l, innerNodeAnnotations = innerNodeAnnotations, innerNodeColors = innerNodeColors, parentIndex = index, lty=lty, lwd=lwd)
  x.l     <- out.l$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- cluster$merge[index,2]
  h.r <- if(index.r<0) 0 else cluster$height[index.r]
  out.r   <- colorSubTreeForAnnotations(cluster = cluster, index = index.r, innerNodeAnnotations = innerNodeAnnotations, innerNodeColors = innerNodeColors, parentIndex = index, lty=lty, lwd=lwd)
  x.r     <- out.r$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  x.m  <- (x.r + x.l) / 2  
  
  segments(
    x0  = c(x.l, x.l, x.r),
    x1  = c(x.l, x.r, x.r),
    y0  = c(h.l, h.m, h.r),
    y1  = c(h.m, h.m, h.m),
    col = col,
    lty = lty,
    lwd = lwd
  )
  
  list(x=x.m)
}
#########################################################################################
## annotate and process matrix
readClusterData <- function(file, progress = FALSE){
  if(progress)  setProgress(value = 0, detail = "Parsing")
  dataFrame <- read.csv(file = file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  if(progress)  incProgress(amount = 0.1, detail = "Preprocessing")
  columnNames <- unlist(as.matrix(dataFrame[3, ]))
  colnames(dataFrame) <- columnNames
  #endOfAnnotation <- tail(x = which(is.na(suppressWarnings(as.numeric(columnNames)))), n = 1)[[1]]
  endOfAnnotation <- max(which(dataFrame[1, ] == ""))
  tagsSector <- dataFrame[2, 1:endOfAnnotation]
  measurementColumns <- c(14, endOfAnnotation)
  numberOfMeasurementColumns <- measurementColumns[[2]] - measurementColumns[[1]] + 1
  fragmentColumns <- c(endOfAnnotation + 1, ncol(dataFrame))
  fragmentMasses <- as.numeric(columnNames[fragmentColumns[[1]] : fragmentColumns[[2]]])
  fragmentFrequency <- as.numeric(unlist(as.matrix(dataFrame[1, fragmentColumns[[1]] : fragmentColumns[[2]]])))
  fragmentAbundance <- as.numeric(unlist(as.matrix(dataFrame[2, fragmentColumns[[1]] : fragmentColumns[[2]]])))
  numberOfFragments <- length(fragmentMasses)
  dataFrameOriginal <- dataFrame
  dataFrameHeader <- dataFrame[1:3, ]
  dataFrame <- dataFrame[4:nrow(dataFrame), ]
  idColumns <- which(tagsSector == "ID")
  #precursorLabels <- paste(dataFrame[, 1], dataFrame[, 2], sep = " / ")
  precursorLabels <- apply(X = dataFrame[, idColumns], MARGIN = 1, FUN = paste, collapse = " / ")
  
  dupplicated <- which(duplicated(precursorLabels))
  if(length(dupplicated) > 0){
    precursorLabels <- precursorLabels[-dupplicated]
    dataFrame         <- dataFrame[-dupplicated, ]
    dataFrameOriginal <- dataFrameOriginal[-(dupplicated + 3)]
  }
  
  rownames(dataFrame) <- precursorLabels
  rownames(dataFrameOriginal) <- c("HeaderForFragmentCounts", "HeaderForGroupsAndFragmentIntensities", "Header", precursorLabels)
  numberOfPrecursors <- nrow(dataFrame)
  
  #dataColumns <- which(!is.na(suppressWarnings(as.numeric(tagsSector))))
  dataColumns <- which(tagsSector != "ID" & tagsSector != "")
  groups <- unique(unlist(tagsSector[dataColumns]))
  groupsStartEnd <- list()
  for(groupIdx in 1:length(groups))
    groupsStartEnd[[groupIdx]] <- c(min(which(tagsSector == groups[[groupIdx]])), max(which(tagsSector == groups[[groupIdx]])))
  groupsStartEndMatrix <- t(matrix(data = unlist(groupsStartEnd), nrow = 2))
  rownames(groupsStartEndMatrix) <- groups
  colnames(groupsStartEndMatrix) <- c("Start", "End")
  
  ####################
  ## measurement data to colors
  if(progress)  incProgress(amount = 0.1, detail = "Coloring")
  if(progress)  incProgress(amount = 0, detail = "Coloring init")
  
  dataFrameMeasurements <- data.frame(matrix(nrow = nrow(dataFrame), ncol = 0))
  rownames(dataFrameMeasurements) <- rownames(dataFrame)
  
  ## mean of groups
  if(progress)  incProgress(amount = 0, detail = "Coloring gather data")
  dataColumnsNameFunctionFromIndex <- function(groupIdx){
    paste(groups[[groupIdx]], "_", columnNames[groupsStartEndMatrix[groupIdx, 1]:groupsStartEndMatrix[groupIdx, 2]], sep = "")
  }
  dataColumnsNameFunctionFromName <- function(group){
    dataColumnsNameFunctionFromIndex(match(x = group, table = groups))
  }
  dataColumnsNameFunctionFromNames <- function(groups){
    unlist(lapply(X = groups, FUN = dataColumnsNameFunctionFromName))
  }
  groupNameFunctionFromDataColumnName <- function(dataColumnName){
    groupIdx <- which(unlist(lapply(X = groups, FUN = function(x){
      dataColumnNames <- dataColumnsNameFunctionFromName(x)
      any(dataColumnNames == dataColumnName)
    })))
    groups[[groupIdx]]
  }
  for(groupIdx in 1:length(groups)){
    dataColumnNames <- dataColumnsNameFunctionFromIndex(groupIdx)
    dataFrameMeasurements[, dataColumnNames] <- data.matrix(dataFrame[, groupsStartEndMatrix[[groupIdx, 1]]:groupsStartEndMatrix[[groupIdx, 2]]])
  }
  dataMeanColumnNameFunctionFromIndex  <- function(groupIdx){
    return(dataMeanColumnNameFunctionFromName(groups[[groupIdx]]))
  }
  dataMeanColumnNameFunctionFromName  <- function(group){
    return(paste(group, "_mean", sep = ""))
  }
  dataMeanColumnNames <- list()
  for(groupIdx in 1:length(groups)){
    dataColumnName <- dataMeanColumnNameFunctionFromIndex(groupIdx)
    dataMeanColumnNames[[groupIdx]] <- dataColumnName
    dataFrameMeasurements[, dataColumnName] <- apply(X = data.matrix(dataFrame[, groupsStartEndMatrix[[groupIdx, 1]]:groupsStartEndMatrix[[groupIdx, 2]]]), MARGIN = 1, FUN = mean)
  }
  dataMeanColumnNames <- unlist(dataMeanColumnNames)
  
  if(progress)  incProgress(amount = 0, detail = "Coloring LFC")
  ## log fold change between groups
  lfcColumnNameFunctionFromIndex <- function(groupIdxOne, groupIdxTwo){
    lfcColumnNameFunctionFromName(groups[[groupIdxOne]], groups[[groupIdxTwo]])
  }
  lfcColumnNameFunctionFromName <- function(groupOne, groupTwo){
    return(paste("LFC", groupOne, "vs", groupTwo, sep = "_"))
  }
  lfcColumnNames <- list()
  for(groupIdx1 in 1:length(groups))
    for(groupIdx2 in 1:length(groups)){
      lfcColumnName <- lfcColumnNameFunctionFromIndex(groupIdx1, groupIdx2)
      lfcColumnNames[[length(lfcColumnNames) + 1]] <- lfcColumnName
      dataFrameMeasurements[, lfcColumnName] <- log(
        x = dataFrameMeasurements[, dataMeanColumnNameFunctionFromIndex(groupIdx2)] / dataFrameMeasurements[, dataMeanColumnNameFunctionFromIndex(groupIdx1)], 
        base = 2
      )
      
      ## tackle zero values
      dataFrameMeasurements[is.na(dataFrameMeasurements[, lfcColumnName]), lfcColumnName] <- 0
      dataFrameMeasurements[is.infinite(dataFrameMeasurements[, lfcColumnName]), lfcColumnName] <- 0
    }
  lfcColumnNames <- unlist(lfcColumnNames)
  
  groupNameFromGroupIndex <- function(groupIdx){
    return(groups[[groupIdx]])
  }
  groupIdxFromGroupName <- function(group){
    return(match(x = group, table = groups))
  }
  
  ## map to colors
  if(progress)  incProgress(amount = 0, detail = "Coloring matrix")
  matrixDataFrame <- data.matrix(dataFrameMeasurements)
  
  absMax <- max(matrixDataFrame[, dataMeanColumnNames])
  logFoldChangeMinMax <- c(min(dataFrameMeasurements[, lfcColumnNames]), max(dataFrameMeasurements[, lfcColumnNames]))
  logFoldChangeMax <- max(abs(logFoldChangeMinMax))
  colorMapAbsoluteData  <- makecmap(
    x = c(0, absMax), n = 100, 
    colFn = colorRampPalette(c('white', 'black'))
  )
  colorMapLogFoldChange <- makecmap(
    x = c(-logFoldChangeMax, logFoldChangeMax), n = 100, 
    colFn = colorRampPalette(c('blue', 'white', 'red'))
  )
  
  columnGroupLabels <- sapply(X = groups, FUN = function(x){ rep(x = x, times = length(dataColumnsNameFunctionFromName(x))) })
  columnGroupOrgLabels <- columnNames[min(groupsStartEndMatrix):max(groupsStartEndMatrix)]
  
  ## translate and box colors
  if(progress)  incProgress(amount = 0, detail = "Coloring box")
  colorDataFrame <- dataFrameMeasurements
  colorDataFrame[, dataMeanColumnNames] <- cmap(matrixDataFrame[, dataMeanColumnNames], colorMapAbsoluteData)
  colorDataFrame[, lfcColumnNames]      <- cmap(matrixDataFrame[, lfcColumnNames     ], colorMapLogFoldChange)
  colorMatrixDataFrame <- as.matrix(colorDataFrame)
  
  ##########################################
  ## compute features
  
  if(progress)  incProgress(amount = 0.1, detail = "Features")
  ## get features
  
  dataFrameData <- as.matrix(dataFrame[, fragmentColumns[[1]] : fragmentColumns[[2]]])
  
  ## real features
  featureMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfFragments)
  rownames(featureMatrix) <- precursorLabels
  colnames(featureMatrix) <- fragmentMasses
  for(i in 1:numberOfPrecursors){
    if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Features ", i, " / ", numberOfPrecursors, sep = ""))
    featureVector <- as.numeric(unlist(dataFrameData[i, ]))
    featureMatrix[i, ] <- featureVector
    featureMatrix[i, is.na(featureMatrix[i, ])] <- 0
  }
  
  ## binary features
  if(progress)  incProgress(amount = 0, detail = "Features binary")
  featureIndeces <- list()
  #featureMatrixBinary <- matrix(nrow = numberOfPrecursors, ncol = numberOfFragments)
  for(i in 1:numberOfPrecursors){
    if(progress)  incProgress(amount = 0.1 / numberOfPrecursors, detail = paste("Feature indeces ", i, " / ", numberOfPrecursors, sep = ""))
    featureVectorBinary <- featureMatrix[i, ] != 0
    #featureMatrixBinary[i, ] <- featureVectorBinary
    featureIndeces[[i]] <- which(featureVectorBinary)
  }
  featureIndexMatrix <- matrix(nrow = numberOfPrecursors, ncol = max(sapply(X = featureIndeces, FUN = length)))
  for(i in 1:numberOfPrecursors)
    featureIndexMatrix[i, 1:length(featureIndeces[[i]])] <- featureIndeces[[i]]
  
  ## remove unused columns
  fragmentThere <- apply(X = featureMatrix, MARGIN = 2, FUN = function(x){any(x != 0)})
  minimumMass <- min(fragmentMasses[fragmentThere])
  maximumMass <- max(fragmentMasses[fragmentThere])
  
  ## annotations
  annoArrayOfLists    <- vector(mode='list', length=numberOfPrecursors)
  annoArrayIsArtifact <- vector(mode='logical', length=numberOfPrecursors)
  annoPresentAnnotationsList <- list()
  annoPresentColorsList <- list()
  
  annoPresentAnnotationsList[[1]] <- "Ignore"
  annoPresentColorsList[[1]] <- "red"
  
  ##########################################
  ## box
  if(progress)  incProgress(amount = 0.1, detail = "Boxing")
  dataList <- list()
  ## data
  dataList$dataFrameOriginal <- dataFrameOriginal
  dataList$dataFrameHeader <- dataFrameHeader
  dataList$dataFrame <- dataFrame
  dataList$numberOfPrecursors <- numberOfPrecursors
  dataList$groups <- groups
  dataList$columnGroupLabels <- columnGroupLabels
  dataList$columnGroupOrgLabels <- columnGroupOrgLabels
  #dataList$groupsStartEndMatrix <- groupsStartEndMatrix
  ## data: fragments
  dataList$fragmentMasses <- fragmentMasses
  dataList$fragmentFrequency <- fragmentFrequency
  dataList$fragmentAbundance <- fragmentAbundance
  dataList$minimumMass <- minimumMass
  dataList$maximumMass <- maximumMass
  dataList$precursorLabels <- precursorLabels
  ## data: abundancies
  dataList$dataFrameMeasurements <- dataFrameMeasurements
  dataList$logFoldChangeMax <- logFoldChangeMax
  dataList$absMax <- absMax
  dataList$colorMatrixDataFrame <- colorMatrixDataFrame
  dataList$colorMapAbsoluteData <- colorMapAbsoluteData
  dataList$colorMapLogFoldChange <- colorMapLogFoldChange
  dataList$dataColumnsNameFunctionFromName <- dataColumnsNameFunctionFromName
  dataList$dataColumnsNameFunctionFromIndex <- dataColumnsNameFunctionFromIndex
  dataList$dataColumnsNameFunctionFromNames <- dataColumnsNameFunctionFromNames
  dataList$groupNameFunctionFromDataColumnName <- groupNameFunctionFromDataColumnName
  dataList$dataMeanColumnNameFunctionFromName <- dataMeanColumnNameFunctionFromName
  dataList$dataMeanColumnNameFunctionFromIndex <- dataMeanColumnNameFunctionFromIndex
  dataList$lfcColumnNameFunctionFromName <- lfcColumnNameFunctionFromName
  dataList$lfcColumnNameFunctionFromIndex <- lfcColumnNameFunctionFromIndex
  dataList$groupNameFromGroupIndex <- groupNameFromGroupIndex
  dataList$groupIdxFromGroupName <- groupIdxFromGroupName
  
  ## features
  #dataList$distanceMatrix <- distanceMatrix
  dataList$featureMatrix <- featureMatrix
  #dataList$featureMatrixBinary <- featureMatrixBinary
  dataList$featureIndeces <- featureIndeces
  dataList$featureIndexMatrix <- featureIndexMatrix
  ## annotations
  dataList$annoArrayOfLists <- annoArrayOfLists
  dataList$annoArrayIsArtifact <- annoArrayIsArtifact
  dataList$annoPresentAnnotationsList <- annoPresentAnnotationsList
  dataList$annoPresentColorsList <- annoPresentColorsList
  
  if(progress)  setProgress(1)
  
  return(dataList)
}
filterData <- function(dataList, groups, filter_average, filter_lfc, filterList_ms2_masses, filter_ms2_ppm, includeIgnoredPrecursors, progress = FALSE){
  ##########################################
  ## filter
  filter <- rep(x = TRUE, times = dataList$numberOfPrecursors)
  
  ## filter_average
  if(!is.null(filter_average))
    filter <- filter & apply(dataList$dataFrameMeasurements[, sapply(X = groups, FUN = dataList$dataMeanColumnNameFunctionFromName)], 1, mean) >= filter_average
  
  ## filter_lfc
  if(!is.null(filter_lfc)){
    if(filter_lfc != 0){
      if(length(groups) != 2){  stop("The number of groups for LFC is not equal to two!") }
      if(filter_lfc > 0)
        filter <- filter & dataList$dataFrameMeasurements[, dataList$lfcColumnNameFunctionFromName(groups[[1]], groups[[2]])] >= filter_lfc
      else
        filter <- filter & dataList$dataFrameMeasurements[, dataList$lfcColumnNameFunctionFromName(groups[[1]], groups[[2]])] <= filter_lfc
    }
  }
  
  ## filter_ms2_masses, filter_ms2_ppm
  if(!is.null(filterList_ms2_masses) & !is.null(filter_ms2_ppm) & length(filterList_ms2_masses) > 0){
    error <- abs(dataList$fragmentMasses) * filter_ms2_ppm / 1E6
    filterAll <- rep(x = FALSE, times = dataList$numberOfPrecursors)
    for(filter_ms2_masses in filterList_ms2_masses){
      fragmentColumns <- vector(mode = "numeric")
      
      for(fragmentIndex in 1:length(filter_ms2_masses)){
        distances <- abs(dataList$fragmentMasses - filter_ms2_masses[[fragmentIndex]])
        if(any(distances <= error, na.rm = TRUE)){
          minColumn <- which.min(distances)
          fragmentColumns[[length(fragmentColumns) + 1]] <- minColumn
        }
      }
      filterTemp <- rep(x = TRUE, times = dataList$numberOfPrecursors)
      if(length(fragmentColumns) > 0){
        for(columns in fragmentColumns)
          filterTemp <- filterTemp & dataList$featureMatrix[, columns] != 0
      } else {
        filterTemp <- rep(x = FALSE, times = dataList$numberOfPrecursors)
      }
      filterAll <- filterAll | filterTemp
    }
    filter <- filter & filterAll
  }
  
  ## include ignored precursors
  if(!includeIgnoredPrecursors)
    filter <- filter & !dataList$annoArrayIsArtifact
  
  filter <- which(filter)
  
  resultObj <<- list()
  resultObj$filter <<- filter
  resultObj$numberOfPrecursorsFiltered <<- length(filter)
  resultObj$groups <<- groups
  resultObj$filter_average <<- filter_average
  resultObj$filter_lfc <<- filter_lfc
  resultObj$filterList_ms2_masses <<- filterList_ms2_masses
  resultObj$filter_ms2_ppm <<- filter_ms2_ppm
  resultObj$includeIgnoredPrecursors <<- includeIgnoredPrecursors
  
  return (resultObj)
}
calculateDistanceMatrix <- function(dataList, filter, distance = "Jaccard", progress = FALSE){
  
  filter <- filter$filter
  numberOfPrecursors <- length(filter)
  
  ## compute distance matrix:
  distanceMatrix <- NULL
  switch(distance,
    "Jaccard"={
      featureIndeces <- dataList$featureIndeces[filter]
      
      distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            distanceMatrix[i, j] <- 0
            next
          }
          
          #distanceMatrix[i, j] <- 1 - length(which(featureMatrixBinary[i, ] & featureMatrixBinary[j, ])) / length(which(featureMatrixBinary[i, ] | featureMatrixBinary[j, ]))
          #distanceMatrix[i, j] <- 1 - sum(featureMatrixBinary[i, ] & featureMatrixBinary[j, ]) / sum(featureMatrixBinary[i, ] | featureMatrixBinary[j, ])
          
          intersectionCount <- sum(featureIndeces[[i]] %in% featureIndeces[[j]])
          distanceMatrix[i, j] <- 1 - intersectionCount / (length(featureIndeces[[i]]) + length(featureIndeces[[j]]) - intersectionCount)
        }
      }
    },
    "Jaccard (intensity-weighted)"={
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      
      distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            distanceMatrix[i, j] <- 0
            next
          }
          
          intersection <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
          
          onlyI <- setdiff(x = featureIndeces[[i]], y = intersection)
          onlyJ <- setdiff(x = featureIndeces[[j]], y = intersection)
          sumOnlyI <- sum(featureMatrix[i, onlyI])
          sumOnlyJ <- sum(featureMatrix[j, onlyJ])
          
          relevance  <- apply(X = as.matrix(featureMatrix[, intersection]), MARGIN = 2, FUN = function(x) {max(x[i], x[j])})
          similarity <- 1 - abs(featureMatrix[i, intersection] - featureMatrix[j, intersection]) / relevance
          intersectionSum <- sum(relevance * similarity)
          unionSum <- intersectionSum + sumOnlyI + sumOnlyJ
          
          distanceMatrix[i, j] <- 1 - intersectionSum / unionSum
        }
      }
    },
    "Jaccard (intensity-weighted map)"={
      ## intersection i and j: featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
      ## diff i - j: featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]]
      ## diff j - i: featureIndeces[[j]][!featureIndeces[[j]] %in% featureIndeces[[i]]]
      ## union i or j: c(featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]], featureIndeces[[j]])
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      intensityMapping <- function(x){
        if(x == 0){
          ## x = 0
          return(0)
        } else if(x < 0.2){
          ## 0 <= x < 0.2
          return(0.01)
        } else if(x < 0.4){
          ## 0.2 <= x < 0.4
          return(0.2)
        } else {
          ## 0.4 <= x <= Inf
          return(1)
        }
      }
      intensityMapping <- Vectorize(FUN = intensityMapping, vectorize.args = "x")
      
      distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            distanceMatrix[i, j] <- 0
            next
          }
          
          intersection <- featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
          
          if(length(intersection) == 0){
            distance <- 1#sumOnlyI + sumOnlyJ
          } else {
            onlyI <- featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]]
            onlyJ <- featureIndeces[[j]][!featureIndeces[[j]] %in% featureIndeces[[i]]]
            
            if(length(onlyI) == 0){
              sumOnlyI <- 0
            } else {
              sumOnlyI <- sum(sapply(X = featureMatrix[i, onlyI], FUN = intensityMapping))
            }
            if(length(onlyJ) == 0){
              sumOnlyJ <- 0
            } else {
              sumOnlyJ <- sum(sapply(X = featureMatrix[j, onlyJ], FUN = intensityMapping))
            }
            
            maxIntensity <- apply(X = as.matrix(featureMatrix[c(i, j), intersection]), MARGIN = 2, FUN = max)
            intersectionSum <- sum(sapply(X = maxIntensity, FUN = intensityMapping))
            unionSum <- intersectionSum + sumOnlyI + sumOnlyJ
            distance <- 1 - intersectionSum / unionSum
          }
          
          distanceMatrix[i, j] <- distance
        }
      }
    },
    "Similarity (intensity-weighted)"={
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      
      similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            similarityMatrix[i, j] <- 0
            next
          }
          
          intersection <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
          relevance  <- apply(X = as.matrix(featureMatrix[, intersection]), MARGIN = 2, FUN = function(x) {max(x[i], x[j])})
          similarity <- 1 - abs(featureMatrix[i, intersection] - featureMatrix[j, intersection]) / relevance
          intersectionSum <- sum(relevance * similarity)
          
          similarityMatrix[i, j] <- intersectionSum
        }
      }
      distanceMatrix <- max(similarityMatrix) - similarityMatrix
    },
    "Jaccard (intensity-fragment-count-weighted)"={
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      fragmentFrequency <- dataList$fragmentFrequency
      
      distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            distanceMatrix[i, j] <- 0
            next
          }
          
          intersection <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
          
          onlyI <- setdiff(x = featureIndeces[[i]], y = intersection)
          onlyJ <- setdiff(x = featureIndeces[[j]], y = intersection)
          sumOnlyI <- sum(featureMatrix[i, onlyI] * fragmentFrequency[onlyI])
          sumOnlyJ <- sum(featureMatrix[j, onlyJ] * fragmentFrequency[onlyJ])
          
          relevance  <- apply(X = as.matrix(featureMatrix[c(i, j), intersection]), MARGIN = 2, FUN = max)
          similarity <- 1 - abs(featureMatrix[i, intersection] - featureMatrix[j, intersection]) / relevance
          relevance2 <- fragmentFrequency[intersection]
          intersectionSum <- sum(relevance * relevance2 * similarity)
          unionSum <- intersectionSum + sumOnlyI + sumOnlyJ
          
          distanceMatrix[i, j] <- 1 - intersectionSum / unionSum
        }
      }
    },
    "Similarity (intensity-fragment-count-weighted)"={
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      fragmentFrequency <- dataList$fragmentFrequency
      
      similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            similarityMatrix[i, j] <- 0
            next
          }
          
          intersection <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
          relevance  <- apply(X = as.matrix(featureMatrix[, intersection]), MARGIN = 2, FUN = function(x) {max(x[i], x[j])})
          similarity <- 1 - abs(featureMatrix[i, intersection] - featureMatrix[j, intersection]) / relevance
          relevance2 <- fragmentFrequency[intersection]
          intersectionSum <- sum(relevance * relevance2 * similarity)
          
          similarityMatrix[i, j] <- intersectionSum
        }
      }
      distanceMatrix <- max(similarityMatrix) - similarityMatrix
    },
    "Jaccard (fragment-count-weighted)"={
      featureIndexMatrix <- dataList$featureIndexMatrix[filter, ]
      fragmentFrequency <- dataList$fragmentFrequency
      
      counter <<- 0
      distanceMatrix <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(x)
        {  
          counter <<- counter + 1
          if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", counter, " / ", numberOfPrecursors, sep = ""))
          intersectionSum <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(y){  sum(fragmentFrequency[x[x %in% y]], na.rm = TRUE)  })
          unionSum        <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(y){  sum(fragmentFrequency[c(x[!x %in% y], y)], na.rm = TRUE)  })
          
          1 - intersectionSum / unionSum
        }
      )
    },
    "Manhatten"={
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      
      ## Rasmussen 2008
      distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            distanceMatrix[i, j] <- 0
            next
          }
          
          intersection <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
          intersectionSum <- sum(abs(featureMatrix[i, intersection] - featureMatrix[j, intersection]))
          
          onlyI <- setdiff(x = featureIndeces[[i]], y = intersection)
          onlyJ <- setdiff(x = featureIndeces[[j]], y = intersection)
          sumOnlyI <- sum(featureMatrix[i, onlyI])
          sumOnlyJ <- sum(featureMatrix[j, onlyJ])
          
          distanceMatrix[i, j] <- intersectionSum + sumOnlyI + sumOnlyJ
        }
      }
    },
    "NDP"={
      featureIndeces <- dataList$featureIndeces[filter]
      featureMatrix <- dataList$featureMatrix[filter, ]
      
      ## Gaquerel 2015: standard normalized dot product (NDP) / cosine correlation
      similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
      for(i in 1:numberOfPrecursors){
        if(progress)  incProgress(amount = 1 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
        for(j in 1:numberOfPrecursors){
          if(i == j){
            similarityMatrix[i, j] <- 0
            next
          }
          
          intersection <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
          intersectionSum <- sum(featureMatrix[i, intersection] * featureMatrix[j, intersection])^2
          iSum <- sum(featureMatrix[i, featureIndeces[[i]]] * featureMatrix[i, featureIndeces[[i]]])
          jSum <- sum(featureMatrix[j, featureIndeces[[j]]] * featureMatrix[j, featureIndeces[[j]]])
          
          similarityMatrix[i, j] <- intersectionSum / (iSum * jSum)
        }
      }
      distanceMatrix <- max(similarityMatrix) - similarityMatrix
    },
    stop(paste("Unknown distance (", distance, ")!", sep = ""))
  )
  
  return(distanceMatrix)
}
calculateCluster <- function(dataList, filter, distanceMatrix, method, progress = FALSE){
  numberOfPrecursorsFiltered <- length(filter)
  ##########################################
  ## compute gui stuff
  
  if(progress)  incProgress(amount = 0, detail = "Clustering")
  ## compute and annotate cluster
  dist <- stats::as.dist(m = distanceMatrix)
  cluster <- hclust(d = dist, method = method)
  cluster$labels <- dataList$precursorLabels[filter]
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  leafHeightSpacing <- 0.04
  
  #opt <- order.optimal(dist = dist, merge = cluster$merge)
  #cluster$merge <- opt$merge
  #cluster$order <- opt$order
  
  ## compute (transitive) cluster members, cluster positions, and leaf heights
  if(progress)  incProgress(amount = 0.5, detail = "Analyze cluster")
  
  resultObj <- analyzeTreeFromRoot(dataList, cluster = cluster, filter)
  innerNodeMembersTree <- resultObj$innerNodeMembersTree
  innerNodeMembers <- resultObj$innerNodeMembers
  #innerNodeFeaturesBinary <- resultObj$innerNodeFeaturesBinary
  innerNodeFeaturesIntersection <- resultObj$innerNodeFeaturesIntersection
  innerNodeFeaturesUnion <- resultObj$innerNodeFeaturesUnion
  innerNodeFeaturesCounts <- resultObj$innerNodeFeaturesCounts
  innerNodeFeaturesMeanAbundance <- resultObj$innerNodeFeaturesMeanAbundance
  innerNodeFeaturesIntersectionCounter <- resultObj$innerNodeFeaturesIntersectionCounter
  innerNodeFeaturesUnionCounter <- resultObj$innerNodeFeaturesUnionCounter
  innerNodePosition <- resultObj$innerNodePosition
  leafHeights <- resultObj$leafHeights
  
  ## dendrogram leaf ends for normal plot
  #leafHeights <- leafHeights - leafHeightSpacing
  leafHeights <- rep(x = 0, times = length(leafHeights))
  
  ## compute x- and y-coordinates and point-labels
  coordinatesX <- unlist(c(innerNodePosition, match(x = 1:numberOfPrecursorsFiltered, table = cluster$order)))
  coordinatesY <- unlist(c(cluster$height, leafHeights))
  labels <- unlist(c(1:numberOfInnerNodes, -(1:numberOfPrecursorsFiltered)))
  text <- unlist(c(innerNodeFeaturesIntersectionCounter, apply(X = dataList$featureMatrix[filter, ], MARGIN = 1, FUN = function(x) sum(x != 0))))
  
  ##########################################
  ## box
  if(progress)  incProgress(amount = 0.5, detail = "Boxing")
  filterList <- list()
  ## filter
  filterList$filter <- filter
  filterList$numberOfPrecursorsFiltered <- numberOfPrecursorsFiltered
  ## cluster
  
  filterList$innerNodeMembersTree <- innerNodeMembersTree
  filterList$innerNodeMembers <- innerNodeMembers
  #filterList$innerNodeFeaturesBinary <- innerNodeFeaturesBinary
  filterList$innerNodeFeaturesIntersection <- innerNodeFeaturesIntersection
  filterList$innerNodeFeaturesUnion <- innerNodeFeaturesUnion
  filterList$innerNodeFeaturesCounts <- innerNodeFeaturesCounts
  filterList$innerNodeFeaturesMeanAbundance <- innerNodeFeaturesMeanAbundance
  filterList$innerNodeFeaturesIntersectionCounter <- innerNodeFeaturesIntersectionCounter
  filterList$innerNodeFeaturesUnionCounter <- innerNodeFeaturesUnionCounter
  filterList$cluster <- cluster
  ## poi
  filterList$poiCoordinatesX <- coordinatesX
  filterList$poiCoordinatesY <- coordinatesY
  filterList$poiText <- text
  filterList$poiLabels <- labels
  
  if(progress)  setProgress(1)
  
  return(filterList)
}
#########################################################################################
## data fetching
getMS2spectrum <- function(dataList, clusterDataList, clusterLabel){
  if(clusterLabel < 0){
    ###############################################
    ## leaf
    precursorIndex <- clusterDataList$filter[[-clusterLabel]]
    precursorMass  <- as.numeric(dataList$dataFrame$"m/z"[[precursorIndex]])
    adduct <- dataList$dataFrame$"Adduct ion name"[[precursorIndex]]
    neutralMassCorrection <- NA
    ionMode <- NA
    switch(adduct,
           "[M-H]-"={
             neutralMassCorrection <- 1.008
             ionMode <- -1
           },
           "[M+H]+"={
             neutralMassCorrection <- -1.008
             ionMode <- 1
           },
           stop(paste("Unknown adduct (", adduct, ")!", sep = ""))
    )
    neutralMass <- precursorMass + neutralMassCorrection
    features <- dataList$featureIndeces[[precursorIndex]]
    fragmentsX <- dataList$fragmentMasses[features]
    fragmentsY <- as.numeric(dataList$featureMatrix[precursorIndex, features])
    
    fragmentsPositive <- fragmentsX > 0
    fragmentsPositiveX <- fragmentsX[fragmentsPositive]
    fragmentsPositiveY <- fragmentsY[fragmentsPositive]
    fragmentStrings <- paste(fragmentsPositiveX, fragmentsPositiveY, sep = " ", collapse = "; ")
    
    ##http://msbi.ipb-halle.de/MetFragBeta/LandingPage.jspx?limit=1000&ionmode=-1&database=pubchem&mzppm=7&mzabs=0.005&mass=448.468&formula=C16H20N2O9S2&mzabs=0.05&peaks=130.0655 288214.8119 ; 207.0589 422771.0127 ; 208.0622  87002.3217 ; 210.1334   2674.1707 ; 351.1016  27580.9393 ; 369.1115 739357.5045 ; 370.1148 143864.9611 ; 385.1094   5971.8328 ; 391.0937 337133.4536 ; 392.1025  40126.6888 ; 407.0678   3095.0322 ; 449.0690  37952.2515 
    landingPageUrl <- paste(sep = "",
      "http://msbi.ipb-halle.de/MetFragBeta/LandingPage.jspx?",
      "mass=", neutralMass, "&",
      "formula=", "", "&",
      "ionmode=", ionMode, "&",
      #"limit=", "1000", "&",
      "database=", "pubchem", "&",
      #"mzppm=", "7", "&"
      #"mzabs=", "0.005", "&",
      "peaks=", fragmentStrings
    )
    
    infoText <- paste("Precursor ''", clusterDataList$cluster$labels[[-clusterLabel]], "'' has ", length(fragmentsX), " fragments", sep = "")
    #cat(paste(paste(sort(fragmentsX), collapse = ", "), "\n", sep = ""))
    if(!is.na(neutralMass)){
      #writeClipboard(landingPageUrl, format = 1)
      landingPageUrlForLink <- landingPageUrl
    }
    else
      landingPageUrlForLink <- NULL
    
    columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
    dataFrame     <- data.frame(dataList$dataFrameMeasurements[precursorIndex, columnNames])
    colnames(dataFrame) <- columnNames
    
    featureMatrix <- data.frame(dataList$featureMatrix[precursorIndex, features])
    featureMatrix <- t(featureMatrix)
    
    precursorSet <- precursorIndex
    dataFrame     <- cbind(dataFrame, featureMatrix)
  } else {
    ###############################################
    ## inner node
    clusterIndex <- clusterLabel
    clusterMembers <- sort(clusterDataList$innerNodeMembers[[clusterIndex]])
    
    #features <- clusterDataList$innerNodeFeaturesBinary[[clusterIndex]]
    featuresIntersection <- clusterDataList$innerNodeFeaturesIntersection[[clusterIndex]]
    featuresUnion <- clusterDataList$innerNodeFeaturesUnion[[clusterIndex]]
    fragmentsX <- dataList$fragmentMasses[featuresIntersection]
    fragmentsY <- apply(X = data.matrix(dataList$featureMatrix[clusterMembers, featuresIntersection]), MARGIN = 2, FUN = mean)
    
    infoText <- paste("Cluster ''", paste(dataList$precursorLabels[clusterMembers], sep = ", ", collapse = ", "), "'' has ", length(fragmentsX), " fragments in common", sep = "")
    #cat(paste(paste(sort(fragmentsX), collapse = ", "), "\n", sep = ""))
    landingPageUrlForLink <- NULL
    
    ###############################################
    ## table data
    ## TODO
    columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
    dataFrame     <- dataList$dataFrameMeasurements[clusterMembers, columnNames]
    
    maximumNumberOfFeatures <- 100
    #featureMatrix <- dataList$featureMatrix[clusterMembers, featuresUnion]
    featureMatrix <- data.frame(dataList$featureMatrix[clusterMembers, featuresIntersection])
    colnames(featureMatrix) <- dataList$fragmentMasses[featuresIntersection]
    if(ncol(featureMatrix) > maximumNumberOfFeatures)
      featureMatrix <- featureMatrix[, 1:maximumNumberOfFeatures]
    
    precursorSet <- clusterMembers
    dataFrame <- cbind(dataFrame, featureMatrix)
  }
  
  ## selected data
  order <- order(fragmentsX)
  fragmentsX <- fragmentsX[order]
  fragmentsY <- fragmentsY[order]
  
  resultObj <- list()
  resultObj$fragmentsX <- fragmentsX
  resultObj$fragmentsY <- fragmentsY
  resultObj$infoText <- infoText
  resultObj$landingPageUrl <- landingPageUrlForLink
  resultObj$precursorSet <- precursorSet
  resultObj$dataFrame <- dataFrame
  
  return(resultObj)
}
getMS2spectrumOfPrecursor <- function(dataList, precursorIndex){
  featureIndeces <- dataList$featureIndeces[precursorIndex, ]
  featureMasses  <- dataList$fragmentMasses[featureIndeces]
  featureValues  <- dataList$featureMatrix [precursorIndex, featureIndeces]
  
  resultObj <- list()
  resultObj$fragmentMasses <- featureMasses
  resultObj$fragmentAbundances <- featureValues
  
  return(resultObj)
}
#########################################################################################
## plotting
colorLabels <- function(labels, clusterMembers, color, labelsToRemove = NULL){
  colLab <- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
      
      ## $label "New Hampshire"
      ## $members 1
      ## $height 0
      ## $leaf TRUE
      ## $class "dendrogram"
      
      nodeLabelHere <- a$label
      nodeIndexHere <- match(x = nodeLabelHere, table = labels)
      
      if(length(na.omit(match(x = nodeIndexHere, table = clusterMembers))) > 0)
        attr(n, "nodePar") <- c(a$nodePar, lab.col = color) # change the node color
      
      if(!is.null(labelsToRemove))
        if(length(na.omit(match(x = nodeLabelHere, table = labelsToRemove))) > 0)
          attr(n, "label") <- "" # clear label
    }
    return(n)
  }
  
  return(colLab)
}
calcClusterPlotDendrogram <- function(dataList, filter, clusterDataList, annoPresentAnnotationsList, annoPresentColorsList, distanceMeasure, nodeIndex = NULL, nodeIndeces = NULL, xInterval = NULL){
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## cluster
  par(mar=c(7,4,2,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  dend <- as.dendrogram(clusterDataList$cluster)
  
  ## remove labels below the y-axis
  rightMostInvisibleLabelIndex <- floor(xInterval[[1]] - (xInterval[[2]] - xInterval[[1]]) * 0.04)
  if(rightMostInvisibleLabelIndex > 0){
    labelsToRemove <- clusterDataList$cluster$labels[clusterDataList$cluster$order][1:rightMostInvisibleLabelIndex]
    
    colLab <- colorLabels(clusterDataList$cluster$labels, NULL, NULL, labelsToRemove)
    dend <- dendrapply(dend, colLab)
  }
  
  ## color labels for multiple sub-roots
  if(!is.null(nodeIndeces)){
    clusterMembers <- c(unlist(clusterDataList$innerNodeMembersTree[nodeIndeces[nodeIndeces > 0]]), -nodeIndeces[nodeIndeces < 0])
    
    colLab <- colorLabels(clusterDataList$cluster$labels, clusterMembers, 'green')
    dend <- dendrapply(dend, colLab)
  }
  ## color labels for single sub-root
  if(!is.null(nodeIndex)){
    if(nodeIndex > 0){
      clusterMembers <- clusterDataList$innerNodeMembersTree[[nodeIndex]]
    } else {
      clusterMembers <- c(-nodeIndex)
    }
    
    colLab <- colorLabels(clusterDataList$cluster$labels, clusterMembers, 'blue')
    dend <- dendrapply(dend, colLab)
  }
  
  ## plot
  plot(x = dend, xlab = "", ylab = distanceMeasure, main = "Precursor cluster dendrogram", sub = "", xlim = xInterval)
  
  #xIntervalSize <<- xInterval[[2]] - xInterval[[1]]
  
  ## color tree for annotations
  resultObj <- analyzeTreeFromRoot(dataList, cluster = clusterDataList$cluster, filter)
  innerNodeFeaturesAnnotations <- resultObj$innerNodeFeaturesAnnotations
  rootIndex <- length(clusterDataList$cluster$height)
  innerNodeColors <- lapply(X = innerNodeFeaturesAnnotations, FUN = function(x){ unlist(lapply(X = x, FUN = function(y){ annoPresentColorsList[[match(x = y, table = annoPresentAnnotationsList)]] })) })
  
  for(idx in 1:length(resultObj$innerNodeFeaturesIntersection)){
    if(length(innerNodeFeaturesAnnotations) < idx)
      print(paste(idx, "", "", length(resultObj$innerNodeFeaturesIntersection[[idx]])))
    else
      print(paste(idx, length(innerNodeFeaturesAnnotations[[idx]]), length(innerNodeColors[[idx]]), length(resultObj$innerNodeFeaturesIntersection[[idx]])))
  }
  
  colorSubTreeForAnnotations(cluster = clusterDataList$cluster, index = rootIndex, innerNodeAnnotations = innerNodeFeaturesAnnotations, innerNodeColors = innerNodeColors, parentIndex = NULL)
  
  if(FALSE){
  ## color tree for multiple sub-roots
  if(!is.null(nodeIndeces)){
    for(nodeIdx in nodeIndeces[nodeIndeces > 0]){
      ## sub tree coloring
      clusterMembers <- clusterDataList$innerNodeMembersTree[[nodeIdx]]
      a2r_counter <<- min(match(x = clusterMembers, table = clusterDataList$cluster$order)) - 1
      colorSubTree(cluster = clusterDataList$cluster, index = nodeIdx, col = "green")
    }
  }
  ## color tree for single sub-root
  if(!is.null(nodeIndex)){
    if(nodeIndex > 0){
      ## sub tree coloring
      clusterMembers <- clusterDataList$innerNodeMembersTree[[nodeIndex]]
      a2r_counter <<- min(match(x = clusterMembers, table = clusterDataList$cluster$order)) - 1
      colorSubTree(cluster = clusterDataList$cluster, index = nodeIndex, col = "blue")
    }
  }
  }
  
  ## POI's
  pointSizes <- rep(x = 1., times = length(clusterDataList$poiCoordinatesX))
  pointSizes2 <- rep(x = 2/3, times = length(clusterDataList$poiCoordinatesX))
  pointColors <- rep(x = "black", times = length(clusterDataList$poiCoordinatesX))
  pointColors2 <- rep(x = "gray", times = length(clusterDataList$poiCoordinatesX))
  if(!is.null(nodeIndex)){
    if(nodeIndex < 0)
      idx <- as.integer(length(pointSizes) / 2) - nodeIndex
    else
      idx <- nodeIndex
    
    pointSizes[[idx]] <- 2.
    pointSizes2[[idx]] <- 4/3
    pointColors[[idx]] <- "blue"
    pointColors2[[idx]] <- "blue"
  }
  
  points(x = clusterDataList$poiCoordinatesX, y = clusterDataList$poiCoordinatesY, col = pointColors, pch=19, cex=pointSizes)
  points(x = clusterDataList$poiCoordinatesX, y = clusterDataList$poiCoordinatesY, col = pointColors2, pch=19, cex=pointSizes2)
  graphics::text(  x = clusterDataList$poiCoordinatesX, y = clusterDataList$poiCoordinatesY + 0.02, labels = clusterDataList$poiText, pos = 4)
}
calcClusterPlotHeatmap <- function(dataList, filter, clusterDataList, xInterval = NULL){
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## heatmap
  columnsOfInterest <- c(
    dataList$dataMeanColumnNameFunctionFromName(filter$groups[[1]]), dataList$dataMeanColumnNameFunctionFromName(filter$groups[[2]]), 
    dataList$lfcColumnNameFunctionFromName(filter$groups[[1]], filter$groups[[2]])
  )
  
  par(mar=c(0,4,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right) ## c(title, axis, label)
  
  if(length(columnsOfInterest) == 0){
    plot.new()
  } else {
    colorOne <- dataList$colorMatrixDataFrame[clusterDataList$filter, columnsOfInterest[[1]]][clusterDataList$cluster$order]
    colorTwo <- dataList$colorMatrixDataFrame[clusterDataList$filter, columnsOfInterest[[2]]][clusterDataList$cluster$order]
    colorLFC <- dataList$colorMatrixDataFrame[clusterDataList$filter, columnsOfInterest[[3]]][clusterDataList$cluster$order]
    
    plot(x = c(1, clusterDataList$numberOfPrecursorsFiltered), y = c(0, 3), type= "n", xlab = "", ylab = "", axes = FALSE, xlim = xInterval, ylim = c(0, 3))
    for(i in 1:clusterDataList$numberOfPrecursorsFiltered){
      rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 2, ytop = 3, col = colorLFC[[i]], border = NA)
      rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 1, ytop = 2, col = colorOne[[i]], border = NA)
      rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 0, ytop = 1, col = colorTwo[[i]], border = NA)
    }
    axis(side = 2, at = c(0.5, 1.5, 2.5), labels = c(filter$groups[[2]], filter$groups[[1]], "LFC"), las = 2, tick = TRUE)
  }
}
calcClusterPlotHeatmapLegend <- function(dataList){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  legend_imageAbs <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = dataList$absMax,           to =  0,                         length.out = 100), map = dataList$colorMapAbsoluteData ), nrow=1)))
  legend_imageLFC <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = dataList$logFoldChangeMax, to = -dataList$logFoldChangeMax, length.out = 100), map = dataList$colorMapLogFoldChange), nrow=1)))
  
  plot.new()
  plot.window(xlim = c(0, 3), ylim = c(0, 1))
  epsilon <- 0.05
  rasterImage(image = legend_imageAbs, xleft = 0, ybottom = 0, xright = 1, ytop = 0.5-epsilon)
  graphics::text(x = 2, y = seq(0,0.5-epsilon,l=5), labels = format(x = seq(0, dataList$absMax, l=5), scientific = TRUE, digits = 0))
  
  rasterImage(image = legend_imageLFC, xleft = 0, ybottom = 0.5+epsilon, xright = 1, ytop = 1)
  graphics::text(x = 2, y = seq(0.5+epsilon,1,l=5), labels = format(x = seq(-dataList$logFoldChangeMax, dataList$logFoldChangeMax, l=5), digits = 0))
}
calcClusterPlotMS2 <- function(dataList, fragmentsX = c(), fragmentsY = c(), fragmentsX_02 = NULL, fragmentsY_02 = NULL, xInterval = NULL, selectedFragmentIndex = NULL){
  ####################
  ## fragment spectrum
  if(is.null(xInterval))
    xInterval <- c(dataList$minimumMass, dataList$maximumMass)
  
  ## abundances greater one
  fragmentsY[fragmentsY > 1] <- 1
  
  if(is.null(fragmentsX_02)){
    yInterval <- c(0, 1)
    nodeColors <- rep(x = "black", times = length(fragmentsX))
    dataX <- fragmentsX
    dataY <- fragmentsY
    yTickPositions <- c(0, 0.25, 0.5, 0.75, 1)
    yTickLabels <- c(0, "", 0.5, "", 1)
  } else {
    nodeColors <- rep(x = "black", times = length(fragmentsX) + length(fragmentsX_02))
    
    if(is.null(fragmentsX)){
      yInterval <- c(0, 1)
      dataX <- fragmentsX_02
      dataY <- fragmentsY_02
      yTickPositions <- c(0, 0.25, 0.5, 0.75, 1)
      yTickLabels <- c(0, "", 0.5, "", 1)
    } else {
      yInterval <- c(-1, 1)
      dataX <- c(fragmentsX, fragmentsX_02)
      dataY <- c(fragmentsY, -fragmentsY_02)
      yTickPositions <- c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
      yTickLabels <- c(1, "", 0.5, "", 0, "", 0.5, "", 1)
    }
  }
  
  ## node selection
  if(!is.null(selectedFragmentIndex))
    nodeColors[[selectedFragmentIndex]] <- "green"
  
  par(mar=c(5,4,3,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  plot(x = dataX, y = dataX, ylab = "Relative abundance", xlab = "", xlim = xInterval, ylim = yInterval, xaxt='n', yaxt='n', col = nodeColors)
  axis(side = 2, at = yTickPositions, labels = yTickLabels)
  axis(side = 3)
  mtext(side = 3, "m/z", line = 2)
  
  if(!is.null(fragmentsX_02)){
    xIntervalSize <- xInterval[[2]] - xInterval[[1]]
    xl <- xInterval[[1]] - xIntervalSize
    xr <- xInterval[[2]] + xIntervalSize
    segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  }
    
  tickPositions <- dataX
  #minimumTickLabelShift <- (dataList$maximumMass - dataList$minimumMass) / 75
  #if(length(tickPositions) > 1)
  #  for(i in 2:length(tickPositions))
  #    if(tickPositions[[i]] - tickPositions[[i - 1]] < minimumTickLabelShift)
  #      tickPositions[[i]] <- tickPositions[[i - 1]] + minimumTickLabelShift
  if(length(dataX) > 0){
    ## axis
    axis(side = 1, at = dataX, labels = FALSE, las = 2)
    axis(side = 1, at = tickPositions, labels = dataX, las = 2, tick = FALSE)
    
    ## points
    points(x = dataX, y = dataY, col = nodeColors, type = "h", lwd=4)
    
    pointSizes <- rep(x = 1., times = length(dataX))
    if(!is.null(selectedFragmentIndex))
      pointSizes[[selectedFragmentIndex]] <- 2.
    pointSizes2 <- rep(x = 2/3, times = length(dataX))
    if(!is.null(selectedFragmentIndex))
      pointSizes2[[selectedFragmentIndex]] <- 4/3
    pointColors <- rep(x = "black", times = length(dataX))
    if(!is.null(selectedFragmentIndex))
      pointColors[[selectedFragmentIndex]] <- "green"
    pointColors2 <- rep(x = "gray", times = length(dataX))
    if(!is.null(selectedFragmentIndex))
      pointColors2[[selectedFragmentIndex]] <- "green"
    
    points(x = dataX, y = dataY, col = pointColors, pch=19, cex=pointSizes)
    points(x = dataX, y = dataY, col = pointColors2, pch=19, cex=pointSizes2)
    
    #points(x = dataX, y = dataY, col = nodeColors, pch=19, cex=0.7)
    #for(idx in 1:length(dataX)){
    #  points(x = dataX[[idx]], y = dataY[[idx]], col = "red", type = "h", lwd=4)
    #  points(x = dataX[[idx]], y = dataY[[idx]], col = "red", pch=19, cex=0.7)
    #  #text(x = dataX[[idx]], y = dataY[[idx]], labels = dataX[[idx]], pos = 4)
    #}
  }
}
#calcPCA <- function(dataList, filter, scaling, meanCentered, logTransform, method){
calcPCA <- function(dataList, filter, scaling, logTransform){
  dataFrame <- dataList$dataFrameMeasurements[filter$filter, dataList$dataColumnsNameFunctionFromNames(filter$groups)]
  dataFrame <- t(dataFrame)
  
  if(logTransform){
    dataFrame[dataFrame < 1] <- 1
    dataFrame <- log2(dataFrame)
  }
  
  switch(scaling,
    "None"={
      ## do nothing
    },
    "Mean center"={
      ## subtract mean
      dataFrame <- t(apply(X = dataFrame, MARGIN = 1, FUN = function(x){ x - mean(x) }))
    },
    "Autoscaling (unit variance)"={
      ## subtract mean and divide by variance
      dataFrame <- t(apply(X = dataFrame, MARGIN = 1, FUN = function(x){ (x - mean(x)) / var(x = x) }))
    },
    "Pareto"={
      ## subtract mean and divide by sqrt of variance
      dataFrame <- t(apply(X = dataFrame, MARGIN = 1, FUN = function(x){ (x - mean(x)) / sqrt(var(x = x)) }))
    },
    stop(paste("Unknown scaling (", scaling, ")!", sep = ""))
  )
  
  ## "FactoMineR" package
  #pca = PCA(dataFrame, graph = FALSE, scale.unit = unitVariance)
  ## scores   : pca3$ind$coord
  ## loadings : pca3$var$coord
  ## variances: pca3$eig
  
  ## "pcaMethods" package
  #pca <- pca(object = dataFrame, method = method, nPcs = 5, scale = scaling, center = meanCentered)
  #pca <- pca(object = dataFrame, nPcs = 5, scale = "none")
  #
  #returnObj <- list()
  #returnObj$scores   <- pca@scores
  #returnObj$loadings <- pca@loadings
  #returnObj$variance <- pca@sDev
  
  ## "stats" package
  pca <- prcomp(x = dataFrame, retx = TRUE, center = FALSE, scale. = FALSE)
  returnObj <- list()
  returnObj$scores   <- pca$x
  returnObj$loadings <- pca$rotation
  returnObj$variance <- pca$sdev
  
  return(returnObj)
}
plotPCAscores <- function(pcaObj, dataList, filter, pcaDimensionOne, pcaDimensionTwo, xInterval = NULL, yInterval = NULL){
  palette <- colorPalette()
  colors <- palette[unlist(lapply(X = filter$groups, FUN = function(x){ 
    groupIdx <- dataList$groupIdxFromGroupName(x)
    rep(x = groupIdx, times = length(dataList$dataColumnsNameFunctionFromName(x)))
  }))]
  
  #dataDimOne <- pca$ind$coord[, pcaDimensionOne]
  #dataDimTwo <- pca$ind$coord[, pcaDimensionTwo]
  dataDimOne <- pcaObj$scores[, pcaDimensionOne]
  dataDimTwo <- pcaObj$scores[, pcaDimensionTwo]
  xAxisLabel <- paste("t_", pcaDimensionOne, sep = "")
  yAxisLabel <- paste("t_", pcaDimensionTwo, sep = "")
  
  xMin <- min(dataDimOne)
  xMax <- max(dataDimOne)
  yMin <- min(dataDimTwo)
  yMax <- max(dataDimTwo)
  
  if(is.null(xInterval))
    xInterval <- c(xMin, xMax)
  if(is.null(yInterval))
    yInterval <- c(yMin, yMax)
  
  par(mar=c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  plot(x = dataDimOne, y = dataDimTwo, xlim = xInterval, ylim = yInterval, xlab = xAxisLabel, ylab = yAxisLabel, main = "Scores", col = colors, pch=19, cex=1.)
  
  ## axis
  xInt <- xMax - xMin
  yInt <- yMax - yMin
  xl <- xMin - xInt
  xr <- xMax + xInt
  yl <- yMin - yInt
  yr <- yMax + yInt
  segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  segments(x0 = 0, x1 = 0, y0 = yl, y1 = yr, col = "black", lwd = 1)
  
  labels <- dataList$dataColumnsNameFunctionFromNames(filter$groups)
  graphics::text(x = dataDimOne, y = dataDimTwo, labels = labels, pos = 4)
}
plotPCAloadings <- function(pcaObj, dataList, filter, pcaDimensionOne, pcaDimensionTwo, pcaLoading = NULL, pcaLoadingSet = NULL, xInterval = NULL, yInterval = NULL){
  #rPlot(Dim.1 ~ Dim.2, data = pca$var$coord[, c(1, 2)])
  
  #varianceOne <- format(x = pca$eig[[pcaDimensionOne, 2]], digits = 3)
  #varianceTwo <- format(x = pca$eig[[pcaDimensionTwo, 2]], digits = 3)
  varianceOne <- format(x = pcaObj$variance[[pcaDimensionOne]], digits = 3)
  varianceTwo <- format(x = pcaObj$variance[[pcaDimensionTwo]], digits = 3)
  
  xAxisLabel  <- paste("p_", pcaDimensionOne, " (", varianceOne, "%)", sep = "")
  yAxisLabel  <- paste("p_", pcaDimensionTwo, " (", varianceTwo, "%)", sep = "")
  
  #dataDimOne <- pca$var$coord[, pcaDimensionOne]
  #dataDimTwo <- pca$var$coord[, pcaDimensionTwo]
  dataDimOne <- pcaObj$loadings[, pcaDimensionOne]
  dataDimTwo <- pcaObj$loadings[, pcaDimensionTwo]
  
  xMin <- min(dataDimOne)
  xMax <- max(dataDimOne)
  yMin <- min(dataDimTwo)
  yMax <- max(dataDimTwo)
  
  if(is.null(xInterval))
    xInterval <- c(xMin, xMax)
  if(is.null(yInterval))
    yInterval <- c(yMin, yMax)
  
  pointColors <- rep(x = "black", times = length(dataDimOne))
  if(sum(pcaLoadingSet) > 0)
    pointColors[pcaLoadingSet] <- "green"
  if(!is.null(pcaLoading))
    pointColors[[pcaLoading]] <- "blue"
  pointColors2 <- rep(x = "gray", times = length(dataDimOne))
  if(sum(pcaLoadingSet) > 0)
    pointColors2[pcaLoadingSet] <- "green"
  if(!is.null(pcaLoading))
    pointColors2[[pcaLoading]] <- "blue"
  
  pointSizes <- rep(x = 1., times = length(dataDimOne))
  if(!is.null(pcaLoading))
    pointSizes[[pcaLoading]] <- 2.
  pointSizes2 <- rep(x = 2/3, times = length(dataDimOne))
  if(!is.null(pcaLoading))
    pointSizes2[[pcaLoading]] <- 4/3
  
  if(sum(pcaLoadingSet) > 0){
    dataDimOne <- c(dataDimOne[-pcaLoadingSet], dataDimOne[pcaLoadingSet])
    dataDimTwo <- c(dataDimTwo[-pcaLoadingSet], dataDimTwo[pcaLoadingSet])
    pointColors <- c(pointColors[-pcaLoadingSet], pointColors[pcaLoadingSet])
    pointColors2 <- c(pointColors2[-pcaLoadingSet], pointColors2[pcaLoadingSet])
    pointSizes <- c(pointSizes[-pcaLoadingSet], pointSizes[pcaLoadingSet])
    pointSizes2 <- c(pointSizes2[-pcaLoadingSet], pointSizes2[pcaLoadingSet])
  }
  if(!is.null(pcaLoading)){
    if(sum(pcaLoadingSet) == 0)
      newIdx <- pcaLoading
    else{
      newIdx <- match(x = pcaLoading, table = which(pcaLoadingSet)) + (length(dataDimOne) - sum(pcaLoadingSet))
    }
    
    dataDimOne <- c(dataDimOne[-newIdx], dataDimOne[newIdx])
    dataDimTwo <- c(dataDimTwo[-newIdx], dataDimTwo[newIdx])
    pointColors <- c(pointColors[-newIdx], pointColors[newIdx])
    pointColors2 <- c(pointColors2[-newIdx], pointColors2[newIdx])
    pointSizes <- c(pointSizes[-newIdx], pointSizes[newIdx])
    pointSizes2 <- c(pointSizes2[-newIdx], pointSizes2[newIdx])
  }
  
  par(mar=c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  #plot(x = dataDimOne, y = dataDimTwo, xlim = xInterval, ylim = yInterval, xlab = xAxisLabel, ylab = yAxisLabel, main = "Loadings", pch=19, cex=0.7, col = nodeColors)
  plot(x = NULL, y = NULL, xlim = xInterval, ylim = yInterval, xlab = xAxisLabel, ylab = yAxisLabel, main = "Loadings")
  points(x = dataDimOne, y = dataDimTwo, col = pointColors, pch=19, cex=pointSizes)
  points(x = dataDimOne, y = dataDimTwo, col = pointColors2, pch=19, cex=pointSizes2)
  
  ## axis
  xInt <- xMax - xMin
  yInt <- yMax - yMin
  xl <- xMin - xInt
  xr <- xMax + xInt
  yl <- yMin - yInt
  yr <- yMax + yInt
  segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  segments(x0 = 0, x1 = 0, y0 = yl, y1 = yr, col = "black", lwd = 1)
}
colorPalette <- function(){
  ## http://tools.medialab.sciences-po.fr/iwanthue/
  ## 20 colors
  ## or
  ## https://github.com/johnbaums/hues/blob/master/R/iwanthue.R
  ## library(colorBrewer)
  ## colorRampPalette(c("blue", "red"))( 4)
  ## palette(rainbow(6))
  
  #palette <- palette(c(
  palette <- c(
    rgb(184,88,184, maxColorValue=255),
    rgb(102,178,48, maxColorValue=255),
    rgb(220,64,59, maxColorValue=255),
    rgb(78,167,149, maxColorValue=255),
    rgb(117,79,33, maxColorValue=255),
    rgb(131,59,93, maxColorValue=255),
    rgb(116,159,202, maxColorValue=255),
    rgb(176,158,56, maxColorValue=255),
    rgb(114,119,221, maxColorValue=255),
    rgb(219,119,40, maxColorValue=255),
    rgb(72,101,46, maxColorValue=255),
    rgb(213,66,135, maxColorValue=255),
    rgb(70,91,112, maxColorValue=255),
    rgb(213,115,121, maxColorValue=255),
    rgb(91,76,141, maxColorValue=255),
    rgb(198,137,190, maxColorValue=255),
    rgb(98,175,100, maxColorValue=255),
    rgb(146,56,45, maxColorValue=255),
    rgb(207,78,219, maxColorValue=255),
    rgb(206,136,82, maxColorValue=255)
  )
  #))
  return(palette)
}
