
#source("https://bioconductor.org/biocLite.R")
#biocLite("mzR")
library("mzR")
#install.packages("squash")
library("squash")

#install.packages("FactoMineR")
library("FactoMineR")
#require(devtools)
#install_github('rCharts', 'ramnathv')
#library("rCharts")
#source("http://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")
library("pcaMethods")
#install.packages("cba")
#library("cba")
#install.packages("matrixStats")
library("matrixStats")
library("Matrix")
#install.packages("plotrix")
library("plotrix")
library("tools")

#########################################################################################
## constants
minimumProportionOfLeafs <- 0.75
minimumProportionToShowFragment <- 0.5

clusterNodePointSize0 <- 2/3
clusterNodePointSize1 <- 4/3
clusterNodePointSize2 <- 6/3
clusterNodePointSize3 <- 8/3

ms2StickPointSizeInitial <- 1.
ms2StickPointSizeInitialSmall <- 2/3.
ms2StickPointSizeEmph <- 1.5
ms2StickPointSizeEmphSmall <- 3/3.

#########################################################################################
## tree helper
analyzeTreeFromRoot <- function(dataList, cluster, filter){
  numberOfPrecursorsFiltered <- length(filter)
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  rootIndex <- length(cluster$height)
  
  ## create fields and compute stuff
  innerNodeHeightIncreasesHere <<- list()
  innerNodeMembersTreeHere <<- list()
  innerNodeMembersHere <<- list()
  innerNodeMembersInnerHere <<- list()
  innerNodeFeaturesIntersectionHere <<- list()
  innerNodeFeaturesAnnotations <<- list()
  innerNodeFeaturesUnionHere <<- list()
  #innerNodeFeaturesCountsHere <<- list()
  #innerNodeFeaturesCountsMatrixHere <<- matrix(nrow = numberOfInnerNodes, ncol = length(dataList$fragmentMasses))
  innerNodeFeaturesCountsMatrixHere <<- sparseMatrix(i = numberOfInnerNodes, j = length(dataList$fragmentMasses), x = 0)
  innerNodeFeaturesPresentHere <<- list()
  innerNodeFeaturesMeanAbundanceHere <<- list()
  innerNodeFeaturesCounterIntersectionHere <<- list()
  innerNodeFeaturesCounterUnionHere <<- list()
  innerNodePositionHere <<- vector(mode = "numeric", length = numberOfInnerNodes)
  leafHeightsHere <<- vector(mode = "numeric", length = numberOfPrecursorsFiltered)
  innerNodeMembersHere[1:numberOfPrecursorsFiltered] <<- NA
  
  analyzeTree(dataList, cluster, filter, rootIndex)
  
  ## box
  resultObj <- list()
  resultObj$innerNodeHeightIncreases <- innerNodeHeightIncreasesHere
  resultObj$innerNodeMembersTree <- innerNodeMembersTreeHere
  resultObj$innerNodeMembers <- innerNodeMembersHere
  resultObj$innerNodeMembersInner <- innerNodeMembersInnerHere
  resultObj$innerNodeFeaturesIntersection <- innerNodeFeaturesIntersectionHere
  resultObj$innerNodeFeaturesAnnotations <- innerNodeFeaturesAnnotations
  resultObj$innerNodeFeaturesUnion <- innerNodeFeaturesUnionHere
  resultObj$innerNodeFeaturesCountsMatrix <- innerNodeFeaturesCountsMatrixHere
  resultObj$innerNodeFeaturesPresent <- innerNodeFeaturesPresentHere
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
    resultObj$membersInner <- NULL
    resultObj$featuresBinaryIntersection <- featuresBinaryIntersection
    resultObj$featuresBinaryUnion <- featuresBinaryUnion
    #resultObj$features <- which(featuresBinary)
    #resultObj$featuresCounts <- rep(x = 1, times = resultObj$featuresCounter)
    resultObj$featuresCounts <- featureCounts
    #resultObj$featuresMeanAbundance <- featureValues
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
    membersInner <- c(resultObj.l$membersInner, resultObj.r$membersInner, nodeIdx)
    featuresBinaryIntersection <- resultObj.l$featuresBinaryIntersection & resultObj.r$featuresBinaryIntersection
    featuresBinaryUnion        <- resultObj.l$featuresBinaryUnion        | resultObj.r$featuresBinaryUnion
    #features <- which(featuresBinary)
    featuresCounts <- resultObj.l$featuresCounts + resultObj.r$featuresCounts
    #featuresMeanAbundance <- (resultObj.l$featuresMeanAbundance * resultObj.l$featuresCounts + resultObj.r$featuresMeanAbundance * resultObj.r$featuresCounts) / featuresCounts
    featuresCounterIntersection <- sum(featuresBinaryIntersection)
    featuresCounterUnion        <- sum(featuresBinaryUnion)
    position <- (resultObj.l$position + resultObj.r$position) / 2
    featuresAnnotations <- intersect(x = unlist(resultObj.l$featuresAnnotations), y = unlist(resultObj.r$featuresAnnotations))
    
    ## box
    resultObj <- list()
    resultObj$membersTree <- membersTree
    resultObj$members <- members
    resultObj$membersInner <- membersInner
    resultObj$featuresBinaryIntersection <- featuresBinaryIntersection
    resultObj$featuresBinaryUnion <- featuresBinaryUnion
    #resultObj$features <- union(x = resultObj.l$features, y = resultObj.r$features)
    #resultObj$features <- features
    resultObj$featuresCounts <- featuresCounts
    #resultObj$featuresMeanAbundance <- featuresMeanAbundance
    resultObj$featuresCounterIntersection <- featuresCounterIntersection
    resultObj$featuresCounterUnion <- featuresCounterUnion
    resultObj$position <- position
    resultObj$featuresAnnotations <- featuresAnnotations
    
    ## set values
    innerNodeMembersTreeHere[[nodeIdx]] <<- resultObj$membersTree
    innerNodeMembersInnerHere[[nodeIdx]] <<- resultObj$membersInner
    innerNodeMembersHere[[nodeIdx]] <<- resultObj$members
    innerNodePositionHere[[nodeIdx]] <<- resultObj$position
    innerNodeFeaturesIntersectionHere[[nodeIdx]] <<- which(resultObj$featuresBinaryIntersection)
    if(!is.null(featuresAnnotations))
      innerNodeFeaturesAnnotations[[nodeIdx]] <<- resultObj$featuresAnnotations
    innerNodeFeaturesUnionHere[[nodeIdx]] <<- which(resultObj$featuresBinaryUnion)
    #innerNodeFeaturesBinaryHere[[nodeIdx]] <<- resultObj$featuresBinary
    innerNodeFeaturesCountsMatrixHere[nodeIdx, ] <<- resultObj$featuresCounts
    innerNodeFeaturesPresentHere[[nodeIdx]] <<- sum(resultObj$featuresCounts / length(resultObj$members) >= minimumProportionOfLeafs)
    #innerNodeFeaturesMeanAbundanceHere[[nodeIdx]] <<- resultObj$featuresMeanAbundance
    innerNodeFeaturesCounterIntersectionHere[[nodeIdx]] <<- resultObj$featuresCounterIntersection
    innerNodeFeaturesCounterUnionHere[[nodeIdx]] <<- resultObj$featuresCounterUnion
    
    ## leaf heights in case of leafs and check wether to show a poi or not
    height <- cluster$height[[nodeIdx]]
    idxLeft <- cluster$merge[[nodeIdx, 1]]
    idxRight <- cluster$merge[[nodeIdx, 2]]
    increaseBelow <- FALSE
    
    ## left and right node height
    if(idxLeft < 0){
      height.l <- height
      leafHeightsHere[[-idxLeft]] <<- height.l
    } else {
      height.l <- cluster$height[[idxLeft]]
      increaseBelow <- increaseBelow | innerNodeHeightIncreasesHere[[idxLeft]]
    }
    if(idxRight < 0){
      height.r <- height
      leafHeightsHere[[-idxRight]] <<- height.r
    } else {
      height.r <- cluster$height[[idxRight]]
      increaseBelow <- increaseBelow | innerNodeHeightIncreasesHere[[idxRight]]
    }
    
    ## is only two leafs below?
    increasesHere <- any(height > height.l, height > height.r, increaseBelow)
    innerNodeHeightIncreasesHere[[nodeIdx]] <<- any(height > height.l, height > height.r, increaseBelow)
    
    if(increasesHere){
      if(idxLeft >= 0)  innerNodeHeightIncreasesHere[[idxLeft]] <<- TRUE
      if(idxRight >= 0) innerNodeHeightIncreasesHere[[idxRight]] <<- TRUE
    }
    
    #print(paste(innerNodeHeightIncreasesHere[[nodeIdx]], height, height > height.l, height.l, height > height.r, height.r, idxLeft < 0 & idxRight < 0))
    
    return(resultObj)
  }
}
analyzeTreeFromRootForAnnotations <- function(dataList, cluster, filter){
  numberOfPrecursorsFiltered <- length(filter)
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  rootIndex <- length(cluster$height)
  
  ## create fields and compute stuff
  innerNodeFeaturesAnnotations <<- list()
  
  analyzeTreeForAnnotations(dataList, cluster, filter, rootIndex)
  
  ## box
  resultObj <- list()
  resultObj$innerNodeFeaturesAnnotations <- innerNodeFeaturesAnnotations
  
  return(resultObj)
}
analyzeTreeForAnnotations <- function(dataList, cluster, filter, nodeIdx){
  if(nodeIdx < 0){
    ###################################
    ## leaf
    leafIdx <- -nodeIdx
    precursorIndex <- filter[leafIdx]
    featuresAnnotations <- dataList$annoArrayOfLists[[precursorIndex]]
    
    #print(paste(nodeIdx, featuresAnnotations, class(featuresAnnotations), length(featuresAnnotations)))
    
    ## box
    resultObj <- list()
    resultObj$featuresAnnotations <- featuresAnnotations
    
    return(resultObj)
  } else {
    ###################################
    ## inner node
    resultObj.l <- analyzeTreeForAnnotations(dataList, cluster, filter, cluster$merge[nodeIdx, 1])
    resultObj.r <- analyzeTreeForAnnotations(dataList, cluster, filter, cluster$merge[nodeIdx, 2])
    
    if(length(resultObj.l$featuresAnnotations) > 0 & length(resultObj.r$featuresAnnotations) > 0)
      featuresAnnotations <- intersect(x = unlist(resultObj.l$featuresAnnotations), y = unlist(resultObj.r$featuresAnnotations))
    else
      featuresAnnotations <- ""
    
    #print(paste(nodeIdx, cluster$merge[nodeIdx, 1], cluster$merge[nodeIdx, 2]))
    #print(paste(featuresAnnotations, length(resultObj.l$featuresAnnotations), length(resultObj.r$featuresAnnotations)))
    #print(paste(featuresAnnotations, resultObj.l$featuresAnnotations, resultObj.r$featuresAnnotations))
    #print(paste(featuresAnnotations, resultObj.l$featuresAnnotations[[1]], resultObj.r$featuresAnnotations[[1]]))
    #print(paste(class(featuresAnnotations), class(resultObj.l$featuresAnnotations), class(resultObj.r$featuresAnnotations)))
    
    resultObj <- list()
    if(any(length(featuresAnnotations) == 0, featuresAnnotations == "")){
      featuresAnnotations <- NULL
      resultObj$featuresAnnotations <- ""
    } else {
      resultObj$featuresAnnotations <- featuresAnnotations
      innerNodeFeaturesAnnotations[[nodeIdx]] <<- resultObj$featuresAnnotations
    }
    
    return(resultObj)
  }
}
getSetOfSubTreesFromRootForMass <- function(dataList, fragmentMass, filter, clusterDataList){
  yesNoFunction <- function(precursorIndex){
    features <- dataList$featureIndeces[[precursorIndex]]
    fragmentMasses <- dataList$fragmentMasses[features]
    comprisesFragmentMass <- any(fragmentMasses == fragmentMass)
    return(comprisesFragmentMass)
  }
  getSetOfSubTreesFromRoot(
    filter = filter, 
    clusterDataList = clusterDataList, 
    yesNoFunction = yesNoFunction
  )
}
getSetOfSubTreesFromRootForPrecursorSet <- function(dataList, precursorSet, filter, clusterDataList){
  yesNoFunction <- function(precursorIndex){
    take <- precursorIndex %in% precursorSet
    return(take)
  }
  getSetOfSubTreesFromRoot(
    filter = filter, 
    clusterDataList = clusterDataList, 
    yesNoFunction = yesNoFunction
  )
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
    precursorIndex <- filter[[leafIndex]]
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
colorSubTreeForAnnotations <- function(cluster, index, innerNodeAnnotations, setOfColorSets, parentIndex, parentAnnotation, parentColor, lwd = 1, lty = 1){
  #########################################################################################
  ## leaf case
  if(index<0){ # it is a leaf
    a2r_counter <<- a2r_counter + 1
    return(list(
      x = a2r_counter
    ))       
  }
  
  #########################################################################################
  ## determine color by annotations
  
  ## parent annotations
  if(any(is.null(parentIndex), length(innerNodeAnnotations) < parentIndex))
    ## no parent or no annotation
    parentAnnotations <- NULL
  else
    ## there is annotation
    parentAnnotations <- innerNodeAnnotations[[parentIndex]]
  
  ## current annotations
  if(length(innerNodeAnnotations) < index)
    currentAnnotations <- NULL
  else
    currentAnnotations <- innerNodeAnnotations[[index]]
  
  ## current color
  newAnnotations <- setdiff(x = currentAnnotations, y = parentAnnotations)
  
  if(length(newAnnotations) == 0){
    if(length(parentAnnotations) == 0){
      ## no annotations
      annotation <- "Unknown"
      color      <- "black"
    } else {
      ## there are annotations, but no new ones
      #color <- setOfColorSets[[parentIndex]][[1]]
      annotation <- parentAnnotation
      color      <- parentColor
    }
  } else{
    ## there are new annotations -> take the first new annotation
    #color <- setOfColorSets[[index]][[1]]
    #color <- setOfColorSets[[index]][[length(setOfColorSets[[index]])]]
    newAnnotationsIndeces <- match(x = newAnnotations, table = currentAnnotations)
    
    newAnnotations <- innerNodeAnnotations[[index]][newAnnotationsIndeces]
    newColors      <- setOfColorSets[[index]][newAnnotationsIndeces]
    
    annotation <- newAnnotations[[1]]
    color      <- newColors[[1]]
  }
  
  #print("hihi")
  ##print(paste(setOfColorSets[[index]]))
  #print(paste(is.null(newAnnotations), "; ", is.null(parentAnnotations), "; ", is.null(currentAnnotations), sep = ""))
  #print(paste(length(newAnnotations), "; ", length(parentAnnotations), "; ", length(currentAnnotations), sep = ""))
  #print(paste(class(newAnnotations), "; ", class(parentAnnotations), "; ", class(currentAnnotations), sep = ""))
  #print(paste(newAnnotations, "; ", parentAnnotations, "; ", currentAnnotations), sep = "")
  #print(paste(unlist(newAnnotations), "; ", unlist(parentAnnotations), "; ", unlist(currentAnnotations), sep = ""))
  #print(paste(length(unlist(newAnnotations)), "; ", length(unlist(parentAnnotations)), "; ", length(unlist(currentAnnotations)), sep = ""))
  #print(paste(index, color))
  innerNodeAnnotations[[index]] <<- annotation
  innerNodeColors[[index]] <<- color
  
  #########################################################################################
  ## draw recursively
  h.m   <- cluster$height[index]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- cluster$merge[index,1]
  
  h.l <- if(index.l<0) 0 else cluster$height[index.l]
  
  out.l   <- colorSubTreeForAnnotations(cluster = cluster, index = index.l, innerNodeAnnotations = innerNodeAnnotations, setOfColorSets = setOfColorSets, parentIndex = index, parentAnnotation = annotation, parentColor = color, lty=lty, lwd=lwd)
  x.l     <- out.l$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- cluster$merge[index,2]
  h.r <- if(index.r<0) 0 else cluster$height[index.r]
  out.r   <- colorSubTreeForAnnotations(cluster = cluster, index = index.r, innerNodeAnnotations = innerNodeAnnotations, setOfColorSets = setOfColorSets, parentIndex = index, parentAnnotation = annotation, parentColor = color, lty=lty, lwd=lwd)
  x.r     <- out.r$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  x.m  <- (x.r + x.l) / 2  
  
  #if(color != "black")
  #  print(paste("i", index, "pi", parentIndex, "il", index.l, "ir", index.r, "a", paste(currentAnnotations, collapse = "-"), "na", paste(newAnnotations, collapse = "-"), "pa", paste(parentAnnotations, collapse = "-"), "c", color))
  segments(
    x0  = c(x.l, x.l, x.r),
    x1  = c(x.l, x.r, x.r),
    y0  = c(h.l, h.m, h.r),
    y1  = c(h.m, h.m, h.m),
    col = color,
    lty = lty,
    lwd = lwd
  )
  
  list(x=x.m)
}
#########################################################################################
## annotate and process matrix
readClusterDataFromProjectFile <- function(file, progress = FALSE){
  if(progress)  setProgress(value = 0, detail = "Parsing") else print("Parsing")
  extension <- file_ext(file)
  if(extension == "gz")
    file <- gzfile(file, "r")
  
  dataFrame <- read.csv(file = file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  if(progress)  incProgress(amount = 0.1, detail = "Preprocessing") else print("Preprocessing")
  
  dataList <- readProjectData(dataFrame, progress)
  return(dataList)
}
readProjectData <- function(dataFrame, progress = FALSE){
  #########################################################################################
  ## read and parse
  columnNames <- unlist(as.matrix(dataFrame[3, ]))
  
  ## insert annotation column if not there XXX remove?
  annotationColorsName <- "AnnotationColors"
  annotationColorsMapInitValue <- paste(annotationColorsName, "={}", sep = "")
  annotationColumnName <- "Annotation"
  if(!any(columnNames == annotationColumnName, na.rm = TRUE)){
    ## insert if not there
    target <- 2
    
    if(target == 0 | target == ncol(dataFrame))
      stop("Cannot insert column!")
    
    dataFrame <- cbind(dataFrame[,1:target,drop=F], as.data.frame(x = rep(x = "", times = nrow(dataFrame)), stringsAsFactors = FALSE), dataFrame[,(target+1):ncol(dataFrame),drop=F])
    dataFrame[2, target + 1] <- annotationColorsMapInitValue
    dataFrame[3, target + 1] <- annotationColumnName
    
    columnNames <- c(columnNames[1:target], annotationColumnName, columnNames[(target+1):length(columnNames)])
  }
  colnames(dataFrame) <- columnNames
  annotationColumnIndex <- which(columnNames == annotationColumnName)
  
  ## find columns
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
  #idColumns <- which(tagsSector == "ID")
  #idColumns <- idColumns[]
  
  ## divide
  dataFrameHeader <- dataFrame[1:3, ]
  dataFrameMS1Header <- dataFrame[1:3, 1:endOfAnnotation]
  dataFrameInfos <- dataFrame[4:nrow(dataFrame), 1:endOfAnnotation]
  dataFrame <- dataFrame[4:nrow(dataFrame), ]
  
  ## mz/rt is aligned by .
  #precursorLabels <- paste(dataFrame[, 1], dataFrame[, 2], sep = " / ")
  #precursorLabels <- apply(X = dataFrame[, idColumns], MARGIN = 1, FUN = paste, collapse = " / ")
  numberOfPrecursors <- nrow(dataFrame)
  mzs <- dataFrame[, "m/z"]
  rts <- dataFrame[, "RT"]
  for(i in 1:numberOfPrecursors)
    if(length(grep(x = mzs[[i]], pattern = ".*\\..*")) == 0)
      mzs[[i]] <- paste(mzs[[i]], ".0", sep = "")
  for(i in 1:numberOfPrecursors)
    if(length(grep(x = rts[[i]], pattern = ".*\\..*")) == 0)
      rts[[i]] <- paste(rts[[i]], ".0", sep = "")
  
  tmpMatrixMz <- nchar(matrix(unlist(strsplit(x = mzs, split = "\\.")), nrow = 2))
  maxMz1 <- max(tmpMatrixMz[1, ])
  maxMz2 <- max(tmpMatrixMz[2, ])
  minMz1 <- min(tmpMatrixMz[1, ])
  minMz2 <- min(tmpMatrixMz[2, ])
  tmpMatrixRt <- nchar(matrix(unlist(strsplit(x = rts, split = "\\.")), nrow = 2))
  maxRt1 <- max(tmpMatrixRt[1, ])
  maxRt2 <- max(tmpMatrixRt[2, ])
  minRt1 <- min(tmpMatrixRt[1, ])
  minRt2 <- min(tmpMatrixRt[2, ])
  
  mzLabels <- vector(mode = "character", length = numberOfPrecursors)
  for(i in 1:numberOfPrecursors)
    mzLabels[[i]] <- paste(
      paste(rep(x = "  ", times = maxMz1 - tmpMatrixMz[1, i]), collapse = ""),
      #substr(x = mzs[[i]], start = 1, stop = tmpMatrixMz[1, i] + 1 + minMz2),
      mzs[[i]],
      paste(rep(x = "0", times = maxMz2 - tmpMatrixMz[2, i]), collapse = ""),
      sep = ""
    )
  rtLabels <- vector(mode = "character", length = numberOfPrecursors)
  for(i in 1:numberOfPrecursors)
    rtLabels[[i]] <- paste(
      paste(rep(x = "  ", times = maxRt1 - tmpMatrixRt[1, i]), collapse = ""),
      #substr(x = rts[[i]], start = 1, stop = tmpMatrixRt[1, i] + 1 + minRt2),
      rts[[i]],
      paste(rep(x = "0", times = maxRt2 - tmpMatrixRt[2, i]), collapse = ""),
      sep = ""
    )
  
  precursorLabels <- paste(mzLabels, rtLabels, sep = " / ")
  
  ## remove duplicate rows
  dupplicated <- which(duplicated(precursorLabels))
  if(length(dupplicated) > 0){
    #dupplicated2 <- dupplicated + 3
    
    precursorLabels <- precursorLabels[-dupplicated]
    dataFrame         <- dataFrame[-dupplicated, ]
    #dataFrameOriginal <- dataFrameOriginal[-dupplicated2, ]
    dataFrameInfos <- dataFrameInfos[-dupplicated, ]
  }
  
  ## annotate
  rownames(dataFrame) <- precursorLabels
  #rownames(dataFrameOriginal) <- precursorLabelsWithHeader
  numberOfPrecursors <- nrow(dataFrame)
  colnames(dataFrameInfos) <- columnNames[1:endOfAnnotation]
  rownames(dataFrameInfos) <- precursorLabels
  
  headerLabels <- c("HeaderForFragmentCounts", "HeaderForGroupsAndFragmentIntensities", "Header")
  rownames(dataFrameHeader) <- headerLabels
  colnames(dataFrameHeader) <- columnNames
  
  rownames(dataFrameMS1Header) <- headerLabels
  colnames(dataFrameMS1Header) <- columnNames[1:endOfAnnotation]
  
  ## MS1 measurement data columns
  #dataColumns <- which(!is.na(suppressWarnings(as.numeric(tagsSector))))
  dataColumns <- which(tagsSector != "ID" & tagsSector != "")
  dataColumns <- dataColumns[-which(dataColumns == annotationColumnIndex)]
  groups <- unique(unlist(tagsSector[dataColumns]))
  groupsStartEnd <- list()
  for(groupIdx in 1:length(groups))
    groupsStartEnd[[groupIdx]] <- c(min(which(tagsSector == groups[[groupIdx]])), max(which(tagsSector == groups[[groupIdx]])))
  groupsStartEndMatrix <- t(matrix(data = unlist(groupsStartEnd), nrow = 2))
  rownames(groupsStartEndMatrix) <- groups
  colnames(groupsStartEndMatrix) <- c("Start", "End")
  
  ####################
  ## MS1 measurement data: mean and LFC
  if(progress)  incProgress(amount = 0.1, detail = "Coloring") else print("Coloring")
  if(progress)  incProgress(amount = 0, detail = "Coloring init") else print("Coloring init")
  
  dataFrameMeasurements <- data.frame(matrix(nrow = nrow(dataFrame), ncol = 0))
  rownames(dataFrameMeasurements) <- rownames(dataFrame)
  
  ## column name functions
  if(progress)  incProgress(amount = 0, detail = "Coloring naming functions") else print("Coloring naming functions")
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
    dataFrameMeasurements[, dataColumnNames] <- data.matrix(dataFrame[, groupsStartEndMatrix[[groupIdx, 1]]:groupsStartEndMatrix[[groupIdx, 2]], drop=FALSE])
  }
  dataMeanColumnNameFunctionFromIndex  <- function(groupIdx){
    return(dataMeanColumnNameFunctionFromName(groups[[groupIdx]]))
  }
  dataMeanColumnNameFunctionFromName  <- function(group){
    return(paste(group, "_mean", sep = ""))
  }
  
  lfcColumnNameFunctionFromIndex <- function(groupIdxOne, groupIdxTwo){
    lfcColumnNameFunctionFromName(groups[[groupIdxOne]], groups[[groupIdxTwo]])
  }
  lfcColumnNameFunctionFromName <- function(groupOne, groupTwo){
    return(paste("LFC", groupOne, "vs", groupTwo, sep = "_"))
  }
  
  groupNameFromGroupIndex <- function(groupIdx){
    return(groups[[groupIdx]])
  }
  groupIdxFromGroupName <- function(group){
    return(match(x = group, table = groups))
  }
  
  if(progress)  incProgress(amount = 0, detail = "Coloring gather data") else print("Coloring gather data")
  ## mean data columns
  dataMeanColumnNames <- list()
  for(groupIdx in 1:length(groups)){
    dataColumnName <- dataMeanColumnNameFunctionFromIndex(groupIdx)
    dataMeanColumnNames[[groupIdx]] <- dataColumnName
    if(class(unlist(dataFrame[, groupsStartEndMatrix[[groupIdx, 1]]:groupsStartEndMatrix[[groupIdx, 2]]])) == "character")
      for(colIdx in groupsStartEndMatrix[[groupIdx, 1]]:groupsStartEndMatrix[[groupIdx, 2]])
        dataFrame[, colIdx] <- as.numeric(dataFrame[, colIdx])
    
    dataFrameMeasurements[, dataColumnName] <- apply(X = data.matrix(dataFrame[, groupsStartEndMatrix[[groupIdx, 1]]:groupsStartEndMatrix[[groupIdx, 2]]]), MARGIN = 1, FUN = mean)
    dataFrameMeasurements[is.na(dataFrameMeasurements[, dataColumnName]), dataColumnName] <- 0
    #if(any(dataFrameMeasurements[, dataColumnName] < 0))
    #  print(paste(groupIdx, dataColumnName, paste(which(dataFrameMeasurements[, dataColumnName]), collapse = "; ")))
  }
  dataMeanColumnNames <- unlist(dataMeanColumnNames)
  
  ## all replicates mean
  dataFrameMeasurements[, "meanAll"] <- apply(
    X = data.matrix(dataFrame[, 
                              unlist(apply(X = groupsStartEndMatrix, MARGIN = 1, FUN = function(x) {seq(from = x[[1]], to = x[[2]])})),
                              drop=FALSE]), 
    MARGIN = 1, FUN = mean
  )
  
  meanAllMax <- max(dataFrameMeasurements[, "meanAll"])
  dataFrameMeasurements[, "meanAll"] <- dataFrameMeasurements[, "meanAll"] / meanAllMax
  
  #print(dataFrameMeasurements[, "logMeanAll"])
  
  ## log fold change between groups
  lfcColumnNames <- list()
  for(groupIdx1 in 1:length(groups))
    for(groupIdx2 in 1:length(groups)){
      lfcColumnName <- lfcColumnNameFunctionFromIndex(groupIdx1, groupIdx2)
      lfcColumnNames[[length(lfcColumnNames) + 1]] <- lfcColumnName
      dataFrameMeasurements[, lfcColumnName] <- log(
        x = dataFrameMeasurements[, dataMeanColumnNameFunctionFromIndex(groupIdx1)] / dataFrameMeasurements[, dataMeanColumnNameFunctionFromIndex(groupIdx2)], 
        base = 2
      )
      
      ## tackle zero values
      dataFrameMeasurements[is.na(dataFrameMeasurements[, lfcColumnName]), lfcColumnName] <- 0
      dataFrameMeasurements[is.infinite(dataFrameMeasurements[, lfcColumnName]), lfcColumnName] <- 0
    }
  lfcColumnNames <- unlist(lfcColumnNames)
  
  ####################
  ## MS1 measurement data to colors
  if(progress)  incProgress(amount = 0, detail = "Coloring matrix") else print("Coloring matrix")
  matrixDataFrame <- data.matrix(dataFrameMeasurements)
  matrixDataFrame[, dataMeanColumnNames] <- log10(matrixDataFrame[, dataMeanColumnNames])
  matrixDataFrame[is.infinite(matrixDataFrame)] <- 0
  #matrixDataFrame[matrixDataFrame < 0] <- 0
  
  ## min / max
  logAbsMax <- max(matrixDataFrame[, dataMeanColumnNames])
  logFoldChangeMinMax <- c(min(matrixDataFrame[, lfcColumnNames]), max(matrixDataFrame[, lfcColumnNames]))
  logFoldChangeMax <- max(abs(logFoldChangeMinMax))
  
  ## maps
  colorMapAbsoluteData  <- makecmap(
    x = c(0, logAbsMax), n = 100, 
    #colFn = colorRampPalette(c('white', 'black'))
    colFn = colorRampPalette(rainbow(18)[10:1])
  )
  colorMapLogFoldChange <- makecmap(
    x = c(-logFoldChangeMax, logFoldChangeMax), n = 100, 
    colFn = colorRampPalette(c('blue', 'white', 'red'))
  )
  
  columnGroupLabels <- sapply(X = groups, FUN = function(x){ rep(x = x, times = length(dataColumnsNameFunctionFromName(x))) })
  columnGroupOrgLabels <- columnNames[min(groupsStartEndMatrix):max(groupsStartEndMatrix)]
  
  ## translate and box colors
  if(progress)  incProgress(amount = 0, detail = "Coloring box") else print("Coloring box")
  colorDataFrame <- dataFrameMeasurements
  colorDataFrame[, dataMeanColumnNames] <- cmap(x = matrixDataFrame[, dataMeanColumnNames], map = colorMapAbsoluteData)
  colorDataFrame[, lfcColumnNames]      <- cmap(x = matrixDataFrame[, lfcColumnNames     ], map = colorMapLogFoldChange)
  colorMatrixDataFrame <- as.matrix(colorDataFrame)
  
  ##########################################
  ## compute features
  if(progress)  incProgress(amount = 0.1, detail = "Features") else print("Features")
  
  ## get features
  dataFrameData <- as.matrix(dataFrame[, fragmentColumns[[1]] : fragmentColumns[[2]]])
  
  ## featureMatrix
  ## featureIndeces
  featureIndeces <- list()
  #featureMatrixBinary <- matrix(nrow = numberOfPrecursors, ncol = length(fragmentMasses))
  featureCount <- vector(mode = "numeric", length = numberOfPrecursors)
  matrixRows <- vector(mode = "numeric")
  matrixCols <- vector(mode = "numeric")
  matrixVals <- vector(mode = "numeric")
  itemIndex <- 1
  fragmentMassPresent <- rep(x = FALSE, times = length(fragmentMasses))
  for(i in 1:numberOfPrecursors){
    if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
      if(progress)  incProgress(amount = 0.3 / 10, detail = paste("Features ", i, " / ", numberOfPrecursors, sep = ""))
    ## data
    #featureVector <- as.numeric(unlist(dataFrame[i, fragmentColumns[[1]] : fragmentColumns[[2]]]))
    featureVector <- as.numeric(unlist(dataFrameData[i, ]))
    featureVector[is.na(featureVector)] <- 0
    featureVectorBinary <- featureVector != 0
    featureIndecesHere <- which(featureVectorBinary)
    
    ## box
    numberOfFeatures <- length(featureIndecesHere)
    rows <- rep(x = i, times = numberOfFeatures)
    cols <- featureIndecesHere
    vals <- featureVector[featureIndecesHere]
    
    matrixRows[itemIndex:(itemIndex + numberOfFeatures - 1)] <- rows
    matrixCols[itemIndex:(itemIndex + numberOfFeatures - 1)] <- cols
    matrixVals[itemIndex:(itemIndex + numberOfFeatures - 1)] <- vals
    itemIndex <- itemIndex + numberOfFeatures
    
    featureIndeces[[i]] <- featureIndecesHere
    #featureMatrixBinary[i, ] <- featureVectorBinary
    featureCount[[i]] <- numberOfFeatures
    fragmentMassPresent[cols] <- TRUE
  }
  
  if(progress)  incProgress(amount = 0.1, detail = "Feature postprocessing") else print("Feature postprocessing")
  ## featureMatrix and annotation
  featureMatrix <- sparseMatrix(i = matrixRows, j = matrixCols, x = matrixVals)
  
  fragmentMasses <- fragmentMasses[1:ncol(featureMatrix)]
  #fragmentMasses <- fragmentMasses[fragmentMassPresent]
  rownames(featureMatrix) <- precursorLabels
  colnames(featureMatrix) <- fragmentMasses
  
  ## featureIndexMatrix
  featureIndexMatrix <- matrix(nrow = numberOfPrecursors, ncol = max(sapply(X = featureIndeces, FUN = length)))
  rownames(featureIndexMatrix) <- precursorLabels
  for(i in 1:numberOfPrecursors)
    featureIndexMatrix[i, 1:length(featureIndeces[[i]])] <- featureIndeces[[i]]
  
  ## remove columns without data
  fragmentThere <- apply(X = featureMatrix, MARGIN = 2, FUN = function(x){any(x != 0)})
  minimumMass <- min(fragmentMasses[fragmentThere])
  maximumMass <- max(fragmentMasses[fragmentThere])
  
  #########################################################################################
  ## precursor annotation fields
  if(progress)  incProgress(amount = 0.1, detail = "Feature annotations") else print("Feature annotations")
  annotationValueIgnore <- "Ignore"
  annotationColorIgnore <- "red"
  
  ## present annotations
  annotations    <- vector(mode='list', length=numberOfPrecursors)
  #annotations[1:numberOfPrecursors] <- dataFrame[, annotationColumnName]
  annoVals <- dataFrame[, annotationColumnName]
  for(i in 1:numberOfPrecursors){
    #print(paste(i, annoVals[[i]], nchar(annoVals[[i]]), class(annoVals[[i]])))
    if(nchar(annoVals[[i]]) > 0){
      annotations[[i]] <- as.list(unlist(strsplit(x = annoVals[[i]], split = ", ")))
      #print(paste("a1", i, annotations[[i]], length(annotations[[i]]), class(annotations[[i]])))
    }
    else{
      annotations[[i]] <- list()
      #print(paste("a2", i, annotations[[i]], length(annotations[[i]]), class(annotations[[i]])))
    }
  }
  
  annoArrayOfLists    <- vector(mode='list', length=numberOfPrecursors)
  annoArrayIsArtifact <- vector(mode='logical', length=numberOfPrecursors)
  for(i in 1:numberOfPrecursors){
    ignoreCheck <- annotations[[i]] == annotationValueIgnore
    ignoreThere <- any(ignoreCheck)
    
    if(ignoreThere){
      idx <- which(ignoreCheck)
      annotations[[i]] <- annotations[[i]][-idx]
    }
    annoArrayOfLists[[i]]    <- annotations[[i]]
    #print(paste("b", i, annoArrayOfLists[[i]], length(annoArrayOfLists[[i]]), class(annoArrayOfLists[[i]])))
    
    annoArrayIsArtifact[[i]] <- ignoreThere
  }
  
  ## present annos with colors
  annotationColorsValue <- tagsSector[annotationColumnIndex]
  annotationColorsMapValue <- substr(
    x = annotationColorsValue, 
    start = nchar(paste(annotationColorsName, "={", sep = "")) + 1, 
    stop = nchar(annotationColorsValue) - 1
  )
  
  if(nchar(annotationColorsMapValue) > 0){
    annotationColorsMapValuePairs <- unlist(strsplit(x = annotationColorsMapValue, split = ", "))
    annotationColorsMapValues <- unlist(strsplit(x = annotationColorsMapValuePairs, split = "="))
    annotationColorsMapKeys   <- annotationColorsMapValues[seq(from = 1, to = length(annotationColorsMapValues), by = 2)]
    annotationColorsMapValues <- annotationColorsMapValues[seq(from = 2, to = length(annotationColorsMapValues), by = 2)]
  } else {
    annotationColorsMapKeys <- NULL
    annotationColorsMapValues <- NULL
  }
  
  annoPresentAnnotationsList <- list()
  annoPresentColorsList <- list()
  annoPresentAnnotationsList[[1]] <- annotationValueIgnore
  annoPresentColorsList[[1]] <- annotationColorIgnore
  if(length(annotationColorsMapKeys) > 0)
    for(i in 1:length(annotationColorsMapKeys)){
      annoPresentAnnotationsList[[1 + i]] <- annotationColorsMapKeys[[i]]
      annoPresentColorsList     [[1 + i]] <- annotationColorsMapValues[[i]]
    }
  
  ##########################################
  ## box
  if(progress)  incProgress(amount = 0.1, detail = "Boxing") else print("Boxing")
  dataList <- list()
  ## data
  #dataList$dataFrameOriginal <- dataFrameOriginal
  dataList$dataFrameHeader <- dataFrameHeader
  dataList$dataFrameMS1Header <- dataFrameMS1Header
  dataList$dataFrameInfos <- dataFrameInfos
  #dataList$dataFrame <- dataFrame
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
  dataList$meanAllMax <- meanAllMax
  dataList$logFoldChangeMax <- logFoldChangeMax
  dataList$logAbsMax <- logAbsMax
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
  dataList$featureMatrix <- featureMatrix
  #dataList$featureMatrixBinary <- featureMatrixBinary
  dataList$featureIndeces <- featureIndeces
  dataList$featureCount <- featureCount
  dataList$featureIndexMatrix <- featureIndexMatrix
  ## annotations
  dataList$annotationColumnName <- annotationColumnName
  dataList$annotationColorsName <- annotationColorsName
  dataList$annotationColumnIndex <- annotationColumnIndex
  dataList$annotationValueIgnore <- annotationValueIgnore
  dataList$annotationColorIgnore <- annotationColorIgnore
  dataList$annoArrayOfLists <- annoArrayOfLists
  dataList$annoArrayIsArtifact <- annoArrayIsArtifact
  dataList$annoPresentAnnotationsList <- annoPresentAnnotationsList
  dataList$annoPresentColorsList <- annoPresentColorsList
  
  resultObj <- getMS2plotData(dataList = dataList)
  
  maxNumberOfFragments <- max(resultObj$numberOfFragments)
  colorMapFragmentData  <- makecmap(
    x = c(0, maxNumberOfFragments), n = 100, 
    colFn = colorRampPalette(c('grey', 'black'))
  )
  #fragmentCountColors <- ?cmap(0:maxNumberOfFragments, colorMapFragmentData)
  
  dataList$numberOfFragments <- resultObj$numberOfFragments
  dataList$averageAbundance  <- resultObj$averageAbundance
  dataList$masses            <- resultObj$masses
  dataList$colorMapFragmentData <- colorMapFragmentData
  #dataList$fragmentCountColors  <- resultObj$fragmentCountColors
  
  colorMapPropabilities  <- makecmap(
    x = c(0, 1), n = 100, 
    colFn = colorRampPalette(c('black', 'grey'))
  )
  dataList$colorMapPropabilities <- colorMapPropabilities
  
  if(progress)  setProgress(1) else print("Ready")
  
  ## 950 932 688
  ## 634 336 248
  ## 321 972 296
  ##   9 090 088
  ##  11 753 432
  #print(sort( sapply(ls(),function(x){object.size(get(x))})))
  #memory.profile()
  
  return(dataList)
}
filterData <- function(dataList, groups, filter_average, filter_lfc, filterList_ms2_masses, filter_ms2_ppm, filter_ms1_mass, filter_ms1_ppm, includeIgnoredPrecursors, progress = FALSE){
  ##########################################
  ## filter
  filter <- rep(x = TRUE, times = dataList$numberOfPrecursors)
  
  ## filter_average
  if(!is.null(filter_average))
    filter <- filter & apply(X = as.data.frame(dataList$dataFrameMeasurements[, sapply(X = as.vector(groups), FUN = dataList$dataMeanColumnNameFunctionFromName)]), MARGIN = 1, FUN = mean) >= filter_average
  
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
  
  ## filter_ms1_mass, filter_ms1_ppm
  if(!is.null(filter_ms1_mass) & !is.null(filter_ms1_ppm) & length(filter_ms1_mass) > 0){
    precursorMasses <- as.numeric(dataList$dataFrameInfos$"m/z")
    error <- precursorMasses * filter_ms1_ppm / 1E6
    distances <- abs(precursorMasses - filter_ms1_mass)
    filterMS1mass <- distances <= error
    
    filter <- filter & filterMS1mass
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
calculateDistanceMatrix <- function(dataList, filter, distanceMeasure = "Jaccard", progress = FALSE){
  numberOfPrecursors <- length(filter)
  
  ## compute distance matrix:
  distanceMatrix <- NULL
  switch(distanceMeasure,
         "Jaccard"={
           featureIndeces <- dataList$featureIndeces[filter]
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Jaccard (intensity-weighted pure)"={
           featureIndeces <- dataList$featureIndeces[filter]
           featureMatrix <- dataList$featureMatrix[filter, ]
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Jaccard (intensity-weighted map) exp"={
           ## intersection i and j: featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
           ## diff i - j: featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]]
           ## diff j - i: featureIndeces[[j]][!featureIndeces[[j]] %in% featureIndeces[[i]]]
           ## union i or j: c(featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]], featureIndeces[[j]])
           featureIndeces <- dataList$featureIndeces[filter]
           featureMatrix <- dataList$featureMatrix[filter, ]
           featureMatrix[featureMatrix <  0.2 & featureMatrix >= 0.01] <- 0.01
           featureMatrix[featureMatrix <  0.4 & featureMatrix >= 0.20] <- 0.2
           featureMatrix[featureMatrix >= 0.4] <- 1
           featureMatrix <- as.matrix(featureMatrix)
           featureMatrixBinary <- dataList$featureMatrixBinary[filter, ]
           
          "%ni%" = Negate("%in%")
          
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
             else          print(paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
             
             featureMatrixBinaryHere <- t(featureMatrixBinary)
             
             ## intersecting features
             intersections <- featureMatrixBinary[i, ] & featureMatrixBinaryHere
             onlyIs <- featureMatrixBinary[i, ] & (!featureMatrixBinaryHere)
             onlyJs <- featureMatrixBinaryHere  & (!featureMatrixBinary[i, ])
             
             onlyIs[is.na(onlyIs)] <- FALSE
             onlyJs[is.na(onlyJs)] <- FALSE
             
             sumOnlyIs  <- apply(X = onlyIs, MARGIN = 2, FUN = function(x){ sum(featureMatrix[i, ] & x) })
             #sumOnlyIs <- apply(X = onlyIs, MARGIN = 2, FUN = function(x){ sum(featureMatrix[i, x]) })
             featureMatrix2 <- featureMatrix
             featureMatrix2[!t(onlyJs)] <- 0
             sumOnlyJs <- apply(X = featureMatrix2, MARGIN = 1, FUN = sum )
             
             featureMatrix2 <- featureMatrix
             featureMatrix2[!t(intersections)] <- 0
             intersectionSums <- apply(X = featureMatrix2, MARGIN = 1, FUN = sum )
             
             unionSums <- intersectionSums + sumOnlyIs + sumOnlyJs
             distances <- 1 - intersectionSums / unionSums
             
             distanceMatrix[i, ] <- distances
           }
         },
         "Jaccard (intensity-weighted)"={
           ## map the fragment intensities to three intensity categories and compute distance from intensity categories of intersecting / all fragments
           ## 
           ## intersection i and j: featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
           ## diff i - j: featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]]
           ## diff j - i: featureIndeces[[j]][!featureIndeces[[j]] %in% featureIndeces[[i]]]
           ## union i or j: c(featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]], featureIndeces[[j]])
           featureIndeces <- dataList$featureIndeces[filter]
           featureMatrix <- dataList$featureMatrix[filter, ]
           featureMatrix[featureMatrix <  0.2 & featureMatrix >= 0.01] <- 0.01
           featureMatrix[featureMatrix <  0.4 & featureMatrix >= 0.20] <- 0.2
           featureMatrix[featureMatrix >= 0.4] <- 1
           featureMatrix <- as.matrix(featureMatrix)
          # intensityMapping <- function(x){
          #   if(x == 0){
          #     ## x = 0
          #     return(0)
          #   } else if(x < 0.2){
          #     ## 0 <= x < 0.2
          #     return(0.01)
          #   } else if(x < 0.4){
          #     ## 0.2 <= x < 0.4
          #     return(0.2)
          #   } else {
          #     ## 0.4 <= x <= Inf
          #     return(1)
          #   }
          # }
          # intensityMapping <- Vectorize(FUN = intensityMapping, vectorize.args = "x")
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
               else          print(paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
             for(j in 1:numberOfPrecursors){
               if(i == j){
                 distanceMatrix[i, j] <- 0
                 next
               }
               
               ## intersecting features
               intersection <- featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
               
               if(length(intersection) == 0){
                 ## max distance
                 distance <- 1#sumOnlyI + sumOnlyJ
               } else {
                 onlyI <- featureIndeces[[i]][!featureIndeces[[i]] %in% featureIndeces[[j]]]
                 onlyJ <- featureIndeces[[j]][!featureIndeces[[j]] %in% featureIndeces[[i]]]
                 
                 if(length(onlyI) == 0){
                   sumOnlyI <- 0
                 } else {
                   sumOnlyI <- sum(featureMatrix[i, onlyI, drop = FALSE])
                 }
                 if(length(onlyJ) == 0){
                   sumOnlyJ <- 0
                 } else {
                   sumOnlyJ <- sum(featureMatrix[j, onlyJ, drop = FALSE])
                 }
                 
                 #maxIntensity <- apply(X = as.matrix(featureMatrix[c(i, j), intersection, drop = FALSE]), MARGIN = 2, FUN = max)
                 maxIntensity <- apply(X = featureMatrix[c(i, j), intersection, drop = FALSE], MARGIN = 2, FUN = max)
                 #maxIntensity <- featureMatrix[i, intersection]
                 intersectionSum <- sum(maxIntensity)
                 unionSum <- intersectionSum + sumOnlyI + sumOnlyJ
                 distance <- 1 - intersectionSum / unionSum
               }
               
               distanceMatrix[i, j] <- distance
             }
           }
         },
         "Jaccard (intensity-weighted map) saveOld"={
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
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
             else          print(paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
           featureMatrix <- as.matrix(dataList$featureMatrix[filter, ])
           
           similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
           featureMatrix <- as.matrix(dataList$featureMatrix[filter, ])
           fragmentFrequency <- dataList$fragmentFrequency
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
           featureMatrix <- as.matrix(dataList$featureMatrix[filter, ])
           fragmentFrequency <- dataList$fragmentFrequency
           
           similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", counter, " / ", numberOfPrecursors, sep = ""))
             intersectionSum <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(y){  sum(fragmentFrequency[x[x %in% y]], na.rm = TRUE)  })
             unionSum        <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(y){  sum(fragmentFrequency[c(x[!x %in% y], y)], na.rm = TRUE)  })
             
             1 - intersectionSum / unionSum
           }
           )
         },
         "Manhatten"={
           featureIndeces <- dataList$featureIndeces[filter]
           featureMatrix <- as.matrix(dataList$featureMatrix[filter, ])
           
           ## Rasmussen 2008
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "NDP (Normalized dot product)"={
           featureIndeces <- dataList$featureIndeces[filter]
           featureMatrix <- dataList$featureMatrix[filter, ]
           
           ## Gaquerel 2015: standard normalized dot product (NDP) / cosine correlation
           similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(numberOfPrecursors >= 10 & ((i %% (as.integer(numberOfPrecursors/10))) == 0))
               if(progress)  incProgress(amount = 1 / 10, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
             for(j in 1:numberOfPrecursors){
               if(i == j){
                 similarityMatrix[i, j] <- 0
                 next
               }
               
               intersection    <- intersect(x = featureIndeces[[i]], y = featureIndeces[[j]])
               intersectionSum <- sum(sqrt(featureMatrix[i, intersection]) * dataList$fragmentMasses[intersection]^2 * sqrt(featureMatrix[j, intersection]) * dataList$fragmentMasses[intersection]^2)^2
               iSum            <- sum((sqrt(featureMatrix[i, featureIndeces[[i]]]) * dataList$fragmentMasses[featureIndeces[[i]]]^2)^2)
               jSum            <- sum((sqrt(featureMatrix[j, featureIndeces[[j]]]) * dataList$fragmentMasses[featureIndeces[[j]]]^2)^2)
               
               similarityMatrix[i, j] <- intersectionSum / (iSum * jSum)
             }
           }
           distanceMatrix <- max(similarityMatrix) - similarityMatrix
         },
         stop(paste("Unknown distance (", distance, ")!", sep = ""))
  )
  
  rownames(distanceMatrix) <- dataList$precursorLabels[filter]
  colnames(distanceMatrix) <- dataList$precursorLabels[filter]
  
  returnObj <- list(
    distanceMatrix = distanceMatrix,
    filter = filter,
    distanceMeasure = distanceMeasure
  )
  
  return(returnObj)
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
  
  ## optimal leaf ordering
  #opt <- order.optimal(dist = dist, merge = cluster$merge)
  #cluster$merge <- opt$merge
  #cluster$order <- opt$order
  
  ## compute (transitive) cluster members, cluster positions, and leaf heights
  if(progress)  incProgress(amount = 0.5, detail = "Analyze cluster")
  
  resultObj <- analyzeTreeFromRoot(dataList, cluster = cluster, filter)
  
  innerNodeHeightIncreases <- unlist(resultObj$innerNodeHeightIncreases)
  innerNodeMembersInner <- resultObj$innerNodeMembersInner
  innerNodeMembersTree  <- resultObj$innerNodeMembersTree
  innerNodeMembers <- resultObj$innerNodeMembers
  #innerNodeFeaturesBinary <- resultObj$innerNodeFeaturesBinary
  innerNodeFeaturesIntersection <- resultObj$innerNodeFeaturesIntersection
  innerNodeFeaturesUnion <- resultObj$innerNodeFeaturesUnion
  innerNodeFeaturesCountsMatrix <- resultObj$innerNodeFeaturesCountsMatrix
  innerNodeFeaturesPresent <- resultObj$innerNodeFeaturesPresent
  #innerNodeFeaturesMeanAbundance <- resultObj$innerNodeFeaturesMeanAbundance
  innerNodeFeaturesIntersectionCounter <- resultObj$innerNodeFeaturesIntersectionCounter
  innerNodeFeaturesUnionCounter <- resultObj$innerNodeFeaturesUnionCounter
  innerNodePosition <- resultObj$innerNodePosition
  leafHeights <- resultObj$leafHeights
  
  #innerNodeFeaturesCountsMatrix <- matrix(data = unlist(innerNodeFeaturesCounts), nrow = length(innerNodeFeaturesCounts))
  
  ## dendrogram leaf ends for normal plot
  #leafHeights <- leafHeights - leafHeightSpacing
  leafHeights <- rep(x = 0, times = length(leafHeights))
  
  ## compute x- and y-coordinates and point-labels
  coordinatesX <- unlist(c(innerNodePosition, match(x = 1:numberOfPrecursorsFiltered, table = cluster$order)))
  coordinatesY <- unlist(c(cluster$height, leafHeights))
  
  precursorFeatureCount <- dataList$featureCount[filter]
  innerNodeUnionlabels <- as.character(innerNodeFeaturesUnionCounter)
  innerNodeUnionlabels[!innerNodeHeightIncreases] <- ""
  innerNodeIntersectionlabels <- as.character(innerNodeFeaturesIntersectionCounter)
  innerNodeIntersectionlabels[!innerNodeHeightIncreases] <- ""
  
  drawPoi <- unlist(c(
    innerNodeHeightIncreases, 
    rep(x = TRUE, times = length(filter))
  ))
  union <- unlist(c(
    innerNodeUnionlabels, 
    precursorFeatureCount
  ))
  intersection <- unlist(c(
    innerNodeIntersectionlabels, 
    precursorFeatureCount
    #vector(mode = "character", length = length(filter))
  ))
  intersectionSmooth <- c(innerNodeFeaturesPresent, precursorFeatureCount)
  poiLabels <- unlist(c(1:numberOfInnerNodes, -(1:numberOfPrecursorsFiltered)))
  
  ##########################################
  ## box
  if(progress)  incProgress(amount = 0.5, detail = "Boxing")
  filterList <- list()
  ## filter
  filterList$filter <- filter
  #filterList$distanceMatrix <- distanceMatrix
  filterList$numberOfPrecursorsFiltered <- numberOfPrecursorsFiltered
  ## cluster
  
  filterList$innerNodeMembersTree <- innerNodeMembersTree
  filterList$innerNodeMembersInner <- innerNodeMembersInner
  filterList$innerNodeMembers <- innerNodeMembers
  #filterList$innerNodeFeaturesBinary <- innerNodeFeaturesBinary
  filterList$innerNodeFeaturesIntersection <- innerNodeFeaturesIntersection
  filterList$innerNodeFeaturesUnion <- innerNodeFeaturesUnion
  filterList$innerNodeFeaturesCountsMatrix <- innerNodeFeaturesCountsMatrix
  filterList$innerNodeFeaturesPresent <- innerNodeFeaturesPresent
  #filterList$innerNodeFeaturesMeanAbundance <- innerNodeFeaturesMeanAbundance
  filterList$innerNodeFeaturesIntersectionCounter <- innerNodeFeaturesIntersectionCounter
  filterList$innerNodeFeaturesUnionCounter <- innerNodeFeaturesUnionCounter
  filterList$cluster <- cluster
  ## poi
  filterList$numberOfPois    <- length(coordinatesX)
  filterList$poiCoordinatesX <- coordinatesX
  filterList$poiCoordinatesY <- coordinatesY
  filterList$poiUnion        <- union
  filterList$poiIntersection <- intersection
  filterList$poiIntersectionSmooth <- intersectionSmooth
  filterList$poiLabels       <- poiLabels
  filterList$drawPoi         <- drawPoi
  
  if(progress)  setProgress(1)
  
  ## 1 167 983 544
  ##   127 968 176
  ##   127 922 016
  ##   128 042 056
  ##    13 761 512
  #print(sort( sapply(ls(),function(x){object.size(get(x))})))
  
  return(filterList)
}
calculatePCA <- function(dataList, filterObj, scaling, logTransform){
  dataFrame <- dataList$dataFrameMeasurements[filterObj$filter, dataList$dataColumnsNameFunctionFromNames(filterObj$groups)]
  dataFrame <- t(dataFrame)
  
  if(logTransform){
    dataFrame[dataFrame < 1] <- 1
    dataFrame <- log2(dataFrame)
  }
  
  switch(scaling,
         "None"={
           ## do nothing
           dataFrame2 <- dataFrame
         },
         "Mean center"={
           ## subtract mean
           dataFrame2 <- as.data.frame(apply(X = dataFrame, MARGIN = 2, FUN = function(x){
             x - mean(x = x)
           }))
         },
         "Autoscaling (unit variance)"={
           ## subtract mean and divide by variance
           dataFrame2 <- as.data.frame(apply(X = dataFrame, MARGIN = 2, FUN = function(x){
             (x - mean(x = x)) / sd(x = x)
           }))
         },
         "Pareto"={
           ## subtract mean and divide by sqrt of variance
           dataFrame2 <- as.data.frame(apply(X = dataFrame, MARGIN = 2, FUN = function(x){
             (x - mean(x = x)) / sqrt(sd(x = x))
           }))
         },
         stop(paste("Unknown scaling (", scaling, ")!", sep = ""))
  )
  
  dataFrame2[is.na(dataFrame2)] <- dataFrame[is.na(dataFrame2)]
  
  numberOfComponents <- 5
  returnObj <- list()
  pcaLibrary <- c("stats", "FactoMineR", "pcaMethods")[[2]]
  switch(pcaLibrary,
         "stats"={
           ## "stats" package
           pca <- prcomp(x = dataFrame2, retx = TRUE, center = FALSE, scale. = FALSE)
           returnObj$scores   <- pca$x
           returnObj$loadings <- pca$rotation
           returnObj$variance <- pca$sdev
         },
         "FactoMineR"={
           ## "FactoMineR" package
           pca = PCA(X = dataFrame2, graph = FALSE, scale.unit = FALSE, ncp = numberOfComponents)
           returnObj$scores   <- pca$ind$coord
           returnObj$loadings <- pca$var$coord
           returnObj$variance <- pca$eig$"percentage of variance"
         },
         "pcaMethods"={
           ## "pcaMethods" package
           pca <- pca(object = dataFrame2, method = "svd", nPcs = numberOfComponents, scale = "none", center = FALSE)
           returnObj$scores   <- pca@scores
           returnObj$loadings <- pca@loadings
           returnObj$variance <- pca@sDev
         },
         stop(paste("Unknown PCA library (", pcaLibrary, ")!", sep = ""))
  )
  
  ## valid data
  #$ scores  : num [1:19, 1:5] -355 807 808 -2720 748 ...
  #..- attr(*, "dimnames")=List of 2
  #.. ..$ : chr [1:19] "A1_1" "A1_2" "A1_3" "A1_4" ...
  #.. ..$ : chr [1:5] "Dim.1" "Dim.2" "Dim.3" "Dim.4" ...
  #$ loadings: num [1:358, 1:5] 1.07 1.11 12.3 -2.21 36.6 ...
  #..- attr(*, "dimnames")=List of 2
  #.. ..$ : chr [1:358] "  105.0700 /   161.85" "  122.0968 /   162.04" "  125.9873 /     48.28" "  129.0557 / 1019.34" ...
  #.. ..$ : chr [1:5] "Dim.1" "Dim.2" "Dim.3" "Dim.4" ...
  #$ variance: num [1:18] 35.96 21.67 11.66 8.15 4.56 ...
  #
  ## artificial data two groups
  #$ scores  : num [1:2, 1] 0 0
  #..- attr(*, "dimnames")=List of 2
  #.. ..$ : chr [1:2] "A_1" "B_2"
  #.. ..$ : chr "Dim.1"
  #$ loadings: Named num [1:455] 0 0 0 0 0 0 0 0 0 0 ...
  #..- attr(*, "names")= chr [1:455] "  52.1 / 18.74" "  53.0 /   9.13" "  55.1 / 19.60" "  56.1 / 12.91" ...
  #$ variance: num NaN
  #
  ## artificial data one group
  #$ scores  : num [1:455, 1] 2.89e-15 2.89e-15 2.89e-15 2.89e-15 2.89e-15 ...
  #..- attr(*, "dimnames")=List of 2
  #.. ..$ : chr [1:455] "  52.1 / 18.74" "  53.0 /   9.13" "  55.1 / 19.60" "  56.1 / 12.91" ...
  #.. ..$ : chr "Dim.1"
  #$ loadings: num 2.89e-15
  #$ variance: num 100
  #
  #str(returnObj)
  
  if(any(is.na(returnObj$variance), length(returnObj$variance) == 1)){
    ## in case of artificial data
    numberOfPrecursors <- dataList$numberOfPrecursors
    numberOfSamples    <- nrow(returnObj$scores)
    returnObj$scores   <- matrix(nrow = numberOfSamples, ncol = numberOfComponents)
    returnObj$loadings <- matrix(nrow = numberOfPrecursors, ncol = numberOfComponents)
    returnObj$variance <- vector(mode = "numeric", length = numberOfComponents)
  }
  
  return(returnObj)
}
#########################################################################################
## data fetching
getMetFragLink <- function(dataList, precursorIndex){
  features <- dataList$featureIndeces[[precursorIndex]]
  fragmentsX <- dataList$fragmentMasses[features]
  fragmentsY <- as.numeric(dataList$featureMatrix[precursorIndex, features])
  fragmentsY[fragmentsY > 1] <- 1
  
  precursorMass  <- as.numeric(dataList$dataFrameInfos$"m/z"[[precursorIndex]])
  if("Adduct ion name" %in% colnames(dataList$dataFrameInfos))
    adduct <- dataList$dataFrameInfos$"Adduct ion name"[[precursorIndex]]
  else
    adduct <- dataList$dataFrameInfos$"Adduct.ion.name"[[precursorIndex]]
  neutralMassCorrection <- NA
  ionMode <- NA
  print(adduct)
  #generateLink <- TRUE
  switch(adduct,
         "[M-H]-"={
           neutralMassCorrection <- 1.008
           ionMode <- -1
         },
         "[M+H]+"={
           neutralMassCorrection <- -1.008
           ionMode <- 1
         },
         "Unknown"={
           #generateLink <- FALSE
           neutralMassCorrection <- NA
           ionMode <- NA
         },{
           #stop(paste("Unknown adduct (", adduct, ")!", sep = ""))
           print(paste("###### Unknown adduct (", adduct, ")!", sep = ""))
           #generateLink <- FALSE
           neutralMassCorrection <- NA
           ionMode <- NA
         }
  )
  neutralMass <- precursorMass + neutralMassCorrection
  
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
  
  if(!is.na(neutralMass)){
    #writeClipboard(landingPageUrl, format = 1)
    landingPageUrlForLink <- landingPageUrl
  }
  else
    landingPageUrlForLink <- NULL
  
  return(landingPageUrlForLink)
}
getMS2spectrum <- function(dataList, clusterDataList, clusterLabel){
  if(clusterLabel < 0){
    ###############################################
    ## leaf
    precursorIndex <- clusterDataList$filter[[-clusterLabel]]
    features <- dataList$featureIndeces[[precursorIndex]]
    
    fragmentsX <- dataList$fragmentMasses[features]
    fragmentsY <- as.numeric(dataList$featureMatrix[precursorIndex, features])
    fragmentsY[fragmentsY > 1] <- 1
    fragmentsColor <- rep(x = "black", times = length(fragmentsY))
    
    infoText <- paste("The MS/MS spectrum of MS feature '", trimws(gsub(x = clusterDataList$cluster$labels[[-clusterLabel]], pattern = " +", replacement = " ")), "' comprises ", length(fragmentsX), " fragments.", sep = "")
    
    landingPageUrlForLink <- getMetFragLink(dataList, precursorIndex)
    
    precursorSet <- precursorIndex
  } else {
    ###############################################
    ## inner node
    clusterIndex <- clusterLabel
    clusterMembers <- sort(clusterDataList$innerNodeMembers[[clusterIndex]])
    
    featuresIntersection <- clusterDataList$innerNodeFeaturesIntersection[[clusterIndex]]
    featuresUnion <- clusterDataList$innerNodeFeaturesUnion[[clusterIndex]]
    #fragmentsX <- dataList$fragmentMasses[featuresIntersection]
    #fragmentsY <- apply(X = data.matrix(dataList$featureMatrix[clusterMembers, featuresIntersection]), MARGIN = 2, FUN = mean)
    fragmentsX <- dataList$fragmentMasses[featuresUnion]
    fragmentsY <- apply(X = data.matrix(dataList$featureMatrix[clusterMembers, featuresUnion]), MARGIN = 2, FUN = mean)
    #fragmentsColor <- rep(x = "black", times = length(fragmentsY))
    
    props <- clusterDataList$innerNodeFeaturesCounts[clusterIndex, featuresUnion] / length(clusterDataList$innerNodeMembers[[clusterIndex]])
    #fragmentsColor <- cmap(x = props, map = dataList$colorMapPropabilities)
    fragmentsColor <- vector(length = length(fragmentsX))
    fragmentsColor[props >= minimumProportionOfLeafs] <- "black"
    fragmentsColor[props < minimumProportionOfLeafs] <- "grey"
    
    fragmentsX <- fragmentsX[props > minimumProportionToShowFragment]
    fragmentsY <- fragmentsY[props > minimumProportionToShowFragment]
    fragmentsColor <- fragmentsColor[props > minimumProportionToShowFragment]
    
    infoText <- paste("This cluster comprises ", length(clusterMembers), " MS features which have ", length(fragmentsX), " fragment(s) in common.", sep = "")
    landingPageUrlForLink <- NULL
    
    precursorSet <- clusterMembers
  }
  
  ## selected data
  order <- order(fragmentsX)
  fragmentsX <- fragmentsX[order]
  fragmentsY <- fragmentsY[order]
  fragmentsColor <- fragmentsColor[order]
  
  numberOfPrecursors <- length(precursorSet)
  
  resultObj <- list()
  resultObj$fragmentMasses <- fragmentsX
  resultObj$fragmentAbundances <- fragmentsY
  resultObj$fragmentColor <- fragmentsColor
  resultObj$infoText <- infoText
  resultObj$landingPageUrl <- landingPageUrlForLink
  resultObj$precursorSet <- precursorSet
  resultObj$numberOfPrecursors <- numberOfPrecursors
  
  return(resultObj)
}
## XXX adapt getTableFromTreeSelection
getTableFromPrecursorSet <- function(dataList, precursorSet){
  ###############################################
  ## table data
  numberOfPrecursors <- length(precursorSet)
  
  ## measurements
  columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
  dataFrameMeasurements     <- data.frame(dataList$dataFrameMeasurements[precursorSet, columnNames, drop=FALSE])
  colnames(dataFrameMeasurements) <- columnNames
  rownames(dataFrameMeasurements) <- dataList$precursorLabels[precursorSet]
  
  ## MS2 fragments
  props <- apply(
    X = dataList$featureMatrix[precursorSet, , drop = FALSE], 
    MARGIN = 2, 
    FUN = function(x){ sum(x != 0) / numberOfPrecursors }
  )
  featureIndeces <- which(props > minimumProportionToShowFragment)
  featureMatrix <- data.frame(as.matrix(dataList$featureMatrix[precursorSet, featureIndeces, drop = FALSE]))
  rownames(featureMatrix) <- dataList$precursorLabels[precursorSet]
  colnames(featureMatrix) <- colnames(dataList$featureMatrix)[featureIndeces]
  
  ## annotations
  setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
  setOfAnnotations <- unlist(lapply(X = setOfAnnotationSets, FUN = function(x){
    paste(x, collapse = ", ")
  }))
  annotationDataFrame <- data.frame(Annotations = setOfAnnotations, row.names = rownames(dataFrameMeasurements))
  
  ## box
  precursorLabels <- rownames(dataFrameMeasurements)
  ms1abundanceDataFrame <- dataFrameMeasurements
  ms2fragmentDataFrame <- featureMatrix
  
  resultObj <- list(
    precursorSet = precursorSet,
    featureIndeces = featureIndeces,
    ms1abundanceDataFrame = ms1abundanceDataFrame,
    ms2fragmentDataFrame = ms2fragmentDataFrame,
    annotationDataFrame = annotationDataFrame
  )
  
  return(resultObj)
}
getTableFromTreeSelection2 <- function(dataList, precursorSet, featuresIntersection){
  ###############################################
  ## table data
  
  ## measurements
  columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
  dataFrameMeasurements     <- data.frame(dataList$dataFrameMeasurements[precursorSet, columnNames])
  colnames(dataFrameMeasurements) <- columnNames
  rownames(dataFrameMeasurements) <- dataList$precursorLabels[precursorSet]
  
  ## MS2 fragments
  featureMatrix <- data.frame(as.matrix(dataList$featureMatrix[precursorSet, featuresIntersection, drop = FALSE]))
  rownames(featureMatrix) <- dataList$precursorLabels[precursorSet]
  
  ## annotations
  setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
  setOfAnnotations <- unlist(lapply(X = setOfAnnotationSets, FUN = function(x){
    paste(x, collapse = ", ")
  }))
  annotationDataFrame <- data.frame(Annotations = setOfAnnotations, row.names = rownames(dataFrameMeasurements))
  
  ## box
  precursorLabels <- rownames(dataFrameMeasurements)
  ms1abundanceDataFrame <- dataFrameMeasurements
  ms2fragmentDataFrame <- featureMatrix
  
  resultObj <- list()
  resultObj$precursorSet <- precursorSet
  resultObj$precursorLabels <- precursorLabels
  resultObj$ms1abundanceDataFrame <- ms1abundanceDataFrame
  resultObj$ms2fragmentDataFrame <- ms2fragmentDataFrame
  resultObj$annotationDataFrame <- annotationDataFrame
  
  return(resultObj)
}
getTableFromTreeSelection <- function(dataList, clusterDataList, clusterLabel){
  if(clusterLabel < 0){
    ###############################################
    ## leaf
    precursorIndex <- clusterDataList$filter[[-clusterLabel]]
    precursorSet <- precursorIndex
    features <- dataList$featureIndeces[[precursorIndex]]
    
    ###############################################
    ## table data
    
    ## measurements
    columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
    dataFrameMeasurements     <- data.frame(dataList$dataFrameMeasurements[precursorIndex, columnNames])
    colnames(dataFrameMeasurements) <- columnNames
    
    ## MS2 fragments
    featureMatrix <- data.frame(as.matrix(dataList$featureMatrix[precursorIndex, features, drop = FALSE]))
    rownames(featureMatrix) <- rownames(dataList$dataFrameMeasurements)[[precursorIndex]]
    
    ## annotations
    setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
    setOfAnnotations <- unlist(lapply(X = setOfAnnotationSets, FUN = function(x){
      paste(x, collapse = ", ")
    }))
    annotationDataFrame <- data.frame(Annotations = setOfAnnotations, row.names = rownames(dataFrameMeasurements))
    
    ## box
    precursorLabels <- rownames(dataFrameMeasurements)
    ms1abundanceDataFrame <- dataFrameMeasurements
    ms2fragmentDataFrame <- featureMatrix
    featuresIntersection <- features
  } else {
    ###############################################
    ## inner node
    clusterIndex <- clusterLabel
    precursorSet <- sort(clusterDataList$innerNodeMembers[[clusterIndex]])
    featuresIntersection <- clusterDataList$innerNodeFeaturesIntersection[[clusterIndex]]
    
    ###############################################
    ## table data
    
    ## measurements
    columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
    dataFrameMeasurements   <- dataList$dataFrameMeasurements[precursorSet, columnNames]
    
    ## MS2 fragments
    maximumNumberOfFeatures <- 100
    featureMatrix <- data.frame(as.matrix(dataList$featureMatrix[precursorSet, featuresIntersection]))
    colnames(featureMatrix) <- dataList$fragmentMasses[featuresIntersection]
    if(ncol(featureMatrix) > maximumNumberOfFeatures)
      featureMatrix <- featureMatrix[, 1:maximumNumberOfFeatures]
    
    ## annotations
    setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
    setOfAnnotations <- unlist(lapply(X = setOfAnnotationSets, FUN = function(x){
      paste(x, collapse = ", ")
    }))
    annotationDataFrame <- data.frame(Annotations = setOfAnnotations, row.names = rownames(dataFrameMeasurements))
    
    ## box
    precursorLabels <- rownames(dataFrameMeasurements)
    ms1abundanceDataFrame <- dataFrameMeasurements
    ms2fragmentDataFrame <- featureMatrix
  }
  
  resultObj <- list()
  resultObj$precursorSet <- precursorSet
  resultObj$precursorLabels <- precursorLabels
  resultObj$featuresIntersection <- featuresIntersection
  resultObj$ms1abundanceDataFrame <- ms1abundanceDataFrame
  resultObj$ms2fragmentDataFrame <- ms2fragmentDataFrame
  resultObj$annotationDataFrame <- annotationDataFrame
  
  return(resultObj)
}
getPrecursorSetFromTreeSelections <- function(clusterDataList, clusterLabels){
  precursorSet <- NULL
  for(clusterLabel in clusterLabels)
    precursorSet <- c(precursorSet, getPrecursorSetFromTreeSelection(clusterDataList, clusterLabel))
  return(precursorSet)
}
getPrecursorSetFromTreeSelection <- function(clusterDataList, clusterLabel){
  if(clusterLabel < 0){
    ###############################################
    ## leaf
    precursorIndex <- clusterDataList$filter[[-clusterLabel]]
    precursorSet <- precursorIndex
  } else {
    ###############################################
    ## inner node
    clusterIndex <- clusterLabel
    precursorSet <- sort(clusterDataList$innerNodeMembers[[clusterIndex]])
  }
  return(precursorSet)
}
getMS2spectrumOfPrecursor <- function(dataList, precursorIndex){
  featureIndeces <- dataList$featureIndeces[[precursorIndex]]
  featureMasses  <- dataList$fragmentMasses[featureIndeces]
  featureValues  <- dataList$featureMatrix [precursorIndex, featureIndeces]
  featureValues[featureValues > 1] <- 1
  featureColor   <- rep(x = "black", times = length(featureMasses))
  
  infoText <- paste("The MS/MS spectrum of MS feature '", trimws(gsub(x = dataList$precursorLabels[[precursorIndex]], pattern = " +", replacement = " ")), "' comprises ", length(featureMasses), " fragments.", sep = "")
  
  resultObj <- list()
  resultObj$fragmentMasses <- featureMasses
  resultObj$fragmentAbundances <- featureValues
  resultObj$fragmentColor <- featureColor
  resultObj$infoText <- infoText
  
  return(resultObj)
}
getMS2plotData <- function(dataList){
  numberOfFragments <- apply(X = dataList$featureMatrix, MARGIN = 2, FUN = function(x){sum(x != 0)})
  sumOfAbundances <- apply(X = dataList$featureMatrix, MARGIN = 2, FUN = sum)
  averageAbundance <- sumOfAbundances / numberOfFragments
  averageAbundance[numberOfFragments == 0] <- 0
  masses <- dataList$fragmentMasses
  
  presentFragments <- numberOfFragments > 0
  
  numberOfFragments <- numberOfFragments[presentFragments]
  averageAbundance  <- averageAbundance[presentFragments]
  masses            <- masses[presentFragments]
  
  resultObj <- list()
  resultObj$numberOfFragments <- numberOfFragments
  resultObj$averageAbundance  <- averageAbundance
  resultObj$masses            <- masses
  
  return(resultObj)
}
plotFragments <- function(dataList, xInterval = NULL){
  if(is.null(xInterval)){
    xMin <- min(dataList$masses)
    xMax <- max(dataList$masses)
    xInterval <- c(xMin, xMax)
  } else {
    xMin <- xInterval[[1]]
    xMax <- xInterval[[2]]
  }
  
  massIntervalSelection <- dataList$masses >= xMin & dataList$masses <= xMax
  numberOfFragments <- dataList$numberOfFragments[massIntervalSelection]
  masses            <- dataList$masses[massIntervalSelection]
  
  minimumNumberOfFragments <- 5
  selection <- massIntervalSelection & (dataList$numberOfFragments >= minimumNumberOfFragments)
  
  numberOfFragments <- dataList$numberOfFragments[selection]
  masses            <- dataList$masses[selection]
  
  yMin <- 0
  if(sum(selection) == 0)
    yMax <- 1
  else
    yMax <- max(numberOfFragments)
  yInterval <- c(yMin, yMax)
  
  #colors <- rep(x = "black", times = length(dataList$masses))
  colors <- cmap(x = numberOfFragments, map = dataList$colorMapFragmentData)
  
  #########################################################################################################
  ## plot
  par(mar=c(5,3,2,0), mgp = c(2, 1, 0))  ## c(bottom, left, top, right)
  #plot(x = dataList$masses, y = dataList$numberOfFragments, xlab = "Fragment mass", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = "Fragment plot", col = colors, pch=19, cex=1., xaxt='n')
  #plot(x = masses, y = numberOfFragments, xlab = "", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  plot(x = NULL, y = NULL, xlab = "m/z", ylab = "MS/MS spectrum count", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  axis(side = 3)
  
  ## axis
  dataOrder <- order(numberOfFragments)
  #axis(side = 1, at = dataList$masses, labels = dataList$masses, las = 2, tick = TRUE, col = colors)
  #for(i in 1:length(dataList$masses))
  for(i in dataOrder)
    axis(side = 1, at = masses[[i]], labels = format(x = masses[[i]], digits = 0, nsmall = 4), las = 2, tick = TRUE, col.axis = colors[[i]])
  
  points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], type = "h", lwd=4)
  #points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], pch=19, cex=1.)
  
  resultObj <- list()
  resultObj$poiFragmentX <- masses
  resultObj$poiFragmentY <- numberOfFragments
  
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
calcPlotDendrogram <- function(dataList, filter, clusterDataList, annoPresentAnnotationsList, annoPresentColorsList, distanceMeasure, selectionFragmentTreeNodeSet = NULL, selectionAnalysisTreeNodeSet = NULL, selectionSearchTreeNodeSet = NULL, showClusterLabels, xInterval = NULL){
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## cluster
  par(mar=c(7.25,4,2,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  dend <- as.dendrogram(clusterDataList$cluster)
  
  ## remove labels left of the y-axis
  rightMostInvisibleLabelIndex <- floor(xInterval[[1]] - (xInterval[[2]] - xInterval[[1]]) * 0.04)
  if(rightMostInvisibleLabelIndex > 0){
    labelsToRemove <- clusterDataList$cluster$labels[clusterDataList$cluster$order][1:rightMostInvisibleLabelIndex]
    
    colLab <- colorLabels(clusterDataList$cluster$labels, NULL, NULL, labelsToRemove)
    dend <- dendrapply(dend, colLab)
  }
  
  ## color labels for search sub-roots
  if(!is.null(selectionSearchTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTree[selectionSearchTreeNodeSet[selectionSearchTreeNodeSet > 0]]), 
      -selectionSearchTreeNodeSet[selectionSearchTreeNodeSet < 0]
    )
    
    colLab <- colorLabels(clusterDataList$cluster$labels, clusterMembers, 'red')
    dend <- dendrapply(dend, colLab)
  }
  ## color labels for fragment sub-roots
  if(!is.null(selectionFragmentTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTree[selectionFragmentTreeNodeSet[selectionFragmentTreeNodeSet > 0]]), 
      -selectionFragmentTreeNodeSet[selectionFragmentTreeNodeSet < 0]
    )
    
    colLab <- colorLabels(clusterDataList$cluster$labels, clusterMembers, 'green')
    dend <- dendrapply(dend, colLab)
  }
  ## color labels for analysis sub-root
  if(!is.null(selectionAnalysisTreeNodeSet)){
    for(selectionAnalysisTreeNode in selectionAnalysisTreeNodeSet){
      clusterMembers <- NULL
      if(selectionAnalysisTreeNode > 0){
        clusterMembers <- c(clusterMembers, clusterDataList$innerNodeMembersTree[[selectionAnalysisTreeNode]])
      } else {
        clusterMembers <- c(clusterMembers, -selectionAnalysisTreeNode)
      }
    }
    
    colLab <- colorLabels(clusterDataList$cluster$labels, clusterMembers, 'blue')
    dend <- dendrapply(dend, colLab)
  }
  ## plot
  plot(x = dend, xlab = "", ylab = distanceMeasure, main = "Hierarchical cluster dendrogram", sub = "", xlim = xInterval)
  
  ## color tree for annotations
  resultObjTree <- analyzeTreeFromRootForAnnotations(dataList, cluster = clusterDataList$cluster, filter)
  innerNodeFeaturesAnnotations <- resultObjTree$innerNodeFeaturesAnnotations
  
  rootIndex <- length(clusterDataList$cluster$height)
  setOfColorSets <- setOfAnnotationSetsToSetOfColorSets(dataList, innerNodeFeaturesAnnotations)
  innerNodeMembersInner <- clusterDataList$innerNodeMembersInner
  innerNodeMembersTree <- clusterDataList$innerNodeMembersTree
  
  poisX <- clusterDataList$poiCoordinatesX
  poisY <- clusterDataList$poiCoordinatesY
  numberOfPois <- clusterDataList$numberOfPois
  numberOfPoisDrawn <- sum(clusterDataList$drawPoi)
  
  poisX <- poisX[clusterDataList$drawPoi]
  poisY <- poisY[clusterDataList$drawPoi]
  
  a2r_counter <<- 0
  numberOfInnerNodes <- as.integer(numberOfPois / 2)
  resultObjAnno <- getPrecursorColors(dataList = dataList, precursorSet = filter)
  leafColors <- resultObjAnno$setOfColors
  
  innerNodeColors      <<- vector(length = numberOfInnerNodes)
  innerNodeAnnotations <<- vector(length = numberOfInnerNodes)
  colorSubTreeForAnnotations(cluster = clusterDataList$cluster, index = rootIndex, innerNodeAnnotations = innerNodeFeaturesAnnotations, setOfColorSets = setOfColorSets, parentIndex = NULL, parentAnnotation = "Unknown", parentColor = "black")
  
  ## coloring of nodes by annotation
  pointSizesAnno  <- rep(x = clusterNodePointSize0, times = numberOfPoisDrawn)
  pointColorsAnno <- unlist(c(innerNodeColors, leafColors)[clusterDataList$drawPoi])
  
  ## selections
  pointsAnalysis <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionAnalysisTreeNodeSet)){
    indeces <- NULL
    for(selectionAnalysisTreeNode in selectionAnalysisTreeNodeSet){
      if(selectionAnalysisTreeNode < 0) indeces <- c(indeces, numberOfInnerNodes + abs(selectionAnalysisTreeNode))  ## leaf
      else                              indeces <- c(indeces, unlist(c(innerNodeMembersInner[[selectionAnalysisTreeNode]], numberOfInnerNodes + innerNodeMembersTree[[selectionAnalysisTreeNode]])))  ## inner node
    }
    pointsAnalysis[indeces] <- TRUE
  }
  pointsAnalysis <- pointsAnalysis[clusterDataList$drawPoi]
  
  pointsFragment <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionFragmentTreeNodeSet)){
    for(selectionFragmentTreeNode in selectionFragmentTreeNodeSet){
      if(selectionFragmentTreeNode < 0) idx <- numberOfInnerNodes + abs(selectionFragmentTreeNode)  ## leaf
      else                              idx <- unlist(c(innerNodeMembersInner[[selectionFragmentTreeNode]], numberOfInnerNodes + innerNodeMembersTree[[selectionFragmentTreeNode]]))  ## inner node
      pointsFragment[idx] <- TRUE
    }
  }
  pointsFragment <- pointsFragment[clusterDataList$drawPoi]
  
  pointsSearch <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionSearchTreeNodeSet)){
    for(selectionSearchTreeNode in selectionSearchTreeNodeSet){
      if(selectionSearchTreeNode < 0) idx <- numberOfInnerNodes + abs(selectionSearchTreeNode)  ## leaf
      else                              idx <- unlist(c(innerNodeMembersInner[[selectionSearchTreeNode]], numberOfInnerNodes + innerNodeMembersTree[[selectionSearchTreeNode]]))  ## inner node
      pointsSearch[idx] <- TRUE
    }
  }
  pointsSearch <- pointsSearch[clusterDataList$drawPoi]
  ## calc points
  resultObjPoints <- generatePoints(poisX, poisY, pointSizesAnno, pointColorsAnno, pointsAnalysis, pointsFragment, pointsSearch)
  pointSizes  <- resultObjPoints$pointSizes
  pointColors <- resultObjPoints$pointColors
  poisXpoints <- resultObjPoints$poisXpoints
  poisYpoints <- resultObjPoints$poisYpoints
  
  ## draw points
  points(x = poisXpoints, y = poisYpoints, col = pointColors, pch=19, cex=pointSizes)
  
  ## point labels
  if(showClusterLabels){
    poisXlabels <- poisX
    poisYlabels <- poisY
    
    plotRange <- par("usr") ## (x1, x2, y1, y2)
    yIntervalSize <- plotRange[[4]] - plotRange[[3]]
    poiLabelShift <- yIntervalSize / 50
    poiLabels <- clusterDataList$poiIntersectionSmooth[clusterDataList$drawPoi]
    #poiLabels <- clusterDataList$innerNodeFeaturesPresent[clusterDataList$drawPoi]
    graphics::text(x = poisXlabels - 0.0, y = poisYlabels + poiLabelShift, labels = poiLabels, pos = 4) ## pos: 1 below, 2 left, 3 above, 4 right
  }
  
  allAnnotations <- c(resultObjAnno$setOfAnnotations, innerNodeAnnotations)
  allColors      <- c(resultObjAnno$setOfColors, innerNodeColors)
  
  uniqueIndeces     <- which(!duplicated(allAnnotations))
  uniqueAnnotations <- allAnnotations[uniqueIndeces]
  uniqueColors      <- allColors[uniqueIndeces]
  
  resultList <- list(
    setOfAnnotations = uniqueAnnotations,
    setOfColors      = uniqueColors
  )
  
  return(resultList)
}
calcPlotHeatmap <- function(dataList, filterObj, clusterDataList, xInterval = NULL){
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## heatmap
  columnsOfInterest <- c(
    dataList$dataMeanColumnNameFunctionFromName(filterObj$groups[[1]]), dataList$dataMeanColumnNameFunctionFromName(filterObj$groups[[2]]), 
    dataList$lfcColumnNameFunctionFromName(filterObj$groups[[1]], filterObj$groups[[2]])
  )
  
  par(mar=c(0,4,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right) ## c(title, axis, label)
  
  if(length(columnsOfInterest) == 0){
    plot.new()
  } else {
    colorOne <- dataList$colorMatrixDataFrame[filterObj$filter, columnsOfInterest[[1]]][clusterDataList$cluster$order]
    colorTwo <- dataList$colorMatrixDataFrame[filterObj$filter, columnsOfInterest[[2]]][clusterDataList$cluster$order]
    colorLFC <- dataList$colorMatrixDataFrame[filterObj$filter, columnsOfInterest[[3]]][clusterDataList$cluster$order]
    
    plot(x = c(1, clusterDataList$numberOfPrecursorsFiltered), y = c(0, 3), type= "n", xlab = "", ylab = "", axes = FALSE, xlim = xInterval, ylim = c(0, 3))
    for(i in 1:clusterDataList$numberOfPrecursorsFiltered){
      rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 2, ytop = 3, col = colorLFC[[i]], border = NA)
      rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 1, ytop = 2, col = colorOne[[i]], border = NA)
      rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 0, ytop = 1, col = colorTwo[[i]], border = NA)
    }
    axis(side = 2, at = c(0.5, 1.5, 2.5), labels = c(filterObj$groups[[2]], filterObj$groups[[1]], "LFC"), las = 2, tick = TRUE)
  }
}
reorderAnnotationsForLegend <- function(annoLabels, annoColors){
  ignoreThere  <- any(annoLabels == "Ignore")
  unknownThere <- any(annoLabels == "Unknown")
  numberOfRealAnnotations <- length(annoLabels)
  if(ignoreThere){
    idx <- which(annoLabels == "Ignore")
    annoLabels <- c(annoLabels[-idx], annoLabels[idx])
    annoColors <- c(annoColors[-idx], annoColors[idx])
    numberOfRealAnnotations <- numberOfRealAnnotations - 1
  }
  if(unknownThere){
    idx <- which(annoLabels == "Unknown")
    annoLabels <- c(annoLabels[-idx], annoLabels[idx])
    annoColors <- c(annoColors[-idx], annoColors[idx])
    numberOfRealAnnotations <- numberOfRealAnnotations - 1
  }
  
  if(numberOfRealAnnotations > 0){
    order <- order(annoLabels[1:numberOfRealAnnotations])
    annoLabels[1:numberOfRealAnnotations] <- annoLabels[1:numberOfRealAnnotations][order]
    annoColors[1:numberOfRealAnnotations] <- annoColors[1:numberOfRealAnnotations][order]
  }
  
  resultList <- list(
    annoLabels = annoLabels,
    annoColors = annoColors
  )
}
calcPlotAnnoLegend <- function(annoLabels, annoColors){
  ## get and reorder annotations
  resultObj <- reorderAnnotationsForLegend(annoLabels, annoColors)
  annoLabels <- resultObj$annoLabels
  annoColors <- resultObj$annoColors
  
  calcPlotLegend(annoLabels, annoColors, "Annotations")
}
calcPlotScoresGroupsLegend <- function(groups, colors){
  ## get and reorder annotations
  calcPlotLegend(groups, colors, "Scores")
}
calcPlotLegend <- function(annoLabels, annoColors, title){
  ## layout
  numberOfLines <- length(annoLabels) + 1
  ySpacing <- 1 / (numberOfLines + 1)
  xSpacing <- 0.1
  
  labels <- c(paste(title, ":", sep = ""), annoLabels)
  xPositions <- c(-0.05, rep(x = xSpacing, times = length(annoLabels)))
  #yPositions <- seq(from = 1, to = 0, by = 1/(-length(xPositions)))
  yPositions <- seq(from = 1 - ySpacing, to = ySpacing, length.out = numberOfLines)
  
  symbolXPositions <- rep(x = xSpacing * 0.75, times = length(annoLabels))
  symbolYPositions <- yPositions[2:length(yPositions)]
  
  ##################################################################################################
  ## plot
  plotLegendWithBalls(labels, xPositions, yPositions, symbolXPositions, symbolYPositions, annoColors, xSpacing*0.075)
}
calcPlotAnnoLegendForImage <- function(annoLabels, annoColors, maximumNumberOfLines){
  ## get and reorder annotations
  resultObj <- reorderAnnotationsForLegend(annoLabels, annoColors)
  annoLabels <- resultObj$annoLabels
  annoColors <- resultObj$annoColors
  
  calcPlotLegendForImage(annoLabels, annoColors, "Annotations", maximumNumberOfLines)
}
calcPlotScoresGroupsLegendForImage <- function(groups, colors, maximumNumberOfLines){
  ## get and reorder annotations
  calcPlotLegendForImage(groups, colors, "Scores", maximumNumberOfLines)
}
calcPlotLegendForImage <- function(annoLabels, annoColors, title, maximumNumberOfLines){
  ## layout
  numberOfLines <- length(annoLabels) + 1
  ySpacing <- 1 / maximumNumberOfLines
  xSpacing <- 0.1
  
  labels <- c(paste(title, ":", sep = ""), annoLabels)
  xPositions <- c(-0.05, rep(x = xSpacing, times = length(annoLabels)))
  #yPositions <- seq(from = 1, to = 0, by = 1/(-length(xPositions)))
  yPositions <- seq(from = 1 - ySpacing, to = ySpacing, length.out = maximumNumberOfLines)
  yPositions <- yPositions[1:numberOfLines]
  
  symbolXPositions <- rep(x = xSpacing * 0.75, times = length(annoLabels))
  symbolYPositions <- yPositions[2:length(yPositions)]
  symbolYPositions <- symbolYPositions[1:length(symbolXPositions)]
  
  if(numberOfLines > maximumNumberOfLines){
    labels <- labels[1:(maximumNumberOfLines)]
    annoColors <- annoColors[1:(maximumNumberOfLines - 1)]
    xPositions <- xPositions[1:maximumNumberOfLines]
    yPositions <- yPositions[1:maximumNumberOfLines]
    
    labels[length(labels)] <- paste("...", (numberOfLines - maximumNumberOfLines + 1), " more", sep = "")
    annoColors <- annoColors[1:(length(annoColors) - 1)]
    symbolXPositions <- symbolXPositions[1:(length(symbolXPositions) - 1)]
    symbolYPositions <- symbolYPositions[1:(length(symbolYPositions) - 1)]
  }
  
  ##################################################################################################
  ## plot
  plotLegendWithBalls(labels, xPositions, yPositions, symbolXPositions, symbolYPositions, annoColors, xSpacing*0.075)
}
plotLegendWithBalls <- function(labels, xPositions, yPositions, circleXPositions, circleYPositions, annoColors, radius){
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  ## circles
  for(i in 1:length(annoColors))
    draw.circle(x = circleXPositions[[i]], y = circleYPositions[[i]], radius = radius, nv=50, border=annoColors[[i]], col = annoColors[[i]], lty=1, lwd=5)
  
  ## labels
  graphics::text(x = xPositions, y = yPositions, labels = labels, pos = 4)
}
calcPlotMS2Legend <- function(dataList){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  xSpacing <- 0.1
  stickLabels <- c(
    paste("Presence >= ", minimumProportionOfLeafs, sep = ""), 
    paste("Presence > ", minimumProportionToShowFragment, sep = ""), 
    "Selected fragment"
  )
  annoColorsStick <- c("black", "grey", "green")
  annoColorsBallSmall <- c("grey", "grey", "grey")
  annoColorsBallBig <- c("black", "black", "green")
  labels <- c("Fragment stick colors", stickLabels)
  xPositions <- c(-0.05, rep(x = xSpacing, times = length(stickLabels)))
  #yPositions <- seq(from = 1, to = 0, by = 1/(-length(xPositions)))
  yPositions <- seq(from = 0.9, to = -0.1, by = -0.25)[1:length(xPositions)]
  
  symbolXPositions <- rep(x = xSpacing * 0.75, times = length(stickLabels))
  symbolYPositions <- yPositions[2:length(yPositions)]
  
  ##################################################################################################
  ## plot
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  ## points
  pointSizes <- rep(x = ms2StickPointSizeInitial, times = length(stickLabels))
  pointSizesSmall <- rep(x = ms2StickPointSizeInitialSmall, times = length(stickLabels))
  
  delta <- 0.05
  segments(x0 = symbolXPositions, x1 = symbolXPositions, y0 = symbolYPositions - delta, y1 = symbolYPositions + delta, col = annoColorsStick, lwd = 3)
  points(x = symbolXPositions, y = symbolYPositions + delta, col = annoColorsBallBig, pch=19, cex=pointSizes)
  points(x = symbolXPositions, y = symbolYPositions + delta, col = annoColorsBallSmall, pch=19, cex=pointSizesSmall)
  
  ## labels
  graphics::text(x = xPositions, y = yPositions, labels = labels, pos = 4)
}
calcPlotDendrogramLegend <- function(){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  xSpacing <- 0.1
  stickLabels <- c(
    "Selection by analysis", 
    "Selection by fragment", 
    "Selection by search"
  )
  selectionColors <- c("blue", "green", "red")
  labels <- c("Selection marks", stickLabels)
  xPositions <- c(-0.05, rep(x = xSpacing, times = length(stickLabels)))
  #yPositions <- seq(from = 1, to = 0, by = 1/(-length(xPositions)))
  yPositions <- seq(from = 0.9, to = -0.1, by = -0.25)[1:length(xPositions)]
  
  symbolXPositions <- rep(x = xSpacing * 0.75, times = length(stickLabels))
  symbolYPositions <- yPositions[2:length(yPositions)]
  
  ##################################################################################################
  ## plot
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  ## points
  pointSizes      <- rep(x = clusterNodePointSize1, times = length(stickLabels))
  pointSizesSmall <- rep(x = clusterNodePointSize0, times = length(stickLabels))
  
  points(x = symbolXPositions, y = symbolYPositions, col = selectionColors, pch=19, cex=pointSizes)
  points(x = symbolXPositions, y = symbolYPositions, col = rep(x = "white", times = length(selectionColors)), pch=19, cex=pointSizesSmall)
  
  ## labels
  graphics::text(x = xPositions, y = yPositions, labels = labels, pos = 4)
}
calcPlotHeatmapLegend <- function(dataList){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  legend_imageAbs <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = dataList$logAbsMax,           to =  0,                         length.out = 100), map = dataList$colorMapAbsoluteData ), nrow=1)))
  legend_imageLFC <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = dataList$logFoldChangeMax, to = -dataList$logFoldChangeMax, length.out = 100), map = dataList$colorMapLogFoldChange), nrow=1)))
  
  minY <- 0.0
  maxY <- 0.9
  epsilon <- 0.075
  middle <- ((maxY - minY) / 2)
  #absMaxY <- 0.5-epsilon
  absMaxY <- middle - epsilon
  lfcMinY <- middle + epsilon
  
  ## abs
  maxExp <- as.integer(dataList$logAbsMax)
  exps <- 0:maxExp
  absLegendLabels    = paste("10^", exps, sep = "")
  #absLegendLabels    = vector(length = length(exps))
  #for(i in 1:length(exps)){
  #  print(paste('10', exps[[i]]))
  #  print(expression('10'^exps[[i]]))
  #  print(expression('10'^eval(exps[[i]])))
  #  print(bquote(10^.(exps[[i]])))
  #  print(expression(bquote(10^.(exps[[i]]))))
  #  #absLegendLabels[[i]] <- expression('10'^exps[[i]])
  #  absLegendLabels[[i]] <- expression('10'^"g")
  #}
  maxYlabelPosition <- absMaxY * (maxExp / dataList$logAbsMax)
  if(is.na(maxYlabelPosition))  ## data range is zero
    maxYlabelPosition <- absMaxY
  yLabelPositionIncrement <- maxYlabelPosition / (length(exps) - 1)
  if(is.na(yLabelPositionIncrement) | is.infinite(yLabelPositionIncrement)) ## data range is zero
    yLabelPositionIncrement <- maxYlabelPosition
  
  absLegendPositions <- seq(from = minY, to = maxYlabelPosition, by = yLabelPositionIncrement)
  
  ## lfc
  lfcLegendPositions <- seq(lfcMinY,maxY,l=5)
  lfcLegendLabels <- format(x = seq(-dataList$logFoldChangeMax, dataList$logFoldChangeMax, l=5), digits = 0)
  
  ##################################################################################################
  ## plot
  plot.new()
  plot.window(xlim = c(0, 3), ylim = c(0, 1))
  
  ## abs legend
  rasterImage(image = legend_imageAbs, xleft = 0, ybottom = minY, xright = 1, ytop = absMaxY)
  graphics::text(x = 2, y = absLegendPositions, labels = absLegendLabels)
  segments(
    x0  = rep(x = 0, times = length(absLegendPositions)),
    x1  = rep(x = 1.2, times = length(absLegendPositions)),
    y0  = absLegendPositions,
    y1  = absLegendPositions
  )
  segments(## lower hori; upper hori; left vert; right vert
    x0  = c(0, 0, 0, 1),
    x1  = c(1, 1, 0, 1),
    y0  = c(minY, absMaxY, minY   , minY   ),
    y1  = c(minY, absMaxY, absMaxY, absMaxY)
  )
  graphics::text(x = -0.2, y = absMaxY + 0.05, labels = "log10(MS abundance)", pos = 4)
  
  ## lfc legend
  rasterImage(image = legend_imageLFC, xleft = 0, ybottom = lfcMinY, xright = 1, ytop = maxY)
  graphics::text(x = 2, y = seq(lfcMinY,maxY,l=5), labels = lfcLegendLabels)
  segments(## lower hori; upper hori; left vert; right vert
    x0  = c(0, 0, 0, 1),
    x1  = c(1, 1, 0, 1),
    y0  = c(lfcMinY, maxY, lfcMinY, lfcMinY),
    y1  = c(lfcMinY, maxY, maxY   , maxY   )
  )
  segments(
    x0  = rep(x = 0, times = length(lfcLegendPositions)),
    x1  = rep(x = 1.2, times = length(lfcLegendPositions)),
    y0  = lfcLegendPositions,
    y1  = lfcLegendPositions
  )
  graphics::text(x = -0.2, y = maxY + 0.05, labels = "log2(MS fold change)", pos = 4)
}
calcPlotMS2 <- function(dataList, fragmentsX = c(), fragmentsY = c(), fragmentsColor = c(), fragmentsX_02 = NULL, fragmentsY_02 = NULL, fragmentsColor_02 = NULL, xInterval = NULL, selectedFragmentIndex = NULL){
  ####################
  ## fragment spectrum
  if(is.null(xInterval))
    xInterval <- c(dataList$minimumMass, dataList$maximumMass)
  
  ## abundances greater one
  fragmentsY[fragmentsY > 1] <- 1
  
  ## y-axis
  if(is.null(fragmentsX_02)){
    yInterval <- c(0, 1)
    #nodeColors <- rep(x = "black", times = length(fragmentsX))
    nodeColors <- fragmentsColor
    dataX <- fragmentsX
    dataY <- fragmentsY
    yTickPositions <- c(0, 0.25, 0.5, 0.75, 1)
    yTickLabels <- c(0, "", 0.5, "", 1)
  } else {
    #nodeColors <- rep(x = "black", times = length(fragmentsX) + length(fragmentsX_02))
    nodeColors <- c(fragmentsColor, fragmentsColor_02)
    
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
  plot(x = dataX, y = dataX, ylab = "Relative abundance", xlab = "m/z", xlim = xInterval, ylim = yInterval, xaxt='n', yaxt='n', col = nodeColors)
  axis(side = 2, at = yTickPositions, labels = yTickLabels)
  axis(side = 3)
  title("Fragment plot", line = 2)
  #mtext(side = 3, "m/z", line = 2)
  
  ## a-axis line
  if(!is.null(fragmentsX_02)){
    xIntervalSize <- xInterval[[2]] - xInterval[[1]]
    xl <- xInterval[[1]] - xIntervalSize
    xr <- xInterval[[2]] + xIntervalSize
    segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  }
  
  tickPositions <- dataX
  if(length(dataX) > 0){
    ## axis
    axis(side = 1, at = dataX, labels = FALSE, las = 2)
    axis(side = 1, at = tickPositions, labels = format(x = dataX, digits = 0, nsmall = 4), las = 2, tick = FALSE)
    
    ## points
    points(x = dataX, y = dataY, col = nodeColors, type = "h", lwd=4)
    
    pointSizes <- rep(x = ms2StickPointSizeInitial, times = length(dataX))
    pointSizesSmall <- rep(x = ms2StickPointSizeInitialSmall, times = length(dataX))
    pointColors <- rep(x = "black", times = length(dataX))
    pointColorsSmall <- rep(x = "gray", times = length(dataX))
    if(!is.null(selectedFragmentIndex)){
      pointSizes[[selectedFragmentIndex]] <- ms2StickPointSizeEmph
      pointSizesSmall[[selectedFragmentIndex]] <- ms2StickPointSizeEmphSmall
      pointColors[[selectedFragmentIndex]] <- "green"
      #pointColorsSmall[[selectedFragmentIndex]] <- "green"
    }
    
    points(x = dataX, y = dataY, col = pointColors, pch=19, cex=pointSizes)
    points(x = dataX, y = dataY, col = pointColorsSmall, pch=19, cex=pointSizesSmall)
  }
}
calcPlotPCAscores <- function(pcaObj, dataList, filterObj, pcaDimensionOne, pcaDimensionTwo, showScoresLabels, xInterval = NULL, yInterval = NULL){
  palette <- colorPalette()
  colorsForReplicates <- palette[unlist(lapply(X = filterObj$groups, FUN = function(x){ 
    groupIdx <- dataList$groupIdxFromGroupName(x)
    rep(x = groupIdx, times = length(dataList$dataColumnsNameFunctionFromName(x)))
  }))]
  colorsForGroups <- palette[unlist(lapply(X = filterObj$groups, FUN = dataList$groupIdxFromGroupName))]
  
  dataDimOne <- pcaObj$scores[, pcaDimensionOne]
  dataDimTwo <- pcaObj$scores[, pcaDimensionTwo]
  
  varianceOne <- format(x = pcaObj$variance[[pcaDimensionOne]], digits = 3)
  varianceTwo <- format(x = pcaObj$variance[[pcaDimensionTwo]], digits = 3)
  
  xAxisLabel  <- paste("t_", pcaDimensionOne, " (", varianceOne, "%)", sep = "")
  yAxisLabel  <- paste("t_", pcaDimensionTwo, " (", varianceTwo, "%)", sep = "")
  
  xMin <- min(dataDimOne)
  xMax <- max(dataDimOne)
  yMin <- min(dataDimTwo)
  yMax <- max(dataDimTwo)
  
  if(any(is.na(c(xMin, xMax, yMin, yMax)))){
    xMin <- -1
    xMax <- 1
    yMin <- -1
    yMax <- 1
  }
  
  if(is.null(xInterval))
    xInterval <- c(xMin, xMax)
  if(is.null(yInterval))
    yInterval <- c(yMin, yMax)
  
  par(mar=c(3 + 0.35, 3, 2, 1), mgp = c(2, 1, 0))  ## c(bottom, left, top, right)
  plot(x = dataDimOne, y = dataDimTwo, xlim = xInterval, ylim = yInterval, xlab = xAxisLabel, ylab = yAxisLabel, main = "Scores", col = colorsForReplicates, pch=19, cex=1.)
  
  ## axis
  xInt <- xMax - xMin
  yInt <- yMax - yMin
  xl <- xMin - xInt
  xr <- xMax + xInt
  yl <- yMin - yInt
  yr <- yMax + yInt
  segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  segments(x0 = 0, x1 = 0, y0 = yl, y1 = yr, col = "black", lwd = 1)
  
  if(showScoresLabels){
    labels <- dataList$dataColumnsNameFunctionFromNames(filterObj$groups)
    graphics::text(x = dataDimOne, y = dataDimTwo, labels = labels, pos = 4)
  }
  
  resultList <- list(
    groups = filterObj$groups,
    colors = colorsForGroups
  )
  return(resultList)
}
calcPlotPCAloadings <- function(pcaObj, dataList, filter, pcaDimensionOne, pcaDimensionTwo, selectionFragmentPcaLoadingSet = NULL, selectionAnalysisPcaLoadingSet = NULL, selectionSearchPcaLoadingSet = NULL, xInterval = NULL, yInterval = NULL, showLoadingsLabels = FALSE, showLoadingsAbundance = FALSE){
  varianceOne <- format(x = pcaObj$variance[[pcaDimensionOne]], digits = 3)
  varianceTwo <- format(x = pcaObj$variance[[pcaDimensionTwo]], digits = 3)
  
  xAxisLabel  <- paste("p_", pcaDimensionOne, " (", varianceOne, "%)", sep = "")
  yAxisLabel  <- paste("p_", pcaDimensionTwo, " (", varianceTwo, "%)", sep = "")
  
  dataDimOne <- pcaObj$loadings[, pcaDimensionOne]
  dataDimTwo <- pcaObj$loadings[, pcaDimensionTwo]
  
  xMin <- min(dataDimOne)
  xMax <- max(dataDimOne)
  yMin <- min(dataDimTwo)
  yMax <- max(dataDimTwo)
  
  if(any(is.na(c(xMin, xMax, yMin, yMax)))){
    xMin <- -1
    xMax <- 1
    yMin <- -1
    yMax <- 1
  }
  
  if(is.null(xInterval))
    xInterval <- c(xMin, xMax)
  if(is.null(yInterval))
    yInterval <- c(yMin, yMax)
  
  numberOfPrecursors <- length(dataDimOne)
  poisX <- dataDimOne
  poisY <- dataDimTwo
  
  pointSizesAnno  <- rep(x = clusterNodePointSize0, times = numberOfPrecursors)
  resultObjAnno <- getPrecursorColors(dataList, filter)
  pointColorsAnno <- resultObjAnno$setOfColors
  
  pointsAnalysis <- vector(mode = "logical", length = numberOfPrecursors)
  pointsAnalysis[selectionAnalysisPcaLoadingSet] <- TRUE
  pointsFragment <- vector(mode = "logical", length = numberOfPrecursors)
  pointsFragment[selectionFragmentPcaLoadingSet] <- TRUE
  pointsSearch <- vector(mode = "logical", length = numberOfPrecursors)
  pointsSearch[selectionSearchPcaLoadingSet] <- TRUE
  
  resultObjPoints <- generatePoints(poisX, poisY, pointSizesAnno, pointColorsAnno, pointsAnalysis, pointsFragment, pointsSearch)
  pointSizes    <- resultObjPoints$pointSizes
  pointColors   <- resultObjPoints$pointColors
  poisXpoints   <- resultObjPoints$poisXpoints
  poisYpoints   <- resultObjPoints$poisYpoints
  mappingToData <- resultObjPoints$mappingToData
  
  ## points
  #points(x = poisXpoints, y = poisYpoints, col = pointColors, pch=19, cex=pointSizes)
  if(showLoadingsAbundance){
    precursorMeansNorm <- dataList$dataFrameMeasurements[filter, "meanAll"]
    precursorMeansNorm <- precursorMeansNorm[mappingToData]
    pointSizes <- pointSizes * 2 * precursorMeansNorm
  }
  
  par(mar=c(3 + 0.35, 3, 2, 1), mgp = c(2, 1, 0))  ## c(bottom, left, top, right)
  #plot(x = dataDimOne, y = dataDimTwo, xlim = xInterval, ylim = yInterval, xlab = xAxisLabel, ylab = yAxisLabel, main = "Loadings", pch=19, cex=0.7, col = nodeColors)
  plot(x = NULL, y = NULL, xlim = xInterval, ylim = yInterval, xlab = xAxisLabel, ylab = yAxisLabel, main = "Loadings")
  points(x = poisXpoints, y = poisYpoints, col = pointColors, pch=19, cex=pointSizes)
  
  ## axis
  xInt <- xMax - xMin
  yInt <- yMax - yMin
  xl <- xMin - xInt
  xr <- xMax + xInt
  yl <- yMin - yInt
  yr <- yMax + yInt
  segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  segments(x0 = 0, x1 = 0, y0 = yl, y1 = yr, col = "black", lwd = 1)
  
  if(showLoadingsLabels){
    labels <- dataList$precursorLabels[filter]
    graphics::text(  x = poisX - 0.0, y = poisY + 0.0, labels = labels, pos = 4)
  }
  
  uniqueIndeces     <- which(!duplicated(resultObjAnno$setOfAnnotations))
  uniqueAnnotations <- resultObjAnno$setOfAnnotations[uniqueIndeces]
  uniqueColors      <- resultObjAnno$setOfColors[uniqueIndeces]
  
  resultList <- list(
    setOfAnnotations = uniqueAnnotations,
    setOfColors      = uniqueColors
  )
  return(resultList)
}
generatePoints <- function(poisX, poisY, pointSizesAnno, pointColorsAnno, pointsAnalysis, pointsFragment, pointsSearch){
  numberOfPoisDrawn <- length(poisX)
  
  ## analysis
  pointSizesAnalysis  <- vector(mode = "numeric", length = numberOfPoisDrawn)
  pointColorsAnalysis <- vector(length = numberOfPoisDrawn)
  pointSizesAnalysis[pointsAnalysis] <- clusterNodePointSize1
  pointColorsAnalysis[pointsAnalysis] <- "blue"
  
  ## fragment
  pointSizesFragment  <- vector(mode = "numeric", length = numberOfPoisDrawn)
  pointColorsFragment <- vector(length = numberOfPoisDrawn)
  intersection <- pointsAnalysis & pointsFragment
  difference   <- pointsFragment & (!intersection)
  pointSizesFragment[intersection] <- clusterNodePointSize2
  pointSizesFragment[difference] <- clusterNodePointSize1
  pointColorsFragment[pointsFragment] <- "green"
  
  ## search
  pointSizesSearch  <- vector(mode = "numeric", length = numberOfPoisDrawn)
  pointColorsSearch <- vector(length = numberOfPoisDrawn)
  intersection  <- pointsAnalysis & pointsFragment & pointsSearch
  intersection2 <- (pointsSearch & pointsAnalysis | pointsSearch & pointsFragment) & (!intersection)
  difference   <- pointsSearch & (!intersection) & (!intersection2)
  pointSizesSearch[intersection]  <- clusterNodePointSize3
  pointSizesSearch[intersection2] <- clusterNodePointSize2
  pointSizesSearch[difference] <- clusterNodePointSize1
  pointColorsSearch[pointsSearch] <- "red"
  
  pointSizes  <- c(pointSizesSearch[pointsSearch], pointSizesFragment[pointsFragment], pointSizesAnalysis[pointsAnalysis], pointSizesAnno)
  pointColors <- c(pointColorsSearch[pointsSearch], pointColorsFragment[pointsFragment], pointColorsAnalysis[pointsAnalysis], pointColorsAnno)
  poisXpoints <- c(poisX[pointsSearch], poisX[pointsFragment], poisX[pointsAnalysis], poisX)
  poisYpoints <- c(poisY[pointsSearch], poisY[pointsFragment], poisY[pointsAnalysis], poisY)
  
  mappingToData <- c(which(pointsSearch), which(pointsFragment), which(pointsAnalysis), 1:length(poisY))
  
  resultObj <- list(
    mappingToData  = mappingToData,
    pointSizes  = pointSizes,
    pointColors = pointColors,
    poisXpoints = poisXpoints,
    poisYpoints = poisYpoints
  )
  
  return(resultObj)
}
colorPalette2 <- function(){
  ## http://tools.medialab.sciences-po.fr/iwanthue/
  ## 20 colors
  ## or
  ## https://github.com/johnbaums/hues/blob/master/R/iwanthue.R
  ## library(colorBrewer)
  ## colorRampPalette(c("blue", "red"))( 4)
  ## palette(rainbow(6))
  
  #palette <- palette(c(
  palette <- c(
    rgb(184, 88,184, maxColorValue=255),
    rgb(102,178, 48, maxColorValue=255),
    rgb(220, 64, 59, maxColorValue=255),
    rgb( 78,167,149, maxColorValue=255),
    rgb(117, 79, 33, maxColorValue=255),
    rgb(131, 59, 93, maxColorValue=255),
    rgb(116,159,202, maxColorValue=255),
    rgb(176,158, 56, maxColorValue=255),
    rgb(114,119,221, maxColorValue=255),
    rgb(219,119, 40, maxColorValue=255),
    rgb( 72,101, 46, maxColorValue=255),
    rgb(213, 66,135, maxColorValue=255),
    rgb( 70, 91,112, maxColorValue=255),
    rgb(213,115,121, maxColorValue=255),
    rgb( 91, 76,141, maxColorValue=255),
    rgb(198,137,190, maxColorValue=255),
    rgb( 98,175,100, maxColorValue=255),
    rgb(146, 56, 45, maxColorValue=255),
    rgb(207, 78,219, maxColorValue=255),
    rgb(206,136, 82, maxColorValue=255)
  )
  #))
  return(palette)
}
colorPalette <- function(){
  palette <- c(
    "blue",
    "red",
    "yellow",
    "green",
    "brown",
    "deepskyblue",
    "orange",
    "deeppink",
    "aquamarine",##
    "burlywood", 
    "cadetblue",
    "coral",
    "cornflowerblue",
    "cyan",##
    "darkblue",
    "firebrick",
    "goldenrod",
    "indianred",
    "khaki",##
    "magenta",
    "maroon",
    "beige",
    "moccasin",
    "olivedrab",
    "orangered",
    "orchid",
    "paleturquoise3",##
    "rosybrown",
    "salmon",
    "seagreen3",
    "skyblue",
    "steelblue"
  )
  return(palette)
}
if(FALSE){
  colorPalette3 <- function(){
    ## http://tools.medialab.sciences-po.fr/iwanthue/
    ## 10 colors fancy light background
    ## or
    ## https://github.com/johnbaums/hues/blob/master/R/iwanthue.R
    ## library(colorBrewer)
    ## colorRampPalette(c("blue", "red"))( 4)
    ## palette(rainbow(6))
    
    palette <- rainbow(10)
    #palette <- c(
    #  rgb(230,196, 93, maxColorValue=255), # 8
    #  rgb(209,227,101, maxColorValue=255), # 1
    #  rgb(147,227,129, maxColorValue=255), # 6
    #  rgb(193,215,140, maxColorValue=255), # 9
    #  rgb(225,178,215, maxColorValue=255), # 2
    #  rgb( 96,222,221, maxColorValue=255), # 3
    #  rgb(237,168,130, maxColorValue=255), # 4
    #  rgb(108,227,175, maxColorValue=255), # 5
    #  rgb(165,205,225, maxColorValue=255), # 7
    #  rgb(160,214,181, maxColorValue=255)  # 0
    #)
    return(palette)
  }
}
## annotation stuff
precursorSetToSetOfAnnotationSets <- function(dataList, precursorSet){
  setOfAnnotationSets <- lapply(X = precursorSet, FUN = function(x){
    annotationSet <- dataList$annoArrayOfLists[[x]]
    if(dataList$annoArrayIsArtifact[[x]])
      annotationSet <- c(annotationSet, "Ignore")
    return(unlist(annotationSet))
  })
  return(setOfAnnotationSets)
}
setOfAnnotationSetsToSetOfColorSets <- function(dataList, setOfAnnotationSets){
  setOfColorSets <- lapply(X = setOfAnnotationSets, FUN = function(x){
    if(is.null(x))
      ## no annotation
      colors <- "black"
    else
      ## at least one annotation
      colors <- unlist(lapply(X = x, FUN = function(y){
        dataList$annoPresentColorsList[[match(x = y, table = dataList$annoPresentAnnotationsList)]] }
      ))
    
    return(colors)
  })
  return(setOfColorSets)
}
getPrecursorColors <- function(dataList, precursorSet){
  setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
  setOfColorSets      <- setOfAnnotationSetsToSetOfColorSets(dataList, setOfAnnotationSets)
  setOfColors <- lapply(X = setOfColorSets, FUN = function(x){
    if(any(x == "red"))
      ## at least one annotation is ignore --> take ignore
      color <- "red"
    else{
      switch(as.character(length(x)),
        "0"={## no annotation
          color <- "black"
        },
        "1"={## one annotation
          color <- x
        },
        {## multiple annotations --> take the one which is primary
          color <- x[[1]]
        }
      )## end switch
    }## end else
  })## end lapply
  setOfAnnotations <- lapply(X = setOfAnnotationSets, FUN = function(x){
    if(any(x == "Ignore"))
      ## at least one annotation is ignore --> take ignore
      annotation <- "Ignore"
    else{
      switch(as.character(length(x)),
        "0"={## no annotation
          annotation <- "Unknown"
        },
        "1"={## one annotation
          annotation <- x
        },
        {## multiple annotations --> take the one which is primary
          annotation <- x[[1]]
        }
      )## end switch
    }## end else
  })## end lapply
  
  resultList <- list(
    setOfColors      = unlist(setOfColors),
    setOfAnnotations = unlist(setOfAnnotations)
  )
  
  #return(unlist(setOfColors))
  return(resultList)
}
