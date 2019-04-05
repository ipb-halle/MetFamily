
#########################################################################################
## tree helper
analyzeTreeFromRoot <- function(dataList, cluster, filter){
  numberOfPrecursorsFiltered <- length(filter)
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  rootIndex <- length(cluster$height)
  
  ## create fields and compute stuff
  innerNodeHeightIncreasesHere <<- vector(mode = "logical", length = numberOfInnerNodes)
  innerNodeMembersTreeLeavesHere <<- list()
  innerNodeMembersPrecursorsHere <<- list()
  innerNodeMembersTreeClustersHere <<- list()
  innerNodeFeaturesIntersectionHere <<- list()
  innerNodeFeaturesUnionHere <<- list()
  innerNodeFeaturesCountsMatrixHere <<- sparseMatrix(i = numberOfInnerNodes, j = length(dataList$fragmentMasses), x = 0)
  innerNodeFeaturesPresentHere <<- list()
  #innerNodeFeaturesIntersectionCounterHere <<- vector(mode = "numeric", length = numberOfInnerNodes)
  #innerNodeFeaturesUnionCounterHere <<- vector(mode = "logical", length = numberOfInnerNodes)
  innerNodePositionHere <<- vector(mode = "numeric", length = numberOfInnerNodes)
  leafHeightsHere <<- vector(mode = "numeric", length = numberOfPrecursorsFiltered)
  innerNodeMembersPrecursorsHere[seq_len(numberOfPrecursorsFiltered)] <<- NA
  
  analyzeTree(dataList, cluster, filter, rootIndex)
  
  ## box
  resultObj <- list()
  resultObj$innerNodeHeightIncreases <- innerNodeHeightIncreasesHere
  resultObj$innerNodeMembersTreeLeaves <- innerNodeMembersTreeLeavesHere
  resultObj$innerNodeMembersPrecursors <- innerNodeMembersPrecursorsHere
  resultObj$innerNodeMembersTreeClusters <- innerNodeMembersTreeClustersHere
  resultObj$innerNodeFeaturesIntersection <- innerNodeFeaturesIntersectionHere
  resultObj$innerNodeFeaturesUnion <- innerNodeFeaturesUnionHere
  resultObj$innerNodeFeaturesCountsMatrix <- innerNodeFeaturesCountsMatrixHere
  resultObj$innerNodeFeaturesPresent <- innerNodeFeaturesPresentHere
  #resultObj$innerNodeFeaturesIntersectionCounter <- innerNodeFeaturesIntersectionCounterHere
  #resultObj$innerNodeFeaturesUnionCounter <- innerNodeFeaturesUnionCounterHere
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
    #featuresCounterIntersection <- featuresCounter
    #featuresCounterUnion        <- featuresCounter
    position <- match(x = leafIdx, table = cluster$order)
    featuresAnnotations <- dataList$annoArrayOfLists[[precursorIndex]]
    
    ## box
    resultObj <- list()
    resultObj$membersTreeLeaves <- leafIdx
    resultObj$membersPrecursors <- precursorIndex
    resultObj$membersTreeClusters <- NULL
    resultObj$featuresBinaryIntersection <- featuresBinaryIntersection
    resultObj$featuresBinaryUnion <- featuresBinaryUnion
    resultObj$featuresCounts <- featureCounts
    #resultObj$featuresCounterIntersection <- featuresCounterIntersection
    #resultObj$featuresCounterUnion        <- featuresCounterUnion
    resultObj$position <- position
    resultObj$featuresAnnotations <- featuresAnnotations
    
    return(resultObj)
  } else {
    ###################################
    ## inner node
    resultObj.l <- analyzeTree(dataList, cluster, filter, cluster$merge[nodeIdx, 1])
    resultObj.r <- analyzeTree(dataList, cluster, filter, cluster$merge[nodeIdx, 2])
    
    membersTreeLeaves <- c(resultObj.l$membersTreeLeaves, resultObj.r$membersTreeLeaves)
    membersPrecursors <- c(resultObj.l$membersPrecursors, resultObj.r$membersPrecursors)
    membersTreeClusters <- c(resultObj.l$membersTreeClusters, resultObj.r$membersTreeClusters, nodeIdx)
    featuresBinaryIntersection <- resultObj.l$featuresBinaryIntersection & resultObj.r$featuresBinaryIntersection
    featuresBinaryUnion        <- resultObj.l$featuresBinaryUnion        | resultObj.r$featuresBinaryUnion
    featuresCounts <- resultObj.l$featuresCounts + resultObj.r$featuresCounts
    #featuresCounterIntersection <- sum(featuresBinaryIntersection)
    #featuresCounterUnion        <- sum(featuresBinaryUnion)
    position <- (resultObj.l$position + resultObj.r$position) / 2
    featuresAnnotations <- intersect(x = unlist(resultObj.l$featuresAnnotations), y = unlist(resultObj.r$featuresAnnotations))
    
    ## box
    resultObj <- list()
    resultObj$membersTreeLeaves <- membersTreeLeaves
    resultObj$membersPrecursors <- membersPrecursors
    resultObj$membersTreeClusters <- membersTreeClusters
    resultObj$featuresBinaryIntersection <- featuresBinaryIntersection
    resultObj$featuresBinaryUnion <- featuresBinaryUnion
    resultObj$featuresCounts <- featuresCounts
    #resultObj$featuresCounterIntersection <- featuresCounterIntersection
    #resultObj$featuresCounterUnion <- featuresCounterUnion
    resultObj$position <- position
    resultObj$featuresAnnotations <- featuresAnnotations
    
    ## set values
    innerNodeMembersTreeLeavesHere[[nodeIdx]] <<- resultObj$membersTreeLeaves
    innerNodeMembersTreeClustersHere[[nodeIdx]] <<- resultObj$membersTreeClusters
    innerNodeMembersPrecursorsHere[[nodeIdx]] <<- resultObj$membersPrecursors
    innerNodePositionHere[[nodeIdx]] <<- resultObj$position
    innerNodeFeaturesIntersectionHere[[nodeIdx]] <<- which(resultObj$featuresBinaryIntersection)
    innerNodeFeaturesUnionHere[[nodeIdx]] <<- which(resultObj$featuresBinaryUnion)
    innerNodeFeaturesCountsMatrixHere[nodeIdx, ] <<- resultObj$featuresCounts
    innerNodeFeaturesPresentHere[[nodeIdx]] <<- sum(resultObj$featuresCounts / length(resultObj$membersPrecursors) >= minimumProportionOfLeafs)
    #innerNodeFeaturesIntersectionCounterHere[[nodeIdx]] <<- resultObj$featuresCounterIntersection
    #innerNodeFeaturesUnionCounterHere[[nodeIdx]] <<- resultObj$featuresCounterUnion
    
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
    
    return(resultObj)
  }
}
analyzeTreeFromRootForAnnotations <- function(dataList, cluster, filter){
  numberOfPrecursorsFiltered <- length(filter)
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  rootIndex <- length(cluster$height)
  
  ## create fields and compute stuff
  innerNodeFeaturesAnnotationsHere <<- list()
  
  analyzeTreeForAnnotations(dataList, cluster, filter, rootIndex)
  
  ## box
  resultObj <- list()
  resultObj$innerNodeFeaturesAnnotations <- innerNodeFeaturesAnnotationsHere
  
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
    
    resultObj <- list()
    if(any(length(featuresAnnotations) == 0, featuresAnnotations == "")){
      featuresAnnotations <- NULL
      resultObj$featuresAnnotations <- ""
    } else {
      resultObj$featuresAnnotations <- featuresAnnotations
      innerNodeFeaturesAnnotationsHere[[nodeIdx]] <<- resultObj$featuresAnnotations
    }
    
    return(resultObj)
  }
}
analyzeTreeForFrequentFragments <- function(clusterDataList, nodeIdx, minimumNumberOfChildren){
  if(nodeIdx < 0){
    ###################################
    ## leaf
    #leafIdx <- -nodeIdx
    #precursorIndex <- filter[leafIdx]
    #getMS2spectrumInfoForPrecursorLeaf()
    
    ## box
    return(data.frame("NodeIdx" = nodeIdx, "numberOfChildren" = 0, "numberOfFrequentFeatures" = NA))
  } else {
    ###################################
    ## inner node
    df.l <- analyzeTreeForFrequentFragments(clusterDataList, clusterDataList$cluster$merge[nodeIdx, 1], minimumNumberOfChildren)
    df.r <- analyzeTreeForFrequentFragments(clusterDataList, clusterDataList$cluster$merge[nodeIdx, 2], minimumNumberOfChildren)
    df <- rbind(df.l, df.r)
    numberOfChildren <- sum(df$"NodeIdx" < 0)
    
    if(numberOfChildren >= minimumNumberOfChildren & any(clusterDataList$innerNodeFeaturesCountsMatrix[nodeIdx, ] > (numberOfChildren * minimumProportionOfLeafs))){
      if(any(df$numberOfFrequentFeatures > 0, na.rm = TRUE))
        df <- df[-which(df$numberOfFrequentFeatures > 0), ]
      df[nrow(df)+1, ] <- c(nodeIdx, numberOfChildren, sum(clusterDataList$innerNodeFeaturesCountsMatrix[nodeIdx, ] > (numberOfChildren * minimumProportionOfLeafs)))
    }
    
    return(df)
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
  result <- list()
  if(index<0){ # it is a leaf
    leafIndex <- -index
    precursorIndex <- filter[[leafIndex]]
    criterionFulfilled <- yesNoFunction(precursorIndex)
    
    result$results <- NULL
    result$criterionFulfilled <- criterionFulfilled
  } else {
    result.l  <- getSetOfSubTrees(filter, clusterDataList, yesNoFunction, clusterDataList$cluster$merge[index, 1])
    result.r  <- getSetOfSubTrees(filter, clusterDataList, yesNoFunction, clusterDataList$cluster$merge[index, 2])
    
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
  }
  
  #print(paste(index, result$criterionFulfilled, paste(result$results, collapse = ";")))
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
    ## no parent or no annotation entry
    parentAnnotations <- NULL
  else
    ## there is annotation entry - possibly NULL
    parentAnnotations <- innerNodeAnnotations[[parentIndex]]
  
  ## current annotations
  if(length(innerNodeAnnotations) < index)
    ## no annotation entry
    currentAnnotations <- NULL
  else
    ## there is annotation entry - possibly NULL
    currentAnnotations <- innerNodeAnnotations[[index]]
  
  ## calculate the current color
  newAnnotations <- setdiff(x = currentAnnotations, y = parentAnnotations)
  
  if(length(newAnnotations) == 0){
    if(length(parentAnnotations) == 0){
      ## no annotations at all
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
  
  x0  <- c(out.l$x0,  c(x.l,   x.l,   x.r  ), out.r$x0 )
  x1  <- c(out.l$x1,  c(x.l,   x.r,   x.r  ), out.r$x1 )
  y0  <- c(out.l$y0,  c(h.l,   h.m,   h.r  ), out.r$y0 )
  y1  <- c(out.l$y1,  c(h.m,   h.m,   h.m  ), out.r$y1 )
  col <- c(out.l$col, c(color, color, color), out.r$col)
  
  #segments(
  #  x0  = c(x.l, x.l, x.r),
  #  x1  = c(x.l, x.r, x.r),
  #  y0  = c(h.l, h.m, h.r),
  #  y1  = c(h.m, h.m, h.m),
  #  col = color,
  #  lty = lty,
  #  lwd = lwd
  #)
  
  list(
    x   = x.m,
    x0  = x0,
    x1  = x1,
    y0  = y0,
    y1  = y1,
    col = col
  )
}
drawDendrogram <- function(cluster, index, lwd = 1, lty = 1){
  #########################################################################################
  ## leaf case
  if(index<0){ # it is a leaf
    a2r_counter <<- a2r_counter + 1
    return(list(
      x = a2r_counter
    ))       
  }
  
  #########################################################################################
  ## draw recursively
  h.m   <- cluster$height[index]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- cluster$merge[index,1]
  
  h.l <- ifelse(index.l<0, 0, cluster$height[index.l])
  
  out.l   <- drawDendrogram(cluster = cluster, index = index.l, lty=lty, lwd=lwd)
  x.l     <- out.l$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- cluster$merge[index,2]
  h.r <- ifelse(index.r<0, 0, cluster$height[index.r])
  out.r   <- drawDendrogram(cluster = cluster, index = index.r, lty=lty, lwd=lwd)
  x.r     <- out.r$x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  x.m  <- (x.r + x.l) / 2  
  
  x0  <- c(out.l$x0,  c(x.l,   x.l,   x.r  ), out.r$x0 )
  x1  <- c(out.l$x1,  c(x.l,   x.r,   x.r  ), out.r$x1 )
  y0  <- c(out.l$y0,  c(h.l,   h.m,   h.r  ), out.r$y0 )
  y1  <- c(out.l$y1,  c(h.m,   h.m,   h.m  ), out.r$y1 )
  #col <- c(out.l$col, c(color, color, color), out.r$col)
  
  #segments(
  #  x0  = c(x.l, x.l, x.r),
  #  x1  = c(x.l, x.r, x.r),
  #  y0  = c(h.l, h.m, h.r),
  #  y1  = c(h.m, h.m, h.m),
  #  col = color,
  #  lty = lty,
  #  lwd = lwd
  #)
  
  list(
    x   = x.m,
    x0  = x0,
    x1  = x1,
    y0  = y0,
    y1  = y1#,
    #col = col
  )
}
colorSubTreeForAnnotations2 <- function(cluster, index, dataList, filter, lwd = 1, lty = 1){
  #########################################################################################
  ## leaf case
  if(index<0){ # it is a leaf
    a2r_counter <<- a2r_counter + 1
    annos  <- unlist(dataList$annoArrayOfLists[filter][[-index]])
    #annos  <- dataList$annoArrayOfLists[[cluster$order[[-index]]]]
    #annos  <- dataList$annoArrayOfLists[[which(cluster$order == -index)]]
    colors <- dataList$annoPresentColorsList[unlist(dataList$annoPresentAnnotationsList) %in% annos]
    
    #print(paste(index, paste(colors, collapse = "-"), paste(annos, collapse = "-")))
    #print(paste("leaf", a2r_counter, index, paste(colors, collapse = "-"), cluster$order[[-index]], which(cluster$order == -index)))
    
    return(list(
      x = a2r_counter,
      x0  = vector(length = 0, mode = "numeric"),
      x1  = vector(length = 0, mode = "numeric"),
      y0  = vector(length = 0, mode = "numeric"),
      y1  = vector(length = 0, mode = "numeric"),
      col = vector(length = 0, mode = "character"),
      annos  = annos,
    #  annos  = "Unknown",
      colors = colors
    #  colors = "black"
    ))
  }
  
  #########################################################################################
  ## draw recursively
  h.m   <- cluster$height[index]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- cluster$merge[index,1]
  
  h.l <- if(index.l<0) 0 else cluster$height[index.l]
  
  out.l     <- colorSubTreeForAnnotations2(cluster = cluster, index = index.l, dataList = dataList, filter = filter, lty=lty, lwd=lwd)
  x.l       <- out.l$x
  annos.l   <- out.l$annos
  colors.l  <- out.l$colors
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- cluster$merge[index,2]
  h.r <- if(index.r<0) 0 else cluster$height[index.r]
  out.r     <- colorSubTreeForAnnotations2(cluster = cluster, index = index.r, dataList = dataList, filter = filter, lty=lty, lwd=lwd)
  x.r       <- out.r$x
  annos.r   <- out.r$annos
  colors.r  <- out.r$colors
  
  
  #########################################################################################
  ## determine color by annotations
  annos <- intersect(annos.r, annos.l)
  
  if(length(annos) == 0){
    ## no common annotations
    annos  <- "Unknown"
    colors <- "black"
    color  <- "black"
  } else {
    colors <- colors.l[annos.l %in% annos]
    color  <- colors[[1]]
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  x.m  <- (x.r + x.l) / 2  
  
  ## left from bottom to top; on top from left to right; right from bottom to top
  x0  <- c(out.l$x0,  c(x.l,   x.l,   x.r  ), out.r$x0 )
  x1  <- c(out.l$x1,  c(x.l,   x.r,   x.r  ), out.r$x1 )
  y0  <- c(out.l$y0,  c(h.l,   h.m,   h.r  ), out.r$y0 )
  y1  <- c(out.l$y1,  c(h.m,   h.m,   h.m  ), out.r$y1 )
  #col <- c(out.l$col, color, color, color, out.r$col)
  col <- c(out.l$col, color, out.r$col)
  
  #print(paste(
  #  index, index.l, index.r, color
  #  #length(out.l$col), length(out.r$col), length(col), 
  #  #paste(out.l$col, collapse = "-"), paste(out.r$col, collapse = "-"), 
  #  #paste(col, collapse = "-")
  #))
  
  #print(paste("cluster", index, "(", index.l, index.r, ")", color, x.l, x.m, x.r, h.l, h.m, h.r))
  
  #segments(
  #  x0  = c(x.l, x.l, x.r),
  #  x1  = c(x.l, x.r, x.r),
  #  y0  = c(h.l, h.m, h.r),
  #  y1  = c(h.m, h.m, h.m),
  #  col = color,
  #  lty = lty,
  #  lwd = lwd
  #)
  
  list(
    x   = x.m,
    x0  = x0,
    x1  = x1,
    y0  = y0,
    y1  = y1,
    col = col,
    annos  = annos,
    colors = colors
  )
}
