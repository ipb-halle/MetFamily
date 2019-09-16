
#########################################################################################
## constants
minimumProportionOfLeafs <- 0.75
minimumProportionToShowFragment <- 0.5

clusterNodePointSize0 <- 2/3
clusterNodePointSize1 <- 3/3
clusterNodePointSize2 <- 4/3
clusterNodePointSize3 <- 5/3

ms2StickPointSizeInitial <- 1.
ms2StickPointSizeInitialSmall <- 2/3.
#ms2StickPointSizeEmph <- 1.5
#ms2StickPointSizeEmphSmall <- 3/3.
ms2StickPointSizeEmph <- 1
ms2StickPointSizeEmphSmall <- 2/3.
ms2StickPointSizeMaximumMultiplier <- 0.75
dendrogramClusterPointSizeMaximumMultiplier <- 0.75
dendrogramHeatmapLeftMargin <- 10#6#4

#########################################################################################
## plotting
colorLabels <- function(labels, clusterMembers, color, labelsToRemove = NULL, newLabels = NULL){
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
      if(!is.null(newLabels))
        attr(n, "label") <- newLabels[[nodeIndexHere]] # new label
    }
    return(n)
  }
  
  return(colLab)
}

calcPlotDendrogram <- function(dataList, filter, clusterDataList, annoPresentAnnotationsList, annoPresentColorsList, distanceMeasure, selectionFragmentTreeNodeSet = NULL, selectionAnalysisTreeNodeSet = NULL, selectionSearchTreeNodeSet = NULL, showClusterLabels, hcaPrecursorLabels, xInterval = NULL){
  if(FALSE){
    dataList_ <<- dataList
    filter__ <<- filter
    clusterDataList_ <<- clusterDataList
    annoPresentAnnotationsList_ <<- annoPresentAnnotationsList
    annoPresentColorsList_ <<- annoPresentColorsList
    distanceMeasure_ <<- distanceMeasure
    selectionFragmentTreeNodeSet_ <<- selectionFragmentTreeNodeSet
    selectionAnalysisTreeNodeSet_ <<- selectionAnalysisTreeNodeSet
    selectionSearchTreeNodeSet_ <<- selectionSearchTreeNodeSet
    showClusterLabels_ <<- showClusterLabels
    hcaPrecursorLabels_ <<- hcaPrecursorLabels
    xInterval_ <<- xInterval
  }
  if(FALSE){
    dataList <<- dataList_
    filter <<- filter__
    clusterDataList <<- clusterDataList_
    annoPresentAnnotationsList <<- annoPresentAnnotationsList_
    annoPresentColorsList <<- annoPresentColorsList_
    distanceMeasure <<- distanceMeasure_
    selectionFragmentTreeNodeSet <<- selectionFragmentTreeNodeSet_
    selectionAnalysisTreeNodeSet <<- selectionAnalysisTreeNodeSet_
    selectionSearchTreeNodeSet <<- selectionSearchTreeNodeSet_
    showClusterLabels <<- showClusterLabels_
    hcaPrecursorLabels <<- hcaPrecursorLabels_
    xInterval <<- xInterval_
  }
  
  
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## hcaPrecursorLabels
  precursorLabels <- NULL
  switch(as.character(hcaPrecursorLabels),
         "m/z / RT"={
           precursorLabels <- clusterDataList$cluster$labels
         },
         "Metabolite name"={
           precursorLabels <- dataList$dataFrameInfos[filter, "Metabolite name"]
           
           maximumNumberOfCharacters <- 17
           tooLong <- nchar(x = precursorLabels) > maximumNumberOfCharacters
           precursorLabels[tooLong] <- paste(substring(text = precursorLabels[tooLong], first = 1, last = maximumNumberOfCharacters - 3), rep(x = "...", times = sum(tooLong)), sep = "")
         },
         "Metabolite family"={
           precursorLabels <- unlist(lapply(X = dataList$annoArrayOfLists[filter], FUN = function(x){
             if(length(x) == 0)
               return("Unknown")
             else
               return(x[[1]])
           }))
           
           maximumNumberOfCharacters <- 17
           tooLong <- nchar(x = precursorLabels) > maximumNumberOfCharacters
           precursorLabels[tooLong] <- paste(substring(text = precursorLabels[tooLong], first = 1, last = maximumNumberOfCharacters - 3), rep(x = "...", times = sum(tooLong)), sep = "")
         },
         {## unknown state
           stop(paste("Unknown hcaPrecursorLabels value", hcaPrecursorLabels))
         }
  )## end switch
  
  precursorLabelsWithIdx <- paste(precursorLabels, "_", seq_len(length(precursorLabels)), sep = "")
  
  ## remove labels left of the y-axis
  rightMostInvisibleLabelIndex <- floor(xInterval[[1]] - (xInterval[[2]] - xInterval[[1]]) * 0.04)
  if(rightMostInvisibleLabelIndex > 0){
    #labelsToRemove <- precursorLabelsWithIdx[clusterDataList$cluster$order][1:rightMostInvisibleLabelIndex]
    #length(na.omit(match(x = precursorLabelsWithIdx, table = labelsToRemove))) > 0
    labelIndecesToRemove <- clusterDataList$cluster$order[seq_len(rightMostInvisibleLabelIndex)]
    precursorLabelsWithIdx[labelIndecesToRemove] <- ""
    precursorLabels[labelIndecesToRemove] <- ""
  }
  
  clusterDataList$cluster$labels <- precursorLabelsWithIdx
  
  ####################
  ## cluster
  par(mar=c(7.25,dendrogramHeatmapLeftMargin,2,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  dend <- as.dendrogram(clusterDataList$cluster)
  
  ## color labels for search sub-roots
  if(!is.null(selectionSearchTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTreeLeaves[selectionSearchTreeNodeSet[selectionSearchTreeNodeSet > 0]]), 
      -selectionSearchTreeNodeSet[selectionSearchTreeNodeSet < 0]
    )
    
    colLab <- colorLabels(precursorLabelsWithIdx, clusterMembers, 'red')
    dend <- dendrapply(dend, colLab)
  }
  ## color labels for fragment sub-roots
  if(!is.null(selectionFragmentTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTreeLeaves[selectionFragmentTreeNodeSet[selectionFragmentTreeNodeSet > 0]]), 
      -selectionFragmentTreeNodeSet[selectionFragmentTreeNodeSet < 0]
    )
    
    colLab <- colorLabels(precursorLabelsWithIdx, clusterMembers, 'green')
    dend <- dendrapply(dend, colLab)
  }
  ## color labels for analysis sub-root
  if(!is.null(selectionAnalysisTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTreeLeaves[selectionAnalysisTreeNodeSet[selectionAnalysisTreeNodeSet > 0]]), 
      -selectionAnalysisTreeNodeSet[selectionAnalysisTreeNodeSet < 0]
    )
    
    colLab <- colorLabels(precursorLabelsWithIdx, clusterMembers, 'blue')
    dend <- dendrapply(dend, colLab)
  }
  #if(!is.null(selectionAnalysisTreeNodeSet)){
  #  for(selectionAnalysisTreeNode in selectionAnalysisTreeNodeSet){
  #    clusterMembers <- NULL
  #    if(selectionAnalysisTreeNode > 0){
  #      clusterMembers <- c(clusterMembers, clusterDataList$innerNodeMembersTreeLeaves[[selectionAnalysisTreeNode]])
  #    } else {
  #      clusterMembers <- c(clusterMembers, -selectionAnalysisTreeNode)
  #    }
  #  }
  #  
  #  colLab <- colorLabels(precursorLabelsWithIdx, clusterMembers, 'blue')
  #  dend <- dendrapply(dend, colLab)
  #}
  
  ### remove labels left of the y-axis
  #rightMostInvisibleLabelIndex <- floor(xInterval[[1]] - (xInterval[[2]] - xInterval[[1]]) * 0.04)
  #if(rightMostInvisibleLabelIndex > 0){
  #  labelsToRemove <- precursorLabelsWithIdx[clusterDataList$cluster$order][1:rightMostInvisibleLabelIndex]
  #  
  #  colLab <- colorLabels(precursorLabelsWithIdx, NULL, NULL, labelsToRemove)
  #  dend <- dendrapply(dend, colLab)
  #}
  
  ## remove precursorLabel indeces
  colLab <- colorLabels(precursorLabelsWithIdx, NULL, NULL, NULL, precursorLabels)
  dend <- dendrapply(dend, colLab)
  
  ## plot
  plot(x = dend, xlab = "", ylab = distanceMeasure, main = "Hierarchical cluster dendrogram", sub = "", xlim = xInterval)
  
  ## color tree for annotations
  resultObjTree <- analyzeTreeFromRootForAnnotations(dataList, cluster = clusterDataList$cluster, filter)
  innerNodeFeaturesAnnotations <- resultObjTree$innerNodeFeaturesAnnotations
  
  rootIndex <- length(clusterDataList$cluster$height)
  setOfColorSets <- setOfAnnotationSetsToSetOfColorSets(dataList, innerNodeFeaturesAnnotations)
  innerNodeMembersTreeClusters <- clusterDataList$innerNodeMembersTreeClusters
  innerNodeMembersTreeLeaves <- clusterDataList$innerNodeMembersTreeLeaves
  
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
      else                              indeces <- c(indeces, unlist(c(innerNodeMembersTreeClusters[[selectionAnalysisTreeNode]], numberOfInnerNodes + innerNodeMembersTreeLeaves[[selectionAnalysisTreeNode]])))  ## inner node
    }
    pointsAnalysis[indeces] <- TRUE
  }
  pointsAnalysis <- pointsAnalysis[clusterDataList$drawPoi]
  
  pointsFragment <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionFragmentTreeNodeSet)){
    for(selectionFragmentTreeNode in selectionFragmentTreeNodeSet){
      if(selectionFragmentTreeNode < 0) idx <- numberOfInnerNodes + abs(selectionFragmentTreeNode)  ## leaf
      else                              idx <- unlist(c(innerNodeMembersTreeClusters[[selectionFragmentTreeNode]], numberOfInnerNodes + innerNodeMembersTreeLeaves[[selectionFragmentTreeNode]]))  ## inner node
      pointsFragment[idx] <- TRUE
    }
  }
  pointsFragment <- pointsFragment[clusterDataList$drawPoi]
  
  pointsSearch <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionSearchTreeNodeSet)){
    for(selectionSearchTreeNode in selectionSearchTreeNodeSet){
      if(selectionSearchTreeNode < 0) idx <- numberOfInnerNodes + abs(selectionSearchTreeNode)  ## leaf
      else                              idx <- unlist(c(innerNodeMembersTreeClusters[[selectionSearchTreeNode]], numberOfInnerNodes + innerNodeMembersTreeLeaves[[selectionSearchTreeNode]]))  ## inner node
      pointsSearch[idx] <- TRUE
    }
  }
  pointsSearch <- pointsSearch[clusterDataList$drawPoi]
  
  ## cluster discriminativity
  clusterDiscriminativity <- clusterDataList$clusterDiscriminativity[clusterDataList$drawPoi]
  
  ## calc points
  resultObjPoints <- generatePoints(
    poisX = poisX, poisY = poisY, 
    pointSizesAnno = pointSizesAnno, pointColorsAnno = pointColorsAnno, 
    pointsAnalysis = pointsAnalysis, pointsFragment = pointsFragment, pointsSearch = pointsSearch,
    pointSizeModifier = clusterDiscriminativity
  )
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
calcPlotDendrogram_plotly <- function(
  dataList, filterObj, clusterDataList, 
  #annoPresentAnnotationsList, annoPresentColorsList, 
  distanceMeasure, 
  showClusterLabels, hcaPrecursorLabels, 
  selectionFragmentTreeNodeSet = NULL, selectionAnalysisTreeNodeSet = NULL, selectionSearchTreeNodeSet = NULL, 
  selectedSelection, heatmapContent, heatmapOrdering, heatmapProportion
){
  
  if(FALSE){
    dataList_ <<- dataList
    filterObj_1 <<- filterObj
    clusterDataList_ <<- clusterDataList
    #annoPresentAnnotationsList_ <<- annoPresentAnnotationsList
    #annoPresentColorsList_ <<- annoPresentColorsList
    distanceMeasure_ <<- distanceMeasure
    showClusterLabels_ <<- showClusterLabels
    hcaPrecursorLabels_ <<- hcaPrecursorLabels
    selectionFragmentTreeNodeSet_ <<- selectionFragmentTreeNodeSet
    selectionAnalysisTreeNodeSet_ <<- selectionAnalysisTreeNodeSet
    selectionSearchTreeNodeSet_ <<- selectionSearchTreeNodeSet
    selectedSelection_ <<- selectedSelection
    heatmapContent_ <<- heatmapContent
    heatmapOrdering_ <<- heatmapOrdering
    heatmapProportion_ <<- heatmapProportion
    #xInterval_ <<- xInterval
    #stop()
  }
  if(FALSE){
    dataList <- dataList_
    filterObj <- filterObj_1
    clusterDataList <- clusterDataList_
    #annoPresentAnnotationsList <- annoPresentAnnotationsList_
    #annoPresentColorsList <- annoPresentColorsList_
    distanceMeasure <- distanceMeasure_
    showClusterLabels <- showClusterLabels_
    hcaPrecursorLabels <- hcaPrecursorLabels_
    selectionFragmentTreeNodeSet <- selectionFragmentTreeNodeSet_
    selectionAnalysisTreeNodeSet <- selectionAnalysisTreeNodeSet_
    selectionSearchTreeNodeSet <- selectionSearchTreeNodeSet_
    selectedSelection <- selectedSelection_
    heatmapContent <- heatmapContent_
    heatmapOrdering <- heatmapOrdering_
    heatmapProportion <- heatmapProportion_
    #xInterval <- xInterval_
  }
  
  #stop("hihi")
  
  
  #if(is.null(xInterval))
  #  xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## hcaPrecursorLabels
  precursorLabels <- NULL
  switch(as.character(hcaPrecursorLabels),
         "m/z / RT"={
           precursorLabels <- clusterDataList$cluster$labels
         },
         "Metabolite name"={
           precursorLabels <- dataList$dataFrameInfos[filterObj$filter, "Metabolite name"]
           
           maximumNumberOfCharacters <- 17
           tooLong <- nchar(x = precursorLabels) > maximumNumberOfCharacters
           precursorLabels[tooLong] <- paste(substring(text = precursorLabels[tooLong], first = 1, last = maximumNumberOfCharacters - 3), rep(x = "...", times = sum(tooLong)), sep = "")
         },
         "Metabolite family"={
           precursorLabels <- unlist(lapply(X = dataList$annoArrayOfLists[filterObj$filter], FUN = function(x){
             if(length(x) == 0)
               return("Unknown")
             else
               return(x[[1]])
           }))
           
           maximumNumberOfCharacters <- 17
           tooLong <- nchar(x = precursorLabels) > maximumNumberOfCharacters
           precursorLabels[tooLong] <- paste(substring(text = precursorLabels[tooLong], first = 1, last = maximumNumberOfCharacters - 3), rep(x = "...", times = sum(tooLong)), sep = "")
         },
         {## unknown state
           stop(paste("Unknown hcaPrecursorLabels value", hcaPrecursorLabels))
         }
  )## end switch
  
  if(FALSE){
  precursorLabelsWithIdx <- paste(precursorLabels, "_", seq_len(length(precursorLabels)), sep = "")
  
  ## remove labels left of the y-axis
  rightMostInvisibleLabelIndex <- 0#floor(xInterval[[1]] - (xInterval[[2]] - xInterval[[1]]) * 0.04)
  if(rightMostInvisibleLabelIndex > 0){
    #labelsToRemove <- precursorLabelsWithIdx[clusterDataList$cluster$order][1:rightMostInvisibleLabelIndex]
    #length(na.omit(match(x = precursorLabelsWithIdx, table = labelsToRemove))) > 0
    labelIndecesToRemove <- clusterDataList$cluster$order[seq_len(rightMostInvisibleLabelIndex)]
    precursorLabelsWithIdx[labelIndecesToRemove] <- ""
    precursorLabels[labelIndecesToRemove] <- ""
  }
  
  clusterDataList$cluster$labels <- precursorLabelsWithIdx
  }
  ## poi labels
  poiLabels <- as.character(clusterDataList$poiIntersectionSmooth[clusterDataList$drawPoi])
  poiHovers <- vector(mode = "character", length = length(poiLabels))
  
  poiLabelsHere <- clusterDataList$poiLabels[clusterDataList$drawPoi]
  
  for(poiIdx in seq_along(poiLabels)){
    nodeLabel <- poiLabelsHere[[poiIdx]]
    resultObj <- getMS2spectrum(dataList = dataList, clusterDataList = clusterDataList, treeLabel = nodeLabel)
    
    if(nodeLabel < 0){ ## leaf
      featureFamilies <- dataList$annoArrayOfLists[[resultObj$precursorSet]]
      featureFamilies <- ifelse(
        test = length(featureFamilies) == 0, 
        yes = "None", 
        no = paste(unlist(featureFamilies), collapse = ", ")
      )
      
      poiHovers[[poiIdx]] <- paste(
        "Precursor: ",           dataList$precursorLabels[[resultObj$precursorSet]],           "<br>", 
        "Number of fragments: ", length(resultObj$fragmentMasses),                             "<br>", 
        "Name: ",                dataList$dataFrameInfos[[resultObj$precursorSet, "Metabolite name"]], "<br>", 
        "Metabolite families: ", featureFamilies,
        sep = "")
    } else { ## cluster
      
      featureFamilies <- Reduce(intersect, dataList$annoArrayOfLists[resultObj$precursorSet])
      featureFamilies <- ifelse(
        test = length(featureFamilies) == 0, 
        yes = "None", 
        no = paste(unlist(featureFamilies), collapse = ", ")
      )
      
      poiHovers[[poiIdx]] <- paste(
        "Cluster discriminating power: ", resultObj$clusterDiscriminativity,      "<br>", 
        "Precursors in cluster: ",        resultObj$numberOfPrecursors,           "<br>", 
        "Frequent fragments: ",           length(resultObj$fragmentMasses),       "<br>", 
        "Characteristic fragments: ",     sum(resultObj$fragmentColor == "black"),"<br>",
        "Common metabolite families: ",   featureFamilies,
        sep = "")
    }
    
    ## cluster
    #$fragmentMasses
    #$fragmentAbundances
    #$fragmentColor
    #$fragmentDiscriminativity
    #$clusterDiscriminativity
    #$infoText
    #$precursorSet
    #$numberOfPrecursors
    
    ## leaf
    #$fragmentMasses
    #$fragmentAbundances
    #$fragmentColor
    #$fragmentDiscriminativity
    #$infoText
    #$landingPageUrl
    #$precursorSet
    #$numberOfPrecursors
  }
  
  
  ####################
  ## cluster
  
  if(FALSE){
  par(mar=c(7.25,4,2,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  dend <- as.dendrogram(clusterDataList$cluster)
  }
  
  leafLabelColors <- rep(x = "black", times = length(precursorLabels))
  
  ## color labels for search sub-roots
  if(!is.null(selectionSearchTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTreeLeaves[selectionSearchTreeNodeSet[selectionSearchTreeNodeSet > 0]]), 
      -selectionSearchTreeNodeSet[selectionSearchTreeNodeSet < 0]
    )
    leafLabelColors[clusterMembers] <- 'red'
  }
  ## color labels for fragment sub-roots
  if(!is.null(selectionFragmentTreeNodeSet)){
    clusterMembers <- c(
      unlist(clusterDataList$innerNodeMembersTreeLeaves[selectionFragmentTreeNodeSet[selectionFragmentTreeNodeSet > 0]]), 
      -selectionFragmentTreeNodeSet[selectionFragmentTreeNodeSet < 0]
    )
    leafLabelColors[clusterMembers] <- 'green'
  }
  ## color labels for analysis sub-root
  if(!is.null(selectionAnalysisTreeNodeSet)){
    clusterMembers <- NULL
    for(selectionAnalysisTreeNode in selectionAnalysisTreeNodeSet){
      if(selectionAnalysisTreeNode > 0){
        clusterMembers <- c(clusterMembers, clusterDataList$innerNodeMembersTreeLeaves[[selectionAnalysisTreeNode]])
      } else {
        clusterMembers <- c(clusterMembers, -selectionAnalysisTreeNode)
      }
    }
    leafLabelColors[clusterMembers] <- 'blue'
  }
  
  ### remove labels left of the y-axis
  #rightMostInvisibleLabelIndex <- floor(xInterval[[1]] - (xInterval[[2]] - xInterval[[1]]) * 0.04)
  #if(rightMostInvisibleLabelIndex > 0){
  #  labelsToRemove <- precursorLabelsWithIdx[clusterDataList$cluster$order][1:rightMostInvisibleLabelIndex]
  #  
  #  colLab <- colorLabels(precursorLabelsWithIdx, NULL, NULL, labelsToRemove)
  #  dend <- dendrapply(dend, colLab)
  #}
  
  ## remove precursorLabel indeces
  #colLab <- colorLabels(precursorLabelsWithIdx, NULL, NULL, NULL, precursorLabels)
  #dend <- dendrapply(dend, colLab)
  
  ## color tree for annotations
  resultObjTree <- analyzeTreeFromRootForAnnotations(dataList, cluster = clusterDataList$cluster, filterObj$filter)
  innerNodeFeaturesAnnotations <- resultObjTree$innerNodeFeaturesAnnotations
  
  rootIndex <- length(clusterDataList$cluster$height)
  setOfColorSets <- setOfAnnotationSetsToSetOfColorSets(dataList, innerNodeFeaturesAnnotations)
  innerNodeMembersTreeClusters <- clusterDataList$innerNodeMembersTreeClusters
  innerNodeMembersTreeLeaves <- clusterDataList$innerNodeMembersTreeLeaves
  
  poisX <- clusterDataList$poiCoordinatesX
  poisY <- clusterDataList$poiCoordinatesY
  numberOfPois <- clusterDataList$numberOfPois
  numberOfPoisDrawn <- sum(clusterDataList$drawPoi)
  
  poisX <- poisX[clusterDataList$drawPoi]
  poisY <- poisY[clusterDataList$drawPoi]
  
  a2r_counter <<- 0
  numberOfInnerNodes <- as.integer(numberOfPois / 2)
  resultObjAnno <- getPrecursorColors(dataList = dataList, precursorSet = filterObj$filter)
  leafColors <- resultObjAnno$setOfColors
  
  innerNodeColors      <<- vector(length = numberOfInnerNodes)
  innerNodeAnnotations <<- vector(length = numberOfInnerNodes)
  segmentListAnno <- colorSubTreeForAnnotations(cluster = clusterDataList$cluster, index = rootIndex, innerNodeAnnotations = innerNodeFeaturesAnnotations, setOfColorSets = setOfColorSets, parentIndex = NULL, parentAnnotation = "Unknown", parentColor = "black")
  
  ## coloring of nodes by annotation
  numberOfPoisDrawn <- sum(clusterDataList$drawPoi)
  pointSizesAnno  <- rep(x = clusterNodePointSize0, times = numberOfPoisDrawn)
  pointColorsAnno <- unlist(c(innerNodeColors, leafColors)[clusterDataList$drawPoi])
  
  ############################
  ## dendrogram selections
  pointsAnalysis <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionAnalysisTreeNodeSet)){
    indeces <- NULL
    for(selectionAnalysisTreeNode in selectionAnalysisTreeNodeSet){
      if(selectionAnalysisTreeNode < 0) indeces <- c(indeces, numberOfInnerNodes + abs(selectionAnalysisTreeNode))  ## leaf
      else                              indeces <- c(indeces, unlist(c(innerNodeMembersTreeClusters[[selectionAnalysisTreeNode]], numberOfInnerNodes + innerNodeMembersTreeLeaves[[selectionAnalysisTreeNode]])))  ## inner node
    }
    pointsAnalysis[indeces] <- TRUE
  }
  pointsAnalysis <- pointsAnalysis[clusterDataList$drawPoi]
  
  pointsFragment <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionFragmentTreeNodeSet)){
    for(selectionFragmentTreeNode in selectionFragmentTreeNodeSet){
      if(selectionFragmentTreeNode < 0) idx <- numberOfInnerNodes + abs(selectionFragmentTreeNode)  ## leaf
      else                              idx <- unlist(c(innerNodeMembersTreeClusters[[selectionFragmentTreeNode]], numberOfInnerNodes + innerNodeMembersTreeLeaves[[selectionFragmentTreeNode]]))  ## inner node
      pointsFragment[idx] <- TRUE
    }
  }
  pointsFragment <- pointsFragment[clusterDataList$drawPoi]
  
  pointsSearch <- vector(mode = "logical", length = numberOfPois)
  if(!is.null(selectionSearchTreeNodeSet)){
    for(selectionSearchTreeNode in selectionSearchTreeNodeSet){
      if(selectionSearchTreeNode < 0) idx <- numberOfInnerNodes + abs(selectionSearchTreeNode)  ## leaf
      else                              idx <- unlist(c(innerNodeMembersTreeClusters[[selectionSearchTreeNode]], numberOfInnerNodes + innerNodeMembersTreeLeaves[[selectionSearchTreeNode]]))  ## inner node
      pointsSearch[idx] <- TRUE
    }
  }
  pointsSearch <- pointsSearch[clusterDataList$drawPoi]
  
  ## cluster discriminativity
  clusterDiscriminativity <- clusterDataList$clusterDiscriminativity[clusterDataList$drawPoi]
  
  ## calc points
  resultObjPoints <- generatePoints(
    poisX = poisX, poisY = poisY, 
    pointSizesAnno = pointSizesAnno, pointColorsAnno = pointColorsAnno, 
    pointsAnalysis = pointsAnalysis, pointsFragment = pointsFragment, pointsSearch = pointsSearch,
    pointSizeModifier = clusterDiscriminativity
  )
  pointSizes  <- resultObjPoints$pointSizes
  pointColors <- resultObjPoints$pointColors
  poisXpoints <- resultObjPoints$poisXpoints
  poisYpoints <- resultObjPoints$poisYpoints
  
  ####################
  ## heatmap
  
  ## heatmap data
  groups <- dataList$groups
  switch(heatmapContent,
         "Log-fold-change"={## log-fold-change
           columnsOfInterest <- c(
             dataList$dataMeanColumnNameFunctionFromName(filterObj$groups[[1]]), 
             dataList$dataMeanColumnNameFunctionFromName(filterObj$groups[[2]]), 
             dataList$lfcColumnNameFunctionFromName(filterObj$groups[[1]], filterObj$groups[[2]])
           )
           columnsOfInterestAbs <- columnsOfInterest[1:2]
           columnsOfInterestLFC <- columnsOfInterest[3]
           labels = c(filterObj$groups[[1]], filterObj$groups[[2]], "LFC")
           labelsAbs <- labels[1:2]
           labelsLFC <- labels[3]
         },
         "Abundance by group"={## groups
           columnsOfInterest <- unlist(lapply(X = groups, FUN = function(x){
             dataList$dataMeanColumnNameFunctionFromName(x)
           }))
           columnsOfInterest <- rev(columnsOfInterest) ## plot is bottom to top
           columnsOfInterestAbs <- columnsOfInterest
           columnsOfInterestLFC <- NULL
           labels <- groups
           labelsAbs <- labels[1:2]
           labelsLFC <- NULL
         },
         "Abundance by sample"={## samples
           columnsOfInterest <- dataList$dataColumnsNameFunctionFromGroupNames(groups = groups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
           columnsOfInterest <- dataList$orderColumnNames(groupSampleDataFrame = dataList$groupSampleDataFrame, columnNames = columnsOfInterest)
           columnsOfInterest <- rev(columnsOfInterest) ## plot is bottom to top
           columnsOfInterestAbs <- columnsOfInterest
           columnsOfInterestLFC <- NULL
           labels <- columnsOfInterest
           labelsAbs <- labels[1:2]
           labelsLFC <- NULL
         },
         {## unknown state
           stop(paste("Unknown heatmapContent value", heatmapContent))
         }
  )## end switch
  numberOfGroups <- length(columnsOfInterest)
  
  ## heatmap ordering
  switch(heatmapOrdering,
         "Specified order"={## no clustering
           doCalculateCluster <- FALSE
         },
         "MS1 clustering"={## do clustering
           if(heatmapContent == "Log-fold-change"){
             doCalculateCluster <- FALSE
           } else {
             doCalculateCluster <- TRUE
           }
         },
         {## unknown state
           stop(paste("Unknown heatmapOrdering value", heatmapOrdering))
         }
  )## end switch
  
  ## selections and selection colors
  #selectedTreeNodeSet <- NULL
  #frameColor <- NULL
  if(!is.null(selectedSelection)){
    analysisSelection <- {
      selectedTreeNodeSet <- selectionAnalysisTreeNodeSet
      frameColor <- "blue"
    }
    fragmentSelection <- {
      selectedTreeNodeSet <- selectionFragmentTreeNodeSet
      frameColor <- "green"
    }
    searchSelection <- {
      selectedTreeNodeSet <- selectionSearchTreeNodeSet
      frameColor <- "red"
    }
    switch(selectedSelection,
           "Analysis_HCA" = analysisSelection,
           "Analysis_PCA" = analysisSelection,
           "Fragment_HCA" = fragmentSelection,
           "Fragment_PCA" = fragmentSelection,
           "Search_HCA"   = searchSelection,
           "Search_PCA"   = searchSelection,
           {## unknown state
             stop(paste("Unknown selectedSelection value", selectedSelection))
           }
    )## end switch
  }
  ## in case of tree selections compute leaf intervals
  #if(!is.null(selectedTreeNodeSet)){
    clusterMemberList <- lapply(X = selectedTreeNodeSet, FUN = function(x){
      ifelse(test = x < 0, yes = -x, no = clusterDataList$innerNodeMembersTreeLeaves[x])[[1]]
    })
    positionsList <- lapply(X = clusterMemberList, FUN = function(x){
      match(x = x, table = clusterDataList$cluster$order)
    })
    intervals <- lapply(X = clusterMemberList, FUN = function(x){
      positions <- match(x = x, table = clusterDataList$cluster$order)
      c(min(positions), max(positions))
    })
    intervalMatrix <- as.matrix(as.data.frame(intervals))
    
    ## merge adjacent intervals TODO
    #intervalMatrix <- matrix(data=unlist(intervals), nrow=2,ncol=length(intervals))
    #if(ncol(intervalMatrix) > 1){
    #  intervalMatrix <- intervalMatrix[, order(intervalMatrix[1, ])]
    #  
    #  sapply(X = seq_len(ncol(intervalMatrix) - 1), FUN = function(x){
    #    intervalMatrix[2, x] + 1 == intervalMatrix[1, x+1]
    #  })
    #}
    
  #}
  
  ## MS1 clustering and colors
  if(doCalculateCluster){
    data_s_p <- t(as.matrix(as.data.frame(lapply(X = columnsOfInterest, FUN = function(x){
      values <- dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
      if(is.null(selectedTreeNodeSet))
        return(values)
      
      unlist(apply(X = intervalMatrix, MARGIN = 2, FUN = function(x){
        values[seq(from=x[[1]], to = x[[2]], by = 1)]
      }))
    }))))
    
    dist <- dist(data_s_p)
    cluster <- hclust(d = dist, method = "ward.D")
    
    ## optimal leaf ordering
    if(numberOfGroups > 2){
      opt <- order.optimal(dist = dist, merge = cluster$merge)
      cluster$merge <- opt$merge
      cluster$order <- opt$order
    }
    
    labels <- labels[cluster$order]
    columnsOfInterest <- columnsOfInterest[cluster$order]
    
    colors <- lapply(X = columnsOfInterest[cluster$order], FUN = function(x){
      dataList$colorMatrixDataFrame[filterObj$filter, x][clusterDataList$cluster$order]
    })
    values <- lapply(X = columnsOfInterest[cluster$order], FUN = function(x){
      dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
    })
    valuesAbs <- values
    valuesLFC <- NULL
  } else {
    #colors <- lapply(X = columnsOfInterest, FUN = function(x){
    #  dataList$colorMatrixDataFrame[filterObj$filter, x][clusterDataList$cluster$order]
    #})
    colors <- sapply(columnsOfInterest,FUN = function(x){
      dataList$colorMatrixDataFrame[filterObj$filter, x][clusterDataList$cluster$order]
    },simplify = FALSE,USE.NAMES = TRUE)
    values <- lapply(X = columnsOfInterest, FUN = function(x){
      dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
    })
    valuesAbs <- lapply(X = columnsOfInterestAbs, FUN = function(x){
      dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
    })
    valuesLFC <- lapply(X = columnsOfInterestLFC, FUN = function(x){
      dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
    })
  }
  
  ## assemble heatmap values and colors
  numberOfGroupsAbs <- length(valuesAbs)
  numberOfGroupsLFC <- length(valuesLFC)
  
  matrixAbs <- t(matrix(data = unlist(valuesAbs), nrow = clusterDataList$numberOfPrecursorsFiltered, ncol = numberOfGroupsAbs))
  matrixAbsLog <- unlist(log10(matrixAbs))
  matrixAbsLog[is.infinite(matrixAbsLog)] <- 0
  
  #paletteAbs <- colorRampPalette(rainbow(18)[10:1])
  paletteAbs <- rainbow(18)[10:1]
  paletteLFC <- c('blue', 'white', 'red')
  
  if(numberOfGroupsLFC > 0){
    matrixLFC <- t(matrix(data = unlist(valuesLFC), nrow = clusterDataList$numberOfPrecursorsFiltered, ncol = numberOfGroupsLFC))
    matrixLFC[is.infinite(matrixLFC)] <- 0
  } else {
    matrixLFC <- matrix(nrow = 0, ncol = clusterDataList$numberOfPrecursorsFiltered)
  }
  
  if(FALSE){
  quantiles <- c(-dataList$logFoldChangeMax, 0, dataList$logFoldChangeMax, dataList$logFoldChangeMax*1.01, 10^(seq(from=dataList$logAbsMax/10, to=dataList$logAbsMax, length.out=9)))
  quantiles <- approx(x = quantiles, n = 10000)$y
  colors    <- c('blue', 'white', 'red', rainbow(18)[10:1])
  #colors    <- c('blue', 'white', 'red', rainbow(50)[17:8])
  colors    <- colorRampPalette(colors)(10000)
  
  matrixAbsTmp <- matrixAbs
  matrixAbsTmp[matrixAbsTmp<=dataList$logFoldChangeMax] <- dataList$logFoldChangeMax*1.01
  valueMatrix <- rbind(matrixLFC, matrixAbsTmp)
  valueMatrix <- matrix(data = #quantiles[
    cut(as.vector(valueMatrix), breaks = quantiles, right = TRUE, labels = FALSE)
    #]
  , nrow=numberOfGroups, ncol=clusterDataList$numberOfPrecursorsFiltered)
  #valueMatrix <- valueMatrix[rev(seq_len(nrow(valueMatrix))), ]
  
  m <- rbind(matrixLFC, matrixAbsTmp)
  m <- valueMatrix
  colors2    <- c('blue', 'white', 'red', rainbow(18)[10:1])
  vals <- unlist(m)
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric(palette = colors2, domain = c(-dataList$logFoldChangeMax,10^dataList$logAbsMax))(vals)
  colz <- setNames(object = data.frame(vals[o], cols[o]), nm = NULL)
  
  
  colors    <- c('blue', 'white', 'red', rainbow(18)[10:1])
  
  minLFCnew <- 3/13*dataList$logAbsMax
  matrixLFCTmp <- (matrixLFC - dataList$logFoldChangeMax) / (2 * dataList$logFoldChangeMax)  * minLFCnew
  valueMatrix <- rbind(matrixLFCTmp, matrixAbsLog)
  #valueMatrix <- valueMatrix + minLFCnew
  }
  
  #################################
  ## dendrogram elements
  #subDivide <- function(items){
  #  ## add NA at every fourth position
  #  chunks <- split(items, cut(seq_along(items), ceiling(length(items) / 3), labels = FALSE))
  #  chunks <- lapply(X = chunks, FUN = function(x){
  #    c(x, NA)
  #  })
  #  unlist(chunks)
  #}
  
  #a2r_counter <<- 0
  #segmentListDend <- drawDendrogram(cluster = clusterDataList$cluster, index = rootIndex)
  #segmentListDend$x0 <- subDivide(segmentListDend$x0)
  #segmentListDend$x1 <- subDivide(segmentListDend$x1)
  #segmentListDend$y0 <- subDivide(segmentListDend$y0)
  #segmentListDend$y1 <- subDivide(segmentListDend$y1)
  
  a2r_counter <<- 0
  segmentListAnno <- colorSubTreeForAnnotations2(cluster = clusterDataList$cluster, index = rootIndex, dataList = dataList, filter = filterObj$filter)
  
  #segmentListAnno$x0 <- subDivide(segmentListAnno$x0)
  #segmentListAnno$x1 <- subDivide(segmentListAnno$x1)
  #segmentListAnno$y0 <- subDivide(segmentListAnno$y0)
  #segmentListAnno$y1 <- subDivide(segmentListAnno$y1)
  
  presentColors <- unique(segmentListAnno$col)
  #if("black" %in% presentColors)
  #  presentColors <- presentColors[-which(presentColors == "black")]
  
  ## list of lists
  annoByAnno <- sapply(X = presentColors, simplify = FALSE, FUN = function(color){
    indeces <- which(segmentListAnno$col == color)
    indeces <- unlist(sapply(X = indeces, simplify = FALSE, FUN = function(index){
      ((index - 1) * 3 + 1) : ((index - 1) * 3 + 3)
    }))
    list(
      x0=segmentListAnno$x0[indeces], 
      x1=segmentListAnno$x1[indeces], 
      y0=segmentListAnno$y0[indeces], 
      y1=segmentListAnno$y1[indeces]
    )
  })
  
  ##################################################################################
  ## plot
  curveNumberToCurveName <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(curveNumberToCurveName) <- c("curveNumber", "name")
  curveNumber <- -1 ## first curveNumber is 0
  
  ## create plot
  ## https://stackoverflow.com/questions/45861236/segments-prevent-the-display-of-tooltips-of-markers-in-r-plotly
  
  hoverlabel <- list(
    bgcolor = "grey",
    bordercolor = "black"
  )
  
  #p <- plot_ly(type = 'scatter', mode = 'lines')
  p <- plot_ly(source = "dendLabelsHeatmap")
  ## heatmap
  if(FALSE && numberOfGroupsLFC > 0){
    p <- add_heatmap(p, x = xCoordinates, y = labelsLFC, 
                     z = matrixLFC, name = 'heatmap_lfc', #type = 'heatmap', mode = NULL, 
                     #colors = palette(50),
                     colors = paletteLFC, zauto = FALSE, zmin = -dataList$logFoldChangeMax, zmax = dataList$logFoldChangeMax,
                     xaxis="x", yaxis="y4", showlegend = FALSE, showscale=FALSE
                     
    )
    curveNumberToCurveName[nrow(curveNumberToCurveName)+1, ] <- c(curveNumber <- curveNumber+1, "heatmap_lfc")
  }
  
  matrixForHeatmap <- rbind(
    matrixAbsLog,
    matrix(data=rep(x = 0, times = numberOfGroupsLFC * clusterDataList$numberOfPrecursorsFiltered), nrow=numberOfGroupsLFC, ncol=clusterDataList$numberOfPrecursorsFiltered)
  )
  p <- add_heatmap(p, x = xCoordinates, y = labels, 
                   z = matrixForHeatmap, name = 'heatmap_abs', #type = 'heatmap', mode = NULL, 
                   #colors = palette(50),
                   colors = paletteAbs, zauto = FALSE, zmin = 0, zmax = dataList$logAbsMax,
                   xaxis="x", yaxis="y3", showlegend = FALSE, showscale=FALSE
  )
  curveNumberToCurveName[nrow(curveNumberToCurveName)+1, ] <- c(curveNumber <- curveNumber+1, "heatmap_abs")
  if(FALSE){
    p <- add_heatmap(p, x = xCoordinates, y = labelsAbs, 
                     z = matrixAbs, name = 'abs', #type = 'heatmap', mode = NULL, 
                     #colors = palette(50),
                     colors = paletteAbs, zauto = FALSE, zmin = 0, zmax = 10^dataList$logAbsMax,
                     xaxis="x", yaxis="y3", showlegend = FALSE, showscale=FALSE
    )
    
    #p <- add_heatmap(p, x = xCoordinates, y = labels, 
    #                 z = valueMatrix, name = 'abs', #type = 'heatmap', mode = NULL, 
    #                 #colors = palette(50),
    #                 colors = colors, 
    #                 zauto = FALSE, zmin = 1, zmax = 10000,
    #                 #zauto = FALSE, zmin = -dataList$logFoldChangeMax, zmax = 10^dataList$logAbsMax,
    #                 xaxis="x", yaxis="y3", showlegend = FALSE#, showscale=FALSE
    #)
    p <- add_heatmap(p, x = xCoordinates, y = labels, 
                     z = m, name = 'abs', #type = 'heatmap', mode = NULL, 
                     #colors = palette(50),
                     color=colors,
                     #colorscale = colz, 
                     #zauto = FALSE, zmin = 1, zmax = 10000,
                     #zauto = FALSE, zmin = -dataList$logFoldChangeMax, zmax = 10^dataList$logAbsMax,
                     xaxis="x", yaxis="y3", showlegend = FALSE#, showscale=FALSE
    )
    #p <- add_heatmap(p, data=data.frame(x=xCoordinates, y=rep("bla",times=ncol(valueMatrix)),z=valueMatrix[1,],col=colors[valueMatrix[1,]], stringsAsFactors = F), x = ~x, y = ~y, 
    #                 z = ~z, name = 'abs', #type = 'heatmap', mode = NULL, 
    #                 #colors = palette(50),
    #                 colors = ~col, #zauto = FALSE, zmin = -dataList$logFoldChangeMax, zmax = 10^dataList$logAbsMax,
    #                 xaxis="x", yaxis="y3", showlegend = FALSE#, showscale=FALSE
    #)
    
    add_heatmap(p, x = xCoordinates, y = labels, 
                z = valueMatrix, name = 'abs', #type = 'heatmap', mode = NULL, 
                #colors = palette(50),
                colors = colors, #zauto = FALSE, zmin = -minLFCnew, zmax = 10^dataList$logAbsMax,
                xaxis="x", yaxis="y3", showlegend = FALSE#, showscale=FALSE
    )
  }
  
  
  ## draw dendrogram
  #print(paste("1", length(xDend), length(yDend), min(xDend, na.rm = TRUE), max(xDend, na.rm = TRUE)))
  #p <- add_trace(p, x = ~xDend, y = ~yDend, name = 'dend', mode = 'lines', 
  #               hoverinfo="none", line = list(color = "black")
  #)
  ## color dendrogram
  for(colorIdx in seq_along(annoByAnno)){
    curveName <- paste('dendAnno', presentColors[[colorIdx]], colorIdx, sep = "_")
    p <- add_segments(p, x = annoByAnno[[colorIdx]]$x0, y = annoByAnno[[colorIdx]]$y0, xend = annoByAnno[[colorIdx]]$x1, yend = annoByAnno[[colorIdx]]$y1,
                      name = curveName, 
                      mode = 'lines',
                      #mode = 'lines+markers',
                      #hoverinfo="skip", #"none", 
                      hoverinfo="text", hovertext=rep(x = "error", times = length(annoByAnno[[colorIdx]]$x0)), 
                      #hoverinfo="text", hovertext=seq_along(annoByAnno[[colorIdx]]$x0), 
                      #hoverinfo="text", hovertext=tooltip2, 
                      line = list(
                        color = presentColors[[colorIdx]]#,
                        
                        ## TODO for selections
                        #dash = "dash", ## ("solid", "dot", "dash", "longdash", "dashdot", or "longdashdot")
                        
                        #layer="below traces"
                        #hovermode = FALSE,
                        #layer="below"
                      ),
                      #captureevents = FALSE,
                      hoveron="fills",
                      xaxis="x", yaxis="y", showlegend = FALSE
    )
    curveNumberToCurveName[nrow(curveNumberToCurveName)+1, ] <- c(curveNumber <- curveNumber+1, curveName)
  }
  ## node labels
  if(FALSE){
  #print(paste("3", length(poisX), length(poisY), length(poiLabels), min(poisX), max(poisX)))
  xOffset <- 0.3#(xInterval[[2]] - xInterval[[1]] + 1) / 50
  yOffset <- max(clusterDataList$cluster$height) / 30
  p <- add_text(p, x = poisX + xOffset, y = poisY + yOffset, name = 'node labels', mode = 'text', 
                 hoverinfo="none", text=poiLabels, textposition="top right",#"middle center",#"middle right" #
                 xaxis="x", yaxis="y"
  )
  curveNumberToCurveName[nrow(curveNumberToCurveName)+1, ] <- c(curveNumber <- curveNumber+1, 'node labels')
  }
  
  ## nodes
  poiLabelsPoints <- c(rep(x = "", times = length(poisXpoints) - length(poiLabels)), poiHovers)
  #print(paste("4", length(poisXpoints), length(poisYpoints), length(pointColors), length(pointSizes), min(poisXpoints), max(poisXpoints)))
  p <- add_markers(p, x = poisXpoints, y = poisYpoints, name = 'nodes', #mode = 'markers', 
                 marker=list(
                   size = pointSizes*7., 
                   #sizeref = 1.,
                   color = pointColors,
                   opacity=1#,
                   #layer="above traces"
                   #layer="above"
                 ),
                 hoverinfo="text", hovertext=poiLabelsPoints, hoveron="points+fills", #, hoverlabel = list(bordercolor="blue")
                 hoverlabel = hoverlabel,
                 xaxis="x", yaxis="y", showlegend = FALSE
  )
  curveNumberToCurveName[nrow(curveNumberToCurveName)+1, ] <- c(curveNumber <- curveNumber+1, 'nodes')
  
  ## leaf labels
  if(FALSE){
  ## TODO repair
  xCoordinates <- seq_len(clusterDataList$numberOfPrecursorsFiltered)
  uniqueLeafLabelColors <- unique(leafLabelColors)
  leafLabelsList <- sapply(X = uniqueLeafLabelColors, FUN = function(x){
    which(leafLabelColors == x)
  })
  for(idx in seq_along(leafLabelsList)){
    leafLabelIndeces <- leafLabelsList[[idx]]
    leafLabelsName <- paste('leaf_labels', uniqueLeafLabelColors[[idx]], idx, sep = "_")
    p <- add_annotations(p, x = xCoordinates[leafLabelIndeces], y = 0, name = leafLabelsName, text = precursorLabels[leafLabelIndeces],
                    showarrow = FALSE, textangle = 270, font = list(color = uniqueLeafLabelColors[[idx]]),#bgcolor = uniqueLeafLabelColors[[idx]], 
                    xref = "x", yref = "y2"
    )
    curveNumberToCurveName[nrow(curveNumberToCurveName)+1, ] <- c(curveNumber <- curveNumber+1, leafLabelsName)
  }
  }
  
  ## plot config
  ## h 2.1 / l 3.4 / d 11
  ## 25 .. 25 .. 11 .. 11
  ## 1  ..  3 .. 10 ..  n
  labelsToDendrogramProportion <- 0.25
  dendrogramLabelsProportion <- 1 - heatmapProportion
  labelsProportion <- dendrogramLabelsProportion * labelsToDendrogramProportion
  dendrogramProportion <- dendrogramLabelsProportion * (1-labelsToDendrogramProportion)
  
  if(numberOfGroupsLFC > 0){
    heatmapProportionLFC <- heatmapProportion * numberOfGroupsLFC / numberOfGroups
    heatmapProportionAbs <- heatmapProportion - heatmapProportionLFC
  } else {
    heatmapProportionLFC <- 0
    heatmapProportionAbs <- heatmapProportion
  }
  
  p <- config(p, displayModeBar = FALSE)
  p <- plotly::layout(p=p,                        # all of layout's properties: /r/reference/#layout
              #showlegend = FALSE,
              #hovermode = F,#"closest",
              title = "Hierarchical cluster dendrogram",
              xaxis = list(
                visible = FALSE
                #scaleanchor = "TODO",
                #domain = xInterval
              ),
              yaxis = list(
                title = distanceMeasure,
                showgrid = FALSE,
                zeroline=FALSE,
                type = "linear",
                domain = c(heatmapProportion + labelsProportion, 1),
                rangemode="nonnegative"
              ),
              yaxis2 = list(
                #title = "TODO",
                #fixedrange = TRUE,
                showgrid = FALSE,
                zeroline=FALSE,
                type = "category",
                domain = c(heatmapProportion, heatmapProportion + labelsProportion)
              ),
              yaxis3 = list(
                title = "TODO",
                #fixedrange = TRUE,
                showgrid = FALSE,
                zeroline=FALSE,
                type = "category",
                domain = c(0, heatmapProportion)
              ),
              #yaxis4 = list(
              #  title = "TODO2",
              #  #fixedrange = TRUE,
              #  showgrid = FALSE,
              #  zeroline=FALSE,
              #  type = "category",
              #  domain = c(heatmapProportionAbs, heatmapProportion)
              #),
              ## heatmap selections
              #if(!is.null(selectedTreeNodeSet))
              #  apply(X = intervalMatrix, MARGIN = 2, FUN = function(x){
              #    rect(xleft = x[[1]] - 0.5, xright = x[[2]] + 0.5, ybottom = 0, ytop = numberOfGroups, border = frameColor, lwd = 2)
              #  })
              shapes = #list(
                c(
                ## LFC rects?
                {if(numberOfGroupsLFC > 0){
                  colorsHere <- colors[[columnsOfInterestLFC]]
                  lapply(X = seq_along(colorsHere), FUN = function(x){
                    list(type = "rect",
                         fillcolor = colorsHere[[x]], line = list(color="none", width=0), opacity = 1.,
                         x0 = x-0.5, x1 = x+0.5, xref = "x",
                         y0 = -0.5 + numberOfGroupsAbs, y1 = 0.5 + numberOfGroupsAbs, yref = "y3"
                    )
                  })
                  } else NULL
                },
                ## selections in the heatmap
                unname(apply(X = intervalMatrix, MARGIN = 2, FUN = function(x){
                  if(length(x) == 0)
                    return(list())
                  else
                    list(type = "rect",
                         fillcolor = "none", line = list(color=frameColor, width=2), opacity = 1.,
                         x0 = x[[1]]-0.5, x1 = x[[2]]+0.5, xref = "x",
                         y0 = -0.5, y1 = numberOfGroups - 0.5, yref = "y3"
                    )
                }))
                )
                #list(type = "rect",
                #     fillcolor = "none", line = list(color=frameColor, width=2), opacity = 1.,
                #     x0 = 3, x1 = 5, xref = "x",
                #     y0 = -0.5, y1 = numberOfGroupsAbs - 0.5, yref = "y2"
                #)
              #)
              #,
              #annotations = list(
              #  x = xCoordinates,
              #  y = m$mpg,
              #  text = rownames(m),
              #  xref = "x",
              #  yref = "y",
              #  showarrow = TRUE,
              #  arrowhead = 7,
              #  ax = 20,
              #  ay = -40
              #)
  )
  
  resultObj <- list()
  resultObj$curveNumberToCurveName <- curveNumberToCurveName
  resultObj$plotlyPlot             <- p
  resultObj$columnsOfInterest      <- columnsOfInterest
  
  return(resultObj)
  
  ## https://plot.ly/r/plotlyproxy/
  ## https://plot.ly/r/get-requests/
  ## https://cran.r-project.org/web/packages/plotly/plotly.pdf
  ## https://plot.ly/r/shinyapp-linked-click/
  ## https://plot.ly/r/reference/#Layout_and_layout_style_objects
  ## https://plot.ly/r/reference/#layout
  ## https://plot.ly/r/shinyapp-plotly-events/
  
  
  ## p.bars.text = p.bars %>%
  ## plotly::layout(annotations = list(x = 30, y = 0,  text = "Expected", textangle=270, showarrow=F, xanchor="center"))
  
  ## for PCA:
  ## zeroline = TRUE
  
}
calcPlotHeatmap <- function(dataList, filterObj, clusterDataList, selectedTreeNodeSet, frameColor, heatmapContent, heatmapOrdering, xInterval = NULL){
  
  if(FALSE){
    dataList_ <<- dataList
    filterObj_ <<- filterObj
    clusterDataList_ <<- clusterDataList
    selectedTreeNodeSet_ <<- selectedTreeNodeSet
    frameColor_ <<- frameColor
    heatmapContent_ <<- heatmapContent
    heatmapOrdering_ <<- heatmapOrdering
    xInterval_ <<- xInterval
  }
  if(FALSE){
    dataList <- dataList_
    filterObj <- filterObj_
    clusterDataList <- clusterDataList_
    selectedTreeNodeSet <- selectedTreeNodeSet_
    frameColor <- frameColor_
    heatmapContent <- heatmapContent_
    heatmapOrdering <- heatmapOrdering_
    xInterval <- xInterval_
  }
  
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## heatmap
  groups <- dataList$groups
  switch(heatmapContent,
         "Log-fold-change"={## log-fold-change
           columnsOfInterest <- c(
             dataList$dataMeanColumnNameFunctionFromName(filterObj$groups[[1]]), 
             dataList$dataMeanColumnNameFunctionFromName(filterObj$groups[[2]]), 
             dataList$lfcColumnNameFunctionFromName(filterObj$groups[[1]], filterObj$groups[[2]])
           )
           labels = c(filterObj$groups[[1]], filterObj$groups[[2]], "LFC")
         },
         "Abundance by group"={## groups
           columnsOfInterest <- unlist(lapply(X = groups, FUN = function(x){
             dataList$dataMeanColumnNameFunctionFromName(x)
           }))
           #columnsOfInterest <- rev(columnsOfInterest) ## plot is bottom to top
           labels <- groups
         },
         "Abundance by sample"={## samples
           columnsOfInterest <- dataList$dataColumnsNameFunctionFromGroupNames(groups = groups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
           columnsOfInterest <- dataList$orderColumnNames(groupSampleDataFrame = dataList$groupSampleDataFrame, columnNames = columnsOfInterest)
           #columnsOfInterest <- rev(columnsOfInterest) ## plot is bottom to top
           labels <- columnsOfInterest
         },
         {## unknown state
           stop(paste("Unknown heatmapContent value", heatmapContent))
         }
  )## end switch
  numberOfGroups <- length(columnsOfInterest)
  
  ## heatmap ordering
  switch(heatmapOrdering,
         "Specified order"={## no clustering
           doCalculateCluster <- FALSE
         },
         "MS1 clustering"={## do clustering
           if(heatmapContent == "Log-fold-change"){
             doCalculateCluster <- FALSE
           } else {
             doCalculateCluster <- TRUE
           }
         },
         {## unknown state
           stop(paste("Unknown heatmapOrdering value", heatmapOrdering))
         }
  )## end switch
  
  ## in case of tree selections compute leaf intervals
  if(!is.null(selectedTreeNodeSet)){
    clusterMemberList <- lapply(X = selectedTreeNodeSet, FUN = function(x){
      ifelse(test = x < 0, yes = -x, no = clusterDataList$innerNodeMembersTreeLeaves[x])[[1]]
    })
    positionsList <- lapply(X = clusterMemberList, FUN = function(x){
      match(x = x, table = clusterDataList$cluster$order)
    })
    intervals <- lapply(X = clusterMemberList, FUN = function(x){
      positions <- match(x = x, table = clusterDataList$cluster$order)
      c(min(positions), max(positions))
    })
    intervalMatrix <- as.matrix(as.data.frame(intervals))
  }
  
  ## MS1 clustering and colors
  if(doCalculateCluster){
    data_s_p <- t(as.matrix(as.data.frame(lapply(X = columnsOfInterest, FUN = function(x){
      values <- dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
      if(is.null(selectedTreeNodeSet))
        return(values)
      
      unlist(apply(X = intervalMatrix, MARGIN = 2, FUN = function(x){
        values[seq(from=x[[1]], to = x[[2]], by = 1)]
      }))
    }))))
    
    dist <- dist(data_s_p)
    cluster <- hclust(d = dist, method = "ward.D")
    
    ## optimal leaf ordering
    if(numberOfGroups > 2){
      opt <- cba::order.optimal(dist = dist, merge = cluster$merge)
      cluster$merge <- opt$merge
      cluster$order <- opt$order
    }
    
    labels            <- labels[cluster$order]
    columnsOfInterest <- columnsOfInterest[cluster$order]
    
    #colors <- lapply(X = columnsOfInterest[cluster$order], FUN = function(x){
    colors <- lapply(X = columnsOfInterest, FUN = function(x){
      dataList$colorMatrixDataFrame [filterObj$filter, x][clusterDataList$cluster$order]
    })
    #values <- lapply(X = columnsOfInterest[cluster$order], FUN = function(x){
    values <- lapply(X = columnsOfInterest, FUN = function(x){
      dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
    })
    
    columnOrder <- cluster$order
  } else {
    colors <- lapply(X = columnsOfInterest, FUN = function(x){
      dataList$colorMatrixDataFrame [filterObj$filter, x][clusterDataList$cluster$order]
    })
    values <- lapply(X = columnsOfInterest, FUN = function(x){
      dataList$dataFrameMeasurements[filterObj$filter, x][clusterDataList$cluster$order]
    })
    columnOrder <- seq_along(columnsOfInterest)
  }
  
  ## plot
  par(mar=c(0,dendrogramHeatmapLeftMargin,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right) ## c(title, axis, label)
  
  if(numberOfGroups == 0){
    ## nothing there
    print("### calcPlotHeatmap: no columnsOfInterest")
    plot.new()
    return()
  }
  
  if(TRUE){
  plot(x = c(1, clusterDataList$numberOfPrecursorsFiltered), y = c(0, numberOfGroups), type= "n", xlab = "", ylab = "", axes = FALSE, xlim = xInterval, ylim = c(0, numberOfGroups))
  rect(
    xleft   = rep(x = seq_len(clusterDataList$numberOfPrecursorsFiltered) - 0.5, times = numberOfGroups), 
    xright  = rep(x = seq_len(clusterDataList$numberOfPrecursorsFiltered) + 0.5, times = numberOfGroups), 
    ybottom = unlist(lapply(X = seq_len(numberOfGroups) - 1, FUN = function(x){rep(x = x, times = clusterDataList$numberOfPrecursorsFiltered)})),
    ytop    = unlist(lapply(X = seq_len(numberOfGroups)    , FUN = function(x){rep(x = x, times = clusterDataList$numberOfPrecursorsFiltered)})),
    col     = unlist(colors),
    border = NA
  )
  if(!is.null(selectedTreeNodeSet))
    apply(X = intervalMatrix, MARGIN = 2, FUN = function(x){
      rect(xleft = x[[1]] - 0.5, xright = x[[2]] + 0.5, ybottom = 0, ytop = numberOfGroups, border = frameColor, lwd = 2)
    })
  
  axis(side = 2, at = seq(from = 0.5, by = 1, length.out = numberOfGroups), labels = labels, las = 2, tick = TRUE)
  }
  
  if(FALSE){
    m <- matrix(data = unlist(values), nrow = clusterDataList$numberOfPrecursorsFiltered, ncol = numberOfGroups)
    col <- m[1, ]
    names(col) <- labels
    sort(col)
    colUnSorted <- col
  }
  
  if(FALSE){
    plot_ly(
      x = NULL, y = labels,
      z = t(matrix(data = unlist(colors), nrow = clusterDataList$numberOfPrecursorsFiltered, ncol = numberOfGroups)), 
      type = "heatmap"
    )
    
    m <- t(matrix(data = rgb(t(col2rgb(unlist(colors))), maxColorValue=255), nrow = numberOfPrecursorsFiltered, ncol = numberOfGroups))
    m <- matrix(data = c(1:6), nrow = 2, ncol = 3)
    
    
    palette <- colorRampPalette(rainbow(18)[10:1])
    vals <- unique(scales::rescale(c(m)))
  }
  
  if(FALSE){
  
  m <- log10(t(matrix(data = unlist(values), nrow = clusterDataList$numberOfPrecursorsFiltered, ncol = numberOfGroups)))
  m[is.infinite(m)] <- 0
  
  min <- min(m, na.rm = TRUE)
  max <- max(m, na.rm = TRUE)
  
  vals <- unlist(m)
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric(palette = rainbow(18)[10:1], domain = c(min,max))(vals)
  colz <- setNames(object = data.frame(vals[o], cols[o]), nm = NULL)
  
  plot <- plot_ly(
    x = NULL, y = labels,
    z = m, colorscale = colz,# colors = colorRamp(colors = c("red", "green")),#palette,
    type = "heatmap"
  ) %>%
    plotly::layout(                        # all of layout's properties: /r/reference/#layout
      xaxis = list(
        visible = FALSE,
        showgrid = FALSE#,
        #scaleanchor = "TODO",
        #scaleratio = 1
        #domain = xInterval
      ),
      yaxis = list(
        fixedrange = TRUE,
        showgrid = FALSE,
        type = "category"#,
        #tickmode = "auto", #"array"
        #tickvals = 1:numberOfGroups,
        #nticks = numberOfGroups,
        #ticktext = TODO
        #ticklen = -1
        #zeroline = TRUE ## for pca
      )
      #title = "Unemployment", # layout's title: /r/reference/#layout-title
      #xaxis = list(           # layout's xaxis is a named list. List of valid keys: /r/reference/#layout-xaxis
      #  title = "Time",      # xaxis's title: /r/reference/#layout-xaxis-title
      #  showgrid = F),       # xaxis's showgrid: /r/reference/#layout-xaxis-showgrid
      #yaxis = list(           # layout's yaxis is a named list. List of valid keys: /r/reference/#layout-yaxis
      #  title = "uidx")     # yaxis's title: /r/reference/#layout-yaxis-title
    )
  
  ## https://plot.ly/r/plotlyproxy/
  ## https://plot.ly/r/get-requests/
  ## https://cran.r-project.org/web/packages/plotly/plotly.pdf
  ## https://plot.ly/r/shinyapp-linked-click/
  ## https://plot.ly/r/reference/#Layout_and_layout_style_objects
  ## https://plot.ly/r/reference/#layout
  ## https://plot.ly/r/shinyapp-plotly-events/
  return(list(plot=plot, columnsOfInterest=columnsOfInterest))
  }
  
  returnObj <- list()
  returnObj$columnsOfInterest <- columnsOfInterest
  returnObj$columnOrder <- columnOrder
  
  return(returnObj)
}
calcPlotHeatmapOld <- function(dataList, filterObj, clusterDataList, xInterval = NULL){
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
    rect(
      xleft   = rep(x = seq_len(clusterDataList$numberOfPrecursorsFiltered) - 0.5, times = 3), 
      xright  = rep(x = seq_len(clusterDataList$numberOfPrecursorsFiltered) + 0.5, times = 3), 
      ybottom = unlist(lapply(X = 0:2, FUN = function(x){rep(x = x, times = clusterDataList$numberOfPrecursorsFiltered)})),
      ytop    = unlist(lapply(X = 1:3, FUN = function(x){rep(x = x, times = clusterDataList$numberOfPrecursorsFiltered)})),
      col     = c(colorTwo, colorOne, colorLFC),
      border = NA
    )
    
    #for(i in seq_len(clusterDataList$numberOfPrecursorsFiltered)){
    #  rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 2, ytop = 3, col = colorLFC[[i]], border = NA)
    #  rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 1, ytop = 2, col = colorOne[[i]], border = NA)
    #  rect(xleft = i - 0.5, xright = i + 0.5, ybottom = 0, ytop = 1, col = colorTwo[[i]], border = NA)
    #}
    axis(side = 2, at = c(0.5, 1.5, 2.5), labels = c(filterObj$groups[[2]], filterObj$groups[[1]], "LFC"), las = 2, tick = TRUE)
  }
  
  return(columnsOfInterest)
}
reorderAnnotationsForLegend <- function(annoLabels, annoColors){
  if(FALSE){
    annoLabels_ <<- annoLabels
    annoColors_ <<- annoColors
  }
  if(FALSE){
    annoLabels <<- annoLabels_
    annoColors <<- annoColors_
  }
  
  
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
    order <- order(annoLabels[seq_len(numberOfRealAnnotations)])
    annoLabels[seq_len(numberOfRealAnnotations)] <- annoLabels[seq_len(numberOfRealAnnotations)][order]
    annoColors[seq_len(numberOfRealAnnotations)] <- annoColors[seq_len(numberOfRealAnnotations)][order]
  }
  
  resultList <- list(
    annoLabels = annoLabels,
    annoColors = annoColors
  )
}
calcPlotAnnoLegend <- function(annoLabels, annoColors){
  if(is.null(annoLabels)){
    annoLabels <- vector(mode = "character")
    annoColors <- vector(mode = "character")
  }
  
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
calcPlotAnnoLegendForImage <- function(annoLabels, annoColors, maximumNumberOfLines=10){
  ## get and reorder annotations
  resultObj <- reorderAnnotationsForLegend(annoLabels, annoColors)
  annoLabels <- resultObj$annoLabels
  annoColors <- resultObj$annoColors
  
  calcPlotLegendForImage(annoLabels, annoColors, "Annotations", maximumNumberOfLines)
}
calcPlotScoresGroupsLegendForImage <- function(groups, colors, maximumNumberOfLines=10){
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
  yPositions <- yPositions[seq_len(numberOfLines)]
  
  symbolXPositions <- rep(x = xSpacing * 0.75, times = length(annoLabels))
  symbolYPositions <- yPositions[2:length(yPositions)]
  symbolYPositions <- symbolYPositions[seq_len(length(symbolXPositions))]
  
  if(numberOfLines > maximumNumberOfLines){
    labels <- labels[seq_len(maximumNumberOfLines)]
    annoColors <- annoColors[seq_len((maximumNumberOfLines - 1))]
    xPositions <- xPositions[seq_len(maximumNumberOfLines)]
    yPositions <- yPositions[seq_len(maximumNumberOfLines)]
    
    labels[length(labels)] <- paste("...", (numberOfLines - maximumNumberOfLines + 1), " more", sep = "")
    annoColors <- annoColors[seq_len(length(annoColors) - 1)]
    symbolXPositions <- symbolXPositions[seq_len(length(symbolXPositions) - 1)]
    symbolYPositions <- symbolYPositions[seq_len(length(symbolYPositions) - 1)]
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
  for(i in seq_len(length(annoColors)))
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
  yPositions <- seq(from = 0.9, to = -0.1, by = -1/length(labels))[seq_len(length(xPositions))]
  
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
    "Selection by HCA/PCA", 
    "Selection by fragment", 
    "Selection by search"
  )
  selectionColors <- c("blue", "green", "red")
  labels <- c("Selection marks", stickLabels)
  xPositions <- c(-0.05, rep(x = xSpacing, times = length(stickLabels)))
  #yPositions <- seq(from = 1, to = 0, by = 1/(-length(xPositions)))
  yPositions <- seq(from = 0.9, to = -0.1, by = -1/length(labels))[seq_len(length(xPositions))]
  
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
calcPlotDiscriminativityLegend <- function(){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  xSpacing <- 0.1
  stickLabels <- c(
    "0%", 
    "25%", 
    "50%", 
    "75%", 
    "100%"
  )
  labels <- c("Cluster-discriminating power", stickLabels)
  xPositions <- c(-0.05, rep(x = xSpacing, times = length(stickLabels)))
  #yPositions <- seq(from = 1, to = 0, by = 1/(-length(xPositions)))
  yPositions <- seq(from = 0.9, to = -0.1, by = -1/length(labels))[seq_len(length(xPositions))]
  
  symbolXPositions <- rep(x = xSpacing * 0.75, times = length(stickLabels))
  symbolYPositions <- yPositions[2:length(yPositions)]
  
  ##################################################################################################
  ## plot
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  ## points
  pointSizes      <- ms2StickPointSizeInitialSmall + c(0, 0.25, 0.5, 0.75, 1) * ms2StickPointSizeMaximumMultiplier
  
  #dendrogramClusterPointSizeMaximumMultiplier <- 0.75
  points(x = symbolXPositions, y = symbolYPositions, col = "black", pch=19, cex=pointSizes)
  
  ## labels
  graphics::text(x = xPositions, y = yPositions, labels = labels, pos = 4)
  #graphics::text(x = xPositions[[1]], y = yPositions[[1]], labels = labels[[1]], pos = 4)
  #graphics::text(x = xPositions[2:length(xPositions)], y = yPositions[2:length(yPositions)], labels = labels[2:length(labels)], pos = 4, adj = 1)
}
calcPlotHeatmapLegend <- function(dataList){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  legend_imageAbs <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = dataList$logAbsMax,        to =  0,                         length.out = 100), map = dataList$colorMapAbsoluteData ), nrow=1)))
  legend_imageLFC <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = dataList$logFoldChangeMax, to = -dataList$logFoldChangeMax, length.out = 100), map = dataList$colorMapLogFoldChange), nrow=1)))
  
  maximumNumberOfLabelsAbs <- 5
  maximumNumberOfLabelsLFC <- 2.5
  minY <- 0.0
  maxY <- 0.875
  epsilon <- 0.075
  middle <- ((maxY - minY) / 2)
  absMaxY <- middle - epsilon
  lfcMinY <- middle + epsilon
  
  ##############################################
  ## abs
  returnObj <- createTickLabels(maximumNumberOfLabelsAbs, dataList$logAbsMax, "10^")
  absLegendPositions <- returnObj$legendPositions
  absLegendLabels    <- returnObj$legendLabels
  
  absLegendPositions <- minY + absLegendPositions * (absMaxY - minY)
  
  ##############################################
  ## lfc
  returnObj <- createTickLabels(maximumNumberOfLabelsLFC, dataList$logFoldChangeMax, "")
  lfcLegendPositions <- returnObj$legendPositions
  lfcLegendLabels    <- returnObj$legendLabels
  
  ## double
  if(lfcLegendLabels[[1]] == "0"){
    ## do not copy zero
    lfcLegendPositions <- c(rev(lfcLegendPositions), -lfcLegendPositions[2:length(lfcLegendPositions)])
    revLabels <- lfcLegendLabels[2:length(lfcLegendLabels)]
    revLabels[revLabels != ""] <- paste("-", revLabels[revLabels != ""])
    lfcLegendLabels    <- c(rev(lfcLegendLabels   ), revLabels)
  } else {
    lfcLegendPositions <- c(rev(lfcLegendPositions), -lfcLegendPositions)
    revLabels <- lfcLegendLabels
    revLabels[revLabels != ""] <- paste("-", revLabels[revLabels != ""])
    lfcLegendLabels    <- c(rev(lfcLegendLabels   ), revLabels)
  }
  
  lfcLegendPositions <- lfcLegendPositions - min(lfcLegendPositions)
  lfcLegendPositions <- lfcLegendPositions / max(lfcLegendPositions)
  lfcLegendPositions <- lfcMinY + lfcLegendPositions * (maxY - lfcMinY)
  
  ##################################################################################################
  ## plot
  plot.new()
  plot.window(xlim = c(0, 3), ylim = c(0, 1))
  
  ## abs legend
  rasterImage(image = legend_imageAbs, xleft = 0, ybottom = minY, xright = 1, ytop = absMaxY)
  graphics::text(x = 1.2, y = absLegendPositions, labels = absLegendLabels, pos = 4)
  
  ## axis marks
  segments(
    x0  = rep(x = 0, times = length(absLegendPositions)),
    x1  = rep(x = 1.2, times = length(absLegendPositions)),
    y0  = absLegendPositions,
    y1  = absLegendPositions
  )
  ## frame
  segments(## lower hori; upper hori; left vert; right vert
    x0  = c(0, 0, 0, 1),
    x1  = c(1, 1, 0, 1),
    y0  = c(minY, absMaxY, minY   , minY   ),
    y1  = c(minY, absMaxY, absMaxY, absMaxY)
  )
  graphics::text(x = -0.2, y = absMaxY + 0.065, labels = "MS\u00B9 abundance", pos = 4)
  
  ## lfc legend
  rasterImage(image = legend_imageLFC, xleft = 0, ybottom = lfcMinY, xright = 1, ytop = maxY)
  graphics::text(x = 1.2, y = lfcLegendPositions, labels = lfcLegendLabels, pos = 4)
  #graphics::text(x = 2, y = seq(lfcMinY,maxY,l=5), labels = lfcLegendLabels)
  ## frame
  segments(## lower hori; upper hori; left vert; right vert
    x0  = c(0, 0, 0, 1),
    x1  = c(1, 1, 0, 1),
    y0  = c(lfcMinY, maxY, lfcMinY, lfcMinY),
    y1  = c(lfcMinY, maxY, maxY   , maxY   )
  )
  ## axis marks
  segments(
    x0  = rep(x = 0, times = length(lfcLegendPositions)),
    x1  = rep(x = 1.2, times = length(lfcLegendPositions)),
    y0  = lfcLegendPositions,
    y1  = lfcLegendPositions
  )
  graphics::text(x = -0.2, y = maxY + 0.09, labels = "log2(MS\u00B9 fold change)", pos = 4)
}
createTickLabels <- function(maximumNumberOfLabels, max, labelPrefix){
  if(max < 1)
    max <- 1
  maxInteger <- as.integer(max)
  numbers <- 0:maxInteger
  
  maxYlabelPosition <- maxInteger / max
  if(is.na(maxYlabelPosition))  ## data range is zero
    maxYlabelPosition <- 1
  legendPositions <- seq(from = 0, to = maxYlabelPosition, length.out = length(numbers))
  
  # if(length(legendPositions) > maximumNumberOfTicks){
  #   indeces <- rev(seq(from = length(numbers), to = 1, by = -ceiling( length(numbers) / maximumNumberOfTicks )))
  #   
  #   legendPositions <- legendPositions[indeces]
  #   numbers         <- numbers        [indeces]
  # }
  
  legendLabels    = paste(labelPrefix, numbers, sep = "")
  if(length(numbers) > maximumNumberOfLabels){
    indeces <- rev(seq(from = length(numbers), to = 1, by = -ceiling( length(numbers) / maximumNumberOfLabels )))
    
    legendLabels    <- legendLabels[indeces]
    legendPositions <- legendPositions[indeces]
  }
  
  returnObj <- list()
  returnObj$legendPositions <- legendPositions
  returnObj$legendLabels    <- legendLabels
  
  return(returnObj)
}
calcPlotMS2 <- function(dataList, fragmentsX = NULL, fragmentsY = NULL, fragmentsColor = NULL, fragmentsDiscriminativity = NULL, fragmentsX_02 = NULL, fragmentsY_02 = NULL, fragmentsColor_02 = NULL, xInterval = NULL, selectedFragmentIndex = NULL, dendrogramFragmentStatistics = FALSE){
  
  if(FALSE){
    dataList_ <<- dataList
    fragmentsX_ <<- fragmentsX
    fragmentsY_ <<- fragmentsY
    fragmentsColor_ <<- fragmentsColor
    fragmentsDiscriminativity_ <<- fragmentsDiscriminativity
    fragmentsX_02_ <<- fragmentsX_02
    fragmentsY_02_ <<- fragmentsY_02
    fragmentsColor_02_ <<- fragmentsColor_02
    xInterval_ <<- xInterval
    selectedFragmentIndex_ <<- selectedFragmentIndex
  }
  if(FALSE){
    dataList <<- dataList_
    fragmentsX <<- fragmentsX_
    fragmentsY <<- fragmentsY_
    fragmentsColor <<- fragmentsColor_
    fragmentsDiscriminativity <<- fragmentsDiscriminativity_
    fragmentsX_02 <<- fragmentsX_02_
    fragmentsY_02 <<- fragmentsY_02_
    fragmentsColor_02 <<- fragmentsColor_02_
    xInterval <<- xInterval_
    selectedFragmentIndex <<- selectedFragmentIndex_
  }
  
  
  ####################
  ## fragment spectrum
  if(is.null(xInterval))
    xInterval <- c(dataList$minimumMass, dataList$maximumMass)
  
  ## abundances greater one
  fragmentsY[fragmentsY > 1] <- 1
  
  ## y-axis
  if(is.null(fragmentsX_02)){
    ## nothing hovered; maybe something selected
    yInterval <- c(0, 1)
    #nodeColors <- rep(x = "black", times = length(fragmentsX))
    nodeColors <- fragmentsColor
    dataX <- fragmentsX
    dataY <- fragmentsY
    dataX2 <- fragmentsX
    dataY2 <- fragmentsY
    yTickPositions <- c(0, 0.25, 0.5, 0.75, 1)
    yTickLabels <- c(0, "", 0.5, "", 1)
  } else {
    ## something hovered; maybe something selected
    #nodeColors <- rep(x = "black", times = length(fragmentsX) + length(fragmentsX_02))
    nodeColors <- c(fragmentsColor, fragmentsColor_02)
    
    if(is.null(fragmentsX)){
      ## something hovered; nothing selected
      yInterval <- c(0, 1)
      dataX <- fragmentsX_02
      dataY <- fragmentsY_02
      dataX2 <- NULL
      dataY2 <- NULL
      yTickPositions <- c(0, 0.25, 0.5, 0.75, 1)
      yTickLabels <- c(0, "", 0.5, "", 1)
    } else {
      ## something hovered; something selected
      yInterval <- c(-1, 1)
      dataX <- c(fragmentsX, fragmentsX_02)
      dataY <- c(fragmentsY, -fragmentsY_02)
      dataX2 <- fragmentsX
      dataY2 <- fragmentsY
      yTickPositions <- c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
      yTickLabels <- c(1, "", 0.5, "", 0, "", 0.5, "", 1)
    }
  }
  
  ## node selection
  if(!is.null(selectedFragmentIndex))
    nodeColors[[selectedFragmentIndex]] <- "green"
  
  par(mar=c(6,4,3,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  plot(x = dataX, y = dataX, ylab = "Relative abundance", xlab = "m/z", xlim = xInterval, ylim = yInterval, xaxt='n', yaxt='n', col = nodeColors)
  axis(side = 2, at = yTickPositions, labels = yTickLabels)
  axis(side = 3)
  title("Fragment plot", line = 2)
  #mtext(side = 3, "m/z", line = 2)
  
  ## x-axis line
  if(!is.null(fragmentsX) & !is.null(fragmentsX_02)){
    xIntervalSize <- xInterval[[2]] - xInterval[[1]]
    xl <- xInterval[[1]] - xIntervalSize
    xr <- xInterval[[2]] + xIntervalSize
    segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  }
  
  ## label axis, fragment sticks with clickable points
  tickPositions <- dataX
  if(length(dataX) > 0){
    ## axis with the individual fragment m/z's (ticks, labels)
    axis(side = 1, at = dataX, labels = FALSE, las = 2)
    axis(side = 1, at = tickPositions, labels = format(x = dataX, digits = 0, nsmall = 4), las = 2, tick = FALSE)
    
    ## sticks
    points(x = dataX, y = dataY, col = nodeColors, type = "h", lwd=4)
    
    if(!is.null(dataX2)){
      ## clickable points
      pointSizes <- rep(x = ms2StickPointSizeInitial, times = length(dataX2))
      pointSizesSmall <- rep(x = ms2StickPointSizeInitialSmall, times = length(dataX2))
      pointColors <- rep(x = "black", times = length(dataX2))
      pointColorsSmall <- rep(x = "gray", times = length(dataX2))
      if(!is.null(selectedFragmentIndex)){
        pointSizes[[selectedFragmentIndex]] <- ms2StickPointSizeEmph
        pointSizesSmall[[selectedFragmentIndex]] <- ms2StickPointSizeEmphSmall
        pointColors[[selectedFragmentIndex]] <- "green"
        #pointColorsSmall[[selectedFragmentIndex]] <- "green"
      }
      pointSizeMultiplier <- fragmentsDiscriminativity * ms2StickPointSizeMaximumMultiplier
      pointSizes      <- pointSizes      + pointSizeMultiplier
      pointSizesSmall <- pointSizesSmall + pointSizeMultiplier
      
      points(x = dataX2, y = dataY2, col = pointColors,      pch=19, cex=pointSizes)
      points(x = dataX2, y = dataY2, col = pointColorsSmall, pch=19, cex=pointSizesSmall)
    }
  }
  
  if(!is.null(fragmentsX) & !is.null(fragmentsX_02)){
    graphics::text(labels = "Fragments from selection", x = xInterval[[2]], y = 0.9, pos = 2, adj = c(0,0))
    graphics::text(labels = "Fragments from mouse hover", x = xInterval[[2]], y = -0.9, pos = 2, adj = c(0,0))
  }
  if(!is.null(fragmentsX) & is.null(fragmentsX_02)){
    if(!dendrogramFragmentStatistics)
      graphics::text(labels = "Fragments from selection", x = xInterval[[2]], y = 0.95, pos = 2, adj = c(0,0))
    else
      graphics::text(labels = "Dendrogram statistics",    x = xInterval[[2]], y = 0.95, pos = 2, adj = c(0,0))
  }
  if(is.null(fragmentsX) & !is.null(fragmentsX_02)){
    graphics::text(labels = "Fragments from mouse hover", x = xInterval[[2]], y = 0.95, pos = 2, adj = c(0,0))
  }
}
getPcaPerformanceIndicator <- function(pcaObj, isScores){
  variableName <- ifelse(test = isScores, yes = "t", no = "p")
  variableName_dim <- paste(variableName, "_", pcaDimensionOne, sep = "")
  
  if("accurracyContribution" %in% names(pcaObj)){
    ## R2 and Q2 or accurracyContribution
    accOne <- format(x = pcaObj$accurracyContribution[[pcaDimensionOne]]*100, digits = 3)
    accTwo <- format(x = pcaObj$accurracyContribution[[pcaDimensionTwo]]*100, digits = 3)
    
    xAxisLabel  <- paste(variableName_dim, " (acc = ", accOne, "%)", sep = "")
    yAxisLabel  <- paste(variableName_dim, " (acc = ", accTwo, "%)", sep = "")
  } else {
    if(all(c("R2", "Q2") %in% names(pcaObj))){
      ## R2 and Q2
      r2One <- format(x = pcaObj$R2[[pcaDimensionOne]]*100, digits = 3)
      r2Two <- format(x = pcaObj$R2[[pcaDimensionTwo]]*100, digits = 3)
      q2One <- format(x = pcaObj$Q2[[pcaDimensionOne]]*100, digits = 3)
      q2Two <- format(x = pcaObj$Q2[[pcaDimensionTwo]]*100, digits = 3)
      
      xAxisLabel  <- paste(variableName_dim, " (R^2 = ", r2One, "%; Q^2 = ", q2One, "%)", sep = "")
      yAxisLabel  <- paste(variableName_dim, " (R^2 = ", r2Two, "%; Q^2 = ", q2Two, "%)", sep = "")
    } else {
      if("variance" %in% names(pcaObj)){
        ## explained variance
        varianceOne <- format(x = pcaObj$variance[[pcaDimensionOne]]*100, digits = 3)
        varianceTwo <- format(x = pcaObj$variance[[pcaDimensionTwo]]*100, digits = 3)
        xAxisLabel  <- paste(variableName_dim, " (var = ", varianceOne, "%)", sep = "")
        yAxisLabel  <- paste(variableName_dim, " (var = ", varianceTwo, "%)", sep = "")
      } else {
        ## no performace indicator there
        xAxisLabel  <- paste(variableName_dim, sep = "")
        yAxisLabel  <- paste(variableName_dim, sep = "")
      }
    }
  }
  
  resultObj <- list()
  resultObj$xAxisLabel <- xAxisLabel
  resultObj$yAxisLabel <- yAxisLabel
  return(resultObj)
}
calcPlotPCAscores <- function(pcaObj, dataList, filterObj, pcaDimensionOne, pcaDimensionTwo, showScoresLabels, xInterval = NULL, yInterval = NULL){
  palette <- colorPaletteScores()
  
  if(filterObj$filterBySamples){
    colorsForReplicates <- palette[unlist(lapply(X = filterObj$groups, FUN = function(x){ 
      groupIdx <- dataList$groupIdxFromGroupName(x)
      samples <- dataList$dataColumnsNameFunctionFromGroupName(x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
      samples <- intersect(samples, filterObj$sampleSet)
      rep(x = groupIdx, times = length(samples))
    }))]
  } else {
    colorsForReplicates <- palette[unlist(lapply(X = filterObj$groups, FUN = function(x){ 
      groupIdx <- dataList$groupIdxFromGroupName(x)
      rep(x = groupIdx, times = length(dataList$dataColumnsNameFunctionFromGroupName(x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))))
    }))]
  }
  
  ## data
  dataDimOne <- pcaObj$scores[, pcaDimensionOne]
  dataDimTwo <- pcaObj$scores[, pcaDimensionTwo]
  
  ## performance
  resultObj <- getPcaPerformanceIndicator(pcaObj = pcaObj, isScores = TRUE)
  xAxisLabel <- resultObj$xAxisLabel
  yAxisLabel <- resultObj$yAxisLabel
  
  ## xlim / ylim
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
  plot(
    x = dataDimOne, y = dataDimTwo, 
    xlim = xInterval, ylim = yInterval, 
    xlab = xAxisLabel, ylab = yAxisLabel, main = "Scores", 
    col = colorsForReplicates, pch=19, cex=1.
  )
  
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
    
    if(filterObj$filterBySamples){
      labels <- dataList$dataColumnsNameFunctionFromGroupNames(filterObj$groups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
      labels <- intersect(labels, filterObj$sampleSet)
    } else {
      labels <- dataList$dataColumnsNameFunctionFromGroupNames(filterObj$groups, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
    }
    graphics::text(x = dataDimOne, y = dataDimTwo, labels = labels, pos = 4)
  }
}
calcPlotPCAloadings <- function(
  pcaObj, dataList, filter, 
  pcaDimensionOne, pcaDimensionTwo, 
  selectionFragmentPcaLoadingSet = NULL, selectionAnalysisPcaLoadingSet = NULL, selectionSearchPcaLoadingSet = NULL, 
  xInterval = NULL, yInterval = NULL, 
  loadingsLabels = "None", showLoadingsAbundance = FALSE,
  showLoadingsFeaturesAnnotated = TRUE, showLoadingsFeaturesUnannotated = TRUE, showLoadingsFeaturesSelected = TRUE, showLoadingsFeaturesUnselected = TRUE
){
  if(FALSE){
    pcaObj_ <<- pcaObj
    dataList_ <<- dataList
    filter__ <<- filter
    pcaDimensionOne_ <<- pcaDimensionOne
    pcaDimensionTwo_ <<- pcaDimensionTwo
    selectionFragmentPcaLoadingSet_ <<- selectionFragmentPcaLoadingSet
    selectionAnalysisPcaLoadingSet_ <<- selectionAnalysisPcaLoadingSet
    selectionSearchPcaLoadingSet_ <<- selectionSearchPcaLoadingSet
    xInterval_ <<- xInterval
    yInterval_ <<- yInterval
    loadingsLabels_ <<- loadingsLabels
    showLoadingsAbundance_ <<- showLoadingsAbundance
    showLoadingsFeaturesAnnotated_ <<- showLoadingsFeaturesAnnotated
    showLoadingsFeaturesUnannotated_ <<- showLoadingsFeaturesUnannotated
    showLoadingsFeaturesSelected_ <<- showLoadingsFeaturesSelected
    showLoadingsFeaturesUnselected_ <<- showLoadingsFeaturesUnselected
  }
  if(FALSE){
    pcaObj <- pcaObj_
    dataList <- dataList_
    filter <- filter__
    pcaDimensionOne <- pcaDimensionOne_
    pcaDimensionTwo <- pcaDimensionTwo_
    selectionFragmentPcaLoadingSet <- selectionFragmentPcaLoadingSet_
    selectionAnalysisPcaLoadingSet <- selectionAnalysisPcaLoadingSet_
    selectionSearchPcaLoadingSet <- selectionSearchPcaLoadingSet_
    xInterval <- xInterval_
    yInterval <- yInterval_
    loadingsLabels <- loadingsLabels_
    showLoadingsAbundance <- showLoadingsAbundance_
    showLoadingsFeaturesAnnotated <- showLoadingsFeaturesAnnotated_
    showLoadingsFeaturesUnannotated <- showLoadingsFeaturesUnannotated_
    showLoadingsFeaturesSelected <- showLoadingsFeaturesSelected_
    showLoadingsFeaturesUnselected <- showLoadingsFeaturesUnselected_
  }
  
  resultObjAnno <- getPrecursorColors(dataList, filter)
  
  ## shown loadings features
  allFeatures <- seq_len(dataList$numberOfPrecursors)
  annotatedFeatures <- which(resultObjAnno$setOfColors != "black")
  selectedLoadingsFeatures <- union(union(selectionAnalysisPcaLoadingSet, selectionFragmentPcaLoadingSet), selectionSearchPcaLoadingSet)
  filter2 <- NULL
  if(showLoadingsFeaturesAnnotated)
    ## annotated features
    filter2 <- c(filter2, annotatedFeatures)
  if(showLoadingsFeaturesUnannotated)
    ## unannotated features
    filter2 <- c(filter2, setdiff(allFeatures, annotatedFeatures))
  if(showLoadingsFeaturesSelected)
    ## selected features
    filter2 <- c(filter2, selectedLoadingsFeatures)
  if(showLoadingsFeaturesUnselected)
    ## unselected features
    filter2 <- c(filter2, setdiff(allFeatures, selectedLoadingsFeatures))
  
  filter <- intersect(filter, unique(filter2))
  
  resultObjAnno <- getPrecursorColors(dataList, filter)
  #resultObjAnno$setOfAnnotations <- resultObjAnno$setOfAnnotations[filter]
  #resultObjAnno$setOfColors      <- resultObjAnno$setOfColors[filter]
  
  ## data
  dataDimOne <- pcaObj$loadings[, pcaDimensionOne]
  dataDimTwo <- pcaObj$loadings[, pcaDimensionTwo]
  
  ## performance
  resultObj <- getPcaPerformanceIndicator(pcaObj = pcaObj, isScores = TRUE)
  xAxisLabel <- resultObj$xAxisLabel
  yAxisLabel <- resultObj$yAxisLabel
  
  ## xlim / ylim
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
  
  dataDimOne <- dataDimOne[filter]
  dataDimTwo <- dataDimTwo[filter]
  selectionFragmentPcaLoadingSet <- intersect(selectionFragmentPcaLoadingSet, filter)
  selectionAnalysisPcaLoadingSet <- intersect(selectionAnalysisPcaLoadingSet, filter)
  selectionSearchPcaLoadingSet   <- intersect(selectionSearchPcaLoadingSet  , filter)
  
  numberOfPrecursors <- length(dataDimOne)
  poisX <- dataDimOne
  poisY <- dataDimTwo
  
  pointSizesAnno  <- rep(x = clusterNodePointSize0, times = numberOfPrecursors)
  pointColorsAnno <- resultObjAnno$setOfColors
  
  pointsAnalysis <- vector(mode = "logical", length = numberOfPrecursors)
  pointsAnalysis[match(x = selectionAnalysisPcaLoadingSet, table = filter)] <- TRUE
  pointsFragment <- vector(mode = "logical", length = numberOfPrecursors)
  pointsFragment[match(x = selectionFragmentPcaLoadingSet, table = filter)] <- TRUE
  pointsSearch <- vector(mode = "logical", length = numberOfPrecursors)
  pointsSearch[match(x = selectionSearchPcaLoadingSet, table = filter)] <- TRUE
  
  annotatedPoints <- pointColorsAnno != "black"
  selectedPoints  <- pointsAnalysis | pointsFragment | pointsSearch
  lv1points <-   annotatedPoints  &   selectedPoints
  lv2points <- (!annotatedPoints) &   selectedPoints
  lv3points <-   annotatedPoints  & (!selectedPoints)
  lv4points <- (!annotatedPoints) & (!selectedPoints)
  
  #poisX <- c(poisX[!annotatedPoints], poisX[annotatedPoints])
  #poisY <- c(poisY[!annotatedPoints], poisY[annotatedPoints])
  #pointSizesAnno <- c(pointSizesAnno[!annotatedPoints], pointSizesAnno[annotatedPoints])
  #pointColorsAnno <- c(pointColorsAnno[!annotatedPoints], pointColorsAnno[annotatedPoints])
  #pointsAnalysis <- c(pointsAnalysis[!annotatedPoints], pointsAnalysis[annotatedPoints])
  #pointsFragment <- c(pointsFragment[!annotatedPoints], pointsFragment[annotatedPoints])
  #pointsSearch <- c(pointsSearch[!annotatedPoints], pointsSearch[annotatedPoints])
  
  poisX           <- c(poisX          [lv4points], poisX          [lv3points], poisX          [lv2points], poisX          [lv1points])
  poisY           <- c(poisY          [lv4points], poisY          [lv3points], poisY          [lv2points], poisY          [lv1points])
  pointSizesAnno  <- c(pointSizesAnno [lv4points], pointSizesAnno [lv3points], pointSizesAnno [lv2points], pointSizesAnno [lv1points])
  pointColorsAnno <- c(pointColorsAnno[lv4points], pointColorsAnno[lv3points], pointColorsAnno[lv2points], pointColorsAnno[lv1points])
  pointsAnalysis  <- c(pointsAnalysis [lv4points], pointsAnalysis [lv3points], pointsAnalysis [lv2points], pointsAnalysis [lv1points])
  pointsFragment  <- c(pointsFragment [lv4points], pointsFragment [lv3points], pointsFragment [lv2points], pointsFragment [lv1points])
  pointsSearch    <- c(pointsSearch   [lv4points], pointsSearch   [lv3points], pointsSearch   [lv2points], pointsSearch   [lv1points])
  
  resultObjPoints <- generatePoints(
    poisX = poisX, poisY = poisY, 
    pointSizesAnno = pointSizesAnno, pointColorsAnno = pointColorsAnno, 
    pointsAnalysis = pointsAnalysis, pointsFragment = pointsFragment, pointsSearch = pointsSearch,
    pointSizeModifier = NULL
  )
  pointSizes    <- resultObjPoints$pointSizes
  pointColors   <- resultObjPoints$pointColors
  poisXpoints   <- resultObjPoints$poisXpoints
  poisYpoints   <- resultObjPoints$poisYpoints
  mappingToData <- resultObjPoints$mappingToData
  
  switch(loadingsLabels,
         "None"={## no labels
           labels <- NULL
         },
         "m/z / RT"={## mz/rt
           labels <- dataList$precursorLabels[filter]
           labels <- c(labels[lv4points], labels[lv3points], labels[lv2points], labels[lv1points])
           #labels <- c(labels[!annotatedPoints], labels[annotatedPoints])
         },
         "Metabolite name"={## name
           labels <- dataList$dataFrameInfos[filter, "Metabolite name"]
           labels <- c(labels[lv4points], labels[lv3points], labels[lv2points], labels[lv1points])
           #labels <- c(labels[!annotatedPoints], labels[annotatedPoints])
         },
         "Metabolite family"={## family
           featureFamilies <- dataList$annoArrayOfLists[filter]
           labels <- unlist(lapply(X = featureFamilies, FUN = function(x){
             ifelse(
               test = length(x) == 0, 
               yes = "-", 
               no = paste(unlist(x), collapse = ", ")
             )
           }))
           labels <- c(labels[lv4points], labels[lv3points], labels[lv2points], labels[lv1points])
         },
         {## unknown state
           stop(paste("Unknown loadingsLabels value", loadingsLabels))
         }
  )## end switch
  
  ## points
  #points(x = poisXpoints, y = poisYpoints, col = pointColors, pch=19, cex=pointSizes)
  if(showLoadingsAbundance){
    precursorMeansNorm <- dataList$dataFrameMeasurements[filter, "meanAllNormed"]
    precursorMeansNorm <- c(precursorMeansNorm[lv4points], precursorMeansNorm[lv3points], precursorMeansNorm[lv2points], precursorMeansNorm[lv1points])
    #precursorMeansNorm <- c(precursorMeansNorm[!annotatedPoints], precursorMeansNorm[annotatedPoints])
    precursorMeansNorm <- precursorMeansNorm[mappingToData]
    pointSizes <- pointSizes * 2 * precursorMeansNorm
  }
  
  ############################################################################################
  ## plot
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
  
  if(all(!is.null(labels), length(labels) > 0))
    graphics::text(  x = poisX - 0.0, y = poisY + 0.0, labels = labels, pos = 4)
  
  uniqueIndeces     <- which(!duplicated(resultObjAnno$setOfAnnotations))
  uniqueAnnotations <- resultObjAnno$setOfAnnotations[uniqueIndeces]
  uniqueColors      <- resultObjAnno$setOfColors[uniqueIndeces]
  
  resultList <- list(
    setOfAnnotations = uniqueAnnotations,
    setOfColors      = uniqueColors
  )
  return(resultList)
}
generatePoints <- function(poisX, poisY, pointSizesAnno, pointColorsAnno, pointsAnalysis, pointsFragment, pointsSearch, pointSizeModifier){
  numberOfPoisDrawn <- length(poisX)
  #pointIndeces <- seq_len(numberOfPoisDrawn)
  
  ## analysis
  pointSizesAnalysis  <- vector(mode = "numeric", length = numberOfPoisDrawn)
  pointColorsAnalysis <- vector(length = numberOfPoisDrawn)
  pointSizesAnalysis [pointsAnalysis] <- clusterNodePointSize1
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
  
  if(!is.null(pointSizeModifier)){
    pointSizesSearch   <- pointSizesSearch   + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier
    pointSizesFragment <- pointSizesFragment + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier
    pointSizesAnalysis <- pointSizesAnalysis + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier
    pointSizesAnno     <- pointSizesAnno     + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier
    #pointSizesSearch   <- pointSizesSearch   * (1 + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier)
    #pointSizesFragment <- pointSizesFragment * (1 + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier)
    #pointSizesAnalysis <- pointSizesAnalysis * (1 + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier)
    #pointSizesAnno     <- pointSizesAnno     * (1 + pointSizeModifier * dendrogramClusterPointSizeMaximumMultiplier)
  }
  
  pointSizes        <- c(pointSizesSearch[pointsSearch], pointSizesFragment[pointsFragment], pointSizesAnalysis[pointsAnalysis], pointSizesAnno)
  pointColors       <- c(pointColorsSearch[pointsSearch], pointColorsFragment[pointsFragment], pointColorsAnalysis[pointsAnalysis], pointColorsAnno)
  poisXpoints       <- c(poisX[pointsSearch], poisX[pointsFragment], poisX[pointsAnalysis], poisX)
  poisYpoints       <- c(poisY[pointsSearch], poisY[pointsFragment], poisY[pointsAnalysis], poisY)
  #pointIndecesPoints<- c(pointIndeces[pointsSearch], pointIndeces[pointsFragment], pointIndeces[pointsAnalysis], pointIndeces)
  
  mappingToData <- c(which(pointsSearch), which(pointsFragment), which(pointsAnalysis), seq_len(length(poisY)))
  
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
colorPaletteScores <- function(){
  palette <- colorPalette()
  palette <- c(
    palette[ 8: 1],
    palette[16: 9],
    palette[24:17],
    palette[32:25]
  )
  return(palette)
}
plotFragmentsFromDataList <- function(dataList, xInterval = NULL, yInterval = NULL, relative = FALSE){
  if(is.null(xInterval)){
    xMin <- min(dataList$fragmentMasses)
    xMax <- max(dataList$fragmentMasses)
    xInterval <- c(xMin, xMax)
  } else {
    xMin <- xInterval[[1]]
    xMax <- xInterval[[2]]
  }
  
  massIntervalSelection <- dataList$ms2_masses >= xMin & dataList$ms2_masses <= xMax
  #numberOfFragments <- dataList$ms2_numberOfFragments[massIntervalSelection]
  #masses            <- dataList$ms2_masses[massIntervalSelection]
  
  minimumNumberOfFragments <- 5
  selection <- massIntervalSelection & (dataList$ms2_numberOfFragments >= minimumNumberOfFragments)
  
  numberOfFragments <- dataList$ms2_numberOfFragments[selection]
  masses            <- dataList$ms2_masses[selection]
  
  colors <- cmap(x = numberOfFragments, map = dataList$colorMapFragmentData)
  
  #colors <- rep(x = "black", times = length(dataList$ms2_masses))
  colors <- cmap(x = numberOfFragments, map = dataList$colorMapFragmentData)
  
  if(relative)
    numberOfFragments <- numberOfFragments / dataList$numberOfPrecursors
  
  plotFragments(masses=masses, numberOfFragments=numberOfFragments, colors = colors, numberOfPrecursors=dataList$numberOfPrecursors, xInterval=xInterval, yInterval=yInterval)
}  
plotFragments <- function(masses, numberOfFragments, colors = NULL, numberOfPrecursors, title = NULL, xInterval = NULL, yInterval = NULL){
  if(is.null(xInterval)){
    xMin <- min(masses)
    xMax <- max(masses)
    xInterval <- c(xMin, xMax)
  } else {
    xMin <- xInterval[[1]]
    xMax <- xInterval[[2]]
  }
  
  yMin <- 0
  if(length(masses) == 0)
    yMax <- 1
  else
    yMax <- max(numberOfFragments)
  if(is.null(yInterval))
    yInterval <- c(yMin, yMax)
  
  if(is.null(colors))
    colors <- rep(x = "black", times = length(masses))
  
  #########################################################################################################
  ## plot
  par(mar=c(5,3,ifelse(test = is.null(title), yes = 2, no = 4),3), mgp = c(2, 1, 0))  ## c(bottom, left, top, right)
  #plot(x = dataList$ms2_masses, y = dataList$ms2_numberOfFragments, xlab = "Fragment mass", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = "Fragment plot", col = colors, pch=19, cex=1., xaxt='n')
  #plot(x = masses, y = numberOfFragments, xlab = "", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  #plot(x = NULL, y = NULL, xlab = "m/z", ylab = "Number of spectra", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  plot(x = NULL, y = NULL, xlab = "m/z", ylab = "Number of spectra", xlim = xInterval, ylim = yInterval, main = NULL, pch=19, cex=1., xaxt='n')
  if(!is.null(title))
    title(title, line = 3)
  axis(side = 3)
  
  ## axis
  dataOrder <- order(numberOfFragments)
  #axis(side = 1, at = dataList$ms2_masses, labels = dataList$ms2_masses, las = 2, tick = TRUE, col = colors)
  #for(i in 1:length(dataList$ms2_masses))
  for(i in dataOrder)
    axis(side = 1, at = masses[[i]], labels = format(x = masses[[i]], digits = 0, nsmall = 4), las = 2, tick = TRUE, col.axis = colors[[i]])
  
  points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], type = "h", lwd=4)
  #points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], pch=19, cex=1.)
  
  numberOfFragmentsInPercent <- numberOfFragments / numberOfPrecursors * 100
  yInterval <- c(0, yInterval[[2]] / numberOfPrecursors * 100)
  par(new = TRUE)
  plot(x = c(0, masses), y = c(0, numberOfFragmentsInPercent), xlim = xInterval, ylim = yInterval, axes=FALSE, type="n", xlab = "", ylab = "")
  axis(side = 4, line = NA, at = as.integer(pretty(c(0, numberOfFragmentsInPercent))))
  mtext(side = 4, line = 2, text = "Number of spectra in %")
  
  resultObj <- list()
  resultObj$poiFragmentX <- masses
  resultObj$poiFragmentY <- numberOfFragments
  
  return(resultObj)
}
plotFragments2 <- function(masses, numberOfFragments, numberOfPrecursors, xInterval = NULL, yInterval = NULL){
  if(FALSE){
    masses_ <<- masses
    numberOfFragments_ <<- numberOfFragments
    numberOfPrecursors_ <<- numberOfPrecursors
  }
  if(FALSE){
    masses <- masses_
    numberOfFragments <- numberOfFragments_
    numberOfPrecursors <- numberOfPrecursors_
  }
  
  
  if(is.null(xInterval)){
    xMin <- min(masses)
    xMax <- max(masses)
    xInterval <- c(xMin, xMax)
  } else {
    xMin <- xInterval[[1]]
    xMax <- xInterval[[2]]
  }
  
  massIntervalSelection <- masses >= xMin & masses <= xMax
  #numberOfFragments <- numberOfFragments[massIntervalSelection]
  #masses            <- masses[massIntervalSelection]
  
  minimumNumberOfFragments <- 1
  selection <- massIntervalSelection & (numberOfFragments >= minimumNumberOfFragments)
  
  numberOfFragments <- numberOfFragments[selection]
  masses            <- masses[selection]
  
  ms2PlotDataColorMapFragmentData  <- makecmap(
    x = c(0, max(numberOfFragments)), n = 100, 
    colFn = colorRampPalette(c('grey', 'black'))
  )
  #colors <- rep(x = "black", times = length(masses))
  colors <- cmap(x = numberOfFragments, map = ms2PlotDataColorMapFragmentData)
  
  yMin <- 0
  if(sum(selection) == 0)
    yMax <- 1
  else
    yMax <- max(numberOfFragments)
  if(is.null(yInterval))
    yInterval <- c(yMin, yMax)
  
  
  #########################################################################################################
  ## plot
  par(mar=c(5,3,2,3), mgp = c(2, 1, 0))  ## c(bottom, left, top, right)
  #plot(x = masses, y = numberOfFragments, xlab = "Fragment mass", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = "Fragment plot", col = colors, pch=19, cex=1., xaxt='n')
  #plot(x = masses, y = numberOfFragments, xlab = "", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  plot(x = NULL, y = NULL, xlab = "m/z", ylab = "Number of spectra", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  axis(side = 3)
  
  ## axis
  dataOrder <- order(numberOfFragments)
  #axis(side = 1, at = masses, labels = masses, las = 2, tick = TRUE, col = colors)
  #for(i in 1:length(masses))
  for(i in dataOrder)
    axis(side = 1, at = masses[[i]], labels = format(x = masses[[i]], digits = 0, nsmall = 4), las = 2, tick = TRUE, col.axis = colors[[i]])
  
  points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], type = "h", lwd=4)
  #points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], pch=19, cex=1.)
  
  ## 2nd y-axis TODO
  numberOfFragmentsInPercent <- numberOfFragments / numberOfPrecursors * 100
  par(new = TRUE)
  plot(x = c(0, masses), y = c(0, numberOfFragmentsInPercent), xlim = xInterval, axes=FALSE, type="n", xlab = "", ylab = "")
  axis(side = 4, line = NA, at = as.integer(pretty(c(0, numberOfFragmentsInPercent))))
  mtext(side = 4, line = 2, text = "Number of spectra in %")
  
  resultObj <- list()
  resultObj$poiFragmentX <- masses
  resultObj$poiFragmentY <- numberOfFragments
  
  return(resultObj)
}
calcPlotSpectrumVsClass_small_old <- function(dataX_spec, dataY_spec, dataX_class, dataY_class, xInterval){
  yInterval <- c(-1, 1)
  
  colors_spec <- rep(x = "grey", times = length(dataX_spec))
  colors_spec[dataX_spec %in% dataX_class] <- "black"
  
  colors_class <- rep(x = "black", times = length(dataY_spec))
  
  dataY_class <- -dataY_class
  
  par(mar=c(0,0,0,0))  ## c(bottom, left, top, right)
  plot(NA, ylab = "", xlab = "", xlim = xInterval, ylim = yInterval, xaxt='n', yaxt='n')
  
  segments(x0 = xInterval[[1]], x1 = xInterval[[2]], y0 = 0, y1 = 0, col = "black", lwd = 1)
  
  points(x = dataX_spec,  y = dataY_spec,  col = colors_spec,  type = "h", lwd=4)
  points(x = dataX_class, y = dataY_class, col = "black", type = "h", lwd=4)
  
}
calcPlotSpectrumVsClass_small <- function(masses_spec, intensity_spec, colors_spec, masses_class, frequency_class, colors_class, xInterval){
  yInterval <- c(-1, 1)
  
  intensity_spec[intensity_spec > 1] <- 1
  frequency_class <- -frequency_class
  
  ## plot
  par(mar=c(0,0,0,0))  ## c(bottom, left, top, right)
  plot(NA, ylab = "", xlab = "", xlim = xInterval, ylim = yInterval, xaxt='n', yaxt='n')
  
  ## x-axis line
  xIntervalSize <- xInterval[[2]] - xInterval[[1]]
  xl <- xInterval[[1]] - xIntervalSize
  xr <- xInterval[[2]] + xIntervalSize
  segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  
  ## sticks
  points(x = masses_spec,  y = intensity_spec,  col = colors_spec,  type = "h", lwd=4)
  points(x = masses_class, y = frequency_class, col = colors_class, type = "h", lwd=4)
}
calcPlotSpectrumVsClass_big <- function(masses_spec, intensity_spec, colors_spec, masses_class, frequency_class, colors_class, singleSpec, xInterval){
  if(TRUE){
    masses_spec_ <<- masses_spec
    intensity_spec_ <<- intensity_spec
    colors_spec_ <<- colors_spec
    masses_class_ <<- masses_class
    frequency_class_ <<- frequency_class
    colors_class_ <<- colors_class
    singleSpec_ <<- singleSpec
    xInterval_ <<- xInterval
  }
  if(FALSE){
    masses_spec <- masses_spec_
    intensity_spec <- intensity_spec_
    colors_spec <- colors_spec_
    masses_class <- masses_class_
    frequency_class <- frequency_class_
    colors_class <- colors_class_
    singleSpec <<- singleSpec_
    xInterval <- xInterval_
  }
  
  onlyClass <- is.null(masses_class)
  yInterval <- c(ifelse(test = is.null(masses_class), yes = 0, no = -1), 1)
  
  ## abundances greater one
  intensity_spec[intensity_spec > 1] <- 1
  
  if(!is.null(frequency_class))
    frequency_class <- -frequency_class
  
  yTickPositions <- c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
  yTickLabels <- c(1, "", 0.5, "", 0, "", 0.5, "", 1)
  
  pointSizes <- rep(x = ms2StickPointSizeInitial, times = length(masses_spec))
  pointSizesSmall <- rep(x = ms2StickPointSizeInitialSmall, times = length(masses_spec))
  pointColors <- rep(x = "black", times = length(masses_spec))
  pointColorsSmall <- rep(x = "gray", times = length(masses_spec))
  
  these <- which(colors_spec == "black")
  if(length(these) == 0)
    these <- length(masses_spec)+1
  
  plotTitle <-  ifelse(test = singleSpec, yes = "Spectrum versus class", no = 
                ifelse(test = onlyClass,  yes = "Metabolite family spectrum", no = "Metabolite family versus class"))
  yAxisLabel <- ifelse(test = singleSpec, yes = "Frequency / Intensity", no = 
                ifelse(test = onlyClass,  yes = "Frequency", no = "Frequency / Frequency"))
  
  ## plot
  par(mar=c(6,4.1,4,0.1), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  plot(NA, ylab = yAxisLabel, xlab = "m/z", xlim = xInterval, ylim = yInterval, xaxt='n', yaxt='n')
  axis(side = 2, at = yTickPositions, labels = yTickLabels)
  axis(side = 3)
  title(plotTitle, line = 2.5)
  #mtext(side = 3, "m/z", line = 2)
  
  ## x-axis line
  xIntervalSize <- xInterval[[2]] - xInterval[[1]]
  xl <- xInterval[[1]] - xIntervalSize
  xr <- xInterval[[2]] + xIntervalSize
  segments(x0 = xl, x1 = xr, y0 = 0, y1 = 0, col = "black", lwd = 1)
  
  ## axis with the individual fragment m/z's (ticks, labels)
  axis(side = 1, at = masses_spec, labels = FALSE, las = 2)
  axis(side = 1, at = masses_spec[-these], labels = format(x = masses_spec[-these], digits = 0, nsmall = 4), las = 2, tick = FALSE, col.axis = "grey")
  if(length(these) == 1){
    if(these != length(masses_spec)+1)
      axis(side = 1, at = masses_spec[ these], labels = format(x = masses_spec[ these], digits = 0, nsmall = 4), las = 2, tick = FALSE, col.axis = "black")
  } else {
    axis(side = 1, at = masses_spec[ these], labels = format(x = masses_spec[ these], digits = 0, nsmall = 4), las = 2, tick = FALSE, col.axis = "black")
  }
  
  ## sticks
  points(x = masses_spec[-these],  y = intensity_spec[-these],  col = colors_spec[-these],  type = "h", lwd=4)
  if(!is.null(masses_class)){
    these2 <- colors_class == "black"
    points(x = masses_class[!these2], y = frequency_class[!these2], col = colors_class[!these2], type = "h", lwd=4)
    points(x = masses_class[ these2], y = frequency_class[ these2], col = colors_class[ these2], type = "h", lwd=4)
  }
  
  ## points
  points(x = masses_spec, y = intensity_spec, col = pointColors,      pch=19, cex=pointSizes)
  points(x = masses_spec, y = intensity_spec, col = pointColorsSmall, pch=19, cex=pointSizesSmall)
  
  ## black ones
  points(x = masses_spec[these], y = intensity_spec[these], col = colors_spec[these], type = "h", lwd=4)
  points(x = masses_spec[these], y = intensity_spec[these], col = pointColors[these], pch=19, cex=pointSizes[these])
  points(x = masses_spec[these], y = intensity_spec[these], col = pointColorsSmall[these], pch=19, cex=pointSizesSmall[these])
  
  ## plot labels
  graphics::text(labels = ifelse(test = singleSpec, yes = "Fragments from spectrum", no = "Fragments from spectra"), x = xInterval[[2]], y = 0.9, pos = 2, adj = c(0,0))
  graphics::text(labels = "Fragments from class",    x = xInterval[[2]], y = -0.9, pos = 2, adj = c(0,0))
}
