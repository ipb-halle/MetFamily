
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
           precursorLabels[tooLong] <- paste(substring(text = precursorLabels[tooLong], first = 1, last = maximumNumberOfCharacters - 3), rep(x = "...", times = length(tooLong)), sep = "")
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
           precursorLabels[tooLong] <- paste(substring(text = precursorLabels[tooLong], first = 1, last = maximumNumberOfCharacters - 3), rep(x = "...", times = length(tooLong)), sep = "")
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
  par(mar=c(7.25,4,2,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
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
    for(selectionAnalysisTreeNode in selectionAnalysisTreeNodeSet){
      clusterMembers <- NULL
      if(selectionAnalysisTreeNode > 0){
        clusterMembers <- c(clusterMembers, clusterDataList$innerNodeMembersTreeLeaves[[selectionAnalysisTreeNode]])
      } else {
        clusterMembers <- c(clusterMembers, -selectionAnalysisTreeNode)
      }
    }
    
    colLab <- colorLabels(precursorLabelsWithIdx, clusterMembers, 'blue')
    dend <- dendrapply(dend, colLab)
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
calcPlotHeatmap <- function(dataList, filterObj, clusterDataList, selectedTreeNodeSet, frameColor, xInterval = NULL){
  if(is.null(xInterval))
    xInterval <- c(1, clusterDataList$numberOfPrecursorsFiltered)
  
  ####################
  ## heatmap
  groups <- rev(dataList$groups)
  columnsOfInterest <- unlist(lapply(X = groups, FUN = function(x){
    dataList$dataMeanColumnNameFunctionFromName(x)
  }))
  numberOfGroups <- length(columnsOfInterest)
  labels <- groups
  
  par(mar=c(0,4,0,0), mgp = c(3, 1, 0))  ## c(bottom, left, top, right) ## c(title, axis, label)
  
  
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
    
    #print("")
    #print(selectedTreeNodeSet)
    #print(clusterMemberList)
    #print(positionsList)
    #print(intervals)
    #print(intervalMatrix)
    #print("")
  }
  
  if(numberOfGroups == 0){
    print("### calcPlotHeatmap: no columnsOfInterest")
    plot.new()
    return()
  }
  
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
  
  colors <- lapply(X = columnsOfInterest[cluster$order], FUN = function(x){
    dataList$colorMatrixDataFrame[filterObj$filter, x][clusterDataList$cluster$order]
  })
  
  plot(x = c(1, clusterDataList$numberOfPrecursorsFiltered), y = c(0, numberOfGroups), type= "n", xlab = "", ylab = "", axes = FALSE, xlim = xInterval, ylim = c(0, numberOfGroups))
  rect(
    xleft   = rep(x = seq_len(clusterDataList$numberOfPrecursorsFiltered) - 0.5, times = numberOfGroups), 
    xright  = rep(x = seq_len(clusterDataList$numberOfPrecursorsFiltered) + 0.5, times = numberOfGroups), 
    ybottom = unlist(lapply(X = 0:(numberOfGroups - 1), FUN = function(x){rep(x = x, times = clusterDataList$numberOfPrecursorsFiltered)})),
    ytop    = unlist(lapply(X = seq_len(numberOfGroups),       FUN = function(x){rep(x = x, times = clusterDataList$numberOfPrecursorsFiltered)})),
    col     = unlist(colors),
    border = NA
  )
  if(!is.null(selectedTreeNodeSet))
    apply(X = intervalMatrix, MARGIN = 2, FUN = function(x){
      rect(xleft = x[[1]] - 0.5, xright = x[[2]] + 0.5, ybottom = 0, ytop = numberOfGroups, border = frameColor, lwd = 2)
    })
  
  axis(side = 2, at = seq(from = 0.5, by = 1, length.out = numberOfGroups), labels = labels, las = 2, tick = TRUE)
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
calcPlotMS2 <- function(dataList, fragmentsX = NULL, fragmentsY = NULL, fragmentsColor = NULL, fragmentsDiscriminativity = NULL, fragmentsX_02 = NULL, fragmentsY_02 = NULL, fragmentsColor_02 = NULL, fragmentsDiscriminativity_02 = NULL, xInterval = NULL, selectedFragmentIndex = NULL){
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
      pointSizeMultiplier <- c(fragmentsDiscriminativity, fragmentsDiscriminativity_02) * ms2StickPointSizeMaximumMultiplier
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
    graphics::text(labels = "Fragments from selection", x = xInterval[[2]], y = 0.95, pos = 2, adj = c(0,0))
  }
  if(is.null(fragmentsX) & !is.null(fragmentsX_02)){
    graphics::text(labels = "Fragments from mouse hover", x = xInterval[[2]], y = 0.95, pos = 2, adj = c(0,0))
  }
}
calcPlotPCAscores <- function(pcaObj, dataList, filterObj, pcaDimensionOne, pcaDimensionTwo, showScoresLabels, xInterval = NULL, yInterval = NULL){
  palette <- colorPaletteScores()
  colorsForReplicates <- palette[unlist(lapply(X = filterObj$groups, FUN = function(x){ 
    groupIdx <- dataList$groupIdxFromGroupName(x)
    rep(x = groupIdx, times = length(dataList$dataColumnsNameFunctionFromName(x)))
  }))]
  
  dataDimOne <- pcaObj$scores[, pcaDimensionOne]
  dataDimTwo <- pcaObj$scores[, pcaDimensionTwo]
  
  #varianceOne <- format(x = pcaObj$variance[[pcaDimensionOne]], digits = 3)
  #varianceTwo <- format(x = pcaObj$variance[[pcaDimensionTwo]], digits = 3)
  r2One <- format(x = pcaObj$R2[[pcaDimensionOne]], digits = 3)
  r2Two <- format(x = pcaObj$R2[[pcaDimensionTwo]], digits = 3)
  q2One <- format(x = pcaObj$Q2[[pcaDimensionOne]], digits = 3)
  q2Two <- format(x = pcaObj$Q2[[pcaDimensionTwo]], digits = 3)
  
  #xAxisLabel  <- paste("t_", pcaDimensionOne, " (", varianceOne, "%)", sep = "")
  #yAxisLabel  <- paste("t_", pcaDimensionTwo, " (", varianceTwo, "%)", sep = "")
  xAxisLabel  <- paste("t_", pcaDimensionOne, " (R^2 = ", r2One, "; Q^2 = ", q2One, ")", sep = "")
  yAxisLabel  <- paste("t_", pcaDimensionTwo, " (R^2 = ", r2Two, "; Q^2 = ", q2Two, ")", sep = "")
  
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
    labels <- dataList$dataColumnsNameFunctionFromNames(filterObj$groups)
    graphics::text(x = dataDimOne, y = dataDimTwo, labels = labels, pos = 4)
  }
}
calcPlotPCAloadings <- function(pcaObj, dataList, filter, pcaDimensionOne, pcaDimensionTwo, selectionFragmentPcaLoadingSet = NULL, selectionAnalysisPcaLoadingSet = NULL, selectionSearchPcaLoadingSet = NULL, xInterval = NULL, yInterval = NULL, loadingsLabels = "None", showLoadingsFeatures = "All", showLoadingsAbundance = FALSE){
  #varianceOne <- format(x = pcaObj$variance[[pcaDimensionOne]], digits = 3)
  #varianceTwo <- format(x = pcaObj$variance[[pcaDimensionTwo]], digits = 3)
  r2One <- format(x = pcaObj$R2[[pcaDimensionOne]], digits = 3)
  r2Two <- format(x = pcaObj$R2[[pcaDimensionTwo]], digits = 3)
  q2One <- format(x = pcaObj$Q2[[pcaDimensionOne]], digits = 3)
  q2Two <- format(x = pcaObj$Q2[[pcaDimensionTwo]], digits = 3)
  
  #xAxisLabel  <- paste("p_", pcaDimensionOne, " (", varianceOne, "%)", sep = "")
  #yAxisLabel  <- paste("p_", pcaDimensionTwo, " (", varianceTwo, "%)", sep = "")
  xAxisLabel  <- paste("p_", pcaDimensionOne, " (R^2 = ", r2One, "; Q^2 = ", q2One, ")", sep = "")
  yAxisLabel  <- paste("p_", pcaDimensionTwo, " (R^2 = ", r2Two, "; Q^2 = ", q2Two, ")", sep = "")
  
  filter2 <- NULL
  switch(showLoadingsFeatures,
         "All"={## all features
           filter2 <- seq_len(dataList$numberOfPrecursors)
         },
         "Only selected"={## selected features
           filter2 <- union(union(selectionAnalysisPcaLoadingSet, selectionFragmentPcaLoadingSet), selectionSearchPcaLoadingSet)
         },
         "Only unselected"={## unselected features
           filter2 <- setdiff(seq_len(dataList$numberOfPrecursors), union(union(selectionAnalysisPcaLoadingSet, selectionFragmentPcaLoadingSet), selectionSearchPcaLoadingSet))
         },
         {## unknown state
           stop(paste("Unknown showLoadingsFeatures value", showLoadingsFeatures))
         }
  )## end switch
  
  filter <- intersect(filter, filter2)
  
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
  
  dataDimOne <- dataDimOne[filter2]
  dataDimTwo <- dataDimTwo[filter2]
  selectionFragmentPcaLoadingSet <- selectionFragmentPcaLoadingSet[filter2]
  selectionAnalysisPcaLoadingSet <- selectionAnalysisPcaLoadingSet[filter2]
  selectionSearchPcaLoadingSet   <- selectionSearchPcaLoadingSet  [filter2]
  
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
  
  pointSizes  <- c(pointSizesSearch[pointsSearch], pointSizesFragment[pointsFragment], pointSizesAnalysis[pointsAnalysis], pointSizesAnno)
  pointColors <- c(pointColorsSearch[pointsSearch], pointColorsFragment[pointsFragment], pointColorsAnalysis[pointsAnalysis], pointColorsAnno)
  poisXpoints <- c(poisX[pointsSearch], poisX[pointsFragment], poisX[pointsAnalysis], poisX)
  poisYpoints <- c(poisY[pointsSearch], poisY[pointsFragment], poisY[pointsAnalysis], poisY)
  
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
plotFragments <- function(dataList, xInterval = NULL, yInterval = NULL, relative = FALSE){
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
  
  #colors <- rep(x = "black", times = length(dataList$masses))
  colors <- cmap(x = numberOfFragments, map = dataList$colorMapFragmentData)
  
  if(relative)
    numberOfFragments <- numberOfFragments / dataList$numberOfPrecursors
  
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
  #plot(x = dataList$masses, y = dataList$numberOfFragments, xlab = "Fragment mass", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = "Fragment plot", col = colors, pch=19, cex=1., xaxt='n')
  #plot(x = masses, y = numberOfFragments, xlab = "", ylab = "Precursors", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  plot(x = NULL, y = NULL, xlab = "m/z", ylab = "Number of spectra", xlim = xInterval, ylim = yInterval, main = NULL, col = colors, pch=19, cex=1., xaxt='n')
  axis(side = 3)
  
  ## axis
  dataOrder <- order(numberOfFragments)
  #axis(side = 1, at = dataList$masses, labels = dataList$masses, las = 2, tick = TRUE, col = colors)
  #for(i in 1:length(dataList$masses))
  for(i in dataOrder)
    axis(side = 1, at = masses[[i]], labels = format(x = masses[[i]], digits = 0, nsmall = 4), las = 2, tick = TRUE, col.axis = colors[[i]])
  
  points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], type = "h", lwd=4)
  #points(x = masses[dataOrder], y = numberOfFragments[dataOrder], col = colors[dataOrder], pch=19, cex=1.)
  
  ## 2nd y-axis TODO
  numberOfFragmentsInPercent <- numberOfFragments / dataList$numberOfPrecursors * 100
  par(new = TRUE)
  plot(x = c(0, masses), y = c(0, numberOfFragmentsInPercent), xlim = xInterval, axes=FALSE, type="n", xlab = "", ylab = "")
  axis(side = 4, line = NA, at = as.integer(pretty(c(0, numberOfFragmentsInPercent))))
  mtext(side = 4, line = 2, text = "Number of spectra in %")
  
  resultObj <- list()
  resultObj$poiFragmentX <- masses
  resultObj$poiFragmentY <- numberOfFragments
  
  return(resultObj)
}