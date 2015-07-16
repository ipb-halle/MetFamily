
#install.packages("squash")
library("squash")

#########################################################################################
## annotate and process matrix
readClusterData <- function(file, distance = "Jaccard", progress = FALSE){
  if(progress)  setProgress(value = 0, detail = "Parsing")
  dataFrame <- read.csv(file = file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  #### TODO remove
  #dataFrame <- dataFrame[c(1:3, 16:20), ]
  #dataFrame <- dataFrame[c(1:50), ]
  ####
  
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
  
  dataFrame <- dataFrame[4:nrow(dataFrame), ]
  idColumns <- which(tagsSector == "ID")
  #precursorLabels <- paste(dataFrame[, 1], dataFrame[, 2], sep = " / ")
  precursorLabels <- apply(X = dataFrame[, idColumns], MARGIN = 1, FUN = paste, collapse = " / ")
  rownames(dataFrame) <- precursorLabels
  numberOfPrecursors <- nrow(dataFrame)
  
  dataColumns <- which(!is.na(suppressWarnings(as.numeric(tagsSector))))
  groups <- unique(as.numeric(tagsSector[dataColumns]))
  groupsStartEnd <- list()
  for(groupIdx in 1:length(groups))
    groupsStartEnd[[groupIdx]] <- c(min(which(tagsSector == groups[[groupIdx]])), max(which(tagsSector == groups[[groupIdx]])))
  
  ####################
  ## measurement data to colors
  if(progress)  incProgress(amount = 0.1, detail = "Coloring")
  
  dataFrameMeasurements <- data.frame(matrix(nrow = nrow(dataFrame), ncol = 0))
  rownames(dataFrameMeasurements) <- rownames(dataFrame)
  
  ## mean of groups
  dataColumnNameFunction <- function(group){
    as.character(groups[[groupIdx]])
  }
  dataColumns <- list()
  for(groupIdx in 1:length(groups)){
    dataColumnName <- dataColumnNameFunction(groups[[groupIdx]])
    dataColumns[[groupIdx]] <- dataColumnName
    dataFrameMeasurements[, dataColumnName] <- apply(X = data.matrix(dataFrame[, groupsStartEnd[[groupIdx]][[1]]:groupsStartEnd[[groupIdx]][[2]]]), MARGIN = 1, FUN = mean)
  }
  dataColumns <- unlist(dataColumns)
  
  ## log fold change between groups
  lfcColumnNameFunction <- function(groupOne, groupTwo){
    paste("LFC", as.character(groupOne), "vs", as.character(groupTwo), sep = "_")
  }
  lfcColumnNames <- list()
  for(groupIdx1 in 1:length(groups))
    for(groupIdx2 in 1:length(groups)){
      lfcColumnName <- lfcColumnNameFunction(groupOne = groups[[groupIdx1]], groupTwo = groups[[groupIdx2]])
      lfcColumnNames[[length(lfcColumnNames) + 1]] <- lfcColumnName
      dataFrameMeasurements[, lfcColumnName] <- log(dataFrameMeasurements[, as.character(groups[[groupIdx1]])] / dataFrameMeasurements[, as.character(groups[[groupIdx2]])])
      
      ## tackle zero values
      dataFrameMeasurements[is.na(dataFrameMeasurements[, lfcColumnName]), lfcColumnName] <- 0
      dataFrameMeasurements[is.infinite(dataFrameMeasurements[, lfcColumnName]), lfcColumnName] <- 0
    }
  lfcColumnNames <- unlist(lfcColumnNames)
  
  ## map to colors
  matrixDataFrame <- data.matrix(dataFrameMeasurements)
  
  absMax <- max(matrixDataFrame[, dataColumns])
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
  
  ## translate and box colors
  colorDataFrame <- dataFrameMeasurements
  colorDataFrame[, dataColumns]    <- cmap(matrixDataFrame[, dataColumns   ], colorMapAbsoluteData)
  colorDataFrame[, lfcColumnNames] <- cmap(matrixDataFrame[, lfcColumnNames], colorMapLogFoldChange)
  colorMatrixDataFrame <- as.matrix(colorDataFrame)
  
  ##########################################
  ## compute distances
  
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
  featureMatrixBinary <- matrix(nrow = numberOfPrecursors, ncol = numberOfFragments)
  for(i in 1:numberOfPrecursors){
    if(progress)  incProgress(amount = 0.1 / numberOfPrecursors, detail = paste("Features binary ", i, " / ", numberOfPrecursors, sep = ""))
    featureVectorBinary <- featureMatrix[i, ] != 0
    featureMatrixBinary[i, ] <- featureVectorBinary
    featureIndeces[[i]] <- which(featureVectorBinary)
  }
  featureIndexMatrix <- matrix(nrow = numberOfPrecursors, ncol = max(sapply(X = featureIndeces, FUN = length)))
  for(i in 1:numberOfPrecursors)
    featureIndexMatrix[i, 1:length(featureIndeces[[i]])] <- featureIndeces[[i]]
  
  ## remove unused columns
  fragmentThere <- apply(X = featureMatrixBinary, MARGIN = 2, FUN = any)
  minimumMass <- min(fragmentMasses[fragmentThere])
  maximumMass <- max(fragmentMasses[fragmentThere])
  
  ## compute distance matrix: 1 - jaccard similarity
  if(progress)  incProgress(amount = 0.1, detail = "Distances")
  
  distanceMatrix <- NULL
  switch(distance,
         "Jaccard"={
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Jaccard (weighted)"={
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Similarity (weighted)"={
           similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Jaccard (weighted2)"={
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Similarity (weighted2)"={
           similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
         "Jaccard (weighted3)"={
           counter <<- 0
           distanceMatrix <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(x)
             {  
               counter <<- counter + 1
               if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", counter, " / ", numberOfPrecursors, sep = ""))
               unionSum        <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(y){  sum(fragmentFrequency[union(    x = x, y = y)], na.rm = TRUE)  })
               intersectionSum <- apply(X = featureIndexMatrix, MARGIN = 1, FUN = function(y){  sum(fragmentFrequency[intersect(x = x, y = y)], na.rm = TRUE)  })
               
               1 - intersectionSum / unionSum
             }
           )
         },
         "Manhatten"={
           ## Rasmussen 2008
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
           ## Gaquerel 2015: standard normalized dot product (NDP) / cosine correlation
           similarityMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in 1:numberOfPrecursors){
             if(progress)  incProgress(amount = 0.2 / numberOfPrecursors, detail = paste("Distances ", i, " / ", numberOfPrecursors, sep = ""))
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
  
  ##########################################
  ## box
  if(progress)  incProgress(amount = 0.1, detail = "Boxing")
  dataList <- list()
  ## data
  dataList$dataFrame <- dataFrame
  dataList$numberOfPrecursors <- numberOfPrecursors
  dataList$groups <- groups
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
  dataList$lfcColumnNameFunction <- lfcColumnNameFunction
  dataList$dataColumnNameFunction <- dataColumnNameFunction
  ## features
  dataList$distanceMatrix <- distanceMatrix
  dataList$featureMatrix <- featureMatrix
  dataList$featureMatrixBinary <- featureMatrixBinary
  dataList$featureIndeces <- featureIndeces
  dataList$featureIndexMatrix <- featureIndexMatrix
  
  if(progress)  setProgress(1)
  
  return(dataList)
}

filterData <- function(dataList, groupOne, groupTwo, filter_average, filter_lfc, filter_ms2_masses, filter_ms2_ppm, progress = FALSE){
  ##########################################
  ## filter
  filter <- rep(x = TRUE, times = dataList$numberOfPrecursors)
  
  ## filter_average
  if(!is.null(filter_average))
    filter <- filter & apply(dataList$dataFrameMeasurements[, c(dataList$dataColumnNameFunction(groupOne), dataList$dataColumnNameFunction(groupTwo))], 1, max) >= filter_average
  
  ## filter_lfc
  if(!is.null(filter_lfc)){
    if(filter_lfc != 0){
      if(filter_lfc > 0)
        filter <- filter & dataList$dataFrameMeasurements[, dataList$lfcColumnNameFunction(groupOne, groupTwo)] >= filter_lfc
      else
        filter <- filter & dataList$dataFrameMeasurements[, dataList$lfcColumnNameFunction(groupOne, groupTwo)] <= filter_lfc
    }
  }
  
  ## filter_ms2_masses, filter_ms2_ppm
  if(!is.null(filter_ms2_masses) & !is.null(filter_ms2_ppm) & length(filter_ms2_masses) > 0){
    error <- abs(dataList$fragmentMasses) * filter_ms2_ppm / 1E6
    #print(paste("ff", error))
    numberOfSearchedFragmentMasses <- length(filter_ms2_masses)
    fragmentColumns <- vector(mode = "numeric")
    
    for(fragmentIndex in 1:numberOfSearchedFragmentMasses){
      columns <- unlist(which(abs(dataList$fragmentMasses - filter_ms2_masses[[fragmentIndex]]) <= error))
      if(length(columns) > 0)
        fragmentColumns[[length(fragmentColumns) + 1]] <- unlist(which(abs(dataList$fragmentMasses - filter_ms2_masses[[fragmentIndex]]) <= error))
    }
    if(length(fragmentColumns) > 0){
      for(columns in fragmentColumns)
        filter <- filter & dataList$featureMatrix[, columns] != 0
    } else {
      filter <- rep(x = FALSE, times = dataList$numberOfPrecursors)
    }
  }
  
  filter <- which(filter)
  return (filter)
}
calculateCluster <- function(dataList, filter, groupOne, groupTwo, filter_average, filter_lfc, filter_ms2_masses, filter_ms2_ppm, progress = FALSE){
  numberOfPrecursorsFiltered <- length(filter)
  ##########################################
  ## compute gui stuff
  
  if(progress)  incProgress(amount = 0.5, detail = "Clustering")
  ## compute and annotate cluster
  cluster <- hclust(d = stats::as.dist(m = dataList$distanceMatrix[filter, filter]), method = "complete")
  cluster$labels <- dataList$precursorLabels[filter]
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  leafHeightSpacing <- 0.04
  
  ## compute (transitive) cluster members, cluster positions, and leaf heights
  innerNodeMembers <- list()
  innerNodeFeaturesBinary <- list()
  innerNodeFeaturesCount <- list()
  innerNodePosition <- vector(mode = "numeric", length = numberOfInnerNodes)
  leafHeights <- vector(mode = "numeric", length = numberOfPrecursorsFiltered)
  for(i in 1:numberOfInnerNodes)
    innerNodeMembers[[i]] <- NA
  
  goAhead <- TRUE
  while(goAhead){
    goAhead <- FALSE
    ## iterate each uncomputed inner node
    for(i in which(is.na(innerNodeMembers))){
      ## check if all cluster-children have been computed
      if(!any(is.na(innerNodeMembers[unlist(cluster$merge[i, ])[unlist(cluster$merge[i, ]) > 0]]))){
        leafs <- -unlist(cluster$merge[i, ])[unlist(cluster$merge[i, ]) < 0]
        nodes <-  unlist(cluster$merge[i, ])[unlist(cluster$merge[i, ]) > 0]
        ## cluster members
        innerNodeMembers[[i]] <- unlist(c(
          leafs, 
          unlist(innerNodeMembers[nodes])
        ))
        ## cluster position
        innerNodePosition[[i]] <- mean(unlist(c(
          match(x = leafs, table = cluster$order), 
          unlist(innerNodePosition[nodes])
        )))
        ## leaf heights
        leafHeights[leafs] <- cluster$height[[i]]
        
        ## features
        innerNodeFeaturesBinary[[i]] <- apply(
          X = dataList$featureMatrixBinary[filter[innerNodeMembers[[i]]], ], MARGIN = 2, FUN = all
        )
        innerNodeFeaturesCount[[i]] <- length(which(innerNodeFeaturesBinary[[i]]))
      }
      goAhead <- TRUE
    }
  }
  
  ## dendrogram leaf ends for normal plot
  #leafHeights <- leafHeights - leafHeightSpacing
  leafHeights <- rep(x = 0, times = length(leafHeights))
  
  ## compute x- and y-coordinates and point-labels
  coordinatesX <- unlist(c(innerNodePosition, match(x = 1:numberOfPrecursorsFiltered, table = cluster$order)))
  coordinatesY <- unlist(c(cluster$height, leafHeights))
  labels <- unlist(c(-(1:numberOfInnerNodes), 1:numberOfPrecursorsFiltered))
  text <- unlist(c(innerNodeFeaturesCount, apply(X = dataList$featureMatrixBinary[filter, ], MARGIN = 1, FUN = function(X) length(which(X)))))
  
  ##########################################
  ## box
  if(progress)  incProgress(amount = 0.1, detail = "Boxing")
  filterList <- list()
  ## parameter
  filterList$groupOne <- groupOne
  filterList$groupTwo <- groupTwo
  filterList$filter_average <- filter_average
  filterList$filter_lfc <- filter_lfc
  filterList$filter_ms2_masses <- filter_ms2_masses
  filterList$filter_ms2_ppm <- filter_ms2_ppm
  ## filter
  filterList$filter <- filter
  filterList$numberOfPrecursorsFiltered <- numberOfPrecursorsFiltered
  ## cluster
  filterList$innerNodeMembers <- innerNodeMembers
  filterList$innerNodeFeaturesBinary <- innerNodeFeaturesBinary
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
getMS2spectrum <- function(dataList, filterList, label){
  if(label > 0){
    ## leaf
    precursorIndex <- filterList$filter[[label]]
    precursorMass  <- as.numeric(dataList$dataFrame$"m/z"[[precursorIndex]])
    adduct <- dataList$dataFrame$"Adduct.ion.name"[[precursorIndex]]
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
    
    features <- dataList$featureMatrixBinary[precursorIndex, ]
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
    
    infoText <- paste("Precursor ''", filterList$cluster$labels[[label]], "'' has ", length(fragmentsX), " fragments", sep = "")
    #cat(paste(paste(sort(fragmentsX), collapse = ", "), "\n", sep = ""))
    if(!is.na(neutralMass)){
      #writeClipboard(landingPageUrl, format = 1)
      landingPageUrlForLink <- landingPageUrl
    }
    else
      landingPageUrlForLink <- NULL
  } else {
    ## inner node
    clusterIndex <- -label
    clusterMembers <- sort(filterList$filter[filterList$innerNodeMembers[[clusterIndex]]])
    
    features <- filterList$innerNodeFeaturesBinary[[clusterIndex]]
    fragmentsX <- dataList$fragmentMasses[features]
    fragmentsY <- apply(X = data.matrix(dataList$featureMatrix[clusterMembers, features]), MARGIN = 2, FUN = mean)
    
    infoText <- paste("Cluster ''", paste(dataList$precursorLabels[clusterMembers], sep = ", ", collapse = ", "), "'' has ", length(fragmentsX), " fragments in common", sep = "")
    #cat(paste(paste(sort(fragmentsX), collapse = ", "), "\n", sep = ""))
    landingPageUrlForLink <- NULL
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
  
  return(resultObj)
}

#########################################################################################
## plotting helper
colorSubTree <- function(
  cluster,
  index, # index of the current node
  lwd = 1,
  lty = 1,
  col = "black"
){
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
#########################################################################################
## plotting
calcClusterPlotDendrogram <- function(filterList, nodeIndex = NULL){
  ####################
  ## cluster
  par(mar=c(7,4,2,1), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  dend <- as.dendrogram(filterList$cluster)
  if(!is.null(nodeIndex)){
    if(nodeIndex < 0){
      clusterMembers <- filterList$innerNodeMembers[[-nodeIndex]]
    } else {
      clusterMembers <- c(nodeIndex)
    }
    
    colLab <- function(n) {
      if(is.leaf(n)) {
        a <- attributes(n)
        nodeLabelHere <- a$label
        nodeIndexHere <- match(x = nodeLabelHere, table = filterList$cluster$labels)
        
        if(length(na.omit(match(x = nodeIndexHere, table = clusterMembers))) > 0)
          attr(n, "nodePar") <- c(a$nodePar, lab.col = 'blue') #   change the node color
      }
      return(n)
    }
    dend <- dendrapply(dend, colLab)
  }
  
  plot(x = dend, xlab = "", ylab = "Jaccard distance", main = "Precursor cluster dendrogram", sub = "")
  
  if(!is.null(nodeIndex)){
    if(nodeIndex < 0){
      ## sub tree coloring
      a2r_counter <<- min(match(x = clusterMembers, table = filterList$cluster$order)) - 1
      colorSubTree(cluster = filterList$cluster, index = -nodeIndex, col = "blue")
    }
  }
  
  for(i in 1:length(filterList$poiCoordinatesX)){
    points(x = filterList$poiCoordinatesX[[i]], y = filterList$poiCoordinatesY[[i]], col = "red", pch=19, cex=1.)
    text(  x = filterList$poiCoordinatesX[[i]], y = filterList$poiCoordinatesY[[i]] + 0.02, labels = filterList$poiText[[i]], pos = 4)
  }
}
calcClusterPlotHeatmap <- function(dataList, filterList){
  ####################
  ## heatmap
  if(is.na(filterList$groupOne) | is.na(filterList$groupTwo)){
    columnsOfInterest <- c()
  } else {
    columnsOfInterest <- c(
      dataList$dataColumnNameFunction(filterList$groupOne), dataList$dataColumnNameFunction(filterList$groupTwo), 
      dataList$lfcColumnNameFunction(groupOne = filterList$groupOne, groupTwo = filterList$groupTwo)
    )
  }
  
  ## horizontal offset to match the cluster dendrogram
  offset <- 1
  offset <- grconvertX(offset, from = "user", to = "device")
  offset <- offset / (filterList$numberOfPrecursorsFiltered * 13)
  
  par(mar=c(0,4 + offset,1,1 + offset), mgp = c(3, 1, 0) + offset)  ## c(bottom, left, top, right)
  
  if(length(columnsOfInterest) == 0){
    plot.new()
    #cimage(zcol = dataList$colorMatrixDataFrame[, columnsOfInterest], axes = FALSE, ylab = "Groups")
  } else {
    cimage(zcol = dataList$colorMatrixDataFrame[filterList$filter, columnsOfInterest], axes = FALSE, ylab = "Groups")
    axis(side = 2, at = c(1, 2, 3), labels = c(filterList$groupOne, filterList$groupTwo, "LFC"), las = 2, tick = TRUE)
  }
}
calcClusterPlotHeatmapLegend <- function(dataList){
  ####################
  ## heatmap legend
  par(mar=c(0,0,0,2), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  
  legend_imageAbs <- as.raster(t(matrix(cmap(seq(from = 0, to = dataList$absMax, length.out = 100), dataList$colorMapAbsoluteData), nrow=1)))
  legend_imageLFC <- as.raster(t(matrix(cmap(seq(from = -dataList$logFoldChangeMax, to = dataList$logFoldChangeMax, length.out = 100), dataList$colorMapLogFoldChange), nrow=1)))
  
  plot.new()
  epsilon <- 0.05
  rasterImage(image = legend_imageAbs, xleft = 0, ybottom = 0, xright = 0.4, ytop = 0.5-epsilon)
  text(x=0.8, y = seq(0,0.5-epsilon,l=5), labels = format(seq(0, dataList$absMax, l=5), scientific = TRUE, digits = 0))
  
  rasterImage(image = legend_imageLFC, xleft = 0, ybottom = 0.5+epsilon, xright = 0.4, ytop = 1)
  text(x=0.8, y = seq(0.5+epsilon,1,l=5), labels = format(x = seq(-dataList$logFoldChangeMax, dataList$logFoldChangeMax, l=5), digits = 0))
}
calcClusterPlotMS2 <- function(dataList, fragmentsX = c(), fragmentsY = c()){
  ####################
  ## fragment spectrum
  if(is.null(fragmentsX) | is.null(fragmentsY)){
    fragmentsX = c()
    fragmentsY = c()
  }
  
  par(mar=c(5,4,1,1), mgp = c(3, 1, 0))  ## c(bottom, left, top, right)
  plot(x = fragmentsX, y = fragmentsY, ylab = "Abundance", xlab = "m/z", xlim = c(dataList$minimumMass, dataList$maximumMass), ylim = c(0, 1.1))
  
  tickPositions <- fragmentsX
  minimumTickLabelShift <- (dataList$maximumMass - dataList$minimumMass) / 100
  if(length(tickPositions) > 1)
    for(i in 2:length(tickPositions))
      if(tickPositions[[i]] - tickPositions[[i - 1]] < minimumTickLabelShift)
        tickPositions[[i]] <- tickPositions[[i - 1]] + minimumTickLabelShift
  if(length(fragmentsX) > 0){
    axis(side = 1, at = fragmentsX, labels = FALSE, las = 2)
    axis(side = 1, at = tickPositions, labels = fragmentsX, las = 2, tick = FALSE)
    for(idx in 1:length(fragmentsX)){
      points(x = fragmentsX[[idx]], y = fragmentsY[[idx]], col = "red", type = "h", lwd=4)
      points(x = fragmentsX[[idx]], y = fragmentsY[[idx]], col = "red", pch=19, cex=0.7)
      #text(x = fragmentsX[[idx]], y = fragmentsY[[idx]], labels = fragmentsX[[idx]], pos = 4)
    }
  }
}
