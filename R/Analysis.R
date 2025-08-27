
#' Filter dataset based on intensity mean values
#' 
#' @param dataList list object
#' @param filter_average numeric
#' @param sampleClasses vector of sample classes
#'
#' @returns boolean vector
#' @export
filterThreshold <- function(dataList, filter_average, sampleClasses = dataList$sampleClasses) {
  
  # need to fix code if false
  stopifnot(identical(sampleClasses, as.vector(sampleClasses)))
  
  mean_colnames <- sapply(X = sampleClasses,
                          FUN = dataList$dataMeanColumnNameFunctionFromName)
  
  unname(
    apply(X = dataList$dataFrameMeasurements[, mean_colnames],
          MARGIN = 1, FUN = mean) >= filter_average
  )

}


#' Filter dataset based on LFC of two groups
#'
#' @param dataList list object
#' @param filter_lfc numeric e.g. 2 or -2
#' @param sampleClasses vector of names for group1 and group2
#'
#' @returns boolean vector
#' @export
filterLFC <- function(dataList, filter_lfc, sampleClasses = dataList$sampleClasses) {

  stopifnot("The number of sampleClasses for LFC is not equal to two!" = length(sampleClasses) == 2)
  
  test_lfc <- if(filter_lfc > 0) filter_lfc else -filter_lfc
  
  dataList$dataFrameMeasurements[
    , dataList$lfcColumnNameFunctionFromName(sampleClasses[[1]], sampleClasses[[2]])] >= test_lfc
  
}

#' Create a filterObject for console use
#' 
#' Mixes code from filterData() and doPerformFiltering_impl() to provide
#' a simple way to create a compatible filterObj in a script.
#' 
#' doPerformFiltering_impl() relies on a DataList in the global environment, 
#' so it can't be used in a script.
#'
#' @param dataList list object
#' @param filter_average numeric
#' @param filter_lfc numeric
#' @param sampleClasses character vector
#'
#' @returns filterObject
#' @export
makeFilterObj <- function(dataList, filter_average = NULL, filter_lfc = NULL, sampleClasses = NULL) {
  
  if(is.null(sampleClasses)) sampleClasses <- dataList$sampleClasses
  sampleSet <- dataList$groupSampleDataFrame$Sample
  
  filterHere <- filterData(
    dataList, sampleClasses = dataList$sampleClasses, 
    sampleSet = sampleSet, filterBySamples = FALSE,
    filter_average = filter_average, filter_lfc = filter_lfc,
    filterList_ms2_masses = NULL, filter_ms2_ppm = NULL, 
    filter_ms1_masses = NULL, filter_ms1_ppm = NULL, 
    includeIgnoredPrecursors = FALSE, progress = FALSE)
  
  # gp: the rest is brute force to reproduce saved filterObj
  filterHere$filter <- filterHere$filter %>% purrr::set_names(dataList$precursorLabels[.])
  
  filterHere$groupSetOriginal <- sampleClasses
  filterHere$sampleSetOriginal <- sampleSet
  filterHere$filterBySamplesOriginal <- FALSE
  
  filterHere$filter_averageOriginal <- if(is.null(filter_average)) 0 else as.character(filter_average)
  filterHere$filter_lfcOriginal <- if (is.null(filter_lfc)) 0 else as.character(filter_lfc)
  
  filterHere$filter_ms2_masses1Original <- list()
  filterHere$filter_ms2_masses2Original <- list()
  filterHere$filter_ms2_masses3Original <- list()
  filterHere$filter_ms2_ppmOriginal <- ""
  filterHere$filter_ms1_massesOriginal <- list()
  filterHere$filter_ms1_ppmOriginal <- ""
  
  filterHere$includeIgnoredPrecursorsOriginal <- FALSE
  
  filterHere
  
}


#' Filter a dataList
#'
#' @param dataList List object
#' @param sampleClasses sample class
#' @param sampleSet ?
#' @param filterBySamples boolean
#' @param filter_average numeric Features that have a mean intensity below that threshold are filtered out.
#' @param filter_lfc numeric Log-fold-change threshold between two treatments. Only implement when exactly two groups are present.
#' @param filterList_ms2_masses ?
#' @param filter_ms2_ppm ?
#' @param filter_ms1_masses ?
#' @param filter_ms1_ppm ?
#' @param includeIgnoredPrecursors boolean
#' @param progress boolean
#'
#' @returns A filtered dataList object
#' @export
filterData <- function(dataList, sampleClasses, sampleSet = NULL, filterBySamples = FALSE, filter_average = NULL,
                       filter_lfc = NULL, filterList_ms2_masses = NULL, filter_ms2_ppm = NULL, 
                       filter_ms1_masses = NULL, filter_ms1_ppm = NULL, 
                       includeIgnoredPrecursors = FALSE, progress = FALSE){
  
  ## filter
  filter <- rep(x = TRUE, times = dataList$numberOfPrecursors)
  
  ## filter_average
  if(!is.null(filter_average)){
    # CHECK gp: Does it need to consider filterBySamples?
    filter <- filter & filterThreshold(dataList, filter_average, sampleClasses)
  }
  
  ## filter_lfc
  if(!is.null(filter_lfc)){
    # CHECK gp: strange condition, leaving it for now
    if(filter_lfc != 0){
      filter <- filter & filterLFC(dataList, filter_lfc, sampleClasses)
    }
  }
  
  ## filter_ms2_masses, filter_ms2_ppm
  if(!is.null(filterList_ms2_masses) && !is.null(filter_ms2_ppm) && length(filterList_ms2_masses) > 0){
    error <- abs(dataList$fragmentMasses) * filter_ms2_ppm / 1E6
    filterFragmentLists <- rep(x = FALSE, times = dataList$numberOfPrecursors)
    # BUG filter_ms2_masses not defined anywhere!!!
    for(filter_ms2_masses in filterList_ms2_masses){
      filterTempFragmentList <- rep(x = TRUE, times = dataList$numberOfPrecursors)
      for(fragmentIndex in seq_len(length(filter_ms2_masses))){
        distances <- abs(dataList$fragmentMasses - filter_ms2_masses[[fragmentIndex]])
        filterTempFragment <- rep(x = FALSE, times = dataList$numberOfPrecursors)
        
        ## set of fragment columns which are conjunct by or
        columns <- which(distances <= error)
        #columns <- which.min(distances <= error)
        for(column in columns)
          filterTempFragment <- filterTempFragment | dataList$featureMatrix[, column] != 0
        
        ## set of fragemnt masses which are conjunct by and
        filterTempFragmentList <- filterTempFragmentList & filterTempFragment
      }
      
      ## set of mass lists which are conjunct by or
      filterFragmentLists <- filterFragmentLists | filterTempFragmentList
    }
    filter <- filter & filterFragmentLists
  }
  
  ## filter_ms1_masses, filter_ms1_ppm
  if(!is.null(filter_ms1_masses) && !is.null(filter_ms1_ppm) && length(filter_ms1_masses) > 0){
    precursorMasses <- as.numeric(dataList$dataFrameInfos$"m/z")
    error <- precursorMasses * filter_ms1_ppm / 1E6
    
    filterMS1masses <- rep(x = FALSE, times = dataList$numberOfPrecursors)
    for(precursorMassIndex in seq_len(length(filter_ms1_masses))){
      distances <- abs(precursorMasses - filter_ms1_masses[[precursorMassIndex]])
      filterPart <- distances <= error
      filterPart[is.na(filterPart)] <- FALSE
      filterMS1masses <- filterMS1masses | filterPart
    }
    filter <- filter & filterMS1masses
  }
  
  ## include ignored precursors
  if(!includeIgnoredPrecursors) {
    filter <- filter & !dataList$annoArrayIsArtifact
  }
  
  filter <- which(filter)
  
  resultObj <- list()
  resultObj$filter <- filter
  resultObj$numberOfPrecursors <- dataList$numberOfPrecursors
  resultObj$numberOfPrecursorsFiltered <- length(filter)
  
  if(is.null(sampleClasses)){
    resultObj$sampleClasses    <- list()
    resultObj$sampleSet <- list()
    resultObj$filterBySamples <- NA
  } else {
    resultObj$sampleClasses    <- sampleClasses
    resultObj$sampleSet <- sampleSet
    resultObj$filterBySamples <- filterBySamples
  }
  
  #resultObj$sampleClasses                   <- ifelse(test = is.null(sampleClasses),                   yes = NA, no = sampleClasses)
  resultObj$filter_average           <- ifelse(test = is.null(filter_average),           yes = 0, no = filter_average)
  resultObj$filter_lfc               <- ifelse(test = is.null(filter_lfc),               yes = 0, no = filter_lfc)
  if(is.null(filterList_ms2_masses)){
    resultObj$filterList_ms2_masses    <- list()
  } else {
    resultObj$filterList_ms2_masses    <- filterList_ms2_masses
  }
  #resultObj$filterList_ms2_masses    <- ifelse(test = is.null(filterList_ms2_masses),    yes = list(), no = filterList_ms2_masses)
  resultObj$filter_ms2_ppm           <- ifelse(test = is.null(filter_ms2_ppm),           yes = "", no = filter_ms2_ppm)
  if(is.null(filter_ms1_masses)){
    resultObj$filter_ms1_masses    <- list()
  } else {
    resultObj$filter_ms1_masses    <- filter_ms1_masses
  }
  #resultObj$filter_ms1_masses        <- ifelse(test = is.null(filter_ms1_masses),        yes = "", no = filter_ms1_masses)
  resultObj$filter_ms1_ppm           <- ifelse(test = is.null(filter_ms1_ppm),           yes = "", no = filter_ms1_ppm)
  resultObj$includeIgnoredPrecursors <- ifelse(test = is.null(includeIgnoredPrecursors), yes = NA, no = includeIgnoredPrecursors)
  
  return (resultObj)
}

#' Calculate a distance matrix between MS/MS spectra in the MetFamily feature Matrix
#'
#' Computes a distance matrix between MS/MS spectra in the feature matrix representation using
#' various distance or similarity measures.
#' The function supports multiple distance metrics, of which the following are also available 
#' in the MetFamily GUI: "Jaccard", "Jaccard (intensity-weighted)", "Jaccard (fragment-count-weighted)"
#' and "NDP (Normalized dot product)". The following list defines the metrics and their characteristics:
#'
#' - "Jaccard": Jaccard distance \eqn{1 - |A \cap B| / |A \cup B|}: the fraction of matching fragments among all fragments,
#'   where A and B are MS/MS features of two different precursors.
#'
#' - "Jaccard (intensity-weighted)": Jaccard distance weighted by the intensity of 
#'   features. First, intensities are discretized: [0.01-0.2[ are set down to 0.01, [0.2, 0.4[ down to 0.2, 
#'   and >= 0.4 increased to 1. For matching fragments the higher intensity is used. Then, Jacard becomes the sum of matching intensities among the sum of all intensities in A and B.
#'   \eqn{1 - \sum_{i\in matches} max(A_i, B_i) / (sum(A) + sum(B))}
#'
#' - "Jaccard (fragment-count-weighted)": Jaccard distance weighted by relative occurance of the fragments among the precursors after filtering
#'   counts: \eqn{1 - \sum_{i\in matches} freq(A_i) / (sum_{\notin matches} freq(A) + sum_{\notin matches}(B))}
#'
#' - "NDP (Normalized dot product)": Normalized dot product similarity: \eqn{NDP = \frac{\left( \sum_{i}^{\text{S1\&S2}} W_{\text{S1},i} W_{\text{S2},i} \right)^2}{\sum_i W_{\text{S1},i}^2 \sum_i W_{\text{S2},i}^2}}
#'   as described in Gaquerel et al. 2015. (10.1073/pnas.1610218113)[https://www.pnas.org/doi/10.1073/pnas.1610218113#sec-4-5]
#'
#' @param dataList List object containing precursor, feature and MS/MS data.
#' @param filter Logical or integer vector indicating for which precursors to include the MS/MS spectra.
#' @param distanceMeasure Character string specifying the distance metric to use.
#'        Supported values include "Jaccard", "Manhatten", "NDP (Normalized dot product)", and others.
#' @param progress Logical or NA. If TRUE, progress is reported via incProgress().
#'
#' @return A list with elements:
#'   \describe{
#'     \item{distanceMatrix}{A numeric matrix of pairwise distances.}
#'     \item{filter}{The filter vector used, same as the input parameter.}
#'     \item{distanceMeasure}{The distance metric used, same as the input parameter.}
#'   }
#' @export
calculateDistanceMatrix <- function(dataList, filter, distanceMeasure = "Jaccard", progress = FALSE,
                                    minNumberFragments = 5, removePrecursor = TRUE) {
  
  stopifnot(is.logical(progress))
  
  numberOfPrecursors <- length(filter)
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Distances 0 / ", numberOfPrecursors, sep = ""))
  ## compute distance matrix:
  lastOut <- proc.time()["user.self"]
  lastPrecursor <- 1
  
  distanceMatrix <- NULL
  switch(distanceMeasure,
         "Jaccard"={
           featureIndeces <- dataList$featureIndeces[filter]
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }

             for(j in seq_len(numberOfPrecursors)){
               if(i == j){
                 distanceMatrix[i, j] <- 0
                 next
               }
               
               intersectionCount <- sum(featureIndeces[[i]] %in% featureIndeces[[j]])
               distanceMatrix[i, j] <- 1 - intersectionCount / (length(featureIndeces[[i]]) + length(featureIndeces[[j]]) - intersectionCount)
             }
           }
         },
         "Jaccard (intensity-weighted pure)"={
           featureIndeces <- dataList$featureIndeces[filter]
           featureMatrix <- dataList$featureMatrix[filter, ]
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
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
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             
             featureMatrixBinaryHere <- t(featureMatrixBinary)
             
             ## intersecting features
             intersections <- featureMatrixBinary[i, ] & featureMatrixBinaryHere
             onlyIs <- featureMatrixBinary[i, ] & (!featureMatrixBinaryHere)
             onlyJs <- featureMatrixBinaryHere  & (!featureMatrixBinary[i, ])
             
             onlyIs[is.na(onlyIs)] <- FALSE
             onlyJs[is.na(onlyJs)] <- FALSE
             
             sumOnlyIs  <- apply(X = onlyIs, MARGIN = 2, FUN = function(x){ sum(featureMatrix[i, ] & x) })

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
           
           distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
               if(i == j){
                 distanceMatrix[i, j] <- 0
                 next
               }
               
               ## intersecting features
               intersection <- featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
               
               if(length(intersection) == 0){
                 ## max distance
                 distance <- 1 # sumOnlyI + sumOnlyJ
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
                 
                 maxIntensity <- apply(X = featureMatrix[c(i, j), intersection, drop = FALSE], MARGIN = 2, FUN = max)
                 intersectionSum <- sum(maxIntensity)
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
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
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
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
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
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
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
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (counter - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- counter
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", counter, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", counter, " / ", numberOfPrecursors, sep = ""))
             }
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
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
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
           for(i in seq_len(numberOfPrecursors)){
             time <- proc.time()["user.self"]
             if(time - lastOut > 1){
               lastOut <- time
               precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
               lastPrecursor <- i
               if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress,     detail = paste("Distance ", i, " / ", numberOfPrecursors, sep = "")) else print(paste("Distance ", i, " / ", numberOfPrecursors, sep = ""))
             }
             for(j in seq_len(numberOfPrecursors)){
               if(i == j){
                 similarityMatrix[i, j] <- 0
                 next
               }
               
               intersection    <- featureIndeces[[i]][featureIndeces[[i]] %in% featureIndeces[[j]]]
               intersectionSum <- sum(sqrt(featureMatrix[i, intersection]) * dataList$fragmentMasses[intersection]^2 * sqrt(featureMatrix[j, intersection]) * dataList$fragmentMasses[intersection]^2)^2
               iSum            <- sum((sqrt(featureMatrix[i, featureIndeces[[i]]]) * dataList$fragmentMasses[featureIndeces[[i]]]^2)^2)
               jSum            <- sum((sqrt(featureMatrix[j, featureIndeces[[j]]]) * dataList$fragmentMasses[featureIndeces[[j]]]^2)^2)
               
               similarityMatrix[i, j] <- intersectionSum / (iSum * jSum)
             }
           }
           distanceMatrix <- max(similarityMatrix) - similarityMatrix
         },
         # (Modified) Cosine (with NL)
         "Cosine" = {
           distanceMatrix <- impl_cosine_similarity(
             dataList, filter, progress = progress, allow_shift = FALSE, nl = FALSE,
             rm_precursor = removePrecursor, num_filter = minNumberFragments)
         },
         "Cosine (with NL)" = {
           distanceMatrix <- impl_cosine_similarity(
             dataList, filter, progress = progress, allow_shift = FALSE, nl = TRUE,
             rm_precursor = removePrecursor, num_filter = minNumberFragments)
         },
         "Modified Cosine" = {
           distanceMatrix <- impl_cosine_similarity(
             dataList, filter, progress = progress, allow_shift = TRUE, nl = FALSE,
             rm_precursor = removePrecursor, num_filter = minNumberFragments)
         },
         "Modified Cosine (with NL)" = {
           distanceMatrix <- impl_cosine_similarity(
             dataList, filter, progress = progress, allow_shift = TRUE, nl = TRUE,
             rm_precursor = removePrecursor, num_filter = minNumberFragments)
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


#' Implement Cosine Similarity
#'
#' @param dataList dataset
#' @param filter filter vector
#' @param progress boolena
#' @param allow_shift boolean, modified cosine when true
#' @param nl boolean, include neutral losses when true
#' @param num_filter integer, minimal number of fragments required to perform comparison
#' @param rm_precursor boolean, remove precursor ions from ms2 spectra
#'
#' @returns distance matrix
#' @export
impl_cosine_similarity <- function(dataList, filter, progress = FALSE,
                        allow_shift = FALSE, nl = FALSE, num_filter=5, rm_precursor=T) {
  
  numberOfPrecursors <- length(filter)
  
  featureIndeces <- dataList$featureIndeces[filter]
  featureMatrix <- dataList$featureMatrix[filter, ]
  
  precursorMasses <- as.numeric(dataList$dataFrameInfos$`m/z`)[filter]
  
  # Convert feature indeces to fragment mzs
  featureMasses <- lapply(featureIndeces, function(idxs) {
    dataList$fragmentMasses[idxs]
  })
  
  # Remove precursorMasses from featureMasses if needed
  if (rm_precursor){
    featureMasses <- mapply(function(fm, pm) {
      fm[abs(fm - pm) > 1] #precursor tolerance window set to 1
    }, featureMasses, precursorMasses, SIMPLIFY = FALSE)
  }
  
  lastOut <- proc.time()["user.self"]
  lastPrecursor <- 1
  
  distanceMatrix <- matrix(nrow = numberOfPrecursors, ncol = numberOfPrecursors)
  
  # Calculate cosine similarity when spectra have at least # of fragments
  keep_idx <- which(sapply(featureMasses, length) >= num_filter)
  
  for (i in 1:(numberOfPrecursors-1)) {
    time <- proc.time()["user.self"]
    if(time - lastOut > 1){
      lastOut <- time
      precursorProgress <- (i - lastPrecursor) / numberOfPrecursors
      lastPrecursor <- i
      if(!is.na(progress))  if(progress)  incProgress(amount = precursorProgress, detail = paste("Distance ", i, " / ", numberOfPrecursors)) else print(paste("Distance ", i, " / ", numberOfPrecursors))
    }
    
    # Skip the row if spectra1 doesn't meet fragment # filter
    if (!(i %in% keep_idx)){
      distanceMatrix[i, ] <- 1
      next
    }
    
    # Get fragments mz and normalized intensities for spectra1
    mz1 <- featureMasses[[i]]
    intensity1 <- featureMatrix[i,]
    intensity1 <- intensity1[names(intensity1) %in% mz1]
    precursor_mz1 <- precursorMasses[i]
    
    for (j in (i+1):numberOfPrecursors) {
      # If spectra2 doesn't meet fragment # filter, put spectra1_spectra2 distance as 1
      if (!j %in% keep_idx){
        distanceMatrix[i, j] <- 1
        next
      }
      
      # Get fragments mz and normalized intensities for spectra
      mz2 <- featureMasses[[j]]
      intensity2 <- featureMatrix[j,]
      intensity2 <- intensity2[names(intensity2) %in% mz2]
      precursor_mz2 <- precursorMasses[j]
      
      # Calculate cosine score: 
      cosine_score <- 
        cosine_similarity(mz1, intensity1, precursor_mz1, 
                          mz2, intensity2, precursor_mz2, 
                          fragment_mz_tolerance = 0.01, allow_shift = allow_shift, 
                          nl = nl, normalize = FALSE, method = 2)
      
      # Convert similarity to distance
      distanceMatrix[i, j] <- 1 - cosine_score
    }
  }
  
  # Fill diagonal as 0 (identical), mirror upper half to lower half
  diag(distanceMatrix) <- 0
  distanceMatrix[lower.tri(distanceMatrix)] <- 
    t(distanceMatrix)[lower.tri(distanceMatrix)]
  
  distanceMatrix
  
}


#' Calculate cosine/modified cosine score
#'
#' @param mz1 numeric vector
#' @param intensity1 numeric vector
#' @param precursor_mz1 numeric
#' @param mz2 numeric vector
#' @param intensity2 numeric vector
#' @param precursor_mz2 numeric
#' @param fragment_mz_tolerance numeric
#' @param allow_shift boolean
#' @param nl boolean
#' @param normalize boolean
#' @param method 1 = normalise by euclidean norm, 2 = sqrt normalised intensity
#' 
#' @importFrom clue solve_LSAP
#'
#' @returns score
#' @export
cosine_similarity <- function(
    mz1, intensity1, precursor_mz1, 
    mz2, intensity2, precursor_mz2, 
    fragment_mz_tolerance = 0.01,
    allow_shift = TRUE,
    nl = FALSE,
    normalize = FALSE, 
    method = 2) {
  
  # Determine if including neutral loss
  if (!nl){
    fragidx1 <- which(mz1>0)
    mz1 <- mz1[fragidx1]
    intensity1 <- intensity1[fragidx1]
    
    fragidx2 <- which(mz2>0) 
    mz2 <- mz2[fragidx2]
    intensity2 <- intensity2[fragidx2]
  }
  
  # Normalize intensity vectors to unit norm (if needed)
  if (normalize & method==1){
    norm1 <- sqrt(sum(intensity1^2))
    norm2 <- sqrt(sum(intensity2^2))
    if (norm1 == 0 || norm2 == 0) return(0)
    intensity1 <- intensity1 / norm1
    intensity2 <- intensity2 / norm2
  }
  
  # Determine if shifting (modified cosine) is used
  mass_diff <- 0
  precursor_mass_diff <- precursor_mz2 - precursor_mz1
  if (allow_shift && abs(precursor_mass_diff) >= fragment_mz_tolerance) {
    mass_diff <- c(0, precursor_mass_diff) 
  }
  
  # Build cost matrix
  cost_matrix <- matrix(0, nrow = length(mz1), ncol = length(mz2))
  
  for (i in seq_along(mz1)) {
    for (z in seq_along(mass_diff)) {
      mz2_shifted <- mz2 - mass_diff[z]
      delta <- abs(mz1[i] - mz2_shifted)
      matched <- which(delta <= fragment_mz_tolerance)
      if (length(matched) >= 0) {
        if (method == 1){
          cost_matrix[i, matched] <- intensity1[i] * intensity2[matched]
        } else if (method == 2){
          # Method 2 normalization: sqrt of normalized intensity
          cost_matrix[i, matched] <- sqrt(intensity1[i])/sqrt(sum(intensity1)) * sqrt(intensity2[matched])/sqrt(sum(intensity2))
        }
      }
    }
  }
  
  # Solve linear sum assignment
  if(nrow(cost_matrix) > ncol(cost_matrix)) cost_matrix <- t(cost_matrix)
  assignment <- clue::solve_LSAP(cost_matrix, maximum = TRUE)
  matched_scores <- cost_matrix[cbind(seq_len(nrow(cost_matrix)), assignment)]
  score <- sum(matched_scores[matched_scores > 0])
  
  return(score)
}


#' Create cluster object
#'
#' @param dataList 
#' @param filterObj 
#' @param distanceMatrix 
#' @param method 
#' @param distanceMeasure 
#' @param progress 
#'
#' @returns clusterDataList object
#' @export
calculateCluster <- function(dataList, filterObj, distanceMatrix, method, distanceMeasure, progress = FALSE){
  
  numberOfPrecursorsFiltered <- length(filterObj$filter)
  ##########################################
  ## compute gui stuff
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "Clustering")
  ## compute and annotate cluster
  dist <- stats::as.dist(m = distanceMatrix)
  cluster <- hclust(d = dist, method = method)
  cluster$labels <- dataList$precursorLabels[filterObj$filter]
  numberOfInnerNodes <- numberOfPrecursorsFiltered - 1
  leafHeightSpacing <- 0.04
  
  ## optimal leaf ordering
  #opt <- order.optimal(dist = dist, merge = cluster$merge)
  #cluster$merge <- opt$merge
  #cluster$order <- opt$order
  
  ## compute (transitive) cluster members, cluster positions, and leaf heights
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.3, detail = "Analyze cluster")
  
  resultObj <- analyzeTreeFromRoot(dataList, cluster = cluster, filterObj$filter)
  
  innerNodeHeightIncreases <- unlist(resultObj$innerNodeHeightIncreases)
  innerNodeMembersTreeClusters <- resultObj$innerNodeMembersTreeClusters
  innerNodeMembersTreeLeaves  <- resultObj$innerNodeMembersTreeLeaves
  innerNodeMembersPrecursors <- resultObj$innerNodeMembersPrecursors
  #innerNodeFeaturesBinary <- resultObj$innerNodeFeaturesBinary
  innerNodeFeaturesIntersection <- resultObj$innerNodeFeaturesIntersection
  innerNodeFeaturesUnion <- resultObj$innerNodeFeaturesUnion
  innerNodeFeaturesCountsMatrix <- resultObj$innerNodeFeaturesCountsMatrix
  innerNodeFeaturesPresent <- resultObj$innerNodeFeaturesPresent
  #innerNodeFeaturesMeanAbundance <- resultObj$innerNodeFeaturesMeanAbundance
  #innerNodeFeaturesIntersectionCounter <- resultObj$innerNodeFeaturesIntersectionCounter
  #innerNodeFeaturesUnionCounter <- resultObj$innerNodeFeaturesUnionCounter
  innerNodePosition <- resultObj$innerNodePosition
  leafHeights <- resultObj$leafHeights
  
  #innerNodeFeaturesCountsMatrix <- matrix(data = unlist(innerNodeFeaturesCounts), nrow = length(innerNodeFeaturesCounts))
  
  ## dendrogram leaf ends for normal plot
  #leafHeights <- leafHeights - leafHeightSpacing
  leafHeights <- rep(x = 0, times = length(leafHeights))
  
  ## compute x- and y-coordinates and point-labels
  poiCoordinatesX <- unlist(c(innerNodePosition, match(x = seq_len(numberOfPrecursorsFiltered), table = cluster$order)))
  poiCoordinatesY <- unlist(c(cluster$height, leafHeights))
  
  precursorFeatureCount <- dataList$featureCount[filterObj$filter]
  #innerNodeUnionlabels <- as.character(innerNodeFeaturesUnionCounter)
  #innerNodeUnionlabels[!innerNodeHeightIncreases] <- ""
  #innerNodeIntersectionlabels <- as.character(innerNodeFeaturesIntersectionCounter)
  #innerNodeIntersectionlabels[!innerNodeHeightIncreases] <- ""
  
  drawPoi <- unlist(c(
    innerNodeHeightIncreases, 
    rep(x = TRUE, times = length(filterObj$filter))
  ))
  #poiUnion <- unlist(c(
  #  innerNodeUnionlabels, 
  #  precursorFeatureCount
  #))
  #poiIntersection <- unlist(c(
  #  innerNodeIntersectionlabels, 
  #  precursorFeatureCount
  #  #vector(mode = "character", length = length(filterObj$filter))
  #))
  poiIntersectionSmooth <- c(innerNodeFeaturesPresent, precursorFeatureCount)
  poiLabels <- unlist(c(seq_len(numberOfInnerNodes), -(seq_len(numberOfPrecursorsFiltered))))
  
  ##########################################
  ## box
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.5, detail = "Boxing")
  clusterDataList <- list()
  ## filter
  clusterDataList$filterObj <- filterObj
  clusterDataList$method <- method
  clusterDataList$distanceMeasure <- distanceMeasure
  #clusterDataList$distanceMatrix <- distanceMatrix
  clusterDataList$numberOfPrecursorsFiltered <- numberOfPrecursorsFiltered
  ## cluster
  
  clusterDataList$innerNodeMembersTreeLeaves <- innerNodeMembersTreeLeaves
  clusterDataList$innerNodeMembersTreeClusters <- innerNodeMembersTreeClusters
  clusterDataList$innerNodeMembersPrecursors <- innerNodeMembersPrecursors
  #clusterDataList$innerNodeFeaturesBinary <- innerNodeFeaturesBinary
  clusterDataList$innerNodeFeaturesIntersection <- innerNodeFeaturesIntersection
  clusterDataList$innerNodeFeaturesUnion <- innerNodeFeaturesUnion
  clusterDataList$innerNodeFeaturesCountsMatrix <- innerNodeFeaturesCountsMatrix
  #clusterDataList$innerNodeFeaturesPresent <- innerNodeFeaturesPresent
  #clusterDataList$innerNodeFeaturesMeanAbundance <- innerNodeFeaturesMeanAbundance
  #clusterDataList$innerNodeFeaturesIntersectionCounter <- innerNodeFeaturesIntersectionCounter
  #clusterDataList$innerNodeFeaturesUnionCounter <- innerNodeFeaturesUnionCounter
  clusterDataList$cluster <- cluster
  ## poi
  clusterDataList$numberOfPois    <- length(poiCoordinatesX)
  clusterDataList$poiCoordinatesX <- poiCoordinatesX
  clusterDataList$poiCoordinatesY <- poiCoordinatesY
  #clusterDataList$poiUnion        <- poiUnion
  #clusterDataList$poiIntersection <- poiIntersection
  clusterDataList$poiIntersectionSmooth <- poiIntersectionSmooth
  clusterDataList$poiLabels       <- poiLabels
  clusterDataList$drawPoi         <- drawPoi
  
  ##########################################
  ## calculate ms2 spectra
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Calculate spectrum information for leaves")
  ms2spectrumInfoForLeaves <- list()
  for(leafIdx in seq_len(numberOfPrecursorsFiltered)) {
    ms2spectrumInfoForLeaves[[leafIdx]] <- getMS2spectrumInfoForPrecursorLeaf(dataList, clusterDataList, treeLabel = -leafIdx, outAdductWarning = FALSE)
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Calculate consensus spectra for clusters")
  ms2spectrumInfoForClusters <- list()
  for(clusterIdx in seq_len(numberOfInnerNodes)) {
    ms2spectrumInfoForClusters[[clusterIdx]] <- getMS2spectrumInfoForCluster(dataList, clusterDataList, treeLabel = clusterIdx)
  }
  
  ## calculate cluster discriminativity
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Calculate cluster discriminativity")
  #clusterDiscriminativity <- vector(mode = "numeric", length = numberOfInnerNodes)
  clusterDiscriminativity <- unlist(suppressWarnings(lapply(X = ms2spectrumInfoForClusters, FUN = function(x){x$clusterDiscriminativity})))
  clusterDiscriminativity <- c(clusterDiscriminativity, rep(x = 0, times = numberOfPrecursorsFiltered))
  
  clusterDataList$ms2spectrumInfoForLeaves    <- ms2spectrumInfoForLeaves
  clusterDataList$ms2spectrumInfoForClusters  <- ms2spectrumInfoForClusters
  clusterDataList$clusterDiscriminativity     <- clusterDiscriminativity
  
  if(!is.na(progress))  if(progress)  setProgress(1)
  
  ## 1 167 983 544
  ##   127 968 176
  ##   127 922 016
  ##   128 042 056
  ##    13 761 512
  #print(sort( sapply(ls(),function(x){object.size(get(x))})))
  
  return(clusterDataList)
}


preprocessDataForPca <- function(dataFrame, scaling, logTransform){
  
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
             #print(x)
             (x - mean(x = x)) / sqrt(sd(x = x))
           }))
           
           dataFrame2[is.na(dataFrame2)] <- 0
         },
         stop(paste("Unknown scaling (", scaling, ")!", sep = ""))
  )
  
  #dataFrame2[is.na(dataFrame2)] <- dataFrame[is.na(dataFrame2)]
  dataFrame2[is.na(dataFrame2)] <- 0
  
  return(dataFrame2)
}

#' Run PCA
#'
#' @param dataList List object 
#' @param dataFrame2 
#' @param ms1AnalysisMethod 
#'
#' @returns A PCA list object
#' @export
performPca <- function(dataList, dataFrame2, ms1AnalysisMethod){
  
  print("######################################################################################")
  print(ms1AnalysisMethod)
  
 
  ## TODO pcaMethods confidence intervals analog to MetaboAnalyst: pcaMethods:::simpleEllipse
  minimumNumberOfComponents <- 5
  numberOfComponents <- min(minimumNumberOfComponents, nrow(dataFrame2))
  
  returnObj <- list()
  returnObj$ms1AnalysisMethod = ms1AnalysisMethod
  
  if(ncol(dataFrame2) < 3){
    ###########################################################
    ## no data or only one sample
    numberOfPrecursors <- dataList$numberOfPrecursors
    numberOfSamples    <- ncol(dataFrame2)
    returnObj$scores   <- matrix(nrow = numberOfSamples,    ncol = numberOfComponents, data = rep(x = 0, times = numberOfSamples    * numberOfComponents))
    returnObj$loadings <- matrix(nrow = numberOfPrecursors, ncol = numberOfComponents, data = rep(x = 0, times = numberOfPrecursors * numberOfComponents))
    returnObj$variance <- vector(mode = "numeric", length = numberOfComponents)
    return(returnObj)
  }
  
  ## choose library
  #ms1AnalysisMethod <- c(
  #  "stats",           # 1
  #  "FactoMineR",      # 2
  #  "pcaMethods",      # 3
  #  "mixOmics_pca",    # 4
  #  "mixOmics_spca",   # 5
  #  "mixOmics_plsda",  # 6
  #  "mixOmics_splsda"  # 7
  #)[[5]]
  
  switch(ms1AnalysisMethod,
         "PCA (Principal Component Analysis)"={
           #ms1AnalysisMethod <- "mixOmics_pca"
           ms1AnalysisMethod <- "pcaMethods"
         },
         "sPCA (Sparse Principal Component Analysis)"={
           ms1AnalysisMethod <- "mixOmics_spca"
         },
         "PLS-DA (Partial Least Squares Discriminant Analysis)"={
           ms1AnalysisMethod <- "mixOmics_plsda"
           #ms1AnalysisMethod <- "caret_plsda"
         },
         "sPLS-DA (Sparse Partial Least Squares Discriminant Analysis)"={
           ms1AnalysisMethod <- "mixOmics_splsda"
         },
         stop(paste("Unknown analysis method (", ms1AnalysisMethod, ")!", sep = ""))
  )
  
  print(paste("Analysis", ms1AnalysisMethod, sep = ": "))
  switch(ms1AnalysisMethod,
         "stats"={
           ## pca from "stats" package
           pca <- stats::prcomp(x = dataFrame2, retx = TRUE, center = FALSE, scale. = FALSE)
           returnObj$scores   <- pca$x
           returnObj$loadings <- pca$rotation
           returnObj$variance <- pca$sdev
         },
         "FactoMineR"={
           ## pca from "FactoMineR" package
           pca = FactoMineR::PCA(X = dataFrame2, graph = FALSE, scale.unit = FALSE, ncp = numberOfComponents)
           returnObj$scores   <- pca$ind$coord
           returnObj$loadings <- pca$var$coord
           returnObj$variance <- pca$eig$"percentage of variance"
         },
         "pcaMethods"={
           ## pca from "pcaMethods" package
           #pca <- pca(object = dataFrame2, method = "robustPca", nPcs = 2, scale = "none", center = FALSE, cv = "q2")
           pca <- pcaMethods::pca(object = dataFrame2, method = "svd", nPcs = numberOfComponents, scale = "none", center = FALSE)
           returnObj$scores   <- pca@scores
           returnObj$loadings <- pca@loadings
           returnObj$variance <- pca@sDev
           
           returnObj$R2 <- pca@R2
           #returnObj$Q2 <- pcaMethods::Q2(object = pca, fold=2)
           returnObj$Q2 <- tryCatch(
             {
               #Q2(object = pca, fold=2)
               pcaMethods::Q2(object = pca, verbose = FALSE)
             },
             error=function(cond) {
               rep(x = "N/A", times = numberOfComponents)
             }
           )
         },
         "mixOmics_pca"={
           ## pca from "mixOmics" package
           pca = mixOmics::pca(X = dataFrame2, ncomp = numberOfComponents, center = FALSE, scale = FALSE)
           returnObj$scores   <- pca$variates[[1]]
           returnObj$loadings <- pca$loadings[[1]]
           returnObj$variance <- pca$explained_variance
           #returnObj$R2 <- 
           #returnObj$Q2 <- 
         },
         "mixOmics_spca"={
           ## pca from "mixOmics" package
           
           # FIX error: 
           # Error in mixOmics::spca: (converted from warning) data
           # contain missing values which will be set to zero for calculations.
           # Consider using center = TRUE for better performance, or impute
           # missing values using 'impute.nipals' function.
           pca = mixOmics::spca(X = dataFrame2, ncomp = numberOfComponents, center = FALSE, scale = FALSE)
           returnObj$scores   <- pca$variates[[1]]
           returnObj$loadings <- pca$loadings[[1]]
           returnObj$variance <- pca$explained_variance
           #returnObj$R2 <- 
           #returnObj$Q2 <- 
           
           #performance <- mixOmics::perf(pca, validation = "loo", progressBar = FALSE)
           #val <- mixOmics::perf(pca, criterion = c("R2", "Q2"))
           
         },
         "mixOmics_plsda"={
           ## plsda "mixOmics" package
           groupLabels  <- unlist(lapply(X = rownames(dataFrame2), FUN = function(x){dataList$groupNameFunctionFromDataColumnName(dataColumnName = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))}))
           pca = mixOmics::plsda(X = dataFrame2, Y = groupLabels, ncomp = numberOfComponents, scale = FALSE)
           
           if(any(pca$explained_variance$X < 0.01)){
             maxComp <- max(which(pca$explained_variance$X >= 0.01))
             numberOfComponents <- maxComp
             
             pca = mixOmics::plsda(X = dataFrame2, Y = groupLabels, ncomp = numberOfComponents, scale = FALSE)
           }
           
           returnObj$scores   <- pca$variates[[1]]
           returnObj$loadings <- pca$loadings[[1]]
           returnObj$variance <- pca$explained_variance$X
           
           #performance <- mixOmics::perf(pca, validation = "loo", progressBar = FALSE)
           #performance$choice.ncomp
           
           if(FALSE){## R2 and Q2?
             performance <- mixOmics::perf(pca, validation = "Mfold", folds = 2, progressBar = FALSE, tol = 1e-20)
             performance <- mixOmics::perf(pca, validation = "loo", progressBar = FALSE, tol = 1e-20)
             
             pca = mixOmics::pca(X = dataFrame2, ncomp = numberOfComponents, center = FALSE, scale = FALSE)
             loadings <- pca$loadings[[1]]
             sumOfLoadings <- apply(X = loadings, MARGIN = 1, FUN = sum)
             toRemove <- which(abs(sumOfLoadings) < 0.0001)
             
             dataFrame3 <- dataFrame2[-toRemove, ]
             pca = mixOmics::plsda(X = dataFrame3, Y = groupLabels, ncomp = numberOfComponents, scale = FALSE)
             
             
             
             randomMatrix <- replicate(n = ncol(dataFrame2), runif(n = nrow(dataFrame2), min = -10000, max = 10000))
             dataFrame3 <- dataFrame2 + randomMatrix
             
             pca = mixOmics::plsda(X = dataFrame3, Y = groupLabels, ncomp = numberOfComponents, scale = FALSE)
             returnObj$scores   <- pca$variates[[1]]
             returnObj$loadings <- pca$loadings[[1]]
             returnObj$variance <- pca$explained_variance
             
             performance <- mixOmics::perf(pca, validation = "Mfold", folds = 2, progressBar = FALSE, tol = 1e-20)
             performance <- mixOmics::perf(pca, validation = "loo", progressBar = FALSE)
             mixOmics::perf(pca, validation = "loo", progressBar = FALSE)
           }
           
           #returnObj$R2 <- performance$R2
           #returnObj$Q2 <- performance$Q2
         },
         "mixOmics_splsda"={
           ## splsda from "mixOmics" package TODO
           groupLabels  <- unlist(lapply(X = rownames(dataFrame2), FUN = function(x){dataList$groupNameFunctionFromDataColumnName(dataColumnName = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))}))
           caret_splsda = mixOmics::splsda(X = dataFrame2, Y = groupLabels, ncomp = numberOfComponents, scale = FALSE)
           returnObj$scores   <- caret_splsda$variates[[1]]
           returnObj$loadings <- caret_splsda$loadings[[1]]
           returnObj$variance <- caret_splsda$explained_variance$X
         },
         "caret_plsda"={
           groupLabels  <- unlist(lapply(X = rownames(dataFrame2), FUN = function(x){dataList$groupNameFunctionFromDataColumnName(dataColumnName = x, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))}))
           caret_plsda <- caret::plsda(x = dataFrame2, y = as.factor(groupLabels), ncomp = numberOfComponents, probMethod = "softmax")
           
           returnObj$scores   <- caret_plsda$scores
           returnObj$loadings <- caret_plsda$loadings
           returnObj$variance <- caret_plsda$Xvar / sum(caret_plsda$Xvar)
           returnObj$accurracyContribution <- leaveOneOutCrossValidation_plsda(dataFrame2, groupLabels, numberOfComponents)
         },
         stop(paste("Unknown analysis method (", ms1AnalysisMethod, ")!", sep = ""))
  )
  
  #str(returnObj)
  
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
  ## artificial data two sampleClasses
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
    #numberOfSamples    <- nrow(returnObj$scores)
    numberOfSamples    <- ncol(dataFrame2)
    returnObj$scores   <- matrix(nrow = numberOfSamples, ncol = numberOfComponents)
    returnObj$loadings <- matrix(nrow = numberOfPrecursors, ncol = numberOfComponents)
    returnObj$variance <- vector(mode = "numeric", length = numberOfComponents)
  }
  return(returnObj)
}
leaveOneOutCrossValidation_plsda <- function(dataFrame2, groupLabels, numberOfComponents){
  ## leave-one-out-cross-validation
  maxNumberOfComponents_train <- min(numberOfComponents, nrow(dataFrame2) - 1 - 1)
  numberOfPositivePredictions_comp <- numeric(length = maxNumberOfComponents_train)
  for(sampleIdx in seq_along(groupLabels)){
    dataFrame2_train <- dataFrame2[-sampleIdx, ]
    groupLabels_train <- groupLabels[-sampleIdx]
    dataFrame2_test <- dataFrame2[sampleIdx, , drop=FALSE]
    groupLabels_test <- groupLabels[sampleIdx]
    
    for(numberOfComponents_train in seq_len(maxNumberOfComponents_train)){
      caret_plsda_train <- caret::plsda(
        x = dataFrame2_train, 
        y = as.factor(groupLabels_train), 
        ncomp = numberOfComponents_train, 
        probMethod = "softmax"
      )
      
      groupLabels_predict <- as.character(predict(caret_plsda_train, dataFrame2_test))
      numberOfPositivePredictions_comp[[numberOfComponents_train]] <- numberOfPositivePredictions_comp[[numberOfComponents_train]] + (groupLabels_predict == groupLabels_test)
    }
  }
  accurracy_comp <- numberOfPositivePredictions_comp / nrow(dataFrame2)
  accurracyContribution_comp <- numeric(length = maxNumberOfComponents_train + 1)
  accurracyContribution_comp[[1]] <- accurracy_comp[[1]]
  accurracyContribution_comp[2:maxNumberOfComponents_train] <- accurracy_comp[2:maxNumberOfComponents_train] - accurracy_comp[1:(maxNumberOfComponents_train-1)]
  accurracyContribution_comp[[maxNumberOfComponents_train + 1]] <- NA
  
  return(accurracyContribution_comp)
}

#' Calculate PCA
#'
#' @param dataList listObject
#' @param filterObj filterObject
#' @param ms1AnalysisMethod 
#' @param scaling 
#' @param logTransform 
#'
#' @returns PCA object
#' @export
calculatePCA <- function(dataList, filterObj, ms1AnalysisMethod, scaling, logTransform){

  ## data selection
  if(filterObj$filterBySamples){
    dataFrame <- dataList$dataFrameMeasurements[filterObj$filter, filterObj$sampleSet]
  } else {
    dataFrame <- dataList$dataFrameMeasurements[filterObj$filter, dataList$dataColumnsNameFunctionFromGroupNames(sampleClasses = filterObj$sampleClasses, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))]
  }
  dataFrame <- t(dataFrame)
  
  ## data scaling
  dataFrame2 <- preprocessDataForPca(dataFrame, scaling, logTransform)
  
  ## data analysis
  returnObj <- performPca(dataList, dataFrame2, ms1AnalysisMethod)
  returnObj$filterObj = filterObj
  returnObj$scaling = scaling
  returnObj$logTransform = logTransform
  
  return(returnObj)
}
