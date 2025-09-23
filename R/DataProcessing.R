# Utility functions used for data processing.
# Main processing functions were moved to import_export_project.R


#' Convert data.frame columns to numeric
#'
#' The data.numericmatrix() function works similar to base::data.matrix()
#' before R-4.0.0 converting character columns to numeric without converting 
#' to factor first, thus returning the actual numeric values.
#' 
#' @param x The data.frame to convert
#'
#' @return A matrix with all columns converted to numeric
#' @export
#'
#' @examples
#' data.numericmatrix(data.frame(a = c("1", "2", "3"), 
#'                               b = c("4", "5", "6")))
#' 
data.numericmatrix <- function(x) {
  for (i in 1:ncol(x)) {
    if (is.character(x[, i])) {
      x[, i] <- as.numeric(as.character(x[, i]))
    }
  }
  as.matrix(x)
}


#' @export
serializeSampleSelectionAndOrder <- function(groupSampleDataFrame)
{
  ## wrap columns
  columnsSerialized <- sapply(X = seq_len(ncol(groupSampleDataFrame)), FUN = function(colIdx){
    cellContent <- paste(groupSampleDataFrame[, colIdx], collapse = "; ")
    paste(colnames(groupSampleDataFrame)[[colIdx]], "=", "(", cellContent, ")", sep = "")
  })
  ## box
  groupSampleDataFrameName <- "SampleSelectionAndOrder"
  groupSampleDataFrameValue <- paste(columnsSerialized, collapse = "|")
  groupSampleDataFrameFieldValue <- paste(groupSampleDataFrameName, "=", "{", groupSampleDataFrameValue, "}", sep = "")
  
  return(groupSampleDataFrameFieldValue)
}

deserializeSampleSelectionAndOrder <- function(groupSampleDataFrameFieldValue){
  ## unbox
  groupSampleDataFrameName <- "SampleSelectionAndOrder"
  groupSampleDataFrameValue <- substr(
    x = groupSampleDataFrameFieldValue, 
    start = nchar(paste(groupSampleDataFrameName, "={", sep = "")) + 1, 
    stop = nchar(groupSampleDataFrameFieldValue) - nchar("}")
  )
  
  ## unwrap
  columnsSerialized <- strsplit(x = groupSampleDataFrameValue, split = "\\|")[[1]]
  columnNames = unlist(lapply(X = strsplit(x = columnsSerialized, split = "="), FUN = function(x){
    x[[1]]
  }))
  groupSampleDataFrame <- as.data.frame(stringsAsFactors = FALSE, lapply(X = strsplit(x = columnsSerialized, split = "="), FUN = function(x){
    cellContent <- x[[2]]
    cellContent <- substr(
      x = cellContent, 
      start = 1 + nchar("("), 
      stop = nchar(cellContent) - nchar(")")
    )
    cellContent <- strsplit(x = cellContent, split = "; ")
  }))
  colnames(groupSampleDataFrame) <- columnNames
  
  groupSampleDataFrame[, "Order"]   <- as.integer(groupSampleDataFrame[, "Order"])
  groupSampleDataFrame[, "Exclude"] <- as.logical(groupSampleDataFrame[, "Exclude"])
  
  return(groupSampleDataFrame)
}


#' @export
serializeParameterSetFile <- function(importParameterSet, toolName, toolVersion){
  ## wrap
  importParametersValue <- paste(names(importParameterSet), importParameterSet, sep = "=", collapse = "\n")
  ## box
  comment <- paste(
    "# This is the set of parameters which have been used for the initial data import to ",
    importParameterSet$toolVersion,
    "\n",
    "# Exported with ",
    toolName, " ", toolVersion,
    sep = ""
  )
  importParametersFileValue <- paste(comment, importParametersValue, sep = "\n")
  return(importParametersFileValue)
}


#' @export
deserializeParameterSetFile <- function(importParametersFileContent){
  ## remove comments
  importParametersValuePairs <- importParametersFileContent[-grep(pattern = "#.*", x = importParametersFileContent)]
  ## deserialize
  importParameterSet <- deserializeParameterSetKeyValuePairs(importParametersValuePairs)
  return(importParameterSet)
}

serializeParameterSet <- function(importParameterSet){
  ## wrap
  importParametersValue <- paste(names(importParameterSet), importParameterSet, sep = "=", collapse = "; ")
  ## box
  importParametersName <- "ImportParameters"
  importParametersFieldValue <- paste(importParametersName, "={", importParametersValue, "}", sep = "")
  return(importParametersFieldValue)
}

deserializeParameterSet <- function(importParametersFieldValue){
  ## unbox
  importParametersName <- "ImportParameters"
  importParametersValue <- substr(
    x = importParametersFieldValue, 
    start = nchar(paste(importParametersName, "={", sep = "")) + 1, 
    stop = nchar(importParametersFieldValue) - nchar("}")
  )
  
  ## unwrap
  importParametersValuePairs <- unlist(strsplit(x = importParametersValue, split = "; "))
  importParameterSet <- deserializeParameterSetKeyValuePairs(importParametersValuePairs)
  return(importParameterSet)
}

deserializeParameterSetKeyValuePairs <- function(importParametersValuePairs){
  ## unwrap
  importParametersValuePairsList <- strsplit(x = importParametersValuePairs, split = "=")
  
  ## catch empty parameter values
  for(i in seq_len(length(importParametersValuePairsList)))
    if(length(importParametersValuePairsList[[i]]) == 1)
      importParametersValuePairsList[[i]][[2]] <- ""
  
  ## split
  importParametersValues <- unlist(importParametersValuePairsList)
  importParametersKeys   <- importParametersValues[seq(from = 1, to = length(importParametersValues), by = 2)]
  importParametersValues <- importParametersValues[seq(from = 2, to = length(importParametersValues), by = 2)]
  
  ## box to list
  importParameterSet        <- as.list(importParametersValues)
  names(importParameterSet) <- importParametersKeys
  
  ## cast logical's and numeric's
  importParameterSet <- castListEntries(importParameterSet)
  return(importParameterSet)
}

#' Cast logical's and numeric's in a list or data.frame
#'
#' Tries to cast a list entry (or column in a data.frame) to logical's, 
#' if that does not create any missing values, it is assumed 
#' to be a logical will be replaced by `as.logical()` conversion.
#' Similarly for numeric entries (or columns). Everything else remains strings
#' 
#' @param list list object
#'
#' @return list of the same lenght with logical's and numeric's casted
#' @export
castListEntries <- function(list){
  ## cast logical's and numeric's
  suppressWarnings(
    for(idx in seq_len(length(list))){
      if(!is.na(as.logical(list[[idx]]))){
        ## logical
        list[[idx]] <- as.logical(list[[idx]])
      } else if(!is.na(as.numeric(list[[idx]]))){
        ## numeric
        list[[idx]] <- as.numeric(list[[idx]])
      } else {
        ## string
      }
    }
  )
  
  return(list)
}

#########################################################################################
## annotation stuff
precursorSetToSetOfAnnotationSets <- function(dataList, precursorSet){
  setOfAnnotationSets <- lapply(X = precursorSet, FUN = function(x){
    annotationSet <- dataList$annoArrayOfLists[[x]]
    if(dataList$annoArrayIsArtifact[[x]]) {
      annotationSet <- c(annotationSet, "Ignore")
    }
    return(unlist(annotationSet))
  })
  return(setOfAnnotationSets)
}

#' Find color for annotations
#'
#' @param dataList dataList object
#' @param setOfAnnotationSets list of features with annotation class
#'
#' @returns list of annotations and colors
#' @noRd
setOfAnnotationSetsToSetOfColorSets <- function(dataList, setOfAnnotationSets){
  # setOfAnnotationSets %>% unlist %>% table  
  
  makeColorSet <- function(x){
    if(is.null(x)) {
      ## no annotation
      colors <- "black"
    } else {
      ## at least one annotation
      matchColor <- function(y) {
        dataList$annoPresentColorsList[[match(x = y, table = dataList$annoPresentAnnotationsList)]]
      }
      colors <- unlist(lapply(x, matchColor))
    }
    return(colors)
  }
  
  setOfColorSets <- lapply(X = setOfAnnotationSets, FUN = makeColorSet)
  return(setOfColorSets)
}


#' Create annotation colors
#'
#' @param dataList list object
#' @param precursorSet vector
#'
#' @returns list with colors and annotations
#' @export
getPrecursorColors <- function(dataList, precursorSet){
  setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
  setOfColorSets      <- setOfAnnotationSetsToSetOfColorSets(dataList, setOfAnnotationSets)
  setOfColors <- lapply(X = setOfColorSets, FUN = function(x){
    if(any(x == "red")) {
      ## at least one annotation is ignore --> take ignore
      color <- "red"
    } else {
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
    if(any(x == "Ignore")) {
      ## at least one annotation is ignore --> take ignore
      annotation <- "Ignore"
    } else {
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

#########################################################################################
## data fetching

#' @export
getMetFragLink <- function(dataList, precursorIndex, outAdductWarning = TRUE){
  features <- dataList$featureIndeces[[precursorIndex]]
  fragmentsX <- dataList$fragmentMasses[features]
  fragmentsY <- as.numeric(dataList$featureMatrix[precursorIndex, features])
  fragmentsY[fragmentsY > 1] <- 1
  
  precursorMass  <- as.numeric(dataList$dataFrameInfos$"m/z"[[precursorIndex]])
  adduct <- "Unknown"
  if("Adduct ion name" %in% colnames(dataList$dataFrameInfos))
    adduct <- dataList$dataFrameInfos$"Adduct ion name"[[precursorIndex]]
  if("Adduct.ion.name" %in% colnames(dataList$dataFrameInfos))
    adduct <- dataList$dataFrameInfos$"Adduct.ion.name"[[precursorIndex]]
  neutralMassCorrection <- NA
  ionMode <- NA
  #generateLink <- TRUE
  
  error <- NULL
  
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
           error <- paste("This MS\u00B9 feature cannot be send to MetFrag, because the adduct is unknown.")
         },{
           #stop(paste("Unknown adduct (", adduct, ")!", sep = ""))
           if(outAdductWarning) print(paste("###### Unknown adduct (", adduct, ")!", sep = ""))
           #generateLink <- FALSE
           neutralMassCorrection <- NA
           ionMode <- NA
           error <- paste("This MS\u00B9 feature cannot be send to MetFrag, because the adduct '", adduct, "' is not supported.")
         }
  )
  neutralMass <- precursorMass + neutralMassCorrection
  
  fragmentsPositive <- fragmentsX > 0
  fragmentsPositiveX <- fragmentsX[fragmentsPositive]
  fragmentsPositiveY <- fragmentsY[fragmentsPositive]
  fragmentStrings <- paste(fragmentsPositiveX, fragmentsPositiveY, sep = " ", collapse = "; ")
  
  metFragOld <- FALSE
  if(metFragOld){
    # http://msbi.ipb-halle.de/MetFragBeta/LandingPage.jspx?limit=1000&ionmode=-1&database=pubchem&mzppm=7&mzabs=0.005&mass=448.468&formula=C16H20N2O9S2&mzabs=0.05&peaks=130.0655 288214.8119 ; 207.0589 422771.0127 ; 208.0622  87002.3217 ; 210.1334   2674.1707 ; 351.1016  27580.9393 ; 369.1115 739357.5045 ; 370.1148 143864.9611 ; 385.1094   5971.8328 ; 391.0937 337133.4536 ; 392.1025  40126.6888 ; 407.0678   3095.0322 ; 449.0690  37952.2515 
    landingPageUrl <- paste(sep = "",
                            "http://msbi.ipb-halle.de/MetFrag/LandingPage.jspx?",
                            "mass=", neutralMass, "&",
                            "formula=", "", "&",
                            "ionmode=", ionMode, "&",
                            #"limit=", "1000", "&",
                            "database=", "pubchem", "&",
                            #"mzppm=", "7", "&"
                            #"mzabs=", "0.005", "&",
                            "peaks=", fragmentStrings
    )
  } else {
    
    fragmentStrings <- paste(fragmentsPositiveX, fragmentsPositiveY, sep = "_", collapse = ";")
    
    ## https://msbi.ipb-halle.de/MetFragBeta/landing.xhtml?FragmentPeakMatchAbsoluteMassDeviation=0.01&FragmentPeakMatchRelativeMassDeviation=10&DatabaseSearchRelativeMassDeviation=10&PeakList=110_100;210_100&IonizedPrecursorMass=200.101&MetFragDatabaseType=PubChem
    #FragmentPeakMatchAbsoluteMassDeviation
    #FragmentPeakMatchRelativeMassDeviation
    #DatabaseSearchRelativeMassDeviation
    #PrecursorCompoundIDs
    #IonizedPrecursorMass
    #NeutralPrecursorMolecularFormula
    #PrecursorIonMode
    #IonizedPrecursorMass
    landingPageUrl <- paste(sep = "",
                            "https://msbi.ipb-halle.de/MetFragBeta/landing.xhtml", "?",
                            #"https://msbi.ipb-halle.de/MetFrag/landing.xhtml", "?",
                            "NeutralPrecursorMass", "=", neutralMass, "&",
                            "PrecursorIonMode", "=", ionMode, "&",
                            "MetFragDatabaseType", "=", "PubChem", "&",
                            "PeakList", "=", fragmentStrings
    )
  }
  
  
  #writeClipboard(landingPageUrl, format = 1)
  returObj <- list(
    error = error,
    landingPageUrl = landingPageUrl,
    precursorMass = precursorMass,
    neutralMass = neutralMass,
    adduct = adduct,
    fragmentMass = fragmentsPositiveX,
    fragmentIntensities = fragmentsPositiveY
  )
  
  return(returObj)
}

#' @export
getMS2spectrum <- function(dataList, clusterDataList, treeLabel){
  if(treeLabel < 0){
    ###############################################
    ## leaf
    #return(getMS2spectrumInfoForPrecursorLeaf(dataList, clusterDataList, treeLabel))
    return(clusterDataList$ms2spectrumInfoForLeaves[[-treeLabel]])
  } else {
    ###############################################
    ## inner node
    #return(getMS2spectrumInfoForCluster(dataList, clusterDataList, treeLabel))
    return(clusterDataList$ms2spectrumInfoForClusters[[treeLabel]])
  }
}

getMS2spectrumInfoForPrecursorLeaf <- function(dataList, clusterDataList, treeLabel, outAdductWarning = TRUE){
  if(treeLabel >= 0)
    return(NULL)
  ###############################################
  ## leaf
  precursorIndex <- clusterDataList$filterObj$filter[[-treeLabel]]
  return(getMS2spectrumInfoForPrecursor(dataList, precursorIndex, outAdductWarning))
}

#' @export
getMS2spectrumInfoForPrecursor <- function(dataList, precursorIndex, outAdductWarning = TRUE){
  
  precursorSet <- precursorIndex
  numberOfPrecursors <- length(precursorSet)
  
  ## fragments
  features <- dataList$featureIndeces[[precursorIndex]]
  fragmentsX <- dataList$fragmentMasses[features]
  fragmentsY <- as.numeric(dataList$featureMatrix[precursorIndex, features])
  fragmentsY[fragmentsY > 1] <- 1
  fragmentsColor <- rep(x = "black", times = length(fragmentsY))
  
  ## fragment discriminativity
  fragmentDiscriminativity <- rep(x = 0, times = length(features))
  
  
  ## info and MetFrag link
  featureID <- trimws(gsub(x = dataList$precursorLabels[[precursorIndex]], pattern = " +", replacement = " "))
  #featureID <- trimws(gsub(x = clusterDataList$cluster$labels[[-treeLabel]], pattern = " +", replacement = " "))
  featureFamilies <- dataList$annoArrayOfLists[[precursorIndex]]
  featureFamilies <- ifelse(
    test = length(featureFamilies) == 0, 
    yes = "None", 
    no = paste(unlist(featureFamilies), collapse = ", ")
  )
  featureName <- dataList$dataFrameInfos[[precursorIndex, "Metabolite name"]]
  
  infoText <- paste(
    "The MS/MS spectrum of MS\u00B9 feature '", 
    featureID, "'", 
    " comprises ", length(fragmentsX), " fragments. Families: ", featureFamilies, ". Name: ", featureName,
    sep = ""
  )
  metFragLinkList <- getMetFragLink(dataList, precursorIndex, outAdductWarning = FALSE)
  
  
  ## order data
  order <- order(fragmentsX)
  fragmentsX <- fragmentsX[order]
  fragmentsY <- fragmentsY[order]
  fragmentsColor <- fragmentsColor[order]
  fragmentDiscriminativity <- fragmentDiscriminativity[order]
  
  ## box
  resultObj <- list()
  resultObj$fragmentMasses <- fragmentsX
  resultObj$fragmentAbundances <- fragmentsY
  resultObj$fragmentColor <- fragmentsColor
  resultObj$fragmentDiscriminativity <- fragmentDiscriminativity
  resultObj$infoText <- infoText
  resultObj$infoFeatureLabel <- featureID
  resultObj$infoFragmentCount <- length(fragmentsX)
  resultObj$infoFamilies <- featureFamilies
  resultObj$infoName <- featureName
  resultObj$metFragLinkList <- metFragLinkList
  resultObj$precursorSet <- precursorSet
  resultObj$numberOfPrecursors <- numberOfPrecursors
  
  return(resultObj)
}
getMS2spectrumInfoForCluster <- function(dataList, clusterDataList, treeLabel,
                                         minimumProportionOfLeafs = 0.75, minimumProportionToShowFragment = 0.5){
  if(treeLabel < 0)
    return(NULL)
  ###############################################
  ## inner node
  clusterIndex <- treeLabel
  clusterMembersPrecursors <- sort(clusterDataList$innerNodeMembersPrecursors[[clusterIndex]])
  precursorSet <- clusterMembersPrecursors
  numberOfPrecursors <- length(precursorSet)
  numberOfClusterMembers <- length(clusterMembersPrecursors)
  
  ## fragments
  featuresIntersection <- clusterDataList$innerNodeFeaturesIntersection[[clusterIndex]]
  featuresUnion <- clusterDataList$innerNodeFeaturesUnion[[clusterIndex]]
  #fragmentsX <- dataList$fragmentMasses[featuresIntersection]
  #fragmentsY <- apply(X = data.numericmatrix(dataList$featureMatrix[clusterMembersPrecursors, featuresIntersection]), MARGIN = 2, FUN = mean)
  fragmentsX <- dataList$fragmentMasses[featuresUnion]
  #fragmentsY <- apply(X = data.numericmatrix(dataList$featureMatrix[clusterMembersPrecursors, featuresUnion]), MARGIN = 2, FUN = mean)
  fragmentsY <- apply( X = data.matrix(dataList$featureMatrix[clusterMembersPrecursors, featuresUnion, drop = FALSE]), MARGIN = 2, FUN = mean )
  
  selectedPositive <- clusterDataList$innerNodeFeaturesCountsMatrix[clusterIndex, featuresUnion]
  coverageSelected <- selectedPositive / numberOfClusterMembers
  #fragmentsColor <- rep(x = "black", times = length(fragmentsY))
  fragmentsColor <- vector(length = length(fragmentsX))
  fragmentsColor[coverageSelected >= minimumProportionOfLeafs] <- "black"
  fragmentsColor[coverageSelected < minimumProportionOfLeafs] <- "grey"
  
  ## fragment discriminativity
  rootIndex <- length(clusterDataList$cluster$height)
  numberOfLeaves <- length(clusterDataList$innerNodeMembersPrecursors[[rootIndex]])
  numberOfNotClusterMembers <- numberOfLeaves - numberOfClusterMembers
  positives <- clusterDataList$innerNodeFeaturesCountsMatrix[rootIndex, featuresUnion]
  notSelectedPositive <- positives - selectedPositive
  
  #selectedNegative <- numberOfClusterMembers - featuresCountsMatrixSelected
  #notSelectedNegative <- numberOfNotClusterMembers - featuresCountsMatrixNotSelected
  #coverageNotSelected <- notSelectedPositive / numberOfNotClusterMembers
  #coverageAll <- positives / numberOfLeaves
  #relativePositives <- (selectedPositive - notSelectedPositive) / numberOfClusterMembers
  #fragmentDiscriminativity <- coverageSelected - coverageNotSelected
  #fragmentDiscriminativity <- coverageSelected * (1 - coverageNotSelected)
  #fragmentDiscriminativity <- ((selectedPositive - notSelectedPositive) / numberOfClusterMembers) * coverageSelected
  fragmentDiscriminativity <- (selectedPositive - notSelectedPositive) / numberOfClusterMembers
  fragmentDiscriminativity[fragmentDiscriminativity < 0] <- 0
  
  ## reduce to fragments above minimumProportionToShowFragment
  fragmentsX               <- fragmentsX              [coverageSelected > minimumProportionToShowFragment]
  fragmentsY               <- fragmentsY              [coverageSelected > minimumProportionToShowFragment]
  fragmentsColor           <- fragmentsColor          [coverageSelected > minimumProportionToShowFragment]
  fragmentDiscriminativity <- fragmentDiscriminativity[coverageSelected > minimumProportionToShowFragment]
  
  if(length(fragmentDiscriminativity) > 0)
    clusterDiscriminativity <- max(fragmentDiscriminativity)
  else
    clusterDiscriminativity <- 0
  
  ## annotations
  minimumProportionOfMembership <- 0.5
  featureFamilies_all <- unlist(unique(dataList$annoArrayOfLists[precursorSet])) ## all families
  proportionOfMembership <- sapply(X = featureFamilies_all, FUN = function(featureFamily){
    numberOfMembersHere <- sum(unlist(lapply(X = dataList$annoArrayOfLists[precursorSet], function(families){featureFamily %in% families})))
    return(numberOfMembersHere / length(precursorSet))
  })
  frequentFamilies <- featureFamilies_all[proportionOfMembership >= minimumProportionOfMembership]
  
  featureFamilies <- unlist(unique(dataList$annoArrayOfLists[precursorSet]))
  featureFamilies <- ifelse(
    test = length(featureFamilies) == 0, 
    yes = "None", 
    no = paste(unlist(featureFamilies), collapse = ", ")
  )
  
  ## info
  infoText <- paste(
    "This cluster has a cluster discriminativity of ", format(x = clusterDiscriminativity*100, digits = 3, nsmall = 2), "%",
    " and comprises ", length(clusterMembersPrecursors), " MS\u00B9 features",
    " which have ", length(fragmentsX), " fragment(s) in common.", 
    sep = ""
  )
  
  ## order data
  order <- order(fragmentsX)
  fragmentsX <- fragmentsX[order]
  fragmentsY <- fragmentsY[order]
  fragmentsColor <- fragmentsColor[order]
  fragmentDiscriminativity <- fragmentDiscriminativity[order]
  
  ## box
  resultObj <- list()
  resultObj$fragmentMasses <- fragmentsX
  resultObj$fragmentAbundances <- fragmentsY
  resultObj$fragmentColor <- fragmentsColor
  resultObj$fragmentDiscriminativity <- fragmentDiscriminativity
  resultObj$clusterDiscriminativity <- clusterDiscriminativity
  resultObj$frequentFamilies <- frequentFamilies
  resultObj$infoText <- infoText
  resultObj$metFragLinkList <- NULL
  resultObj$precursorSet <- precursorSet
  resultObj$numberOfPrecursors <- numberOfPrecursors
  
  return(resultObj)
}
## XXX adapt getTableFromTreeSelection

#' @export
getTableFromPrecursorSet <- function(dataList, precursorSet, minimumProportionToShowFragment = 0.5){
  ###############################################
  ## table data
  numberOfPrecursors <- length(precursorSet)
  
  ## measurements
  columnNames <- unlist(lapply(X = dataList$sampleClasses, FUN = dataList$dataMeanColumnNameFunctionFromName))
  dataFrameMeasurements     <- data.frame(dataList$dataFrameMeasurements[precursorSet, columnNames, drop=FALSE])
  colnames(dataFrameMeasurements) <- columnNames
  rownames(dataFrameMeasurements) <- dataList$precursorLabels[precursorSet]
  dataFrameMeasurements <- format(x = dataFrameMeasurements, nsmall = 4)
  
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
  featureMatrix <- format(x = featureMatrix, digits = 1, nsmall = 4)
  if(length(featureIndeces) > 0){
    featureMatrix[featureMatrix=="0.0000"] <- "-"
    featureMatrix[featureMatrix=="1.0000"] <- "1"
  }
  
  ## annotations
  setOfAnnotationSets <- precursorSetToSetOfAnnotationSets(dataList, precursorSet)
  setOfAnnotations <- unlist(lapply(X = setOfAnnotationSets, FUN = function(x){
    paste(x, collapse = ", ")
  }))
  annotationDataFrame <- data.frame("Annotation" = setOfAnnotations, row.names = rownames(dataFrameMeasurements))
  
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

#' @export
getPrecursorSetFromTreeSelections <- function(clusterDataList, clusterLabels){
  precursorSet <- NULL
  for(clusterLabel in clusterLabels)
    precursorSet <- c(precursorSet, getPrecursorSetFromTreeSelection(clusterDataList, clusterLabel))
  return(precursorSet)
}

#' @export
getPrecursorSetFromTreeSelection <- function(clusterDataList, clusterLabel){
  if(clusterLabel < 0){
    ###############################################
    ## leaf
    precursorIndex <- clusterDataList$filterObj$filter[[-clusterLabel]]
    precursorSet <- precursorIndex
  } else {
    ###############################################
    ## inner node
    clusterIndex <- clusterLabel
    precursorSet <- sort(clusterDataList$innerNodeMembersPrecursors[[clusterIndex]])
  }
  return(precursorSet)
}

#' @export
getSpectrumStatistics <- function(dataList, precursorSet){
  
  fragmentCounts <- Matrix::colSums(x = dataList$featureMatrix[precursorSet, , drop=FALSE] != 0)
  theseFragments <- fragmentCounts > 0
  fragmentCounts <- fragmentCounts[theseFragments]
  fragmentMasses <- dataList$ms2_masses[theseFragments]
  return(list(
    fragmentMasses = fragmentMasses,
    fragmentCounts = fragmentCounts
  ))
}
getMS2plotData <- function(matrixRows, matrixCols, matrixVals, fragmentMasses){
  numberOfFragments <- length(fragmentMasses)
  meanIntensity <- vector(mode = "numeric", length = numberOfFragments)
  fragmentCount <- vector(mode = "numeric", length = numberOfFragments)
  for(colIdx in seq_len(numberOfFragments)){
    intensities <- matrixVals[matrixCols == colIdx]
    fragmentCount[[colIdx]] <- length(intensities)
    meanIntensity[[colIdx]] <- mean(x = intensities)
  }
  
  presentFragments <- fragmentCount > 0
  
  resultObj <- list()
  resultObj$numberOfFragments <- fragmentCount
  resultObj$averageAbundance  <- meanIntensity
  resultObj$masses            <- fragmentMasses
  
  return(resultObj)
}
regExExtraction = function(pattern, x, ...) {
  args = list(...)
  args[['perl']] = T
  
  re = do.call(gregexpr, c(list(pattern, x), args))
  
  mapply(function(re, x){
    
    cap = sapply(attr(re, 'capture.names'), function(n, re, x){
      start = attr(re, 'capture.start')[, n]
      len   = attr(re, 'capture.length')[, n]
      end   = start + len - 1
      tok   = substr(rep(x, length(start)), start, end)
      
      return(tok)
    }, re, x, simplify=F, USE.NAMES=T)
    
    return(cap)
  }, re, x, SIMPLIFY=F)
}

#' @export
numericVectorToStringForEval <- function(vec){
  return(paste("c(", paste(vec, collapse = ","), ")", sep = ""))
}

#' @export
colorVectorToStringForEval <- function(vec){
  return(paste("c('", paste(vec, collapse = "','"), "')", sep = ""))
}

#' Swap a vector's names and values
#' 
#' Similar to `searchable::invert()`
#'
#' @param x vector
#'
#' @returns string
swap_names_and_values <- function(x) {
  
  if( is.null( names(x) ) ) stop( "vector does not have names.")
  
  v <- names(x)
  names(v) <- as.character(x)
  
  return(v)
  
}
