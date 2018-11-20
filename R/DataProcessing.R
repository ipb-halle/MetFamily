
#########################################################################################
## annotate and process matrix
sparseMatrixToString <- function(matrixRows, matrixCols, matrixVals, parameterSet){
  matrixRows <- c(matrixRows, 1)
  matrixCols <- c(matrixCols, 1)
  matrixVals <- c(matrixVals, serializeParameterSet(parameterSet))
  
  ## TODO performance 25s
  ## convert matrix to dataframe
  numberOfRows    <- max(matrixRows)
  numberOfColumns <- max(matrixCols)
  
  lines <- vector(mode = "character", length = numberOfRows)
  for(rowIdx in seq_len(numberOfRows)){
    indeces <- matrixRows == rowIdx
    tokens  <- vector(mode = "character", length = numberOfColumns)
    tokens[matrixCols[indeces]] <- matrixVals[indeces]
    lines[[rowIdx]] <- paste(tokens, collapse = "\t")
  }
  
  return(lines)
}

readClusterDataFromProjectFile <- function(file, progress = FALSE){
  if(progress)  setProgress(value = 0, detail = "Parsing") else print("Parsing")
  extension <- file_ext(file)
  if(extension == "gz"){
    file <- gzfile(file, "r")
  } else {
    file <- file(file, "r")
  }
  
  suppressWarnings(
    fileLines <- readLines(con = file)
  )
  base::close(con = file)
  
  dataList <- readProjectData(fileLines = fileLines, progress = progress)
  fileLines <- NULL
  
  return(dataList)
}
readProjectData <- function(fileLines, progress = FALSE){
  allowedTags <- c("ID")
  allowedTagPrefixes <- c("AnnotationColors=")
  
  ##################################################################################################
  ## parse data
  if(progress)  incProgress(amount = 0.1, detail = "Preprocessing") else print("Preprocessing")
  
  numberOfRows <- length(fileLines)
  numberOfMS1features <- as.integer(numberOfRows - 3)
  
  ## header
  line1Tokens <- strsplit(x = fileLines[[1]], split = "\t")[[1]]
  line2Tokens <- strsplit(x = fileLines[[2]], split = "\t")[[1]]
  line3Tokens <- strsplit(x = fileLines[[3]], split = "\t")[[1]]
  
  ## metabolite profile vs fragmentMatrix
  numberOfColumns <- length(line1Tokens)
  line1TokensOffset <- 2
  fragmentMatrixStart <- min(which(line1Tokens[(1 + line1TokensOffset):numberOfColumns] != "")) + line1TokensOffset
  numberOfMetaboliteProfileColumns <- fragmentMatrixStart - 1
  numberOfFragmentGroups <- numberOfColumns - numberOfMetaboliteProfileColumns
  
  ## extract infos from header
  importParameters <- line1Tokens[[1]]
  if(nchar(importParameters) == 0){
    ## import parameterSet not there: backward compatibility - add if not there
    importParameters <- "ImportParameters={projectName=MetFamily project; projectDescription=; toolVersion=MetFamily 1.0; minimumIntensityOfMaximalMS2peak=2000; minimumProportionOfMS2peaks=0.05; mzDeviationAbsolute_grouping=0.01; mzDeviationInPPM_grouping=10; doPrecursorDeisotoping=TRUE; mzDeviationAbsolute_precursorDeisotoping=0.001; mzDeviationInPPM_precursorDeisotoping=10; maximumRtDifference=0.02; doMs2PeakGroupDeisotoping=FALSE; mzDeviationAbsolute_ms2PeakGroupDeisotoping=0.01; mzDeviationInPPM_ms2PeakGroupDeisotoping=10; proportionOfMatchingPeaks_ms2PeakGroupDeisotoping=0.9; mzDeviationAbsolute_mapping=0.01; minimumNumberOfMS2PeaksPerGroup=1; neutralLossesPrecursorToFragments=TRUE; neutralLossesFragmentsToFragments=FALSE}"
  }
  groupSampleDataFrameFieldValue <- line1Tokens[[2]]
  
  fragmentGroupsNumberOfFramgents <- as.integer(line1Tokens[fragmentMatrixStart:numberOfColumns])
  line1Tokens <- NULL
  
  tagsSector <- line2Tokens[seq_len(numberOfMetaboliteProfileColumns)]
  fragmentGroupsAverageIntensity <- as.numeric(line2Tokens[fragmentMatrixStart:numberOfColumns])
  line2Tokens <- NULL
  
  metaboliteProfileColumnNames <- line3Tokens[seq_len(numberOfMetaboliteProfileColumns)]
  fragmentGroupsAverageMass <- as.numeric(line3Tokens[fragmentMatrixStart:numberOfColumns])
  line3Tokens <- NULL
  
  if(any(duplicated(metaboliteProfileColumnNames)))
    stop(paste("Duplicated column names in the metabolite profile: ", paste(sort(unique(metaboliteProfileColumnNames[duplicated(metaboliteProfileColumnNames)])), collapse = "; ")))
  
  #########################################################################
  ## extract metabolite profile and fragment matrix
  metaboliteProfile <- as.data.frame(matrix(nrow = numberOfMS1features, ncol = numberOfMetaboliteProfileColumns))
  colnames(metaboliteProfile) <- metaboliteProfileColumnNames
  
  listMatrixRows <- list()
  listMatrixCols <- list()
  listMatrixVals <- list()
  
  lastOut <- proc.time()["user.self"]
  lastRow <- 1
  for(rowIdx in seq_len(numberOfMS1features)){
    time <- proc.time()["user.self"]
    if(time - lastOut > 1){
      lastOut <- time
      rowProgress <- (rowIdx - lastRow) / numberOfMS1features
      lastRow <- rowIdx
      if(progress)  incProgress(amount = rowProgress*0.2,     detail = paste("Preprocessing ", rowIdx, " / ", numberOfMS1features, sep = "")) else print(paste("Preprocessing ", rowIdx, " / ", numberOfMS1features, sep = ""))
    }
    
    lineIdx <- rowIdx + 3
    tokens <- str_split(string = fileLines[[lineIdx]], pattern = "\t")[[1]]
    
    ## metabolite profile
    metaboliteProfile[rowIdx, ] <- tokens[seq_len(numberOfMetaboliteProfileColumns)]
    
    ## fragment matrix
    tokens <- tokens[fragmentMatrixStart:numberOfColumns]
    nonEmpty <- tokens != ""
    indeces <- which(nonEmpty)
    numberOfEntries <- length(indeces)
    listMatrixRows[[rowIdx]] <- rep(x = rowIdx, times = numberOfEntries)
    listMatrixCols[[rowIdx]] <- indeces
    listMatrixVals[[rowIdx]] <- tokens[nonEmpty]
  }
  matrixRows <- as.integer(unlist(listMatrixRows))
  matrixCols <- as.integer(unlist(listMatrixCols))
  matrixVals <- as.numeric(unlist(listMatrixVals))
  listMatrixRows <- NULL
  listMatrixCols <- NULL
  listMatrixVals <- NULL
  
  ## header
  dataFrameHeader <- cbind(
    data.frame(rbind(
      c(importParameters, rep(x = "", times = numberOfMetaboliteProfileColumns - 1)),
      tagsSector,
      metaboliteProfileColumnNames
    ), stringsAsFactors = FALSE),
    data.frame(rbind(
      fragmentGroupsNumberOfFramgents,
      fragmentGroupsAverageIntensity,
      fragmentGroupsAverageMass
    ), stringsAsFactors = FALSE)
  )
  headerLabels <- c("HeaderForFragmentCounts", "HeaderForGroupsAndFragmentIntensities", "Header")
  rownames(dataFrameHeader) <- headerLabels
  headerColumnNames <- c(metaboliteProfileColumnNames, fragmentGroupsAverageMass)
  colnames(dataFrameHeader) <- headerColumnNames
  
  ## import parameterSet
  importParameterSet <- deserializeParameterSet(importParameters)
  
  ## insert annotation column if not there
  annotationColorsName <- "AnnotationColors"
  annotationColorsMapInitValue <- paste(annotationColorsName, "={}", sep = "")
  annotationColumnName <- "Annotation"
  if(!any(metaboliteProfileColumnNames == annotationColumnName, na.rm = TRUE)){
    ## annotation column backward compatibility - insert if not there
    target <- 2
    
    if(target == 0 | target == numberOfMetaboliteProfileColumns)
      stop("Cannot insert column!")
    
    metaboliteProfile <- cbind(
      metaboliteProfile[,seq_len(target),drop=F], 
      as.data.frame(x = rep(x = "", times = numberOfMS1features), stringsAsFactors = FALSE), 
      metaboliteProfile[, (target+1):numberOfMetaboliteProfileColumns, drop=FALSE]
    )
    dataFrameHeader <- cbind(
      dataFrameHeader[,seq_len(target),drop=F], 
      as.data.frame(x = rep(x = "", times = numberOfMS1features), stringsAsFactors = FALSE), 
      dataFrameHeader[, (target+1):numberOfColumns, drop=FALSE]
    )
    numberOfMetaboliteProfileColumns <- numberOfMetaboliteProfileColumns + 1
    metaboliteProfileColumnNames <- c(metaboliteProfileColumnNames[seq_len(target)], annotationColumnName, metaboliteProfileColumnNames[(target+1):numberOfMetaboliteProfileColumns])
    colnames(metaboliteProfile) <- metaboliteProfileColumnNames
    headerColumnNames <- c(metaboliteProfileColumnNames, fragmentGroupsAverageMass)
    colnames(dataFrameHeader) <- headerColumnNames
    
    dataFrameHeader[2, target + 1] <- annotationColorsMapInitValue
    dataFrameHeader[3, target + 1] <- annotationColumnName
  }
  
  annotationColumnIndex <- which(metaboliteProfileColumnNames == annotationColumnName)
  annotationColorsValue <- dataFrameHeader[2, annotationColumnIndex]
  
  dataFrameMS1Header <- dataFrameHeader[, seq_len(numberOfMetaboliteProfileColumns)]
  
  ##################################################################################################
  ## MS1 feature IDs
  
  ## mz/rt is aligned by '.'
  mzs <- metaboliteProfile[, "m/z"]
  rts <- metaboliteProfile[, "RT"]
  
  ## add .0 if necessary
  for(i in seq_len(numberOfMS1features))
    if(length(grep(x = mzs[[i]], pattern = ".*\\..*")) == 0)
      mzs[[i]] <- paste(mzs[[i]], ".0", sep = "")
  for(i in seq_len(numberOfMS1features))
    if(length(grep(x = rts[[i]], pattern = ".*\\..*")) == 0)
      rts[[i]] <- paste(rts[[i]], ".0", sep = "")
  
  regexResult <- gregexpr(pattern = "^(?<before>\\d+)\\.(?<after>\\d+)$", text = mzs, perl = TRUE)
  mzStartsBefore  <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.start")[[1]]}))
  mzStartsAfter   <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.start")[[2]]}))
  mzLengthsBefore <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.length")[[1]]}))
  mzLengthsAfter  <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.length")[[2]]}))
  
  mzMaxBefore <- max(mzLengthsBefore)
  mzMaxAfter  <- max(mzLengthsAfter )
  
  regexResult <- gregexpr(pattern = "^(?<before>\\d+)\\.(?<after>\\d+)$", text = rts, perl = TRUE)
  rtStartsBefore  <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.start")[[1]]}))
  rtStartsAfter   <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.start")[[2]]}))
  rtLengthsBefore <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.length")[[1]]}))
  rtLengthsAfter  <- unlist(lapply(X = regexResult, FUN = function(x){attr(x = x, which = "capture.length")[[2]]}))
  
  rtMaxBefore <- max(rtLengthsBefore)
  rtMaxAfter  <- max(rtLengthsAfter )
  
  
  maximumNumberOfDecimalPlacesForMz <- 3
  maximumNumberOfDecimalPlacesForRt <- 2
  
  for(idx in seq_along(mzs)){
    mzStartBefore  <- mzStartsBefore [[idx]]
    mzStartAfter   <- mzStartsAfter  [[idx]]
    mzLengthBefore <- mzLengthsBefore[[idx]]
    mzLengthAfter  <- mzLengthsAfter [[idx]]
    
    rtStartBefore  <- rtStartsBefore [[idx]]
    rtStartAfter   <- rtStartsAfter  [[idx]]
    rtLengthBefore <- rtLengthsBefore[[idx]]
    rtLengthAfter  <- rtLengthsAfter [[idx]]
    
    mzBefore <- substr(start = mzStartBefore, stop = mzStartBefore + mzLengthBefore - 1, x = mzs[[idx]])
    mzAfter  <- substr(start = mzStartAfter,  stop = mzStartAfter  + mzLengthAfter  - 1, x = mzs[[idx]])
    
    rtBefore <- substr(start = rtStartBefore, stop = rtStartBefore + rtLengthBefore - 1, x = rts[[idx]])
    rtAfter  <- substr(start = rtStartAfter,  stop = rtStartAfter  + rtLengthAfter  - 1, x = rts[[idx]])
    
    if(nchar(mzBefore) < mzMaxBefore) 
      mzBefore <- paste(
        paste(rep(x = "  ", times = mzMaxBefore - nchar(mzBefore)), collapse = ""),
        mzBefore,
        sep = ""
      )
    if(nchar(mzAfter) > maximumNumberOfDecimalPlacesForMz)
      mzAfter <- substr(x = mzAfter, start = 1, stop = maximumNumberOfDecimalPlacesForMz)
    if(nchar(mzAfter) < maximumNumberOfDecimalPlacesForMz)
      mzAfter <- paste(
        mzAfter,
        paste(rep(x = "0", times = maximumNumberOfDecimalPlacesForMz - nchar(mzAfter)), collapse = ""),
        #paste(rep(x = "  ", times = maximumNumberOfDecimalPlacesForMz - nchar(mzAfter)), collapse = ""),
        sep = ""
      )
    
    if(nchar(rtBefore) < rtMaxBefore) 
      rtBefore <- paste(
        paste(rep(x = "  ", times = rtMaxBefore - nchar(rtBefore)), collapse = ""),
        rtBefore,
        sep = ""
      )
    if(nchar(rtAfter) > maximumNumberOfDecimalPlacesForRt)
      rtAfter <- substr(x = rtAfter, start = 1, stop = maximumNumberOfDecimalPlacesForRt)
    if(nchar(rtAfter) < maximumNumberOfDecimalPlacesForRt)
      rtAfter <- paste(
        rtAfter,
        #paste(rep(x = "  ", times = maximumNumberOfDecimalPlacesForRt - nchar(rtAfter)), collapse = ""),
        paste(rep(x = "0", times = maximumNumberOfDecimalPlacesForRt - nchar(rtAfter)), collapse = ""),
        sep = ""
      )
    
    mzs[[idx]] <- paste(mzBefore, mzAfter, sep = ".")
    rts[[idx]] <- paste(rtBefore, rtAfter, sep = ".")
  }
  precursorLabels <- paste(mzs, rts, sep = " / ")
  
  ## remove duplicated MS1 features
  duplicated <- which(duplicated(precursorLabels))
  numberOfDuplicated <- length(duplicated)
  if(numberOfDuplicated > 0){
    precursorLabels <- precursorLabels[-duplicated]
    metaboliteProfile <- metaboliteProfile[-duplicated, ]
    numberOfMS1features <- numberOfMS1features - numberOfDuplicated
    for(duplicatedRowIdxIdx in seq_along(duplicated)){
      duplicatedRowIdx <- duplicated[[duplicatedRowIdxIdx]]
      
      ## remove row from matrix
      indeces1 <- which(matrixRows == duplicatedRowIdx)
      if(length(indeces1) == 0)
        next
      
      matrixRows <- matrixRows[-indeces1]
      matrixCols <- matrixCols[-indeces1]
      matrixVals <- matrixVals[-indeces1]
      
      ## update subsequent matrix rows
      indeces2 <- which(matrixRows > duplicatedRowIdx)
      matrixRows[indeces2] <- matrixRows[indeces2] - 1
      
      ## update subsequent duplicated rows
      indeces3 <- which(duplicated > duplicatedRowIdx)
      duplicated[indeces3] <- duplicated[indeces3] - 1
    }
  }
  rownames(metaboliteProfile) <- precursorLabels
  
  #############################################################################################
  ## process features
  if(progress)  incProgress(amount = 0.1, detail = "Features") else print("Features")
  
  ## get features
  featureIndeces <- list()
  featureCount <- vector(mode = "numeric", length = numberOfMS1features)
  #fragmentMassPresent <- rep(x = FALSE, times = length(fragmentGroupsAverageMass))
  for(i in seq_len(numberOfMS1features)){
    # if(numberOfMS1features >= 10 & ((i %% (as.integer(numberOfMS1features/10))) == 0))
    #   if(progress)  incProgress(amount = 0.3 / 10, detail = paste("Features ", i, " / ", numberOfMS1features, sep = ""))
    ## data
    indecesHere <- which(matrixRows == i)
    featureIndecesHere <- matrixCols[indecesHere]
    numberOfFeatures <- length(featureIndecesHere)
    
    featureIndeces[[i]] <- featureIndecesHere
    featureCount[[i]] <- numberOfFeatures
    #fragmentMassPresent[featureIndecesHere] <- TRUE
  }
  
  if(progress)  incProgress(amount = 0.1, detail = "Feature postprocessing") else print("Feature postprocessing")
  
  ## ms2 plot data
  # resultObj <- getMS2plotData(matrixRows, matrixCols, matrixVals, fragmentMasses = fragmentGroupsAverageMass)
  # ms2PlotDataNumberOfFragments <- resultObj$numberOfFragments
  # ms2PlotDataAverageAbundance  <- resultObj$averageAbundance
  # ms2PlotDataFragmentMasses    <- resultObj$masses
  ms2PlotDataNumberOfFragments <- fragmentGroupsNumberOfFramgents
  ms2PlotDataAverageAbundance  <- fragmentGroupsAverageIntensity
  ms2PlotDataFragmentMasses    <- fragmentGroupsAverageMass
  maxNumberOfFragments <- max(ms2PlotDataNumberOfFragments)
  ms2PlotDataColorMapFragmentData  <- makecmap(
    x = c(0, maxNumberOfFragments), n = 100, 
    colFn = colorRampPalette(c('grey', 'black'))
  )
  
  ## featureMatrix and annotation
  featureMatrix <- sparseMatrix(i = matrixRows, j = matrixCols, x = matrixVals)
  matrixRows <- NULL
  matrixCols <- NULL
  matrixVals <- NULL
  
  #fragmentGroupsAverageMass <- fragmentGroupsAverageMass[1:ncol(featureMatrix)]
  #fragmentGroupsAverageMass <- fragmentGroupsAverageMass[fragmentMassPresent]
  rownames(featureMatrix) <- precursorLabels
  colnames(featureMatrix) <- fragmentGroupsAverageMass
  
  ## featureIndexMatrix
  featureIndexMatrix <- matrix(nrow = numberOfMS1features, ncol = max(sapply(X = featureIndeces, FUN = length)))
  rownames(featureIndexMatrix) <- precursorLabels
  for(i in seq_len(numberOfMS1features))
    featureIndexMatrix[i, seq_len(length(featureIndeces[[i]]))] <- featureIndeces[[i]]
  
  # ## remove columns without data
  # fragmentThere <- apply(X = featureMatrix, MARGIN = 2, FUN = function(x){any(x != 0)})
  # minimumMass <- min(fragmentGroupsAverageMass[fragmentThere])
  # maximumMass <- max(fragmentGroupsAverageMass[fragmentThere])
  minimumMass <- min(fragmentGroupsAverageMass)
  maximumMass <- max(fragmentGroupsAverageMass)
  
  ##################################################################################################
  ## process sample measurements
  
  ## sample columns
  sampleColumns <- tagsSector != ""
  for(allowedTag in allowedTags)
    sampleColumns[grep(x = tagsSector, pattern = paste("^", allowedTag, "$", sep = ""))] <- FALSE
  for(allowedTagPrefix in allowedTagPrefixes)
    sampleColumns[grep(x = tagsSector, pattern = paste("^", allowedTagPrefix, sep = ""))] <- FALSE
  sampleColumns <- which(sampleColumns)
  sampleColumnsStartEnd <- c(min(sampleColumns), max(sampleColumns))
  
  groups <- unique(tagsSector[sampleColumns])
  numberOfGroups <- length(groups)
  
  sampleNamesToExclude <- NULL
  
  
  dataColumnIndecesFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude = NULL){
    which(tagsSector == groups[[groupIdx]] & !(metaboliteProfileColumnNames %in% sampleNamesToExclude))
  }
  dataColumnsNameFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude = NULL){
    #sampleNames = paste(groups[[groupIdx]], "_", metaboliteProfileColumnNames[dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)], sep = "")
    sampleNames = metaboliteProfileColumnNames[dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)]
    #sampleNames = sampleNames[!(sampleNames %in% sampleNamesToExclude)]
    return(sampleNames)
  }
  dataColumnsNameFunctionFromGroupName <- function(group, sampleNamesToExclude = NULL){
    dataColumnsNameFunctionFromGroupIndex(groupIdx = match(x = group, table = groups), sampleNamesToExclude = sampleNamesToExclude)
  }
  dataColumnsNameFunctionFromGroupNames <- function(groups, sampleNamesToExclude = NULL){
    unlist(lapply(X = groups, FUN = function(x){dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)}))
  }
  groupNameFunctionFromDataColumnName <- function(dataColumnName, sampleNamesToExclude = NULL){
    groupIdx <- which(unlist(lapply(X = groups, FUN = function(x){
      dataColumnNames <- dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)
      any(dataColumnNames == dataColumnName)
    })))
    groups[[groupIdx]]
  }
  lfcColumnNameFunctionFromString <- function(columnName){
    tokens <- strsplit(x = columnName, split = "_vs_")[[1]]
    groupOne <- strsplit(x = tokens[[1]], split = "LFC_")[[1]][[2]]
    groupTwo <- tokens[length(tokens)]
    return(c(groupOne, groupTwo))
  }
  dataMeanColumnNameFunctionFromString <- function(columnName){
    group <- substr(x = columnName, start = 1, stop = nchar(columnName) - nchar("_mean"))
    return(group)
  }
  
  ## manage group and samples: order and exclusion
  groupNames  <- tagsSector[sampleColumns]
  sampleNames <- metaboliteProfileColumnNames[sampleColumns]
  
  if(nchar(groupSampleDataFrameFieldValue) == 0){
    ## not there: backward compatibility - add if not there
    groupSampleDataFrame <- data.frame(stringsAsFactors = FALSE, 
                                       "Group"   = groupNames,
                                       "Sample"  = sampleNames,
                                       "Order"   = seq_along(sampleNames),
                                       "Exclude" = rep(x = FALSE, times = length(sampleNames))
    )
  } else {
    groupSampleDataFrame <- deserializeSampleSelectionAndOrder(groupSampleDataFrameFieldValue)
  }
  dataFrameMS1Header[[1,2]] <- serializeSampleSelectionAndOrder(groupSampleDataFrame)
  
  returnObj <- processMS1data(
    sampleNamesToExclude=sampleNamesToExclude, numberOfMS1features=numberOfMS1features, precursorLabels=precursorLabels, 
    groups=groups, metaboliteProfileColumnNames=metaboliteProfileColumnNames, tagsSector = tagsSector, 
    dataColumnIndecesFunctionFromGroupIndex=dataColumnIndecesFunctionFromGroupIndex, dataColumnsNameFunctionFromGroupIndex=dataColumnsNameFunctionFromGroupIndex, dataColumnsNameFunctionFromGroupName=dataColumnsNameFunctionFromGroupName, dataColumnsNameFunctionFromGroupNames=dataColumnsNameFunctionFromGroupNames, groupNameFunctionFromDataColumnName=groupNameFunctionFromDataColumnName,
    metaboliteProfile=metaboliteProfile, progress=progress
  )
  
  ## name functions
  dataMeanColumnNameFunctionFromIndex <- returnObj$dataMeanColumnNameFunctionFromIndex
  dataMeanColumnNameFunctionFromName  <- returnObj$dataMeanColumnNameFunctionFromName
  lfcColumnNameFunctionFromIndex      <- returnObj$lfcColumnNameFunctionFromIndex
  lfcColumnNameFunctionFromName       <- returnObj$lfcColumnNameFunctionFromName
  groupNameFromGroupIndex             <- returnObj$groupNameFromGroupIndex
  groupIdxFromGroupName               <- returnObj$groupIdxFromGroupName
  ## data and names
  dataFrameMeasurements               <- returnObj$dataFrameMeasurements
  ## colors
  colorMatrixDataFrame                <- returnObj$colorMatrixDataFrame
  colorMapAbsoluteData                <- returnObj$colorMapAbsoluteData
  colorMapLogFoldChange               <- returnObj$colorMapLogFoldChange
  columnGroupLabels                   <- returnObj$columnGroupLabels
  ## constants
  meanAllMax       <- returnObj$meanAllMax
  logFoldChangeMax <- returnObj$logFoldChangeMax
  logAbsMax        <- returnObj$logAbsMax
  
  #########################################################################################
  ## precursor annotation fields
  if(progress)  incProgress(amount = 0.1, detail = "Feature annotations") else print("Feature annotations")
  annotationValueIgnore <- "Ignore"
  annotationColorIgnore <- "red"
  
  ## present annotations
  annotations    <- vector(mode='list', length=numberOfMS1features)
  #annotations[1:numberOfMS1features] <- dataFrame[, annotationColumnName]
  annoVals <- metaboliteProfile[, annotationColumnName]
  for(i in seq_len(numberOfMS1features)){
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
  
  annoArrayOfLists    <- vector(mode='list', length=numberOfMS1features)
  annoArrayIsArtifact <- vector(mode='logical', length=numberOfMS1features)
  for(i in seq_len(numberOfMS1features)){
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
  annotationColorsMapValue <- substr(
    x = annotationColorsValue, 
    start = nchar(paste(annotationColorsName, "={", sep = "")) + 1, 
    stop = nchar(annotationColorsValue) - nchar("}")
  )
  
  if(nchar(annotationColorsMapValue) > 0){
    #annotationColorsMapValuePairsTmp <- unlist(strsplit(x = annotationColorsMapValue, split = "="))
    #annotationColorsMapValues <- sapply(X = strsplit(x = annotationColorsMapValuePairsTmp[2:length(annotationColorsMapValuePairsTmp)], split = ", "), FUN = function(token){
    #  token[[1]]
    #})
    #if(length(annotationColorsMapValuePairsTmp) < 3){
    #  annotationColorsMapKeys <- annotationColorsMapValuePairsTmp[[1]]
    #}else{
    #  annotationColorsMapKeys <- c(annotationColorsMapValuePairsTmp[[1]], substr(
    #    x = annotationColorsMapValuePairsTmp[2:(length(annotationColorsMapValuePairsTmp) - 1)], 
    #    start = nchar(annotationColorsMapValues) + nchar(", ") + 1, 
    #    stop = nchar(annotationColorsMapValuePairsTmp[2:length(annotationColorsMapValuePairsTmp)])
    #  ))
    #}
    
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
    for(i in seq_len(length(annotationColorsMapKeys))){
      annoPresentAnnotationsList[[1 + i]] <- annotationColorsMapKeys[[i]]
      annoPresentColorsList     [[1 + i]] <- annotationColorsMapValues[[i]]
    }
  
  ## check consistency
  if(!all(unique(unlist(annoArrayOfLists)) %in% unlist(annoPresentAnnotationsList))){
    missing <- unique(unlist(annoArrayOfLists))[!(unique(unlist(annoArrayOfLists)) %in% unlist(annoPresentAnnotationsList))]
    stop(paste("Annotation(s)", paste(missing, collapse = "; "), "missing in present annotations list"))
  }
  if(!all(unlist(annoPresentAnnotationsList) %in% unique(c(annotationValueIgnore, unlist(annoArrayOfLists))))){
    missing <- unlist(annoPresentAnnotationsList)[!(unlist(annoPresentAnnotationsList) %in% unique(c(annotationValueIgnore, unlist(annoArrayOfLists))))]
    stop(paste("Present annotation(s)", paste(missing, collapse = "; "), "missing in annotations"))
  }
  
  ##################################################################################################
  ## box
  if(progress)  incProgress(amount = 0.1, detail = "Boxing") else print("Boxing")
  dataList <- list()
  ## data
  dataList$dataFrameHeader <- dataFrameHeader
  dataList$dataFrameMS1Header <- dataFrameMS1Header
  dataList$dataFrameInfos <- metaboliteProfile
  dataList$importParameterSet <- importParameterSet
  dataList$numberOfPrecursors <- numberOfMS1features
  dataList$numberOfDuplicatedPrecursors <- numberOfDuplicated
  dataList$groups <- groups
  dataList$columnGroupLabels <- columnGroupLabels
  dataList$groupSampleDataFrame <- groupSampleDataFrame
  dataList$metaboliteProfileColumnNames <- metaboliteProfileColumnNames
  dataList$tagsSector <- tagsSector
  ## data: fragments
  dataList$fragmentMasses <- fragmentGroupsAverageMass
  dataList$fragmentFrequency <- fragmentGroupsNumberOfFramgents
  dataList$fragmentAbundance <- fragmentGroupsAverageIntensity
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
  ## data: column name functions
  dataList$dataColumnsNameFunctionFromGroupName <- dataColumnsNameFunctionFromGroupName
  dataList$dataColumnsNameFunctionFromGroupIndex <- dataColumnsNameFunctionFromGroupIndex
  dataList$dataColumnsNameFunctionFromGroupNames <- dataColumnsNameFunctionFromGroupNames
  dataList$groupNameFunctionFromDataColumnName <- groupNameFunctionFromDataColumnName
  dataList$lfcColumnNameFunctionFromString <- lfcColumnNameFunctionFromString
  dataList$dataMeanColumnNameFunctionFromString <- dataMeanColumnNameFunctionFromString
  dataList$dataColumnIndecesFunctionFromGroupIndex <- dataColumnIndecesFunctionFromGroupIndex
  dataList$dataMeanColumnNameFunctionFromName <- dataMeanColumnNameFunctionFromName
  dataList$dataMeanColumnNameFunctionFromIndex <- dataMeanColumnNameFunctionFromIndex
  dataList$lfcColumnNameFunctionFromName <- lfcColumnNameFunctionFromName
  dataList$lfcColumnNameFunctionFromIndex <- lfcColumnNameFunctionFromIndex
  dataList$groupNameFromGroupIndex <- groupNameFromGroupIndex
  dataList$groupIdxFromGroupName <- groupIdxFromGroupName
  ## features
  dataList$featureMatrix <- featureMatrix
  dataList$featureIndeces <- featureIndeces
  dataList$featureCount <- featureCount
  dataList$featureIndexMatrix <- featureIndexMatrix
  ## ms2 plot data
  dataList$ms2_numberOfFragments <- ms2PlotDataNumberOfFragments
  dataList$ms2_averageAbundance  <- ms2PlotDataAverageAbundance
  dataList$ms2_masses            <- ms2PlotDataFragmentMasses
  dataList$colorMapFragmentData <- ms2PlotDataColorMapFragmentData
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
  
  if(progress)  setProgress(1) else print("Ready")
  
  
  ## redefine MS1 column functions
  dataColumnIndecesFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude){
    which(dataList$tagsSector == dataList$groups[[groupIdx]] & !(dataList$metaboliteProfileColumnNames %in% sampleNamesToExclude))
  }
  dataList$dataColumnIndecesFunctionFromGroupIndex <- dataColumnIndecesFunctionFromGroupIndex
  dataColumnsNameFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude){
    #sampleNames = paste(dataList$groups[[groupIdx]], "_", metaboliteProfileColumnNames[dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)], sep = "")
    dataList$metaboliteProfileColumnNames[dataList$dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)]
    #sampleNames = sampleNames[!(sampleNames %in% sampleNamesToExclude)]
    #return(sampleNames)
  }
  dataList$dataColumnsNameFunctionFromGroupIndex <- dataColumnsNameFunctionFromGroupIndex
  dataColumnsNameFunctionFromGroupName <- function(group, sampleNamesToExclude){
    dataColumns <- dataList$dataColumnsNameFunctionFromGroupIndex(groupIdx = match(x = group, table = dataList$groups), sampleNamesToExclude = sampleNamesToExclude)
  }
  dataList$dataColumnsNameFunctionFromGroupName <- dataColumnsNameFunctionFromGroupName
  dataColumnsNameFunctionFromGroupNames <- function(groups, sampleNamesToExclude){
    unlist(lapply(X = groups, FUN = function(x){
      dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)
    }))
  }
  dataList$dataColumnsNameFunctionFromGroupNames <- dataColumnsNameFunctionFromGroupNames
  groupNameFunctionFromDataColumnName <- function(dataColumnName, sampleNamesToExclude){
    groupIdx <- which(unlist(lapply(X = dataList$groups, FUN = function(x){
      dataColumnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)
      any(dataColumnNames == dataColumnName)
    })))
    dataList$groups[[groupIdx]]
  }
  dataList$groupNameFunctionFromDataColumnName <- groupNameFunctionFromDataColumnName
  
  orderColumnNames <- function(groupSampleDataFrame, columnNames){
    order <- groupSampleDataFrame[groupSampleDataFrame[, "Sample"] %in% columnNames, "Order"]
    list <- list()
    list[order] <- columnNames
    columnNames <- unlist(list)
    return(columnNames)
  }
  dataList$orderColumnNames <- orderColumnNames
  
  ## define sample in-/exclusion functions
  excludedSamples <- function(groupSampleDataFrame, groups = dataList$groups){
    #dataList$groupSampleDataFrame[, "Sample"][ dataList$groupSampleDataFrame[, "Exclude"]]
    samples    =  groupSampleDataFrame[, "Sample"]
    isExcluded =  groupSampleDataFrame[, "Exclude"]
    isGroup    =  groupSampleDataFrame[, "Group"] %in% groups
    return(samples[isExcluded & isGroup])
  }
  dataList$excludedSamples <- excludedSamples
  includedSamples <- function(groupSampleDataFrame, groups = dataList$groups){
    #dataList$groupSampleDataFrame[, "Sample"][!dataList$groupSampleDataFrame[, "Exclude"]]
    samples    =  groupSampleDataFrame[, "Sample"]
    isIncluded = !groupSampleDataFrame[, "Exclude"]
    isGroup    =  groupSampleDataFrame[, "Group"] %in% groups
    return(samples[isIncluded & isGroup])
  }
  dataList$includedSamples <- includedSamples
  
  includedGroups <- function(groupSampleDataFrame, samples = dataList$groupSampleDataFrame[, "Sample"]){
    unique(unlist(lapply(X = intersect(samples, dataList$includedSamples(groupSampleDataFrame)), FUN = function(sampleName){
      dataList$groupNameFunctionFromDataColumnName(dataColumnName = sampleName, sampleNamesToExclude = dataList$excludedSamples(dataList$groupSampleDataFrame))
    })))
  }
  dataList$includedGroups <- includedGroups
  excludedGroups <- function(groupSampleDataFrame, samples = dataList$groupSampleDataFrame[, "Sample"]){
    setdiff(dataList$groups, dataList$includedGroups(groupSampleDataFrame, samples)) 
  }
  dataList$excludedGroups <- excludedGroups
  
  
  ## 950 932 688
  ## 634 336 248
  ## 321 972 296
  ##   9 090 088
  ##  11 753 432
  ##  13 272 240
  #print(sort( sapply(ls(),function(x){object.size(get(x))})))
  #memory.profile()
  
  return(dataList)
}

## tagsSector <- dataFrameMS1Header[2, ]
## dataList$dataFrameInfos <- metaboliteProfile
processMS1data <- function(
  sampleNamesToExclude, numberOfMS1features, precursorLabels, 
  groups, metaboliteProfileColumnNames, 
  dataColumnIndecesFunctionFromGroupIndex, dataColumnsNameFunctionFromGroupIndex, dataColumnsNameFunctionFromGroupName, dataColumnsNameFunctionFromGroupNames, groupNameFunctionFromDataColumnName,
  tagsSector, metaboliteProfile, progress
){
  numberOfGroups <- length(groups)
  
  ####################
  ## MS1 measurement data: mean and LFC
  if(progress)  incProgress(amount = 0.1, detail = "Coloring") else print("Coloring")
  if(progress)  incProgress(amount = 0, detail = "Coloring init") else print("Coloring init")
  
  dataFrameMeasurements <- data.frame(matrix(nrow = numberOfMS1features, ncol = 0))
  rownames(dataFrameMeasurements) <- precursorLabels
  
  ## column name functions
  if(progress)  incProgress(amount = 0, detail = "Coloring naming functions") else print("Coloring naming functions")
  
  ## store data of groups
  dataColumnNames <- list()
  for(groupIdx in seq_len(numberOfGroups)){
    dataColumnNamesHere <- dataColumnsNameFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)
    dataColumnNames <- c(dataColumnNames, dataColumnNamesHere)
    dataFrameMeasurements[, dataColumnNamesHere] <- data.matrix(metaboliteProfile[, dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude), drop = FALSE])
  }
  dataColumnNames <- unlist(dataColumnNames)
  
  dataMeanColumnNameFunctionFromName  <- function(group){
    return(paste(group, "_mean", sep = ""))
  }
  dataMeanColumnNameFunctionFromIndex  <- function(groupIdx){
    return(dataMeanColumnNameFunctionFromName(groups[[groupIdx]]))
  }
  
  lfcColumnNameFunctionFromName <- function(groupOne, groupTwo){
    return(paste("LFC", groupOne, "vs", groupTwo, sep = "_"))
  }
  lfcColumnNameFunctionFromIndex <- function(groupIdxOne, groupIdxTwo){
    lfcColumnNameFunctionFromName(groups[[groupIdxOne]], groups[[groupIdxTwo]])
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
  for(groupIdx in seq_len(numberOfGroups)){
    dataMeanColumnName <- dataMeanColumnNameFunctionFromIndex(groupIdx)
    dataMeanColumnNames[[groupIdx]] <- dataMeanColumnName
    if(class(unlist(metaboliteProfile[, dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)])) == "character")
      for(colIdx in dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude))
        metaboliteProfile[, colIdx] <- as.numeric(metaboliteProfile[, colIdx])
    
    dataFrameMeasurements[, dataMeanColumnName] <- apply(X = data.matrix(metaboliteProfile[, dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude), drop=FALSE]), MARGIN = 1, FUN = mean)
    dataFrameMeasurements[is.na(dataFrameMeasurements[, dataMeanColumnName]), dataMeanColumnName] <- 0
  }
  dataMeanColumnNames <- unlist(dataMeanColumnNames)
  
  ## all replicates mean
  dataFrameMeasurements[, "meanAllNormed"] <- apply(
    X = data.matrix(metaboliteProfile[, 
                                      unlist(lapply(X = seq_len(numberOfGroups), FUN = function(x) {dataColumnIndecesFunctionFromGroupIndex(groupIdx = x, sampleNamesToExclude = sampleNamesToExclude)})),
                                      drop=FALSE]), 
    MARGIN = 1, FUN = mean
  )
  
  meanAllMax <- max(dataFrameMeasurements[, "meanAllNormed"])
  if(meanAllMax != 0)
    dataFrameMeasurements[, "meanAllNormed"] <- dataFrameMeasurements[, "meanAllNormed"] / meanAllMax
  
  ## log fold change between groups
  lfcColumnNames <- list()
  for(groupIdx1 in seq_len(numberOfGroups))
    for(groupIdx2 in seq_len(numberOfGroups)){
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
  
  
  #########################################################################################
  ## MS1 measurement data to colors
  if(progress)  incProgress(amount = 0, detail = "Coloring matrix") else print("Coloring matrix")
  
  matrixDataFrame <- data.matrix(dataFrameMeasurements)
  
  matrixDataFrame[, dataColumnNames    ][matrixDataFrame[, dataColumnNames    ] < 1] <- 1
  matrixDataFrame[, dataMeanColumnNames][matrixDataFrame[, dataMeanColumnNames] < 1] <- 1
  #matrixDataFrame[matrixDataFrame[, dataMeanColumnNames] < 1] <- 1
  
  matrixDataFrame[, dataColumnNames]     <- log10(matrixDataFrame[, dataColumnNames])
  matrixDataFrame[, dataMeanColumnNames] <- log10(matrixDataFrame[, dataMeanColumnNames])
  matrixDataFrame[is.infinite(matrixDataFrame)] <- 0
  #matrixDataFrame[matrixDataFrame < 0] <- 0
  
  ## min / max
  logAbsMin <- min(0, min(matrixDataFrame[, dataMeanColumnNames]))
  #logAbsMax <- max(matrixDataFrame[, dataMeanColumnNames])
  logAbsMax <- max(matrixDataFrame[, c(dataColumnNames, dataMeanColumnNames)])
  #logAbsMax <- max(matrixDataFrame[, dataColumnsNameFunctionFromGroupNames(groups = groups)])
  logFoldChangeMinMax <- c(min(matrixDataFrame[, lfcColumnNames]), max(matrixDataFrame[, lfcColumnNames]))
  logFoldChangeMax <- max(abs(logFoldChangeMinMax))
  if(logFoldChangeMax < 1)
    logFoldChangeMax <- 1
  
  ## maps
  colorMapAbsoluteData  <- makecmap(
    x = c(logAbsMin, logAbsMax), n = 100, 
    #colFn = colorRampPalette(c('white', 'black'))
    colFn = colorRampPalette(rainbow(18)[10:1])
  )
  colorMapLogFoldChange <- makecmap(
    x = c(-logFoldChangeMax, logFoldChangeMax), n = 100, 
    colFn = colorRampPalette(c('blue', 'white', 'red'))
  )
  
  columnGroupLabels <- sapply(X = groups, FUN = function(x){ rep(x = x, times = length(dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude))) })
  
  ## translate and box colors
  if(progress)  incProgress(amount = 0, detail = "Coloring box") else print("Coloring box")
  colorDataFrame <- dataFrameMeasurements
  colorDataFrame[, dataColumnNames    ] <- cmap(x = matrixDataFrame[, dataColumnNames    ], map = colorMapAbsoluteData)
  colorDataFrame[, dataMeanColumnNames] <- cmap(x = matrixDataFrame[, dataMeanColumnNames], map = colorMapAbsoluteData)
  colorDataFrame[, lfcColumnNames]      <- cmap(x = matrixDataFrame[, lfcColumnNames     ], map = colorMapLogFoldChange)
  colorMatrixDataFrame <- as.matrix(colorDataFrame)
  
  returnObj <- list(
    ## name functions
    #dataColumnsNameFunctionFromGroupIndex=dataColumnsNameFunctionFromGroupIndex,
    #dataColumnsNameFunctionFromGroupName=dataColumnsNameFunctionFromGroupName,
    #dataColumnsNameFunctionFromGroupNames=dataColumnsNameFunctionFromGroupNames,
    #groupNameFunctionFromDataColumnName=groupNameFunctionFromDataColumnName,
    dataMeanColumnNameFunctionFromIndex=dataMeanColumnNameFunctionFromIndex,
    dataMeanColumnNameFunctionFromName=dataMeanColumnNameFunctionFromName,
    lfcColumnNameFunctionFromIndex=lfcColumnNameFunctionFromIndex,
    lfcColumnNameFunctionFromName=lfcColumnNameFunctionFromName,
    groupNameFromGroupIndex=groupNameFromGroupIndex,
    groupIdxFromGroupName=groupIdxFromGroupName,
    ## data and names
    dataFrameMeasurements=dataFrameMeasurements,
    #dataMeanColumnNames=dataMeanColumnNames,
    #lfcColumnNames=lfcColumnNames,
    ## colors
    colorMatrixDataFrame=colorMatrixDataFrame,
    #matrixDataFrame=matrixDataFrame,
    colorMapAbsoluteData=colorMapAbsoluteData,
    colorMapLogFoldChange=colorMapLogFoldChange,
    columnGroupLabels=columnGroupLabels,
    ## constants
    meanAllMax=meanAllMax,
    logFoldChangeMax=logFoldChangeMax,
    logAbsMax=logAbsMax
  )
}

serializeSampleSelectionAndOrder <- function(groupSampleDataFrame){
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

#########################################################################################
## data fetching
getMetFragLink <- function(dataList, precursorIndex){
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
           print(paste("###### Unknown adduct (", adduct, ")!", sep = ""))
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
  
  #fragmentStrings <- paste(fragmentsPositiveX, fragmentsPositiveY, sep = "_", collapse = ";")
  
  ## https://msbi.ipb-halle.de/MetFragBeta/landing.xhtml?FragmentPeakMatchAbsoluteMassDeviation=0.01&FragmentPeakMatchRelativeMassDeviation=10&DatabaseSearchRelativeMassDeviation=10&PeakList=110_100;210_100&IonizedPrecursorMass=200.101&MetFragDatabaseType=PubChem
  #FragmentPeakMatchAbsoluteMassDeviation
  #FragmentPeakMatchRelativeMassDeviation
  #DatabaseSearchRelativeMassDeviation
  #PrecursorCompoundIDs
  #IonizedPrecursorMass
  #NeutralPrecursorMolecularFormula
  #PrecursorIonMode
  #IonizedPrecursorMass
  #landingPageUrl <- paste(sep = "",
  #                        "https://msbi.ipb-halle.de/MetFragBeta/landing.xhtml", "?",
  #                        "NeutralPrecursorMass", "=", neutralMass, "&",
  #                        "PrecursorIonMode", "=", ionMode, "&",
  #                        "MetFragDatabaseType", "=", "PubChem", "&",
  #                        "PeakList", "=", fragmentStrings
  #)
  
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
getMS2spectrumInfoForPrecursorLeaf <- function(dataList, clusterDataList, treeLabel){
  if(treeLabel >= 0)
    return(NULL)
  ###############################################
  ## leaf
  precursorIndex <- clusterDataList$filterObj$filter[[-treeLabel]]
  return(getMS2spectrumInfoForPrecursor(dataList, precursorIndex))
}
getMS2spectrumInfoForPrecursor <- function(dataList, precursorIndex){
  
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
  metFragLinkList <- getMetFragLink(dataList, precursorIndex)
  
  
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
getMS2spectrumInfoForCluster <- function(dataList, clusterDataList, treeLabel){
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
  #fragmentsY <- apply(X = data.matrix(dataList$featureMatrix[clusterMembersPrecursors, featuresIntersection]), MARGIN = 2, FUN = mean)
  fragmentsX <- dataList$fragmentMasses[featuresUnion]
  fragmentsY <- apply(X = data.matrix(dataList$featureMatrix[clusterMembersPrecursors, featuresUnion]), MARGIN = 2, FUN = mean)
  
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
getTableFromPrecursorSet <- function(dataList, precursorSet){
  ###############################################
  ## table data
  numberOfPrecursors <- length(precursorSet)
  
  ## measurements
  columnNames <- unlist(lapply(X = dataList$groups, FUN = dataList$dataMeanColumnNameFunctionFromName))
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
  featureMatrix <- format(x = featureMatrix, digits = 0, nsmall = 4)
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
getSpectrumStatistics <- function(dataList, precursorSet){
  if(FALSE){
    dataList_ <<- dataList
    precursorSet_ <<- precursorSet
  }
  if(FALSE){
    dataList <- dataList_
    precursorSet <- precursorSet_
  }
  
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
numericVectorToStringForEval <- function(vec){
  return(paste("c(", paste(vec, collapse = ","), ")", sep = ""))
}
colorVectorToStringForEval <- function(vec){
  return(paste("c('", paste(vec, collapse = "','"), "')", sep = ""))
}
