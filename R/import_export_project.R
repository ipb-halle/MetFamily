# Import and export input and project files
# - Parameter list
# - Project file



#' Create a project dataList
#' 
#' Wrapper for a clunky workflow
#'
#' @param ms1_path file for MS1 intensity .txt file
#' @param ms2_path file for MS2 .msp file
#' @param annot_path file for annotations (.csv or .tsv)
#' @param parameterSet list of parameters to use
#' @param siriusFileColumnName One of "NPC class", "NPC superclass", "NPC pathway", "ClassyFire subclass",
#'  "ClassyFire class", "ClassyFire superclass"
#' @returns dataList
#' @export
projectFromFiles <- function(
    ms1_path, ms2_path, 
    annot_path = NULL,
    siriusFileColumnName = c("ClassyFire superclass"),
    parameterSet = NULL
) {
  
  if (is.null(parameterSet)) {
    parameterSet <- parameterSetDefault()
  }
  
  resultObj <- convertToProjectFile(filePeakMatrixPath = ms1_path, 
                                    fileSpectra = ms2_path, parameterSet) 
  
  lines <- sparseMatrixToString(resultObj$matrixRows, resultObj$matrixCols,
                                resultObj$matrixVals, parameterSet)
  
  dataList <- readProjectData(lines)
  
  if(!is.null(annot_path)) {
    dataList <- add_qfeatures(dataList, qfeatures = resultObj$qfeatures, 
                              fileAnnotation = annot_path, siriusFileColumnName)
  }
  
  dataList
}


#' Default values
#'
#' @returns list parameterSet
#' @export
importParameterSetInit <- function() {
  
  list(
    minimumIntensityOfMaximalMS2peak = 2000,
    minimumNumbersOfFragments = 2,
    minimumAbsoluteMS2peaks = 500,
    minimumProportionOfMS2peaks = 0.05,
    neutralLossesPrecursorToFragments = TRUE,
    neutralLossesFragmentsToFragments = FALSE,
    mzDeviationAbsolute_grouping = 0.01,
    mzDeviationInPPM_grouping = 10,
    showImportParametersAdvanced = FALSE,
    doPrecursorDeisotoping = TRUE,
    mzDeviationAbsolute_precursorDeisotoping = 0.01,
    mzDeviationInPPM_precursorDeisotoping = 10,
    maximumRtDifference = 0.05,
    doMs2PeakGroupDeisotoping = TRUE,
    mzDeviationAbsolute_ms2PeakGroupDeisotoping = 0.01,
    mzDeviationInPPM_ms2PeakGroupDeisotoping = 10
  )
  
}


#' Create a parameterSet object
#'
#' To use in a script
#'
#' @returns a parameterSet list object of default values used by the app
#' @export
parameterSetDefault <- function() {
  
  ipsi <- importParameterSetInit()
  
  c(
    list(
      projectName = paste0("MetFamily project (created ", 
                           gsub(" ", "_", gsub(":", ".", Sys.time())),
                           ")"),
      projectDescription = "",
      toolVersion = "MetFamily 1.0"
    ),
    
    # assigned from app input in ui.R
    ipsi[c("minimumIntensityOfMaximalMS2peak", "minimumNumbersOfFragments",
           "minimumAbsoluteMS2peaks", "minimumProportionOfMS2peaks", 
           "mzDeviationAbsolute_grouping", "mzDeviationInPPM_grouping", 
           "doPrecursorDeisotoping", "mzDeviationAbsolute_precursorDeisotoping", 
           "mzDeviationInPPM_precursorDeisotoping", "maximumRtDifference", 
           "doMs2PeakGroupDeisotoping", "mzDeviationAbsolute_ms2PeakGroupDeisotoping", 
           "mzDeviationInPPM_ms2PeakGroupDeisotoping")],
    
    # fixed parameters, see server_guiTabInput
    list(
      proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9,
      mzDeviationAbsolute_mapping = 0.01
    ),
    
    # more app inputs
    ipsi[c("neutralLossesPrecursorToFragments",
           "neutralLossesFragmentsToFragments")]
  )
}


#' Flatten listObject to string for export
#' 
#' Also used before running readProjectData
#'
#' @param matrixRows 
#' @param matrixCols 
#' @param matrixVals 
#' @param parameterSet 
#'
#' @returns big string object
#' @export
sparseMatrixToString <- function(
    matrixRows, 
    matrixCols, 
    matrixVals, parameterSet
) {
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


#' Read MetFamily Project data saved by the export function
#'
#' Supports reading from plain and gzip'ed files
#' 
#' @param file Path to file to read
#' @param progress Whether to update a shiny Progress bar
#'
#' @return A big dataList. 
#' 
#' @seealso [readProjectData]
#' @export
readProjectFile <- function(file, progress = FALSE) {
  if(!is.na(progress)) {
    if(progress) {
      setProgress(value = 0, detail = "Parsing")
    } else {
      print("Parsing")
    }
  }
  
  extension <- tools::file_ext(file)
  
  if(extension == "gz") {
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


#' Read MetFamily Project data saved by the export function
#'
#' Read in a project described as lines. Main workhorse used when loading a project
#' file or separate input files
#'
#' @param fileLines Character vector with content of a project file
#' @param progress Whether to update a shiny Progress bar
#'
#' @return A big dataList. 
#' 
#' @seealso [processMS1data]
#' @export
#' @importFrom stringr str_split
readProjectData <- function(fileLines, progress = FALSE) {
  allowedTags <- c("ID")
  allowedTagPrefixes <- c("AnnotationColors=")
  
  ##################################################################################################
  ## parse data
  if(!is.na(progress)) {
    if(progress) {  
      incProgress(amount = 0.1, detail = "Preprocessing") 
    }
  } else { 
    print("Preprocessing")
  }
  
  numberOfRows <- length(fileLines)
  numberOfMS1features <- as.integer(numberOfRows - 3)
  
  ## header
  line1Tokens <- strsplit(x = fileLines[[1]], split = "\t")[[1]]
  line2Tokens <- strsplit(x = fileLines[[2]], split = "\t")[[1]]
  line3Tokens <- strsplit(x = fileLines[[3]], split = "\t")[[1]]
  
  # extract version
  versionString <- line1Tokens[3]
  line1Tokens[3] <- ""
  
  # For backwards compatibilty, change , to ; in annotations
  # old project files have no semicolon, with entries separates by commas
  annoField <- line2Tokens[3]
  if(!stringr::str_detect(annoField, ";") &&
     stringr::str_count(annoField, ",") == stringr::str_count(annoField, "=") - 2) {
    line2Tokens[3] <- stringr::str_replace_all(line2Tokens[3], ",", ";")
  }
  
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
    importParameters <- parameterSetDefault() %>% serializeParameterSet
    
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
  
  if(any(duplicated(metaboliteProfileColumnNames))) {
    stop(paste("Duplicated column names in the metabolite profile: ", 
               paste(sort(unique(metaboliteProfileColumnNames[duplicated(metaboliteProfileColumnNames)])), collapse = "; ")))
  }
  
  #########################################################################
  ## extract metabolite profile and fragment matrix
  metaboliteProfile <- as.data.frame(matrix(nrow = numberOfMS1features, 
                                            ncol = numberOfMetaboliteProfileColumns))
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
      if(!is.na(progress))  
        if(progress)  
          incProgress(amount = rowProgress*0.2,     
                      detail = paste("Preprocessing ", rowIdx, " / ", numberOfMS1features, sep = "")) 
      else 
        print(paste("Preprocessing ", rowIdx, " / ", numberOfMS1features, sep = ""))
    }
    
    lineIdx <- rowIdx + 3
    tokens <- stringr::str_split(string = fileLines[[lineIdx]], pattern = "\t")[[1]]
    
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
  
  # gp: Previous location of "Start of importing annotation part1 from two"
  
  listMatrixVals <- NULL
  
  ## header
  dataFrameHeader <- cbind(
    data.frame(rbind(
      c(importParameters, rep(x = "", times = numberOfMetaboliteProfileColumns - 1)),
      tagsSector,
      metaboliteProfileColumnNames), stringsAsFactors = FALSE),
    data.frame(rbind(
      fragmentGroupsNumberOfFramgents,
      fragmentGroupsAverageIntensity,
      fragmentGroupsAverageMass
    ), stringsAsFactors = FALSE)
  )
  headerLabels <- c("HeaderForFragmentCounts", 
                    "HeaderForGroupsAndFragmentIntensities", 
                    "Header")
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
    metaboliteProfileColumnNames <- c(metaboliteProfileColumnNames[seq_len(target)], 
                                      annotationColumnName, 
                                      metaboliteProfileColumnNames[(target+1):numberOfMetaboliteProfileColumns])
    colnames(metaboliteProfile) <- metaboliteProfileColumnNames
    headerColumnNames <- c(metaboliteProfileColumnNames, fragmentGroupsAverageMass)
    colnames(dataFrameHeader) <- headerColumnNames
    
    dataFrameHeader[2, target + 1] <- annotationColorsMapInitValue
    dataFrameHeader[3, target + 1] <- annotationColumnName
  }  # gp: Previous location of "Start of importing annotation part2 from two"
  
  annotationColumnIndex <- which(metaboliteProfileColumnNames == annotationColumnName)
  annotationColorsValue <- dataFrameHeader[2, annotationColumnIndex]
  
  dataFrameMS1Header <- dataFrameHeader[, seq_len(numberOfMetaboliteProfileColumns)]
  
  ##################################################################################################
  ## MS1 feature IDs
  
  ## mz/rt is aligned by '.'
  mzs <- metaboliteProfile[, "m/z"]
  rts <- metaboliteProfile[, "RT"]
  
  ## add .0 if necessary
  for(i in seq_len(numberOfMS1features)) {
    if(length(grep(x = mzs[[i]], pattern = ".*\\..*")) == 0) {
      mzs[[i]] <- paste(mzs[[i]], ".0", sep = "") }}
  
  for(i in seq_len(numberOfMS1features)) {
    if(length(grep(x = rts[[i]], pattern = ".*\\..*")) == 0) {
      rts[[i]] <- paste(rts[[i]], ".0", sep = "") }}
  
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
  
  # TODO document this block
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
    
    if(nchar(mzBefore) < mzMaxBefore) {
      mzBefore <- paste(
        paste(rep(x = "  ", times = mzMaxBefore - nchar(mzBefore)), collapse = ""),
        mzBefore,
        sep = "")
    }
    if(nchar(mzAfter) > maximumNumberOfDecimalPlacesForMz) {
      mzAfter <- substr(x = mzAfter, start = 1, stop = maximumNumberOfDecimalPlacesForMz)
    }
    if(nchar(mzAfter) < maximumNumberOfDecimalPlacesForMz) {
      mzAfter <- paste(
        mzAfter,
        paste(rep(x = "0", times = maximumNumberOfDecimalPlacesForMz - nchar(mzAfter)), collapse = ""),
        sep = "")
    }
    if(nchar(rtBefore) < rtMaxBefore) {
      rtBefore <- paste(
        paste(rep(x = "  ", times = rtMaxBefore - nchar(rtBefore)), collapse = ""),
        rtBefore,
        sep = "")
    }
    if(nchar(rtAfter) > maximumNumberOfDecimalPlacesForRt) {
      rtAfter <- substr(x = rtAfter, start = 1, stop = maximumNumberOfDecimalPlacesForRt)
    }
    if(nchar(rtAfter) < maximumNumberOfDecimalPlacesForRt) {
      rtAfter <- paste(
        rtAfter,
        paste(rep(x = "0", times = maximumNumberOfDecimalPlacesForRt - nchar(rtAfter)), collapse = ""),
        sep = "")
    }
    
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
      if(length(indeces1) == 0) {
        next
      }
      
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
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Features") else print("Features")
  
  ## get features
  featureIndeces <- list()
  featureCount <- vector(mode = "numeric", length = numberOfMS1features)
  
  for(i in seq_len(numberOfMS1features)){
    indecesHere <- which(matrixRows == i)
    featureIndecesHere <- matrixCols[indecesHere]
    numberOfFeatures <- length(featureIndecesHere)
    
    featureIndeces[[i]] <- featureIndecesHere
    featureCount[[i]] <- numberOfFeatures
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Feature postprocessing") else print("Feature postprocessing")
  
  ## ms2 plot data
  ms2PlotDataNumberOfFragments <- fragmentGroupsNumberOfFramgents
  ms2PlotDataAverageAbundance  <- fragmentGroupsAverageIntensity
  ms2PlotDataFragmentMasses    <- fragmentGroupsAverageMass
  maxNumberOfFragments <- max(ms2PlotDataNumberOfFragments)
  ms2PlotDataColorMapFragmentData  <- squash::makecmap(
    x = c(0, maxNumberOfFragments), n = 100, 
    colFn = colorRampPalette(c('grey', 'black'))
  )
  
  ## featureMatrix and annotation
  featureMatrix <- sparseMatrix(i = matrixRows, j = matrixCols, x = matrixVals)
  matrixRows <- NULL
  matrixCols <- NULL
  matrixVals <- NULL
  
  rownames(featureMatrix) <- precursorLabels
  colnames(featureMatrix) <- fragmentGroupsAverageMass
  
  ## featureIndexMatrix
  featureIndexMatrix <- matrix(nrow = numberOfMS1features, ncol = max(sapply(X = featureIndeces, FUN = length)))
  rownames(featureIndexMatrix) <- precursorLabels
  for(i in seq_len(numberOfMS1features)) {
    featureIndexMatrix[i, seq_len(length(featureIndeces[[i]]))] <- featureIndeces[[i]]
  }
  
  minimumMass <- min(fragmentGroupsAverageMass)
  maximumMass <- max(fragmentGroupsAverageMass)
  
  ##################################################################################################
  ## process sample measurements
  
  ## sample columns
  sampleColumns <- tagsSector != ""
  for(allowedTag in allowedTags) {
    sampleColumns[grep(x = tagsSector, pattern = paste("^", allowedTag, "$", sep = ""))] <- FALSE
  }
  for(allowedTagPrefix in allowedTagPrefixes) {
    sampleColumns[grep(x = tagsSector, pattern = paste("^", allowedTagPrefix, sep = ""))] <- FALSE
  }
  sampleColumns <- which(sampleColumns)
  sampleColumnsStartEnd <- c(min(sampleColumns), max(sampleColumns))
  
  sampleClasses <- unique(tagsSector[sampleColumns])
  numberOfGroups <- length(sampleClasses)
  
  sampleNamesToExclude <- NULL
  
  
  dataColumnIndecesFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude = NULL){
    which(tagsSector == sampleClasses[[groupIdx]] & !(metaboliteProfileColumnNames %in% sampleNamesToExclude))
  }
  dataColumnsNameFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude = NULL){
    sampleNames = metaboliteProfileColumnNames[dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)]
    return(sampleNames)
  }
  dataColumnsNameFunctionFromGroupName <- function(group, sampleNamesToExclude = NULL){
    dataColumnsNameFunctionFromGroupIndex(groupIdx = match(x = group, table = sampleClasses), sampleNamesToExclude = sampleNamesToExclude)
  }
  dataColumnsNameFunctionFromGroupNames <- function(sampleClasses, sampleNamesToExclude = NULL){
    unlist(lapply(X = sampleClasses, FUN = function(x){dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)}))
  }
  groupNameFunctionFromDataColumnName <- function(dataColumnName, sampleNamesToExclude = NULL){
    groupIdx <- which(unlist(lapply(X = sampleClasses, FUN = function(x){
      dataColumnNames <- dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)
      any(dataColumnNames == dataColumnName)
    })))
    sampleClasses[[groupIdx]]
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
    sampleClasses=sampleClasses, metaboliteProfileColumnNames=metaboliteProfileColumnNames, tagsSector = tagsSector, 
    dataColumnIndecesFunctionFromGroupIndex=dataColumnIndecesFunctionFromGroupIndex, 
    dataColumnsNameFunctionFromGroupIndex=dataColumnsNameFunctionFromGroupIndex, 
    dataColumnsNameFunctionFromGroupName=dataColumnsNameFunctionFromGroupName, 
    dataColumnsNameFunctionFromGroupNames=dataColumnsNameFunctionFromGroupNames, 
    groupNameFunctionFromDataColumnName=groupNameFunctionFromDataColumnName,
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
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Feature annotations") else print("Feature annotations")
  annotationValueIgnore <- "Ignore"
  annotationColorIgnore <- "red"
  
  
  ## present annotations
  annotations    <- vector(mode='list', length=numberOfMS1features)
  annoVals <- metaboliteProfile[, annotationColumnName]
  
  for(i in seq_len(numberOfMS1features)){
    if(nchar(annoVals[[i]]) > 0){
      annotations[[i]] <- as.list(unlist(strsplit(x = annoVals[[i]], split = "; ")))
    }
    else{
      annotations[[i]] <- list()
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
    annoArrayIsArtifact[[i]] <- ignoreThere
  }
  
  ## present annos with colors
  annotationColorsMapValue <- substr(
    x = annotationColorsValue, 
    start = nchar(paste(annotationColorsName, "={", sep = "")) + 1, 
    stop = nchar(annotationColorsValue) - nchar("}")
  )
  
  if(nchar(annotationColorsMapValue) > 0){
    
    annotationColorsMapValuePairs <- unlist(strsplit(x = annotationColorsMapValue, split = "; "))
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
  if(length(annotationColorsMapKeys) > 0) {
    for(i in seq_along(annotationColorsMapKeys)){
      annoPresentAnnotationsList[[1 + i]] <- annotationColorsMapKeys[[i]]
      annoPresentColorsList     [[1 + i]] <- annotationColorsMapValues[[i]]
    }
  }
  
  ## check consistency
  if(!all(unique(unlist(annoArrayOfLists)) %in% unlist(annoPresentAnnotationsList))){
    missing <- unique(unlist(annoArrayOfLists))[!(unique(unlist(annoArrayOfLists)) %in% unlist(annoPresentAnnotationsList))]
    stop(paste("Annotation(s)", paste(missing, collapse = "; "), "missing in present annotations list"))
  }
  # FIX annotations "" missing
  if(!all(unlist(annoPresentAnnotationsList) %in% unique(c(annotationValueIgnore, unlist(annoArrayOfLists))))){
    missing <- unlist(annoPresentAnnotationsList)[!(unlist(annoPresentAnnotationsList) %in% unique(c(annotationValueIgnore, unlist(annoArrayOfLists))))]
    stop(paste("Present annotation(s)", paste(missing, collapse = "; "), "missing in annotations"))
  }
  
  ##################################################################################################
  ## box
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Boxing") else print("Boxing")
  dataList <- list()
  ## data
  dataList$dataFrameHeader <- dataFrameHeader
  dataList$dataFrameMS1Header <- dataFrameMS1Header
  dataList$dataFrameInfos <- metaboliteProfile
  dataList$importParameterSet <- importParameterSet
  dataList$numberOfPrecursors <- numberOfMS1features
  dataList$numberOfDuplicatedPrecursors <- numberOfDuplicated
  dataList$sampleClasses <- sampleClasses
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
  
  if(!is.na(progress))  if(progress)  setProgress(1) else print("Ready")
  
  
  ## redefine MS1 column functions
  dataColumnIndecesFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude){
    which(dataList$tagsSector == dataList$sampleClasses[[groupIdx]] & !(dataList$metaboliteProfileColumnNames %in% sampleNamesToExclude))
  }
  dataList$dataColumnIndecesFunctionFromGroupIndex <- dataColumnIndecesFunctionFromGroupIndex
  
  dataColumnsNameFunctionFromGroupIndex <- function(groupIdx, sampleNamesToExclude){
    dataList$metaboliteProfileColumnNames[dataList$dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)]
  }
  dataList$dataColumnsNameFunctionFromGroupIndex <- dataColumnsNameFunctionFromGroupIndex
  
  dataColumnsNameFunctionFromGroupName <- function(group, sampleNamesToExclude){
    dataColumns <- dataList$dataColumnsNameFunctionFromGroupIndex(groupIdx = match(x = group, table = dataList$sampleClasses), sampleNamesToExclude = sampleNamesToExclude)
  }
  dataList$dataColumnsNameFunctionFromGroupName <- dataColumnsNameFunctionFromGroupName
  
  dataColumnsNameFunctionFromGroupNames <- function(sampleClasses, sampleNamesToExclude){
    unlist(lapply(X = sampleClasses, FUN = function(x){
      dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)
    }))
  }
  dataList$dataColumnsNameFunctionFromGroupNames <- dataColumnsNameFunctionFromGroupNames
  
  groupNameFunctionFromDataColumnName <- function(dataColumnName, sampleNamesToExclude){
    groupIdx <- which(unlist(lapply(X = dataList$sampleClasses, FUN = function(x){
      dataColumnNames <- dataList$dataColumnsNameFunctionFromGroupName(group = x, sampleNamesToExclude = sampleNamesToExclude)
      any(dataColumnNames == dataColumnName)
    })))
    dataList$sampleClasses[[groupIdx]]
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
  excludedSamples <- function(groupSampleDataFrame, sampleClasses = dataList$sampleClasses){
    samples    =  groupSampleDataFrame[, "Sample"]
    isExcluded =  groupSampleDataFrame[, "Exclude"]
    isGroup    =  groupSampleDataFrame[, "Group"] %in% sampleClasses
    return(samples[isExcluded & isGroup])
  }
  dataList$excludedSamples <- excludedSamples
  
  includedSamples <- function(groupSampleDataFrame, sampleClasses = dataList$sampleClasses){
    samples    =  groupSampleDataFrame[, "Sample"]
    isIncluded = !groupSampleDataFrame[, "Exclude"]
    isGroup    =  groupSampleDataFrame[, "Group"] %in% sampleClasses
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
    setdiff(dataList$sampleClasses, dataList$includedGroups(groupSampleDataFrame, samples)) 
  }
  dataList$excludedGroups <- excludedGroups
  
  return(dataList)
}


#' Add qfeatures to dataList
#' 
#' 
#'
#' @param dataList Output from readProjectData.
#' @param qfeatures qfeature object, can be taken from resultObj$qfeatures.
#' @param fileAnnotation character Path for sirius annotation file.
#' @param siriusFileColumnName One of "NPC class", "NPC superclass", "NPC pathway", "ClassyFire subclass",
#'  "ClassyFire class", "ClassyFire superclass"
#' @returns The dataList object with added sirius annotations.
#' @export
add_qfeatures <- function(
    dataList,
    qfeatures,
    fileAnnotation = NULL,
    siriusFileColumnName = "ClassyFire superclass"
) {
  # This function takes snippets previously in convertToProjectFile and readProjectData
  # to streamline the process and declutter the aforementionned functions.
  
  # The function does not do anything without annotation file
  if (is.null(fileAnnotation)) {
    return(dataList)
  }
  
  qfeatures <- addSiriusAnnotations(qfeatures, siriusFile = fileAnnotation, 
                                    siriusFileColumnName = siriusFileColumnName)
  
  # previously used test, not sure if still needed
  if (is.null(attr(rowData(qfeatures[[1]]), "annotation column"))) {
    stop("No annotation")
  }
  
  #Start of importing  annotation part1 from two
  
  # Extract the relevant data: Alignment ID and the annotation column from qfeatures
  annot_colname <- attr(rowData(qfeatures[[1]]), "annotation column")
  annotation_data <- rowData(qfeatures[[1]])[[annot_colname]]
  alignment_ids <- rowData(qfeatures[[1]])[["Alignment ID"]]
  
  # Find the matching indices between metaboliteProfile and annotation_data
  metaboliteProfile <- dataList$dataFrameInfos
  matching_indices <- match(metaboliteProfile[["Alignment ID"]], alignment_ids)
  
  metaboliteProfile$Annotation[!is.na(matching_indices)] <- annotation_data[matching_indices[!is.na(matching_indices)]] 
  #eliminate NAs replace by "" so nchar(annoVals[[i]]) > 0 works in l. 597
  metaboliteProfile$Annotation[is.na(metaboliteProfile$Annotation)] <- "" 
  
  
  #Start of importing annotation part2 from two
  
  #adding HEX color codes from external annotations to the annotationColorsMapInitValue of dataFrameHeader
  
  # Copy the selected column by user, Remove duplicates and exclude the first row
  
  # Split annotations based on ;
  uniqueAnnotations0 <- unique(unlist(strsplit(metaboliteProfile$Annotation, ";")))
  uniqueAnnotations0 <- unique(metaboliteProfile$Annotation)
  
  # remove empty value
  uniqueAnnotations0 <- uniqueAnnotations0[!uniqueAnnotations0 == ""]
  
  uniqueAnnotations <- paste0(uniqueAnnotations0, "=")
  # Add a random string from the hex color list to each element of uniqueAnnotions
  # strings_list <- c("#000000", "#FFFFFF", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#808080", "#C0C0C0", "#FFA500", "#FFC0CB", "#FFD700", "#A52A2A")
  # uniqueAnnotations <- paste0(uniqueAnnotations, sample(strings_list, length(uniqueAnnotations), replace = TRUE))
  
  # gp: removed red as it is already used for ignore category
  # TODO automate color generation, do not limit it to 40, could always use same colors, no need for sampling
  allowedCols <- c("blue", "yellow", "green", "brown", "deepskyblue", "orange", "deeppink", "aquamarine", "burlywood", "cadetblue", "coral", "cornflowerblue", "cyan", "firebrick", "goldenrod", "indianred", "khaki", "magenta", "maroon", "beige", "moccasin", "olivedrab", "orangered", "orchid", "paleturquoise3", "rosybrown", "salmon", "seagreen3", "skyblue", "steelblue", "#BF360C", "#33691E", "#311B92", "#880E4F", "#1A237E", "#006064", "#004D40", "#FF6F00", "#E65100")
  uniqueAnnotations <- paste0(uniqueAnnotations, sample(allowedCols, length(uniqueAnnotations), replace = TRUE))
  
  # Format uniqueAnnotations into a single line with comma-separated values
  uniqueAnnotations1 <- paste(uniqueAnnotations, collapse = "; ")
  
  #uniqueAnnotationsHexs <- paste("AnnotationColors={", paste(uniqueAnnotations1, collapse = ","), "}")# this line introduces a space after the first Item of the object, therefore, replaced with the following to remove the space
  uniqueAnnotationsHexs <- gsub("AnnotationColors=\\{\\s+", "AnnotationColors={", paste0("AnnotationColors={", uniqueAnnotations1, "}"))
  # Put back in dataFrameHeader
  dataList$dataFrameHeader$Annotation[2] <- uniqueAnnotationsHexs
  
  
  # adjust the rest -----------------------------
  # gp: Some of these are using brute force instead of repeating logic from readProjectData
  # I think this is fine for now, needs to be tested more thoroughly later (or rewrite readProjectData)
  
  dataList$dataFrameInfos <- metaboliteProfile
  
  # dataFrameMS1Header - brute force
  dataList$dataFrameMS1Header$Annotation[2] <- uniqueAnnotationsHexs
  
  # annoArrayOfLists - brute force
  annoArrayOfLists1 <- 
    base::split(metaboliteProfile$Annotation, seq_along(metaboliteProfile$Annotation))[!metaboliteProfile$Annotation == ""] %>% 
    unname
  annoArrayOfLists2 <- base::split(annoArrayOfLists1, seq_along(annoArrayOfLists1))
  dataList$annoArrayOfLists[!metaboliteProfile$Annotation == ""] <- annoArrayOfLists2
  
  # annoPresentAnnotationsList
  dataList$annoPresentAnnotationsList <- c(list("Ignore"), uniqueAnnotations0 %>% as.list)
  
  # annoPresentColorsList
  dataList$annoPresentColorsList <- c(list("red"), uniqueAnnotations %>% strsplit("=") %>% lapply(function(x) x[2]))
  
  # tagsSector
  dataList$tagsSector[3] <- uniqueAnnotationsHexs
  
  dataList
  
}


#' Process MS-Dial-like MS1 data.frame
#' 
#' Processing of MS-Dial-like MS1 data.frame. Includes calculation 
#' of MS1 data mean and log-fold-change (LFC) data
#'
#' @param sampleNamesToExclude 
#' @param numberOfMS1features 
#' @param precursorLabels 
#' @param sampleClasses 
#' @param metaboliteProfileColumnNames 
#' @param dataColumnIndecesFunctionFromGroupIndex 
#' @param dataColumnsNameFunctionFromGroupIndex 
#' @param dataColumnsNameFunctionFromGroupName 
#' @param dataColumnsNameFunctionFromGroupNames 
#' @param groupNameFunctionFromDataColumnName 
#' @param tagsSector 
#' @param metaboliteProfile 
#' @param progress 
#'
#' @return ?
#' @export
#' @importFrom grDevices colorRampPalette rainbow
processMS1data <- function(
    sampleNamesToExclude, 
    numberOfMS1features, 
    precursorLabels, 
    sampleClasses, 
    metaboliteProfileColumnNames, 
    dataColumnIndecesFunctionFromGroupIndex, 
    dataColumnsNameFunctionFromGroupIndex, 
    dataColumnsNameFunctionFromGroupName, 
    dataColumnsNameFunctionFromGroupNames, 
    groupNameFunctionFromDataColumnName,
    tagsSector, 
    metaboliteProfile, 
    progress=FALSE
) {
  numberOfGroups <- length(sampleClasses)
  
  ####################
  ## MS1 measurement data: mean and LFC
  if(!is.na(progress))  
    if(progress)  
      incProgress(amount = 0.1, detail = "Coloring") 
  else 
    print("Coloring")
  
  if(!is.na(progress))  
    if(progress)  
      incProgress(amount = 0, detail = "Coloring init") 
  else 
    print("Coloring init")
  
  dataFrameMeasurements <- data.frame(matrix(nrow = numberOfMS1features, ncol = 0))
  rownames(dataFrameMeasurements) <- precursorLabels
  
  ## column name functions
  if(!is.na(progress))  
    if(progress)  
      incProgress(amount = 0, detail = "Coloring naming functions") 
  else 
    print("Coloring naming functions")
  
  ## store data of sampleClasses
  dataColumnNames <- list()
  for(groupIdx in seq_len(numberOfGroups)){
    dataColumnNamesHere <- dataColumnsNameFunctionFromGroupIndex(groupIdx = groupIdx, 
                                                                 sampleNamesToExclude = sampleNamesToExclude)
    dataColumnNames <- c(dataColumnNames, dataColumnNamesHere)
    dataFrameMeasurements[, dataColumnNamesHere] <- data.numericmatrix(metaboliteProfile[, dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, 
                                                                                                                                   sampleNamesToExclude = sampleNamesToExclude), 
                                                                                         drop = FALSE])
  }
  dataColumnNames <- unlist(dataColumnNames)
  
  dataMeanColumnNameFunctionFromName  <- function(group){
    return(paste(group, "_mean", sep = ""))
  }
  
  dataMeanColumnNameFunctionFromIndex  <- function(groupIdx){
    return(dataMeanColumnNameFunctionFromName(sampleClasses[[groupIdx]]))
  }
  
  lfcColumnNameFunctionFromName <- function(groupOne, groupTwo){
    return(paste("LFC", groupOne, "vs", groupTwo, sep = "_"))
  }
  
  lfcColumnNameFunctionFromIndex <- function(groupIdxOne, groupIdxTwo){
    lfcColumnNameFunctionFromName(sampleClasses[[groupIdxOne]], sampleClasses[[groupIdxTwo]])
  }
  
  groupNameFromGroupIndex <- function(groupIdx){
    return(sampleClasses[[groupIdx]])
  }
  
  groupIdxFromGroupName <- function(group){
    return(match(x = group, table = sampleClasses))
  }
  
  if(!is.na(progress))  
    if(progress)  incProgress(amount = 0, detail = "Coloring gather data") 
  else 
    print("Coloring gather data")
  
  ## mean data columns
  dataMeanColumnNames <- list()
  for(groupIdx in seq_len(numberOfGroups)){
    dataMeanColumnName <- dataMeanColumnNameFunctionFromIndex(groupIdx)
    dataMeanColumnNames[[groupIdx]] <- dataMeanColumnName
    if(is(unlist(metaboliteProfile[, dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude)]),"character"))
      for(colIdx in dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude))
        metaboliteProfile[, colIdx] <- as.numeric(metaboliteProfile[, colIdx])
    
    dataFrameMeasurements[, dataMeanColumnName] <- apply(X = data.numericmatrix(metaboliteProfile[, dataColumnIndecesFunctionFromGroupIndex(groupIdx = groupIdx, sampleNamesToExclude = sampleNamesToExclude), drop=FALSE]), MARGIN = 1, FUN = mean)
    dataFrameMeasurements[is.na(dataFrameMeasurements[, dataMeanColumnName]), dataMeanColumnName] <- 0
  }
  dataMeanColumnNames <- unlist(dataMeanColumnNames)
  
  ## all replicates mean
  dataFrameMeasurements[, "meanAllNormed"] <- apply(
    X = data.numericmatrix(metaboliteProfile[, 
                                             unlist(lapply(X = seq_len(numberOfGroups), FUN = function(x) {dataColumnIndecesFunctionFromGroupIndex(groupIdx = x, sampleNamesToExclude = sampleNamesToExclude)})),
                                             drop=FALSE]), 
    MARGIN = 1, FUN = mean
  )
  
  meanAllMax <- max(dataFrameMeasurements[, "meanAllNormed"])
  if(meanAllMax != 0)
    dataFrameMeasurements[, "meanAllNormed"] <- dataFrameMeasurements[, "meanAllNormed"] / meanAllMax
  
  ## log fold change between sampleClasses
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
  if(!is.na(progress))  
    if(progress)  
      incProgress(amount = 0, detail = "Coloring matrix") 
  else 
    print("Coloring matrix")
  
  matrixDataFrame <- data.numericmatrix(dataFrameMeasurements)
  
  matrixDataFrame[, dataColumnNames    ][matrixDataFrame[, dataColumnNames    ] < 1] <- 1
  matrixDataFrame[, dataMeanColumnNames][matrixDataFrame[, dataMeanColumnNames] < 1] <- 1
  
  matrixDataFrame[, dataColumnNames]     <- log10(matrixDataFrame[, dataColumnNames])
  matrixDataFrame[, dataMeanColumnNames] <- log10(matrixDataFrame[, dataMeanColumnNames])
  matrixDataFrame[is.infinite(matrixDataFrame)] <- 0
  
  ## min / max
  logAbsMin <- min(0, min(matrixDataFrame[, dataMeanColumnNames]))
  logAbsMax <- max(matrixDataFrame[, c(dataColumnNames, dataMeanColumnNames)])
  logFoldChangeMinMax <- c(min(matrixDataFrame[, lfcColumnNames]), max(matrixDataFrame[, lfcColumnNames]))
  logFoldChangeMax <- max(abs(logFoldChangeMinMax))
  if(logFoldChangeMax < 1)
    logFoldChangeMax <- 1
  
  ## maps
  colorMapAbsoluteData  <- squash::makecmap(
    x = c(logAbsMin, logAbsMax), n = 100, 
    colFn = colorRampPalette(rainbow(18)[10:1])
  )
  colorMapLogFoldChange <- squash::makecmap(
    x = c(-logFoldChangeMax, logFoldChangeMax), n = 100, 
    colFn = colorRampPalette(c('blue', 'white', 'red'))
  )
  
  columnGroupLabels <- sapply(X = sampleClasses, 
                              FUN = function(x){ 
                                rep(x = x, 
                                    times = length(dataColumnsNameFunctionFromGroupName(group = x, 
                                                                                        sampleNamesToExclude = sampleNamesToExclude))) 
                              })
  
  ## translate and box colors
  if(!is.na(progress)) {
    if(progress) incProgress(amount = 0, detail = "Coloring box")
  } else {
    print("Coloring box")
  }
  
  # TODO document why all the colors
  colorDataFrame <- dataFrameMeasurements
  colorDataFrame[, dataColumnNames    ] <- squash::cmap(x = matrixDataFrame[, dataColumnNames    ], map = colorMapAbsoluteData)
  colorDataFrame[, dataMeanColumnNames] <- squash::cmap(x = matrixDataFrame[, dataMeanColumnNames], map = colorMapAbsoluteData)
  colorDataFrame[, lfcColumnNames]      <- squash::cmap(x = matrixDataFrame[, lfcColumnNames     ], map = colorMapLogFoldChange)
  colorMatrixDataFrame <- as.matrix(colorDataFrame)
  
  returnObj <- list(
    ## name functions
    dataMeanColumnNameFunctionFromIndex=dataMeanColumnNameFunctionFromIndex,
    dataMeanColumnNameFunctionFromName=dataMeanColumnNameFunctionFromName,
    lfcColumnNameFunctionFromIndex=lfcColumnNameFunctionFromIndex,
    lfcColumnNameFunctionFromName=lfcColumnNameFunctionFromName,
    groupNameFromGroupIndex=groupNameFromGroupIndex,
    groupIdxFromGroupName=groupIdxFromGroupName,
    ## data and names
    dataFrameMeasurements=dataFrameMeasurements,
    ## colors
    colorMatrixDataFrame=colorMatrixDataFrame,
    colorMapAbsoluteData=colorMapAbsoluteData,
    colorMapLogFoldChange=colorMapLogFoldChange,
    columnGroupLabels=columnGroupLabels,
    ## constants
    meanAllMax=meanAllMax,
    logFoldChangeMax=logFoldChangeMax,
    logAbsMax=logAbsMax
  )
}


#' Export MetFamily project
#'
#' Creates a `.gz` file to store the MetFamily project.
#'
#' Implements the logic of prepareAllPrecursors in a generalized manner.
#' Replaces `createExportMatrix()` as a package function.
#'
#' @param dataList main dataset
#' @param file filename
#' @param precursorSet entries to include, defaults to all
#' @param compressed if FALSE, skip compressing. If True, compress
#'
#' @returns filename 
#' @export
writeProjectFile <- function(dataList, precursorSet = NULL, file, compressed=TRUE) {
  
  if (compressed) {
    if(tools::file_ext(file) != "gz") {
      file <- paste0(file, ".gz")
    }
  }
  
  if (is.null(precursorSet)) {
    precursorSet <- seq_len(dataList$numberOfPrecursors)
  }
  
  ## fragment matrix
  fragmentMatrix      <- dataList$featureMatrix[precursorSet, ]
  dgTMatrix <- as(fragmentMatrix, "dgTMatrix")
  matrixRows <- dgTMatrix@i + 1
  matrixCols <- dgTMatrix@j + 1
  matrixVals <- dgTMatrix@x
  
  numberOfColumns <- ncol(fragmentMatrix)
  numberOfRows <- nrow(fragmentMatrix)
  chunkSize <- 1000
  numberOfChunks <- ceiling(numberOfColumns / chunkSize)
  
  fragmentCounts      <- vector(mode = "integer", length = numberOfColumns)
  fragmentIntensities <- vector(mode = "numeric", length = numberOfColumns)
  fragmentMasses      <- dataList$fragmentMasses
  linesMatrix <- matrix(nrow = numberOfRows, ncol = numberOfChunks)
  
  for(chunkIdx in seq_len(numberOfChunks)) {
    colStart <- 1 + (chunkIdx - 1) * chunkSize
    colEnd <- colStart + chunkSize - 1
    if(chunkIdx == numberOfChunks) {
      colEnd <- numberOfColumns
    }
    
    numberOfColumnsHere <- colEnd - colStart + 1
    numberOfRowsHere <- max(matrixRows)
    indeces <- matrixCols >= colStart & matrixCols <= colEnd
    
    fragmentMatrixPart <- matrix(
      data = rep(x = "", times = numberOfRowsHere * numberOfColumnsHere), 
      nrow = numberOfRowsHere,
      ncol = numberOfColumnsHere)
    fragmentMatrixPart[cbind(matrixRows[indeces], matrixCols[indeces] - colStart + 1)] <- matrixVals[indeces]
    
    fragmentCountsPart      <- apply(X = fragmentMatrixPart, MARGIN = 2, FUN = function(x){ sum(x != "") })
    fragmentIntensitiesPart <- apply(X = fragmentMatrixPart, MARGIN = 2, FUN = function(x){ sum(as.numeric(x), na.rm = TRUE) }) / fragmentCountsPart
    
    linesPart <- apply(X = fragmentMatrixPart, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
    
    fragmentCounts[colStart:colEnd] <- fragmentCountsPart
    fragmentIntensities[colStart:colEnd] <- fragmentIntensitiesPart
    linesMatrix[, chunkIdx] <- linesPart
  }
  
  ## assemble
  linesFragmentMatrixWithHeader <- c(
    paste(fragmentCounts, collapse = "\t"),
    paste(fragmentIntensities, collapse = "\t"),
    paste(fragmentMasses, collapse = "\t"),
    apply(X = linesMatrix, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
  )
  
  ################################################################################
  ## MS1 matrix
  dataList$dataFrameMS1Header[[1,2]] <- serializeSampleSelectionAndOrder(dataList$groupSampleDataFrame)
  ms1Matrix     <- rbind(
    dataList$dataFrameMS1Header,
    dataList$dataFrameInfos[precursorSet, ]
  )
  ms1Matrix     <- as.matrix(ms1Matrix)
  
  ###########################################################
  ## export annotations
  
  ## process annotations
  annotations <- dataList$annoArrayOfLists
  for(i in 1:length(annotations)) {
    if(dataList$annoArrayIsArtifact[[i]]) {
      annotations[[i]] <- c(annotations[[i]], dataList$annotationValueIgnore)
    }
  }
  
  annotationStrings <- vector(mode = "character", length = length(annotations))
  # NOTE gp: some annotations contain commas, so use semi-colon to separate
  for(i in 1:length(annotations)) {
    if(length(annotations[[i]]) > 0) {
      annotationStrings[[i]] <- paste(annotations[[i]], collapse = "; ")
    } else {
      annotationStrings[[i]] <- ""
    }
  }
  annotationStrings <- annotationStrings[precursorSet]
  
  ## process annotation-color-map
  annoPresentAnnotations <- dataList$annoPresentAnnotationsList[-1]
  annoPresentColors      <- dataList$annoPresentColorsList[-1]
  
  if(length(annoPresentAnnotations) > 0) {
    annotationColors <- paste(annoPresentAnnotations, annoPresentColors, sep = "=", collapse = "; ")
  } else {
    annotationColors <- ""
  }
  annotationColors <- paste0(dataList$annotationColorsName, "={", annotationColors, "}")
  
  ## box
  annotationColumn <- c("", annotationColors, dataList$annotationColumnName, annotationStrings)
  
  ms1Matrix[, dataList$annotationColumnIndex] <- annotationColumn
  
  # add MetFamily version
  source(system.file("MetFamily/version.R", package = "MetFamily"), local = T)
  versionString <- paste0(
        R.Version()$version.string, 
        "; MetFamily build: ", metFamilyAppVersion, "_", system(command = "hostname", intern = TRUE),
        "; MetFamily package: ", packageVersion
      )
  
  ms1Matrix[1,3] <- versionString
  
  linesMS1MatrixWithHeader <- apply(X = ms1Matrix, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
  lines <- paste(linesMS1MatrixWithHeader, linesFragmentMatrixWithHeader, sep = "\t")
  
  ## save to file
  con <- ifelse (compressed, 
                 gzfile(description = file, open = "w"),
                 file(description = file, open = "w"))
  
  writeLines(text = lines, con = con)
  base::close(con)
  
  invisible(file)
  
}





