####################################################################################
## aligned spectra
parsePeakAbundanceMatrix <- function(filePeakMatrix, 
                                     doPrecursorDeisotoping, 
                                     mzDeviationInPPM_precursorDeisotoping, 
                                     mzDeviationAbsolute_precursorDeisotoping, 
                                     maximumRtDifference, 
                                     progress=FALSE)
{
  ## read file
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = paste("Parsing MS1 file content...", sep = "")) else print(paste("Parsing MS1 file content...", sep = ""))
  
  dataFrameAll <- read.table(filePeakMatrix, header=FALSE, sep = "\t", as.is=TRUE, quote = "\"", check.names = FALSE, comment.char = "")
  dataFrameAll1 <- read.table(filePeakMatrix, header=TRUE, sep = "\t", as.is=TRUE, quote = "\"", check.names = FALSE, comment.char = "")
  
  oldFormat <- max(which(dataFrameAll[1:5, 1] == "")) == 3
  header_rowNumber <- ifelse(test = oldFormat, yes = 4, no = 5)
  dataFrameHeader <- dataFrameAll[1:header_rowNumber, ]
  
  `%notin%` <- Negate(`%in%`)
  TR<- c("Reference RT","Reference m/z","Comment",
         "Manually modified for quantification",
         "Total score","RT similarity","Average","Stdev")
  IN2<-which(unname(dataFrameHeader[4,]) %notin% TR)
  dataFrameHeader1 <-dataFrameHeader[IN2]
  
  dataFrame <- dataFrameAll[(header_rowNumber + 1):nrow(dataFrameAll), ]
  dataFrame1 <- dataFrame[IN2]
  
  colnames(dataFrame) <- dataFrameHeader[header_rowNumber, ]
  colnames(dataFrame1) <- dataFrameHeader1[header_rowNumber, ]
  
  numberOfPrecursors <- nrow(dataFrame1)
  numberOfPrecursorsPrior <- numberOfPrecursors
  
  columnIndexEndOfAnnotation <- max(match(x = "Class", 
                                          table = dataFrameHeader1[1, ]), 
                                    na.rm = TRUE)
  
  if(ncol(dataFrame1) > columnIndexEndOfAnnotation){
    dataColumnStartEndIndeces <- c(columnIndexEndOfAnnotation + 1, ncol(dataFrame1))
    numberOfDataColumns <- dataColumnStartEndIndeces[[2]] - dataColumnStartEndIndeces[[1]] + 1
    dataColumnNames <- colnames(dataFrame1)[dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
    
    sampleClass          <- dataFrameHeader1[1, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
    sampleType           <- dataFrameHeader1[2, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
    sampleInjectionOrder <- dataFrameHeader1[3, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
    batchID              <- NULL
    if(!oldFormat)
      batchID            <- dataFrameHeader1[4, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
  } else {
    dataColumnStartEndIndeces <- NULL
    numberOfDataColumns <- 0
    dataColumnNames <- NULL
    
    sampleClass          <- NULL
    sampleType           <- NULL
    sampleInjectionOrder <- NULL
    batchID              <- NULL
  }
  
  commaNumbers <- sum(grepl(x = dataFrame1$"Average Mz", pattern = "^(\\d+,\\d+$)|(^\\d+$)"))
  decimalSeparatorIsComma <- commaNumbers == nrow(dataFrame1)
  if(decimalSeparatorIsComma){
    if(!is.null(dataFrame1$"Average Rt(min)"))     dataFrame1$"Average Rt(min)"     <- gsub(x = gsub(x = dataFrame1$"Average Rt(min)", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Average Mz"))          dataFrame1$"Average Mz"          <- gsub(x = gsub(x = dataFrame1$"Average Mz", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Fill %"))              dataFrame1$"Fill %"              <- gsub(x = gsub(x = dataFrame1$"Fill %", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Dot product"))         dataFrame1$"Dot product"         <- gsub(x = gsub(x = dataFrame1$"Dot product", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Reverse dot product")) dataFrame1$"Reverse dot product" <- gsub(x = gsub(x = dataFrame1$"Reverse dot product", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Fragment presence %")) dataFrame1$"Fragment presence %" <- gsub(x = gsub(x = dataFrame1$"Fragment presence %", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    
    ## replace -1 by 0
    if(numberOfDataColumns > 0) {
      for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
        dataFrame1[ , colIdx] <- gsub(x = gsub(x = dataFrame1[ , colIdx], pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
      }
    }
  }
  
  ###################
  ## column formats
  if(!is.null(dataFrame1$"Average Rt(min)"))     dataFrame1$"Average Rt(min)"     <- as.numeric(dataFrame1$"Average Rt(min)")
  if(!is.null(dataFrame1$"Average Mz"))          dataFrame1$"Average Mz"          <- as.numeric(dataFrame1$"Average Mz")
  if(!is.null(dataFrame1$"Fill %"))              dataFrame1$"Fill %"              <- as.numeric(dataFrame1$"Fill %")
  if(!is.null(dataFrame1$"MS/MS included"))      dataFrame1$"MS/MS included"      <- as.logical(dataFrame1$"MS/MS included")
  if(!is.null(dataFrame1$"Dot product"))         dataFrame1$"Dot product"         <- as.numeric(dataFrame1$"Dot product")
  if(!is.null(dataFrame1$"Reverse dot product")) dataFrame1$"Reverse dot product" <- as.numeric(dataFrame1$"Reverse dot product")
  if(!is.null(dataFrame1$"Fragment presence %")) dataFrame1$"Fragment presence %" <- as.numeric(dataFrame1$"Fragment presence %")
  
  #####################
  ## sorted by m/z (needed for deisotoping)
  if(!is.null(dataFrame1$"Average Mz"))
    dataFrame1 <- dataFrame1[order(dataFrame1$"Average Mz"), ]
  
  ## replace -1 by 0
  if(numberOfDataColumns > 0){
    for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
      dataFrame1[ , colIdx] <- as.numeric(dataFrame1[ , colIdx])
      if(!is.na(sum(dataFrame1[,colIdx] == -1)))
        dataFrame1[(dataFrame1[,colIdx] == -1),colIdx] <- 0
    }
  }
  
  ## deisotoping
  numberOfRemovedIsotopePeaks <- 0
  if(doPrecursorDeisotoping & !is.null(dataFrame1$"Average Mz")){
    if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Precursor deisotoping...", sep = "")) else print(paste("Precursor deisotoping...", sep = ""))
    distance13Cminus12C <- 1.0033548378
    
    ## mark isotope precursors
    precursorsToRemove <- vector(mode = "logical", length = numberOfPrecursors)
    
    if(numberOfDataColumns > 0){
      intensities <- dataFrame1[ , dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
      medians <- apply(X = as.matrix(intensities), MARGIN = 1, FUN = median)
    }
    
    for(precursorIdx in seq_len(numberOfPrecursors)){
      if((precursorIdx %% (as.integer(numberOfPrecursors/10))) == 0)
        if(!is.na(progress))  if(progress)  incProgress(amount = 0.0, detail = paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = "")) else print(paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = ""))
      
      mzError <- dataFrame1$"Average Mz"[[precursorIdx]] * mzDeviationInPPM_precursorDeisotoping / 1000000
      mzError <- max(mzError, mzDeviationAbsolute_precursorDeisotoping)
      
      ## RT difference <= maximumRtDifference
      validPrecursorsInRt <- abs(dataFrame1$"Average Rt(min)"[[precursorIdx]] - dataFrame1$"Average Rt(min)"[-precursorIdx]) <= maximumRtDifference
      
      ## MZ difference around 1.0033548378 (first isotope) or 1.0033548378 * 2 (second isotope)
      validPrecursorsInMz1 <- abs((dataFrame1$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 1) - dataFrame1$"Average Mz"[-precursorIdx]) <= mzError
      validPrecursorsInMz2 <- abs((dataFrame1$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 2) - dataFrame1$"Average Mz"[-precursorIdx]) <= mzError
      validPrecursorsInMz <- validPrecursorsInMz1 | validPrecursorsInMz2
      
      ## intensity gets smaller in the isotope spectrum
      if(numberOfDataColumns > 0){
        validPrecursorsInIntensity <- (medians[-precursorIdx] - medians[[precursorIdx]]) > 0
      } else {
        validPrecursorsInIntensity <- TRUE
      }
      
      if(any(validPrecursorsInRt & validPrecursorsInMz & validPrecursorsInIntensity))
        precursorsToRemove[[precursorIdx]] <- TRUE
    }
    
    ## remove isotopes
    dataFrame1 <- dataFrame1[!precursorsToRemove, ]
    
    numberOfRemovedIsotopePeaks <- sum(precursorsToRemove)
    numberOfPrecursors <- nrow(dataFrame1)
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Boxing...", sep = "")) else print(paste("Boxing...", sep = ""))
  returnObj <- list()
  returnObj$dataFrame1 <- dataFrame1
  
  ## meta
  returnObj$oldFormat <- oldFormat
  returnObj$numberOfPrecursors <- numberOfPrecursors
  returnObj$dataColumnStartEndIndeces <- dataColumnStartEndIndeces
  returnObj$numberOfDataColumns <- numberOfDataColumns
  
  ## group anno
  returnObj$sampleClass          <- sampleClass
  returnObj$sampleType           <- sampleType
  returnObj$sampleInjectionOrder <- sampleInjectionOrder
  returnObj$batchID              <- batchID
  
  ## misc
  returnObj$numberOfPrecursorsPrior <- numberOfPrecursorsPrior
  returnObj$numberOfRemovedIsotopePeaks <- numberOfRemovedIsotopePeaks
  
  return (returnObj)
}


####################################################################################
## parse MS/MS spectra
parseMSP <- function(fileSpectra, 
                     minimumIntensityOfMaximalMS2peak, 
                     minimumProportionOfMS2peaks, 
                     neutralLossesPrecursorToFragments, 
                     neutralLossesFragmentsToFragments, 
                     progress = FALSE){
  fileLines <- readLines(con = fileSpectra)
  
  returnObj <- parseMSP_chunk(fileLines, 
                              minimumIntensityOfMaximalMS2peak, 
                              minimumProportionOfMS2peaks, 
                              neutralLossesPrecursorToFragments, 
                              neutralLossesFragmentsToFragments, 
                              progress=FALSE)
  returnObj$fileSpectra <- fileSpectra
  
  return(returnObj)
}

parseMSP_big <- function(fileSpectra, 
                         minimumIntensityOfMaximalMS2peak, 
                         minimumProportionOfMS2peaks, 
                         neutralLossesPrecursorToFragments, 
                         neutralLossesFragmentsToFragments, 
                         progress = FALSE){
  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Read file") else print("MS/MS file: Read file")
  
  fileLines <- readLines(con = fileSpectra)
  
  if(progress)  incProgress(amount = 0, detail = "Decompose file to chunks") else print("Decompose file to chunks")
  ## entry lines in file
  isName	    <- grepl(pattern = "^Name:",			 x = fileLines)
  isNAME      <- grepl(pattern = "^NAME:",			 x = fileLines)
  isBeginIons	<- grepl(pattern = "^BEGIN IONS$", x = fileLines)
  
  ## entry line intervals in file
  entryBorders   <- c(which(isName | isNAME | isBeginIons), length(fileLines)+1)
  entryIntervals <- matrix(data = unlist(lapply(X = seq_len(length(entryBorders) - 1), FUN = function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
  
  rm(isName, isNAME, isBeginIons)
  
  ## split in chunks
  maximumNumberOfLinesPerChunk <- as.integer(1E6)
  fileLineChunks <- list()
  while(length(fileLines) > 0){
    numberOfEntries <- max(which(entryIntervals[2, ] <= maximumNumberOfLinesPerChunk))
    numberOfFileLinesInChunk <- entryIntervals[2, numberOfEntries]
    
    fileLinesInChunk <- fileLines[1:numberOfFileLinesInChunk]
    fileLineChunks[[length(fileLineChunks) + 1]] <- fileLinesInChunk
    
    ## remaining stuff
    fileLines      <- fileLines[-(1:numberOfFileLinesInChunk)]
    entryIntervals <- entryIntervals[, -(1:numberOfEntries)]
    if(ncol(entryIntervals) > 0)
      entryIntervals <- entryIntervals - min(entryIntervals) + 1
  }
  
  numberOfChunks <- length(fileLineChunks)
  fileLineCountOfChunks  <- unlist(lapply(X = fileLineChunks, FUN = length))
  fileLineOffsetOfChunks <- c(0, sapply(X = seq_along(fileLineCountOfChunks)[-1], FUN = function(x){sum(fileLineCountOfChunks[1:(x-1)])}))
  
  if(progress)  incProgress(amount = 0, detail = paste("Parse", numberOfChunks, "chunks")) else print(paste("Parse", numberOfChunks, "chunks"))
  ## merge chunks
  returnObj <- list()
  returnObj$fileSpectra <- fileSpectra
  returnObj$spectraList <- list()
  returnObj$numberOfSpectra <- 0
  returnObj$numberOfMS2PeaksOriginal <- 0
  returnObj$numberOfMS2PeaksWithNeutralLosses <- 0
  returnObj$numberOfMS2PeaksAboveThreshold <- 0
  returnObj$numberOfMS2PeaksBelowThreshold <- 0
  returnObj$precursorMz <- vector(mode = "numeric")
  returnObj$precursorRt <- vector(mode = "numeric")
  
  for(chunkIdx in seq_len(numberOfChunks)){
    if(progress)  incProgress(amount = 0, detail = paste("Chunk", chunkIdx, "/", numberOfChunks)) else print(paste("Chunk", chunkIdx, "/", numberOfChunks))
    
    fileLines  <- fileLineChunks[[1]]
    fileLineChunks[[1]] <- NULL
    
    returnObj2 <- parseMSP_chunk(
      fileLines = fileLines, 
      minimumIntensityOfMaximalMS2peak = minimumIntensityOfMaximalMS2peak, 
      minimumProportionOfMS2peaks = minimumProportionOfMS2peaks, 
      neutralLossesPrecursorToFragments = neutralLossesPrecursorToFragments, 
      neutralLossesFragmentsToFragments = neutralLossesFragmentsToFragments, 
      offset = fileLineOffsetOfChunks[[chunkIdx]], progress = NA
    )
    
    returnObj$spectraList                       <- c(returnObj$spectraList,                        returnObj2$spectraList)
    returnObj$numberOfSpectra                   <-   returnObj$numberOfSpectra                   + returnObj2$numberOfSpectra
    returnObj$numberOfSpectraOriginal           <-   returnObj$numberOfSpectraOriginal           + returnObj2$numberOfSpectraOriginal
    returnObj$numberOfMS2PeaksOriginal          <-   returnObj$numberOfMS2PeaksOriginal          + returnObj2$numberOfMS2PeaksOriginal
    returnObj$numberOfMS2PeaksWithNeutralLosses <-   returnObj$numberOfMS2PeaksWithNeutralLosses + returnObj2$numberOfMS2PeaksWithNeutralLosses
    returnObj$numberOfMS2PeaksAboveThreshold    <-   returnObj$numberOfMS2PeaksAboveThreshold    + returnObj2$numberOfMS2PeaksAboveThreshold
    returnObj$numberOfMS2PeaksBelowThreshold    <-   returnObj$numberOfMS2PeaksBelowThreshold    + returnObj2$numberOfMS2PeaksBelowThreshold
    returnObj$numberOfTooHeavyFragments         <-   returnObj$numberOfTooHeavyFragments         + returnObj2$numberOfTooHeavyFragments
    returnObj$numberOfSpectraDiscardedDueToNoPeaks      <-   returnObj$numberOfSpectraDiscardedDueToNoPeaks      + returnObj2$numberOfSpectraDiscardedDueToNoPeaks
    returnObj$numberOfSpectraDiscardedDueToMaxIntensity <-   returnObj$numberOfSpectraDiscardedDueToMaxIntensity + returnObj2$numberOfSpectraDiscardedDueToMaxIntensity
    returnObj$numberOfSpectraDiscardedDueToTooHeavy     <-   returnObj$numberOfSpectraDiscardedDueToTooHeavy     + returnObj2$numberOfSpectraDiscardedDueToTooHeavy
    returnObj$precursorMz                       <- c(returnObj$precursorMz,                        returnObj2$precursorMz)
    returnObj$precursorRt                       <- c(returnObj$precursorRt,                        returnObj2$precursorRt)
  }
  
  return(returnObj)
}

#####################
parseMSP_chunk <- function(fileLines, 
                           minimumIntensityOfMaximalMS2peak, 
                           minimumProportionOfMS2peaks, 
                           neutralLossesPrecursorToFragments, 
                           neutralLossesFragmentsToFragments, 
                           offset = 0, 
                           progress = FALSE)
{
  
  ## LC-MS/MS entry:
  ## NAME: Unknown
  ## RETENTIONTIME: 3.215358
  ## PRECURSORMZ: 78.91963
  ## METABOLITENAME: 
  ## ADDUCTIONNAME: [M-H]-
  ## Num Peaks: 2
  ## 76.97093  754
  ## 76.98951  754
  ## 
  ## GC-MS additional properties:
  ## SCANNUMBER: 518
  ## MODELION: 59
  ## MODELIONHEIGHT: 924
  ## MODELIONAREA: 924
  ## INTEGRATEDHEIGHT: 924
  ## INTEGRATEDAREA: 924
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Read file") else print("MS/MS file: Read file")
  
  #fileLines <- readLines(con = fileSpectra)
  numberOfFileLines <- length(fileLines)
  
  numberOfMS2PeaksOriginal <<- 0
  numberOfMS2PeaksWithNeutralLosses <<- 0
  numberOfTooHeavyFragments <<- 0
  numberOfMS2PeaksAboveThreshold <<- 0
  numberOfMS2PeaksBelowThreshold <<- 0
  numberOfSpectraDiscardedDueToNoPeaks      <<- 0
  numberOfSpectraDiscardedDueToMaxIntensity <<- 0
  numberOfSpectraDiscardedDueToTooHeavy <<- 0
  
  ## start with empty lines or not?
  endOfRecord <- TRUE
  if(numberOfFileLines > 0)
    if(nchar(trimws(fileLines[[1]])) > 0)
      endOfRecord <- FALSE
  
  ## check for pattern
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Parse") else print("MS/MS file: Parse")
  isBI  	<- grepl(pattern = "^BEGIN IONS$",                    		x = fileLines)
  isName	<- grepl(pattern = "(^Name:)|(^NAME:)",               		x = fileLines)
  isNAme	<- grepl(pattern = "^NAME=",                           		x = fileLines)
  isNAME	<- grepl(pattern = "^TITLE=",                          		x = fileLines)
  isRT	  <- grepl(pattern = "^RETENTIONTIME:",						          x = fileLines)
  isRt	  <- grepl(pattern = "^retention time:",			  	          x = fileLines)
  isRts	  <- grepl(pattern = "^RTINSECONDS=",	     		  	          x = fileLines)
  isMZ	  <- grepl(pattern = "^PRECURSORMZ:",							          x = fileLines)
  isMz	  <- grepl(pattern = "^precursor m/z:",						          x = fileLines)
  isTotEm	<- grepl(pattern = "^total exact mass:",					        x = fileLines)
  isEMass	<- grepl(pattern = "(^exact mass:)|(^EXACT_MASS:)",       x = fileLines)
  isPMass	<- grepl(pattern = "^PEPMASS=",   							          x = fileLines)
  isE2Mass<- grepl(pattern = "^EXACTMASS=", 							          x = fileLines)
  isMetN	<- grepl(pattern = "^METABOLITENAME:",						        x = fileLines)
  isAddN	<- grepl(pattern = "(^ADDUCTIONNAME:)|(^Adductionname:)",	x = fileLines)
  isScanN	<- grepl(pattern = "^SCANNUMBER:",							          x = fileLines)
  isMIon	<- grepl(pattern = "^MODELION:",							            x = fileLines)
  isPrety	<- grepl(pattern = "^precursor type:",						        x = fileLines)
  isPretyp	<- grepl(pattern = "^PRECURSORTYPE:",						          x = fileLines)
  isNumP	<- grepl(pattern = "^Num Peaks:",							            x = fileLines)
  isPeak	<- grepl(pattern = "^\\d+(((\\.)|(,))\\d+)?[ \t]\\d+(((\\.)|(,))\\d+)?$",	x = fileLines)
  isCoCl	<- grepl(pattern = "^compound class:",						        x = fileLines)
  isInty	<- grepl(pattern = "^instrument type:",						        x = fileLines)
  isIntyp	<- grepl(pattern = "^INSTRUMENTTYPE:",						        x = fileLines)
  isIntype<- grepl(pattern = "^SOURCE_INSTRUMENT=",					        x = fileLines)
  isInt  	<- grepl(pattern = "^INSTRUMENT:",						            x = fileLines)
  isInchi	<- grepl(pattern = "(^InChI:)|(^InChI=)|(^INCHI=)|(^INCHI:)", x = fileLines)
  isInchiKey	<- grepl(pattern = "(^InChIKey:)|(^INCHIKEY:)|(^InChIKey=)|(^INCHIKEY=)|(^INCHIAUX=)",	      x = fileLines)
  isSmiles  	<- grepl(pattern = "(^SMILES:)|(^SMILES=)",		        x = fileLines)
  
  someStrings <- trimws(c(
    substring(text = fileLines[isMZ], first = nchar("RETENTIONTIME:") + 1), 
    substring(text = fileLines[isMZ], first = nchar("PRECURSORMZ:") + 1), 
    substring(text = fileLines[isMz], first = nchar("precursor m/z:") + 1)
  ))
  decimalDelimiterIsComma <- ifelse(test = length(someStrings) > 0, yes = all(grepl(x = someStrings, pattern = "^(\\d+,\\d+$)|(^\\d+$)")), no = FALSE)
  if(decimalDelimiterIsComma){
    fileLines_02 <- gsub(pattern = ",", replacement = ".", x = gsub(pattern = "\\.", replacement = "", x = fileLines))
  } else {
    fileLines_02 <- fileLines
  }
  
  ## extract
  suppressWarnings({
    parsedName	 <- 						      trimws(substring(text = fileLines, first = nchar("Name:") + 1))
    parsedNAme	 <- 						      trimws(substring(text = fileLines, first = nchar("NAME=") + 1))
    parsedNAME	 <- 						      trimws(substring(text = fileLines, first = nchar("TITLE=") + 1))
    parsedRT  	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("RETENTIONTIME:") + 1)))
    parsedRt     <- as.numeric(unlist(lapply(X = strsplit(x =   trimws(substring(text = fileLines_02, first = nchar("retention time:") + 1)), split = " "), FUN = function(x){if(length(x)==0) return(NA) else return(x[[1]])})))
    parsedRts  	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("RTINSECONDS=") + 1)))
    parsedMZ	   <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("PRECURSORMZ:") + 1)))
    parsedMz	   <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("precursor m/z:") + 1)))
    parsedTotEm  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("total exact mass:") + 1)))
    parsedEMass  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("exact mass:") + 1)))
    parsedE2Mass <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("EXACTMASS=") + 1)))
    parsedPMass  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("PEPMASS=") + 1)))
    parsedMetN	 <- 						      trimws(substring(text = fileLines, first = nchar("METABOLITENAME:") + 1))
    parsedAddN	 <- 						      trimws(substring(text = fileLines, first = nchar("Adductionname:") + 1))
    parsedScanN  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("SCANNUMBER:") + 1)))
    parsedMIon	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("MODELION:") + 1)))
    parsedPrety  <- 						      trimws(substring(text = fileLines, first = nchar("precursor type:") + 1))
    parsedPretyp <- 						      trimws(substring(text = fileLines, first = nchar("PRECURSORTYPE:") + 1))
    parsedNumP	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("Num Peaks:") + 1)))
    
    parsedTokensTmp <- strsplit(x = trimws(fileLines_02), split = "[ \t]")
    #parsedms2Peaks_mz <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){head(x = x, n = 1)})))
    parsedms2Peaks_mz  <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){if(length(x)<1) return(NA) else return(x[[1]])})))
    parsedms2Peaks_int <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){if(length(x)<2) return(NA) else return(x[[2]])})))
    
    parsedCoCl	 <- 						      trimws(substring(text = fileLines, first = nchar("compound class:") + 1))
    parsedInty	 <- 						      trimws(substring(text = fileLines, first = nchar("instrument type:") + 1))
    parsedIntype <- 						      trimws(substring(text = fileLines, first = nchar("SOURCE_INSTRUMENT=") + 1))
    parsedIntyp  <- 						      trimws(substring(text = fileLines, first = nchar("INSTRUMENTTYPE:") + 1))
    parsedInt    <- 						      trimws(substring(text = fileLines, first = nchar("INSTRUMENT:") + 1))
    parsedInchi  <- 						      trimws(substring(text = fileLines, first = nchar("InChI:") + 1))
    parsedInchiKey  <- 						    trimws(substring(text = fileLines, first = nchar("InChIKey:") + 1))
    parsedSmiles    <- 						    trimws(substring(text = fileLines, first = nchar("SMILES:") + 1))
  })
  
  ## entry line intervals in file
  entryBorders   <- c(which(isName | isNAME | isBI), length(fileLines)+1)
  entryIntervals <- matrix(data = unlist(lapply(X = seq_len(length(entryBorders) - 1), FUN = function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
  
  numberOfSpectraOriginal <- length(entryBorders)
  
  ## do it
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Assemble spectra") else print("MS/MS file: Assemble spectra")
  suppressWarnings(
    spectraList <- apply(X = entryIntervals, MARGIN = 2, FUN = function(x){
      #print(x)
      
      fileLines2      <- fileLines 	 		[x[[1]]:x[[2]]]
      fileLines_022   <- fileLines_02 	[x[[1]]:x[[2]]]
      parsedName2	 		<- parsedName	 		[x[[1]]:x[[2]]]
      parsedNAme2	 		<- parsedNAme	 		[x[[1]]:x[[2]]]
      parsedNAME2	 		<- parsedNAME	 		[x[[1]]:x[[2]]]
      parsedRT2  	 		<- parsedRT  	 		[x[[1]]:x[[2]]]
      parsedRt2     	<- parsedRt     	[x[[1]]:x[[2]]]
      parsedRts2     	<- parsedRts     	[x[[1]]:x[[2]]]
      parsedMZ2	 		  <- parsedMZ	 		  [x[[1]]:x[[2]]]
      parsedMz2	 		  <- parsedMz	 		  [x[[1]]:x[[2]]]
      parsedTotEm2  	<- parsedTotEm  	[x[[1]]:x[[2]]]
      parsedEMass2  	<- parsedEMass  	[x[[1]]:x[[2]]]
      parsedE2Mass2  	<- parsedE2Mass  	[x[[1]]:x[[2]]]
      parsedPMass2  	<- parsedPMass  	[x[[1]]:x[[2]]]
      parsedMetN2	 		<- parsedMetN	 		[x[[1]]:x[[2]]]
      parsedAddN2	 		<- parsedAddN	 		[x[[1]]:x[[2]]]
      parsedScanN2  	<- parsedScanN  	[x[[1]]:x[[2]]]
      parsedMIon2	 		<- parsedMIon	 		[x[[1]]:x[[2]]]
      parsedPrety2  	<- parsedPrety  	[x[[1]]:x[[2]]]
      parsedPretyp2  	<- parsedPretyp  	[x[[1]]:x[[2]]]
      parsedNumP2	 		<- parsedNumP	 		[x[[1]]:x[[2]]]
      parsedms2Peaks_mz2	<- parsedms2Peaks_mz	[x[[1]]:x[[2]]]
      parsedms2Peaks_int2	<- parsedms2Peaks_int	[x[[1]]:x[[2]]]
      parsedCoCl2	 		<- parsedCoCl	 		[x[[1]]:x[[2]]]
      parsedInty2	 		<- parsedInty	 		[x[[1]]:x[[2]]]
      parsedIntype2		<- parsedIntype		[x[[1]]:x[[2]]]
      parsedIntyp2		<- parsedIntyp		[x[[1]]:x[[2]]]
      parsedInt2	  	<- parsedInt  		[x[[1]]:x[[2]]]
      parsedInchi2  	<- parsedInchi  	[x[[1]]:x[[2]]]
      parsedInchiKey2 <- parsedInchiKey [x[[1]]:x[[2]]]
      parsedSmiles2  	<- parsedSmiles  	[x[[1]]:x[[2]]]
      
      isName2	 		<- isName	 		[x[[1]]:x[[2]]]
      isNAme2	 		<- isName	 		[x[[1]]:x[[2]]]
      isNAME2	 		<- isNAME	 		[x[[1]]:x[[2]]]
      isRT2  	 		<- isRT  	 		[x[[1]]:x[[2]]]
      isRt2     	<- isRt     	[x[[1]]:x[[2]]]
      isRts2     	<- isRts     	[x[[1]]:x[[2]]]
      isMZ2	 		  <- isMZ	 		  [x[[1]]:x[[2]]]
      isMz2	 		  <- isMz	 		  [x[[1]]:x[[2]]]
      isTotEm2  	<- isTotEm  	[x[[1]]:x[[2]]]
      isEMass2  	<- isEMass  	[x[[1]]:x[[2]]]
      isE2Mass2  	<- isE2Mass  	[x[[1]]:x[[2]]]
      isPMass2  	<- isPMass  	[x[[1]]:x[[2]]]
      isMetN2	 		<- isMetN	 		[x[[1]]:x[[2]]]
      isAddN2	 		<- isAddN	 		[x[[1]]:x[[2]]]
      isScanN2  	<- isScanN  	[x[[1]]:x[[2]]]
      isMIon2	 		<- isMIon	 		[x[[1]]:x[[2]]]
      isPrety2  	<- isPrety  	[x[[1]]:x[[2]]]
      isPretyp2  	<- isPretyp  	[x[[1]]:x[[2]]]
      isNumP2	 		<- isNumP	 		[x[[1]]:x[[2]]]
      isPeak2	    <- isPeak	    [x[[1]]:x[[2]]]
      isCoCl2	 		<- isCoCl	 		[x[[1]]:x[[2]]]
      isInty2	 		<- isInty	 		[x[[1]]:x[[2]]]
      isIntype2		<- isIntype		[x[[1]]:x[[2]]]
      isIntyp2		<- isIntyp		[x[[1]]:x[[2]]]
      isInt2  		<- isInt  		[x[[1]]:x[[2]]]
      isInchi2  	<- isInchi  	[x[[1]]:x[[2]]]
      isInchiKey2 <- isInchiKey [x[[1]]:x[[2]]]
      isSmiles2  	<- isSmiles  	[x[[1]]:x[[2]]]
      
      name <- NULL
      ms1Int <- NA
      rt <- NULL
      mz <- NULL
      metName <- "Unknown"
      adduct <- "Unknown"
      scanNumber <- NA
      quantMass <- NA
      peakNumber <- NA
      ms2Peaks_mz  <- vector(mode = "numeric")
      ms2Peaks_int <- vector(mode = "numeric")
      compoundClass <- "Unknown"
      instrumentType <- "Unknown"
      inchi <- ""
      inchiKey <- ""
      smiles <- ""
      
      if(any(isName2))
        ## name
        name        <- parsedName2	[[which(isName2)[[1]]]]
      if(any(isNAme2))
        ## name
        name        <- parsedNAme2	[[which(isNAme2)[[1]]]]
      if(any(isNAME2))
        ## name
        name        <- parsedNAME2	[[which(isNAME2)[[1]]]]
      if(any(isRT2))
        ## retention time
        rt          <- parsedRT2	[[which(isRT2)[[1]]]]
      if(any(isRt2) & any(is.null(rt), is.na(rt)))
        rt          <- parsedRt2	[[which(isRt2)[[1]]]]
      if(any(isRts2) & any(is.null(rt), is.na(rt)))
        rt          <- parsedRts2	[[which(isRts2)[[1]]]]
      if(any(isMZ2))
        ## precursor m/z
        mz          <- parsedMZ2	[[which(isMZ2)[[1]]]]
      if(any(isMz2)    & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0)))
        mz          <- parsedMz2	[[which(isMz2)[[1]]]]
      if(any(isTotEm2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0)) & any(isPrety2)){
        mzTmp <- parsedTotEm2[[which(isTotEm2)[[1]]]]
        pit   <- parsedPrety2[[which(isPrety2)[[1]]]]
        switch(pit,
               "[M-H]-" = { mz <- mzTmp -  1.008  },
               "[M-H]"  = { mz <- mzTmp -  1.008  },
               "[M+H]+" = { mz <- mzTmp +  1.008  },
               "[M+H]"  = { mz <- mzTmp +  1.008  },
               "[M+Na]+"= { mz <- mzTmp + 22.9898 },
               "[M+Na]" = { mz <- mzTmp + 22.9898 }#,
        )
      }
      if(any(isEMass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0))){
        mz <- parsedEMass2[[which(isEMass2)[[1]]]]
      }
      if(any(isPMass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0))){
        if(!is.na(parsedPMass2[[which(isPMass2)[[1]]]])){
          mz  <- parsedPMass2[[which(isPMass2)[[1]]]]
        } else {
          tmp <- as.numeric(strsplit(x = trimws(substring(text = fileLines_022[[which(isPMass2)[[1]]]], first = nchar("PEPMASS=") + 1)), split = "[\t ]")[[1]])
          mz  <- tmp[[1]]
          if(length(tmp) > 1)
            ms1Int <- tmp[[2]]
        }
      }
      if(any(isE2Mass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0))){
        mz <- parsedE2Mass2[[which(isE2Mass2)[[1]]]]
      }
      if(any(isMetN2))
        metName     <- parsedMetN2	[[which(isMetN2)[[1]]]]
      if(any(isAddN2)  & adduct == "Unknown")
        ## adduct
        adduct      <- parsedAddN2	[[which(isAddN2)[[1]]]]
      if(any(isPrety2) & adduct == "Unknown")
        adduct      <- parsedPrety2[[which(isPrety2)[[1]]]]
      if(any(isPretyp2) & adduct == "Unknown")
        adduct      <- parsedPretyp2[[which(isPretyp2)[[1]]]]
      if(any(isScanN2))
        scanNumber  <- parsedScanN2[[which(isScanN2)[[1]]]]
      if(any(isMIon2))
        quantMass  <- parsedMIon2	[[which(isMIon2)[[1]]]]
      if(any(isNumP2))
        ## #peaks
        peakNumber  <- parsedNumP2	[[which(isNumP2)[[1]]]]
      if(any(isPeak2)){
        ## MS2 peaks: "178.88669\t230"
        ms2Peaks_mz  <- parsedms2Peaks_mz2 [which(isPeak2)]
        ms2Peaks_int <- parsedms2Peaks_int2[which(isPeak2)]
      }
      if(any(isCoCl2))
        ## compound class
        compoundClass  <- parsedCoCl2	[[which(isCoCl2)[[1]]]]
      if(any(isInty2))
        ## instrument type
        instrumentType  <- parsedInty2	[[which(isInty2)[[1]]]]
      if(any(isIntype2))
        ## instrument type
        instrumentType  <- parsedIntype2	[[which(isIntype2)[[1]]]]
      if(any(isIntyp2))
        ## instrument type
        instrumentType  <- parsedIntyp2	[[which(isIntyp2)[[1]]]]
      if(all(any(isInt2), instrumentType %in% c("Unknown", "NA")))
        ## instrument
        instrumentType  <- parsedInt2	[[which(isInt2)[[1]]]]
      if(any(isInchi2))
        ## structure
        inchi  <- parsedInchi2 [[which(isInchi2)[[1]]]]
      if(any(isInchiKey2))
        ## structure
        inchiKey  <- parsedInchiKey2 [[which(isInchiKey2)[[1]]]]
      if(any(isSmiles2))
        ## structure
        smiles  <- parsedSmiles2 [[which(isSmiles2)[[1]]]]
      ## end of parsing
      
      if(is.null(rt))
        rt <- 0
      
      ms2Peaks_mz_original  <- ms2Peaks_mz
      ms2Peaks_int_original <- ms2Peaks_int
      
      numberOfMS2PeaksOriginal <<- numberOfMS2PeaksOriginal + length(ms2Peaks_mz)
      if(length(ms2Peaks_mz) == 0)
        numberOfSpectraDiscardedDueToNoPeaks <<- numberOfSpectraDiscardedDueToNoPeaks + 1
      
      ###################################################################
      ## filter fragments with mass greater than precursor
      numberOfTooHeavyFragmentsHere <- 0
      if(all(!is.null(mz), !is.na(mz))){
        tooHeavy <- ms2Peaks_mz > mz
        ms2Peaks_mz  <- ms2Peaks_mz [!tooHeavy]
        ms2Peaks_int <- ms2Peaks_int[!tooHeavy]
        numberOfTooHeavyFragmentsHere <- sum(tooHeavy)
        
        if(length(ms2Peaks_mz) == 0 & numberOfTooHeavyFragmentsHere > 0)
          numberOfSpectraDiscardedDueToTooHeavy <<- numberOfSpectraDiscardedDueToTooHeavy + 1
      }
      numberOfTooHeavyFragments <<- numberOfTooHeavyFragments + numberOfTooHeavyFragmentsHere
      
      ###################################################################
      ## filter for ms2 peak intensity
      peakNumber <- length(ms2Peaks_mz)
      if(peakNumber > 0){
        ###################################################################
        ## filter for ms2 peak intensity relative to maximum peak intensity in spectrum
        maximumIntensity <- max(ms2Peaks_int)
        if(maximumIntensity >= minimumIntensityOfMaximalMS2peak){
          ## spectrum is considered
          intensityThreshold <- maximumIntensity * minimumProportionOfMS2peaks
          fragmentsAboveThreshold <- ms2Peaks_int >= intensityThreshold
          
          ms2Peaks_mz  <- ms2Peaks_mz [fragmentsAboveThreshold]
          ms2Peaks_int <- ms2Peaks_int[fragmentsAboveThreshold]
          numberOfMS2PeaksAboveThreshold <<- numberOfMS2PeaksAboveThreshold + sum( fragmentsAboveThreshold)
          numberOfMS2PeaksBelowThreshold <<- numberOfMS2PeaksBelowThreshold + sum(!fragmentsAboveThreshold)
        } else {
          ## spectrum is not considered
          numberOfMS2PeaksBelowThreshold <<- numberOfMS2PeaksBelowThreshold + length(ms2Peaks_mz)
          numberOfSpectraDiscardedDueToMaxIntensity <<- numberOfSpectraDiscardedDueToMaxIntensity + 1
          ms2Peaks_mz  <- vector(mode = "numeric")
          ms2Peaks_int <- vector(mode = "numeric")
        }
      }
      
      ###################################################################
      ## normalize ms2 peaks to maximum = 1
      peakNumber <- length(ms2Peaks_mz)
      if(peakNumber > 0){
        max <- max(ms2Peaks_int)
        ms2Peaks_int <- ms2Peaks_int / max
      }
      
      ###################################################################
      ## add neutral losses
      if(peakNumber > 0){
        #################################
        ## neutral losses regarding the precursor
        if(all(!is.null(mz), !is.na(mz), neutralLossesPrecursorToFragments)){
          ms2PeaksNLPF_mz  <- ms2Peaks_mz - as.numeric(mz)
          ms2PeaksNLPF_int <- ms2Peaks_int
        } else {
          ms2PeaksNLPF_mz  <- vector(mode = "numeric")
          ms2PeaksNLPF_int <- vector(mode = "numeric")
        }
        #################################
        ## neutral losses amongst fragments
        if(neutralLossesFragmentsToFragments){
          m_mz  <- outer(X = ms2Peaks_mz,  Y = ms2Peaks_mz,  FUN = function(x,y){x-y})
          m_int <- outer(X = ms2Peaks_int, Y = ms2Peaks_int, FUN = function(x,y){(x+y) / 2})
          upper <- upper.tri(x = m_mz)
          ms2PeaksNLFF_mz  <- m_mz [upper]
          ms2PeaksNLFF_int <- m_int[upper]
        } else {
          ms2PeaksNLFF_mz  <- vector(mode = "numeric")
          ms2PeaksNLFF_int <- vector(mode = "numeric")
        }
        
        ms2Peaks_mz  <- c(ms2Peaks_mz,  ms2PeaksNLPF_mz,  ms2PeaksNLFF_mz)
        ms2Peaks_int <- c(ms2Peaks_int, ms2PeaksNLPF_int, ms2PeaksNLFF_int)
      }
      
      ###################################################################
      ## precursor mz
      #mz <- ifelse(test = !is.null(mz), yes = round(as.numeric(mz), digits = 4), no = ifelse(test = !is.null(scanNumber), yes = scanNumber, no = max(ms2Peaks_mz)))
      if(all(!is.null(mz), !is.na(mz))){
        mz <- round(as.numeric(mz), digits = 4)
      } else {
        if(!is.na(quantMass)){
          mz <- quantMass
        } else {
          if(!is.na(scanNumber)){
            mz <- scanNumber
          } else {
            mz <- max(ms2Peaks_mz)
          }
        }
      }
      
      ###################################################################
      ## string representation of spectrum
      spectrumString <- paste(ms2Peaks_mz_original, ms2Peaks_int_original, sep = " ", collapse = ";")
      
      ###################################################################
      ## built ms set
      spectrumItem <- list(
        name = name,
        ms1Int = ms1Int,
        #rt = round(as.numeric(rt), digits = 2),
        rt = rt,
        mz = mz,
        metName = metName,
        adduct = adduct,
        quantMass = quantMass,
        compoundClass = compoundClass,
        instrumentType = instrumentType,
        inchi = inchi,
        inchiKey = inchiKey,
        smiles = smiles,
        #peakNumber = as.numeric(peakNumber),
        peakNumber = length(ms2Peaks_mz),
        ms2Peaks_mz  = ms2Peaks_mz,
        ms2Peaks_int = ms2Peaks_int,
        spectrumString = spectrumString,
        entryInterval = x + offset
      )
      if(spectrumItem$peakNumber > 0){
        ## add
        numberOfMS2PeaksWithNeutralLosses <<- numberOfMS2PeaksWithNeutralLosses + spectrumItem$peakNumber
        return(spectrumItem)
      } else
        return(NULL)
    })
  )
  
  rm(
    isName,
    isNAME,
    isRT,
    isRt,
    isMZ,
    isMz,
    isTotEm,
    isEMass,
    isE2Mass,
    isPMass,
    isMetN,
    isAddN,
    isScanN,
    isMIon,
    isPrety,
    isNumP,
    isPeak,
    isCoCl,
    isInty,
    isInchi,
    isInchiKey,
    isSmiles,
    parsedName,
    parsedNAME,
    parsedRT,
    parsedRt,
    parsedMZ,
    parsedMz,
    parsedTotEm,
    parsedEMass,
    parsedE2Mass,
    parsedPMass,
    parsedMetN,
    parsedAddN,
    parsedScanN,
    parsedMIon,
    parsedPrety,
    parsedNumP,
    parsedTokensTmp,
    parsedms2Peaks_mz,
    parsedms2Peaks_int,
    parsedCoCl,
    parsedInty,
    parsedInchi,
    parsedInchiKey,
    parsedSmiles
  )
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "MS/MS file: Box") else print("MS/MS file: Box")
  
  ## remove NULL entries?
  spectraList[unlist(lapply(X = spectraList, FUN = is.null))] <- NULL
  
  numberOfSpectra <- length(spectraList)
  
  ## postprocess
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file postprocessing", sep = "")) else print(paste("MS/MS file postprocessing", sep = ""))
  
  precursorMz <- unlist(lapply(X = spectraList, FUN = function(x){
    if(is.null(x))  return(NULL)
    else            return(x$mz)
  }))
  precursorRt <- unlist(lapply(X = spectraList, FUN = function(x){
    if(is.null(x))  return(NULL)
    else            return(x$rt)
  }))
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file boxing", sep = "")) else print(paste("MS/MS file boxing", sep = ""))
  returnObj <- list()
  returnObj$fileSpectra <- NA
  returnObj$spectraList <- spectraList
  returnObj$numberOfSpectra <- numberOfSpectra
  returnObj$numberOfSpectraOriginal <- numberOfSpectraOriginal
  returnObj$numberOfMS2PeaksOriginal <- numberOfMS2PeaksOriginal
  returnObj$numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses
  returnObj$numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold
  returnObj$numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold
  returnObj$numberOfTooHeavyFragments <- numberOfTooHeavyFragments
  returnObj$numberOfSpectraDiscardedDueToNoPeaks <- numberOfSpectraDiscardedDueToNoPeaks
  returnObj$numberOfSpectraDiscardedDueToMaxIntensity <- numberOfSpectraDiscardedDueToMaxIntensity
  returnObj$numberOfSpectraDiscardedDueToTooHeavy <- numberOfSpectraDiscardedDueToTooHeavy
  returnObj$precursorMz <- precursorMz
  returnObj$precursorRt <- precursorRt
  
  return(returnObj)
}

parseMSP_attributes <- function(fileSpectra, progress = FALSE, flexiblePeakList = FALSE, multiplePeaksPerLine = FALSE, includeIDasRecordSeparator=TRUE, includeNAMEasRecordSeparator=TRUE, includeTITLEasRecordSeparator=TRUE, returnEmptySpectra = FALSE){
  fileLines <- readLines(con = fileSpectra)
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Read file") else print("MS/MS file: Read file")
  
  numberOfFileLines <- length(fileLines)
  
  ## start with empty lines or not?
  endOfRecord <- TRUE
  if(numberOfFileLines > 0)
    if(nchar(trimws(fileLines[[1]])) > 0)
      endOfRecord <- FALSE
  
  ## check for pattern
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Parse") else print("MS/MS file: Parse")
  isID  	<- grepl(pattern = "^ID:",                            		x = fileLines)
  isBI  	<- grepl(pattern = "^BEGIN IONS$",                    		x = fileLines)
  isName	<- grepl(pattern = "(^Name:)|(^NAME:)",               		x = fileLines)
  isNAme	<- grepl(pattern = "^NAME=",                           		x = fileLines)
  isTITLE	<- grepl(pattern = "^TITLE=",                          		x = fileLines)
  isAccession	<- grepl(pattern = "^ACCESSION:",                          		x = fileLines)
  
  if(!includeIDasRecordSeparator) isID <- rep(x = F, times = length(isID))
  if(!includeNAMEasRecordSeparator){
    isName <- rep(x = F, times = length(isName))
    isNAme <- rep(x = F, times = length(isNAme))
  }
  if(!includeTITLEasRecordSeparator) isTITLE <- rep(x = F, times = length(isTITLE))
  
  numberSmall <- "\\d+(\\.\\d+)?"
  numberBig   <- "\\d+(\\.\\d+(E\\d+)?)?"
  mzValueRegEx <- "(\\d+(\\.\\d+)?)"
  intensityRegex <- "(\\d+((\\.\\d+)?([eE](-)?\\d+)?)?)"
  annotationRegex <- "\".+\""
  if(flexiblePeakList){
    isPeak	<- grepl(pattern = "^[ \t]*\\d+(\\.\\d+([eE](-)?\\d+)?)?([ \t]\\d+(\\.\\d+([eE](-)?\\d+)?)?)*[ \t]*$",	x = fileLines)
  } else {
    isPeak	<- grepl(pattern = paste("^[ \t]*", mzValueRegEx, "[ \t]", intensityRegex, "([ \t]+((", mzValueRegEx, "[ \t]", intensityRegex, ")|(", annotationRegex, ")))*[ \t]*$", sep = ""),	x = fileLines)
  }
  isEmpty	<- nchar(trimws(fileLines)) == 0
  
  tagVector   <- unlist(lapply(X = str_split(string = fileLines, pattern = "(:)|(=)"), FUN = function(x){x[[1]]}))
  valueVector <- trimws(substr(x = fileLines, start = nchar(tagVector) + 1 + 1, stop = nchar(fileLines)))
  
  ## entry line intervals in file
  entryBorders   <- c(which(isName | isNAme | isTITLE | isBI | isID | isAccession), length(fileLines)+1)
  entryIntervals <- matrix(data = unlist(lapply(X = seq_len(length(entryBorders) - 1), FUN = function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
  
  ## do it
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Assemble spectra") else print("MS/MS file: Assemble spectra")
  spectraList <- apply(X = entryIntervals, MARGIN = 2, FUN = function(x){
    #print(x)
    fileLines2   <- fileLines  [x[[1]]:x[[2]]]
    isPeak2	     <- isPeak	   [x[[1]]:x[[2]]]
    isEmpty2     <- isEmpty    [x[[1]]:x[[2]]]
    
    tagVector2   <- tagVector  [x[[1]]:x[[2]]]
    valueVector2 <- valueVector[x[[1]]:x[[2]]]
    
    ###################################################################
    ## built ms set
    spectrumItem <- list()
    spectrumItem[tagVector2[!isPeak2 & !isEmpty2]] <- trimws(valueVector2[!isPeak2 & !isEmpty2])
    
    if(!is.null(spectrumItem$"Num Peaks"))
      if(spectrumItem$"Num Peaks" == "0" & !returnEmptySpectra)
        return(NULL)
    
    peakLines <- fileLines2[isPeak2]
    peakLines <- trimws(gsub(x = peakLines, pattern = "\".*\"", replacement = ""))
    
    if(multiplePeaksPerLine){
      peakLines <- unlist(strsplit(x = peakLines, split = "[ \t]"))
    } else {
      peakLines <- unlist(lapply(X = strsplit(x = peakLines, split = "[ \t]"), FUN = function(peaktokens){peaktokens[1:2]}))
    }
    spectrumItem["peaks"] <- paste(peakLines, collapse = " ")
    spectrumItem["peaks"] <- trimws(gsub(x = spectrumItem["peaks"], pattern = "  ", replacement = " "))
    
    ## check peaks
    tokens <- strsplit(x = spectrumItem[["peaks"]], split = "[ \t]")[[1]]
    if(length(tokens) == 0){#spectrumItem$"Num Peaks" == "0"){
      mzs  <- character(0)
      ints <- character(0)
      warning(paste(basename(fileSpectra), ": Empty peak list for [", paste(x, collapse = ","), "]", sep = ""))
    } else {
      mzs  <- tokens[seq(from=1, to=length(tokens), by = 2)]
      ints <- tokens[seq(from=2, to=length(tokens), by = 2)]
    }
    if(length(mzs) != length(ints)) stop("error in parsing peaks")
    
    ## handle duplicated tags
    duplicatedTags    <- unique(names(spectrumItem)[duplicated(names(spectrumItem))])
    if(length(duplicatedTags) > 0){
      duplicated <- sapply(X = duplicatedTags, FUN = function(x){
        unlist(sapply(X = seq_along(spectrumItem), FUN = function(y){ if(names(spectrumItem[y])==x) return(y) }))
      }, simplify = F)
      
      indecesToRemove <- vector(mode = "integer", length = 0)
      for(idx in seq_along(duplicated)){
        indeces <- duplicated[[idx]]
        representant <- indeces[[1]]
        indecesToRemoveHere <- indeces[-1]
        
        spectrumItem[[representant]] <- paste(spectrumItem[indeces], sep = "; ")
        indecesToRemove <- c(indecesToRemove, indecesToRemoveHere)
      }
      spectrumItem <- spectrumItem[-indecesToRemove]
    }
    
    return(spectrumItem)
  })
  
  rm(
    isName,
    isNAme,
    isTITLE,
    isBI,
    isID,
    isPeak,
    tagVector,
    valueVector
  )
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "MS/MS file: Box") else print("MS/MS file: Box")
  
  ## remove NULL entries?
  spectraList[unlist(lapply(X = spectraList, FUN = is.null))] <- NULL
  
  numberOfSpectra <- length(spectraList)
  
  ## postprocess
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file postprocessing", sep = "")) else print(paste("MS/MS file postprocessing", sep = ""))
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file boxing", sep = "")) else print(paste("MS/MS file boxing", sep = ""))
  returnObj <- list()
  returnObj$fileSpectra <- fileSpectra
  returnObj$spectraList <- spectraList
  returnObj$numberOfSpectra <- numberOfSpectra
  
  return(returnObj)
}





####################################################################################
## built matrix

builtMatrix <- function(spectraList, 
                        mzDeviationAbsolute_grouping, 
                        mzDeviationInPPM_grouping, 
                        doMs2PeakGroupDeisotoping, 
                        mzDeviationAbsolute_ms2PeakGroupDeisotoping, 
                        mzDeviationInPPM_ms2PeakGroupDeisotoping, 
                        proportionOfMatchingPeaks_ms2PeakGroupDeisotoping, 
                        progress = FALSE)
{
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping preprocessing...", sep = "")) else print(paste("Fragment grouping preprocessing...", sep = ""))
  numberOfSpectra <- length(spectraList)
  
  fragment_mz   <- as.vector(unlist(lapply(X = spectraList, FUN = function(x){x$ms2Peaks_mz })))
  fragment_int  <- as.vector(unlist(lapply(X = spectraList, FUN = function(x){x$ms2Peaks_int})))
  fragment_spec <- as.vector(unlist(lapply(X = spectraList, FUN = function(x){x$peakNumber  })))
  fragment_spec <- rep(x = seq_len(numberOfSpectra), times = fragment_spec)
  numberOfMS2Peaks <- length(fragment_mz)
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping preprocessing ready", sep = "")) else print(paste("Fragment grouping preprocessing ready", sep = ""))
  if(!is.na(progress))  if(progress)  incProgress(amount = 0,     detail = paste("Fragment grouping", sep = "")) else print(paste("Fragment grouping", sep = ""))
  startTime <- Sys.time()
  
  resultObj <- mzClustGeneric(
    p = matrix(data = c(fragment_mz, fragment_spec), nrow = numberOfMS2Peaks, ncol = 2), 
    sampclass = NULL, 
    mzppm = mzDeviationInPPM_grouping, mzabs = mzDeviationAbsolute_grouping, 
    minsamp = 1, minfrac = 0,
    progress
  )
  
  endTime <- Sys.time()
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.2,     detail = paste("Fragment grouping ready (", difftime(time1 = endTime, time2 = startTime, units = "secs"), "s)", sep = "")) else print(paste("Fragment grouping ready (", difftime(time1 = endTime, time2 = startTime, units = "secs"), "s)", sep = ""))
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group postprocessing", sep = "")) else print(paste("Fragment group postprocessing", sep = ""))
  startTime <- Sys.time()
  matrixRows <- vector(mode = "numeric")
  matrixCols <- vector(mode = "numeric")
  matrixVals <- vector(mode = "numeric")
  
  numberOfMS2PeakGroups <- nrow(resultObj$mat)
  fragmentMasses <- resultObj$mat[, "mzmed"]
  for(groupIdx in seq_len(numberOfMS2PeakGroups)){
    if(numberOfMS2PeakGroups > 10)
      if((groupIdx %% (as.integer(numberOfMS2PeakGroups/10))) == 0)
        if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Fragment group postprocessing: ", groupIdx, " / ", numberOfMS2PeakGroups, sep = "")) else print(paste("Fragment group postprocessing: ", groupIdx, " / ", numberOfMS2PeakGroups, sep = ""))
    
    groupMembers <- resultObj$idx[[groupIdx]]
    numberOfFragmentsInGroup <- length(groupMembers)
    
    matrixCols <- c(matrixCols, rep(x = groupIdx, times = numberOfFragmentsInGroup))
    matrixRows <- c(matrixRows, fragment_spec[groupMembers])
    matrixVals <- c(matrixVals, fragment_int[groupMembers])
  }
  
  numberOfCollisions <- sum(duplicated(cbind(matrixRows, matrixCols)))
  endTime <- Sys.time()
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Fragment group postprocessing ready (", difftime(time1 = endTime, time2 = startTime, units = "secs"), "s)", sep = "")) else print(paste("Fragment group postprocessing ready (", difftime(time1 = endTime, time2 = startTime, units = "secs"), "s)", sep = ""))
  
  if(length(matrixRows) == 0){
    ## box results
    if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group boxing", sep = "")) else print(paste("Fragment group boxing", sep = ""))
    returnObj <- list()
    returnObj$matrix <- matrix(nrow = 0, ncol = 0)
    returnObj$numberOfSpectra <- numberOfSpectra
    returnObj$fragmentMasses <- vector(mode = "numeric", length = 0)
    returnObj$numberOfCollisions <- numberOfCollisions
    
    returnObj$numberOfMS2Peaks <- numberOfMS2Peaks
    returnObj$numberOfMS2PeaksPrior <- 0
    returnObj$numberOfRemovedMS2IsotopePeaks <- 0
    returnObj$numberOfMS2PeakGroups <- numberOfMS2PeakGroups
    returnObj$numberOfMS2PeakGroupsPrior <- 0
    returnObj$numberOfRemovedMS2PeakGroupIsotopeColumns <- 0
    
    return(returnObj)
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Boxing to matrix", sep = "")) else print(paste("Boxing to matrix", sep = ""))
  matrix <- sparseMatrix(i = matrixRows, j = matrixCols, x = matrixVals, dims = c(numberOfSpectra, numberOfMS2PeakGroups))
  
  orderTempCol <- order(fragmentMasses)
  matrix <- matrix[, orderTempCol]
  fragmentMasses <- fragmentMasses[orderTempCol]
  
  ## deisotoping
  numberOfMS2PeakGroupsPrior <- numberOfMS2PeakGroups
  numberOfRemovedMS2PeakGroupIsotopeColumns <- 0
  numberOfMS2PeaksPrior <- numberOfMS2Peaks
  numberOfRemovedMS2IsotopePeaks <- 0
  if(doMs2PeakGroupDeisotoping){
    if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group deisotoping", sep = "")) else print(paste("Fragment group deisotoping", sep = ""))
    startTime <- Sys.time()
    
    distance13Cminus12C <- 1.0033548378
    ## mark isotope precursors
    ms2PeakGroupsToRemove <- vector(mode = "logical", length = numberOfMS2PeakGroups)
    for(ms2PeakGroupIdx in seq_len(numberOfMS2PeakGroups)){
      if(numberOfMS2PeakGroups > 10)
        if((ms2PeakGroupIdx %% (as.integer(numberOfMS2PeakGroups/10))) == 0)
          if(!is.na(progress)){
            if(progress)  incProgress(amount = 0.0, detail = paste("Fragment group deisotoping ", ms2PeakGroupIdx, " / ", numberOfMS2PeakGroups, sep = "")) else print(paste("Fragment group deisotoping ", ms2PeakGroupIdx, " / ", numberOfMS2PeakGroups, sep = ""))
          }
      mzError <- abs(fragmentMasses[[ms2PeakGroupIdx]] * mzDeviationInPPM_ms2PeakGroupDeisotoping / 1E6)
      mzError <- max(mzError, mzDeviationAbsolute_ms2PeakGroupDeisotoping)
      
      ## MZ difference around 1.0033548378 (first isotope) or 1.0033548378 * 2 (second isotope)
      if(fragmentMasses[[ms2PeakGroupIdx]] > 0){
        ## fragment
        distances <- (fragmentMasses[[ms2PeakGroupIdx]] - distance13Cminus12C)     - fragmentMasses[-ms2PeakGroupIdx]
      } else {
        ## neutral loss
        distances <- (fragmentMasses[[ms2PeakGroupIdx]] + distance13Cminus12C)     - fragmentMasses[-ms2PeakGroupIdx]
      }
      validInMz <- abs(distances) <= mzError
      
      if(!any(validInMz))
        next
      validInMz <- which(validInMz)
      if(fragmentMasses[[ms2PeakGroupIdx]] < 0)
        validInMz <- validInMz + 1
      
      ## isotopic fragments are mainly in spectra with monoisotopic fragments
      fragmentIntensitiesHere <- matrix[, ms2PeakGroupIdx]
      
      isotopicThere <- fragmentIntensitiesHere != 0
      numberOfFragmentPeaksHere <- sum(isotopicThere)
      
      validInOverlap <- apply(X = matrix(data = matrix[, validInMz], 
                                         nrow = numberOfSpectra), 
                              MARGIN = 2, FUN = function(x){
                                monoisotopicThere <- x != 0
                                precursorInCommon <- isotopicThere & monoisotopicThere
                                validOverlap <- (sum(precursorInCommon) / sum(isotopicThere)) >= proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
                                return(validOverlap)
                              })
      
      if(!any(validInOverlap))
        next
      
      monoisotopicFragmentColumn <- min(validInMz[validInOverlap])
      
      
      ## intensity gets smaller in the isotope spectrum
      monoisotopicThere <- matrix[, monoisotopicFragmentColumn] != 0
      numberOfMonoisotopicFragmentPeaks <- sum(monoisotopicThere)
      precursorInCommon <- isotopicThere & monoisotopicThere
      validToRemove <- (matrix[, monoisotopicFragmentColumn] > matrix[, ms2PeakGroupIdx]) & isotopicThere
      numberOfRemovedPeaks <- sum(validToRemove)
      matrix[validToRemove, ms2PeakGroupIdx] <- rep(x = 0, times = numberOfRemovedPeaks)
      numberOfRemovedMS2IsotopePeaks <- numberOfRemovedMS2IsotopePeaks + numberOfRemovedPeaks
      
      if(sum(matrix[, ms2PeakGroupIdx] != 0) == 0)
        ms2PeakGroupsToRemove[[ms2PeakGroupIdx]] <- TRUE
    }
    
    ## remove
    matrix <- matrix[, !ms2PeakGroupsToRemove]
    
    numberOfRemovedMS2PeakGroupIsotopeColumns <- sum(ms2PeakGroupsToRemove)
    fragmentMasses <- fragmentMasses[!ms2PeakGroupsToRemove]
    numberOfMS2PeakGroups <- length(fragmentMasses)
    numberOfMS2Peaks <- numberOfMS2PeaksPrior - numberOfRemovedMS2IsotopePeaks
    
    endTime <- Sys.time()
    if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group deisotoping ready (", difftime(time1 = endTime, time2 = startTime, units = "secs"), "s)", sep = "")) else print(paste("Fragment group deisotoping ready (", difftime(time1 = endTime, time2 = startTime, units = "secs"), "s)", sep = ""))
  }
  
  precursorMz <- unlist(lapply(X = spectraList, FUN = function(x){
    if(is.null(x))  return(NULL)
    else            return(x$mz)
  }))
  precursorRt <- unlist(lapply(X = spectraList, FUN = function(x){
    if(is.null(x))  return(NULL)
    else            return(x$rt)
  }))
  matrix@Dimnames[[1]] <- paste(precursorMz, precursorRt, sep = " / ") #TODO: mz / rt / x for uniqueness...?
  matrix@Dimnames[[2]] <- fragmentMasses
  
  ## box results
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group boxing", sep = "")) else print(paste("Fragment group boxing", sep = ""))
  returnObj <- list()
  returnObj$matrix <- matrix
  returnObj$numberOfSpectra <- numberOfSpectra
  returnObj$fragmentMasses <- fragmentMasses
  returnObj$numberOfCollisions <- numberOfCollisions
  
  returnObj$numberOfMS2Peaks <- numberOfMS2Peaks
  returnObj$numberOfMS2PeaksPrior <- numberOfMS2PeaksPrior
  returnObj$numberOfRemovedMS2IsotopePeaks <- numberOfRemovedMS2IsotopePeaks
  returnObj$numberOfMS2PeakGroups <- numberOfMS2PeakGroups
  returnObj$numberOfMS2PeakGroupsPrior <- numberOfMS2PeakGroupsPrior
  returnObj$numberOfRemovedMS2PeakGroupIsotopeColumns <- numberOfRemovedMS2PeakGroupIsotopeColumns
  
  return(returnObj)
}

#' Alignment of high resolution mass spectra
#'
#' adapted from R package xcms: xcms_1.44.0, package path: R/mzClust.R
#' Reference: 
#' Alignment of high resolution mass spectra: development of a heuristic approach for metabolomics
#' Metabolomics June 2006, Volume 2, Issue 2, pp 75-83
#' http://link.springer.com/article/10.1007%2Fs11306-006-0021-7
#' 
#' @param p 
#' @param sampclass 
#' @param mzppm 
#' @param mzabs 
#' @param minsamp 
#' @param minfrac 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples
mzClustGeneric <- function(p, 
                           sampclass=NULL, 
                           mzppm = 20, 
                           mzabs = 0, 
                           minsamp = 1, 
                           minfrac=0.5, 
                           progress = FALSE)
{
  makeBin <- function(pos){
    if(pos > numpeaks)
      return(list(pos=pos,bin=c(-1)))
    
    bin <- pord[pos]
    pos <- pos+1
    basepeak <- p[bin[1],1]
    
    error_range <- c(basepeak, abs(basepeak)*error_window+basepeak+2*mzabs)
    while(pos < numpeaks && p[pord[pos],1] <= error_range[2]) {
      bin <- c(bin,pord[pos])
      pos <- pos + 1
    }
    
    lst <- list(pos=pos,bin=bin)
    lst
  }
  meanDeviationOverLimit <- function(bin){
    bin_mz <- p[bin,1]
    m <- mean(bin_mz)
    error_range <- c(m-ppm_error*abs(m)-mzabs, ppm_error*abs(m)+m+mzabs)
    if(length(bin_mz[(bin_mz > error_range[2]) |
                     (bin_mz < error_range[1])]) > 0 ) {
      return(TRUE)
    } else { FALSE }
  }
  bin2output <- function(bin){
    gcount <- integer(length(classnum))
    if(length(gcount) != 0){
      for(i in seq(along = bin)){
        class_idx <- sampclass[p[bin[i],2]]
        gcount[class_idx] <- gcount[class_idx] + 1
      }
    }
    ## make sure, that if no classes given, 'any' is false
    if(length(bin) < minsamp || (!any(gcount >= classnum*minfrac) && length(gcount)>0))
      return(list())
    groupvec <- c(rep(NA,4+length(gcount)))
    groupvec[1] <- mean(p[bin,1])
    groupvec[2:3] <- range(p[bin,1])
    groupvec[4] <- length(bin)
    sorted <- order(p[bin,1])
    grp_members <- bin[sorted]
    groupvec[4+seq(along = gcount)] <- gcount
    lst <- list(stat=groupvec,members=grp_members)
    lst
  }
  ppm_error <- mzppm/1000000
  error_window <- 2*ppm_error
  
  ## numeric version of classlabel
  if(is.null(sampclass)){
    classnum <- integer(0)
    classnames <- seq(along=classnum)
  } else {
    classnames <- levels(sampclass)
    sampclass <- as.vector(unclass(sampclass))
    
    classnum <- integer(max(sampclass))
  }
  
  for(i in seq(along = classnum))
    classnum[i] <- sum(sampclass == i)
  
  numpeaks <- nrow(p)
  
  groupmat <- matrix(nrow = 512, ncol = 4 + length(classnum))
  groupindex <- vector("list", 512)
  
  pord <- order(p[,1])
  pos <- c(1)
  binNumber <- 1
  newbin <- makeBin(pos)
  binA <- newbin$bin
  pos <- newbin$pos
  
  lastOut <- proc.time()["user.self"]
  lastPos <- 1
  
  loopCounter <- 0
  while(TRUE){
    loopCounter <- loopCounter + 1
    
    if(binNumber +4 > nrow(groupmat)){
      groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
      groupindex <- c(groupindex, vector("list", length(groupindex)))
    }
    
    ## progress output
    time <- proc.time()["user.self"]
    if(time - lastOut > 1){
      lastOut <- time
      peakProgress <- (pos - lastPos) / numpeaks
      lastPos <- pos
      if(!is.na(progress))  if(progress)  incProgress(amount = peakProgress * 0.2,     detail = paste("Fragment grouping ", pos, " / ", numpeaks, sep = "")) else {
        print(paste("Fragment grouping ", pos, " / ", numpeaks, sep = ""))
      }
    }
    
    newbin <- makeBin(pos)
    binB <- newbin$bin
    pos <- newbin$pos
    
    if(binB[1] < 0){
      ## cancel
      out <- bin2output(binA)
      if(length(out) != 0){
        groupmat[binNumber,] <- out$stat
        groupindex[[binNumber]] <- out$members
        binNumber <- binNumber + 1
      }
      break
    }
    max_binA <- max(p[binA,1])
    min_binB <- min(p[binB,1])
    
    binclust <- 0
    if(   max_binA + abs(max_binA)*error_window+2*mzabs >= min_binB 
          && min_binB - abs(min_binB)*error_window - 2*mzabs <= max_binA){
      binC <- c(binA,binB)
      binclust <- 1
    } else {
      if(meanDeviationOverLimit(binA)){
        binC <- binA
        binclust <- 2
      }
    }
    
    ## case: not in range or mean deviation over limit
    ## perform hierarchical clustering
    if(binclust != 0){
      
      if(length(unique(p[binC,1])) > 10000){
        ## debugging
        stop(paste("Too many fragments for clustering:", 
                   ppm_error, mzabs, length(p[binC,1]), length(unique(p[binC,1]))))
      }
      
      bin <- p[binC,1]
      uniqueBin <- unique(bin)
      
      ## xcms:::mzClust_hclust is using 
      ## the fast C implementation: .C("R_mzClust_hclust"
      mzFragmentGroups <- xcms:::mzClust_hclust(uniqueBin,ppm_error,mzabs)
      mzFragmentGroups2 <- vector(mode = "integer", length = length(bin))
      for(idx in seq_along(uniqueBin))
        mzFragmentGroups2[bin==uniqueBin[[idx]]] <- mzFragmentGroups[[idx]]
      
      mzFragmentGroups <- mzFragmentGroups2
      
      last_group <- mzFragmentGroups[which.max(p[binC,1])]
      binA <- binC[which(mzFragmentGroups == last_group)]
      
      ## bug fix where there were not enough empty rows in the matrix 
      ## (in case of more than four new mzFragmentGroups)
      if(binNumber + last_group > nrow(groupmat)){
        groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
        groupindex <- c(groupindex, vector("list", length(groupindex)))
      }
      
      if(max(mzFragmentGroups) >1){
        for(c in 1:max(mzFragmentGroups)){
          if(c == last_group){
            next
          }
          tmp_grp <- which(mzFragmentGroups == c)
          tmp_c <- binC[tmp_grp]
          out <- bin2output(tmp_c)
          if(length(out) != 0){
            groupmat[binNumber,] <- out$stat
            groupindex[[binNumber]] <- out$members
            binNumber <- binNumber + 1
          }
        }
      }
    }
    
    if(binclust != 1){
      out <- bin2output(binA)
      if(length(out) != 0){
        groupmat[binNumber,] <- out$stat
        groupindex[[binNumber]] <- out$members
        binNumber <- binNumber + 1
      }
      binA <- binB
    }
  }
  colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "npeaks", classnames)
  
  binNumber <- binNumber - 1
  groupmat <- groupmat[seq(length = binNumber), ]
  groupindex <- groupindex[seq(length = binNumber)]
  cat("\n")
  flush.console()
  return(list(mat=groupmat,idx=groupindex))
}

convertToProjectFile <- function(filePeakMatrixPath, 
                                 fileSpectra,
                                 fileAnnotation,
                                 parameterSet, 
                                 progress = FALSE){
  ####################################################################################
  ## parse MS/MS spectra
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file...", sep = "")) else print(paste("Parsing MS/MS file...", sep = ""))
  
  returnObj <- parseMSP(
    fileSpectra = fileSpectra, 
    minimumIntensityOfMaximalMS2peak = parameterSet$minimumIntensityOfMaximalMS2peak, 
    minimumProportionOfMS2peaks = parameterSet$minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments = parameterSet$neutralLossesPrecursorToFragments,
    neutralLossesFragmentsToFragments = parameterSet$neutralLossesFragmentsToFragments,
    progress = progress
  )
  spectraList <- returnObj$spectraList
  numberOfSpectra <- returnObj$numberOfSpectra
  numberOfSpectraOriginal <- returnObj$numberOfSpectraOriginal
  precursorMz     <- returnObj$precursorMz
  precursorRt     <- returnObj$precursorRt
  numberOfMS2PeaksOriginal <- returnObj$numberOfMS2PeaksOriginal
  numberOfMS2PeaksWithNeutralLosses <- returnObj$numberOfMS2PeaksWithNeutralLosses
  numberOfMS2PeaksAboveThreshold <- returnObj$numberOfMS2PeaksAboveThreshold
  numberOfMS2PeaksBelowThreshold <- returnObj$numberOfMS2PeaksBelowThreshold
  numberOfTooHeavyFragments <- returnObj$numberOfTooHeavyFragments
  numberOfSpectraDiscardedDueToNoPeaks <- returnObj$numberOfSpectraDiscardedDueToNoPeaks
  numberOfSpectraDiscardedDueToMaxIntensity <- returnObj$numberOfSpectraDiscardedDueToMaxIntensity
  numberOfSpectraDiscardedDueToTooHeavy <- returnObj$numberOfSpectraDiscardedDueToTooHeavy
  
  if(numberOfSpectra == 0)
    return("Number of spectra is zero")
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file ready", sep = "")) else print(paste("Parsing MS/MS file ready", sep = ""))
  
  rm(returnObj)
  
  
  filePeakMatrixQF <- readMSDial(filePeakMatrixPath)
  if (!is.null(fileAnnotation)){
    # TODO: determine colums to merge by
    filePeakMatrixQF <- addSiriusAnnotations(filePeakMatrixQF,fileAnnotation)
  }
  
  returnObj <- convertToProjectFile2(
    filePeakMatrixQF = filePeakMatrixQF, 
    spectraList = spectraList, precursorMz = precursorMz, precursorRt = precursorRt, 
    metaboliteFamilies = rep(x = "", times = numberOfSpectra), uniqueMetaboliteFamilies = NULL, metaboliteFamilyColors = NULL, 
    parameterSet = parameterSet, 
    progress = progress
  )
  returnObj$numberOfSpectraOriginal <- numberOfSpectraOriginal
  returnObj$numberOfMS2PeaksOriginal <- numberOfMS2PeaksOriginal
  returnObj$numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses
  returnObj$numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold
  returnObj$numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold
  returnObj$numberOfTooHeavyFragments <- numberOfTooHeavyFragments
  returnObj$numberOfSpectraDiscardedDueToNoPeaks <- numberOfSpectraDiscardedDueToNoPeaks
  returnObj$numberOfSpectraDiscardedDueToMaxIntensity <- numberOfSpectraDiscardedDueToMaxIntensity
  returnObj$numberOfSpectraDiscardedDueToTooHeavy <- numberOfSpectraDiscardedDueToTooHeavy
  
  return(returnObj)
}

convertToProjectFile2 <- function(filePeakMatrixQF, 
                                  spectraList, 
                                  precursorMz, 
                                  precursorRt, 
                                  metaboliteFamilies, 
                                  uniqueMetaboliteFamilies, 
                                  metaboliteFamilyColors, 
                                  furtherProperties = list(), 
                                  parameterSet, 
                                  progress = FALSE)
{
  numberOfSpectraParsed <- length(spectraList)
  
  ####################################################################################
  ## metabolite profile
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = paste("Parsing MS1 file...", sep = "")) else print(paste("Parsing MS1 file...", sep = ""))
  
  if(!is.null(filePeakMatrixQF)){
    returnObj <- parsePeakAbundanceMatrixQF(
      qfeatures = filePeakMatrixQF, 
      doPrecursorDeisotoping = parameterSet$doPrecursorDeisotoping, 
      mzDeviationInPPM_precursorDeisotoping = parameterSet$mzDeviationInPPM_precursorDeisotoping, 
      mzDeviationAbsolute_precursorDeisotoping = parameterSet$mzDeviationAbsolute_precursorDeisotoping, 
      maximumRtDifference = parameterSet$maximumRtDifference,
      progress = progress
    )
    dataFrame <- returnObj$dataFrame
    oldFormat <- returnObj$oldFormat
    numberOfPrecursors <- returnObj$numberOfPrecursors
    dataColumnStartEndIndeces <- returnObj$dataColumnStartEndIndeces
    numberOfDataColumns <- returnObj$numberOfDataColumns
    
    groupLabels          <- returnObj$sampleClass
    sampleType           <- returnObj$sampleType
    sampleInjectionOrder <- returnObj$sampleInjectionOrder
    batchID              <- returnObj$batchID
    
    numberOfParsedMs1Features <- returnObj$numberOfPrecursorsPrior
    numberOfRemovedPrecursorIsotopePeaks <- returnObj$numberOfRemovedIsotopePeaks
    
    rm(returnObj)
  } else {
    propList <- list(
      "Average Rt(min)" = precursorRt,
      "Average Mz" = precursorMz,
      "Metabolite name" = unlist(lapply(X = spectraList, FUN = function(x){x$name})),
      "Adduct ion name" = unlist(lapply(X = spectraList, FUN = function(x){x$adduct}))
    )
    propList <- c(propList, furtherProperties)
    propList <- c(propList, list("MySample" = rep(x = 0, times = numberOfSpectraParsed)))
    dataFrame <- data.frame(propList, check.names = FALSE, stringsAsFactors = FALSE)
    oldFormat <- FALSE
    
    numberOfPrecursors <- numberOfSpectraParsed
    dataColumnStartEndIndeces <- c(5,5) + length(furtherProperties)
    numberOfDataColumns <- 1
    
    groupLabels <- c("Unknown")
    sampleType <- c("Sample")
    sampleInjectionOrder <- -1
    batchID              <- -1
    
    numberOfPrecursorsPrior <- numberOfPrecursors
    numberOfRemovedPrecursorIsotopePeaks <- 0
    
    numberOfParsedMs1Features <- -1
  }
  
  isGC <- "EI spectrum" %in% colnames(dataFrame)
  if(isGC){
    ## in case of GC-MS data set this value according to retention time
    #assignment <- unlist(lapply(X = precursorRt, FUN = function(x){
    #  precursorMz[[which.min(abs(x-dataFrame$"Average Rt(min)"))]]
    #}))
    
    #diffAll <- abs(outer(X = precursorRt, Y = dataFrame$"Average Rt(min)", FUN = function(x, y){abs(x-y)}))
    #allHits <- apply(X = diffAll, MARGIN = 2, FUN = function(x){which(x == min(x[x < parameterSet$mzDeviationAbsolute_mapping], Inf))})
    #assignment <- unlist(lapply(X = allHits, FUN = function(x){
    #  if(length(x) == 0)
    #    return(-1)
    #  return(precursorMz[[x[[1]]]])
    #}))
    
    dataFrame$"Average Mz" <- as.numeric(dataFrame$"Quant mass")
    
    if(nrow(dataFrame) == length(spectraList)){
      ## replace rounded rts in the MS1 data by the unrounded ones from the MS/MS data
      dataFrame$"Average Rt(min)" <- unlist(lapply(X = spectraList, FUN = function(x){x$rt}))
    }
    
    ## remove empty spectra
    nullSpectra <- unlist(lapply(X = spectraList, FUN = is.null))
    if(any(nullSpectra)){
      spectraList[nullSpectra] <- NULL
      dataFrame <- dataFrame[-which(nullSpectra),]
      metaboliteFamilies <- metaboliteFamilies[-which(nullSpectra)] ## unique metabolite families and colors...?
      furtherProperties <- lapply(X = furtherProperties, FUN = function(props){props[-which(nullSpectra)]})
      numberOfPrecursors <- length(spectraList)
      precursorMz <- dataFrame$"Average Mz"
      precursorRt <- dataFrame$"Average Rt(min)"
      numberOfSpectra <- length(spectraList)
    }
  }
  
  ## remove redundant MS1 features
  precursorLabels <- paste(dataFrame$"Average Mz", dataFrame$"Average Rt(min)", sep = " / ")
  dupplicated <- which(duplicated(precursorLabels))
  numberOfDupplicated <- length(dupplicated)
  if(numberOfDupplicated > 0){
    precursorLabels <- precursorLabels[-dupplicated]
    dataFrame <- dataFrame[-dupplicated, ]
  }
  numberOfPrecursors <- numberOfPrecursors - numberOfDupplicated
  
  ####################################################################################
  ## map MS1 and MS/MS to each other
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Postprocessing matrix...") else print("Postprocessing matrix...")
  
  ## order rows by precursor m/z and order columns by fragment group m/z (needed for deisotoping)
  orderMS1features <- order(precursorMz)
  
  precursorMz <- precursorMz[orderMS1features]
  precursorRt <- precursorRt[orderMS1features]
  spectraList <- spectraList[orderMS1features]
  metaboliteFamilies <- metaboliteFamilies[orderMS1features]
  furtherProperties <- lapply(X = furtherProperties, FUN = function(props){props[orderMS1features]})
  
  #TODO: resolve ?
  #temporary fix
  #filePeakMatrix <- NULL
  
  if(!is.null(filePeakMatrixQF)){
    ## allHits: dataFrame$"Average Mz" --> precursorMz; allHits indexes the spectraList
    diffAll <- abs(outer(X = precursorMz, Y = dataFrame$"Average Mz", FUN = function(x, y){abs(x-y)}))
    allHits <- apply(X = diffAll, MARGIN = 2, FUN = function(x){which(x == min(x[x < parameterSet$mzDeviationAbsolute_mapping], Inf))})
    rm(diffAll)
  } else {
    allHits <- lapply(X = dataFrame$"Average Mz", FUN = function(x){which(precursorMz == x)})
  }
  if(is.array(allHits))
    allHits <- as.list(as.data.frame(allHits))
  
  numberOfUnmappedPrecursorsMz <- 0
  numberOfUnmappedPrecursorsRt <- 0
  
  for(i in seq_len(numberOfPrecursors)){
    numberOfItems <- length(allHits[[i]])
    if(numberOfItems == 1)
      if(is.na(allHits[[i]]))
        numberOfItems <- 0
    
    if(numberOfItems == 0){    ## no hit
      allHits[[i]] <- NA
      numberOfUnmappedPrecursorsMz <- numberOfUnmappedPrecursorsMz + 1
    }
    else{    ## take hit with minimum absolute RT difference
      absoluteDifferences <- abs(dataFrame$"Average Rt(min)"[i] - precursorRt[allHits[[i]]])
      bestIdx <- which.min(absoluteDifferences)
      if(length(bestIdx) > 1)
        bestIdx <- bestIdx[[1]]
      if(absoluteDifferences[[bestIdx]] <= parameterSet$maximumRtDifference){
        allHits[[i]] <- allHits[[i]][[bestIdx]]
      } else {
        allHits[[i]] <- NA
        numberOfUnmappedPrecursorsRt <- numberOfUnmappedPrecursorsRt + 1
      }
    }
  }
  ## apply
  hitsTempAll     <- unlist(allHits)
  naRows <- is.na(hitsTempAll)
  hitsTempAll     <- hitsTempAll[!naRows]
  ## sub set matrices
  numberOfUnmappedSpectra <- length(spectraList) - length(spectraList[hitsTempAll])
  spectraList      <- spectraList[hitsTempAll]
  numberOfSpectra <- length(spectraList)
  ## sub set 
  metaboliteFamiliesAll <- metaboliteFamilies[hitsTempAll]
  furtherPropertiesAll <- lapply(X = furtherProperties, FUN = function(props){props[hitsTempAll]})
  
  dataFrame <- dataFrame[!naRows, ]
  precursorMzAll <- dataFrame$"Average Mz"
  precursorRtAll <- dataFrame$"Average Rt(min)"
  numberOfPrecursors <- length(precursorMzAll)
  numberOfUnmappedPrecursors <- sum(naRows)
  
  ####################################################################################
  ## built matrix
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Building fragment mzFragmentGroups...", sep = "")) else print(paste("Building fragment mzFragmentGroups...", sep = ""))
  if(isGC){
    ## round GC mzs to integer
    for(spectrumIdx in seq_len(length(spectraList))){
      spectraList[[spectrumIdx]]$ms2Peaks_mz <- round(x = spectraList[[spectrumIdx]]$ms2Peaks_mz, digits = 1)
    }
  }
  
  returnObj <- builtMatrix(
    spectraList = spectraList, 
    mzDeviationAbsolute_grouping = parameterSet$mzDeviationAbsolute_grouping, 
    mzDeviationInPPM_grouping = parameterSet$mzDeviationInPPM_grouping, 
    doMs2PeakGroupDeisotoping = parameterSet$doMs2PeakGroupDeisotoping, 
    mzDeviationAbsolute_ms2PeakGroupDeisotoping = parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping,
    mzDeviationInPPM_ms2PeakGroupDeisotoping = parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping,
    proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping, 
    progress = progress
  )
  matrix <- returnObj$matrix
  fragmentMasses <- returnObj$fragmentMasses
  numberOfMS2PeakGroups <- returnObj$numberOfMS2PeakGroups
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Building fragment mzFragmentGroups ready", sep = "")) else print(paste("Building fragment mzFragmentGroups ready", sep = ""))
  
  orderTempCol <- order(fragmentMasses)
  matrix <- matrix[, orderTempCol]
  fragmentMasses <- fragmentMasses[orderTempCol]
  
  rm(returnObj)
  
  ########################################
  ## filter empty fragment mzFragmentGroups
  numberOfFragments <- length(fragmentMasses)
  fragmentGroupNonEmpty <- vector(mode = "logical", length = numberOfFragments)
  for(colIdx in seq_len(numberOfFragments))
    fragmentGroupNonEmpty[[colIdx]] <- sum(matrix[, colIdx] != 0) > 0
  matrix <- matrix[, fragmentGroupNonEmpty]
  fragmentMasses <- fragmentMasses[fragmentGroupNonEmpty]
  numberOfFragments <- length(fragmentMasses)
  
  ########################################
  ## filter small fragment masses
  minimumFragmentMass <- 5
  tmp <- abs(fragmentMasses) < minimumFragmentMass
  fragmentMasses <- fragmentMasses[!tmp]
  matrix     <- matrix[, !tmp]
  numberOfMS2PeakGroupsAll <- ncol(matrix)
  
  ########################################
  ## convert sparse matrix to three-vector format
  dgTMatrix <- as(matrix, "dgTMatrix")
  rm(matrix)
  gc()
  
  matrixRows <- dgTMatrix@i
  matrixCols <- dgTMatrix@j
  matrixVals <- dgTMatrix@x
  rm(dgTMatrix)
  gc()
  
  ## minimum index from 0 to 1
  matrixRows <- matrixRows + 1
  matrixCols <- matrixCols + 1
  
  numberOfItems <- length(matrixRows)
  
  ########################################
  ## first fragments, then neutral losses
  colIdx <- min(which(fragmentMasses > 0))
  numberOfColumns <- length(fragmentMasses)
  
  matrixCols[matrixCols < colIdx] <- matrixCols[matrixCols < colIdx] + numberOfColumns
  matrixCols <- matrixCols - (colIdx - 1)
  
  fragmentMasses <- unlist(c(fragmentMasses[colIdx:numberOfColumns], fragmentMasses[seq_len(colIdx - 1)]))
  
  ## row ordering
  numberOfAnnotationColumns <- ncol(dataFrame)
  rowOrder <- order(precursorMzAll)
  
  ####################################################################################
  ## matrix assembly with additional columns and column head
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Boxing...") else print("Boxing...")
  
  ## additional rows stuff
  matrixRows <- as.integer(matrixRows)
  matrixCols <- as.integer(matrixCols)
  
  ## TODO performance 30s
  numberOfFragments <- length(fragmentMasses)
  meanIntensity <- vector(mode = "numeric", length = numberOfFragments)
  fragmentCount <- vector(mode = "numeric", length = numberOfFragments)
  for(colIdx in seq_len(numberOfFragments)){
    intensities <- matrixVals[matrixCols == colIdx]
    fragmentCount[[colIdx]] <- length(intensities)
    meanIntensity[[colIdx]] <- mean(x = intensities)
  }
  
  ## add stuff
  numberOfRows <- length(precursorMzAll)
  numberOfPrimaryAnnotationColumns <- 3
  columnOffset <- numberOfPrimaryAnnotationColumns + numberOfAnnotationColumns# + length(furtherPropertiesAll)
  matrixCols <- matrixCols + columnOffset
  
  ######################################
  ## additional columns
  
  ## precursor m/z
  matrixRows <- c(matrixRows, seq_len(numberOfRows))
  matrixCols <- c(matrixCols, rep(x = 1, times = numberOfRows))
  matrixVals <- c(matrixVals, precursorMzAll)
  
  ## precursor RT
  matrixRows <- c(matrixRows, seq_len(numberOfRows))
  matrixCols <- c(matrixCols, rep(x = 2, times = numberOfRows))
  matrixVals <- c(matrixVals, precursorRtAll)
  
  ## annotation column
  matrixRows <- c(matrixRows, seq_len(numberOfRows))
  matrixCols <- c(matrixCols, rep(x = 3, times = numberOfRows))
  matrixVals <- c(matrixVals, metaboliteFamiliesAll)
  
  ## all present stuff
  for(colIdx in seq_len(numberOfAnnotationColumns)){
    matrixRows <- c(matrixRows, seq_len(numberOfRows))
    matrixCols <- c(matrixCols, rep(x = numberOfPrimaryAnnotationColumns + colIdx, times = numberOfRows))
    matrixVals <- c(matrixVals, dataFrame[, colIdx])
  }
  
  ## sort rows
  rowOrderReverse <- vector(mode = "numeric", length = numberOfRows)
  for(i in seq_len(numberOfRows))
    rowOrderReverse[[i]] <- which(rowOrder == i)
  matrixRows <- rowOrderReverse[matrixRows]
  
  ######################################
  ## additional rows
  numberOfColumns <- numberOfFragments + columnOffset
  
  ##################
  ## additional row: head
  headRow <- c("m/z", "RT", "Annotation", names(dataFrame)[seq_len(numberOfAnnotationColumns)], fragmentMasses)
  matrixRows <- matrixRows + 1
  
  matrixRows <- c(matrixRows, rep(x = 1, times = numberOfColumns))
  matrixCols <- c(matrixCols, seq_len(numberOfColumns))
  matrixVals <- c(matrixVals, headRow)
  
  ##################
  ## additional row: average fragment intensity
  intensityRow <- c(rep(x = "", times = numberOfPrimaryAnnotationColumns + numberOfAnnotationColumns), meanIntensity)
  
  matrixRows <- matrixRows + 1
  
  matrixRows <- c(matrixRows, rep(x = 1, times = numberOfColumns))
  matrixCols <- c(matrixCols, seq_len(numberOfColumns))
  matrixVals <- c(matrixVals, intensityRow)
  
  ## columns tags
  if(numberOfDataColumns > 0){
    dataColumnStartIdx <- numberOfPrimaryAnnotationColumns + dataColumnStartEndIndeces[[1]]
    dataColumns <- dataColumnStartIdx:(dataColumnStartIdx + numberOfDataColumns - 1)
  } else {
    dataColumns <- NULL
  }
  
  ## AnnotationColors={AS=#0000FF, SQT-glucosides=#FF0000}
  if(!is.null(uniqueMetaboliteFamilies) & !is.null(metaboliteFamilyColors)){
    annotationColorsValue <- paste(uniqueMetaboliteFamilies, metaboliteFamilyColors, sep = "=", collapse = ", ")
  } else {
    annotationColorsValue <- ""
  }
  annotationColorsFieldValue <- paste("AnnotationColors={", annotationColorsValue, "}", sep = "")
  
  matrixRows <- c(matrixRows, rep(x = 1, times = 3 + numberOfDataColumns))
  matrixCols <- c(matrixCols, 1, 2, 3, dataColumns)
  matrixVals <- c(matrixVals, "ID", "ID", annotationColorsFieldValue, groupLabels)
  
  ##################
  ## additional row: number of present fragment entries
  fragmentCountRow <- c(rep(x = "", times = numberOfPrimaryAnnotationColumns + numberOfAnnotationColumns), fragmentCount)
  
  matrixRows <- matrixRows + 1
  
  matrixRows <- c(matrixRows, rep(x = 1, times = numberOfColumns))
  matrixCols <- c(matrixCols, seq_len(numberOfColumns))
  matrixVals <- c(matrixVals, fragmentCountRow)
  
  ##################
  ## return
  resultObj <- list(
    matrixRows = as.numeric(unlist(matrixRows)),
    matrixCols = as.numeric(unlist(matrixCols)),
    matrixVals = unlist(matrixVals),
    numberOfPrecursors = numberOfPrecursors,
    numberOfParsedSpectra = numberOfSpectraParsed,
    numberOfParsedMs1Features = numberOfParsedMs1Features,
    numberOfRemovedPrecursorIsotopePeaks = numberOfRemovedPrecursorIsotopePeaks,
    numberOfDuplicatedPrecursors = numberOfDupplicated,
    numberOfUnmappedSpectra = numberOfUnmappedSpectra,
    numberOfUnmappedPrecursors = numberOfUnmappedPrecursors,
    numberOfUnmappedPrecursorsMz = numberOfUnmappedPrecursorsMz,
    numberOfUnmappedPrecursorsRt = numberOfUnmappedPrecursorsRt
  )
  
  if(!is.na(progress))  if(progress)  setProgress(1) else print("Ready")
  return(resultObj)
}
