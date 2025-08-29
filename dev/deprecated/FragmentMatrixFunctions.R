
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