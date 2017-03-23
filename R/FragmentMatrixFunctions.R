
library("xcms")
library("Matrix")

####################################################################################
## aligned spectra
parsePeakAbundanceMatrix <- function(filePeakMatrix, doPrecursorDeisotoping, mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, maximumRtDifference, progress){
  ## known columns
  assumedColumns <- c(
    ## in file            ## parsed
    "Alignment ID",       "Alignment.ID",
    "Average Rt(min)",    "Average.Rt.min.",
    "Average Mz",         "Average.Mz",
    "Metabolite name",    "Metabolite.name",
    "Adduct ion name",    "Adduct.ion.name",
    "Fill %",             "Fill..",
    "MS/MS included",     "MS.MS.included",
    "INCHIKEY",           "INCHIKEY",
    "SMILES",             "SMILES",
    "LINK",               "LINK",
    "SCORE",              "SCORE",
    "Dot product",        "Dot.product",
    "Reverse dot product","Reverse.dot.product",
    "Fragment presence %","Fragment.presence.%",
    "Spectrum reference file name", "Spectrum.reference.file.name"
  )
  
  ## read file
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = paste("Parsing MS¹ file content...", sep = "")) else print(paste("Parsing MS¹ file content...", sep = ""))
  
  dataFrameAll <- read.table(filePeakMatrix, header=FALSE, sep = "\t", as.is=TRUE, quote = "", check.names = FALSE, comment.char = "")
  dataFrameHeader <- dataFrameAll[1:4, ]
  dataFrame <- dataFrameAll[5:nrow(dataFrameAll), ]
  colnames(dataFrame) <- dataFrameHeader[4, ]
  
  numberOfPrecursors <- nrow(dataFrame)
  numberOfPrecursorsPrior <- numberOfPrecursors
  
  ## columns
  columnIndexEndOfAnnotation <- max(match(x = "Class", table = dataFrameHeader[1, ]), na.rm = TRUE)
  
  if(ncol(dataFrame) > columnIndexEndOfAnnotation){
    dataColumnStartEndIndeces <- c(columnIndexEndOfAnnotation + 1, ncol(dataFrame))
    numberOfDataColumns <- dataColumnStartEndIndeces[[2]] - dataColumnStartEndIndeces[[1]] + 1
    dataColumnNames <- colnames(dataFrame)[dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
    
    sampleClass <- dataFrameHeader[1, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
    sampleType  <- dataFrameHeader[2, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
    sampleInjectionOrder <- dataFrameHeader[3, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
  } else {
    dataColumnStartEndIndeces <- NULL
    numberOfDataColumns <- 0
    dataColumnNames <- NULL
    
    sampleClass <- NULL
    sampleType  <- NULL
    sampleInjectionOrder <- NULL
  }
  
  ## sorted by m/z (needed for deisotoping)
  if(!is.null(dataFrame$"Average Mz"))
    dataFrame <- dataFrame[order(dataFrame$"Average Mz"), ]
  
  ## replace -1 by 0
  if(numberOfDataColumns > 0){
    for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
      dataFrame[ , colIdx] <- as.numeric(dataFrame[ , colIdx])
      if(!is.na(sum(dataFrame[,colIdx] == -1)))
        dataFrame[(dataFrame[,colIdx] == -1),colIdx] <- 0
    }
  }
  
  ########################################
  ## column formats
  if(!is.null(dataFrame$"Average Rt(min)"))     dataFrame$"Average Rt(min)"     <- as.numeric(dataFrame$"Average Rt(min)")
  if(!is.null(dataFrame$"Average.Rt.min."))     dataFrame$"Average.Rt.min."     <- as.numeric(dataFrame$"Average.Rt.min.")
  if(!is.null(dataFrame$"Average Mz"))          dataFrame$"Average Mz"          <- as.numeric(dataFrame$"Average Mz")
  if(!is.null(dataFrame$"Average.Mz"))          dataFrame$"Average.Mz"          <- as.numeric(dataFrame$"Average.Mz")
  if(!is.null(dataFrame$"Fill %"))              dataFrame$"Fill %"              <- as.numeric(dataFrame$"Fill %")
  if(!is.null(dataFrame$"Fill.."))              dataFrame$"Fill.."              <- as.numeric(dataFrame$"Fill..")
  if(!is.null(dataFrame$"MS/MS included"))      dataFrame$"MS/MS included"      <- as.logical(dataFrame$"MS/MS included")
  if(!is.null(dataFrame$"MS.MS.included"))      dataFrame$"MS.MS.included"      <- as.logical(dataFrame$"MS.MS.included")
  if(!is.null(dataFrame$"Dot product"))         dataFrame$"Dot product"         <- as.numeric(dataFrame$"Dot product")
  if(!is.null(dataFrame$"Dot.product"))         dataFrame$"Dot.product"         <- as.numeric(dataFrame$"Dot.product")
  if(!is.null(dataFrame$"Reverse dot product")) dataFrame$"Reverse dot product" <- as.numeric(dataFrame$"Reverse dot product")
  if(!is.null(dataFrame$"Reverse.dot.product")) dataFrame$"Reverse.dot.product" <- as.numeric(dataFrame$"Reverse.dot.product")
  if(!is.null(dataFrame$"Fragment presence %")) dataFrame$"Fragment presence %" <- as.numeric(dataFrame$"Fragment presence %")
  if(!is.null(dataFrame$"Fragment.presence.%")) dataFrame$"Fragment.presence.%" <- as.numeric(dataFrame$"Fragment.presence.%")
  
  ########################################
  ## deisotoping
  numberOfRemovedIsotopePeaks <- 0
  if(doPrecursorDeisotoping & !is.null(dataFrame$"Average Mz")){
    if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Precursor deisotoping...", sep = "")) else print(paste("Precursor deisotoping...", sep = ""))
    
    distance13Cminus12C <- 1.0033548378
    ## mark isotope precursors
    precursorsToRemove <- vector(mode = "logical", length = numberOfPrecursors)
    
    if(numberOfDataColumns > 0){
      intensities <- dataFrame[ , dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
      medians <- apply(X = as.matrix(intensities), MARGIN = 1, FUN = median)
    }
    for(precursorIdx in seq_len(numberOfPrecursors)){
      if((precursorIdx %% (as.integer(numberOfPrecursors/10))) == 0)
        if(!is.na(progress))  if(progress)  incProgress(amount = 0.0, detail = paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = "")) else print(paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = ""))
      
      mzError <- dataFrame$"Average Mz"[[precursorIdx]] * mzDeviationInPPM_precursorDeisotoping / 1000000
      mzError <- max(mzError, mzDeviationAbsolute_precursorDeisotoping)
      
      ## RT difference <= maximumRtDifference
      validPrecursorsInRt <- abs(dataFrame$"Average Rt(min)"[[precursorIdx]] - dataFrame$"Average Rt(min)"[-precursorIdx]) <= maximumRtDifference
      ## MZ difference around 1.0033548378 (first isotope) or 1.0033548378 * 2 (second isotope)
      validPrecursorsInMz1 <- abs((dataFrame$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 1) - dataFrame$"Average Mz"[-precursorIdx]) <= mzError
      validPrecursorsInMz2 <- abs((dataFrame$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 2) - dataFrame$"Average Mz"[-precursorIdx]) <= mzError
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
    dataFrame <- dataFrame[!precursorsToRemove, ]
    
    numberOfRemovedIsotopePeaks <- sum(precursorsToRemove)
    numberOfPrecursors <- nrow(dataFrame)
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Boxing...", sep = "")) else print(paste("Boxing...", sep = ""))
  returnObj <- list()
  returnObj$dataFrame <- dataFrame
  ## meta
  returnObj$numberOfPrecursors <- numberOfPrecursors
  returnObj$dataColumnStartEndIndeces <- dataColumnStartEndIndeces
  returnObj$numberOfDataColumns <- numberOfDataColumns
  ## group anno
  returnObj$sampleClass <- sampleClass
  returnObj$sampleType <- sampleType
  returnObj$sampleInjectionOrder <- sampleInjectionOrder
  ## misc
  returnObj$numberOfPrecursorsPrior <- numberOfPrecursorsPrior
  returnObj$numberOfRemovedIsotopePeaks <- numberOfRemovedIsotopePeaks
  
  return (returnObj)
}

####################################################################################
## parse MS/MS spectra
## 
## NAME: Unknown
## RETENTIONTIME: 3.215358
## PRECURSORMZ: 78.91963
## METABOLITENAME: 
## ADDUCTIONNAME: [M-H]-
## Num Peaks: 2
## 76.97093  754
## 76.98951  754
## 
parseMSP <- function(fileSpectra, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, neutralLossesPrecursorToFragments, neutralLossesFragmentsToFragments, progress = FALSE){
  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Read file") else print("MS/MS file: Read file")
  
  fileLines <- readLines(con = fileSpectra)
  numberOfFileLines <- length(fileLines)
  
  ## GC-MS entry:
  #SCANNUMBER: 518
  #MODELION: 59
  #MODELIONHEIGHT: 924
  #MODELIONAREA: 924
  #INTEGRATEDHEIGHT: 924
  #INTEGRATEDAREA: 924
  
  numberOfMS2Peaks <- 0
  numberOfMS2PeaksWithNeutralLosses <- 0
  numberOfMS2PeaksAboveThreshold <- 0
  numberOfMS2PeaksBelowThreshold <- 0
  
  ## start with empty lines or not?
  endOfRecord <- TRUE
  if(numberOfFileLines > 0)
    if(nchar(trimws(fileLines[[1]])) > 0)
      endOfRecord <- FALSE
  
  ## check for pattern
  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Parse") else print("MS/MS file: Parse")
  isName	<- grepl(pattern = "^Name:",						    	            x = fileLines)
  isNAME	<- grepl(pattern = "^NAME:",						    	            x = fileLines)
  isRT	  <- grepl(pattern = "^RETENTIONTIME:",						          x = fileLines)
  isRt	  <- grepl(pattern = "^retention time:",			  	          x = fileLines)
  isMZ	  <- grepl(pattern = "^PRECURSORMZ:",							          x = fileLines)
  isMz	  <- grepl(pattern = "^precursor m/z:",						          x = fileLines)
  isTotEm	<- grepl(pattern = "^total exact mass:",					        x = fileLines)
  isEMass	<- grepl(pattern = "^exact mass:",							          x = fileLines)
  isMetN	<- grepl(pattern = "^METABOLITENAME:",						        x = fileLines)
  isADDN	<- grepl(pattern = "^ADDUCTIONNAME:",						          x = fileLines)
  isAddN	<- grepl(pattern = "^Adductionname:",						          x = fileLines)
  isScanN	<- grepl(pattern = "^SCANNUMBER:",							          x = fileLines)
  isMIon	<- grepl(pattern = "^MODELION:",							            x = fileLines)
  isPrety	<- grepl(pattern = "^precursor type:",						        x = fileLines)
  isNumP	<- grepl(pattern = "^Num Peaks:",							            x = fileLines)
  isPeak	<- grepl(pattern = "^\\d+(\\.\\d+)?[ \t]\\d+(\\.\\d+)?",	x = fileLines)
  isCoCl	<- grepl(pattern = "^compound class:",						        x = fileLines)
  isInty	<- grepl(pattern = "^instrument type:",						        x = fileLines)
  isInchi	<- grepl(pattern = "^InChI:",								              x = fileLines)
  
  ## extract
  suppressWarnings({
    parsedNAME	 <- 						      trimws(substring(text = fileLines, first = nchar("NAME:") + 1))
    parsedName	 <- 						      trimws(substring(text = fileLines, first = nchar("Name:") + 1))
    parsedRT  	 <- as.numeric(				trimws(substring(text = fileLines, first = nchar("RETENTIONTIME:") + 1)))
    parsedRt     <- as.numeric(unlist(lapply(X = strsplit(x =   trimws(substring(text = fileLines, first = nchar("retention time:") + 1)), split = " "), FUN = function(x){if(length(x)==0) return(NA) else return(x[[1]])})))
    parsedMZ	   <- as.numeric(				trimws(substring(text = fileLines, first = nchar("PRECURSORMZ:") + 1)))
    parsedMz	   <- as.numeric(				trimws(substring(text = fileLines, first = nchar("precursor m/z:") + 1)))
    parsedTotEm  <- as.numeric(				trimws(substring(text = fileLines, first = nchar("total exact mass:") + 1)))
    parsedEMass  <- as.numeric(				trimws(substring(text = fileLines, first = nchar("exact mass:") + 1)))
    parsedMetN	 <- 						      trimws(substring(text = fileLines, first = nchar("METABOLITENAME:") + 1))
    parsedADDN	 <- 						      trimws(substring(text = fileLines, first = nchar("ADDUCTIONNAME:") + 1))
    parsedAddN	 <- 						      trimws(substring(text = fileLines, first = nchar("Adductionname:") + 1))
    parsedScanN  <- as.numeric(				trimws(substring(text = fileLines, first = nchar("SCANNUMBER:") + 1)))
    parsedMIon	 <- as.numeric(				trimws(substring(text = fileLines, first = nchar("MODELION:") + 1)))
    parsedPrety  <- 						      trimws(substring(text = fileLines, first = nchar("precursor type:") + 1))
    parsedNumP	 <- as.numeric(				trimws(substring(text = fileLines, first = nchar("Num Peaks:") + 1)))
    
    parsedTokensTmp <- strsplit(x = trimws(fileLines), split = "[ \t]")
    #parsedms2Peaks_mz <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){head(x = x, n = 1)})))
    parsedms2Peaks_mz  <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){if(length(x)<1) return(NA) else return(x[[1]])})))
    parsedms2Peaks_int <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){if(length(x)<2) return(NA) else return(x[[2]])})))
    
    parsedCoCl	 <- 						      trimws(substring(text = fileLines, first = nchar("compound class:") + 1))
    parsedInty	 <- 						      trimws(substring(text = fileLines, first = nchar("instrument type:") + 1))
    parsedInchi  <- 						      trimws(substring(text = fileLines, first = nchar("InChI:") + 1))
  })
  
  ## entry line intervals in file
  entryBorders   <- c(which(isName | isNAME), length(fileLines))
  entryIntervals <- matrix(data = unlist(lapply(X = seq_len(length(entryBorders) - 1), FUN = function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
  
  ## do it
  if(progress)  incProgress(amount = 0, detail = "MS/MS file: Assemble spectra") else print("MS/MS file: Assemble spectra")
  spectraList <- apply(X = entryIntervals, MARGIN = 2, FUN = function(x){
    parsedNAME2	 		<- parsedNAME	 		[x[[1]]:x[[2]]]
    parsedName2	 		<- parsedName	 		[x[[1]]:x[[2]]]
    parsedRT2  	 		<- parsedRT  	 		[x[[1]]:x[[2]]]
    parsedRt2     		<- parsedRt     		[x[[1]]:x[[2]]]
    parsedMZ2	 		<- parsedMZ	 		    [x[[1]]:x[[2]]]
    parsedMz2	 		<- parsedMz	 		    [x[[1]]:x[[2]]]
    parsedTotEm2  		<- parsedTotEm  		[x[[1]]:x[[2]]]
    parsedEMass2  		<- parsedEMass  		[x[[1]]:x[[2]]]
    parsedMetN2	 		<- parsedMetN	 		[x[[1]]:x[[2]]]
    parsedADDN2	 		<- parsedADDN	 		[x[[1]]:x[[2]]]
    parsedAddN2	 		<- parsedAddN	 		[x[[1]]:x[[2]]]
    parsedScanN2  		<- parsedScanN  		[x[[1]]:x[[2]]]
    parsedMIon2	 		<- parsedMIon	 		[x[[1]]:x[[2]]]
    parsedPrety2  		<- parsedPrety  		[x[[1]]:x[[2]]]
    parsedNumP2	 		<- parsedNumP	 		[x[[1]]:x[[2]]]
    parsedms2Peaks_mz2	<- parsedms2Peaks_mz	[x[[1]]:x[[2]]]
    parsedms2Peaks_int2	<- parsedms2Peaks_int	[x[[1]]:x[[2]]]
    parsedCoCl2	 		<- parsedCoCl	 		[x[[1]]:x[[2]]]
    parsedInty2	 		<- parsedInty	 		[x[[1]]:x[[2]]]
    parsedInchi2  		<- parsedInchi  		[x[[1]]:x[[2]]]
    
    isNAME2	 		<- isNAME	 		[x[[1]]:x[[2]]]
    isName2	 		<- isName	 		[x[[1]]:x[[2]]]
    isRT2  	 		<- isRT  	 		[x[[1]]:x[[2]]]
    isRt2     		<- isRt     		[x[[1]]:x[[2]]]
    isMZ2	 		<- isMZ	 		    [x[[1]]:x[[2]]]
    isMz2	 		<- isMz	 		    [x[[1]]:x[[2]]]
    isTotEm2  		<- isTotEm  		[x[[1]]:x[[2]]]
    isEMass2  		<- isEMass  		[x[[1]]:x[[2]]]
    isMetN2	 		<- isMetN	 		[x[[1]]:x[[2]]]
    isADDN2	 		<- isADDN	 		[x[[1]]:x[[2]]]
    isAddN2	 		<- isAddN	 		[x[[1]]:x[[2]]]
    isScanN2  		<- isScanN  		[x[[1]]:x[[2]]]
    isMIon2	 		<- isMIon	 		[x[[1]]:x[[2]]]
    isPrety2  		<- isPrety  		[x[[1]]:x[[2]]]
    isNumP2	 		<- isNumP	 		[x[[1]]:x[[2]]]
    isPeak2	<- isPeak	[x[[1]]:x[[2]]]
    isCoCl2	 		<- isCoCl	 		[x[[1]]:x[[2]]]
    isInty2	 		<- isInty	 		[x[[1]]:x[[2]]]
    isInchi2  		<- isInchi  		[x[[1]]:x[[2]]]
    
    name <- NULL
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
    
    if(any(isName2))
      ## name
      name        <- parsedName2	[[which(isName2)[[1]]]]
    if(any(isNAME2) & !is.null(name))
      ## name
      name        <- parsedNAME2	[[which(isNAME2)[[1]]]]
    if(any(isRT2))
      ## retention time
      rt          <- parsedRT2	[[which(isRT2)[[1]]]]
    if(any(isRt2) & any(is.null(rt), is.na(rt)))
      rt          <- parsedRt2	[[which(isRt2)[[1]]]]
    if(any(isMZ2))
      ## precursor m/z
      mz          <- parsedMZ2	[[which(isMZ2)[[1]]]]
    if(any(isMz2)    & any(is.null(mz), is.na(mz)))
      mz          <- parsedMz2	[[which(isMz2)[[1]]]]
    if(any(isTotEm2) & any(is.null(mz), is.na(mz)) & any(isPrety2)){
      mzTmp <- parsedTotEm2[[which(isTotEm2)[[1]]]]
      pit   <- parsedPrety2[[which(isPrety2)[[1]]]]
      switch(pit,
             "[M-H]-" = { mz <- mzTmp -  1.008  },
             "[M-H]"  = { mz <- mzTmp -  1.008  },
             "[M+H]+" = { mz <- mzTmp +  1.008  },
             "[M+H]"  = { mz <- mzTmp +  1.008  },
             "[M+Na]+"= { mz <- mzTmp + 22.9898 },
             "[M+Na]" = { mz <- mzTmp + 22.9898 }#,
             #stop(paste("Unknown precursor ion type, pit, ")!", sep = ""))
      )
    }
    if(any(isEMass2) & any(is.null(mz), is.na(mz)) & any(isPrety2)){
      mzTmp <- parsedEMass2[[which(isEMass2)[[1]]]]
      pit   <- parsedPrety2[[which(isPrety2)[[1]]]]
      switch(pit,
             "[M-H]-" = { mz <- mzTmp -  1.008  },
             "[M-H]"  = { mz <- mzTmp -  1.008  },
             "[M+H]+" = { mz <- mzTmp +  1.008  },
             "[M+H]"  = { mz <- mzTmp +  1.008  },
             "[M+Na]+"= { mz <- mzTmp + 22.9898 },
             "[M+Na]" = { mz <- mzTmp + 22.9898 }#,
             #stop(paste("Unknown precursor ion type, pit, ")!", sep = ""))
      )
    }
    if(any(isMetN2))
      metName     <- parsedMetN2	[[which(isMetN2)[[1]]]]
    if(any(isADDN2))
      ## adduct
      adduct      <- parsedADDN2	[[which(isADDN2)[[1]]]]
    if(any(isAddN2)  & adduct == "Unknown")
      adduct      <- parsedAddN2	[[which(isAddN2)[[1]]]]
    if(any(isPrety2) & adduct == "Unknown")
      adduct      <- parsedPrety2[[which(isPrety2)[[1]]]]
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
    if(any(isInchi2))
      ## compound class
      inchi  <- parsedInchi2 [[which(isInchi2)[[1]]]]
    
    
    
    ## end of entry
    
    if(is.null(rt))
      rt <- 0
    
    #if(is.null(mz))
    #  ## in case of gc
    #  mz <- max(ms2Peaks_mz)
    
    ###################################################################
    ## filter fragments with mass greater than precursor
    if(!is.null(mz)){
      tooHeavy <- ms2Peaks_mz > mz
      ms2Peaks_mz  <- ms2Peaks_mz [!tooHeavy]
      ms2Peaks_int <- ms2Peaks_int[!tooHeavy]
    }
    
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
        numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold + length(ms2Peaks_mz)
      } else {
        ## spectrum is not considered
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
      if(neutralLossesPrecursorToFragments){
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
    ## built ms set
    spectrumItem <- list(
      name = name,
      rt = round(as.numeric(rt), digits = 2),
      mz = ifelse(test = !is.null(mz), yes = round(as.numeric(mz), digits = 4), no = ifelse(test = !is.null(scanNumber), yes = scanNumber, no = max(ms2Peaks_mz))),
      metName = metName,
      adduct = adduct,
      quantMass = quantMass,
      compoundClass = compoundClass,
      instrumentType = instrumentType,
      inchi = inchi,
      #peakNumber = as.numeric(peakNumber),
      peakNumber = length(ms2Peaks_mz),
      ms2Peaks_mz  = ms2Peaks_mz,
      ms2Peaks_int = ms2Peaks_int
    )
    if(spectrumItem$peakNumber > 0){
      ## add
      numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses + spectrumItem$peakNumber
      return(spectrumItem)
    } else
      return(NULL)
  })
  
  if(progress)  incProgress(amount = 0.1, detail = "MS/MS file: Box") else print("MS/MS file: Box")
  
  ## remove NULL entries?
  spectraList[unlist(lapply(X = spectraList, FUN = is.null))] <- NULL
  
  numberOfSpectra <- length(spectraList)
  
  ## postprocess
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file postprocessing", sep = "")) else print(paste("MS/MS file postprocessing", sep = ""))
  #precursorMz <- vector(mode = "numeric", length = numberOfSpectra)
  #for(spectrumIdx in seq_len(length.out = numberOfSpectra))
  #  precursorMz[[spectrumIdx]] <- spectraList[[spectrumIdx]]$mz
  
  #precursorRt <- vector(mode = "numeric", length = numberOfSpectra)
  #for(spectrumIdx in seq_len(length.out = numberOfSpectra))
  #  precursorRt[[spectrumIdx]] <- spectraList[[spectrumIdx]]$rt
  
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
  returnObj$fileSpectra <- fileSpectra
  returnObj$spectraList <- spectraList
  returnObj$numberOfSpectra <- numberOfSpectra
  returnObj$numberOfMS2Peaks <- numberOfMS2Peaks
  returnObj$numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses
  returnObj$numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold
  returnObj$numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold
  returnObj$precursorMz <- precursorMz
  returnObj$precursorRt <- precursorRt
  
  return(returnObj)
}

####################################################################################
## built matrix

builtMatrix <- function(spectraList, mzDeviationAbsolute_grouping, mzDeviationInPPM_grouping, doMs2PeakGroupDeisotoping, mzDeviationAbsolute_ms2PeakGroupDeisotoping, mzDeviationInPPM_ms2PeakGroupDeisotoping, proportionOfMatchingPeaks_ms2PeakGroupDeisotoping, progress = FALSE){
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping preprocessing...", sep = "")) else print(paste("Fragment grouping preprocessing...", sep = ""))
  numberOfSpectra <- length(spectraList)
  
  # spectrumItem <- list(
  #   name = name,
  #   rt = round(as.numeric(rt), digits = 2),
  #   mz = round(as.numeric(mz), digits = 4),
  #   metName = metName,
  #   adduct = adduct,
  #   #peakNumber = as.numeric(peakNumber),
  #   peakNumber = length(ms2Peaks_mz),
  #   ms2Peaks_mz  = ms2Peaks_mz,
  #   ms2Peaks_int = ms2Peaks_int
  # )
  
  fragment_mz   <- as.vector(unlist(lapply(X = spectraList, FUN = function(x){x$ms2Peaks_mz })))
  fragment_int  <- as.vector(unlist(lapply(X = spectraList, FUN = function(x){x$ms2Peaks_int})))
  fragment_spec <- as.vector(unlist(lapply(X = spectraList, FUN = function(x){x$peakNumber  })))
  fragment_spec <- rep(x = seq_len(numberOfSpectra), times = fragment_spec)
  numberOfMS2Peaks <- length(fragment_mz)
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping preprocessing ready", sep = "")) else print(paste("Fragment grouping preprocessing ready", sep = ""))
  if(!is.na(progress))  if(progress)  incProgress(amount = 0,     detail = paste("Fragment grouping", sep = "")) else print(paste("Fragment grouping", sep = ""))
  
  #resultObj <- xcms:::mzClustGeneric(
  resultObj <- mzClustGeneric(
    p = matrix(data = c(fragment_mz, fragment_spec), nrow = numberOfMS2Peaks, ncol = 2), 
    sampclass = NULL, 
    mzppm = mzDeviationInPPM_grouping, mzabs = mzDeviationAbsolute_grouping, 
    minsamp = 1, minfrac = 0,
    progress
  )
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.2,     detail = paste("Fragment grouping ready", sep = "")) else print(paste("Fragment grouping ready", sep = ""))
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group postprocessing", sep = "")) else print(paste("Fragment group postprocessing", sep = ""))
  matrixRows <- vector(mode = "numeric")
  matrixCols <- vector(mode = "numeric")
  matrixVals <- vector(mode = "numeric")
  
  numberOfMS2PeakGroups <- nrow(resultObj$mat)
  fragmentMasses <- resultObj$mat[, "mzmed"]
  for(groupIdx in seq_len(numberOfMS2PeakGroups)){
    if((groupIdx %% (as.integer(numberOfMS2PeakGroups/10))) == 0)
      if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Fragment group postprocessing: ", groupIdx, " / ", numberOfMS2PeakGroups, sep = "")) else print(paste("Fragment group postprocessing: ", groupIdx, " / ", numberOfMS2PeakGroups, sep = ""))
    
    groupMembers <- resultObj$idx[[groupIdx]]
    numberOfFragmentsInGroup <- length(groupMembers)
    
    matrixCols <- c(matrixCols, rep(x = groupIdx, times = numberOfFragmentsInGroup))
    matrixRows <- c(matrixRows, fragment_spec[groupMembers])
    matrixVals <- c(matrixVals, fragment_int[groupMembers])
  }
  
  #rm(ms2PeakGroupList)
  numberOfCollisions <- sum(duplicated(cbind(matrixRows, matrixCols)))
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Fragment group postprocessing ready", sep = "")) else print(paste("Fragment group postprocessing ready", sep = ""))
  
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
  matrix <- sparseMatrix(i = matrixRows, j = matrixCols, x = matrixVals)
  #rm(matrixRows)
  #rm(matrixCols)
  #rm(matrixVals)
  #gc()
  
  
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
    distance13Cminus12C <- 1.0033548378
    ## mark isotope precursors
    ms2PeakGroupsToRemove <- vector(mode = "logical", length = numberOfMS2PeakGroups)
    for(ms2PeakGroupIdx in seq_len(numberOfMS2PeakGroups)){
      if((ms2PeakGroupIdx %% (as.integer(numberOfMS2PeakGroups/10))) == 0)
        if(!is.na(progress))  if(progress)  incProgress(amount = 0.0, detail = paste("Fragment group deisotoping ", ms2PeakGroupIdx, " / ", numberOfMS2PeakGroups, sep = "")) else print(paste("Fragment group deisotoping ", ms2PeakGroupIdx, " / ", numberOfMS2PeakGroups, sep = ""))
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
      #validInMz1 <- abs((fragmentMasses[[ms2PeakGroupIdx]] - distance13Cminus12C)     - fragmentMasses[-ms2PeakGroupIdx]) <= mzError
      #validInMz2 <- abs((fragmentMasses[[ms2PeakGroupIdx]] - distance13Cminus12C * 2) - fragmentMasses[-ms2PeakGroupIdx]) <= mzError
      #validInMz <- validInMz1 | validInMz2
      if(!any(validInMz))
        next
      validInMz <- which(validInMz)
      if(fragmentMasses[[ms2PeakGroupIdx]] < 0)
        validInMz <- validInMz + 1
      
      ## isotopic fragments are mainly in spectra with monoisotopic fragments
      fragmentIntensitiesHere <- matrix[, ms2PeakGroupIdx]
      isotopicThere <- fragmentIntensitiesHere != 0
      numberOfFragmentPeaksHere <- sum(isotopicThere)
      
      validInOverlap <- apply(X = matrix(data = matrix[, validInMz], nrow = numberOfSpectra), MARGIN = 2, FUN = function(x){
        monoisotopicThere <- x != 0
        precursorInCommon <- isotopicThere & monoisotopicThere
        validOverlap <- (sum(precursorInCommon) / sum(isotopicThere)) >= proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
        return(validOverlap)
      })
      
      if(!any(validInOverlap))
        next
      
      monoisotopicFragmentColumn <- min(validInMz[validInOverlap])
      
      ## intensity gets smaller in the isotope spectrum
      #monoisotopicFragmentIntensities <- matrix[, monoisotopicFragmentColumn]
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
    
    if(!is.na(progress))  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group deisotoping ready", sep = "")) else print(paste("Fragment group deisotoping ready", sep = ""))
  }
  
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

## adapted from R package xcms: xcms_1.44.0, package path: R/mzClust.R
## 
## Reference: 
## Alignment of high resolution mass spectra: development of a heuristic approach for metabolomics
## Metabolomics June 2006, Volume 2, Issue 2, pp 75-83
## http://link.springer.com/article/10.1007%2Fs11306-006-0021-7
mzClustGeneric <- function(p, sampclass=NULL, mzppm = 20, mzabs = 0, minsamp = 1, minfrac=0.5, progress = FALSE){
  makeBin <- function(pos){
    if(pos > numpeaks)
      return(list(pos=pos,bin=c(-1)))
    
    bin <- pord[pos]
    pos <- pos+1
    basepeak <- p[bin[1],1]
    error_range <- c(basepeak, basepeak*error_window+basepeak+2*mzabs)
    while(pos < numpeaks && p[pord[pos],1] <= error_range[2]) {
      bin <- c(bin,pord[pos])
      pos <- pos + 1
    }
    
    #if(pos %% (numpeaks%/%100+1) == 0) {
    #  cat(format(((pos-1)/numpeaks*100),digits=1,nsmall=2)," ")
    #  flush.console()
    #}
    
    lst <- list(pos=pos,bin=bin)
    lst
  }
  meanDeviationOverLimit <- function(bin){
    bin_mz <- p[bin,1]
    m <- mean(bin_mz)
    error_range <- c(m-ppm_error*m-mzabs, ppm_error*m+m+mzabs)
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
  while(TRUE){
    if(binNumber +4 > nrow(groupmat)){
      groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
      groupindex <- c(groupindex, vector("list", length(groupindex)))
    }
    
    time <- proc.time()["user.self"]
    if(time - lastOut > 1){
      lastOut <- time
      peakProgress <- (pos - lastPos) / numpeaks
      lastPos <- pos
      if(!is.na(progress))  if(progress)  incProgress(amount = peakProgress * 0.2,     detail = paste("Fragment grouping ", pos, " / ", numpeaks, sep = "")) else print(paste("Fragment grouping ", pos, " / ", numpeaks, sep = ""))
    }
    
    newbin <- makeBin(pos)
    binB <- newbin$bin
    pos <- newbin$pos
    
    if(binB[1] < 0){
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
    if(max_binA + max_binA*error_window+2*mzabs >= min_binB && min_binB - min_binB*error_window - 2*mzabs <= max_binA){
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
      groups <- xcms:::mzClust_hclust(p[binC,1],ppm_error,mzabs)
      
      last_group <- groups[which.max(p[binC,1])]
      binA <- binC[which(groups == last_group)]
      if(max(groups) >1){
        for(c in 1:max(groups)){
          if(c == last_group){
            next
          }
          tmp_grp <- which(groups == c)
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

convertToProjectFile <- function(filePeakMatrix, fileSpectra, parameterSet, progress = FALSE){
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
  precursorMz     <- returnObj$precursorMz
  precursorRt     <- returnObj$precursorRt
  
  if(numberOfSpectra == 0)
    return("Number of spectra is zero")
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file ready", sep = "")) else print(paste("Parsing MS/MS file ready", sep = ""))
  
  ## out
  #print(paste("parsing file", fileSpectra, "finished"))
  #print(paste("Number of MS2 spectra: ", numberOfSpectra))
  #print(paste("Avg number of MS2 peaks per MS2 spectrum: ", returnObj$numberOfMS2PeaksAboveThreshold / numberOfSpectra))
  #print(paste("Number of MS2 peaks: ", returnObj$numberOfMS2Peaks))
  #print(paste("Number of MS2 peaks with neutral losses: ", returnObj$numberOfMS2PeaksWithNeutralLosses))
  #print(paste("Number of MS2 peaks considered: ", returnObj$numberOfMS2PeaksAboveThreshold))
  #print(paste("Number of MS2 peaks not considered: ", returnObj$numberOfMS2PeaksBelowThreshold))
  
  rm(returnObj)
  return(convertToProjectFile2(
    filePeakMatrix, 
    spectraList, precursorMz, precursorRt, 
    rep(x = "", times = numberOfSpectra), NULL, NULL, 
    parameterSet, 
    progress
  ))
}
## metaboliteFamilies <- rep(x = "", times = numberOfSpectra)
## uniqueMetaboliteFamilies <- NULL
## metaboliteFamilyColors <- NULL
convertToProjectFile2 <- function(filePeakMatrix, spectraList, precursorMz, precursorRt, metaboliteFamilies, uniqueMetaboliteFamilies, metaboliteFamilyColors, parameterSet, progress = FALSE){
  numberOfSpectra <- length(spectraList)
  
  ####################################################################################
  ## metabolite profile
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = paste("Parsing MS¹ file...", sep = "")) else print(paste("Parsing MS¹ file...", sep = ""))
  
  if(!is.null(filePeakMatrix)){
    returnObj <- parsePeakAbundanceMatrix(
      filePeakMatrix = filePeakMatrix, 
      doPrecursorDeisotoping = parameterSet$doPrecursorDeisotoping, 
      mzDeviationInPPM_precursorDeisotoping = parameterSet$mzDeviationInPPM_precursorDeisotoping, 
      mzDeviationAbsolute_precursorDeisotoping = parameterSet$mzDeviationAbsolute_precursorDeisotoping, 
      maximumRtDifference = parameterSet$maximumRtDifference,
      progress = progress
    )
    dataFrame <- returnObj$dataFrame
    numberOfPrecursors <- returnObj$numberOfPrecursors
    dataColumnStartEndIndeces <- returnObj$dataColumnStartEndIndeces
    numberOfDataColumns <- returnObj$numberOfDataColumns
    
    #groupLabels <- returnObj$groupLabels
    groupLabels <- returnObj$sampleClass
    sampleType <- returnObj$sampleType
    sampleInjectionOrder <- returnObj$sampleInjectionOrder
    
    numberOfPrecursorsPrior <- returnObj$numberOfPrecursorsPrior
    numberOfRemovedPrecursorIsotopePeaks <- returnObj$numberOfRemovedIsotopePeaks
    
    rm(returnObj)
  } else {
    dataFrame <- data.frame(
      "Average Rt(min)" = precursorRt,
      "Average Mz" = precursorMz,
      "Metabolite name" = unlist(lapply(X = spectraList, FUN = function(x){x$name})),
      "Adduct ion name" = unlist(lapply(X = spectraList, FUN = function(x){x$adduct})),
      "MySample" = rep(x = 0, times = numberOfSpectra),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    numberOfPrecursors <- numberOfSpectra
    dataColumnStartEndIndeces <- c(5,5)
    numberOfDataColumns <- 1
    
    #groupLabels <- returnObj$groupLabels
    groupLabels <- c("Unknown")
    sampleType <- c("Sample")
    sampleInjectionOrder <- c(-1)
    
    numberOfPrecursorsPrior <- numberOfPrecursors
    numberOfRemovedPrecursorIsotopePeaks <- 0
  }
  
  if(is.null(dataFrame$"Average Mz")){
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
    
    dataFrame$"Average Mz" <- seq_len(nrow(dataFrame))
  }
  
  ## out
  #print(paste("Parsing file ", filePeakMatrix, " finished", sep = ""))
  #print(paste("Number of precursors:", numberOfPrecursors))
  #if(parameterSet$doPrecursorDeisotoping){
  #  print(paste("Number of isotope peaks removed:", numberOfRemovedPrecursorIsotopePeaks, "/", numberOfPrecursorsPrior, "=", round(x = numberOfRemovedPrecursorIsotopePeaks / numberOfPrecursorsPrior * 100, digits = 1), "%"))
  #}
  
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
  
  ########################################
  ## map spectra to significant precursors by m/z and RT
  #assignment <- unlist(lapply(X = seq_along(precursorMz), FUN = function(i){
  #  mz <- precursorMz[[i]]
  #  indeces <- which(abs(mz-dataFrame$"Average Mz") <= parameterSet$mzDeviationAbsolute_mapping)
  #  if(length(indeces) == 0)
  #    return(NA)
  #  if(length(indeces) == 1)
  #    return(indeces)
  #  rt <- precursorRt[[i]]
  #  absoluteDifferences <- abs(dataFrame$"Average Rt(min)"[indeces] - rt)
  #  bestIdx <- which.min(absoluteDifferences)
  #  return(indeces[[bestIdx]])
  #}))
  isGC <- "EI spectrum" %in% colnames(dataFrame)
  
  if(!isGC){
    ## order rows by precursor m/z and order columns by fragment group m/z (needed for deisotoping)
    orderMS1features <- order(precursorMz)
    
    precursorMz <- precursorMz[orderMS1features]
    precursorRt <- precursorRt[orderMS1features]
    spectraList <- spectraList[orderMS1features]
    metaboliteFamilies <- metaboliteFamilies[orderMS1features]
    
    diffAll <- abs(outer(X = precursorMz, Y = dataFrame$"Average Mz", FUN = function(x, y){abs(x-y)}))
    allHits <- apply(X = diffAll, MARGIN = 2, FUN = function(x){which(x == min(x[x < parameterSet$mzDeviationAbsolute_mapping], Inf))})
    if(is.array(allHits))
      allHits <- as.list(as.data.frame(allHits))
    
    #bestHits <- as.integer(rep(x = NA, times = numberOfPrecursors))
    rm(diffAll)
    for(i in seq_len(numberOfPrecursors)){
      numberOfItems <- length(allHits[[i]])
      if(numberOfItems == 1)
        if(is.na(allHits[[i]]))
          numberOfItems <- 0
      
      if(numberOfItems == 0){    ## no hit
        allHits[[i]] <- NA
      }
      else{    ## take hit with minimum absolute RT difference
        # which(lapply(X = allHits, FUN = length)>1)
        absoluteDifferences <- abs(dataFrame$"Average Rt(min)"[i] - precursorRt[allHits[[i]]])
        bestIdx <- which.min(absoluteDifferences)
        if(length(bestIdx) > 1)
          bestIdx <- bestIdx[[1]]
        if(absoluteDifferences[[bestIdx]] <= parameterSet$maximumRtDifference){
          allHits[[i]] <- allHits[[i]][[bestIdx]]
        } else {
          allHits[[i]] <- NA
        }
        #allHits[[i]] <- allHits[[i]][[1]]
      }
    }
    ## apply
    hitsTempAll     <- unlist(allHits)
    naRows <- is.na(hitsTempAll)
    hitsTempAll     <- hitsTempAll[!naRows]
    ## sub set matrices
    spectraList      <- spectraList[hitsTempAll]
    numberOfSpectra <- length(spectraList)
    ## sub set 
    precursorMzAll <- precursorMz[hitsTempAll]
    precursorRtAll <- precursorRt[hitsTempAll]
    metaboliteFamiliesAll <- metaboliteFamilies[hitsTempAll]
    
    dataFrame <- dataFrame[!naRows, ]
    numberOfPrecursors <- length(precursorMzAll)
  } else {
    ## sub set matrices
    if(FALSE){
      spectraList      <- list()
      for(rowIdx in seq_len(numberOfPrecursors)){
        tokens <- strsplit(x = strsplit(x = dataFrame$"EI spectrum"[rowIdx], split = ";")[[1]], split = " ")
        ms2Peaks_mz  <- as.numeric(unlist(lapply(X = tokens, FUN = function(x){x[[1]]})))
        ms2Peaks_int <- as.numeric(unlist(lapply(X = tokens, FUN = function(x){x[[2]]})))
        ms2Peaks_int <- ms2Peaks_int/max(ms2Peaks_int)
        
        spectraList[[length(spectraList) + 1]] <- list(
          name = dataFrame$`Metabolite name`[[rowIdx]],
          rt = dataFrame$`Average Rt(min)`[[rowIdx]],
          mz = rowIdx,
          metName = dataFrame$`Metabolite name`[[rowIdx]],
          adduct = "Unknown",
          quantMass = dataFrame$`Quant mass`[[rowIdx]],
          compoundClass = "Unknown",
          instrumentType = "Unknown",
          inchi = "Unknown",
          #peakNumber = as.numeric(peakNumber),
          peakNumber = length(ms2Peaks_mz),
          ms2Peaks_mz  = ms2Peaks_mz,
          ms2Peaks_int = ms2Peaks_int
        )
      }
    }
    
    nullSpectra <- unlist(lapply(X = spectraList, FUN = is.null))
    if(any(nullSpectra)){
      spectraList[nullSpectra] <- NULL
      dataFrame <- dataFrame[-which(nullSpectra),]
      metaboliteFamiliesAll <- metaboliteFamilies[-which(nullSpectra)] ## unique metabolite families and colors...?
    } else {
      metaboliteFamiliesAll <- metaboliteFamilies
    }
    
    numberOfPrecursors <- length(spectraList)
    numberOfSpectra <- length(spectraList)
    precursorMzAll <- seq_len(numberOfPrecursors)
    precursorRtAll <- dataFrame$`Average Rt(min)`
    
    #dataFrame <- dataFrame[!naRows, ]
    #numberOfPrecursors <- length(precursorMzAll)
  }
  
  ####################################################################################
  ## built matrix
  #print("")
  #print("building matrix...")
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Building fragment groups...", sep = "")) else print(paste("Building fragment groups...", sep = ""))
  if(isGC){
    for(spectrumIdx in seq_len(length(spectraList))){
      spectraList[[spectrumIdx]]$ms2Peaks_mz <- round(x = spectraList[[spectrumIdx]]$ms2Peaks_mz, digits = 0)
    }
  }
  
  returnObj <- builtMatrix(
    #returnObj <- builtMatrixOld(
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
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Building fragment groups ready", sep = "")) else print(paste("Building fragment groups ready", sep = ""))
  
  orderTempCol <- order(fragmentMasses)
  matrix <- matrix[, orderTempCol]
  fragmentMasses <- fragmentMasses[orderTempCol]
  
  ## out
  #print("building matrix finished")
  #print(paste("Number of processed items: ", returnObj$numberOfMS2Peaks))
  #print(paste("Number of collisions within MS2 spectra:", returnObj$numberOfCollisions, "/", returnObj$numberOfMS2Peaks, "=", round(x = returnObj$numberOfCollisions / returnObj$numberOfMS2Peaks * 100, digits = 1), "%"))
  #if(parameterSet$doMs2PeakGroupDeisotoping){
  #  print(paste("Number of removed MS2 peaks:", returnObj$numberOfRemovedMS2IsotopePeaks, "/", returnObj$numberOfMS2PeaksPrior, "=", round(x = returnObj$numberOfRemovedMS2IsotopePeaks / returnObj$numberOfMS2PeaksPrior * 100, digits = 1), "%"))
  #  print(paste("Number of completely removed MS2 peak groups:", returnObj$numberOfRemovedMS2PeakGroupIsotopeColumns, "/", returnObj$numberOfMS2PeakGroupsPrior, "=", round(x = returnObj$numberOfRemovedMS2PeakGroupIsotopeColumns / returnObj$numberOfMS2PeakGroupsPrior * 100, digits = 1), "%"))
  #}
  
  rm(returnObj)
  
  ########################################
  ## filter empty fragment groups
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
  
  #matrix <- cbind(matrix[, colIdx:length(fragmentMzAll)], matrix[, 1:(colIdx - 1)])
  matrixCols[matrixCols < colIdx] <- matrixCols[matrixCols < colIdx] + numberOfColumns
  matrixCols <- matrixCols - (colIdx - 1)
  
  fragmentMasses <- unlist(c(fragmentMasses[colIdx:numberOfColumns], fragmentMasses[seq_len(colIdx - 1)]))
  
  ## row ordering
  numberOfAnnotationColumns <- ncol(dataFrame)
  rowOrder <- order(precursorMzAll)
  
  #print("postprocessing matrix finished")
  
  ####################################################################################
  ## matrix assembly with additional columns and column head
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = "Boxing...") else print("Boxing...")
  #print("")
  #print("Boxing...")
  
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
  #numberOfPrimaryAnnotationColumns <- 4
  numberOfPrimaryAnnotationColumns <- 3
  columnOffset <- numberOfPrimaryAnnotationColumns + numberOfAnnotationColumns
  matrixCols <- matrixCols + columnOffset
  
  ######################################
  ## additional columns
  
  ## precursor m/z
  #matrix <- cbind(precursorMzAll, matrix)
  matrixRows <- c(matrixRows, seq_len(numberOfRows))
  matrixCols <- c(matrixCols, rep(x = 1, times = numberOfRows))
  matrixVals <- c(matrixVals, precursorMzAll)
  
  ## precursor RT
  #matrix <- cbind(precursorRtAll, matrix)
  matrixRows <- c(matrixRows, seq_len(numberOfRows))
  matrixCols <- c(matrixCols, rep(x = 2, times = numberOfRows))
  matrixVals <- c(matrixVals, precursorRtAll)
  
  ## annotation column
  #matrix <- cbind(precursorRtAll, matrix)
  matrixRows <- c(matrixRows, seq_len(numberOfRows))
  matrixCols <- c(matrixCols, rep(x = 3, times = numberOfRows))
  #matrixVals <- c(matrixVals, rep(x = "", times = numberOfPrecursors))
  matrixVals <- c(matrixVals, metaboliteFamiliesAll)
  
  ## all present stuff
  #matrix <- cbind(dataFrameSignificants[, 1:numberOfAnnotationColumns], matrix)
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
  
  # newMatrixRows <- vector(mode = "numeric", length = length(matrixRows))
  # for(i in 1:length(matrixRows))
  #   newMatrixRows[[i]] <- which(rowOrder == matrixRows[[i]])
  # matrixRows <- newMatrixRows
  
  ######################################
  ## additional rows
  numberOfColumns <- numberOfFragments + columnOffset
  
  ##################
  ## additional row: head
  headRow <- c("m/z", "RT", "Annotation", names(dataFrame)[seq_len(numberOfAnnotationColumns)], fragmentMasses)
  #matrix <- rbind(headRow, matrix)
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
  if(!is.null(uniqueMetaboliteFamilies) & !is.null(metaboliteFamilyColors))
    annotationColorsValue <- paste(uniqueMetaboliteFamilies, metaboliteFamilyColors, sep = "=", collapse = ", ")
  else
    annotationColorsValue <- ""
  annotationColorsFieldValue <- paste("AnnotationColors={", annotationColorsValue, "}", sep = "")
  
  matrixRows <- c(matrixRows, rep(x = 1, times = 3 + numberOfDataColumns))
  matrixCols <- c(matrixCols, 1, 2, 3, dataColumns)
  #matrixVals <- c(matrixVals, "ID", "ID", "AnnotationColors={}", groupLabels)
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
    numberOfPrecursors = numberOfPrecursors
  )
  
  if(!is.na(progress))  if(progress)  setProgress(1) else print("Ready")
  return(resultObj)
}
