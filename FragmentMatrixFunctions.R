
####################################################################################
## aligned spectra
parsePeakAbundanceMatrix <- function(filePeakMatrix, doPrecursorDeisotoping, mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, maximumRtDifference){
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
  dataFrameAll <- read.table(filePeakMatrix, header=FALSE, sep = "\t", as.is=TRUE, quote = "")
  dataFrameHeader <- dataFrameAll[1:4, ]
  dataFrame <- dataFrameAll[5:nrow(dataFrameAll), ]
  colnames(dataFrame) <- dataFrameHeader[4, ]
  
  numberOfPrecursors <- nrow(dataFrame)
  numberOfPrecursorsPrior <- numberOfPrecursors
  
  ## columns
  #columnIndexEndOfAnnotation <- max(match(x = assumedColumns, table = colnames(dataFrame)), na.rm = TRUE)
  #dataColumnStartEndIndeces <- c(columnIndexEndOfAnnotation + 1, ncol(dataFrame))
  #numberOfDataColumns <- dataColumnStartEndIndeces[[2]] - dataColumnStartEndIndeces[[1]] + 1
  #dataColumnNames <- colnames(dataFrame)[dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
  #columnIndexEndOfAnnotation <- max(match(x = assumedColumns, table = colnames(dataFrame)), na.rm = TRUE)
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
  
  ## replicate groups
  #distanceMatrix <- stringdistmatrix(a = dataColumnNames, b = dataColumnNames)
  #cluster <- hclust(d = as.dist(distanceMatrix))
  ##plot(x = cluster)
  #mergesWithLeafs <- apply(X = cluster$merge, MARGIN = 1, FUN = function(x) any(x < 0))
  #cutHeight <- max(cluster$height[mergesWithLeafs])
  #groupLabels <- cutree(tree = cluster, h = cutHeight)
  #groupLabels <- sampleClass
  
  ## sorted by m/z (needed for deisotoping)
  dataFrame <- dataFrame[order(dataFrame$"Average Mz"), ]
  
  ## replace -1 by 0
  if(numberOfDataColumns > 0){
    for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
      dataFrame[ , colIdx] <- as.numeric(dataFrame[ , colIdx])
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
  
  # -  "Alignment ID",       "Alignment.ID",
  # +  "Average Rt(min)",    "Average.Rt.min.",
  # +  "Average Mz",         "Average.Mz",
  # -  "Metabolite name",    "Metabolite.name",
  # -  "Adduct ion name",    "Adduct.ion.name",
  # +  "Fill %",             "Fill..",
  # +  "MS/MS included",     "MS.MS.included",
  # -  "INCHIKEY",           "INCHIKEY",
  # -  "SMILES",             "SMILES",
  # -  "LINK",               "LINK",
  # -  "SCORE",              "SCORE",
  # +  "Dot product",        "Dot.product",
  # +  "Reverse dot product","Reverse.dot.product",
  # +  "Fragment presence %","Fragment.presence.%",
  # -  "Spectrum reference file name", "Spectrum.reference.file.name"
  
  ########################################
  ## deisotoping
  numberOfRemovedIsotopePeaks <- 0
  if(doPrecursorDeisotoping){
    distance13Cminus12C <- 1.0033548378
    ## mark isotope precursors
    precursorsToRemove <- vector(mode = "logical", length = numberOfPrecursors)
    
    if(numberOfDataColumns > 0){
      intensities <- dataFrame[ , dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
      medians <- apply(X = as.matrix(intensities), MARGIN = 1, FUN = median)
    }
    for(precursorIdx in 1:numberOfPrecursors){
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
  fileLines <- readLines(con = fileSpectra)
  numberOfFileLines <- length(fileLines)
  spectraList <- list()
  ms2Peaks <- list()
  numberOfMS2Peaks <- 0
  numberOfMS2PeaksWithNeutralLosses <- 0
  numberOfMS2PeaksAboveThreshold <- 0
  numberOfMS2PeaksBelowThreshold <- 0
  for(lineIdx in 1:numberOfFileLines){
    ## progress
    if((lineIdx %% (as.integer(numberOfFileLines/10))) == 0)
      if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file: ", lineIdx, " / ", numberOfFileLines, sep = "")) else print(paste("MS/MS file: ", lineIdx, " / ", numberOfFileLines, sep = ""))
    
    ## current line
    line <- fileLines[[lineIdx]]
    
    if(nchar(trimws(line)) == 0){
      ## end of entry
      
      ###################################################################
      ## filter for ms2 peak intensity
      peakNumber <- length(ms2Peaks)
      if(peakNumber > 0){
        ms2PeaksFiltered <- list()
        
        if(TRUE){
          ###################################################################
          ## filter for ms2 peak intensity relative to maximum peak intensity in spectrum
          maximumIntensity <- 0
          for(ms2PeakIdx in 1:length(ms2Peaks)){
            if(ms2Peaks[[ms2PeakIdx]]$intensity > maximumIntensity){
              maximumIntensity <- ms2Peaks[[ms2PeakIdx]]$intensity
            }
          }
          if(maximumIntensity >= minimumIntensityOfMaximalMS2peak){
            intensityThreshold <- maximumIntensity * minimumProportionOfMS2peaks
            for(ms2PeakIdx in 1:length(ms2Peaks)){
              numberOfMS2Peaks <- numberOfMS2Peaks + 1
              if(ms2Peaks[[ms2PeakIdx]]$intensity >= intensityThreshold){
                ms2PeaksFiltered[[length(ms2PeaksFiltered) + 1]] <- ms2Peaks[[ms2PeakIdx]]
                numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold + 1
              } else
                numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold + 1
            }
          }
        }
        ####
        if(FALSE){
          ###################################################################
          ## filter for ms2 peak intensity regarding absolute threshold
          ms2PeaksFiltered <- list()
          for(ms2PeakIdx in 1:length(ms2Peaks)){
            numberOfMS2Peaks <- numberOfMS2Peaks + 1
            if(ms2Peaks[[ms2PeakIdx]]$intensity >= minimumMS2PeakIntensity){
              ms2PeaksFiltered[[length(ms2PeaksFiltered) + 1]] <- ms2Peaks[[ms2PeakIdx]]
              numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold + 1
            } else
              numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold + 1
          }
        }
        
        ms2Peaks <- ms2PeaksFiltered
      }
      
      ###################################################################
      ## normalize ms2 peaks to maximum = 1
      peakNumber <- length(ms2Peaks)
      if(peakNumber > 0){
        max <- 0
        for(ms2PeakIdx in 1:length(ms2Peaks))
          if(ms2Peaks[[ms2PeakIdx]]$intensity > max)
            max <- ms2Peaks[[ms2PeakIdx]]$intensity
        for(ms2PeakIdx in 1:length(ms2Peaks))
          ms2Peaks[[ms2PeakIdx]]$intensity <- ms2Peaks[[ms2PeakIdx]]$intensity / max
      }
      
      ###################################################################
      ## add neutral losses
      if(peakNumber > 0){
        #################################
        ## neutral losses regarding the precursor
        if(neutralLossesPrecursorToFragments){
          ms2PeaksNL <- list()
          for(ms2PeakIdx in 1:length(ms2Peaks))
            ms2PeaksNL[[length(ms2PeaksNL) + 1]] <- list(
              mz        = ms2Peaks[[ms2PeakIdx]]$mz - as.numeric(mz),
              intensity = ms2Peaks[[ms2PeakIdx]]$intensity
            )
          ms2Peaks <- c(ms2Peaks, ms2PeaksNL)
        }
        #################################
        ## neutral losses amongst fragments
        if(neutralLossesFragmentsToFragments){
          ms2PeaksNL <- list()
          for(ms2PeakIdx1 in 1:(length(ms2Peaks) - 1))
            for(ms2PeakIdx2 in (ms2PeakIdx1 + 1):length(ms2Peaks))
              ms2PeaksNL[[length(ms2PeaksNL) + 1]] <- list(
                mz        = -abs(ms2Peaks[[ms2PeakIdx1]]$mz - ms2Peaks[[ms2PeakIdx2]]$mz),
                intensity = (ms2Peaks[[ms2PeakIdx1]]$intensity + ms2Peaks[[ms2PeakIdx2]]$intensity) / 2
              )
          ms2Peaks <- c(ms2Peaks, ms2PeaksNL)
        }
      }
      
      ###################################################################
      ## built ms set
      msItem <- list(
        name = name,
        rt = round(as.numeric(rt), digits = 2),
        mz = round(as.numeric(mz), digits = 4),
        metName = metName,
        adduct = adduct,
        #peakNumber = as.numeric(peakNumber),
        peakNumber = length(ms2Peaks),
        ms2Peaks = ms2Peaks
      )
      if(msItem$peakNumber > 0){
        ## add
        spectraList[[length(spectraList) + 1]] <- msItem
        numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses + msItem$peakNumber
      }
      
      ###################################################################
      ## init for next entry
      ms2Peaks <- list()
    } else if(grepl(pattern = "^NAME:", x = line)[[1]]){
      #name        <- strsplit(x = line, split = " ")[[1]][[2]]
      name        <- trimws(substring(text = line, first = nchar("NAME:") + 1))
    } else if(grepl(pattern = "^RETENTIONTIME:", x = line)[[1]]){
      rt          <- trimws(substring(text = line, first = nchar("RETENTIONTIME:") + 1))
    } else if(grepl(pattern = "^PRECURSORMZ:", x = line)[[1]]){
      mz          <- trimws(substring(text = line, first = nchar("PRECURSORMZ:") + 1))
    } else if(grepl(pattern = "^METABOLITENAME:", x = line)[[1]]){
      if(identical(x = trimws(line), y = "METABOLITENAME:")){
        metName <- ""
      } else{
        metName     <- trimws(substring(text = line, first = nchar("METABOLITENAME:") + 1))
      }
    } else if(grepl(pattern = "^ADDUCTIONNAME:", x = line)[[1]]){
      adduct      <- trimws(substring(text = line, first = nchar("ADDUCTIONNAME:") + 1))
    } else if(grepl(pattern = "^Num Peaks:", x = line)[[1]]){
      peakNumber  <- trimws(substring(text = line, first = nchar("Num Peaks:") + 1))
    } else{
      ## MS2 peaks: "178.88669\t230"
      tokens <- strsplit(x = line, split = "[ \t]")[[1]]
      ms2item <- list(
        mz        = as.numeric(tokens[[1]]),
        intensity = as.numeric(tokens[[2]])
      )
      ms2Peaks[[length(ms2Peaks) + 1]] <- ms2item
    }
  }
  numberOfSpectra <- length(spectraList)
  
  ## postprocess
  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file postprocessing", sep = "")) else print(paste("MS/MS file postprocessing", sep = ""))
  precursorMz <- vector(mode = "numeric", length = numberOfSpectra)
  for(spectrumIdx in seq_len(length.out = numberOfSpectra))
    precursorMz[[spectrumIdx]] <- spectraList[[spectrumIdx]]$mz
  
  precursorRt <- vector(mode = "numeric", length = numberOfSpectra)
  for(spectrumIdx in seq_len(length.out = numberOfSpectra))
    precursorRt[[spectrumIdx]] <- spectraList[[spectrumIdx]]$rt
  
  if(progress)  incProgress(amount = 0.01, detail = paste("MS/MS file boxing", sep = "")) else print(paste("MS/MS file boxing", sep = ""))
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
  
  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping preprocessing...", sep = "")) else print(paste("Fragment grouping preprocessing...", sep = ""))
  numberOfSpectra <- length(spectraList)
  
  muList <- vector(mode = "numeric", length = 0)
  ms2PeakGroupList <- list()
  newGroupCreated <- FALSE
  
  matrixRows <- vector(mode = "numeric")
  matrixCols <- vector(mode = "numeric")
  matrixVals <- vector(mode = "numeric")
  itemIndex <- 1
  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping preprocessing ready", sep = "")) else print(paste("Fragment grouping preprocessing ready", sep = ""))
  
  #Rprof("/home/htreutle/Code/RprofMetSWATH03.txt")
  for(spectrumIdx in seq_len(length.out = numberOfSpectra)){
    ## progress
    if((spectrumIdx %% (as.integer(numberOfSpectra/10))) == 0)
      if(progress)  incProgress(amount = 1/numberOfSpectra * 0.2, detail = paste("Fragment grouping: ", spectrumIdx, " / ", numberOfSpectra, sep = "")) else print(paste("Fragment grouping: ", spectrumIdx, " / ", numberOfSpectra, sep = ""))
    #if((spectrumIdx %% (as.integer(numberOfSpectra/100))) == 0)
    #  if(progress)  incProgress(amount = 0.005, detail = paste("Fragment grouping: ", spectrumIdx, " / ", numberOfSpectra, sep = "")) else print(paste("Fragment grouping: ", spectrumIdx, " / ", numberOfSpectra, sep = ""))
    
    ## current spectrum
    spectrum <- spectraList[[spectrumIdx]]
    
    ## iterate MS2 peaks
    for(ms2PeakIdx in 1:spectrum$peakNumber){
      ## current MS2 peak
      ms2Peak <- spectrum$ms2Peaks[[ms2PeakIdx]]
      
      ## determine similar MS2 peak group if there
      #groupIdxForMS2Peak <- -1
      error <- max(mzDeviationAbsolute_grouping, ms2Peak$mz * (mzDeviationInPPM_grouping / 1E6))
      deviations <- abs(muList - ms2Peak$mz)
      indeces <- which(deviations <= error)
      if(length(indeces) > 0){
        ## take best matching fragment group
        groupIdxForMS2Peak <- indeces[[which.min(deviations[indeces])]]
        #groupIdxForMS2Peak <- indeces[[1]]
        
        ## update MS2 peak group
        ms2PeakGroupList[[groupIdxForMS2Peak]][[length(ms2PeakGroupList[[groupIdxForMS2Peak]]) + 1]] <- ms2Peak$mz
        muList[[groupIdxForMS2Peak]] <- mean(x = ms2PeakGroupList[[groupIdxForMS2Peak]])
      } else {
        ## create new MS2 peak group
        groupIdxForMS2Peak <- length(ms2PeakGroupList) + 1
        
        ms2PeakGroupList[[groupIdxForMS2Peak]] <- ms2Peak$mz
        
        # ms2PeakGroupList[[groupIdxForMS2Peak]] <- vector(mode = "numeric", length = 1)
        # 
        # ## update MS2 peak group
        # ms2PeakGroupList[[groupIdxForMS2Peak]][[1]] <- ms2Peak$mz
        muList[[groupIdxForMS2Peak]] <- ms2Peak$mz
      }
      
      ## add MS2 peak intensity
      matrixRows[[itemIndex]] <- spectrumIdx
      matrixCols[[itemIndex]] <- groupIdxForMS2Peak
      matrixVals[[itemIndex]] <- ms2Peak$intensity
      
      # ## update MS2 peak group
      # ms2PeakGroupList[[groupIdxForMS2Peak]][[length(ms2PeakGroupList[[groupIdxForMS2Peak]]) + 1]] <- ms2Peak$mz
      # #muList[[groupIdxForMS2Peak]] <- median(x = unlist(ms2PeakGroupList[[groupIdxForMS2Peak]], use.names = FALSE))
      # muList[[groupIdxForMS2Peak]] <- mean(x = unlist(ms2PeakGroupList[[groupIdxForMS2Peak]], use.names = FALSE))
      
      # ## handle unlisted list
      # if(newGroupCreated){
      #   muListUnlisted <- unlist(muList, use.names = FALSE)
      # } else {
      #   muListUnlisted[[groupIdxForMS2Peak]] <- muList[[groupIdxForMS2Peak]]
      # }
      
      itemIndex <- itemIndex + 1
    }
  }
  #Rprof(NULL)
  #summaryRprof(filename = "/home/htreutle/Code/RprofMetSWATH03.txt")
  
  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group postprocessing", sep = "")) else print(paste("Fragment group postprocessing", sep = ""))
  
  #rm(ms2PeakGroupList)
  numberOfCollisions <- sum(duplicated(cbind(matrixRows, matrixCols)))
  numberOfMS2PeakGroups <- length(muList)
  
  if(length(matrixRows) == 0){
    ## box results
    if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group boxing", sep = "")) else print(paste("Fragment group boxing", sep = ""))
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
  fragmentMasses <- unlist(muList)
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
  numberOfMS2PeaksPrior <- itemIndex - 1
  numberOfMS2Peaks <- numberOfMS2PeaksPrior
  numberOfRemovedMS2IsotopePeaks <- 0
  if(doMs2PeakGroupDeisotoping){
    if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group deisotoping", sep = "")) else print(paste("Fragment group deisotoping", sep = ""))
    distance13Cminus12C <- 1.0033548378
    ## mark isotope precursors
    ms2PeakGroupsToRemove <- vector(mode = "logical", length = numberOfMS2PeakGroups)
    for(ms2PeakGroupIdx in 1:numberOfMS2PeakGroups){
      if((ms2PeakGroupIdx %% (as.integer(numberOfMS2PeakGroups/10))) == 0)
        if(progress)  incProgress(amount = 0.0, detail = paste("Fragment group deisotoping ", ms2PeakGroupIdx, " / ", numberOfMS2PeakGroups, sep = "")) else print(paste("Fragment group deisotoping ", ms2PeakGroupIdx, " / ", numberOfMS2PeakGroups, sep = ""))
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
      monoisotopicFragmentIntensities <- matrix[, monoisotopicFragmentColumn]
      monoisotopicThere <- matrix[, monoisotopicFragmentColumn] != 0
      numberOfMonoisotopicFragmentPeaks <- sum(monoisotopicThere)
      precursorInCommon <- isotopicThere & monoisotopicThere
      validToRemove <- (matrix[, monoisotopicFragmentColumn] > matrix[, ms2PeakGroupIdx]) & isotopicThere
      numberOfRemovedPeaks <- sum(validToRemove)
      matrix[validToRemove, ms2PeakGroupIdx] <- rep(x = 0, times = numberOfRemovedPeaks)
      numberOfRemovedMS2IsotopePeaks <- numberOfRemovedMS2IsotopePeaks + numberOfRemovedPeaks
      
      #print(paste(ms2PeakGroupIdx, "(", numberOfFragmentPeaksHere, ")", "vs", monoisotopicFragmentColumn, "(", numberOfMonoisotopicFragmentPeaks, ")", "is", fragmentMasses[[ms2PeakGroupIdx]], "vs", fragmentMasses[[monoisotopicFragmentColumn]]))
      
      if(sum(matrix[, ms2PeakGroupIdx] != 0) == 0)
        ms2PeakGroupsToRemove[[ms2PeakGroupIdx]] <- TRUE
    }
    
    ## remove
    matrix <- matrix[, !ms2PeakGroupsToRemove]
    
    numberOfRemovedMS2PeakGroupIsotopeColumns <- sum(ms2PeakGroupsToRemove)
    fragmentMasses <- fragmentMasses[!ms2PeakGroupsToRemove]
    numberOfMS2PeakGroups <- length(fragmentMasses)
    numberOfMS2Peaks <- numberOfMS2PeaksPrior - numberOfRemovedMS2IsotopePeaks
  }
  
  ## box results
  if(progress)  incProgress(amount = 0.05, detail = paste("Fragment group boxing", sep = "")) else print(paste("Fragment group boxing", sep = ""))
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

convertToProjectFile <- function(filePeakMatrix, fileSpectra, parameterSet, progress = FALSE){
  ####################################################################################
  ## aligned spectra
  if(progress)  incProgress(amount = 0.1, detail = paste("Parsing MS¹ file...", sep = "")) else print(paste("Parsing MS¹ file...", sep = ""))
  
  returnObj <- parsePeakAbundanceMatrix(
    filePeakMatrix = filePeakMatrix, 
    doPrecursorDeisotoping = parameterSet$doPrecursorDeisotoping, 
    mzDeviationInPPM_precursorDeisotoping = parameterSet$mzDeviationInPPM_precursorDeisotoping, 
    mzDeviationAbsolute_precursorDeisotoping = parameterSet$mzDeviationAbsolute_precursorDeisotoping, 
    maximumRtDifference = parameterSet$maximumRtDifference
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
  
  ## out
  #print(paste("Parsing file ", filePeakMatrix, " finished", sep = ""))
  #print(paste("Number of precursors:", numberOfPrecursors))
  #if(parameterSet$doPrecursorDeisotoping){
  #  print(paste("Number of isotope peaks removed:", numberOfRemovedPrecursorIsotopePeaks, "/", numberOfPrecursorsPrior, "=", round(x = numberOfRemovedPrecursorIsotopePeaks / numberOfPrecursorsPrior * 100, digits = 1), "%"))
  #}
  
  rm(returnObj)
  
  ####################################################################################
  ## parse MS/MS spectra
  
  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file...", sep = "")) else print(paste("Parsing MS/MS file...", sep = ""))
  
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
  
  if(length(spectraList) == 0)
    return("Number of spectra is zero")
  
  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file ready", sep = "")) else print(paste("Parsing MS/MS file ready", sep = ""))
  
  ## out
  #print(paste("parsing file", fileSpectra, "finished"))
  #print(paste("Number of MS2 spectra: ", numberOfSpectra))
  #print(paste("Avg number of MS2 peaks per MS2 spectrum: ", returnObj$numberOfMS2PeaksAboveThreshold / numberOfSpectra))
  #print(paste("Number of MS2 peaks: ", returnObj$numberOfMS2Peaks))
  #print(paste("Number of MS2 peaks with neutral losses: ", returnObj$numberOfMS2PeaksWithNeutralLosses))
  #print(paste("Number of MS2 peaks considered: ", returnObj$numberOfMS2PeaksAboveThreshold))
  #print(paste("Number of MS2 peaks not considered: ", returnObj$numberOfMS2PeaksBelowThreshold))
  
  rm(returnObj)
  
  ####################################################################################
  ## built matrix
  #print("")
  #print("building matrix...")
  
  if(progress)  incProgress(amount = 0.01, detail = paste("Building fragment groups...", sep = "")) else print(paste("Building fragment groups...", sep = ""))
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
  if(progress)  incProgress(amount = 0.01, detail = paste("Building fragment groups ready", sep = "")) else print(paste("Building fragment groups ready", sep = ""))
  
  #print("building matrix finished")
  #print(paste("Number of processed items: ", returnObj$numberOfMS2Peaks))
  #print(paste("Number of collisions within MS2 spectra:", returnObj$numberOfCollisions, "/", returnObj$numberOfMS2Peaks, "=", round(x = returnObj$numberOfCollisions / returnObj$numberOfMS2Peaks * 100, digits = 1), "%"))
  #if(parameterSet$doMs2PeakGroupDeisotoping){
  #  print(paste("Number of removed MS2 peaks:", returnObj$numberOfRemovedMS2IsotopePeaks, "/", returnObj$numberOfMS2PeaksPrior, "=", round(x = returnObj$numberOfRemovedMS2IsotopePeaks / returnObj$numberOfMS2PeaksPrior * 100, digits = 1), "%"))
  #  print(paste("Number of completely removed MS2 peak groups:", returnObj$numberOfRemovedMS2PeakGroupIsotopeColumns, "/", returnObj$numberOfMS2PeakGroupsPrior, "=", round(x = returnObj$numberOfRemovedMS2PeakGroupIsotopeColumns / returnObj$numberOfMS2PeakGroupsPrior * 100, digits = 1), "%"))
  #}
  
  rm(returnObj)
  
  ####################################################################################
  ## postprocess matrix
  
  if(progress)  incProgress(amount = 0.1, detail = "Postprocessing matrix...") else print("Postprocessing matrix...")
  
  ## order rows by precursor m/z and order columns by fragment group m/z (needed for deisotoping)
  orderTempRow <- order(precursorMz)
  orderTempCol <- order(fragmentMasses)
  
  matrix <- matrix[orderTempRow, orderTempCol]
  precursorMz <- precursorMz[orderTempRow]
  precursorRt <- precursorRt[orderTempRow]
  fragmentMasses <- fragmentMasses[orderTempCol]
  
  ########################################
  ## map spectra to significant precursors by m/z and RT
  diffAll <- abs(outer(X = precursorMz, Y = dataFrame$"Average Mz", FUN = function(x, y) x-y))
  allHits <- apply(X = diffAll, MARGIN = 2, FUN = function(x) which(x == min(x[x < parameterSet$mzDeviationAbsolute_mapping], Inf)))
  rm(diffAll)
  for(i in seq_len(length.out = numberOfPrecursors)){
    numberOfItems <- length(allHits[[i]])
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
  matrixAll      <- matrix[hitsTempAll, ]
  ## sub set 
  precursorMzAll <- precursorMz[hitsTempAll]
  precursorRtAll <- precursorRt[hitsTempAll]
  
  dataFrame <- dataFrame[!naRows, ]
  numberOfPrecursors <- length(precursorMzAll)
  
  ########################################
  ## filter by the minimum number of fragment abundance
  #numberOfMS2PeaksPerGroupAll <- apply(X = matrixAll, MARGIN = 2, FUN = function(x) length(which(x > 0)))
  numberOfMS2PeaksPerGroupAll <- vector(mode = "numeric", length = numberOfMS2PeakGroups)
  for(ms2Idx in 1:numberOfMS2PeakGroups)
    numberOfMS2PeaksPerGroupAll[[ms2Idx]] <- sum(matrixAll[, ms2Idx] > 0)
  orderFragmentMzAll <- which(numberOfMS2PeaksPerGroupAll >= parameterSet$minimumNumberOfMS2PeaksPerGroup)
  fragmentMzAll <- fragmentMasses[orderFragmentMzAll]
  matrixAll <- matrixAll[, orderFragmentMzAll]
  numberOfMS2PeakGroupsAll <- ncol(matrixAll)
  
  ########################################
  ## filter small fragment masses
  minimumFragmentMass <- 5
  tmp <- abs(fragmentMzAll) < minimumFragmentMass
  fragmentMzAll <- fragmentMzAll[!tmp]
  matrixAll     <- matrixAll[, !tmp]
  numberOfMS2PeakGroupsAll <- ncol(matrixAll)
  
  ########################################
  ## convert sparse matrix to three-vector format
  dgTMatrix <- as(matrixAll, "dgTMatrix")
  rm(matrixAll)
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
  colIdx <- min(which(fragmentMzAll > 0))
  numberOfColumns <- length(fragmentMzAll)
  
  #matrixAll <- cbind(matrixAll[, colIdx:length(fragmentMzAll)], matrixAll[, 1:(colIdx - 1)])
  matrixCols[matrixCols < colIdx] <- matrixCols[matrixCols < colIdx] + numberOfColumns
  matrixCols <- matrixCols - (colIdx - 1)
  
  fragmentMzAll <- unlist(c(fragmentMzAll[colIdx:numberOfColumns], fragmentMzAll[1:(colIdx - 1)]))
  
  ## row ordering
  numberOfAnnotationColumns <- ncol(dataFrame)
  rowOrder <- order(precursorMzAll)
  
  #print("postprocessing matrix finished")
  
  ####################################################################################
  ## matrix assembly with additional columns and column head
  
  if(progress)  incProgress(amount = 0.1, detail = "Boxing...") else print("Boxing...")
  #print("")
  #print("Boxing...")
  
  ## additional rows stuff
  numberOfFragments <- length(fragmentMzAll)
  meanIntensity <- vector(mode = "numeric", length = numberOfFragments)
  fragmentCount <- vector(mode = "numeric", length = numberOfFragments)
  for(colIdx in 1:numberOfFragments){
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
  #matrixAll <- cbind(precursorMzAll, matrixAll)
  matrixRows <- c(matrixRows, 1:numberOfRows)
  matrixCols <- c(matrixCols, rep(x = 1, times = numberOfRows))
  matrixVals <- c(matrixVals, precursorMzAll)
  
  ## precursor RT
  #matrixAll <- cbind(precursorRtAll, matrixAll)
  matrixRows <- c(matrixRows, 1:numberOfRows)
  matrixCols <- c(matrixCols, rep(x = 2, times = numberOfRows))
  matrixVals <- c(matrixVals, precursorRtAll)
  
  ## annotation column
  #matrixAll <- cbind(precursorRtAll, matrixAll)
  matrixRows <- c(matrixRows, 1:numberOfRows)
  matrixCols <- c(matrixCols, rep(x = 3, times = numberOfRows))
  matrixVals <- c(matrixVals, rep(x = "", times = numberOfPrecursors))
  
  ## all present stuff
  #matrixAll <- cbind(dataFrameSignificants[, 1:numberOfAnnotationColumns], matrixAll)
  for(colIdx in 1:numberOfAnnotationColumns){
    matrixRows <- c(matrixRows, 1:numberOfRows)
    matrixCols <- c(matrixCols, rep(x = numberOfPrimaryAnnotationColumns + colIdx, times = numberOfRows))
    matrixVals <- c(matrixVals, dataFrame[, colIdx])
  }
  
  ## sort rows
  newMatrixRows <- vector(mode = "numeric", length = length(matrixRows))
  for(i in 1:length(matrixRows))
    newMatrixRows[[i]] <- which(rowOrder == matrixRows[[i]])
  matrixRows <- newMatrixRows
  
  ######################################
  ## additional rows
  numberOfColumns <- numberOfFragments + columnOffset
  
  ##################
  ## additional row: head
  headRow <- c("m/z", "RT", "Annotation", names(dataFrame)[1:numberOfAnnotationColumns], fragmentMzAll)
  #matrixAll <- rbind(headRow, matrixAll)
  matrixRows <- matrixRows + 1
  
  matrixRows <- c(matrixRows, rep(x = 1, times = numberOfColumns))
  matrixCols <- c(matrixCols, 1:numberOfColumns)
  matrixVals <- c(matrixVals, headRow)
  
  ##################
  ## additional row: average fragment intensity
  intensityRow <- c(rep(x = "", times = numberOfPrimaryAnnotationColumns + numberOfAnnotationColumns), meanIntensity)
  
  matrixRows <- matrixRows + 1
  
  matrixRows <- c(matrixRows, rep(x = 1, times = numberOfColumns))
  matrixCols <- c(matrixCols, 1:numberOfColumns)
  matrixVals <- c(matrixVals, intensityRow)
  
  ## columns tags
  if(numberOfDataColumns > 0){
    dataColumnStartIdx <- numberOfPrimaryAnnotationColumns + dataColumnStartEndIndeces[[1]]
    dataColumns <- dataColumnStartIdx:(dataColumnStartIdx + numberOfDataColumns - 1)
  } else {
    dataColumns <- NULL
  }
  
  matrixRows <- c(matrixRows, rep(x = 1, times = 3 + numberOfDataColumns))
  matrixCols <- c(matrixCols, 1, 2, 3, dataColumns)
  matrixVals <- c(matrixVals, "ID", "ID", "AnnotationColors={}", groupLabels)
  
  ##################
  ## additional row: number of present fragment entries
  fragmentCountRow <- c(rep(x = "", times = numberOfPrimaryAnnotationColumns + numberOfAnnotationColumns), fragmentCount)
  
  matrixRows <- matrixRows + 1
  
  matrixRows <- c(matrixRows, rep(x = 1, times = numberOfColumns))
  matrixCols <- c(matrixCols, 1:numberOfColumns)
  matrixVals <- c(matrixVals, fragmentCountRow)
  
  ##################
  ## return
  resultObj <- list(
    matrixRows = as.numeric(unlist(matrixRows)),
    matrixCols = as.numeric(unlist(matrixCols)),
    matrixVals = unlist(matrixVals)
  )
  
  if(progress)  setProgress(1) else print("Ready")
  return(resultObj)
}
exampleInput <- function(){
  filePeakMatrix <- "/home/htreutle/Downloads/MetSWATH/20158251022_rawmatrix_0_less.txt"
  fileSpectra    <- "/home/htreutle/Downloads/MetSWATH/20158251022_spectra_0_less.msp"
  parameterSet <- list()
  parameterSet$minimumIntensityOfMaximalMS2peak                  <- 2000
  parameterSet$minimumProportionOfMS2peaks                       <- 0.05
  parameterSet$mzDeviationAbsolute_grouping                      <- 0.01
  parameterSet$mzDeviationInPPM_grouping                         <- 10
  parameterSet$doPrecursorDeisotoping                            <- TRUE
  parameterSet$mzDeviationAbsolute_precursorDeisotoping          <- 0.001
  parameterSet$mzDeviationInPPM_precursorDeisotoping             <- 10
  parameterSet$maximumRtDifference                               <- 0.02
  parameterSet$doMs2PeakGroupDeisotoping                         <- TRUE
  parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping       <- 0.01
  parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping          <- 10
  parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- 0.9
  parameterSet$mzDeviationAbsolute_mapping                       <- 0.01
  parameterSet$minimumNumberOfMS2PeaksPerGroup                   <- 1
  parameterSet$neutralLossesPrecursorToFragments                 <- TRUE
  parameterSet$neutralLossesFragmentsToFragments                 <- FALSE
  progress       <- FALSE
}