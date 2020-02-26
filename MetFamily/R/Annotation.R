
# List of 19
# $ classifierName         : chr "library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE"
# $ numberOfSpectra        : int 1355
# $ numberOfPositiveSpectra: int 20
# $ numberOfNegativeSpectra: int 1335
# $ algorithm              :List of 5
# $ class                  : chr "Organic compounds; Alkaloids and derivatives"
# $ fragmentMasses         : num [1:12207] -970 -969 -965 -965 -963 ...
# $ classOfClass           : chr "ChemOnt_SubstanceClass"
# $ classifier             : Named num [1:12207] -0.000749 -0.000749 -0.000749 -0.000749 -0.000749 ...
# ..- attr(*, "names")= chr [1:12207] "-970.4801" "-969.496033333333" "-964.5493" "-964.513" ...
# $ frequentFragments      : Named num [1:29] 0.571 0.286 0.214 0.214 0.214 ...
# ..- attr(*, "names")= chr [1:29] "-15.0237479936883" "-16.0318143323429" "122.036402912621" "-0.000562633278350665" ...
# $ characteristicFragments: Named num [1:28] 0.539 0.275 0.203 0.203 0.143 ...
# ..- attr(*, "names")= chr [1:28] "-15.0237479936883" "-16.0318143323429" "122.036402912621" "-57.0267914509474" ...
# $ quantiles              : num [1:1001] 0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 ...
# $ quantilesValuesPositive: num [1:1001] -0.207 -0.207 -0.207 -0.207 -0.207 ...
# $ quantilesValuesNegative: num [1:1001] -0.374 -0.302 -0.283 -0.256 -0.249 ...
# $ positiveScores         : num [1:60] -0.20664 -0.01071 -0.00964 -0.00964 0.04742 ...
# $ negativeScores         : num [1:4010] -0.374 -0.321 -0.309 -0.302 -0.297 ...
# $ AUC                    : num 0.912
# $ TPR_for_FPR_of_5Percent: num 0.667
# $ TNR_for_FNR_of_5Percent: num 0

## filePath <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/old_colSums/2018-04-16_09:14:12_2018-02-13_09:14:10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Classifier.RData"
## propertiesList <- list($library= "2018-02-13_09:14:10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp", maximumNumberOfScores= "10000", removeRareFragments= "FALSE", mergeDuplicatedSpectra= "TRUE", takeSpectraWithMaximumNumberOfPeaks= "FALSE", classOfClass= "ChemOnt|AllClasses", minimumNumberOfPosSpectraPerClass= "10", minimumNumberOfNegSpectraPerClass= "10", numberOfDataSetDecompositions= "10", proportionTraining= "0.7", unitResolution= "FALSE", methodName= "ColSumsPos", paramsString= "smoothIntensities=FALSE", algoName= "method=ColSumsPos; smoothIntensities=FALSE")
## featureMatrix <- dataList$featureMatrix
## parameterSet <- dataList$importParameterSet
doAnnotation <- function(filePath, propertiesList, featureMatrix, parameterSet, classesWhiteList = NULL, progress = FALSE){
  if(TRUE){
    filePath_ <<- filePath
    propertiesList_ <<- propertiesList
    featureMatrix_ <<- featureMatrix
    parameterSet_ <<- parameterSet
    classesWhiteList_ <<- classesWhiteList
  }
  if(FALSE){
    filePath <- filePath_
    propertiesList <- propertiesList_
    featureMatrix <- featureMatrix_
    parameterSet <- parameterSet_
    classesWhiteList = classesWhiteList_
    progress <- FALSE
  }
  
  if(progress)  incProgress(amount = 0, detail = "Init") else print("Init")
  ################################################
  ## given
  fragmentMasses <- as.numeric(colnames(featureMatrix))
  fragmentMasses_test <- fragmentMasses
  numberOfSpectra <- nrow(featureMatrix)
  
  ################################################
  ## classifier
  fdrThreshold <- 0.05
  
  ## load classifier
  if(progress)  incProgress(amount = 0, detail = "Loading classifier") else print("Loading classifier")
  classifiers_class <- NULL
  load(file = filePath) # --> classifiers_class
  #"classifierName"          "numberOfSpectra"         "numberOfPositiveSpectra" "numberOfNegativeSpectra"
  #"algorithm"               "class"                   "fragmentMasses"          "classOfClass"           
  #"classifier"              "frequentFragments"       "characteristicFragments" "quantiles"              
  #"quantilesValuesPositive" "quantilesValuesNegative" "positiveScores"          "negativeScores"         
  #"AUC"                     "TPR_for_FPR_of_5Percent" "TNR_for_FNR_of_5Percent"
  
  ## get classifier function
  method <- propertiesList$algoName
  algorithms <- getAlgorithms()
  algoNames <- unlist(lapply(X = algorithms, FUN = function(x){x$algoName}))
  trainingAlgorithm <- algorithms[[which(algoNames == method)]]
  #trainingAlgorithm <- classifiers_class[[1]]$algorithm
  smoothIntensities <- trainingAlgorithm$params$smoothIntensities
  
  ##################################################
  ## mapping of fragments
  if(progress)  incProgress(amount = 0.1, detail = "Aligning classifier") else print("Aligning classifier")
  fragmentMasses_classifier <- classifiers_class[[1]]$fragmentMasses
  
  #############################################################
  ## mapping fragmentMasses_test to fragmentMasses_classifier
  unitResolution <- as.logical(propertiesList$unitResolution)
  if(unitResolution){
    ## Unit resolution
    fragmentMassesRounded        <- round(fragmentMasses_test)
    duplicatedFragmentMasses_bool <- duplicated(fragmentMassesRounded)
    duplicatedFragmentMasses <- unique(fragmentMassesRounded[duplicatedFragmentMasses_bool])
    #fragmentMassesNew        <- unique(fragmentMassesRounded)
    fragmentMassesNew        <- fragmentMassesRounded[!duplicatedFragmentMasses_bool]
    
    fragmentIndecesToRemove <- vector(mode = "integer", length = 0)
    
    #uniqueFragmentMasses_Indeces  <- which(!duplicatedFragmentMasses_bool)
    #fragmentMassesMapping <- list()
    #fragmentMassesMappingOverhead <- vector(mode = "numeric", length = length(duplicatedFragmentMasses))
    
    for(idx in seq_along(duplicatedFragmentMasses)){
      #if(all(length(duplicatedFragmentMasses)>=10, (idx %% (as.integer(length(duplicatedFragmentMasses)/10))) == 0, progress))
      #  incProgress(amount = 0., detail = paste("Aligning classifier ", idx, " / ", length(duplicatedFragmentMasses), sep = ""))
      
      indeces <- which(fragmentMassesRounded == duplicatedFragmentMasses[[idx]])
      indexRepresentant <- indeces[[1]]
      indecesToRemove   <- indeces[-1]
      
      #fragmentMassesMapping        [[idx]] <- fragmentMasses[indeces]
      #fragmentMassesMappingOverhead[[idx]] <- length(fragmentIndecesToRemove)
      
      featureMatrix[, indexRepresentant] <- Matrix::rowSums(x = featureMatrix[, indeces])
      fragmentIndecesToRemove <- c(fragmentIndecesToRemove, indecesToRemove)
    }
    
    if(length(fragmentIndecesToRemove) > 0){
      featureMatrix       <- featureMatrix[, -fragmentIndecesToRemove]
      fragmentMasses_test <- fragmentMassesNew
      #numberOfMS2PeakGroups <- length(fragmentMasses_test)
      featureMatrix@Dimnames[[2]] <- fragmentMasses_test
      #spectraCount_fragment <- Matrix::colSums(featureMatrix != 0)
    }
    
    mappedFragmentIndeces_source <- which(fragmentMasses_test       %in% fragmentMasses_classifier)
    mappedFragmentIndeces_target <- which(fragmentMasses_classifier %in% fragmentMasses_test      )
    
    fragmentMassesMapped_bool  <- fragmentMassesRounded %in% fragmentMasses_classifier[mappedFragmentIndeces_target]
    fragmentMassesMappedSource <- fragmentMasses       [fragmentMassesMapped_bool]
    fragmentMassesMappedTarget <- fragmentMassesRounded[fragmentMassesMapped_bool]
    
    ## TODO add index vector for mapping of ClassMasses to aligned ClassMasses
    mappingSpectraToClassDf <- data.frame(
      "SpectraMasses" = fragmentMassesMappedSource,
      "ClassMasses"   = fragmentMassesMappedTarget
    )
  } else {
    ## TODO GC-EI vs unitResolution ???
    mzAbs <- parameterSet$mzDeviationAbsolute_grouping
    mzPPM <- parameterSet$mzDeviationInPPM_grouping
    
    mapping_fm_test <- as.integer(rep(x = NA, times = length(fragmentMasses_test)))
    for(fragmentIdx in seq_along(fragmentMasses_test)){
      ## intentionally equals mzAbs in case of neutral losses (with negative values)
      massError <- max(fragmentMasses_test[[fragmentIdx]] * mzPPM / 1E6, mzAbs)
      
      ## mapping and selection of best hit
      hits <- abs(fragmentMasses_test[[fragmentIdx]] - fragmentMasses_classifier) <= massError
      if(sum(hits) == 0)
        next
      mapping_fm_test[[fragmentIdx]] <- which.min(abs(fragmentMasses_test[[fragmentIdx]] - fragmentMasses_classifier))[[1]]
    }
    
    #numberOfMappedFragments <- sum(!is.na(mapping_fm_test))
    mappedFragmentIndeces_target <- mapping_fm_test[!is.na(mapping_fm_test)]
    mappedFragmentIndeces_source <- which(!is.na(mapping_fm_test))
    
    rm(mapping_fm_test)
    
    mappingSpectraToClassDf <- data.frame(
      "SpectraMasses"      = fragmentMasses           [mappedFragmentIndeces_source],
      "ClassMasses"        = fragmentMasses_classifier[mappedFragmentIndeces_target],
      "SpectraMassIndeces" = mappedFragmentIndeces_source,
      "ClassMassIndeces"   = mappedFragmentIndeces_target#,
      #"fragmentMasses"            = fragmentMasses,
      #"fragmentMasses_classifier" = fragmentMasses_classifier
    )
  }
  if(progress)  incProgress(amount = 0, detail = paste("Aligned ", length(mappingSpectraToClassDf$SpectraMasses), "vs", length(fragmentMasses), "/", length(fragmentMasses_classifier), "fragments / NLs")) else print(paste("Aligned ", length(mappingSpectraToClassDf$SpectraMasses), "vs", length(fragmentMasses), "/", length(fragmentMasses_classifier), "fragments / NLs"))
  
  ## map test spectra fragments to classifier spectra fragments
  featureMatrix2 <- sparseMatrix(dims = c(numberOfSpectra, length(fragmentMasses_classifier)), i={}, j={}) * 1 ## ngCMatrix --> dgCMatrix
  featureMatrix2[, mappedFragmentIndeces_target] <- featureMatrix[, mappedFragmentIndeces_source]
  featureMatrix2@Dimnames[[1]] <- rownames(featureMatrix)
  featureMatrix2@Dimnames[[2]] <- fragmentMasses_classifier[mappedFragmentIndeces_target]
  featureMatrix <- featureMatrix2
  rm(featureMatrix2)
  
  #############################################################
  ## classification
  if(progress)  incProgress(amount = 0.1, detail = "Classification") else print("Classification")
  
  classifierClasses <- unlist(lapply(X = classifiers_class, FUN = function(x){x$class}))
  numberOfClasses <- length(classifierClasses)
  classes <- classifierClasses
  matrix <- featureMatrix
  
  if(!is.null(classesWhiteList)){
    classes <- intersect(x = classes, classesWhiteList)
  }
  if(progress)  incProgress(amount = 0, detail = paste("Searching for", length(classes), "/", length(classifierClasses), "classes")) else print(paste("Searching for", length(classes), "/", length(classifierClasses), "classes"))
  
  if(!smoothIntensities){
    matrix@x  [matrix@x   != 0] <- 1
  }
  #results__class      = list()
  #results__spectrum_class <- array(data = list(), dim = c(numberOfSpectra))
  results__class_spectrum <- array(data = list(), dim = c(numberOfClasses))
  
  lastOut <- proc.time()["user.self"]
  lastIdx <- 1
  
  #classIdx <- 1
  for(classIdx in seq_along(classes)){
    ## progress
    time <- proc.time()["user.self"]
    if(time - lastOut > 1){
      lastOut <- time
      precursorProgress <- (classIdx - lastIdx) / numberOfClasses * 0.7
      lastIdx <- classIdx
      if(progress)  incProgress(amount = precursorProgress,     detail = paste("Classification:", classIdx, "/", numberOfClasses)) else print(paste("Classification:", classIdx, "/", numberOfClasses))
    }
    
    class <- classes[[classIdx]]
    
    ###########################################
    ## classify it
    classifier <- classifiers_class[[which(classifierClasses==class)]]
    
    methodFunctions   <- trainingAlgorithm$method
    #methodName        <- trainingAlgorithm$methodName
    #parameters_train  <- trainingAlgorithm$params
    #paramsString      <- trainingAlgorithm$paramsString
    #algoName          <- trainingAlgorithm$algoName
    
    if(is.null(classifier$fragmentMassSelection)){
      fragmentMassSelection <- rep(x = TRUE, times = ncol(matrix))
    } else {
      fragmentMassSelection <- mappingSpectraToClassDf$ClassMassIndeces %in% classifier$fragmentMassSelection
    }
    
    parameters_test <- list()
    parameters_test$matrix_test <- matrix[, fragmentMassSelection]
    parameters_test$classifier  <- classifier$classifier
    
    ########################## do
    startTime <- Sys.time()
    predicted_scores <- tryCatch(expr = {
        do.call(what = methodFunctions$classify, args = parameters_test)
      }, error = function(e) {
        #error <- e
        print(paste("####################"))
        print(paste("####################", e))
        print(paste("####################"))
        return(NULL)
      }
    )
    endTime <- Sys.time()
    time <- difftime(endTime, startTime, units = "secs")[[1]]
    
    ## select by fdr
    #epsilon     <- .Machine$double.eps ^ 0.5 ## 1.490116e-08 ## .Machine$double.eps = 2.220446e-16
    #quantileIdx <- which(classifier$quantiles == (1 - fdrThreshold))
    #quantileIdx <- which(mapply(FUN = function(x, y) {isTRUE(all.equal(x, y))}, classifier$quantiles, 1 - fdrThreshold))
    quantileIdx <- which(sapply(X = classifier$quantiles, FUN = function(quantile) {isTRUE(all.equal(quantile, 1 - fdrThreshold))}))
    if(length(quantileIdx) == 0)  stop(paste("Quantile for fdrThreshold", fdrThreshold, "not there"))
    if(length(quantileIdx) >  1)  stop(paste("Quantile for fdrThreshold", fdrThreshold, "ambiguous"))
    scoreThreshold <- classifier$quantilesValuesNegative[[quantileIdx]]
    
    spectrumIdsPredicted <- which(predicted_scores >= scoreThreshold)
    results__class_spectrum[[classIdx]] <- sapply(X = spectrumIdsPredicted, FUN = function(spectrumId){
      pValue <- 1 - classifier$quantiles[[max(which(
        classifier$quantilesValuesNegative <= predicted_scores[[spectrumId]]
      ))]]
      return(pValue)
    })
    names(results__class_spectrum[[classIdx]]) <- as.character(spectrumIdsPredicted)
  }## substance class
  
  #plot(sapply(1:numberOfSpectra, function(x){length(results__spectrum_class[[x]])}))
  #plot(sapply(1:numberOfClasses, function(x){length(results__class_spectrum[[x]])}))
  #cbind(substr(classes,1,175),sapply(1:numberOfClasses, function(x){length(results__class_spectrum[[x]])}))
  #
  #plot(sapply(1:numberOfClasses, function(x){min(results__class_spectrum[[x]])}))
  
  classToSpectra_class <- lapply(X = results__class_spectrum, FUN = function(x){
    if(length(x)==0){      return(NULL)
    } else {             return(sort(x))
    }
  })
  names(classToSpectra_class) <- classes
  
  
  # List of 13
  # + $ classifierName         : chr "library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE"
  # + $ numberOfSpectra        : int 1355
  # + $ numberOfPositiveSpectra: int 22
  # + $ numberOfNegativeSpectra: int 1333
  # + $ class                  : chr "Organic compounds; Organic acids and derivatives; Carboxylic acids and derivatives; Amino acids, peptides, and "| __truncated__
  # + $ fragmentMasses         : num [1:12207] -970 -969 -965 -965 -963 ...
  # - $ fragmentMassSelection  : indeces
  # + $ classOfClass           : chr "ChemOnt_SubstanceClass"
  # - $ maximumNumberOfScores               : int 10000
  # - $ removeRareFragments                 : logi FALSE
  # - $ mergeDuplicatedSpectra              : logi TRUE
  # - $ takeSpectraWithMaximumNumberOfPeaks : logi FALSE
  # - $ minimumNumberOfPosSpectraPerClass   : int 10
  # - $ minimumNumberOfNegSpectraPerClass   : int 10
  # - $ numberOfDataSetDecompositions       : int 10
  # - $ proportionTraining                  : num 0.7
  # - $ unitResolution                      : logi FALSE
  # + $ AUC                                 : num 0.956
  # + $ AUC_PR                              : num 0.755               ########################### new ###########################
  # + $ TPR_for_FPR_of_5Percent             : num 0.841
  # - $ TNR_for_FNR_of_5Percent             : num 0.715
  # - $ classifier                          :List of 20
  # - ..$ method      : chr "binda"
  # - ..$ modelInfo   :List of 14
  # - ..$ modelType   : chr "Classification"
  # - ..$ results     :'data.frame':	15 obs. of  9 variables:
  # - ...
  # + $ frequentFragments                   : Named num [1:538] 0.328 0.314 0.299 0.277 0.212 ...
  # - ..- attr(*, "names")= chr [1:538] "151.003240136054" "255.03036741573" "285.041143269231" "284.032276344086" ...
  # + $ characteristicFragments             : Named num [1:492] 0.311 0.305 0.295 0.271 0.203 ...
  # - ..- attr(*, "names")= chr [1:492] "151.003240136054" "255.03036741573" "285.041143269231" "284.032276344086" ...
  # - $ minimumFrequency                    : num 0.01
  # - $ frequentFragmentsMeanIntensity      : Named num [1:538] 0.362 0.309 0.647 0.528 0.243 ...
  # - ..- attr(*, "names")= chr [1:538] "151.003240136054" "255.03036741573" "285.041143269231" "284.032276344086" ...
  # - $ characteristicFragmentsMeanIntensity: Named num [1:492] 0.362 0.309 0.647 0.528 0.243 ...
  # - ..- attr(*, "names")= chr [1:492] "151.003240136054" "255.03036741573" "285.041143269231" "284.032276344086" ...
  # - $ importantFragments                  : Named num [1:35] 100 97.8 94.6 86 81.7 ...
  # - ..- attr(*, "names")= chr [1:35] "151.003240136054" "255.03036741573" "285.041143269231" "284.032276344086" ...
  # - $ quantiles                           : num [1:1001] 0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 ...
  # - $ quantilesValuesPositive             : num [1:1001] 1.8e-06 1.8e-06 1.8e-06 1.8e-06 1.8e-06 5.2e-06 5.2e-06 5.2e-06 9.4e-06 9.4e-06 ...
  # - $ quantilesValuesNegative             : num [1:1001] 9.0e-07 9.0e-07 9.0e-07 1.8e-06 1.8e-06 1.8e-06 1.8e-06 1.8e-06 1.8e-06 1.8e-06 ...
  # - $ positiveScores                      : Named num [1:410] 1.80e-06 5.20e-06 9.40e-06 1.21e-05 1.57e-05 1.71e-05 1.77e-05 2.56e-05 2.59e-05 2.59e-05 ...
  # - ..- attr(*, "names")= chr [1:410] "5992" "6534" "5162" "6605" ...
  # - $ negativeScores                      : Named num [1:3960] 9.0e-07 9.0e-07 9.0e-07 9.0e-07 9.0e-07 9.0e-07 9.0e-07 9.0e-07 1.0e-06 1.8e-06 ...
  # - ..- attr(*, "names")= chr [1:3960] "3127" "5623" "5626" "5659" ...
  # + classifier$alternativeSubstanceClasses <- alternativeSubstanceClasses               ########################### new ###########################
  # + classifier$differentSubstanceClasses   <- differentSubstanceClasses                 ########################### new ###########################
  # + classifier$importantFragments <- importance                                         ########################### new ###########################
  # / $ algorithm
  # - ..$ method      :List of 2
  # + ..$ methodName  : chr "ColSums"
  # - ..$ params      :List of 3
  # + ..$ paramsString: chr "smoothIntensities=FALSE, classWeights=FALSE, modelName=binda"
  # + ..$ algoName    : chr "method=caret; smoothIntensities=FALSE, classWeights=FALSE, modelName=binda"
  
  properties_class <- lapply(X = classifiers_class, FUN = function(x){
    c(
    x[c(
      "classifierName",          #: chr "library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE"
      "numberOfSpectra",         #: int 1355
      "numberOfPositiveSpectra", #: int 20
      "numberOfNegativeSpectra", #: int 1335
      #algorithm               #:List of 5
      "class",                   #: chr "Organic compounds; Alkaloids and derivatives"
      "fragmentMasses",          #: num [1:12207] -970 -969 -965 -965 -963 ...
      "classOfClass",            #: chr "ChemOnt_SubstanceClass"
      "frequentFragments",       #: Named num [1:29] 0.571 0.286 0.214 0.214 0.214 ...
      "characteristicFragments", #: Named num [1:28] 0.539 0.275 0.203 0.203 0.143 ...
      "importantFragments",      #
      "alternativeSubstanceClasses", #: chr [1:12]
      "differentSubstanceClasses",   #: chr [1:8]
      "AUC",                      #: num 0.912
      "AUC_PR",                   #: num 0.xxx
      "TPR_for_FPR_of_5Percent" #: num 0.667
    )],
    x[["algorithm"]]["algoName"],
    x[["algorithm"]]["methodName"],
    x[["algorithm"]]["paramsString"]
    )
  })
  
  remove <- sapply(classToSpectra_class, is.null)
  classToSpectra_class <- classToSpectra_class[!remove]
  properties_class     <- properties_class    [!remove]
  
  ## classToSpectra_class[[10]]
  # Named num [1:3] 0.032 0.041 0.042
  # - attr(*, "names")= chr [1:3] "2228" "2049" "2109"
  
  ## properties_class[[1]]
  # List of 13
  # $ classifierName         : chr "library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE"
  # $ numberOfSpectra        : int 1355
  # $ numberOfPositiveSpectra: int 29
  # $ numberOfNegativeSpectra: int 1326
  # $ class                  : chr "Organic compounds; Benzenoids; Benzene and substituted derivatives; Anilides"
  # $ fragmentMasses         : num [1:1668] -970 -969 -965 -963 -960 -959 -958 -956 -955 -954 ...
  # $ classOfClass           : chr "ChemOnt_SubstanceClass"
  # $ frequentFragments      : Named num [1:298] 0.276 0.241 0.207 0.207 0.207 ...
  # ..- attr(*, "names")= chr [1:298] "121" "144" "146" "134" ...
  # $ characteristicFragments: Named num [1:248] 0.218 0.217 0.193 0.177 0.167 ...
  # ..- attr(*, "names")= chr [1:248] "144" "121" "77" "146" ...
  # $ AUC                    : num 0.832
  # $ algoName               : chr "method=ColSums; smoothIntensities=FALSE"
  # $ methodName             : chr "ColSums"
  # $ paramsString           : chr "smoothIntensities=FALSE"
  
  
  resultObj <- list(
    "classToSpectra_class" = classToSpectra_class,
    "properties_class" = properties_class,
    "mappingSpectraToClassDf" = mappingSpectraToClassDf
  )
  return(resultObj)
}
getAlgorithms <- function(){
  methods <- list(
    ## scores
    "CosinusDistance"  = predict_CosinusDistance,
    "ColSums"          = colSums_classifier,
    "ColSumsPos"       = colSumsPos_classifier,
    "Prod"             = predict_Prod,
    "Jaccard"          = predict_Jaccard,
    "JaccardWeighted"  = predict_JaccardWeighted,
    ## classes
    #"LDA"              = predict_LDA,
    "Correlation"      = predict_Correlation#,
    #"RDA"              = predict_RDA,
    #"SVM"              = predict_SVM#,
    #"SOM"              = predict_SOM,
    #"XYF"              = predict_XYF,
    #"NeuralNet"        = predict_NeuralNet
  )
  
  params <- list(
    ## scores
    "CosinusDistance"  = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    "ColSums"          = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "ColSumsPos"          = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "Prod"             = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    "Jaccard"          = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    "JaccardWeighted"  = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    ## classes
    #"LDA"              = list(
    #  list(ratio=FALSE),
    #  list(ratio=TRUE)
    #),
    "Correlation"      = list(
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=TRUE)
    )#,
    #"RDA"              = list(
    #  list()
    #),
    #"SVM"              = list(
    #  list(ratio=FALSE, smoothIntensities=FALSE),
    #  list(ratio=TRUE,  smoothIntensities=FALSE),
    #  list(ratio=FALSE, smoothIntensities=TRUE),
    #  list(ratio=TRUE,  smoothIntensities=TRUE)
    #)
    #"SOM"              = list(
    #  list()
    #),
    #"XYF"              = list(
    #  list()
    #),
    #"NeuralNet"        = list(
    #  list()
    #)
  )
  
  algorithms <- list()
  for(methodIdx in seq_along(methods)){
    methodName <- names(methods)[[methodIdx]]
    paramsLists <- params[[methodName]]
    for(paramsListIdx in seq_along(paramsLists)){
      paramsList <- paramsLists[[paramsListIdx]]
      paramsString <- paste(names(paramsList), paramsList[names(paramsList)], sep = "=", collapse = ", ")
      algoName <- paste("method=", methodName, "; ", paramsString, "", sep = "")
      
      algorithms[[length(algorithms) + 1]] <- list(
        "method" = methods[[methodIdx]],
        "methodName" = methodName,
        "params" = paramsList,
        "paramsString" = paramsString,
        "algoName" = algoName
      )
    }
  }
  
  return(algorithms)
}
getAvailableClassifiers_old <- function(resultFolderForClassifiers){
  classifierFilePaths <- list.files(path = resultFolderForClassifiers, recursive = FALSE, pattern = "^.*.RData$", include.dirs = FALSE, full.names = TRUE)
  classifierFiles     <- basename(classifierFilePaths)
  
  pattern <- paste(
    "^Classifier_",
    
    "library=", "(?<library>.+)", "_",
    "Class=", "(?<type>.+)", "_",
    "AltSC=", "(?<altSubstClass>TRUE|FALSE)", "_", 
    "method=", "(?<methodName>.+)", "_",
    "(?<parameters>.+)", "_",
    "NumClasses=", "(?<classes>\\d+)", 
    
    ".RData$", 
    sep = ""
  )
  
  result <- regExExtraction(pattern, classifierFiles)
  
  if(length(classifierFiles) > 0){
    numberOfProperties <- length(names(result[[1]]))
    propertyNames <- names(result[[1]])
    
    df <- data.frame(
      t(matrix(data = unlist(result), nrow = numberOfProperties)),
      stringsAsFactors = FALSE
    )
    
    ## annotate
    rownames(df) <- classifierFilePaths
    colnames(df) <- propertyNames
    
    ## reorder columns
    df <- df[, c(1,6,2,3,4,5)]
  } else {
    return(
      data.frame(matrix(nrow=0, ncol=6))
    )
  }
  return(df)
}

getClassifierProperties <- function(propertiesFile){
  lines <- readLines(con = propertiesFile)
  linesSplitted <- strsplit(x = lines, split = " = ")
  
  tags   <- unlist(lapply(X = linesSplitted, FUN = function(x){x[[1]]}))
  values <- unlist(lapply(X = linesSplitted, FUN = function(x){x[[2]]}))
  propertiesList <- as.list(values)
  names(propertiesList) <- tags
  
  return(propertiesList)
}
## resultFolderForClassifiers <- "/home/htreutle/Code/Java/MetFam/inst/data/classifiers"
## resultFolderForClassifiers <- "/home/htreutle/Code/Java/MetFam/inst/data/classifier"
getAvailableClassifiers <- function(resultFolderForClassifiers){
  
  resultFolderForClassifiers <- gsub(x = resultFolderForClassifiers, pattern = "app/data/classifier", replacement = "inst/data/classifier")
  
  classifierFilePaths <- list.files(path = resultFolderForClassifiers, recursive = FALSE, pattern = "^.*_Classifier.RData$", include.dirs = FALSE, full.names = TRUE)
  classifierFiles     <- basename(classifierFilePaths)
  resultFiles         <- gsub(x = classifierFiles, pattern = "_Classifier.RData$", replacement = "_Results.tsv")
  propertiesFiles     <- gsub(x = classifierFiles, pattern = "_Classifier.RData$", replacement = "_Properties.txt")
  
  propertiesListsList <- list()
  theseClassifiers <- rep(x = TRUE, times = length(classifierFilePaths))
  for(idx in seq_along(propertiesFiles)){
    propertiesFile <- paste(resultFolderForClassifiers, "/", propertiesFiles[[idx]], sep = "")
    propertiesList <- getClassifierProperties(propertiesFile)
    propertiesListsList[[length(propertiesListsList)+1]] <- propertiesList
  }
  
  classifierFilePaths <- classifierFilePaths[theseClassifiers]
  classifierFiles     <- classifierFiles    [theseClassifiers]
  resultFiles         <- resultFiles        [theseClassifiers]
  propertiesFiles     <- propertiesFiles    [theseClassifiers]
  
  if(length(classifierFiles) > 0){
    numberOfProperties <- length(propertiesListsList[[1]])
    propertyNames <- names(propertiesListsList[[1]])
    
    df <- data.frame(
      t(matrix(data = unlist(propertiesListsList), nrow = numberOfProperties)),
      stringsAsFactors = FALSE
    )
    
    ## annotate
    rownames(df) <- classifierFilePaths
    colnames(df) <- propertyNames
    
    libraryNames <- gsub(x = df$library, pattern = "^\\d\\d\\d\\d\\-\\d\\d\\-\\d\\d_\\d\\d:\\d\\d:\\d\\d_", replacement = "")
    libraryNames <- gsub(x = libraryNames, pattern = "\\.[a-zA-Z]{3,4}$", replacement = "")
    df[, "library"]    <- libraryNames
    #df[, "Class"]      <- ifelse(test = df$processSubstanceClasses, yes = rep(x = "SubstanceClass", times = nrow(df)), no = rep(x = "Substituent", times = nrow(df)))
    df[, "Class"]      <- df$classOfClass
    df[, "Resolution"] <- ifelse(test = df$unitResolution, yes = rep(x = "Low", times = nrow(df)), no = rep(x = "High", times = nrow(df)))
    df[, "Classifier"] <- df[, "methodName"]
    df[, "FilePath"]   <- classifierFilePaths
    ## methodName
    
    ## reorder columns
    dfShow       <- df[, c("library", "Class", "Resolution", "Classifier")]
    dfProperties <- df
  } else {
    dfShow       <- data.frame(matrix(nrow=0, ncol= 4))
    dfProperties <- data.frame(matrix(nrow=0, ncol=18))
  }
  
  resultObj <- list()
  resultObj$availableClassifiersDf           <- dfShow
  resultObj$availableClassifiersDfProperties <- dfProperties
  
  return(resultObj)
}

evaluatePutativeMetaboliteFamiliesOfPrecursorSet <- function(dataList, precursorSet, classToSpectra_class){
  if(FALSE){
    dataList_ <<- dataList
    precursorSet_ <<- precursorSet
    classToSpectra_class_ <<- classToSpectra_class
  }
  if(FALSE){
    dataList <<- dataList_
    precursorSet <<- precursorSet_
    classToSpectra_class <<- classToSpectra_class_
  }
  
  ## fetch hits of interest
  classesAll  <- character()
  featuresAll <- integer()
  pValuesAll  <- numeric()
  
  for(classIdx in seq_along(classToSpectra_class)){
    featuresHere <- names(classToSpectra_class[[classIdx]])
    indeces <- which(featuresHere %in% precursorSet)
    if(length(indeces) > 0){
      classes  <- rep(x = names(classToSpectra_class)[[classIdx]], times = length(indeces))
      features <- featuresHere[indeces]
      pValues  <- unname(classToSpectra_class[[classIdx]])[indeces]
      classesAll  <- c(classesAll, classes)
      featuresAll <- c(featuresAll, features)
      pValuesAll  <- c(pValuesAll, pValues)
    }
  }
  detailDf <- data.frame("MS1_feature"=featuresAll, "Class"=classesAll, "pValue"=pValuesAll, stringsAsFactors = FALSE)
  
  if(nrow(detailDf) == 0){
    return(list("overviewDf" = data.frame("Class"=character(), "pValue"=numeric(), "ProportionInPercent" = numeric(), stringsAsFactors = FALSE), "detailDf" = detailDf))
  }
  if(length(precursorSet) == 1){
    ## single precursor
    overviewDf <- cbind(detailDf$Class, detailDf$pValue, rep(x = 100, times = nrow(detailDf)))
    colnames(overviewDf) <- c("Class", "pValue", "ProportionInPercent")
    overviewDf <- as.data.frame(x = overviewDf, stringsAsFactors = FALSE)
    #printPutativeMetaboliteFamilies <- paste(detailDf$Class, " (pValue=", detailDf$pValue, ")", sep = "")
  } else {
    ## multiple precursors: do statistics
    potentialMetaboliteFamilies <- unique(detailDf$Class)
    potentialMetaboliteFamiliesWithSuperClasses <- sort(unique(unlist(
      lapply(X = strsplit(x = potentialMetaboliteFamilies, split = "; "), FUN = function(x){
        x <- x[x != "NA"]
        class <- sapply(X = seq_len(length(x)), FUN = function(y){
          paste(x[1:y], collapse = "; ")
        })
        return(class)
      })
    )))
    precursorHits <- unname(sapply(X = potentialMetaboliteFamiliesWithSuperClasses, FUN = function(class){
      classHits <- grepl(pattern = paste("^", class, sep = ""), x = detailDf$Class)
      count <- length(unique(detailDf$MS1_feature[classHits]))
      return(count)
    }))
    retainClass <- precursorHits >= (length(precursorSet)/2)
    #retainClass <- precursorHits >= 0
    
    if(sum(retainClass) == 0)
      #return(data.frame("Class"=character(), "pValue"=numeric(), "ProportionInPercent" = numeric(), stringsAsFactors = FALSE))
      return(list("overviewDf" = data.frame("Class"=character(), "pValue"=numeric(), "ProportionInPercent" = numeric(), stringsAsFactors = FALSE), "detailDf" = detailDf))
    
    potentialMetaboliteFamiliesWithSuperClasses <- potentialMetaboliteFamiliesWithSuperClasses[retainClass]
    precursorHits <- precursorHits[retainClass]
    
    ## remove classes without relevance
    detailDf <- detailDf[detailDf$Class %in% potentialMetaboliteFamiliesWithSuperClasses, ]
    
    proportionPercent <- precursorHits / length(precursorSet) * 100
    medianPValue <- unlist(unname(sapply(X = potentialMetaboliteFamiliesWithSuperClasses, FUN = function(class){
      classHits <- grepl(pattern = paste("^", class, sep = ""), x = detailDf$Class)
      value <- median(detailDf$pValue[classHits])
      return(value)
    })))
    
    overviewDf <- data.frame(
      "Class" = potentialMetaboliteFamiliesWithSuperClasses, 
      "pValue" = medianPValue,
      "ProportionInPercent" = format(x = proportionPercent, digits=3, nsmall=1), 
      stringsAsFactors = FALSE
    )
    
    #printPutativeMetaboliteFamilies <- paste(overviewDf$Proportion, "% ", overviewDf$Class, sep = "")
  }
  
  returnObj <- list(
    "overviewDf" = overviewDf,
    "detailDf"   = detailDf
  )
  #return(printPutativeMetaboliteFamilies)
  return(returnObj)
}
evaluatePutativeMetaboliteFamiliesOfPrecursorSet_old <- function(dataList, precursorSet, classToSpectra_class){
  if(FALSE){
    dataList_ <<- dataList
    precursorSet_ <<- precursorSet
    classToSpectra_class_ <<- classToSpectra_class
  }
  if(FALSE){
    dataList <<- dataList_
    precursorSet <<- precursorSet_
    classToSpectra_class <<- classToSpectra_class_
  }
  
  ## hits: list of 
  hitLists <- list()
  for(precursorIndex in precursorSet){
    hitLists[[length(hitLists)+1]] <- list()
    for(classIdx in seq_along(classToSpectra_class)){
      if(!(precursorIndex %in% names(classToSpectra_class[[classIdx]])))
        next
      
      idx <- which(names(classToSpectra_class[[classIdx]]) == precursorIndex)
      hitLists[[length(hitLists)]][[length(hitLists[[length(hitLists)]])+1]] <- c(
        "Class"  = names(classToSpectra_class)[[classIdx]],
        "pValue" = unname(classToSpectra_class[[classIdx]])[[idx]]
      )
    }
  }
  
  if(is.null(unlist(hitLists))){
    return("")
  }
  if(length(precursorSet) == 1){
    ## single precursor
    df <- data.frame(t(matrix(data = unlist(hitLists[[1]]), nrow = 2)), stringsAsFactors = F)
    colnames(df) <- c("Class", "pValue")
    printPutativeMetaboliteFamilies <- paste(df$Class, " (pValue=", df$pValue, ")", sep = "")
  } else {
    ## multiple precursors: do statistics
    
    potentialMetaboliteFamilies <- unlist(hitLists)
    potentialMetaboliteFamilies <- potentialMetaboliteFamilies[seq(from = 1, to = length(potentialMetaboliteFamilies), by = 2)]
    potentialMetaboliteFamilies <- unique(potentialMetaboliteFamilies)
    potentialMetaboliteFamiliesWithSuperClasses <- sort(unique(unlist(
      lapply(X = strsplit(x = potentialMetaboliteFamilies, split = "; "), FUN = function(x){
        x <- x[x != "NA"]
        sapply(X = seq_len(length(x)), FUN = function(y){
          paste(x[1:y], collapse = "; ")
        })
      })
    )))
    precursorHitHits <- unname(sapply(X = potentialMetaboliteFamiliesWithSuperClasses, FUN = function(class){
      count <- sum(unlist(lapply(X = hitLists, FUN = function(y){
        classes <- unlist(lapply(X = y, FUN = function(z){z[[1]]}))
        any(grepl(pattern = paste("^", class, sep = ""), x = classes))
      })))
    }))
    retainClass <- rep(x = FALSE, times = length(potentialMetaboliteFamiliesWithSuperClasses))
    for(idx in seq_along(potentialMetaboliteFamiliesWithSuperClasses)){
      these <- grepl(pattern = paste("^", potentialMetaboliteFamiliesWithSuperClasses[[idx]], sep = ""), x = potentialMetaboliteFamiliesWithSuperClasses[-idx])
      if(all(precursorHitHits[-idx][these] < precursorHitHits[[idx]]))
        retainClass[[idx]] <- TRUE
    }
    
    ## at least 50 % of the precursors must be hits of each class
    retainClass <- retainClass & precursorHitHits >= (length(precursorSet)/2)
    
    potentialMetaboliteFamiliesWithSuperClasses <- potentialMetaboliteFamiliesWithSuperClasses[retainClass]
    precursorHitHits <- precursorHitHits[retainClass]
    proportionPercent <- precursorHitHits / length(precursorSet) * 100
    
    df <- data.frame(
      "Class" = potentialMetaboliteFamiliesWithSuperClasses, 
      "Proportion" = format(x = proportionPercent, digits=3, nsmall=1), 
      #"pValue" = 
      stringsAsFactors = F
    )
    
    printPutativeMetaboliteFamilies <- paste(df$Proportion, "% ", df$Class, sep = "")
  }
  return(printPutativeMetaboliteFamilies)
}

metaboliteFamilyVersusClass <- function(dataList, precursorSet, classToSpectra_class, properties_class, classifierClass, mappingSpectraToClassDf, addClassifierConsensusSpectrum){
  returnObj <- getSpectrumStatistics(dataList = dataList, precursorSet = precursorSet)
  masses_spec <- returnObj$fragmentMasses
  fragmentCounts_spec <- returnObj$fragmentCounts
  frequency_spec <- fragmentCounts_spec / length(precursorSet)
  
  if(all(addClassifierConsensusSpectrum, classifierClass %in% names(classToSpectra_class))){
    ## Plot spectrum vs consensus spectrum
    classIdx <- which(classifierClass == names(classToSpectra_class))
    classProperties         <- properties_class[[classIdx]]
    frequentFragments       <- classProperties$frequentFragments
    characteristicFragments <- classProperties$characteristicFragments
    
    ## class statistics for class plot
    returnObj <- preprocessClassPlot(frequentFragments, characteristicFragments)
    masses_class    <- returnObj$masses_class
    frequency_class <- returnObj$frequency_class
    #colors_class    <- returnObj$colors_class
    
    ## match spec to class
    returnObj <- preprocessSpectrumVsClassPlot(dataList, precursorSet, masses_class, mappingSpectraToClassDf, "Counts")
    masses_spec <- returnObj$masses_spec
    frequency_spec <- returnObj$intensity_spec
    colors_spec <- returnObj$colors_spec
    numberOfMatchingMasses <- returnObj$numberOfMatchingMasses
    matchingMassRowIndeces <- returnObj$matchingMassRowIndeces
    frequency_spec <- frequency_spec / length(precursorSet)
    
    colors_class    <- rep(x = "grey", times = length(masses_class))
    colors_class[masses_class %in% mappingSpectraToClassDf$ClassMasses[matchingMassRowIndeces]] <- "black"
  } else {
    colors_spec <- rep(x = "black", times = length(masses_spec))
    masses_class    <- NULL
    frequency_class <- NULL
    colors_class    <- NULL
  }
  
  returnObj <- list(
      masses_spec     = masses_spec, 
      intensity_spec  = frequency_spec, 
      colors_spec     = colors_spec, 
      masses_class    = masses_class, 
      frequency_class = frequency_class, 
      colors_class    = colors_class
  )
  return(returnObj)
}
preprocessClassPlot <- function(frequentFragments, characteristicFragments){
  frequentMasses       <- as.numeric(names(frequentFragments)) 
  characteristicMasses <- as.numeric(names(characteristicFragments))
  
  masses_class <- unique(c(frequentMasses, characteristicMasses))
  frequency_class <- rep(x = 0.01, times = length(masses_class))
  frequency_class[match(x = frequentMasses, table = masses_class)] <- unname(frequentFragments)
  
  characteristics_class <- rep(x = 0., times = length(masses_class))
  characteristics_class[match(x = characteristicMasses, table = masses_class)] <- unname(characteristicFragments)
  classDataColorMapFragmentData  <- makecmap(
    x = c(0, 1), n = 100, 
    colFn = colorRampPalette(c('grey', 'black'))
  )
  colors_class <- cmap(x = characteristics_class, map = classDataColorMapFragmentData)
  
  returnObj <- list(
    masses_class = masses_class,
    frequency_class = frequency_class,
    colors_class = colors_class
  )
  return(returnObj)
}
preprocessSpectrumVsClassPlot <- function(dataList, precursorIndeces, masses_class, mappingSpectraToClassDf, yType){
  if(length(precursorIndeces) == 1){
    resultObj   <- getMS2spectrumInfoForPrecursor(dataList = dataList, precursorIndex = precursorIndeces)
    masses_spec <- resultObj$fragmentMasses
    switch(yType, 
           "Counts"={
             intensity_spec <- rep(x = 1, times = length(masses_spec))
           },
           "Intensity"={
             intensity_spec <- resultObj$fragmentAbundances
           },
           {
             stop(paste("Unknown yType", yType))
           }
    )
  } else {
    returnObj   <- getSpectrumStatistics(dataList = dataList, precursorSet = precursorIndeces)
    masses_spec <- returnObj$fragmentMasses
    switch(yType, 
           "Counts"={
             intensity_spec <- returnObj$fragmentCounts
           },
           "Intensity"={
             intensity_spec <- resultObj$fragmentAbundances
           },
           {
             stop(paste("Unknown yType", yType))
           }
    )
  }
  
  tolerance <- .Machine$double.eps ^ 0.5 ## default in function all.equal
  
  
  matchingMassRowIndeces <- which(
    as.character(mappingSpectraToClassDf$SpectraMasses) %in% as.character(masses_spec) &
    as.character(mappingSpectraToClassDf$ClassMasses  ) %in% as.character(masses_class)
  )
  #matchingMassRowIndeces <- which(
  #  apply(X = outer(X = mappingSpectraToClassDf$SpectraMasses, Y = masses_spec , FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)}) &
  #  apply(X = outer(X = mappingSpectraToClassDf$ClassMasses  , Y = masses_class, FUN = "-"), MARGIN = 1, FUN = function(x){any(abs(x) <= tolerance)})
  #)
  
  #matchingMassRowIndeces <- which(
  #  mappingSpectraToClassDf$SpectraMasses %in% masses_spec & 
  #  mappingSpectraToClassDf$ClassMasses   %in% masses_class
  
  specIndeces  <- match(x = mappingSpectraToClassDf$SpectraMasses[matchingMassRowIndeces], table = masses_spec )
  #classIndeces <- match(x = mappingSpectraToClassDf$ClassMasses  [matchingMassRowIndeces], table = masses_class)
  
  colors_spec <- rep(x = "grey", times = length(masses_spec))
  #colors_spec[specIndeces] <- colors_class[classIndeces]
  colors_spec[specIndeces] <- "black"
  
  numberOfMatchingMasses <- length(matchingMassRowIndeces)
  
  returnObj <- list(
    masses_spec = masses_spec,
    intensity_spec = intensity_spec,
    colors_spec = colors_spec,
    numberOfMatchingMasses = numberOfMatchingMasses,
    matchingMassRowIndeces = matchingMassRowIndeces
  )
  return(returnObj)
}
