test_that("trivial", {
    expect_equal(1, 1)
})

test_that("exampledata", {
  filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
  data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping=TRUE,
                                   mzDeviationInPPM_precursorDeisotoping=10, mzDeviationAbsolute_precursorDeisotoping=0.01,
                                   maximumRtDifference=0.05,
                                   progress=FALSE)

  ## Test dimensions  
  expect_equal(nrow(data$dataFrame1), 5403)
  expect_equal(ncol(data$dataFrame1),   20)

  ## Test some values
  expect_true(all(summary(t(data$dataFrame1[1,c("TRI03", "TRI02", "TRI01", "LVS03", "LVS02", "LVS01")])) == c("Min.   : 236.0  ", "1st Qu.: 306.2  ", "Median : 357.5  ", "Mean   : 501.0  ", "3rd Qu.: 554.2  ", "Max.   :1146.0  ")))
  expect_true(all(round(summary(data$dataFrame1[, "TRI01"])) == c(0,     809,    1537,    9818,    3207, 4407926)))
})

test_that("MS-Dial 4.X", {
#   filePeakMatrix <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")
#     data <- parsePeakAbundanceMatrix(filePeakMatrix, doPrecursorDeisotoping, 
# 	mzDeviationInPPM_precursorDeisotoping, mzDeviationAbsolute_precursorDeisotoping, 
# 	maximumRtDifference, 
# 	progress=FALSE)
# 
#     expect_equal(ncol(data$dataFrame), 41)
})

test_that("MSP reading works", {
  
  if(FALSE) {
  fileSpectra <- system.file("extdata/showcase/Metabolite_profile_showcase.msp", package = "MetFamily")
  #load(system.file("extdata/testdata/readMSPreturnObj.Rdata", package = "MetFamily"))
  
  returnObj <- parseMSP(fileSpectra = fileSpectra, 
                        minimumIntensityOfMaximalMS2peak = 2000, 
                        minimumProportionOfMS2peaks = 0.05, 
                        neutralLossesPrecursorToFragments = TRUE,
                        neutralLossesFragmentsToFragments = FALSE,
                        progress = FALSE)

  
  expect_equal(returnObj$numberOfSpectra, 2640)
  expect_equal(returnObj$numberOfSpectraOriginal, 5824)
  expect_equal(returnObj$numberOfMS2PeaksOriginal, 145973)
  expect_equal(returnObj$numberOfMS2PeaksWithNeutralLosses, 68738)
  expect_equal(returnObj$numberOfMS2PeaksAboveThreshold, 34369)
  expect_equal(returnObj$numberOfMS2PeaksBelowThreshold, 109797)
  expect_equal(returnObj$numberOfTooHeavyFragments, 1807)
  expect_equal(returnObj$numberOfSpectraDiscardedDueToNoPeaks, 15)
  expect_equal(returnObj$numberOfSpectraDiscardedDueToMaxIntensity, 3163)
  expect_equal(returnObj$numberOfSpectraDiscardedDueToTooHeavy, 5)
  expect_equal(length(returnObj$precursorMz), 2640)
  expect_equal(returnObj$precursorMz[1], 85)
  expect_equal(length(returnObj$precursorRt), 2640)
  expect_equal(returnObj$precursorRt[1], 10.98)
  
  }
  
  # > str(returnObj, max.level = 1)
  # List of 14
  # $ fileSpectra                              : chr "/tmp/RtmpL82I4E/ff16bb92ac5e87cb202eec2e/0.msp"
  # $ spectraList                              :List of 2640
  # $ numberOfSpectra                          : int 2640
  # $ numberOfSpectraOriginal                  : int 5824
  # $ numberOfMS2PeaksOriginal                 : num 145973
  # $ numberOfMS2PeaksWithNeutralLosses        : num 68738
  # $ numberOfMS2PeaksAboveThreshold           : num 34369
  # $ numberOfMS2PeaksBelowThreshold           : num 109797
  # $ numberOfTooHeavyFragments                : num 1807
  # $ numberOfSpectraDiscardedDueToNoPeaks     : num 15
  # $ numberOfSpectraDiscardedDueToMaxIntensity: num 3163
  # $ numberOfSpectraDiscardedDueToTooHeavy    : num 5
  # $ precursorMz                              : num [1:2640] 85 85 85 85 85 ...
  # $ precursorRt                              : num [1:2640] 10.98 7.42 12.4 16.26 13 ...
  # 
  
  
  
  ##
  ## returnObj contains the following elements:
  ##
  
  # returnObj <- list()
  # returnObj$fileSpectra <- fileSpectra
  # returnObj$spectraList <- list()
  # returnObj$numberOfSpectra <- 0
  # returnObj$numberOfMS2PeaksOriginal <- 0
  # returnObj$numberOfMS2PeaksWithNeutralLosses <- 0
  # returnObj$numberOfMS2PeaksAboveThreshold <- 0
  # returnObj$numberOfMS2PeaksBelowThreshold <- 0
  # returnObj$precursorMz <- vector(mode = "numeric")
  # returnObj$precursorRt <- vector(mode = "numeric")
  
})

test_that("Project reading works", {
    fileName <- system.file("extdata/showcase/Project_file_showcase_annotated.csv.gz", package = "MetFamily")
    readDataList <- readClusterDataFromProjectFile(
      file = fileName,
      progress = FALSE)
    #load("inst/extdata/testdata/readClusterDataFromProjectFile.Rdata")
    
    expect_equal(length(readDataList), 59)
    expect_equal(nrow(readDataList$dataFrameHeader), 3)
    expect_equal(ncol(readDataList$dataFrameHeader), 14570)
    expect_equal(length(readDataList$importParameterSet), 19)

  # 
  # > str(dataList, max.level = 1)
  # List of 59
  # $ dataFrameHeader                        :'data.frame':	3 obs. of  14570 variables:
  #   $ dataFrameMS1Header                     :'data.frame':	3 obs. of  23 variables:
  #   $ dataFrameInfos                         :'data.frame':	2585 obs. of  23 variables:
  #   $ importParameterSet                     :List of 19
  # $ numberOfPrecursors                     : int 2585
  # $ numberOfDuplicatedPrecursors           : int 0
  # $ sampleClasses                              : chr [1:2] "TRI" "LVS"
  # $ columnGroupLabels                      : chr [1:3, 1:2] "TRI" "TRI" "TRI" "LVS" ...
  # ..- attr(*, "dimnames")=List of 2
  # $ groupSampleDataFrame                   :'data.frame':	6 obs. of  4 variables:
  #   $ metaboliteProfileColumnNames           : chr [1:23] "m/z" "RT" "Annotation" "Alignment ID" ...
  # $ tagsSector                             : chr [1:23] "ID" "ID" "AnnotationColors={AS=#0000FF, SQT-glucosides=#FF0000}" "" ...
  # $ fragmentMasses                         : num [1:14547] 69 71 71 71 73 ...
  # $ fragmentFrequency                      : int [1:14547] 19 149 3 4 15 77 1 6 11 1 ...
  # $ fragmentAbundance                      : num [1:14547] 0.1162 0.0899 0.131 0.1033 0.1987 ...
  # $ minimumMass                            : num -1008
  # $ maximumMass                            : num 1198
  # $ precursorLabels                        : chr [1:2585] "    85.005 /   7.42" "    85.005 / 10.98" "    85.005 /   8.58" "    85.005 / 12.40" ...
  # $ dataFrameMeasurements                  :'data.frame':	2585 obs. of  13 variables:
  #   $ meanAllMax                             : num 2441679
  # $ logFoldChangeMax                       : num 10.5
  # $ logAbsMax                              : num 6.7
  # $ colorMatrixDataFrame                   : chr [1:2585, 1:13] "#8AFF00" "#E0FF00" "#9CFF00" "#EBFF00" ...
  # ..- attr(*, "dimnames")=List of 2
  # $ colorMapAbsoluteData                   :List of 6
  # $ colorMapLogFoldChange                  :List of 6
  # $ dataColumnsNameFunctionFromGroupName   :function (group, sampleNamesToExclude)  
  #   ..- attr(*, "srcref")= 'srcref' int [1:8] 666 43 668 3 43 3 3786 3788
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ dataColumnsNameFunctionFromGroupIndex  :function (groupIdx, sampleNamesToExclude)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 662 44 664 3 44 3 3782 3784
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ dataColumnsNameFunctionFromGroupNames  :function (sampleClasses, sampleNamesToExclude)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 670 44 674 3 44 3 3790 3794
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ groupNameFunctionFromDataColumnName    :function (dataColumnName, sampleNamesToExclude)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 676 42 682 3 42 3 3796 3802
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ lfcColumnNameFunctionFromString        :function (columnName)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 464 38 469 3 38 3 3584 3589
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ dataMeanColumnNameFunctionFromString   :function (columnName)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 470 43 473 3 43 3 3590 3593
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ dataColumnIndecesFunctionFromGroupIndex:function (groupIdx, sampleNamesToExclude)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 658 46 660 3 46 3 3778 3780
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ dataMeanColumnNameFunctionFromName     :function (group)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 800 42 802 3 42 3 3920 3922
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ dataMeanColumnNameFunctionFromIndex    :function (groupIdx)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 804 43 806 3 43 3 3924 3926
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ lfcColumnNameFunctionFromName          :function (groupOne, groupTwo)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 808 36 810 3 36 3 3928 3930
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ lfcColumnNameFunctionFromIndex         :function (groupIdxOne, groupIdxTwo)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 812 37 814 3 37 3 3932 3934
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ groupNameFromGroupIndex                :function (groupIdx)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 816 30 818 3 30 3 3936 3938
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ groupIdxFromGroupName                  :function (group)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 820 28 822 3 28 3 3940 3942
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ featureMatrix                          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  # $ featureIndeces                         :List of 2585
  # $ featureCount                           : num [1:2585] 1 1 1 1 1 1 1 1 1 1 ...
  # $ featureIndexMatrix                     : int [1:2585, 1:90] 15 15 15 15 15 15 15 15 15 15 ...
  # ..- attr(*, "dimnames")=List of 2
  # $ ms2_numberOfFragments                  : int [1:14547] 19 149 3 4 15 77 1 6 11 1 ...
  # $ ms2_averageAbundance                   : num [1:14547] 0.1162 0.0899 0.131 0.1033 0.1987 ...
  # $ ms2_masses                             : num [1:14547] 69 71 71 71 73 ...
  # $ colorMapFragmentData                   :List of 6
  # $ annotationColumnName                   : chr "Annotation"
  # $ annotationColorsName                   : chr "AnnotationColors"
  # $ annotationColumnIndex                  : int 3
  # $ annotationValueIgnore                  : chr "Ignore"
  # $ annotationColorIgnore                  : chr "red"
  # $ annoArrayOfLists                       :List of 2585
  # $ annoArrayIsArtifact                    : logi [1:2585] FALSE FALSE FALSE FALSE FALSE FALSE ...
  # $ annoPresentAnnotationsList             :List of 3
  # $ annoPresentColorsList                  :List of 3
  # $ orderColumnNames                       :function (groupSampleDataFrame, columnNames)  
  #   ..- attr(*, "srcref")= 'srcref' int [1:8] 685 23 691 3 23 3 3805 3811
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ excludedSamples                        :function (groupSampleDataFrame, sampleClasses = dataList$sampleClasses)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 695 22 700 3 22 3 3815 3820
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ includedSamples                        :function (groupSampleDataFrame, sampleClasses = dataList$sampleClasses)  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 702 22 707 3 22 3 3822 3827
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ includedGroups                         :function (groupSampleDataFrame, samples = dataList$groupSampleDataFrame[, "Sample"])  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 710 21 714 3 21 3 3830 3834
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   $ excludedGroups                         :function (groupSampleDataFrame, samples = dataList$groupSampleDataFrame[, "Sample"])  
  #     ..- attr(*, "srcref")= 'srcref' int [1:8] 716 21 718 3 21 3 3836 3838
  # .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x555c0c5c55d0> 
  #   
  #   
  #   
})


