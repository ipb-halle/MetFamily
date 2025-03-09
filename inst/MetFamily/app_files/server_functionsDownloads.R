
timeStampForFiles <- function(){
  timeStamp <- gsub(" ", "_", gsub(":", ".", Sys.time()))
  return(timeStamp)
}
createImportParameterSetExportFileName <- function(){
  fileProjectName <- dataList$importParameterSet$projectName
  fileProjectName <- gsub(" ", "_", gsub(":", ".", fileProjectName))
  fileName <- paste(timeStampForFiles(), "_", fileProjectName, "_import_parameters.txt", sep = "")
  return(fileName)
}
createExportMatrixName <- function(){
  fileName <- paste(timeStampForFiles(), "_selectedPrecursorMatrix.csv.gz", sep = "")
  return(fileName)
}
createConsensusSpectrumName <- function(annotation){
  fileName <- paste(timeStampForFiles(), "_ConsensusSpectrum_", gsub(pattern = " ", replacement = "_", x = annotation), ".csv", sep = "")
  return(fileName)
}
createClassifierAnnotationName <- function(annotation){
  fileName <- paste(timeStampForFiles(), "_ClassifierAnnotation_", gsub(pattern = " ", replacement = "_", x = annotation), ".csv", sep = "")
  return(fileName)
}
createMetaboliteFamilyProjectFileName <- function(annotation){
  fileName <- paste(timeStampForFiles(), "_", gsub(pattern = " ", replacement = "_", x = annotation), ".csv.gz", sep = "")
  return(fileName)
}
createExportImageName <- function(item, extension){
  fileName <- paste(timeStampForFiles(), "_", item, ".", extension, sep = "")
  return(fileName)
}
createExportDistanceMatrixName <- function(distanceMeasure){
  fileName <- paste(timeStampForFiles(), "_distanceMatrix_", distanceMeasure, ".csv", sep = "")
  return(fileName)
}
createExportMatrix <- function(precursorSet){
  ################################################################################
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
  
  for(chunkIdx in seq_len(numberOfChunks)){
    colStart <- 1 + (chunkIdx - 1) * chunkSize
    colEnd <- colStart + chunkSize - 1
    if(chunkIdx == numberOfChunks)
      colEnd <- numberOfColumns
    
    numberOfColumnsHere <- colEnd - colStart + 1
    numberOfRowsHere <- max(matrixRows)
    indeces <- matrixCols >= colStart & matrixCols <= colEnd
    
    fragmentMatrixPart <- matrix(data = rep(x = "", times = numberOfRowsHere * numberOfColumnsHere), nrow = numberOfRowsHere, ncol = numberOfColumnsHere)
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
  dataList$dataFrameMS1Header[[1,2]] <<- serializeSampleSelectionAndOrder(dataList$groupSampleDataFrame)
  ms1Matrix     <- rbind(
    dataList$dataFrameMS1Header,
    dataList$dataFrameInfos[precursorSet, ]
  )
  ms1Matrix     <- as.matrix(ms1Matrix)
  
  ###########################################################
  ## export annotations
  
  ## process annotations
  annotations <- dataList$annoArrayOfLists
  for(i in 1:length(annotations))
    if(dataList$annoArrayIsArtifact[[i]])
      annotations[[i]] <- c(annotations[[i]], dataList$annotationValueIgnore)
  
  annotationStrings <- vector(mode = "character", length = length(annotations))
  for(i in 1:length(annotations)){
    if(length(annotations[[i]]) > 0)
      annotationStrings[[i]] <- paste(annotations[[i]], collapse = ", ")
    else
      annotationStrings[[i]] <- ""
  }
  annotationStrings <- annotationStrings[precursorSet]
  
  ## process annotaiotn-color-map
  annoPresentAnnotations <- dataList$annoPresentAnnotationsList[-1]
  annoPresentColors      <- dataList$annoPresentColorsList[-1]
  
  if(length(annoPresentAnnotations) > 0){
    annotationColors <- paste(annoPresentAnnotations, annoPresentColors, sep = "=", collapse = ", ")
  } else {
    annotationColors <- ""
  }
  annotationColors <- paste(dataList$annotationColorsName, "={", annotationColors, "}", sep = "")
  
  ## box
  annotationColumn <- c("", annotationColors, dataList$annotationColumnName, annotationStrings)
  
  ms1Matrix[, dataList$annotationColumnIndex] <- annotationColumn
  
  ################################################################################
  ## assemble
  #dataFrame <- cbind(
  #  ms1Matrix,
  #  ms2Matrix
  #)
  linesMS1MatrixWithHeader <- apply(X = ms1Matrix, MARGIN = 1, FUN = function(x){paste(x, collapse = "\t")})
  lines <- paste(linesMS1MatrixWithHeader, linesFragmentMatrixWithHeader, sep = "\t")
  
  return(lines)
}

writeTable <- function(precursorSet, file){
  show_modal_spinner(spin="scaling-squares", 
                     text="\nMerging project files to csv. This can take several minutes!")
  lines <- createExportMatrix(precursorSet)
  gz1 <- gzfile(description = file, open = "w")
  writeLines(text = lines, con = gz1)
  close(gz1)
  remove_modal_spinner()
}

## individual downloads
output$downloadGlobalMS2filteredPrecursors <- downloadHandler(
  filename = function() {
    createExportMatrixName()
  },
  content = function(file) {
    precursorSet <- filterGlobal$filter
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)
output$downloadHcaFilteredPrecursors <- downloadHandler(
  filename = function() {
    createExportMatrixName()
  },
  content = function(file) {
    precursorSet <- filterHca$filter
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)
output$downloadPcaFilteredPrecursors <- downloadHandler(
  filename = function() {
    createExportMatrixName()
  },
  content = function(file) {
    precursorSet <- filterPca$filter
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)
output$downloadSearchPrecursors <- downloadHandler(
  filename = function() {
    createExportMatrixName()
  },
  content = function(file) {
    precursorSet <- filterPca$filter
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)
output$downloadSelectedPrecursors <- downloadHandler(
  filename = function() {
    createExportMatrixName()
  },
  content = function(file) {
    precursorSet <- selectedPrecursorSet
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)

#Obvserve button for exporting the project
observeEvent(input$prepareAllPrecursors, {
  ExportMatrixName <<- createExportMatrixName()
  precursorSet <- 1:dataList$numberOfPrecursors
  writeTable(precursorSet = precursorSet, file = file.path(tempdir(),ExportMatrixName))
  showModal(modalDialog(title = "Download", footer = NULL, size="s", 
                        fluidRow(column(12, p("Your download is ready."))), 
                        fluidRow(column(3,downloadButton(outputId = "downloadAllpreparedPrecursors", 
                                                         label = "Download project"))),
  ))
})

#Serving the modal with the download button to download the project
output$downloadAllpreparedPrecursors <- downloadHandler(
  filename <- ExportMatrixName, 
  content = function(file){
    file.copy(file.path(tempdir(),ExportMatrixName), file)
    file.remove(file.path(tempdir(),ExportMatrixName))
    removeModal()
  },
  contentType = 'text/csv'
)

## download selected
output$downloadHcaSelectedPrecursors <- downloadHandler(
  filename = function() {
    createExportMatrixName()
  },
  content = function(file) {
    ## get selected precursors
    if(is.null(selectionAnalysisTreeNodeSet)){
      ## all precursors
      precursorSet <- filterHca$filter
    } else {
      precursorSet <- getPrecursorSetFromTreeSelections(clusterDataList = clusterDataList, clusterLabels = selectionAnalysisTreeNodeSet)
    }
    
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)

output$downloadImportParameterSet <- downloadHandler(
  filename = function() {
    createImportParameterSetExportFileName()
  },
  content = function(file) {
    fileLines <- serializeParameterSetFile(dataList$importParameterSet, toolName, toolVersion)
    writeLines(text = fileLines, con = file)
  },
  contentType = 'text/csv'
)

## download images
output$downloadHcaImage <- downloadHandler(
  filename = function() {
    fileType <- input$downloadHcaImageType
    createExportImageName("HCA", fileType)
  },
  content = function(file) {
    fileType <- input$downloadHcaImageType
    plotHCA(file, fileType)
  }#,
  #contentType = 'image/png'
)

plotHCA <- function(file, fileType, plotMS2 = TRUE){
  ## parameters
  widthInInch     <- 10
  heigthInInch    <- ifelse(test = plotMS2, yes = 7.5, no = (5.2-1.5)/5.2 * 7.5)
  resolutionInDPI <- 650
  widthInPixel    <- widthInInch  * resolutionInDPI
  heightInPixel   <- heigthInInch * resolutionInDPI
  
  switch(fileType,
         "png"={
           png(filename = file, 
               width = widthInPixel, 
               height = heightInPixel, 
               res = resolutionInDPI, bg = "white")
         },
         "svg"={
           svg(filename = file)
         },
         "pdf"={
           pdf(file = file, title = "PCA image export from MetFam")
         },
         stop(paste("Unknown file type (", fileType, ")!", sep = ""))
  )
  
  if(plotMS2) {
    graphics::layout(
      mat = matrix(
        data = c(1, 1, 1, 1, 2, 2,3, 4, 5, 5,5, 5), nrow = 6, ncol = 2),
      widths = c(4.58, 1.59), 
      heights = c(0.0052, 0.011, 0.01112, 0.01112,0.001111, 0.0111)
    )
  } else {
    graphics::layout(
      mat = matrix(
        data = c(1, 1, 1, 1, 2, 
                 3, 4, 5, 6, 6), 
        nrow = 5, ncol = 2), 
      widths = c(4, 1), 
      heights = c(0.6, 1.4, 0.6, 1.2, 0.5)
    )
  }
  
  drawDendrogramPlotImpl() #1
  drawHeatmapPlotImpl() #2
  
  drawFragmentDiscriminativityLegendImpl() #4
  drawHeatmapLegendImpl() #3
  
  drawAnnotationLegendForImageHCAimpl() #5

  dev.off()
}

output$downloadPcaImage <- downloadHandler(
  filename = function() {
    fileType <- input$downloadPcaImageType
    createExportImageName("PCA", fileType)
  },
  content = function(file) {
    fileType <- input$downloadPcaImageType
    plotPCA(file, fileType)
  }
)

plotPCA <- function(file, fileType, plotMS2 = TRUE){
  ## parameters
  widthInInch     <- 10
  heigthInInch    <- ifelse(test = plotMS2, yes = 6, no = (5.2-1.5)/5.2 * 6)  
  resolutionInDPI <- 650
  widthInPixel    <- widthInInch  * resolutionInDPI
  heightInPixel   <- heigthInInch * resolutionInDPI
  
  switch(fileType,
         "png"={
           png(filename = file, 
               width = widthInPixel, 
               height = heightInPixel, 
               res = resolutionInDPI, bg = "white")
         },
         "svg"={
           svg(filename = file)
         },
         "pdf"={
           pdf(file = file, title = "PCA image export from MetFam")
         },
         stop(paste("Unknown file type (", fileType, ")!", sep = ""))
  )
  
  if(plotMS2){
    graphics::layout(
      mat = matrix(
        data = c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4,5,5,5,5,6,6,6,6,7,7,7,7),nrow=4,ncol=7),
      widths =  c(0.85, 1.15),
      heights = c(1.02,1.01)
    )
  } else {
    graphics::layout(
      mat = matrix(
        data = c(1, 1, 1, 
                 2, 2, 2, 
                 3, 4, 5), 
        nrow = 3, ncol = 3), 
      widths = c(2, 2, 1), 
      heights = c(0.7, 0.6, 2.4)
    )
  }
  
  drawPcaScoresPlotImpl(downloadLayout=TRUE)
  calcPlotScoresGroupsLegendForImage1(scoresGroups$grouXXXps, scoresGroups$colors, 30)
  drawPcaLoadingsPlotImpl()
  drawAnnotationLegendForImagePCAimpl()
  dev.off()
}

output$downloadDistanceMatrix <- downloadHandler(
  filename = function() {
    createExportDistanceMatrixName(currentDistanceMatrixObj$distanceMeasure)
  },
  content = function(file) {
    write.table(x = currentDistanceMatrixObj$distanceMatrix, file = file, sep = "\t", row.names = TRUE, quote = FALSE, col.names=NA)
  },
  contentType = 'text/csv'
)

## download publication data
output$downloadMsData <- downloadHandler(
  filename = function() {
    return("Metabolite_profile_showcase.txt")
  },
  content = function(file) {
    ## copy data for download
    file.copy(system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily"),
              file)
  },
  contentType = "application/zip"
)

## STN: should this really be a fixed fileName ??
output$downloadMsMsData <- downloadHandler(
  filename = function() {
    return("MSMS_library_showcase.msp")
  },
  content = function(file) {
    ## copy data for download
    file.copy(system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily"),
              file)
  },
  contentType = "application/zip"
)

## STN: should this really be a fixed fileName ??
output$downloadFragmentMatrix <- downloadHandler(
  filename = function() {
    return("Fragment_matrix_showcase.csv")
  },
  content = function(file) {
    ## copy data for download
    file.copy(system.file("Fragment_matrix_showcase.csv", package = "MetFamily"),
      file)
  },
  contentType = "application/zip"
)

output$downloadDocShowcaseProtocol <- downloadHandler(
  filename = function() {
    return("MetFamily_Showcase_protocol.pdf")
  },
  content = function(file) {
    ## copy data for download
    file.copy(system.file("extdata/showcase/MetFamily_Showcase_protocol.pdf", package = "MetFamily"),
      file)
  },
  contentType = "application/pdf"
)

output$downloadDocUserGuide <- downloadHandler(
  filename = function() {
    return("MetFamily_user_guide.pdf")
  },
  content = function(file) {
    ## copy data for download
    file.copy(system.file("extdata/showcase/MetFamily_user_guide.pdf", package = "MetFamily"), 
              file)
  },
  contentType = "application/pdf"
)

output$downloadDocInputSpecification <- downloadHandler(
  filename = function() {
    return("MetFamily_Input_Specification.pdf")
  },
  content = function(file) {
    ## copy data for download
    file.copy(system.file("extdata/showcase/MetFamily_Input_Specification.pdf", package = "MetFamily"), 
      file)
  },
  contentType = "application/pdf"
)

#########################################################################################
## consensus spectrum for metabolite families
output$downloadMetaboliteFamilyConsensusSpectrum <- downloadHandler(
  filename = function() {
    createConsensusSpectrumName(allAnnotationNames[[input$familySelectionTable_rows_selected]])
  },
  content = function(file) {
    annotation <- allAnnotationNames[[input$familySelectionTable_rows_selected]]
    precursorSet <- which(unlist(lapply(X = dataList$annoArrayOfLists, 
                                        FUN = function(y){any(y==annotation)})))
    returnObj <- getSpectrumStatistics(dataList = dataList, precursorSet = precursorSet)
    fragmentMasses = returnObj$fragmentMasses
    fragmentCounts = returnObj$fragmentCounts
    fragmentProportion = fragmentCounts / length(precursorSet)
    consensusSpectrumDf <- data.frame(
      "Mass" = fragmentMasses, 
      "Count" = fragmentCounts, 
      "Frequency" = fragmentProportion
    )
    
    write.table(x = consensusSpectrumDf, file = file, sep = "\t", 
                row.names = FALSE, quote = FALSE)
  },
  contentType = 'text/csv'
)

output$downloadMetaboliteFamilyFilteredPrecursors <- downloadHandler(
  filename = function() {
    createMetaboliteFamilyProjectFileName(allAnnotationNames[[input$familySelectionTable_rows_selected]])
  },
  content = function(file) {
    annotation <- allAnnotationNames[[input$familySelectionTable_rows_selected]]
    precursorSet <- which(unlist(lapply(X = dataList$annoArrayOfLists, 
                                        FUN = function(y){any(y==annotation)})))
    writeTable(precursorSet = precursorSet, file = file)
  },
  contentType = 'text/csv'
)

#########################################################################################
## classifier annotation result export
createAnnotationResultTableAll <- function(){
  precursorIndeces <- sort(as.integer(unique(unlist(lapply(X = classToSpectra_class, FUN = names)))))
  
  rows <- list()
  for(precursorIndex in precursorIndeces){
    precursorLabel     <- dataList$precursorLabels[[precursorIndex]]
    mz                 <- dataList$dataFrameInfos[[precursorIndex, "m/z"]]
    rt                 <- dataList$dataFrameInfos[[precursorIndex, "RT"]]
    metaboliteName     <- dataList$dataFrameInfos[[precursorIndex, "Metabolite name"]]
    
    presentAnnotations <- unlist(dataList$annoArrayOfLists[[precursorIndex]])
    if(length(presentAnnotations) == 0){
      presentAnnotations <- ""
    } else {
      presentAnnotations <- sort(unlist(presentAnnotations))
      presentAnnotations <- paste(presentAnnotations, collapse = "; ")
    }
    
    for(classIdx in seq_along(classToSpectra_class)){
      if(!(precursorIndex %in% as.integer(names(classToSpectra_class[[classIdx]]))))
        next
      
      class <- names(classToSpectra_class)[[classIdx]]
      pValue <- format(unname( classToSpectra_class[[classIdx]][[which( precursorIndex == as.integer(names(classToSpectra_class[[classIdx]])) )]] ), digits=4)
      rows[[length(rows)+1]] <- c(
        "Index" = precursorIndex, 
        "Label" = precursorLabel, 
        "m/z"   = mz, 
        "RT"    = rt, 
        "Metabolite name" = metaboliteName, 
        "Annotation (present)" = presentAnnotations, 
        "Annotation (putative)" = class, 
        "pValue" = pValue
      )
    }
  }
  
  head <- c(
    "Index", 
    "Label", 
    "m/z", 
    "RT", 
    "Metabolite name", 
    "Annotation (present)", 
    "Annotation (putative)", 
    "pValue"
  )
  outputDf <- as.data.frame(t(matrix(data = unlist(rows), nrow = length(head))))
  colnames(outputDf) <- head
  return(outputDf)
}

createAnnotationResultTableForClass <- function(){
  selectedRowIdx <- input$annotationResultTableClass_rows_selected
  
  class          <- names(classToSpectra_class)[[selectedRowIdx]]
  classToSpectra <-       classToSpectra_class [[selectedRowIdx]]
  
  precursorIndeces   <- as.integer(names(classToSpectra))
  pValues            <- format(unname(classToSpectra), digits=4)
  precursorLabels    <- dataList$precursorLabels[precursorIndeces]
  mzs                <- dataList$dataFrameInfos[precursorIndeces, "m/z"]
  rts                <- dataList$dataFrameInfos[precursorIndeces, "RT"]
  metaboliteNames    <- dataList$dataFrameInfos[precursorIndeces, "Metabolite name"]
  presentAnnotations <- unlist(lapply(X = dataList$annoArrayOfLists[precursorIndeces], 
                                      FUN = function(x){
                                        if(length(x) == 0){
                                          return("")
                                        } else {
                                          x <- sort(unlist(x))
                                          s <- paste(x, collapse = "; ")
                                          return(s)
                                        }
                                      }))
  
  outputDf <- data.frame(
    "Index" = precursorIndeces, 
    "Label" = precursorLabels, 
    "m/z"   = mzs, 
    "RT"    = rts, 
    "Metabolite name" = metaboliteNames, 
    "Annotation (present)" = presentAnnotations, 
    "Annotation (putative)" = rep(x = class, times = length(precursorIndeces)), 
    "pValue" = pValues
  )
  
  return(outputDf)
}

output$downloadAllAnnotationResults <- downloadHandler(
  filename = function() {
    createClassifierAnnotationName( "All" )
  },
  content = function(file) {
    print(file)
    outputDf <- createAnnotationResultTableAll()
    write.table(x = outputDf, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
  },
  contentType = 'text/csv'
)

output$downloadMetaboliteFamilyAnnotationResults <- downloadHandler(
  filename = function() {
    createClassifierAnnotationName( names(classToSpectra_class)[[input$annotationResultTableClass_rows_selected]] )
  },
  content = function(file) {
    outputDf <- createAnnotationResultTableForClass()
    write.table(x = outputDf, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
  },
  contentType = 'text/csv'
)

#########################################################################################
## report
output$downloadReport2 <- downloadHandler(
  filename = function() {
    return("MetFamilyReport.pdf")
    #return("MetFamilyReport.html")
  },
  content = function(file) {
    ##########################################################################
    ## files
    
    ## source file and tmp file
    tempReportFile <- file.path(tempdir(), "MetFamilyReport.Rmd")
    reportSourceFile <- "report/Report.Rmd"
    
    ## copy template to tmp dir for reasons of file permissions
    file.copy(from = reportSourceFile, to = tempReportFile, overwrite = TRUE)
    
    ##########################################################################
    ## HCA analyses
    drawHCA = state$showHCAplotPanel
    drawPCA = state$showPCAplotPanel
    
    ##########################################################################
    ## HCA analysis
    if(drawHCA){
      imageFileHCA <- file.path(tempdir(), "HcaTmpFile.png")
      plotHCA(file = imageFileHCA, fileType = "png")
    } else {
      imageFileHCA <- ""
    }
    
    ##########################################################################
    ## PCA analysis
    if(drawPCA){
      imageFilePCA <- file.path(tempdir(), "PcaTmpFile.png")
      plotPCA(file = imageFilePCA, fileType = "png")
    } else {
      imageFilePCA <- ""
    }
    
    # Set up parameters to pass to Rmd document
    params <- list(
      creationTime = date(),
      importParameterSet = dataList$importParameterSet,
      drawHCA = drawHCA,
      drawPCA = drawPCA,
      imageFileHCA = imageFileHCA,
      imageFilePCA = imageFilePCA,
      clusterDataList = clusterDataList,
      pcaDataList = pcaDataList
    )
    
    # Knit the document and eval it in a child of the global environment (this isolates the code in the document from the code in this app)
    rmarkdown::render(input = tempReportFile, output_file = file,
                      output_format = "pdf_document",
                      params = params,
                      envir = new.env(parent = globalenv()),
                      quiet = FALSE
    )
  },
  contentType = "application/pdf"
)
