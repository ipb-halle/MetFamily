#' Read Metaboscape Output File into a QFeatures Object
#'
#' This function reads a metabolite profile output file (.csv) from Metaboscape and 
#' converts it into a QFeatures object.
#' 
#' At the moment, sample classes are not considered.
#'
#' @param file A character string specifying the path to the Metaboscape output file (csv format).
#' @param version A character string specifying the version of Metaboscape used to generate the file.
#'   This parameter is currently not used.
#'
#' @return A QFeatures object containing:
#'   \itemize{
#'     \item An assay named "exampleAssay" with the metabolite counts.
#'     \item Row data (feature metadata) extracted from the input file.
#'     \item Column data (sample metadata) extracted from the sample names, including injection order and sample name.
#'   }
#'   
#' @importFrom QFeatures QFeatures
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Metaboscape output file named "data.csv":
#' qf <- readMetaboscape("data.csv") #TODO: System file 
#'
#' 
#' # Examine the structure of the resulting QFeatures object
#' qf
#' 
#' # Access the assay data
#' assay(qf[["exampleAssay"]])
#' 
#' # Access the row data (feature metadata)
#' rowData(qf[["exampleAssay"]])
#' 
#' # Access the column data (sample metadata)
#' colData(qf)
#' }
#'
#'
#' @seealso 
#' \code{\link[QFeatures]{QFeatures}} for more information on the QFeatures class.
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for details on the underlying data structure.
#' 
#' @references
#' #TODO: Bruker metaboscape site 
#' #TODO: Ordering ?
#' 
readMetaboscape <- function(file, version){
  
  # file <- file.path("../file-formats/Metaboscape-export-version 2025b",
  #                   "UTH-2025-07-29-UTH003_002-conyza-test-samples.csv")
  # file.exists(fileSpectra)
  
  table <- readr::read_csv(file, col_types = readr::cols(
    .default = readr::col_character())) %>% as.data.frame
  
  colMeanInt <- stringr::str_detect(colnames(table), "_MeanIntensity")
  startOfSamples <- rev(which(colMeanInt))[1] + 1
  colIdsSamples <- startOfSamples:length(table)
  
  # expected names
  stopifnot(
    identical(colnames(table[,1:11]),
          c("FEATURE_ID", "RT", "PEPMASS", "CCS", "SIGMA_SCORE",
            "NAME_METABOSCAPE", "MOLECULAR_FORMULA", "ADDUCT", 
            "KEGG", "CAS", "MaxIntensity"))
  )
  
  # match MS-Dial names "narrow" format
  table <- table %>% 
    dplyr::rename("Alignment ID" = "FEATURE_ID",
           "Average Rt(min)" = "RT",
           "Average Mz" = "PEPMASS",
           "Metabolite name" = "NAME_METABOSCAPE",
           "Adduct ion name" = "ADDUCT")
  
  # rt in minutes
  table <- table %>% 
    dplyr::mutate(
      `Average Rt(min)` = as.character(as.numeric(`Average Rt(min)`) / 60),
      # needed to match to MGF spectra
      `Average Mz` = as.character(as.numeric(`Average Mz`) - 1.00727))
  
  # colData
  # TODO how to determine sample classes?
  sampleNames <- colnames(table)[colIdsSamples]
  colData <- data.frame(
    Class = sampleNames,
    Type = "Sample",
    row.names = sampleNames
  )
  
  # Extract ids and counts data
  ids <- table %>% dplyr::pull(1)
  countsRaw     <- table[,colIdsSamples]
  countsNumeric <- apply(countsRaw, 2, as.numeric)
  counts <- as.matrix(countsNumeric)
  colnames(counts) <- sampleNames
  rownames(counts) <- ids
  
  # Extract rowData
  rowData <- table[1:11]
  # which(colMeanInt)[1] - 1
  rownames(rowData) <- ids 
  
  # Create SummarizedExperiment object
  sumExp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
  
  # Create QFeatures object
  qf <- QFeatures::QFeatures(
    list(exampleAssay = sumExp), 
    colData = SummarizedExperiment::colData(sumExp)
  )
  
  qf
}