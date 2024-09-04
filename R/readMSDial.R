#' Read MS-DIAL Output File into a QFeatures Object
#'
#' This function reads the output file from MS-DIAL and 
#' converts it into a QFeatures object.
#'
#' @param file A string with the path to the MS-DIAL output file.
#' @param version A character string specifying the version of MS-DIAL used to generate the file.
#'   This parameter is currently not used.
#'
#' @return A QFeatures object containing:
#'   \itemize{
#'     \item An assay named "exampleAssay" with the metabolite counts.
#'     \item Row data (feature metadata) extracted from the input file.
#'     \item Column data (sample metadata) extracted from the input file.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have an MS-DIAL output file named "Metabolite_profile_showcase.txt" in a "data" directory:
#' qf <- readMSDial("data/Metabolite_profile_showcase.txt")
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
#' @importFrom QFeatures QFeatures
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @seealso 
#' \code{\link[QFeatures]{QFeatures}} for more information on the QFeatures class.
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for details on the underlying data structure.
#'
#' @note 
#'
#'
#' @references
#' 
readMSDial <- function(file, version){
    table <- read.table(file, fill = TRUE, sep = "\t",
                        quote = "", header = FALSE)
    
    # Identify the starting row and column of the data
    startRow <- which(table[, 1] != "")[1]
    startCol <- which(table[1, ] != "")[1]
    ##TODO: version dependent error message if startRow or startCol are not as expected.
    
    # Split the table in parts
    colDataRaw <- table[1:startRow, startCol:ncol(table)]
    rowDataRaw <- table[startRow:nrow(table), 1:(startCol)]
    countsRaw <- table[startRow:nrow(table), startCol:ncol(table)]
    
    # Extract ids and counts data
    ids <- rowDataRaw[-1, 1]
    counts <- as.matrix(countsRaw[-1, -1])
    counts <- matrix(as.numeric(counts), nrow = nrow(counts), ncol = ncol(counts))
    colnames(counts) <- as.character(countsRaw[1, -1])
    rownames(counts) <- ids
    
    # Ensure row names of colData match counts column names
    colData <- data.frame(t(colDataRaw[-nrow(colDataRaw), -1]))
    rownames(colData) <- as.character(colDataRaw[nrow(colDataRaw), -1])
    colnames(colData) <- as.character(colDataRaw[-nrow(colDataRaw), 1])

    # Ensure row names of rowData match counts row names
    rowData <- data.frame(rowDataRaw[-1, ], row.names = ids)
    colnames(rowData) <- as.character(rowDataRaw[1,])

    # Create SummarizedExperiment object
   
    sumExp <- SummarizedExperiment(assays = list(counts = counts),
                                    rowData = rowData,
                                    colData = colData)
    ##TODO: Metadata with data source and version 
    
    # Create QFeatures object
    qf <- QFeatures(list(exampleAssay = sumExp), colData = colData(sumExp))
    qf
    ##TODO: name
  }