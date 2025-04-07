#' Read MS-DIAL Output File into a QFeatures Object
#'
#' This function reads the output file from MS-DIAL and converts it into a
#' QFeatures object.
#'
#' @param file A string with the path to the MS-DIAL output file.
#' @param version A character string specifying the version of MS-DIAL used to
#'   generate the file. This parameter is currently not used.
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
#' # Assuming you have an MS-DIAL output file named
#'  "Metabolite_profile_showcase.txt" in a "data" directory:
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
#' @seealso \code{\link[QFeatures]{QFeatures}} for more information on the
#' QFeatures class. \code{\link[SummarizedExperiment]{SummarizedExperiment}} for
#' details on the underlying data structure.
#' 
readMSDial <- function(file, version){
    table <- read.table(file, fill = TRUE, sep = "\t",
                        quote = "", header = FALSE, colClasses = "character", comment.char = "")
    
    # table <- readr::read_tsv(file, col_names = F, col_types = cols(.default = col_character()))
    
    # Identify the starting row and column of the data
    startRow <- which(table[, 1] != "")[1]
    startCol <- which(table[1, ] != "")[1]
    lastCol <- which(!is.na(table[1,])) %>% .[length(.)]
    
    rowDataNames <- table[startRow, 1:startCol] %>% as.character
    colDataNames <- table[1:(startRow - 1), startCol] %>% as.character
    
    # these are the expected values, based on an 'older' version of MS-Dial
    rowDataDefaultNames2020 <- c("Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", 
                                 "Adduct ion name", "Fill %", "MS/MS included", "INCHIKEY", "SMILES", 
                                 "LINK", "Dot product", "Reverse dot product", "Fragment presence %", 
                                 "Spectrum reference file name")
    
    colDataDefaultNames2020 <- c("Class", "Type", "Injection order")
    
    if (startCol == 14) {
      # older version
      msdVersion <- "narrow"
      
      stopifnot(identical(rowDataNames, rowDataDefaultNames2020))
      stopifnot(identical(colDataNames, colDataDefaultNames2020))
      stopifnot(lastCol == ncol(table))
      
      # continue

    } else if (startCol == 35) {
      # 'newer' version
      msdVersion <- "wide"
      
      rowDataDefaultNames2024 <- c("Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", 
                                   "Adduct type", "Post curation result", "Fill %", "MS/MS assigned", 
                                   "Reference RT", "Reference m/z", "Formula", "Ontology", "INCHIKEY", 
                                   "SMILES", "Annotation tag (VS1.0)", "RT matched", "m/z matched", 
                                   "MS/MS matched", "Comment", "Manually modified for quantification", 
                                   "Manually modified for annotation", "Isotope tracking parent ID", 
                                   "Isotope tracking weight number", "RT similarity", "m/z similarity", 
                                   "Simple dot product", "Weighted dot product", "Reverse dot product", 
                                   "Matched peaks count", "Matched peaks percentage", "Total score", 
                                   "S/N average", "Spectrum reference file name", "MS1 isotopic spectrum", 
                                   "MS/MS spectrum")

      stopifnot(identical(rowDataNames, rowDataDefaultNames2024))

            
      # remove average cols
      if (ncol(table) > lastCol) {
        table <- table[,-((lastCol+1):ncol(table))]
      }
      
            
      # colData
      colDataBaseNames2024 <- c("Class", "File type", "Injection order", "Batch ID")
      
      if (length(colDataNames) > 4) {
        # parameter rows present
        nbParam <- length(colDataNames)-4
        paramRows <- paste0("Parameter", seq_len(nb_param))
        
        colDataDefaultNames2024 <- c("Class", paramRows, "File type", "Injection order", "Batch ID")
        stopifnot(identical(colDataNames, colDataDefaultNames2024))
        
        # remove parameter lines
        table <- table[-(2:(1+nbParam)),]
        
        startRow <- which(table[, 1] != "")[1]
        
      } else {
        stopifnot(identical(colDataNames, colDataBaseNames2024))
      }
      
    } else {
      stop("MS-Dial table format not supported.")
    }

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

    # hack back to match "narrow" format
    if (msdVersion == "wide") {
      
      # remove, rename
      colData <- colData %>% dplyr::select(-`Batch ID`, Type = `File type`)
      stopifnot(identical(names(colData), colDataDefaultNames2020))
      
      # create, rename, only keep narrow format
      rowData <- rowData %>% 
        dplyr::mutate(LINK = "No record",
               "Fragment presence %" = "-1") %>% 
        dplyr::rename("Adduct ion name" = "Adduct type",
                      "MS/MS included" = "MS/MS assigned",
                      "Dot product" = "Simple dot product") %>% 
        dplyr::select(all_of(rowDataDefaultNames2020))
      
      stopifnot("MS-Dial rowData names are not as expected" = identical(names(rowData), rowDataDefaultNames2020))
      
    }
    
    # Create SummarizedExperiment object
    sumExp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
                                   rowData = rowData,
                                   colData = colData)
    
    ##TODO: Metadata with data source and version 

    # Create QFeatures object
    qf <- QFeatures::QFeatures(list(exampleAssay = sumExp), colData = SummarizedExperiment::colData(sumExp))
    qf
    
    ##TODO: name
}
