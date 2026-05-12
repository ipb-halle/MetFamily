
#' Read mzmine Output File into a QFeatures Object
#'
#' Read a metabolite profile output file (.csv) from mzmine and 
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
readMZmine <- function(file, version){
  
  table <- readr::read_csv(
    file, col_types = readr::cols(
      .default = readr::col_character())
    ) %>% as.data.frame
    
    # expected names
    stopifnot(
      identical(colnames(table[,1:11]),
      c("id", "mz", "mz_range:min", "mz_range:max", "rt", 
      "rt_range:min", "rt_range:max", "area", "height", 
      "intensity_range:min", "intensity_range:max"))
    )
    
    # tmp <- names(table[, grepl("datafile:", names(table))]) %>%
    #    stringr::str_remove(".*:") %>%
    #    unique
    
    # match MS-Dial names "narrow" format
    table <- table %>% 
      dplyr::rename(
        "Alignment ID" = "id",
        "Average Rt(min)" = "rt",
        "Average Mz" = "mz",
        "Metabolite name" = "preferred_annotation:compound_name",
        "Adduct ion name" = "preferred_annotation:adduct"
      )
    
    heights <- table[, grepl(":height$", names(table))]
    names(heights) <- names(heights) %>% 
      stringr::str_remove(":height") %>%
      stringr::str_remove("datafile:") %>% 
      make.names
    
    counts <- as.matrix(heights)
    rownames(counts) <- table$id
    
    sampleNames <- colnames(heights)
    
    # colData
    # TODO how to determine sample classes?
    colData <- data.frame(
      Class = sampleNames,
      Type = "Sample",
      row.names = sampleNames
    )
    
    # Extract rowData
    # rowData_cols <- names(table)[!grepl("datafile:",names(table))]
    rowData <- table[1:13]
    rownames(rowData) <- table$id
    
    
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
