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
#' # "Metabolite_profile_showcase.txt" in a "data" directory:
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


#' Convert mzTab-M file to QFeatures object (full-featured with metadata enrichment)
#'
#' Parses an mzTab-M file and returns a QFeatures object containing both the
#' Small Molecule Feature (SMF) and Small Molecule Summary (SML) sections
#' with comprehensive metadata enrichment and assay linking.
#'
#' @param file character(1) Path to the mzTab-M file.
#' @param load_sml logical(1) Include SML section and link it to SMF. Default TRUE.
#' @param enrich_coldata logical(1) If TRUE, enrich colData with sample/ms_run/study_variable 
#'   metadata and harmonize across assays. Default FALSE.
#'
#' @return A QFeatures object with assays:
#' \itemize{
#'   \item \code{"SMF"}: Small Molecule Feature assay with quantification data
#'   \item \code{"SML"}: Small Molecule Summary assay with quantification data (if \code{load_sml = TRUE})
#' }
#' Each assay contains:
#' \itemize{
#'   \item \code{assays$intensity}: A double matrix with quantification data
#'   \item \code{rowData}: Non-quantitative columns from respective sections
#'   \item \code{colData}: Enriched metadata when \code{enrich_coldata = TRUE}
#' }
#' The SML assay is linked to the SMF assay via \code{SMF_ID_REFS} references.
#'
#' @details
#' This function provides mzTab-M to QFeatures conversion:
#' \itemize{
#'   \item Processes both SMF and SML sections (when \code{load_sml = TRUE})
#'   \item Requires assay-level quantification columns named \code{abundance_assay[n]}
#'   \item Derives column metadata from MTD \code{assay[...]} keys
#'   \item Creates assay links between SML and SMF based on \code{SMF_ID_REFS}
#'   \item Supports both 1:1 and many-to-many mappings via Hits objects
#'   \item Orders columns by numeric assay index for consistent mapping
#'   \item Uniquifies row identifiers from \code{SMF_ID}/\code{SML_ID} columns
#' }
#'
#' @section Metadata enrichment:
#' When \code{enrich_coldata = TRUE}, the function enriches colData with:
#' \itemize{
#'   \item \code{assay_name}, \code{assay_description}: Direct assay attributes
#'   \item \code{sample_id}, \code{sample_label}, \code{sample_description}: Sample metadata
#'   \item \code{sample_species}, \code{sample_tissue}: Biological context
#'   \item \code{ms_run_id}, \code{ms_run_location}, \code{ms_run_format}: MS run details
#'   \item \code{study_variable}: Study variable membership
#' }
#' The function harmonizes colData across assays to avoid conflicts in QFeatures.
#'
#' @section Linking mechanism:
#' The function creates links between SML and SMF rows using the \code{SMF_ID_REFS}
#' column in the SML section:
#' \itemize{
#'   \item For 1:1 mappings: Uses \code{addAssayLink()} with \code{varFrom}/\code{varTo}
#'   \item For many-to-many mappings: Creates \code{Hits} object for complex relationships
#'   \item Unknown references are dropped to prevent silent mis-linking
#' }
#'
#' @section Requirements:
#' \itemize{
#'   \item mzTab-M file must contain \code{abundance_assay[...]} columns
#'   \item MTD section must contain \code{assay[...]} metadata entries
#'   \item SMF section must contain \code{SMF_ID} column
#'   \item SML section must contain \code{SML_ID} column (if \code{load_sml = TRUE})
#'   \item SML section should contain \code{SMF_ID_REFS} for linking (optional)
#' }
#'
#' @section Edge cases:
#' \itemize{
#'   \item Duplicate IDs will be made unique automatically
#'   \item Conflicting colData columns are dropped with warning when harmonizing
#' }
#'
#' @examples
#' \dontrun{
#' # Basic conversion (SMF + SML, no metadata enrichment)
#' qf <- readMzTabM("path/to/file.mzTab")
#' 
#' # Full-featured conversion with metadata enrichment
#' qf_enriched <- readMzTabM("path/to/file.mzTab", 
#'                                      load_sml = TRUE, 
#'                                      enrich_coldata = TRUE)
#' 
#' # SMF only with metadata enrichment
#' qf_smf_only <- readMzTabM("path/to/file.mzTab", 
#'                                      load_sml = FALSE, 
#'                                      enrich_coldata = TRUE)
#' 
#' # Access enriched metadata
#' colData(qf_enriched)
#' colData(qf_enriched[["SMF"]])
#' }
#'
#' @export
readMzTabM <- function(file, load_sml = TRUE, enrich_coldata = FALSE) {

  # Read mzTab-M
  mztab_df <- rmzTabM::readMzTab(file)

  smf <- rmzTabM::extractSmallMoleculeFeatures(mztab_df)

  if (isTRUE(load_sml)) {
    sml <- rmzTabM::extractSmallMoleculeSummary(mztab_df)
  } else {
    sml <- NULL
  }

  mtd <- rmzTabM::extractMetadata(mztab_df)
  # Helpers
  normalize_ref_tokens <- function(x) {
    if (!nzchar(x)) {
      return(character(0))
    }
    trimws(
      unlist(
        strsplit(gsub(",", "|", x, fixed = TRUE), "|", fixed = TRUE),
        use.names = FALSE
      )
    )
  }
  get_assay_quant_cols <- function(tab) {
    nms <- names(tab)
    keep <- startsWith(nms, "abundance_assay[") & endsWith(nms, "]")
    nms[keep]
  }

  # Check assay-level quantification columns
  smf_assay_quant_cols <- get_assay_quant_cols(smf)
  if (length(smf_assay_quant_cols) == 0) {
    stop("missing abundance_assay[...] columns in mzTab-M file.")
  }
  
  # Derive colData from MTD assay ids
  keys <- as.character(mtd[["V2"]])
  vals <- as.character(mtd[["V3"]])
  keep <- startsWith(keys, "assay[")
  if (!any(keep)) stop("No assay metadata in MTD")
  k <- keys[keep]
  rb <- regexpr("]", k, fixed = TRUE)
  ids <- substr(k, 1, rb)
  coldata_table <- data.frame(id = unique(ids), stringsAsFactors = FALSE)
  colData <- S4Vectors::DataFrame(coldata_table)
  rownames(colData) <- make.names(colData$id, unique = TRUE)

  if (isTRUE(enrich_coldata)) {
    # Enrich colData with MTD metadata linked via assay/sample/ms_run/study_variable
    cv_label <- function(x) {
      ifelse(
        startsWith(x, "[") & endsWith(x, "]"),
        {
          # Strip surrounding brackets then split by comma and use the 3rd field (term label)
          stripped <- sub("^\\[", "", sub("\\]$", "", x))
          parts <- strsplit(stripped, ",", fixed = TRUE)
          vapply(parts, function(p) {
            p <- trimws(p)
            if (length(p) >= 3) p[3] else paste(p, collapse = ", ")
          }, character(1))
        },
        x
      )
    }

    # Helper to add a column only if any non-empty value exists
    add_if_any <- function(df, name, values) {
      if (length(values) && any(!is.na(values) & nzchar(values))) df[[name]] <- values
      df
    }

    assay_ids <- as.character(colData$id)

    # Direct assay attributes
    assay_name <- vals[match(paste0(assay_ids, "-name"), keys)]
    assay_description <- vals[match(paste0(assay_ids, "-description"), keys)]
    assay_sample_ref <- vals[match(paste0(assay_ids, "-sample_ref"), keys)]
    assay_msrun_ref <- vals[match(paste0(assay_ids, "-ms_run_ref"), keys)]

    # Resolve sample fields
    sample_id <- ifelse(nzchar(assay_sample_ref), assay_sample_ref, NA_character_)
    sample_label <- vals[match(sample_id, keys)]
    sample_description <- vals[match(paste0(sample_id, "-description"), keys)]
    sample_species <- cv_label(vals[match(paste0(sample_id, "-species"), keys)])
    sample_organism <- cv_label(vals[match(paste0(sample_id, "-organism"), keys)])
    if (any(!nzchar(sample_species) & nzchar(sample_organism), na.rm = TRUE)) {
      idx <- which(!nzchar(sample_species) & nzchar(sample_organism))
      sample_species[idx] <- sample_organism[idx]
    }
    sample_tissue <- cv_label(vals[match(paste0(sample_id, "-tissue"), keys)])

    # Resolve ms_run fields
    ms_run_id <- ifelse(nzchar(assay_msrun_ref), assay_msrun_ref, NA_character_)
    ms_run_location <- vals[match(paste0(ms_run_id, "-location"), keys)]
    ms_run_format <- cv_label(vals[match(paste0(ms_run_id, "-format"), keys)])

    # Study variable membership: collect labels for any SV referencing an assay
    sv_keys <- keys[startsWith(keys, "study_variable[")]
    sv_vals <- vals[startsWith(keys, "study_variable[")]
    sv_is_refs <- endsWith(sv_keys, "-assay_refs")
    sv_ref_keys <- sv_keys[sv_is_refs]
    sv_ref_vals <- sv_vals[sv_is_refs]
    sv_ids <- sub("-assay_refs$", "", sv_ref_keys)
    # Map SV id -> label (base key without suffix)
    sv_labels <- vals[match(sv_ids, keys)]
    # Build membership map: assay -> concatenated SV labels
    sv_map <- setNames(vector("list", length(sv_ids)), sv_ids)
    for (i in seq_along(sv_ids)) {
      assays_in_sv <- normalize_ref_tokens(ifelse(is.na(sv_ref_vals[i]), "", sv_ref_vals[i]))
      sv_map[[i]] <- assays_in_sv
    }
    study_variable <- vapply(assay_ids, function(aid) {
      hits <- vapply(seq_along(sv_ids), function(i) if (aid %in% sv_map[[i]]) sv_labels[i] else NA_character_, character(1))
      lbls <- hits[!is.na(hits) & nzchar(hits)]
      if (length(lbls)) paste(unique(lbls), collapse = "; ") else NA_character_
    }, character(1))

    # Assign enriched columns where informative
    colData <- add_if_any(colData, "assay_name", assay_name)
    colData <- add_if_any(colData, "assay_description", assay_description)
    colData <- add_if_any(colData, "sample_id", sample_id)
    colData <- add_if_any(colData, "sample_label", sample_label)
    colData <- add_if_any(colData, "sample_description", sample_description)
    colData <- add_if_any(colData, "sample_species", sample_species)
    colData <- add_if_any(colData, "sample_tissue", sample_tissue)
    colData <- add_if_any(colData, "ms_run_id", ms_run_id)
    colData <- add_if_any(colData, "ms_run_location", ms_run_location)
    colData <- add_if_any(colData, "ms_run_format", ms_run_format)
    colData <- add_if_any(colData, "study_variable", study_variable)
  }
  # Build SummarizedExperiment for a section
  build_summarized_experiment <- function(tab, id_column_name, colData) {
    quant_column_names <- get_assay_quant_cols(tab)
    if (!length(quant_column_names)) stop("No quant columns in section")
    # Order columns by numeric index to ensure stable mapping
    extracted <- sub("^[^\\[]*\\[(\\d+)\\].*$", "\\1", quant_column_names)
    num_ids <- as.integer(extracted)
    ids_order <- order(num_ids)
    quant_column_names <- quant_column_names[ids_order]
    num_ids <- num_ids[ids_order]
    entity_ids <- paste0("assay[", num_ids, "]")
    map_idx <- match(entity_ids, colData$id)
    if (anyNA(map_idx)) stop("MTD ids missing for quant columns")
    # Build double matrix with colnames from colData
    tab[quant_column_names] <- lapply(
      tab[quant_column_names],
      function(x) if (is.double(x)) x else suppressWarnings(as.double(x))
    )
    quant_matrix <- as.matrix(tab[quant_column_names])
    colnames(quant_matrix) <- rownames(colData)[map_idx]
    if (!(id_column_name %in% names(tab))) {
      stop(sprintf("Missing id column %s", id_column_name))
    }
    unique_row_ids <- make.unique(as.character(tab[[id_column_name]]))
    rownames(quant_matrix) <- unique_row_ids
    # Keep non-quant columns as rowData
    row_data <- S4Vectors::DataFrame(
      tab[setdiff(names(tab), quant_column_names)]
    )
    SummarizedExperiment::SummarizedExperiment(
      assays = list(intensity = quant_matrix),
      rowData = row_data,
      colData = colData
    )
  }
  
  se_small_molecule_features <- build_summarized_experiment(
    smf,
    id_column_name = "SMF_ID",
    colData = colData
  )

  if (isTRUE(load_sml)) {
    se_small_molecule_summary <- build_summarized_experiment(
      sml,
      id_column_name = "SML_ID",
      colData = colData
    )
  }

  if (isTRUE(enrich_coldata) && isTRUE(load_sml)) {
    # Harmonize colData across assays to avoid conflicts in QFeatures
    harmonize_coldata <- function(se_list) {
      if (!length(se_list)) return(se_list)
      # Common sample names present in all assays
      common_samples <- Reduce(intersect, lapply(se_list, colnames))
      if (!length(common_samples)) return(se_list)
      # Reorder colData to common order
      se_list <- lapply(se_list, function(se) {
        cd <- SummarizedExperiment::colData(se)
        cd <- cd[common_samples, , drop = FALSE]
        SummarizedExperiment::colData(se) <- cd
        colnames(se) <- common_samples
        se
      })
      # Columns present in all assays' colData
      all_cols <- Reduce(intersect, lapply(se_list, function(se) names(SummarizedExperiment::colData(se))))
      if (!length(all_cols)) return(se_list)
      # Identify columns with identical values across assays
      keep_cols <- character(0)
      for (nm in all_cols) {
        vals0 <- as.character(SummarizedExperiment::colData(se_list[[1]])[[nm]])
        same <- TRUE
        for (j in 2:length(se_list)) {
          valsj <- as.character(SummarizedExperiment::colData(se_list[[j]])[[nm]])
          if (!identical(vals0, valsj)) { same <- FALSE; break }
        }
        if (isTRUE(same)) keep_cols <- c(keep_cols, nm)
      }
      dropped <- setdiff(all_cols, keep_cols)
      if (length(dropped)) {
        message(sprintf("Dropping conflicting colData columns: %s", paste(dropped, collapse = ", ")))
      }
      # Apply keep set to each assay
      se_list <- lapply(se_list, function(se) {
        cd <- SummarizedExperiment::colData(se)
        SummarizedExperiment::colData(se) <- cd[, keep_cols, drop = FALSE]
        se
      })
      se_list
    }

    se_list <- list(SMF = se_small_molecule_features, SML = se_small_molecule_summary)
    se_list <- harmonize_coldata(se_list)
    se_small_molecule_features <- se_list[["SMF"]]
    se_small_molecule_summary <- se_list[["SML"]]
  }
  
  # Construct QFeatures
  qfeatures <- if (isTRUE(load_sml)) {
    QFeatures::QFeatures(
      list(
        SMF = se_small_molecule_features,
        SML = se_small_molecule_summary
      )
    )
  } else {
    QFeatures::QFeatures(list(SMF = se_small_molecule_features))
  }

  # Normalize: strip names() on 'study_variable' in shared and assay-level colData to avoid attribute mismatches
  if (isTRUE(enrich_coldata)) {
    cd_shared <- SummarizedExperiment::colData(qfeatures)
    if (length(cd_shared) && ("study_variable" %in% names(cd_shared))) {
      vec <- cd_shared[["study_variable"]]
      names(vec) <- NULL
      cd_shared[["study_variable"]] <- vec
      SummarizedExperiment::colData(qfeatures) <- cd_shared
    }
    assay_names <- names(qfeatures)
    for (an in assay_names) {
      cd_a <- SummarizedExperiment::colData(qfeatures[[an]])
      if (length(cd_a) && ("study_variable" %in% names(cd_a))) {
        v <- cd_a[["study_variable"]]
        names(v) <- NULL
        cd_a[["study_variable"]] <- v
        SummarizedExperiment::colData(qfeatures[[an]]) <- cd_a
      }
    }
  }

  if (isTRUE(load_sml)) {
    # Link SML to SMF via SMF_ID_REFS
    smf_rd <- SummarizedExperiment::rowData(qfeatures[["SMF"]])
    sml_rd <- SummarizedExperiment::rowData(qfeatures[["SML"]])
    col_smf_id <- if ("SMF_ID" %in% names(smf_rd)) "SMF_ID" else stop("SMF_ID not found")
    col_refs <- if ("SMF_ID_REFS" %in% names(sml_rd)) "SMF_ID_REFS" else NA_character_
    # Normalize refs
    refs_chr <- if (!is.na(col_refs)) {
      trimws(as.character(sml_rd[[col_refs]]))
    } else {
      rep("", nrow(sml_rd))
    }
    refs_list <- lapply(refs_chr, normalize_ref_tokens)
    sml_rd$SMF_ID_REFS <- IRanges::CharacterList(refs_list)
  # Re-apply names on 'study_variable' right before assay replacement to ensure identical attributes
  if (isTRUE(enrich_coldata)) {
    cd_shared_tmp <- SummarizedExperiment::colData(qfeatures)
    if (length(cd_shared_tmp) && ("study_variable" %in% names(cd_shared_tmp))) {
      sv_shared <- cd_shared_tmp[["study_variable"]]
      names(sv_shared) <- rownames(cd_shared_tmp)
      cd_shared_tmp[["study_variable"]] <- sv_shared
      SummarizedExperiment::colData(qfeatures) <- cd_shared_tmp
    }
    cd_sml_tmp <- SummarizedExperiment::colData(qfeatures[["SML"]])
    if (length(cd_sml_tmp) && ("study_variable" %in% names(cd_sml_tmp))) {
      sv_sml <- cd_sml_tmp[["study_variable"]]
      names(sv_sml) <- rownames(cd_sml_tmp)
      cd_sml_tmp[["study_variable"]] <- sv_sml
      SummarizedExperiment::colData(qfeatures[["SML"]]) <- cd_sml_tmp
    }
  }
    SummarizedExperiment::rowData(qfeatures[["SML"]]) <- sml_rd
    if (all(lengths(refs_list) == 1L)) {
      # Fast path: exactly one SMF per SML 
      sml_rd$SMF_ID_REF <- unlist(refs_list, use.names = FALSE)
      SummarizedExperiment::rowData(qfeatures[["SML"]]) <- sml_rd
      qfeatures <- QFeatures::addAssayLink(
        qfeatures,
        from = "SMF",
        to = "SML",
        varFrom = col_smf_id,
        varTo = "SMF_ID_REF"
      )
    } else {
      # Build Hits graph for multi-mapping; drop unknown references
      smf_ids <- trimws(as.character(smf_rd[[col_smf_id]]))
      to_idx <- rep.int(seq_along(refs_list), lengths(refs_list))
      ref_tokens <- unlist(refs_list, use.names = FALSE)
      from_idx <- match(ref_tokens, smf_ids)
      keep_mask <- !is.na(from_idx)
      
      # Warn about unknown references
      if (any(!keep_mask)) {
        unknown_refs <- unique(ref_tokens[!keep_mask])
        refs_to_show <- if (length(unknown_refs) <= 5) {
          paste(unknown_refs, collapse = ", ")
        } else {
          paste(c(unknown_refs[1:5], "..."), collapse = ", ")
        }
        warning("Unknown references: ", refs_to_show)
      }
      
      hits <- S4Vectors::Hits(
        from = from_idx[keep_mask],
        to = to_idx[keep_mask],
        nLnode = nrow(qfeatures[["SMF"]]),
        nRnode = nrow(qfeatures[["SML"]])
      )
      assay_link <- QFeatures::AssayLink(
        name = "SML",
        from = "SMF",
        fcol = "SMF_ID_REFS",
        hits = hits
      )
      # Assign links
      qfeatures@assayLinks <- QFeatures::AssayLinks(assay_link)
    }
  }

  qfeatures
}
