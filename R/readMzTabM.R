#' Convert mzTab-M file to QFeatures object
#'
#' Parses an mzTab-M file and converts it to a QFeatures object with 
#' automatic metadata enrichment and linking between sections.
#'
#' @param file character(1) Path to the mzTab-M file.
#' @param load_sml logical(1) Include SML section and link to SMF. Default TRUE.
#'
#' @return QFeatures object with SMF (and optionally SML) assays. 
#'   Column names use sample labels when unique, otherwise assay labels or IDs.
#'   Metadata from MTD section is automatically merged into colData.
#'
#' @details
#' **Requirements**: mzTab-M file with MTD assay entries, SMF section containing
#' \code{SMF_ID} and \code{abundance_assay[n]} columns. Optional: sample, ms_run,
#' and study_variable metadata for enriched colData.
#' 
#' **Column Name Standardization**: The function automatically maps mzTab-M 
#' column names to MetFamily's expected format:
#' \itemize{
#'   \item colData: \code{study_variable_label} -> \code{Class} (for PCA/HCA grouping), 
#'         \code{Type} = "Sample"
#'   \item rowData: \code{SMF_ID} -> \code{Alignment ID}, 
#'         \code{exp_mass_to_charge} -> \code{Average Mz},
#'         \code{retention_time_in_seconds} -> \code{Average Rt(min)} 
#'         (converted from seconds to minutes),
#'         \code{chemical_name} -> \code{Metabolite name},
#'         \code{adduct_ion} -> \code{Adduct ion name},
#'         \code{inchi_key} -> \code{INCHIKEY},
#'         \code{smiles} -> \code{SMILES}
#' }
#' Optional columns that are never used in analysis (Fill %, Dot product, 
#' Fragment presence %, etc.) are not created. Missing columns are set to NA.
#' 
#' **Linking**: SML rows are linked to SMF via \code{SMF_ID_REFS} using either 
#' AssayLinks (1:1) or Hits objects (many-to-many).
#' 
#' **Error handling**: Strict validation with helpful errors for duplicate IDs,
#' missing required sections, and invalid references.
#'
#' @examples
#' \dontrun{
#' # Basic conversion with comprehensive metadata enrichment
#' qf <- readMzTabM("path/to/file.mzTab")
#' 
#' # SMF only (no SML section)
#' qf_smf_only <- readMzTabM("path/to/file.mzTab", load_sml = FALSE)
#' 
#' # Access enriched metadata
#' colData(qf)
#' colData(qf[["SMF"]])
#' 
#' # Access raw MTD section metadata
#' S4Vectors::metadata(qf)$mzTabM$mtd_raw
#' 
#' # Column names will be sample labels when available and unique,
#' # otherwise assay labels, or finally assay IDs as fallback
#' colnames(qf[["SMF"]])
#' }
#'
#' @export
readMzTabM <- function(file, load_sml = TRUE) {
  
  # Helper Functions ----
  # Internal utility functions for metadata parsing and data reshaping
  normalize_ref_tokens <- function(x) {
    if (!length(x) || is.na(x) || !nzchar(x)) return(character(0))
    trimws(unlist(strsplit(gsub(",", "|", x, fixed = TRUE), "|", fixed = TRUE), use.names = FALSE))
  }
  get_assay_quant_cols <- function(tab) {
    nms <- names(tab)
    keep <- startsWith(nms, "abundance_assay[") & endsWith(nms, "]")
    nms[keep]
  }
  # Extract entity index from mzTab keys (e.g., sample[2]-species[1] -> 2)
  extract_index <- function(keys) as.integer(sub("^[^\\[]*\\[(\\d+)\\].*", "\\1", keys))
  # Cast metadata to wide format while preserving array fields
  cast_wide <- function(df, idvar, field_col="field", value_col="value") {
    # Array-like fields that should be kept as separate columns (not cast to wide)
    array_patterns <- c("^scan_polarity\\[", "^cv\\[", "^userParam\\[", "^instrument_", "^software\\[")
    
    # Split into array fields and regular fields
    is_array_field <- sapply(df[[field_col]], function(field) {
      any(sapply(array_patterns, function(pattern) grepl(pattern, field)))
    })
    
    array_df <- df[is_array_field, , drop = FALSE]
    regular_df <- df[!is_array_field, , drop = FALSE]
    
    # Check for problematic duplicates in regular fields
    if (nrow(regular_df) > 0) {
      key <- paste(regular_df[[idvar]], regular_df[[field_col]], sep = "\r")
      if (anyDuplicated(key)) {
        dupe_keys <- key[duplicated(key)]
        dupe_entries <- gsub(".*\r", "", dupe_keys)
        stop("Duplicate metadata entries found for fields: ", 
             paste(unique(dupe_entries), collapse = ", "), 
             ". These fields should appear only once per entity.")
      }
      
      # Cast regular fields to wide
      regular_wide <- reshape(regular_df, idvar = idvar, timevar = field_col, v.names = value_col, direction = "wide")
    } else {
      # Create empty wide format with just the idvar
      regular_wide <- unique(df[idvar])
    }
    
    # Add array fields as individual columns (keeping all duplicates as separate columns)
    if (nrow(array_df) > 0) {
      # For array fields, create unique column names for duplicates
      array_df$unique_field <- make.unique(array_df[[field_col]])
      array_wide <- reshape(array_df, idvar = idvar, timevar = "unique_field", v.names = value_col, direction = "wide")
      
      # Merge regular and array fields
      result <- merge(regular_wide, array_wide, by = idvar, all = TRUE)
    } else {
      result <- regular_wide
    }
    
    result
  }

  # Read mzTab-M
  mztab_df <- rmzTabM::readMzTab(file)

  smf <- rmzTabM::extractSmallMoleculeFeatures(mztab_df)

  if (isTRUE(load_sml)) {
    sml <- rmzTabM::extractSmallMoleculeSummary(mztab_df)
  } else {
    sml <- NULL
  }

  mtd <- rmzTabM::extractMetadata(mztab_df)
  
  keys <- as.character(mtd[["V2"]])
  vals <- as.character(mtd[["V3"]])
  
  # Parse Assay Metadata ----
  # Extract and process assay-level metadata from MTD section
  # Creates assay_tab with assay_index, labels, sample_refs, ms_run_refs
  k_assay <- startsWith(keys, "assay[")
  if (!any(k_assay)) stop("No assay metadata in MTD.")
  a_idx <- extract_index(keys[k_assay])
  a_key <- sub("^assay\\[\\d+\\]-?", "", keys[k_assay])
  a_key[a_key == ""] <- "label"
  assay_df <- data.frame(assay_index = a_idx, field = a_key, value = vals[k_assay], stringsAsFactors = FALSE)
  assay_tab <- cast_wide(assay_df, idvar = "assay_index")
  names(assay_tab) <- sub("^value\\.", "", names(assay_tab))
  names(assay_tab)[names(assay_tab) == "label"] <- "assay_label"
  assay_tab$assay_id <- paste0("assay[", assay_tab$assay_index, "]")
  
  # Parse Sample Metadata  ----
  # Extract sample-level metadata: species, tissues, diseases, descriptions
  # Creates sample_tab for linking assays to biological samples
  k_sample <- startsWith(keys, "sample[")
  sample_tab <- NULL
  if (any(k_sample)) {
    s_idx <- extract_index(keys[k_sample])
    s_key <- sub("^sample\\[\\d+\\]-?", "", keys[k_sample]); s_key[s_key == ""] <- "label"
    sample_df <- data.frame(sample_index = s_idx, field = s_key, value = vals[k_sample], stringsAsFactors = FALSE)
    sample_tab <- cast_wide(sample_df, idvar = "sample_index")
    names(sample_tab) <- sub("^value\\.", "", names(sample_tab))
    names(sample_tab)[names(sample_tab) == "label"] <- "sample_label"
    sample_tab$sample_id <- paste0("sample[", sample_tab$sample_index, "]")
  }
  
  # Parse MS Run Metadata  ----
  # Extract MS run information: file locations, formats, scan polarities
  # Creates msr_tab for linking assays to raw data files
  k_msr <- startsWith(keys, "ms_run[")
  msr_tab <- NULL
  if (any(k_msr)) {
    r_idx <- extract_index(keys[k_msr])
    r_key <- sub("^ms_run\\[\\d+\\]-?", "", keys[k_msr]); r_key[r_key == ""] <- "label"
    msr_df <- data.frame(ms_run_index = r_idx, field = r_key, value = vals[k_msr], stringsAsFactors = FALSE)
    msr_tab <- cast_wide(msr_df, idvar = "ms_run_index")
    names(msr_tab) <- sub("^value\\.", "", names(msr_tab))
    names(msr_tab)[names(msr_tab) == "label"] <- "ms_run_label"
    msr_tab$ms_run_id <- paste0("ms_run[", msr_tab$ms_run_index, "]")
    # Avoid generic names like "location" clashing later:
    if ("location" %in% names(msr_tab)) names(msr_tab)[names(msr_tab) == "location"] <- "ms_run_location"
    if ("format"   %in% names(msr_tab)) names(msr_tab)[names(msr_tab) == "format"]   <- "ms_run_format"
    if ("id_format" %in% names(msr_tab)) names(msr_tab)[names(msr_tab) == "id_format"] <- "ms_run_id_format"
    # scan_polarity[n] columns will be kept as-is (one col per index)
  }
  
  # Parse Study Variables ----
  # Extract experimental conditions/factors that group assays
  # Creates sv_map to link assays to study design (control/treatment groups)
  k_sv <- startsWith(keys, "study_variable[")
  sv_map <- NULL
  if (any(k_sv)) {
    sv_idx  <- extract_index(keys[k_sv])
    sv_key  <- sub("^study_variable\\[\\d+\\]-?", "", keys[k_sv])
    sv_key[sv_key == ""] <- "label"
    sv_df <- data.frame(study_variable_index = sv_idx, field = sv_key, value = vals[k_sv], stringsAsFactors = FALSE)
    # labels
    sv_lab <- sv_df[sv_df$field == "label", c("study_variable_index", "value"), drop = FALSE]
    names(sv_lab)[2] <- "study_variable_label"
    # assay membership
    sv_refs <- sv_df[sv_df$field == "assay_refs", c("study_variable_index", "value"), drop = FALSE]
    if (nrow(sv_refs)) {
      sv_map <- do.call(rbind, lapply(seq_len(nrow(sv_refs)), function(i) {
        ai <- normalize_ref_tokens(sv_refs$value[i])
        if (!length(ai)) return(NULL)
        data.frame(
          assay_id = ai,
          study_variable_index = sv_refs$study_variable_index[i],
          stringsAsFactors = FALSE
        )
      }))
      if (!is.null(sv_map) && nrow(sv_lab)) {
        sv_map <- merge(sv_map, sv_lab, by = "study_variable_index", all.x = TRUE)
      }
    }
  }
  

  # Check assay-level quantification columns
  smf_assay_quant_cols <- get_assay_quant_cols(smf)
  if (length(smf_assay_quant_cols) == 0) {
    stop("missing abundance_assay[...] columns in mzTab-M file.")
  }
  
  # Build Enriched Column Metadata ----
  # Merge assay, sample, ms_run, and study_variable metadata
  if ("sample_ref" %in% names(assay_tab)) assay_tab$sample_id <- assay_tab$sample_ref
  if ("ms_run_ref" %in% names(assay_tab)) assay_tab$ms_run_id <- assay_tab$ms_run_ref
  if (!is.null(sample_tab) && "sample_id" %in% names(assay_tab))
    assay_tab <- merge(assay_tab, sample_tab[, c("sample_id", setdiff(names(sample_tab), "sample_index"))], by = "sample_id", all.x = TRUE)
  if (!is.null(msr_tab) && "ms_run_id" %in% names(assay_tab))
    assay_tab <- merge(assay_tab, msr_tab[, c("ms_run_id", setdiff(names(msr_tab), "ms_run_index"))], by = "ms_run_id", all.x = TRUE)
  if (!is.null(sv_map))
    assay_tab <- merge(assay_tab, sv_map[, c("assay_id", "study_variable_index", "study_variable_label")], by = "assay_id", all.x = TRUE)
  
  assay_tab <- assay_tab[order(assay_tab$assay_index), ]
  # Column naming hierarchy: sample labels (biological context) > assay labels > assay IDs
  # Fall back if sample labels are duplicated (technical replicates case)
  if ("sample_label" %in% names(assay_tab) && all(nzchar(assay_tab$sample_label)) && !anyDuplicated(assay_tab$sample_label)) {
    colnames_choice <- assay_tab$sample_label
  } else if ("assay_label" %in% names(assay_tab) && all(nzchar(assay_tab$assay_label))) {
    if (anyDuplicated(assay_tab$assay_label)) {
      stop("Duplicate assay labels found in mzTab-M file. Assay labels must be unique for column naming.")
    }
    colnames_choice <- assay_tab$assay_label
  } else {
    colnames_choice <- assay_tab$assay_id
  }
  
  # Final check 
  if (anyDuplicated(colnames_choice)) {
    stop("Duplicate column names detected after processing. This indicates a problem with the mzTab-M file structure.")
  }
  
  keep_cd <- intersect(
    names(assay_tab),
    c("assay_id","assay_label","assay_index",
      "sample_id","sample_label",
      "ms_run_id","ms_run_label","ms_run_location","ms_run_format","ms_run_id_format",
      "study_variable_index","study_variable_label",
      grep("^scan_polarity", names(assay_tab), value = TRUE),
      "external_uri")
  )
  colData <- S4Vectors::DataFrame(assay_tab[, keep_cd, drop = FALSE], row.names = colnames_choice)



  # Build SummarizedExperiment for a section
  build_se <- function(tab, id_col, cd) {
    qcols <- get_assay_quant_cols(tab)
    if (!length(qcols)) stop("No abundance_assay[...] columns in ", id_col, " section.")
    idx <- extract_index(qcols); ord <- order(idx); qcols <- qcols[ord]; idx <- idx[ord]
    map <- match(idx, assay_tab$assay_index)
    if (anyNA(map)) stop("Quant columns refer to assay indices not present in MTD assay list.")
    cn <- rownames(cd)[map]
    tab[qcols] <- lapply(tab[qcols], function(x) if (is.double(x)) x else suppressWarnings(as.double(x)))
    m <- as.matrix(tab[qcols]); colnames(m) <- cn
    if (!(id_col %in% names(tab))) stop("Missing ID column ", id_col, " in section.")
    rid <- as.character(tab[[id_col]])
    if (anyDuplicated(rid)) stop("Duplicate ", id_col, " values found.")
    rownames(m) <- rid
    # correct column subsetting for rowData
    rd <- S4Vectors::DataFrame(tab[, setdiff(names(tab), qcols), drop = FALSE])
    SummarizedExperiment::SummarizedExperiment(assays = list(intensity = m), rowData = rd, colData = cd)
  }
  
  se_small_molecule_features <- build_se(smf, "SMF_ID", colData)
  se_small_molecule_summary <- if (isTRUE(load_sml)) build_se(sml, "SML_ID", colData) else NULL

  
  # Construct QFeatures
  qfeatures <- if (isTRUE(load_sml)) {
    QFeatures::QFeatures(experiments = list(SMF = se_small_molecule_features, SML = se_small_molecule_summary), colData = colData)
  } else {
    QFeatures::QFeatures(experiments = list(SMF = se_small_molecule_features), colData = colData)
  }

  # Standardize Column Names for MetFamily Compatibility ----
  # MetFamily expects specific column names that differ from mzTab-M standard
  
  # Standardize colData: Add Class and Type columns
  cd <- SummarizedExperiment::colData(qfeatures)
  
  # Class: Required for PCA/HCA group assignment
  if ("study_variable_label" %in% names(cd) && !all(is.na(cd$study_variable_label))) {
    cd$Class <- cd$study_variable_label
  } else {
    cd$Class <- "Unknown"
    warning("No study_variable_label found in mzTab-M. Setting Class to 'Unknown'.")
  }
  
  # Type: Expected by parser, default to Sample
  cd$Type <- "Sample"
  
  SummarizedExperiment::colData(qfeatures) <- cd
  
  # Standardize rowData for SMF section
  rd_smf <- SummarizedExperiment::rowData(qfeatures[["SMF"]])
  
  # Critical columns with transformations
  rd_smf$`Alignment ID` <- rd_smf$SMF_ID
  rd_smf$`Average Mz` <- rd_smf$exp_mass_to_charge
  
  # Retention time: Convert seconds to minutes
  if ("retention_time_in_seconds" %in% names(rd_smf)) {
    rd_smf$`Average Rt(min)` <- as.numeric(rd_smf$retention_time_in_seconds) / 60
  } else if ("retention_time_in_seconds_start" %in% names(rd_smf)) {
    # Fallback: use start time if main RT missing
    rd_smf$`Average Rt(min)` <- as.numeric(rd_smf$retention_time_in_seconds_start) / 60
    warning("Using retention_time_in_seconds_start as fallback for Average Rt(min)")
  } else {
    rd_smf$`Average Rt(min)` <- NA_real_
    warning("No retention time found in mzTab-M SMF section")
  }
  
  # Metabolite identification (used in display)
  rd_smf$`Metabolite name` <- if ("chemical_name" %in% names(rd_smf)) {
    rd_smf$chemical_name
  } else {
    NA_character_
  }
  
  # Adduct (used for MetFrag links)
  rd_smf$`Adduct ion name` <- if ("adduct_ion" %in% names(rd_smf)) {
    rd_smf$adduct_ion
  } else {
    NA_character_
  }
  
  # Structure identifiers (if available)
  rd_smf$INCHIKEY <- if ("inchi_key" %in% names(rd_smf)) {
    rd_smf$inchi_key
  } else {
    NA_character_
  }
  
  rd_smf$SMILES <- if ("smiles" %in% names(rd_smf)) {
    rd_smf$smiles
  } else {
    NA_character_
  }
  
  SummarizedExperiment::rowData(qfeatures[["SMF"]]) <- rd_smf
  
  # Standardize SML if present (minimal set)
  if ("SML" %in% names(qfeatures)) {
    rd_sml <- SummarizedExperiment::rowData(qfeatures[["SML"]])
    rd_sml$`Alignment ID` <- rd_sml$SML_ID
    rd_sml$`Average Mz` <- if ("exp_mass_to_charge" %in% names(rd_sml)) {
      rd_sml$exp_mass_to_charge
    } else {
      NA_real_
    }
    rd_sml$`Metabolite name` <- if ("chemical_name" %in% names(rd_sml)) {
      rd_sml$chemical_name
    } else {
      NA_character_
    }
    SummarizedExperiment::rowData(qfeatures[["SML"]]) <- rd_sml
  }

  # Create SML to SMF Linkages ----
  # Process SMF_ID_REFS to link small molecule summaries to features
  # Handles both 1:1 and many-to-many relationships using AssayLinks or Hits 
  if (isTRUE(load_sml)) {
    smf_rd <- SummarizedExperiment::rowData(qfeatures[["SMF"]])
    sml_rd <- SummarizedExperiment::rowData(qfeatures[["SML"]])
    refs_chr <- if ("SMF_ID_REFS" %in% names(sml_rd)) as.character(sml_rd[["SMF_ID_REFS"]]) else rep("", nrow(sml_rd))
    refs_list <- IRanges::CharacterList(lapply(refs_chr, normalize_ref_tokens))
    sml_rd$SMF_ID_REFS <- refs_list
    SummarizedExperiment::rowData(qfeatures[["SML"]]) <- sml_rd
    
    smf_ids <- as.character(smf_rd[["SMF_ID"]])
    to_idx <- rep.int(seq_along(refs_list), lengths(refs_list))
    ref_tokens <- unlist(refs_list, use.names = FALSE)
    from_idx <- match(ref_tokens, smf_ids)
    keep <- !is.na(from_idx)
    
    # Error on unknown references
    if (any(!keep)) {
      unknown_refs <- unique(ref_tokens[!keep])
      stop("Unknown SMF_ID references found in SML section: ", paste(unknown_refs, collapse = ", "), 
           ". All SMF_ID_REFS must correspond to existing SMF_ID values.")
    }
    
    if (all(lengths(refs_list) == 1L)) {
      # 1:1 linking
      sml_rd$SMF_ID_REF <- unlist(refs_list, use.names = FALSE)
      SummarizedExperiment::rowData(qfeatures[["SML"]]) <- sml_rd
      qfeatures <- QFeatures::addAssayLink(qfeatures, from = "SMF", to = "SML", varFrom = "SMF_ID", varTo = "SMF_ID_REF")
    } else {
      # Many-to-many linking
      hits <- S4Vectors::Hits(from = from_idx, to = to_idx,
                              nLnode = nrow(qfeatures[["SMF"]]), nRnode = nrow(qfeatures[["SML"]]))
      link <- QFeatures::AssayLink(name = "SML", from = "SMF", fcol = "SMF_ID_REFS", hits = hits)
      qfeatures@assayLinks <- QFeatures::AssayLinks(link)
    }
  }
  
  # Store Raw Metadata ----
  # Preserve complete MTD section for reference 
  meta <- list()
  meta$mtd_raw <- mtd  # Save complete MTD section 
  
  # set metadata via S4Vectors setter
  md <- S4Vectors::metadata(qfeatures)
  md$mzTabM <- meta
  S4Vectors::metadata(qfeatures) <- md

  qfeatures
}
