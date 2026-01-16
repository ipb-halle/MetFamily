# fileSpectra <- file.path("../file-formats/Metaboscape-export-version 2025b",
#                          "UTH-2025-07-29-UTH003_002-conyza-test-samples.sirius.mgf")
# file.exists(fileSpectra)


# TODO create parseMS2 fct which checks file extension and use correct parse_to_list
# Works for MetaboScape MGF files. MS-Dial format is different, with nested header entries

#' Read MGF files
#'
#' Reader for Mascot Generic Format MS2 files. Works similar to `parseMSP_rewrite()`.
#'
#' @inheritParams parseMSP_rewrite
#'
#' @returns A list including a spectraList element and additional stats.
#' @export
parseMGF <- function(
    fileSpectra, 
    minimumIntensityOfMaximalMS2peak = 2000,
    minimumProportionOfMS2peaks = 0.05,
    neutralLossesPrecursorToFragments = T,
    neutralLossesFragmentsToFragments = F,
    mz_tol = 0.17,
    min_frag_nb = 2,
    min_intensity = 500,
    progress = FALSE
) {
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file...", sep = "")) else print(paste("Parsing MS/MS file...", sep = ""))
  
  mgf_entries_parsed <- parseMGF_to_list(fileSpectra)
  
  # filter, add NL
  # can use the msp filter fct directly
  msp_filter_results <- filter_msp_entries(
    msp_entries = mgf_entries_parsed,
    minimumIntensityOfMaximalMS2peak, 
    minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments,  
    neutralLossesFragmentsToFragments,
    mz_tol = mz_tol,
    min_frag_nb = min_frag_nb,
    min_intensity = min_intensity
  )
  
  # format for backward compatibility
  c(
    fileSpectra = fileSpectra,
    msp_filter_results
  )
  
}


#' Process MGF from file to list of entries
#' 
#' Create list object from MGF entry, keeping all headers and values
#' Rename headers based on standards
#'
#' @param fileSpectra `character(1)` file path
#'
#' @returns list of entries
#' @export
parseMGF_to_list <- function(fileSpectra) {
  
  message("MS/MS file: Read file")
  
  # Note: squish also replaces tab with whitespace
  fileLines <- readLines(con = fileSpectra) %>% stringr::str_squish()
  
  message("MS/MS file: Parse")
  
  # delimiters
  begin_lines <- stringr::str_equal(fileLines, "BEGIN IONS")
  end_lines <- stringr::str_equal(fileLines, "END IONS")
  level_lines <- stringr::str_detect(fileLines, "MSLEVEL=")
  
  
  # subset only ms2 content
  lvl_ms2 <- stringr::str_equal(fileLines[level_lines], "MSLEVEL=2")
  from_ms2 <- which(begin_lines)[lvl_ms2] # include begin lines
  to_ms2 <- which(end_lines)[lvl_ms2] - 1
  stopifnot(length(from_ms2) == length(to_ms2) && all(from_ms2 < to_ms2))
  
  keep_lines <- unlist(purrr::map2(from_ms2, to_ms2, `:`))
  lines_ms2 <- fileLines[keep_lines]
  lines_ms2 <- lines_ms2[-1]
  
  isLineEmpty <- stringr::str_equal(lines_ms2, "BEGIN IONS")
  lines_ms2 <- stringr::str_replace(lines_ms2, "^BEGIN IONS$", "")
  
  ## determine line contents
  isLineHeader <- stringr::str_detect(lines_ms2, "=")
  isLinePeaks <- stringr::str_detect(lines_ms2, "^[\\d\\. ]+$")
  # ~ used to be "^\\d+(((\\.)|(,))\\d+)?[ \t]\\d+(((\\.)|(,))\\d+)?$"
  
  nbLines <- length(lines_ms2)
  
  ## determine entry intervals
  firstLines <- c(1, which(isLineEmpty) + 1)
  lastLines <- c(which(isLineEmpty) - 1, nbLines)
  
  # checks
  stopifnot("The content of some lines in the MGF file could not be determined." = 
              unique(isLineEmpty+isLineHeader+isLinePeaks) == 1)
  stopifnot("It is not allowed to have two consecutive empty lines in the MSP file." = 
              !any(c(isLineEmpty, F) & c(F, isLineEmpty)))
  
  ## create named vector of headers
  
  headers0 <- lines_ms2[isLineHeader]
  head_labels0 <- stringr::str_extract(headers0, "^(.*?) ?[:=] ?(.*)$", group = 1)
  head_values <- stringr::str_extract(headers0, "^(.*?) ?[:=] ?(.*)$", group = 2)
  
  
  # use standardized names
  uni_labels <- unique(head_labels0)
  uni_low_labels <- stringr::str_to_lower(uni_labels)
  
  # is rt in min or sec?
  # change values below
  rt_in_sec <- if ("rtinseconds" %in% uni_low_labels) TRUE else FALSE
    
  mgf_dict_table <- mgf_dict()
  uni_stand_labels <- unname(mgf_dict_table[match(uni_low_labels, names(mgf_dict_table))])
  
  
  # keep original if no standard exists
  uni_stand_labels[is.na(uni_stand_labels)] <- uni_labels[is.na(uni_stand_labels)]
  
  # build back long vector
  head_labels <- uni_stand_labels[match(head_labels0, uni_labels)]
  
  headers_long <- character(nbLines)
  headers_long[isLineHeader] <- head_values
  names(headers_long)[isLineHeader] <- head_labels
  
  
  ## peaks data
  
  peaks0 <- lines_ms2[isLinePeaks]
  
  mz <- stringr::str_extract(peaks0, "^(.*) (.*)$", group = 1) %>% as.numeric
  mz_long <- numeric(nbLines)
  mz_long[isLinePeaks] <- mz
  
  int <- stringr::str_extract(peaks0, "^(.*) (.*)$", group = 2) %>% as.numeric
  int_long <- numeric(nbLines)
  int_long[isLinePeaks] <- int
  
  
  # return split entries
  mgf_entries <- purrr::map(seq_along(firstLines), ~{
    
    fl <- firstLines[.x]
    ll <- lastLines[.x]
    
    lineRange <- seq(fl, ll)
    
    # 3 vectors to pull from
    # is there a more clever/pretty way?
    entr <- tibble::lst(
      header = headers_long[lineRange][isLineHeader[lineRange]],
      mz = mz_long[lineRange][isLinePeaks[lineRange]],
      int = int_long[lineRange][isLinePeaks[lineRange]],
      spectrumString = paste(mz, int, sep = " ", collapse = ";"),
      # this won't match in original file as MS1 spectra were removed
      # FIX get the real line numbers
      entry_range = c(fl, ll),
      # changed for mfg as some entries skipped
      id = as.numeric(header["feature_id"])
    )
    
    entr$header["peakNumber"] <- length(entr$mz)
    
    entr
    
  })
  
  # fix rt if in sec
  if (rt_in_sec) {
    
    mgf_entries <- purrr::map(
      mgf_entries, \(entry) {
        
        entry$header["rtime"] <- as.character(
          as.numeric(entry$header["rtime"]) / 60
        )
        
        entry
      }
    )
    
  }

  mgf_entries
  
}


#' MGF header dictionnary
#' 
#' Name matching for stadardized metadata.
#' Standard names based on `Spectra::spectraVariableMapping(MsBackendMsp::MsBackendMsp())`
#'
#' @returns named vector with standardized names as value
#' @export
mgf_dict <- function() {
  tibble::enframe(
    list(
      rtime = c("rtinseconds", "rtinminutes"), #rt
      precursorMz = c("pepmass"), #mz
      feature_id = "feature_id",
      adduct = c("ion"),
      scans = "scans",
      charge = "charge",
      polarity = "polarity"
    )
  ) %>%
    tidyr::unnest(value) %>% dplyr::select(2,1) %>% tibble::deframe()
}
