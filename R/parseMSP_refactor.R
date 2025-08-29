# Rewrite parseMSP ----
# Compartmentalized, easier to debug, faster

#' Read and filter MSP spectra files
#' 
#' Drop-in replacement for parseMSP, just better
#'
#' @param fileSpectra `character(1)` file path
#' @param minimumIntensityOfMaximalMS2peak `numeric(1)` Spectra entries with a
#'   maximum intensity below this value are removed.
#' @param minimumProportionOfMS2peaks `numeric(1)` Between 0 and 1. Within a
#'   spectra, fragments with an intensity below this proportion of the maximum
#'   intensity fragment are removed.
#' @param neutralLossesPrecursorToFragments boolean Should neutral losses to
#'   precursor be calculated.
#' @param neutralLossesFragmentsToFragments boolean Should neutral losses for
#'   each fragment pair be calculated.
#' @param progress boolean, used for interactive app.
#' @param mz_tol `numeric(1)` tolerance used for matching precursor fragment.
#'
#' @returns A list including a spectraList element and additional stats.
#' @export
parseMSP_rewrite <- function(
    fileSpectra, 
    minimumIntensityOfMaximalMS2peak, 
    minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments, 
    neutralLossesFragmentsToFragments,
    mz_tol = 0.17,
    min_frag_nb = 1,
    min_intensity = 0,
    progress = FALSE
) {
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.01, detail = paste("Parsing MS/MS file...", sep = "")) else print(paste("Parsing MS/MS file...", sep = ""))
  
  msp_entries_parsed <- parseMSP_to_list(fileSpectra)
  
  # filter, add NL
  # main workhorse
  msp_filter_results <- filter_msp_entries(
    msp_entries_parsed,
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


#' Process MSP from file to list of entries
#' 
#' Create list object from MSP entry, keeping all headers and values
#' Rename headers based on standards
#'
#' @param fileSpectra `character(1)` file path
#'
#' @returns list of entries
#' @export
parseMSP_to_list <- function(fileSpectra) {
  
  message("MS/MS file: Read file")
  
  # Note: squish also replaces tab with whitespace
  fileLines <- readLines(con = fileSpectra) %>% stringr::str_squish()
  isLineEmpty <- stringr::str_equal(fileLines, "")
  
  message("MS/MS file: Parse")
  
  ## remove empty lines at beginning and end
  nbLines <- length(fileLines)
  firstLine <- which.max(!isLineEmpty)
  lastLine <- nbLines - which.max(!rev(isLineEmpty)) + 1
  
  fileLines <- fileLines[firstLine:lastLine]
  nbLines <- length(fileLines)
  
  ## determine line contents
  
  isLineEmpty <- stringr::str_equal(fileLines, "")
  isLineHeader <- stringr::str_detect(fileLines, ":|=")
  isLinePeaks <- stringr::str_detect(fileLines, "^[\\d\\. ]+$")
  # ~ used to be "^\\d+(((\\.)|(,))\\d+)?[ \t]\\d+(((\\.)|(,))\\d+)?$"
  
  ## determine entry intervals
  firstLines <- c(1, which(isLineEmpty) + 1)
  lastLines <- c(which(isLineEmpty) - 1, nbLines)
  
  # checks
  stopifnot("The content of some lines in the MSP file could not be determined." = 
              unique(isLineEmpty+isLineHeader+isLinePeaks) == 1)
  stopifnot("It is not allowed to have two consecutive empty lines in the MSP file." = 
              !any(c(isLineEmpty, F) & c(F, isLineEmpty)))
  # similar to
  # all(firstLines[-1] == lastLines[-length(lastLines)] + 2)
  
  ## create named vector of headers
  
  headers0 <- fileLines[isLineHeader]
  head_labels0 <- stringr::str_extract(headers0, "^(.*?) ?[:=] ?(.*)$", group = 1)
  head_values <- stringr::str_extract(headers0, "^(.*?) ?[:=] ?(.*)$", group = 2)
  
  
  # use standardized names
  uni_labels <- unique(head_labels0)
  uni_low_labels <- stringr::str_to_lower(uni_labels)
  
  # not implemented, would need to see example file
  if(any(uni_low_labels %in% c("begin ion", "scannumber", "modelion"))) {
    stop(paste0("Label '", label, "' found in the MSP file is not currently implemented. Please contact us."))
  }
  
  msp_dict_table <- msp_dict()
  uni_stand_labels <- unname(msp_dict_table[match(uni_low_labels, names(msp_dict_table))])
  
  
  # keep original if no standard exists
  uni_stand_labels[is.na(uni_stand_labels)] <- uni_labels[is.na(uni_stand_labels)]
  
  # build back long vector
  head_labels <- uni_stand_labels[match(head_labels0, uni_labels)]
  
  headers_long <- character(nbLines)
  headers_long[isLineHeader] <- head_values
  names(headers_long)[isLineHeader] <- head_labels
  
  
  ## peaks data
  
  peaks0 <- fileLines[isLinePeaks]
  
  mz <- stringr::str_extract(peaks0, "^(.*) (.*)$", group = 1) %>% as.numeric
  mz_long <- numeric(nbLines)
  mz_long[isLinePeaks] <- mz
  
  int <- stringr::str_extract(peaks0, "^(.*) (.*)$", group = 2) %>% as.numeric
  int_long <- numeric(nbLines)
  int_long[isLinePeaks] <- int
  
  
  # return split entries
  msp_entries <- purrr::map(seq_along(firstLines), ~{
    
    fl <- firstLines[.x]
    ll <- lastLines[.x]
    
    lineRange <- seq(fl, ll)
    
    # 3 vectors to pull from
    # is there a more clever/pretty way?
    tibble::lst(
      header = headers_long[lineRange][isLineHeader[lineRange]],
      mz = mz_long[lineRange][isLinePeaks[lineRange]],
      int = int_long[lineRange][isLinePeaks[lineRange]],
      spectrumString = paste(mz, int, sep = " ", collapse = ";"),
      entry_range = c(fl, ll),
      id = .x
    )
  })
  
  msp_entries
  
}


#' MSP header dictionnary
#' 
#' Name matching for stadardized metadata.
#' Standard names based on `Spectra::spectraVariableMapping(MsBackendMsp::MsBackendMsp())`
#'
#' @returns named vector with standardized names as value
#' @export
msp_dict <- function() {
  tibble::enframe(
    list(
      name = c("name", "title"),
      # ms1Int = NA,
      rtime = c("retention time", "retentiontime", "rtinseconds"), #rt
      precursorMz = c("precursormz", "precursor m/z", "pepmass"), #mz
      # totalMass needs post-processing based on adduct type
      totalMass = c("total exact mass", "exactmass", "exact mass", "exact_mass"),
      metName = "metabolitename",
      adduct = c("adductionname", "precursor type", "precursortype"),
      # quantMass = "modelion",
      compoundClass = "compound class",
      instrument = c("instrumenttype", "instrument type", "instrument", "source_instrument"), #intrumentType
      inchi = "inchi",
      inchikey = c("inchikey", "inchiaux"), #inchiKey
      smiles = "smiles",
      peakNumber = "num peaks" # but used to return length(ms2Peaks_mz)
    )
  ) %>%
    tidyr::unnest(value) %>% dplyr::select(2,1) %>% tibble::deframe()
}


#' Filter MSP entries
#' 
#' Should be used after `parseMSP_to_list()` to filter entries and MS2 spectra
#' based on a set of parameters.#' 
#'
#' @param msp_entries list of MSP entries created with `parseMSP_to_list()`
#' @inheritParams parseMSP_rewrite
#'
#' @returns list of spectra object, and stats
#' @export
#'
#' @importFrom purrr map map_dbl
filter_msp_entries <- function(
    msp_entries, 
    minimumIntensityOfMaximalMS2peak, 
    minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments, 
    neutralLossesFragmentsToFragments,
    mz_tol = 0.17,
    min_frag_nb = 1,
    min_intensity = 0
) {
  
  # for debugging
  # msp_entries <- msp_entries_parsed
  
  msp_entries_ori <- msp_entries
  
  message("MS/MS file: Filter")
  
  # remove entries too few peaks
  nb_peaks_ori <- purrr::map(msp_entries, ~.x$header["peakNumber"]) %>% as.numeric
  # ~ same as
  # nbPeaks <- map_dbl(msp_entries, ~length(.x$mz))
  
  spectra_no_peaks <- nb_peaks_ori == 0
  msp_entries <- msp_entries[!spectra_no_peaks]
  
  
  ## remove heavy fragments here ----
  # (better here than in filter_msp_entry)
  precursorMz <- purrr::map(msp_entries, ~.x$header["precursorMz"]) %>% as.numeric
  
  heavy_frags <- purrr::map2(msp_entries, precursorMz, ~.x$mz > (.y + mz_tol))
  
  # rm entry if only too heavy frags
  disc_too_heavy <- purrr::map_lgl(heavy_frags, ~all(.x))
  
  msp_entries <- msp_entries[!disc_too_heavy]
  
  # update heavy_frags index
  heavy_frags_sub <- heavy_frags[!disc_too_heavy]
  
  # remove heavy frags
  msp_entries <- purrr::map2(msp_entries, heavy_frags_sub, ~{
    .x$mz <- .x$mz[!.y]
    .x$int <- .x$int[!.y]
    .x
  })
  
  
  ## discard if max intensity below threshold ----
  max_int <- purrr::map_dbl(msp_entries, ~max(.x$int))
  disc_max_int <- max_int < minimumIntensityOfMaximalMS2peak
  
  nb_frag_disc_max_int <- sum(map_dbl(msp_entries[disc_max_int], ~length(.x$mz)))
  
  msp_entries <- msp_entries[!disc_max_int]
  
  message("MS/MS file: Assemble spectra")
  
  # process spectra separately ----
  proc_entries <- purrr::map(msp_entries, ~ process_msp_entry(
    .x,
    minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments, 
    neutralLossesFragmentsToFragments,
    min_intensity
  )
  )
  
  nb_too_low_frags <- purrr::map_dbl(proc_entries, "nb_too_low_frags")
  nb_nat_frags <- purrr::map_dbl(proc_entries, "nb_nat_frags")
  nb_full_frags <- purrr::map_dbl(proc_entries, ~.x$spectra_item[["peakNumber"]])
  
  
  # remove entries with too few remaining fragments
  # TODO update peaks stats
  too_few_peaks <- nb_nat_frags < min_frag_nb
  proc_entries <- proc_entries[!too_few_peaks]
  
  
  # remove individual stats blocks
  spectra_items <- purrr::map(proc_entries, "spectra_item")
  
  # return filtered_spectra + stats
  list(
    spectraList = spectra_items,
    numberOfSpectra = length(proc_entries),
    numberOfSpectraOriginal = length(msp_entries_ori),
    numberOfMS2PeaksOriginal = sum(nb_peaks_ori), 
    numberOfMS2PeaksWithNeutralLosses = sum(nb_full_frags),
    numberOfMS2PeaksAboveThreshold = sum(nb_nat_frags), # final peaks w/o NL
    numberOfMS2PeaksBelowThreshold = nb_frag_disc_max_int + sum(nb_too_low_frags),
    numberOfTooHeavyFragments = sum(unlist(heavy_frags)),
    numberOfSpectraDiscardedDueToNoPeaks = sum(spectra_no_peaks) + sum(too_few_peaks),
    numberOfSpectraDiscardedDueToMaxIntensity = sum(disc_max_int),
    numberOfSpectraDiscardedDueToTooHeavy = sum(disc_too_heavy),
    precursorMz = map_dbl(proc_entries, ~.x$spectra_item[["mz"]]),
    precursorRt = map_dbl(proc_entries, ~.x$spectra_item[["rt"]])
  )
  
}


#' Filter single MSP entry
#' 
#' nb of peaks and peaks with NL are updated in header
#' peaks are filtered, and neutral losses added if requested
#' stats section added
#' returns backwards-compatible spectraItem as used by MetFamily
#' 
#' @param msp_entry MSP entry object
#' @inheritParams parseMSP_rewrite
#'
#' @returns filtered spectra
#' @export
process_msp_entry <- function(
    msp_entry,
    minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments, 
    neutralLossesFragmentsToFragments,
    min_intensity = 0
) {
  
  # msp_entry <- msp_entries[[232]]
  
  header <- msp_entry$header
  head_names <- names(header)
  # NOTE: it's prettier to work with peak data in a dataframe,
  # but it's take about 10x longer
  mz <- msp_entry$mz
  int <- msp_entry$int
  
  # remove low intensity fragments ----
  
  # minimum intensity abs
  low_int_abs <- int < min_intensity
  nb_low_int_abs <- sum(low_int_abs)
  
  # normalize peak heights
  int <- int / max(int)
  
  # low frags rel
  low_int_rel <- int < minimumProportionOfMS2peaks
  nb_low_int_rel <- sum(low_int_rel)
  
  # subset
  # this logic helps to avoid skyscraper peaks to overshadow smaller
  # but still relevant peaks
  low_both_cases <- low_int_abs & low_int_rel
  mz <- mz[!low_both_cases]
  int <- int[!low_both_cases]
  
  
  ## add neutral losses ----
  
  if(neutralLossesPrecursorToFragments) {
    # This uses precursor mz from MS1
    mz_nl <- mz - as.numeric(header["precursorMz"])
    int_nl <- int
  } else {
    mz_nl <- NULL; int_nl <- NULL
  }
  
  if(neutralLossesFragmentsToFragments) {
    # This uses MS2 peaks only
    m_mz  <- outer(X = mz,  Y = mz,  FUN = function(x,y){x-y})
    m_int <- outer(X = int, Y = int, FUN = function(x,y){(x+y) / 2})
    upper <- upper.tri(x = m_mz)
    mz_nlff <- m_mz[upper]
    int_nlff <- m_int[upper]
  } else {
    mz_nlff <- NULL; int_nlff <- NULL
  }
  
  full_mz <- c(mz, mz_nl, mz_nlff)
  full_int <- c(int, int_nl, int_nlff)
  
  
  # return filtered MS2 + basic stats
  # use default values if non-existent
  spectra_item <- tibble::lst(
    name = if("name" %in% head_names) header[["name"]] else NULL,
    ms1Int = NA, # not implemented, need to see MSP example
    rt = if("rtime" %in% head_names) as.numeric(header[["rtime"]]) else 0,
    mz = if("precursorMz" %in% head_names)  round(as.numeric(header[["precursorMz"]]), digits = 4) else NULL,
    metName = if("metName" %in% head_names) header[["metName"]] else "Unknown",
    adduct = if("adduct" %in% head_names) header[["adduct"]] else "Unknown",
    quantMass = NA, # not implemented, need to see MSP example
    compoundClass = if("compoundClass" %in% head_names) header[["compoundClass"]] else "Unknown",
    instrumentType = if("instrument" %in% head_names) header[["instrument"]] else "Unknown",
    inchi = if("inchi" %in% head_names) header[["inchi"]] else "",
    inchiKey = if("inchikey" %in% head_names) header[["inchikey"]] else "",
    smiles = if("smiles" %in% head_names) header[["smiles"]] else "",
    peakNumber = length(full_mz),
    ms2Peaks_mz = full_mz,
    ms2Peaks_int = full_int,
    spectrumString = msp_entry$spectrumString,
    entryInterval = msp_entry$entry_range, # would need + c(0,1) for backwards compatibility
    id = msp_entry$id
  )
  
  tibble::lst(
    spectra_item,
    nb_too_low_frags = nb_low_int_abs + nb_low_int_rel,
    nb_nat_frags = length(mz)
  )
}





# parseMSP using Spectra [WIP] ----


# fileSpectra <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")


#' 
#' library(Spectra)
#' 
#' #' Drop-in replacement for 
#' #'
#' #' @param fileLines 
#' #' @param minimumIntensityOfMaximalMS2peak 
#' #' @param minimumProportionOfMS2peaks 
#' #' @param neutralLossesPrecursorToFragments 
#' #' @param neutralLossesFragmentsToFragments 
#' #' @param offset 
#' #' @param progress 
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#' parseMSP_spectra <- function(
#'     fileSpectra, 
#'     minimumIntensityOfMaximalMS2peak, 
#'     minimumProportionOfMS2peaks, 
#'     neutralLossesPrecursorToFragments, 
#'     neutralLossesFragmentsToFragments, 
#'     progress = FALSE
#' ) {
#'   
#'   
#'   t_spec <- system.time(
#'     sp <- Spectra::Spectra(fileSpectra, source = MsBackendMsp::MsBackendMsp())
#'   )
#'   
#'   t_mf <- system.time(
#'     sp_mf <- MetFamily:::parseMSP(fileSpectra,
#'                       minimumIntensityOfMaximalMS2peak = 1000, 
#'                       minimumProportionOfMS2peaks = 0.1, 
#'                       neutralLossesPrecursorToFragments = T, 
#'                       neutralLossesFragmentsToFragments = F, 
#'                       progress = FALSE)
#'   )
#'   
#'   Spectra::spectraVariables(sp)
#'   
#'   sp$msLevel
#'   sp$name %>% table
#'   
#'   
#'   sp %>% class
#'   sp %>% length
#'   sp %>% lengths
#'   sp[1]
#'   sp %>% class
#'   
#'   object.size(sp)
#'   format(object.size(sp_mf), units = "auto")
#'   # 
#'   
#'   sp@backend@spectraData@listData %>% View
#'   sp$adduct %>% table
#'   table(sp$ADDUCTIONNAME)
#' 
#'   Spectra::spectraVariableMapping(MsBackendMsp::MsBackendMsp())
#'   
#'   system.time(
#'   rmsp <- MsBackendMsp::readMsp(fileSpectra)
#'   )
#'   
#'   table(collisionEnergy(sp))
#'   
#'   precursorCharge(sp) %>% table
#'   precursorIntensity(sp) %>% table
#'   precursorMz(sp) %>% head
#'   
#'   peaksData(sp[4])[[1]]
#'   peaksData(sp[4])@listData[[1]]
#'   isEmpty(sp) %>% table
#'   
#'   as(sp[1:5], "SimpleList")
#'   
#'   rmsp %>% names
#' 
#'   rmsp$intensity  
#'   
#'   Spectra::intensity(sp)
#'   Spectra::mz(sp)
#'   as_tibble(sp)
#'   dplyr::as_tibble(sp)
#'   Spectra::name
#'   
#'   Spectra::peaksVariables(sp)
#'   Spectra::peaksData(sp)[2:4]
#'   
#'   Spectra::spectraData(sp[1:10], columns = c("msLevel", "rtime", "name"))
#'   dataOrigin(sp) %>% head
#'   dataStorage(sp)[1:3]
#'   
#'   
#' }
#' 
#' vignevignevignette("MsBackendMsp", package = "MsBackendMsp")
#' 
#' 
#' system.time(
#' p1 <- MetFamily:::parseMSP(
#'   fileSpectra,
#'   minimumIntensityOfMaximalMS2peak = 1000, 
#'   minimumProportionOfMS2peaks = 0.01, 
#'   neutralLossesPrecursorToFragments = T, 
#'   neutralLossesFragmentsToFragments = F, 
#'   progress = F
#' ))
#' 
#' p1 %>% length
#' 
#' p1$spectraList[[1]]
#' 
#' p1 %>% length
#' 
#' p1$spectraList[[1]]
#' 
#' 
#' identical(p1, p2)
#' 
#' p1 %>% class
#' p2 %>% class
#' 
#' snowparam <- SnowParam(workers = 5, type = "SOCK")
#' register(snowparam, default = TRUE)
#' registered()