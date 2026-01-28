# Filter MSP entries

Should be used after
[`parseMSP_to_list()`](http://ipb-halle.github.io/MetFamily/reference/parseMSP_to_list.md)
to filter entries and MS2 spectra based on a set of parameters.#'

## Usage

``` r
filter_msp_entries(
  msp_entries,
  minimumIntensityOfMaximalMS2peak = 2000,
  minimumProportionOfMS2peaks = 0.05,
  neutralLossesPrecursorToFragments = T,
  neutralLossesFragmentsToFragments = F,
  mz_tol = 0.17,
  min_frag_nb = 2,
  min_intensity = 500
)
```

## Arguments

- msp_entries:

  list of MSP entries created with
  [`parseMSP_to_list()`](http://ipb-halle.github.io/MetFamily/reference/parseMSP_to_list.md)

- minimumIntensityOfMaximalMS2peak:

  `numeric(1)` Spectra entries with a maximum intensity below this value
  are removed.

- minimumProportionOfMS2peaks:

  `numeric(1)` Between 0 and 1. Within a spectra, fragments with an
  intensity below this proportion of the maximum intensity fragment are
  removed.

- neutralLossesPrecursorToFragments:

  boolean Should neutral losses to precursor be calculated.

- neutralLossesFragmentsToFragments:

  boolean Should neutral losses for each fragment pair be calculated.

- mz_tol:

  `numeric(1)` tolerance used for matching precursor fragment.

- min_frag_nb:

  `numeric(1)` minimum number of fragments needed to include spectrum

- min_intensity:

  `numeric(1)` minimum intensity to include fragment

## Value

list of spectra object, and stats
