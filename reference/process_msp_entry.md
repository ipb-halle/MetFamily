# Filter single MSP entry

nb of peaks and peaks with NL are updated in header peaks are filtered,
and neutral losses added if requested stats section added returns
backwards-compatible spectraItem as used by MetFamily

## Usage

``` r
process_msp_entry(
  msp_entry,
  minimumProportionOfMS2peaks,
  neutralLossesPrecursorToFragments,
  neutralLossesFragmentsToFragments,
  min_intensity = 500
)
```

## Arguments

- msp_entry:

  MSP entry object

- minimumProportionOfMS2peaks:

  `numeric(1)` Between 0 and 1. Within a spectra, fragments with an
  intensity below this proportion of the maximum intensity fragment are
  removed.

- neutralLossesPrecursorToFragments:

  boolean Should neutral losses to precursor be calculated.

- neutralLossesFragmentsToFragments:

  boolean Should neutral losses for each fragment pair be calculated.

- min_intensity:

  `numeric(1)` minimum intensity to include fragment

## Value

filtered spectra
