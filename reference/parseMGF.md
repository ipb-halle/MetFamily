# Read MGF files

Reader for Mascot Generic Format MS2 files. Works similar to
[`parseMSP_rewrite()`](http://ipb-halle.github.io/MetFamily/reference/parseMSP_rewrite.md).

## Usage

``` r
parseMGF(
  fileSpectra,
  minimumIntensityOfMaximalMS2peak = 2000,
  minimumProportionOfMS2peaks = 0.05,
  neutralLossesPrecursorToFragments = T,
  neutralLossesFragmentsToFragments = F,
  mz_tol = 0.17,
  min_frag_nb = 2,
  min_intensity = 500,
  progress = FALSE
)
```

## Arguments

- fileSpectra:

  `character(1)` file path

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

- progress:

  boolean, used for interactive app.

## Value

A list including a spectraList element and additional stats.
