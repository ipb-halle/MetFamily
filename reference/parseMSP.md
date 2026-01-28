# Parse MS/MS spectra

**\[superseded\]** Superseded by parseMSP rewrite 2025-08

## Usage

``` r
parseMSP(
  fileSpectra,
  minimumIntensityOfMaximalMS2peak,
  minimumProportionOfMS2peaks,
  neutralLossesPrecursorToFragments,
  neutralLossesFragmentsToFragments,
  progress = FALSE,
  ...
)
```

## Arguments

- fileSpectra:

  filepath

- minimumIntensityOfMaximalMS2peak:

  numeric

- minimumProportionOfMS2peaks:

  numeric

- neutralLossesPrecursorToFragments:

  boolean

- neutralLossesFragmentsToFragments:

  boolean

- progress:

  boolean

- ...:

  ignored. Used for forward compatibility.

## Value

list object

## Details

Read MS2 spectra from MSP file
