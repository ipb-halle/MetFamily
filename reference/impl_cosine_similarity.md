# Implement Cosine Similarity

Implement Cosine Similarity

## Usage

``` r
impl_cosine_similarity(
  dataList,
  filter,
  progress = FALSE,
  allow_shift = FALSE,
  nl = FALSE,
  num_filter = 5,
  rm_precursor = T
)
```

## Arguments

- dataList:

  dataset

- filter:

  filter vector

- progress:

  boolena

- allow_shift:

  boolean, modified cosine when true

- nl:

  boolean, include neutral losses when true

- num_filter:

  integer, minimal number of fragments required to perform comparison

- rm_precursor:

  boolean, remove precursor ions from ms2 spectra

## Value

distance matrix
