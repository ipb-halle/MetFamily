# Calculate cosine/modified cosine score

Calculate cosine/modified cosine score

## Usage

``` r
cosine_similarity(
  mz1,
  intensity1,
  precursor_mz1,
  mz2,
  intensity2,
  precursor_mz2,
  fragment_mz_tolerance = 0.01,
  allow_shift = TRUE,
  nl = FALSE,
  normalize = FALSE,
  method = 2
)
```

## Arguments

- mz1:

  numeric vector

- intensity1:

  numeric vector

- precursor_mz1:

  numeric

- mz2:

  numeric vector

- intensity2:

  numeric vector

- precursor_mz2:

  numeric

- fragment_mz_tolerance:

  numeric

- allow_shift:

  boolean

- nl:

  boolean

- normalize:

  boolean

- method:

  1 = normalise by euclidean norm, 2 = sqrt normalised intensity

## Value

score
