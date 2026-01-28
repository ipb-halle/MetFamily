# Calculate a distance matrix between MS/MS spectra in the MetFamily feature Matrix

Computes a distance matrix between MS/MS spectra in the feature matrix
representation using various distance or similarity measures. The
function supports multiple distance metrics, of which the following are
also available in the MetFamily GUI: "Jaccard", "Jaccard
(intensity-weighted)", "Jaccard (fragment-count-weighted)" and "NDP
(Normalized dot product)". The following list defines the metrics and
their characteristics:

## Usage

``` r
calculateDistanceMatrix(
  dataList,
  filter,
  distanceMeasure = "Jaccard (intensity-weighted)",
  removePrecursor = TRUE,
  minNumberFragments = 5,
  progress = FALSE
)
```

## Arguments

- dataList:

  List object containing precursor, feature and MS/MS data.

- filter:

  Logical or integer vector indicating for which precursors to include
  the MS/MS spectra.

- distanceMeasure:

  Character string specifying the distance metric to use. Supported
  values include "Jaccard", "Manhatten", "NDP (Normalized dot product)",
  and others.

- removePrecursor:

  Logical. Should the precursor ion be removed from the MS2 spectra?
  Defaults to TRUE. Only implemented for cosine-based methods.

- minNumberFragments:

  Numeric. Spectra with fewer fragments (after precursor removal if
  desired) are excluded. Only implemented for cosine-based methods.

- progress:

  Logical or NA. If TRUE, progress is reported via incProgress().

## Value

A list with elements:

- distanceMatrix:

  A numeric matrix of pairwise distances.

- filter:

  The filter vector used, same as the input parameter.

- distanceMeasure:

  The distance metric used, same as the input parameter.

## Details

- "Jaccard": Jaccard distance \\1 - \|A \cap B\| / \|A \cup B\|\\: the
  fraction of matching fragments among all fragments, where A and B are
  MS/MS features of two different precursors.

- "Jaccard (intensity-weighted)": Jaccard distance weighted by the
  intensity of features. First, intensities are discretized:
  \[0.01-0.2\[ are set down to 0.01, \[0.2, 0.4\[ down to 0.2, and \>=
  0.4 increased to 1. For matching fragments the higher intensity is
  used. Then, Jacard becomes the sum of matching intensities among the
  sum of all intensities in A and B. \\1 - \sum\_{i\in matches} max(A_i,
  B_i) / (sum(A) + sum(B))\\

- "Jaccard (fragment-count-weighted)": Jaccard distance weighted by
  relative occurance of the fragments among the precursors after
  filtering counts: \\1 - \sum\_{i\in matches} freq(A_i) / (sum\_{\notin
  matches} freq(A) + sum\_{\notin matches}(B))\\

- "NDP (Normalized dot product)": Normalized dot product similarity:
  \\NDP = \frac{\left( \sum\_{i}^{\text{S1\\S2}} W\_{\text{S1},i}
  W\_{\text{S2},i} \right)^2}{\sum_i W\_{\text{S1},i}^2 \sum_i
  W\_{\text{S2},i}^2}\\ as described in Gaquerel et al. 2015.
  [10.1073/pnas.1610218113](https://www.pnas.org/doi/10.1073/pnas.1610218113#sec-4-5)
