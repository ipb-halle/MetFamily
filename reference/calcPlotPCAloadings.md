# PCA plot of loadings

Show feature loadings on a PCA plot.

## Usage

``` r
calcPlotPCAloadings(
  pcaObj,
  dataList,
  filterVec = NULL,
  pcaDimensionOne = 1,
  pcaDimensionTwo = 2,
  selectionFragmentPcaLoadingSet = NULL,
  selectionAnalysisPcaLoadingSet = NULL,
  selectionSearchPcaLoadingSet = NULL,
  xInterval = NULL,
  yInterval = NULL,
  loadingsLabels = "None",
  showLoadingsAbundance = FALSE,
  showLoadingsFeaturesAnnotated = TRUE,
  showLoadingsFeaturesUnannotated = TRUE,
  showLoadingsFeaturesSelected = TRUE,
  showLoadingsFeaturesUnselected = TRUE
)
```

## Arguments

- pcaObj:

  output from function calculatePCA

- dataList:

  dataList object

- filterVec:

  deprecated, use filterObj element in calculatePCA() instead

- pcaDimensionOne:

  numeric PCA dimension for x axis

- pcaDimensionTwo:

  numeric PCA dimension for y axis

- selectionFragmentPcaLoadingSet:

  numeric vector of indices

- selectionAnalysisPcaLoadingSet:

  numeric vector of indices

- selectionSearchPcaLoadingSet:

  numeric vector of indices

- xInterval:

  numeric Plot parameter

- yInterval:

  numeric Plot parameter

- loadingsLabels:

  boolean One of "None", "m/z / RT", "Metabolite name", or "Metabolite
  family"

- showLoadingsAbundance:

  logical

- showLoadingsFeaturesAnnotated:

  logical

- showLoadingsFeaturesUnannotated:

  logical

- showLoadingsFeaturesSelected:

  logical

- showLoadingsFeaturesUnselected:

  logical

## Value

a resultList object with annotation colors
