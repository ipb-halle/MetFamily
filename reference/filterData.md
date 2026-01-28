# Filter a dataList

Filter a dataList

## Usage

``` r
filterData(
  dataList,
  sampleClasses,
  sampleSet = NULL,
  filterBySamples = FALSE,
  filter_average = NULL,
  filter_lfc = NULL,
  filterList_ms2_masses = NULL,
  filter_ms2_ppm = NULL,
  filter_ms1_masses = NULL,
  filter_ms1_ppm = NULL,
  includeIgnoredPrecursors = FALSE,
  progress = FALSE
)
```

## Arguments

- dataList:

  List object

- sampleClasses:

  sample class

- sampleSet:

  ?

- filterBySamples:

  boolean

- filter_average:

  numeric Features that have a mean intensity below that threshold are
  filtered out.

- filter_lfc:

  numeric Log-fold-change threshold between two treatments. Only
  implement when exactly two groups are present.

- filterList_ms2_masses:

  ?

- filter_ms2_ppm:

  ?

- filter_ms1_masses:

  ?

- filter_ms1_ppm:

  ?

- includeIgnoredPrecursors:

  boolean

- progress:

  boolean

## Value

A filtered dataList object
