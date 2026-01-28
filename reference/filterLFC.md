# Filter dataset based on LFC of two groups

Filter dataset based on LFC of two groups

## Usage

``` r
filterLFC(dataList, filter_lfc, sampleClasses = dataList$sampleClasses)
```

## Arguments

- dataList:

  list object

- filter_lfc:

  numeric e.g. 2 or -2

- sampleClasses:

  vector of names for group1 and group2

## Value

boolean vector
