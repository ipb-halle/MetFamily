# Create a filterObject for console use

Mixes code from filterData() and doPerformFiltering_impl() to provide a
simple way to create a compatible filterObj in a script.

## Usage

``` r
makeFilterObj(
  dataList,
  filter_average = NULL,
  filter_lfc = NULL,
  sampleClasses = NULL
)
```

## Arguments

- dataList:

  list object

- filter_average:

  numeric

- filter_lfc:

  numeric

- sampleClasses:

  character vector

## Value

filterObject

## Details

doPerformFiltering_impl() relies on a DataList in the global
environment, so it can't be used in a script.
