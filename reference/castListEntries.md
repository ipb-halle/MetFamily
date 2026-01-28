# Cast logical's and numeric's in a list or data.frame

Tries to cast a list entry (or column in a data.frame) to logical's, if
that does not create any missing values, it is assumed to be a logical
will be replaced by
[`as.logical()`](https://rdrr.io/r/base/logical.html) conversion.
Similarly for numeric entries (or columns). Everything else remains
strings

## Usage

``` r
castListEntries(list)
```

## Arguments

- list:

  list object

## Value

list of the same lenght with logical's and numeric's casted
