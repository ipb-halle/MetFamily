# Convert data.frame columns to numeric

The data.numericmatrix() function works similar to base::data.matrix()
before R-4.0.0 converting character columns to numeric without
converting to factor first, thus returning the actual numeric values.

## Usage

``` r
data.numericmatrix(x)
```

## Arguments

- x:

  The data.frame to convert

## Value

A matrix with all columns converted to numeric

## Examples

``` r
data.numericmatrix(data.frame(a = c("1", "2", "3"), 
                              b = c("4", "5", "6")))
#>      a b
#> [1,] 1 4
#> [2,] 2 5
#> [3,] 3 6
```
