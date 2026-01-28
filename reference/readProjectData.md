# Read MetFamily Project data saved by the export function

Read in a project described as lines. Main workhorse used when loading a
project file or separate input files

## Usage

``` r
readProjectData(fileLines, progress = FALSE)
```

## Arguments

- fileLines:

  Character vector with content of a project file

- progress:

  Whether to update a shiny Progress bar

## Value

A big dataList.

## See also

[processMS1data](http://ipb-halle.github.io/MetFamily/reference/processMS1data.md)
