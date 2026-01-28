# Export MetFamily project

Creates a `.gz` file to store the MetFamily project.

## Usage

``` r
writeProjectFile(dataList, precursorSet = NULL, file)
```

## Arguments

- dataList:

  main dataset

- precursorSet:

  entries to include, defaults to all

- file:

  filename

## Value

filename

## Details

Implements the logic of prepareAllPrecursors in a generalized manner.
Replaces `createExportMatrix()` as a package function.
