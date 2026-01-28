# Read MS-DIAL Output File into a QFeatures Object

This function reads the output file from MS-DIAL and converts it into a
QFeatures object.

## Usage

``` r
readMSDial(file, version)
```

## Arguments

- file:

  A string with the path to the MS-DIAL output file.

- version:

  A character string specifying the version of MS-DIAL used to generate
  the file. This parameter is currently not used.

## Value

A QFeatures object containing:

- An assay named "exampleAssay" with the metabolite counts.

- Row data (feature metadata) extracted from the input file.

- Column data (sample metadata) extracted from the input file.

## See also

[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
for more information on the QFeatures class.
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
for details on the underlying data structure.

## Examples
