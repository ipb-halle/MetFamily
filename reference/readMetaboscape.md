# Read Metaboscape Output File into a QFeatures Object

This function reads a metabolite profile output file (.csv) from
Metaboscape and converts it into a QFeatures object.

## Usage

``` r
readMetaboscape(file, version)
```

## Arguments

- file:

  A character string specifying the path to the Metaboscape output file
  (csv format).

- version:

  A character string specifying the version of Metaboscape used to
  generate the file. This parameter is currently not used.

## Value

A QFeatures object containing:

- An assay named "exampleAssay" with the metabolite counts.

- Row data (feature metadata) extracted from the input file.

- Column data (sample metadata) extracted from the sample names,
  including injection order and sample name.

## Details

At the moment, sample classes are not considered.

## References

\#TODO: Bruker metaboscape site \#TODO: Ordering ?

## See also

[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
for more information on the QFeatures class.
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
for details on the underlying data structure.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have a Metaboscape output file named "data.csv":
qf <- readMetaboscape("data.csv") #TODO: System file 


# Examine the structure of the resulting QFeatures object
qf

# Access the assay data
assay(qf[["exampleAssay"]])

# Access the row data (feature metadata)
rowData(qf[["exampleAssay"]])

# Access the column data (sample metadata)
colData(qf)
} # }

```
