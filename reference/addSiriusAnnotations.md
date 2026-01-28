# Add Sirius Annotations

Missing values in the Sirius file create empty character values "" in
the annotation, while missing features have NA values. All are
transformed to "" later in `add_qfeatures`.

## Usage

``` r
addSiriusAnnotations(
  qfeatures,
  siriusFile,
  rowData_col = "Alignment ID",
  sirius_col = "featureId",
  siriusFileColumnName = NULL
)
```

## Arguments

- qfeatures:

  dataset

- siriusFile:

  file path

- rowData_col:

  character Qfeatures column to match

- sirius_col:

  character Sirius column to match

- siriusFileColumnName:

  One of "NPC class", "NPC superclass", "NPC pathway", "ClassyFire
  subclass", "ClassyFire class", "ClassyFire superclass". Defaults to
  "ClassyFire superclass" if not specified.

## Value

qfeatures object
