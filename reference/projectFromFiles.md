# Create a project dataList

Wrapper for a clunky workflow

## Usage

``` r
projectFromFiles(
  ms1_path,
  ms2_path,
  annot_path = NULL,
  siriusFileColumnName = c("ClassyFire superclass"),
  parameterSet = NULL
)
```

## Arguments

- ms1_path:

  file for MS1 intensity .txt file

- ms2_path:

  file for MS2 .msp file

- annot_path:

  file for annotations (.csv or .tsv)

- siriusFileColumnName:

  One of "NPC class", "NPC superclass", "NPC pathway", "ClassyFire
  subclass", "ClassyFire class", "ClassyFire superclass"

- parameterSet:

  list of parameters to use

## Value

dataList
