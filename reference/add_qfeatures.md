# Add qfeatures to dataList

Add qfeatures to dataList

## Usage

``` r
add_qfeatures(
  dataList,
  qfeatures,
  fileAnnotation = NULL,
  siriusFileColumnName = "ClassyFire superclass"
)
```

## Arguments

- dataList:

  Output from readProjectData.

- qfeatures:

  qfeature object, can be taken from resultObj\$qfeatures.

- fileAnnotation:

  character Path for sirius annotation file.

- siriusFileColumnName:

  One of "NPC class", "NPC superclass", "NPC pathway", "ClassyFire
  subclass", "ClassyFire class", "ClassyFire superclass"

## Value

The dataList object with added sirius annotations.
