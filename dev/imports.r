
usethis::use_devtools()

use_package("shiny")
# use_package("egg") # is it used?

use_pipe()
use_import_from("shinyBS", "bsTooltip")
use_import_from("grDevices", c("as.raster", "rainbow", "rgb", "colorRampPalette"))
use_import_from("methods", "as")
use_import_from("graphics", c("axis", "mtext", "par", "plot.new",
                              "plot.window", "points", "rasterImage", "rect", "segments",
                              "title"))
use_import_from("methods", c("as", "is"))
use_import_from("stats", c("as.dendrogram", "cor", "dendrapply", "dist",
                           "hclust", "is.leaf", "median", "na.omit", "predict", "sd"))
use_import_from("utils", c("flush.console", "read.table", "read.delim"))

use_import_from("SummarizedExperiment", 
                c("colData", "rowData","rowData<-", "assay"))
use_import_from("plotly", 
                c("plot_ly", "add_trace", "layout"))
use_import_from("stringr", 
                c("str_squish", "str_split", "str_trim", "str_equal",
                  "str_extract", "str_detect", "str_remove", "str_to_lower"))
use_import_from("QFeatures", c("QFeatures", "ncols"))

use_import_from("S4Vectors", "DataFrame")
use_import_from("openxlsx2", "read_xlsx")
use_import_from("Matrix", "sparseMatrix")

use_import_from("dplyr", "select")
use_import_from("purrr", c("set_names", "map", "map2", "map_chr", "map_dbl",
                "map_lgl", "map2_chr", "map2_dbl", "map2_lgl"))

use_import_from("future", "plan")
use_import_from("promises", "future_promise")

use_import_from("readr", "read_tsv")
use_import_from("tibble", c("tibble", "lst", "enframe", "deframe"))
use_import_from("tidyr", c("nest", "unnest"))

# -------------------------------------------------------------------------

document()
check(vignettes = F)


# explore dependencies ----------------------------------------------------

pak::local_deps_tree("MetFamily")
