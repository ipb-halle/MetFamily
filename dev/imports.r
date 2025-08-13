
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

use_import_from("SummarizedExperiment", c("colData", "rowData","rowData<-", "assay"))
use_import_from("plotly", c("plot_ly", "add_trace", "layout"))
use_import_from("stringr", c("str_squish", "str_split", "str_trim"))
use_import_from("QFeatures", c("QFeatures", "ncols"))

use_import_from("S4Vectors", "DataFrame")
use_import_from("openxlsx2", "read_xlsx")
use_import_from("Matrix", "sparseMatrix")

use_import_from("purrr", "set_names")

use_import_from("future", "plan")
use_import_from("promises", "future_promise")
# -------------------------------------------------------------------------

document()
check(vignettes = F)


# explore dependencies ----------------------------------------------------

pak::local_deps_tree("MetFamily")
