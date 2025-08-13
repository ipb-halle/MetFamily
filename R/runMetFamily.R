
##' Run the web app in browser
##'
##' @title Start the web app in browser
##' @param launch.browser Launch web browser to open MetFamily
##' @param ... pass additional parameters down to runApp(), e.g. host="0.0.0.0", port=3838 
##' @export
##' @author Hendrik Treutler, Steffen Neumann
##' @examples
##' \dontrun{
##' runMetFamily()
##' }
runMetFamily <- function(launch.browser = TRUE, ...){
  shiny::runApp(appDir = system.file("MetFamily", package = "MetFamily"), 
  launch.browser = launch.browser, ...)
}

