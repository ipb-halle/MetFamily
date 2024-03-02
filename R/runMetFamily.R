
##' Run the web app in browser
##'
##' @title Start the web app in browser
##' @export
##' @author Hendrik Treutler, Steffen Neumann
##' @examples
##' \dontrun{
##' runMetFamily()
##' }
runMetFamily <- function(){
  shiny::runApp(appDir = system.file("MetFamily", package = "MetFamily"), 
  launch.browser = TRUE)
}

startMetFamily <- runMetFamily
