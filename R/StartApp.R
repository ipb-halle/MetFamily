
##' Start the web app in browser
##'
##' @title Start the web app in browser
##' @export
##' @author Hendrik Treutler
##' @examples
##' startMetFamily()
startMetFamily <- function(){
  shiny::runApp(appDir = "R", launch.browser = TRUE)
}