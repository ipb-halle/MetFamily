
showErrorDialog <- function(msg){
  showDialog("An error occurred", msg)
}
showInfoDialog <- function(msg){
  showDialog("Information", msg)
}
showDialog <- function(title, msg){
  print("Show dialog")
  #output$infoPopupDialog <- renderUI({
  #  bsModal(id = "modalInfoPopupDialog", title = "Information", trigger = "", size = "large", HTML(msg))
  #})
  #toggleModal(session = session, modalId = "modalInfoPopupDialog", toggle = "open")
  showModal(session = session, ui = modalDialog(title = title, HTML(msg)))
}
showUiDialog <- function(modalDialog){
  print("Show ui dialog")
  showModal(session = session, ui = modalDialog)
}
