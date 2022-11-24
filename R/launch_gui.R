#'Start Shiny GUI
#'@description start Shiny GUI.
#'@usage launch_gui()
#'@return shiny app
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@examples
#'#launch_gui()
#'@export
launch_gui<-function(){
  shiny::runApp(system.file("shiny", package = "metabox2"),
                display.mode = "normal",
                launch.browser = TRUE)
}
