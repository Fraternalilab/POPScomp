#! /usr/bin/R

#===============================================================================
# POPSR package
# Auxiliary function to call Shiny app
# (C) 2019 Jens Kleinjung and Franca Fraternali
#===============================================================================

## https://stackoverflow.com/questions/37830819/developing-shiny-app-as-a-package-and-deploying-it-to-shiny-server#37833866
## To make your new package’s shiny app runnable from the shiny server
##   create a file /srv/shiny-server/myapp/app.R with the following code
##   and your app will be available at http://<your_shiny_server_url>/myapp :
#dir <- system.file("shiny-examples", "myapp", package = "mypackage")
#setwd(dir)
#shiny::shinyAppDir(".")

#' @export
runPOPSR = function() {
  appDir = system.file("popsr", package = "POPSR")
  if (appDir == "") {
    stop("Could not find popsr. Try re-installing `POPSR`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

#===============================================================================

