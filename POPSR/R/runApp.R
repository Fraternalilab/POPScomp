#! /usr/bin/R

#===============================================================================
# POPSR package
# Auxiliary function to call Shiny app
# (C) 2019 Jens Kleinjung and Franca Fraternali
#===============================================================================

#' @export
runPOPSR = function() {
  appDir = system.file("popsr", package = "POPSR")
  if (appDir == "") {
    stop("Could not find popsr. Try re-installing `POPSR`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

