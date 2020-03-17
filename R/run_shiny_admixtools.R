#' @export
run_shiny_admixtools = function() {
  appDir = system.file('shiny_admixtools', package = 'admixtools')
  if (appDir == '') {
    stop('Could not find example directory. Try re-installing `admixtools`.', call. = FALSE)
  }
  shiny::enableBookmarking('server')
  shiny::runApp(appDir, display.mode = 'normal')
}




