
#' Launch ADMIXTOOLS 2 GUI
#'
#' @export
run_shiny_admixtools = function() {
  appDir = system.file('shiny_admixtools', package = 'admixtools')
  if (appDir == '') {
    stop('Could not find example directory. Try re-installing `admixtools`.', call. = FALSE)
  }
  # Keep in sync with library() calls in inst/shiny_admixtools/app.R.
  gui_pkgs = c('DT', 'RColorBrewer', 'htmlwidgets', 'igraph', 'magrittr',
               'plotly', 'shiny', 'shinyBS', 'shinyFiles', 'shinyWidgets',
               'shinyalert', 'shinydashboard', 'shinyjs', 'shinythemes',
               'tidyverse')
  missing = gui_pkgs[!vapply(gui_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install_call = sprintf('install.packages(c(%s))',
                           paste0('"', missing, '"', collapse = ', '))
    msg = sprintf('The ADMIXTOOLS 2 GUI needs these packages, which are not installed: %s',
                  paste(missing, collapse = ', '))
    if (interactive()) {
      message(msg)
      ans = utils::menu(c('Yes', 'No'), title = 'Install them now?')
      if (ans == 1) {
        utils::install.packages(missing)
        still_missing = missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
        if (length(still_missing) > 0) {
          stop(sprintf('Install did not satisfy: %s', paste(still_missing, collapse = ', ')),
               call. = FALSE)
        }
      } else {
        stop(sprintf('%s\nInstall them with:\n  %s', msg, install_call), call. = FALSE)
      }
    } else {
      stop(sprintf('%s\nInstall them with:\n  %s', msg, install_call), call. = FALSE)
    }
  }
  shiny::enableBookmarking('server')
  shiny::runApp(appDir, display.mode = 'normal', launch.browser = TRUE)
}




