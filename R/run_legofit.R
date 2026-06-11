#' Run LEGOFIT on a graph or .lgo file and read the result back
#'
#' Thin convenience wrapper around the external `legofit` binary. Writes
#' the graph to a temporary .lgo (if given a graph), invokes `legofit`
#' against a site-pattern file, and parses the fitted output back into an
#' admixtools edge tibble via [read_legofit_output()].
#'
#' Requires the `legofit` executable on `PATH` (or via the `bin`
#' argument). LEGOFIT is an external tool, not an R dependency; this
#' function errors cleanly when it is not found.
#'
#' @param graph_or_lgo An admixtools edge tibble / igraph, or a path to a
#'   .lgo file.
#' @param patterns Path to a LEGOFIT site-pattern (.opf) file.
#' @param bin Path to the `legofit` binary. Default looks on `PATH`.
#' @param args Character vector of extra flags passed to `legofit`.
#' @param graph When `graph_or_lgo` is a .lgo path, the admixtools graph
#'   to align the fitted parameters to (forwarded to read_legofit_output).
#' @return The tibble returned by [read_legofit_output()].
#' @export
run_legofit <- function(graph_or_lgo, patterns, bin = Sys.which("legofit"),
                        args = c("--threads", "1"), graph = NULL) {
  # `--threads` is effectively required by legofit's CLI
  # (usage: legofit [options] --threads <x> <input.lgo> <data.txt>);
  # default to a single deterministic thread. Callers override via `args`.
  if (!nzchar(bin) || !file.exists(bin))
    rlang::abort(
      c("Could not find the `legofit` executable.",
        "i" = "Install LEGOFIT (https://github.com/alanrogers/legofit) and put it on PATH, or pass `bin=`."),
      class = "legofit_binary_not_found")
  if (!file.exists(patterns))
    rlang::abort(sprintf("Site-pattern file not found: %s", patterns),
                 class = "legofit_invalid_input")

  # Resolve the .lgo: write one if handed a graph, else use the path.
  if (is.character(graph_or_lgo) && length(graph_or_lgo) == 1 &&
      file.exists(graph_or_lgo)) {
    lgo_path <- graph_or_lgo
  } else {
    lgo_path <- tempfile(fileext = ".lgo")
    writeLines(graph_to_lgo(graph_or_lgo), lgo_path)
    on.exit(unlink(lgo_path), add = TRUE)
    if (is.null(graph)) graph <- graph_or_lgo
  }

  out_path <- tempfile(fileext = ".legofit")
  err_path <- tempfile(fileext = ".err")
  on.exit(unlink(c(out_path, err_path)), add = TRUE)
  # With stdout/stderr as file paths, system2() returns the exit code as a
  # plain integer (0 = success); there is NO attr(., "status"). Send stderr
  # to its own file so we can surface it on failure.
  code <- system2(bin, c(args, shQuote(lgo_path), shQuote(patterns)),
                  stdout = out_path, stderr = err_path)
  if (!identical(as.integer(code), 0L)) {
    tail_err <- if (file.exists(err_path))
      paste(utils::tail(readLines(err_path, warn = FALSE), 5), collapse = "\n")
    else ""
    rlang::abort(
      c(sprintf("legofit exited with status %s.", code),
        "i" = tail_err),
      class = "legofit_run_failed")
  }

  read_legofit_output(out_path, graph = graph)
}
