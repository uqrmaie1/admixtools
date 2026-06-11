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
#' @details
#' When `graph_or_lgo` is a graph, it is written out with [graph_to_lgo()];
#' extra arguments in `...` are forwarded there. The [graph_to_lgo()] default
#' `time_handling = "fix_admix"` (and `"init"`) cannot place node times for
#' graphs with admixture or missing drift lengths: the resulting NA branch
#' lengths make [graph_to_lgo()]'s topology walk loop forever. `run_legofit()`
#' detects this case up front and aborts with a `legofit_invalid_input`
#' condition pointing you to `time_handling = "free"`, rather than hanging.
#' Pass `time_handling = "free"` through `...` for such graphs.
#'
#' `--threads` is effectively required by `legofit`, so it is always passed
#' (defaulting to `1`) even when you override `args` to add other flags; an
#' explicit `--threads` in `args` wins. The default `args` run a
#' single-threaded stochastic fit. For a deterministic run (e.g. to reproduce
#' a cached result) pass `args = c("--threads", "1", "-1", "-d", "0")`.
#'
#' @param graph_or_lgo An admixtools edge tibble / igraph, or a path to a
#'   .lgo file.
#' @param patterns Path to a LEGOFIT site-pattern (.opf) file.
#' @param bin Path to the `legofit` binary, or a bare command name to look
#'   up on `PATH`. Defaults to `"legofit"` on `PATH`.
#' @param args Character vector of flags passed to `legofit`. `--threads` is
#'   always ensured (see Details).
#' @param graph When `graph_or_lgo` is a .lgo path, the admixtools graph to
#'   align the fitted parameters to. Defaults to the topology parsed from the
#'   .lgo file with [read_lgo()], so a path input returns an edge tibble just
#'   like a graph input.
#' @param ... Additional arguments forwarded to [graph_to_lgo()] when
#'   `graph_or_lgo` is a graph (for example `time_handling = "free"`).
#' @return An edge tibble: the input edges with the LEGOFIT-fitted values
#'   filled in, as returned by [read_legofit_output()] (carrying the
#'   `node_times`, `fit_convergence`, and `identifiability` attributes).
#' @export
run_legofit <- function(graph_or_lgo, patterns, bin = "legofit",
                        args = c("--threads", "1"), graph = NULL, ...) {
  # `--threads` is effectively required by legofit's CLI
  # (usage: legofit [options] --threads <x> <input.lgo> <data.txt>). Ensure it
  # is always present so a caller who overrides `args` to add other flags (e.g.
  # `-1 -d 0`) need not remember to re-supply it. An explicit `--threads` wins.
  if (!"--threads" %in% args) args <- c("--threads", "1", args)

  # Resolve the binary. A `bin` that contains a path separator is treated as
  # an explicit path and used as-is; a bare command name like "legofit" is
  # always looked up on PATH. We key off the separator rather than
  # file.exists() so the guard agrees with how system2() resolves the
  # command: a bare name runs via PATH even if a like-named file happens to
  # sit in the working directory.
  has_sep <- grepl(.Platform$file.sep, bin, fixed = TRUE)
  bin <- if (nzchar(bin) && has_sep) bin else Sys.which(bin)
  if (!nzchar(bin) || !file.exists(bin))
    rlang::abort(
      c("Could not find the `legofit` executable.",
        "i" = "Install LEGOFIT (https://github.com/alanrogers/legofit) and put it on PATH, or pass `bin=`."),
      class = "legofit_binary_not_found")
  if (!file.exists(patterns))
    rlang::abort(sprintf("Site-pattern file not found: %s", patterns),
                 class = "legofit_invalid_input")

  # Resolve the .lgo. A length-1 character argument is meant to be a path;
  # if it doesn't exist, say so plainly rather than failing deep inside
  # graph_to_lgo() trying to coerce a string to a graph.
  if (is.character(graph_or_lgo) && length(graph_or_lgo) == 1) {
    if (!file.exists(graph_or_lgo))
      rlang::abort(sprintf(".lgo file not found: %s", graph_or_lgo),
                   class = "legofit_invalid_input")
    lgo_path <- graph_or_lgo
    # Align the fit to an edge tibble so a .lgo path returns the same shape as
    # a graph input. With graph = NULL, read_legofit_output() would hand back
    # the raw param table instead; parse the .lgo we were given to recover its
    # topology so both entry points return an edge tibble.
    if (is.null(graph)) graph <- read_lgo(lgo_path)
  } else {
    # Fail fast instead of hanging. graph_to_lgo()'s "fix_admix"/"init" time
    # modes feed branch lengths into a bottom-up topology walk that loops
    # forever when any branch length is NA: the NA-fed node never resolves yet
    # the walk keeps reporting progress. That is exactly the admixture /
    # missing-drift case (e.g. the column-less edge tibble from read_lgo()).
    # Predict the branch lengths the walk would see and abort early if any are
    # NA, rather than spinning. (`time_handling = "free"` skips the walk.)
    dots <- list(...)
    th <- if ("time_handling" %in% names(dots)) dots$time_handling else "fix_admix"
    th <- tryCatch(match.arg(th, c("fix_admix", "init", "free")),
                   error = function(cnd) NA_character_)
    if (!is.na(th) && th %in% c("fix_admix", "init")) {
      edg <- suppressWarnings(validate_edge_tibble(coerce_to_edge_tibble(graph_or_lgo)))
      dtt <- if ("drift_to_time" %in% names(dots)) dots$drift_to_time
             else default_drift_to_time
      bl  <- if ("time" %in% names(edg))
               ifelse(!is.na(edg$time), edg$time, dtt(edg$weight, edg$to))
             else dtt(edg$weight, edg$to)
      if (identical(th, "fix_admix")) bl[edg$type == "admix"] <- 0
      if (anyNA(bl))
        rlang::abort(
          c(sprintf(paste("graph_to_lgo(time_handling = \"%s\") cannot place every",
                          "node: some branch lengths are NA (admixture or missing",
                          "drift lengths)."), th),
            "x" = "Running as-is would hang in graph_to_lgo()'s topology walk.",
            "i" = paste("Pass `time_handling = \"free\"` to run_legofit() to let the",
                        "LEGOFIT optimizer fit these times instead.")),
          class = "legofit_invalid_input")
    }
    lgo_path <- tempfile(fileext = ".lgo")
    on.exit(unlink(lgo_path), add = TRUE)          # register before writing
    writeLines(graph_to_lgo(graph_or_lgo, ...), lgo_path)
    if (is.null(graph)) graph <- graph_or_lgo
  }

  out_path <- tempfile(fileext = ".legofit")
  err_path <- tempfile(fileext = ".err")
  on.exit(unlink(c(out_path, err_path)), add = TRUE)
  # With stdout/stderr as file paths, system2() returns the exit code as a
  # plain integer (0 = success); there is NO attr(., "status"). Send stderr
  # to its own file so we can surface it on failure.
  # system2() shQuote()s `command` itself, so pass `bin` raw (quoting it here
  # would double-quote the path); it does NOT quote args, so we quote those.
  code <- system2(bin, c(args, shQuote(lgo_path), shQuote(patterns)),
                  stdout = out_path, stderr = err_path)
  if (!identical(as.integer(code), 0L)) {
    tail_err <- if (file.exists(err_path))
      paste(utils::tail(readLines(err_path, warn = FALSE), 5), collapse = "\n")
    else ""
    msg <- sprintf("legofit exited with status %s.", code)
    rlang::abort(
      if (nzchar(tail_err)) c(msg, "i" = tail_err) else msg,
      class = "legofit_run_failed")
  }

  read_legofit_output(out_path, graph = graph)
}
