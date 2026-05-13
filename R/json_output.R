#' Serialize an admixtools result to JSON
#'
#' Wraps the output of one of the admixtools fitters
#' ([qpadm()], [qpwave()], [qpgraph()], [qpfstats()], [f4()], [f3()], [f2()],
#' [qpdstat()], [qpadm_sweep()], ...) in a schema-versioned envelope and
#' returns a JSON string (or writes to a file / connection).
#'
#' Non-R orchestrators currently call admixtools from Rscript and do their own
#' JSON marshaling per result type — five or six bespoke serializer paths, one
#' per output shape. `result_to_json()` consolidates that into a single
#' contract: pass the result, the function name, and (optionally) the args you
#' invoked it with, and get back a stable JSON document. Consumers in any
#' language can `json.load()` and treat admixtools like any other CLI tool.
#'
#' The envelope is:
#' ```json
#' {
#'   "schema_version":     1,
#'   "function":           "qpadm",
#'   "admixtools_version": "2.0.10",
#'   "args":               { ... echoed back from the caller ... },
#'   "result":             { ... whatever the fitter returned ... }
#' }
#' ```
#'
#' Tibbles / data.frames are serialized row-wise (each row becomes a JSON
#' object), `NA` becomes `null`, integers/doubles are auto-unboxed, and nested
#' list-columns / lists-of-tibbles are recursed into. 3D arrays (e.g. the raw
#' return of [qpfstats()]) serialize to nested JSON arrays without dimnames —
#' if you need named pop axes, convert to a long tibble first via
#' [cubelyr::as.tbl_cube()] or manual unpacking.
#'
#' @export
#' @param result The result list (or single data frame, or array) returned by
#'   an admixtools fitter.
#' @param fn Character scalar naming the function that produced `result`
#'   (e.g. `"qpadm"`, `"qpgraph"`). Stamped into the JSON envelope so consumers
#'   can dispatch on the shape.
#' @param args Optional named list of the arguments used in the call. Echoed
#'   verbatim into the envelope's `args` field. Use `list()` (the default) to
#'   skip.
#' @param file Path or connection to write to. The default (`""`) returns the
#'   JSON as a character scalar without writing.
#' @param pretty If `TRUE`, indent the JSON for human reading. Default
#'   `FALSE` (one-line; smaller for streaming).
#' @param digits Number of significant digits to use for numeric fields.
#'   Default `NA` preserves full precision (matches `jsonlite::toJSON`'s
#'   sentinel for "do not round").
#' @return If `file = ""`, a character scalar containing the JSON document.
#'   Otherwise the JSON is written to `file` and the path is returned invisibly.
#' @seealso [qpadm()], [qpgraph()], [qpfstats()]
#' @examples
#' left = c("Altai_Neanderthal.DG", "Vindija.DG")
#' right = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG", "Switzerland_Bichon.SG")
#' target = "Denisova.DG"
#' fit = qpadm(example_f2_blocks, left, right, target, verbose = FALSE)
#'
#' # JSON string (one line):
#' json = result_to_json(fit, fn = "qpadm",
#'                       args = list(left = left, right = right, target = target))
#' substr(json, 1, 80)
#'
#' # Pretty-printed to stdout for inspection:
#' cat(result_to_json(fit, fn = "qpadm", pretty = TRUE))
result_to_json = function(result, fn, args = list(), file = "",
                          pretty = FALSE, digits = NA) {
  if(!is.character(fn) || length(fn) != 1 || !nzchar(fn))
    stop("'fn' must be a non-empty character scalar (e.g. 'qpadm')")
  if(!is.list(args)) stop("'args' must be a (possibly empty) list")

  payload = list(
    schema_version     = 1L,
    `function`         = fn,
    admixtools_version = as.character(utils::packageVersion("admixtools")),
    args               = args,
    result             = result)

  json = jsonlite::toJSON(payload, auto_unbox = TRUE, na = "null", null = "null",
                          dataframe = "rows", digits = digits, pretty = pretty)

  if(identical(file, "") || is.null(file)) return(as.character(json))
  writeLines(as.character(json), con = file)
  invisible(file)
}
