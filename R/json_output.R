# Schema version baked into the JSON envelope's `schema_version` field. Bump
# this constant when the envelope's structure changes in a way consumers must
# adapt to (e.g. renaming a top-level key, changing the dataframe orientation,
# or removing a field). Decoupled from `utils::packageVersion("admixtools")`
# so pure bug-fix patch releases don't force consumers to re-version their
# parsers when the envelope hasn't actually changed.
.result_to_json_schema_version = 1L


#' Serialize an admixtools result to JSON
#'
#' Wraps the output of one of the admixtools fitters
#' ([qpadm()], [qpwave()], [qpgraph()], [qpfstats()], [f4()], [f3()], [f2()],
#' [qpdstat()], [qpadm_sweep()], ...) in a schema-versioned envelope and
#' returns a JSON string (or writes to a file / connection).
#'
#' Non-R orchestrators currently call admixtools from Rscript and do their own
#' JSON marshaling per result type: five or six bespoke serializer paths, one
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
#' return of [qpfstats()]) serialize to nested JSON arrays without dimnames;
#' if you need named pop axes, convert to a long tibble first via
#' [cubelyr::as.tbl_cube()] or manual unpacking.
#'
#' @export
#' @param result The result list (or single data frame, or array) returned by
#'   an admixtools fitter.
#' @param fn Character scalar naming the function that produced `result`
#'   (e.g. `"qpadm"`, `"qpgraph"`). Stamped into the JSON envelope so consumers
#'   can dispatch on the shape. Defaults to `NULL`, which renders as
#'   `"unknown"`: useful for ad-hoc serialization where the caller doesn't
#'   want to plumb the producer name through. `NA` (any flavour: untyped,
#'   `NA_character_`, etc.) is treated the same as `NULL`.
#' @param args Optional named list of the arguments used in the call. Echoed
#'   verbatim into the envelope's `args` field. Use `list()` (the default) to
#'   skip. **Values must be primitives** (scalars, atomic vectors, or nested
#'   lists thereof): non-primitives like functions, environments, or external
#'   pointers don't have well-defined JSON representations and may serialize
#'   to `{}` or error from `jsonlite`.
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
result_to_json = function(result, fn = NULL, args = list(), file = "",
                          pretty = FALSE, digits = NA) {
  # NULL and any NA flavour (untyped NA, NA_character_, NA_integer_, ...) all
  # map to "unknown". `is.na()` on length-1 inputs returns TRUE for every NA
  # type; we guard the length first because is.na() on a NULL returns logical(0)
  # which would short-circuit the || incorrectly otherwise.
  if(is.null(fn) || (length(fn) == 1 && is.na(fn))) {
    fn = "unknown"
  }
  if(!is.character(fn) || length(fn) != 1 || !nzchar(fn))
    stop("'fn' must be a non-empty character scalar (e.g. 'qpadm') or NULL")
  if(!is.list(args)) stop("'args' must be a (possibly empty) list")

  payload = list(
    schema_version     = .result_to_json_schema_version,
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


#' Read the cache metadata sidecar written by extract_f2()
#'
#' Reads the `cache_metadata.json` file written by [extract_f2()] into the
#' f2 output directory and returns it as a named R list. This gives wrappers
#' and orchestrators a stable API for post-build telemetry (`n_snps`,
#' `n_blocks`, `pops`, ...) without reaching into the internal `*.rds` layout.
#'
#' The exact set of fields depends on which code path produced the cache:
#'
#' - **Regular path** (`qpfstats = FALSE`): `schema_version`,
#'   `admixtools_version`, `built_at`, `pops`, `n_snps`, `n_blocks`,
#'   `blgsize`, `maxmiss`, `minmaf`, `maxmaf`, `minac2`, `auto_only`,
#'   `transitions`, `transversions`, `adjust_pseudohaploid`, `poly_only`,
#'   `apply_corr`, `afprod`, `fst`, `qpfstats` (= `FALSE`), `cache_id`.
#' - **qpfstats path** (`qpfstats = TRUE`): `schema_version`,
#'   `admixtools_version`, `built_at`, `pops`, `n_snps`, `n_blocks`,
#'   `blgsize`, `maxmiss`, `adjust_pseudohaploid`, `qpfstats` (= `TRUE`),
#'   `cache_id`. Filter arguments that qpfstats handles internally
#'   (`minmaf`, `maxmaf`, `auto_only`, `transitions`, etc.) are omitted
#'   because they don't apply on this code path.
#'
#' `n_snps` is the count of SNPs that actually contributed to the f2 blocks
#' (i.e. `sum(block_lengths_f2)`), not the post-filter row count of the
#' snpfile. With the default `poly_only = c('f2')`, non-polymorphic SNPs
#' survive filtering but are excluded from f2 blocks, so `n_snps` can be
#' smaller than the visible row count in `snpdat.tsv.gz`. Both the regular
#' and qpfstats paths report `n_snps` with this same definition.
#'
#' `cache_id` is `NULL` (rendered as JSON `null`) if the cache-id
#' computation failed at build time; `built_at` is ISO 8601 with
#' millisecond precision in UTC. The `cache_id` value mirrors the
#' `.f2_cache_id` sidecar contents.
#'
#' `cache_metadata.json` is written via tempfile + `file.rename`, which is
#' atomic on POSIX filesystems: a SIGKILL mid-`extract_f2` either leaves the
#' previous version intact or commits the new one - never a truncated file.
#'
#' @export
#' @param outdir Path to the f2 output directory passed to [extract_f2()].
#' @return A named list with the fields described above.
#' @seealso [extract_f2()], [compute_f2_cache_id()], [result_to_json()]
#' @examples
#' \dontrun{
#' extract_f2("my_geno_prefix", "f2_out/")
#' meta = read_f2_cache_metadata("f2_out/")
#' meta$n_snps      # SNPs that contributed to f2 blocks
#' meta$n_blocks    # number of jackknife blocks
#' meta$cache_id    # matches the .f2_cache_id sidecar
#' }
read_f2_cache_metadata = function(outdir) {
  path = file.path(outdir, 'cache_metadata.json')
  if(!file.exists(path))
    stop(paste0("cache_metadata.json not found in '", outdir,
                "'. Run extract_f2() with an outdir to generate it."))
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}


# Atomically write `content` to `target` via a sibling tempfile + file.rename.
#
# On POSIX, rename(2) is atomic within a filesystem: a SIGKILL mid-write
# either leaves the previous `target` intact or commits the new bytes -
# never a truncated file. This matters for cache_metadata.json: a partial
# file would fail compute_f2_cache_id's mode-1 JSON parse, leaving the
# orchestrator with a "malformed" warning instead of a clean "missing
# cache" signal.
#
# On Windows, file.rename fails when `target` already exists; we fall back
# to file.copy + unlink. That fallback is not atomic, but Windows has no
# portable atomic-rename-replace primitive in base R anyway, and the
# common case (fresh outdir) never exercises the fallback.
.write_atomic = function(content, target) {
  tmp = tempfile(tmpdir = dirname(target), fileext = ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(content, tmp)
  if(!file.rename(tmp, target)) {
    if(!file.copy(tmp, target, overwrite = TRUE))
      stop("failed to commit atomic write to ", target)
  }
  invisible(target)
}
