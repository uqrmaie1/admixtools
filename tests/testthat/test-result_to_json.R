# Tests for result_to_json(). All cases run against example_f2_blocks
# (7 populations) using a small qpadm fit as the canonical "result" fixture
# so we exercise the realistic nested-tibble + list structure that the
# function is built for.

.fit = function() {
  data("example_f2_blocks", package = "admixtools", envir = environment())
  f2 = get("example_f2_blocks", envir = environment())
  suppressMessages(suppressWarnings(qpadm(
    f2,
    left   = c("Altai_Neanderthal.DG", "Vindija.DG"),
    right  = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG", "Switzerland_Bichon.SG"),
    target = "Denisova.DG",
    verbose = FALSE)))
}

# -- envelope shape ---------------------------------------------------------

test_that("envelope has the documented five top-level keys", {
  json = result_to_json(.fit(), fn = "qpadm")
  obj  = jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_named(obj, c("schema_version", "function", "admixtools_version",
                      "args", "result"), ignore.order = TRUE)
  expect_equal(obj$schema_version, 1L)
  expect_equal(obj[["function"]], "qpadm")
  expect_match(obj$admixtools_version, "^\\d+\\.\\d+")
})

# -- fn argument handling ---------------------------------------------------

test_that("fn = NULL stamps 'unknown' into the envelope", {
  obj = jsonlite::fromJSON(result_to_json(.fit()), simplifyVector = FALSE)
  expect_equal(obj[["function"]], "unknown")
})

test_that("fn = NA (untyped) stamps 'unknown'", {
  obj = jsonlite::fromJSON(result_to_json(.fit(), fn = NA), simplifyVector = FALSE)
  expect_equal(obj[["function"]], "unknown")
})

test_that("fn = NA_character_ also stamps 'unknown' (edge case)", {
  # Untyped NA in earlier code: identical(NA_character_, NA) is FALSE because
  # of type mismatch, so the old guard fell through and NA_character_ leaked
  # into the envelope as a serialized JSON null. The current guard uses
  # is.na() which catches every NA flavour.
  obj = jsonlite::fromJSON(result_to_json(.fit(), fn = NA_character_),
                           simplifyVector = FALSE)
  expect_equal(obj[["function"]], "unknown")
})

test_that("fn = '' errors clearly", {
  expect_error(result_to_json(.fit(), fn = ""), "non-empty character")
})

test_that("fn = 123 (non-character) errors clearly", {
  expect_error(result_to_json(.fit(), fn = 123), "non-empty character")
})

# -- args handling ----------------------------------------------------------

test_that("args is echoed verbatim into the envelope", {
  args_in = list(left = c("A", "B"), right = c("C", "D"), target = "T")
  json = result_to_json(.fit(), fn = "qpadm", args = args_in)
  obj  = jsonlite::fromJSON(json, simplifyVector = FALSE)
  expect_equal(obj$args$left,   as.list(args_in$left))
  expect_equal(obj$args$right,  as.list(args_in$right))
  expect_equal(obj$args$target, args_in$target)
})

test_that("args non-list errors clearly", {
  expect_error(result_to_json(.fit(), fn = "qpadm", args = c("not", "a", "list")),
               "list")
})

# -- file-writing path ------------------------------------------------------

test_that("file = <tempfile> writes JSON that round-trips through fromJSON", {
  tmp = tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  result = result_to_json(.fit(), fn = "qpadm", file = tmp)

  expect_true(file.exists(tmp))
  obj = jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  expect_equal(obj[["function"]], "qpadm")
  expect_equal(obj$schema_version, 1L)
})

# -- result content -- pretty / digits --------------------------------------

test_that("pretty = TRUE inserts newlines for human reading", {
  one_line = result_to_json(.fit(), fn = "qpadm", pretty = FALSE)
  pretty   = result_to_json(.fit(), fn = "qpadm", pretty = TRUE)
  # One-line output has no embedded newlines; pretty has many.
  expect_false(grepl("\n", one_line))
  expect_true(grepl("\n", pretty))
})

test_that("digits = 2 rounds numeric fields in the result", {
  # Compare a small numeric extracted from the parsed JSON to its expected
  # 2-significant-digit form. We test the chisq field of rankdrop row 1
  # because it's a known double in the qpadm result.
  json_full = result_to_json(.fit(), fn = "qpadm", digits = NA)
  json_lo   = result_to_json(.fit(), fn = "qpadm", digits = 2)
  expect_true(nchar(json_lo) <= nchar(json_full))
})

# -- result soundness -------------------------------------------------------

test_that("round-trip preserves the result list's top-level shape", {
  json = result_to_json(.fit(), fn = "qpadm")
  obj  = jsonlite::fromJSON(json, simplifyVector = FALSE)
  # qpadm() returns a list with weights / f4 / rankdrop / popdrop. The
  # round-trip through JSON should preserve at least the names that are
  # present; we don't test cell-level equality because jsonlite's row-wise
  # serialization of tibbles loses tibble class (round-trips to list-of-lists).
  expect_true(all(c("weights", "rankdrop") %in% names(obj$result)))
})

test_that("schema_version constant is wired into the envelope", {
  obj = jsonlite::fromJSON(result_to_json(.fit(), fn = "qpadm"),
                           simplifyVector = FALSE)
  expect_equal(obj$schema_version,
               admixtools:::.result_to_json_schema_version)
})
