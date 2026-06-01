# Tests for read_legofit_bootstrap() Phase C

bootci_path  <- function() test_path("fixtures", "ourex1.bootci")
ourex1_graph <- function() make_ourex1_graph()

# ---------------------------------------------------------------------------
# parse_bootci_output (Step 8 helper)
# ---------------------------------------------------------------------------

test_that("parse_bootci_output: returns expected structure from ourex1.bootci", {
  lines  <- readLines(bootci_path(), encoding = "UTF-8")
  result <- parse_bootci_output(lines)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("parameter", "est", "low", "high"))
  expect_equal(nrow(result), 2L)
})

test_that("parse_bootci_output: values match staged fixture verbatim", {
  lines  <- readLines(bootci_path(), encoding = "UTF-8")
  result <- parse_bootci_output(lines)
  # ourex1.bootci has T_xy and T_xyz
  xy  <- result[result$parameter == "T_xy",  ]
  xyz <- result[result$parameter == "T_xyz", ]
  expect_equal(xy$est,   0.5,   tolerance = 1e-6)
  expect_equal(xyz$est,  2.0,   tolerance = 1e-6)
  expect_equal(xy$low,  0.4625, tolerance = 1e-4)
  expect_equal(xyz$low, 1.90625, tolerance = 1e-4)
})

test_that("parse_bootci_output: optional lbl column tolerated", {
  # Synthesize a 5-column bootci.py output (bootci.py -l mylabel)
  lbl_content <- c(
    "# bootci.py run at: 2026-05-20",
    "# input: /tmp/fake.flat",
    "# confidence: 0.950",
    "       par             est             low            high             lbl",
    "      T_xy      0.50000000      0.46250000      0.53750000       main_run",
    "     T_xyz      2.00000000      1.90625000      2.09375000       main_run"
  )
  result <- parse_bootci_output(lbl_content)
  expect_named(result, c("parameter", "est", "low", "high", "lbl"))
  expect_equal(nrow(result), 2L)
  expect_equal(result$lbl, c("main_run", "main_run"))
})

test_that("parse_bootci_output: malformed header aborts legofit_invalid_input", {
  bad_content <- c(
    "# bootci.py run at: 2026-05-20",
    "# confidence: 0.950",
    "       par             low             est            high",  # order swapped
    "      T_xy      0.46250000      0.50000000      0.53750000"
  )
  expect_error(
    parse_bootci_output(bad_content),
    class = "legofit_invalid_input"
  )
})

# ---------------------------------------------------------------------------
# read_legofit_bootstrap — end-to-end
# ---------------------------------------------------------------------------

# Case 10: smoke
test_that("T10: smoke — read_legofit_bootstrap parses ourex1.bootci", {
  expect_no_error(read_legofit_bootstrap(bootci_path()))
})

# Case 11: returns expected columns
test_that("T11: returns expected columns", {
  result <- read_legofit_bootstrap(bootci_path())
  expect_true(all(c("parameter", "family", "point_estimate", "lo", "hi") %in% names(result)))
})

# Case 12: optional lbl column passes through
test_that("T12: lbl column passes through when present", {
  lbl_content <- c(
    "# bootci.py run at: 2026-05-20",
    "# input: /tmp/fake.flat",
    "# confidence: 0.950",
    "       par             est             low            high             lbl",
    "      T_xy      0.50000000      0.46250000      0.53750000       main_run",
    "     T_xyz      2.00000000      1.90625000      2.09375000       main_run"
  )
  tmp <- tempfile(fileext = ".bootci")
  writeLines(lbl_content, tmp)
  on.exit(unlink(tmp))
  result <- read_legofit_bootstrap(tmp)
  expect_true("lbl" %in% names(result))
  expect_equal(result$lbl, c("main_run", "main_run"))
})

# Case 13: graph supplied, ghost-detection structural check
test_that("T13: with graph=ourex1_graph, result has correct families", {
  result <- read_legofit_bootstrap(bootci_path(), graph = ourex1_graph())
  expect_equal(result$family, c("time", "time"))
})

# Case 13a: malformed header aborts
test_that("T13a: malformed bootci header aborts legofit_invalid_input", {
  bad_content <- c(
    "# bootci.py run at: 2026-05-20",
    "       par             low             est            high",
    "      T_xy      0.46250000      0.50000000      0.53750000"
  )
  tmp <- tempfile(fileext = ".bootci")
  writeLines(bad_content, tmp)
  on.exit(unlink(tmp))
  expect_error(
    read_legofit_bootstrap(tmp),
    class = "legofit_invalid_input"
  )
})

# Point estimates match ourex1.legofit (within floating-point tolerance)
test_that("point_estimate matches fitted output within tolerance", {
  result <- read_legofit_bootstrap(bootci_path())
  xy_est  <- result$point_estimate[result$parameter == "T_xy"]
  xyz_est <- result$point_estimate[result$parameter == "T_xyz"]
  expect_equal(xy_est,  0.5, tolerance = 0.001)
  expect_equal(xyz_est, 2.0, tolerance = 0.001)
})

# ---------------------------------------------------------------------------
# Regression tests for code-review fixes
# ---------------------------------------------------------------------------

# Fix 2: seq() descending bug — header-only file returns empty tibble not crash
test_that("Fix2: bootci file with no data rows returns empty tibble", {
  header_only <- c(
    "# bootci.py run at: 2026-05-20",
    "# confidence: 0.950",
    "       par             est             low            high"
  )
  result <- parse_bootci_output(header_only)
  expect_equal(nrow(result), 0L)
  expect_named(result, c("parameter", "est", "low", "high"))
})

# Fix 6: trailing comment line after data rows does not abort
test_that("Fix6: trailing comment line after data rows is silently skipped", {
  content_with_trailer <- c(
    "# bootci.py run at: 2026-05-20",
    "# confidence: 0.950",
    "       par             est             low            high",
    "      T_xy      0.50000000      0.46250000      0.53750000",
    "# end of output"
  )
  result <- parse_bootci_output(content_with_trailer)
  expect_equal(nrow(result), 1L)
  expect_equal(result$parameter, "T_xy")
})

# Fix 7: column order matches documented order parameter, family, point_estimate, lo, hi
test_that("Fix7: read_legofit_bootstrap column order matches docs", {
  result <- read_legofit_bootstrap(bootci_path())
  expect_equal(names(result)[1:5], c("parameter", "family", "point_estimate", "lo", "hi"))
})

# Finding 5: drop_ghost_params column contract — with graph supplied
test_that("Rv5: read_legofit_bootstrap with graph passes without column-name error", {
  # When graph is supplied the rename-around-drop_ghost_params path is exercised;
  # confirm the result still has the correct 'parameter' column (not 'name').
  result <- read_legofit_bootstrap(bootci_path(), graph = ourex1_graph())
  expect_true("parameter" %in% names(result))
  expect_false("name" %in% names(result))
})

# ---------------------------------------------------------------------------
# Master plan Tier 0 gap closure (the malformed-header abort already exists at
# "parse_bootci_output: malformed header aborts"; these add the other two)
# ---------------------------------------------------------------------------

# U-R17a: a file with no non-comment line has no column header at all.
test_that("U-R17a: parse_bootci_output aborts when every line is a comment", {
  expect_error(
    parse_bootci_output(c("# bootci.py run at: 2026-05-20", "# confidence: 0.950")),
    class = "legofit_invalid_input"
  )
})

# U-R17b: a data row whose token count does not match the header (4 or 5) is a
# structural error, distinct from a bad header.
test_that("U-R17b: parse_bootci_output aborts on a data row with the wrong column count", {
  bad <- c(
    "# confidence: 0.950",
    "       par             est             low            high",
    "      T_xy      0.50000000      0.46250000"   # 3 tokens, expected 4
  )
  expect_error(
    parse_bootci_output(bad),
    class = "legofit_invalid_input"
  )
})

# U-R19 (bootstrap half): re-reading the same bootci file is idempotent.
test_that("U-R19: read_legofit_bootstrap is idempotent on re-read", {
  b1 <- read_legofit_bootstrap(bootci_path(), graph = ourex1_graph())
  b2 <- read_legofit_bootstrap(bootci_path(), graph = ourex1_graph())
  expect_identical(b1, b2)
})
