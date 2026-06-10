# Tier 1 adversarial corpus for the read path (P-6 .. P-12).
# P series.
#
# Each malformed .legofit / .bootci input must yield a typed condition or a
# graceful, documented result. Never a bare stop(), never a crash, and never a
# malformed numeric token silently becoming a plausible-but-wrong number.
# Inputs are generated at test time (writeBin handles the CRLF and BOM cases
# cleanly) rather than committed as binary-ish fixtures.

# A minimal but valid .legofit (Initial + DiffEv + Fitted), used as the base
# that each case mutates.
valid_legofit <- c(
  "Initial parameter values",
  "Fixed:",
  "       T_z = 0",
  "       one = 1",
  "Free:",
  "     T_xyz = 2",
  "      T_xy = 0.5",
  "DiffEv reached_goal. cost=1e-16 spread=3e-07",
  "Fitted parameter values",
  "Free:",
  "     T_xyz = 2",
  "      T_xy = 0.5"
)

# Write `lines` to a tempfile with a chosen line ending and optional UTF-8 BOM.
write_legofit <- function(lines, eol = "\n", bom = FALSE) {
  f <- tempfile(fileext = ".legofit")
  con <- file(f, "wb")
  on.exit(close(con))
  if (bom) writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(charToRaw(paste0(paste(lines, collapse = eol), eol)), con)
  f
}

# P-6 -------------------------------------------------------------------------
test_that("P-6: truncated .legofit parses gracefully and informs on missing convergence", {
  truncated <- valid_legofit[1:7]   # Initial block only; no DiffEv, no Fitted
  f <- write_legofit(truncated)
  expect_message(
    res <- read_legofit_output(f),
    class = "legofit_fit_status_unknown"
  )
  expect_s3_class(res, "tbl_df")          # graceful, not a crash
  expect_gt(nrow(res), 0)                 # available params still surfaced
})

# P-7 -------------------------------------------------------------------------
test_that("P-7: duplicated Fitted header uses the first block", {
  dup <- c(valid_legofit,
           "Fitted parameter values", "Free:",
           "     T_xyz = 99", "      T_xy = 88")
  res <- suppressMessages(read_legofit_output(write_legofit(dup)))
  v <- setNames(res$value, res$name)
  expect_equal(unname(v[["T_xyz"]]), 2)   # first block, not the 99 second block
  expect_equal(unname(v[["T_xy"]]), 0.5)
})

# P-8 -------------------------------------------------------------------------
test_that("P-8: nan and inf are surfaced, not coerced to a wrong finite number", {
  nanf <- valid_legofit; nanf[12] <- "      T_xy = nan"
  r1 <- suppressMessages(read_legofit_output(write_legofit(nanf)))
  expect_true(is.nan(setNames(r1$value, r1$name)[["T_xy"]]))

  inff <- valid_legofit; inff[11] <- "     T_xyz = inf"
  r2 <- suppressMessages(read_legofit_output(write_legofit(inff)))
  expect_true(is.infinite(setNames(r2$value, r2$name)[["T_xyz"]]))
})

# P-9 -------------------------------------------------------------------------
test_that("P-9: a negative time is parsed faithfully (surfaced, not silently changed)", {
  neg <- valid_legofit; neg[12] <- "      T_xy = -0.5"
  res <- suppressMessages(read_legofit_output(write_legofit(neg)))
  # The reader's job is to read, not to police physicality; the value must
  # survive verbatim so a downstream validator can decide. It must NOT be
  # silently zeroed or dropped.
  expect_equal(unname(setNames(res$value, res$name)[["T_xy"]]), -0.5)
})

# P-10 ------------------------------------------------------------------------
test_that("P-10: scientific notation parses; a locale comma decimal yields NA under C locale", {
  withr::local_locale(c(LC_NUMERIC = "C"))
  m <- valid_legofit; m[11] <- "     T_xyz = 2.0e0"; m[12] <- "      T_xy = 0,5"
  res <- suppressMessages(read_legofit_output(write_legofit(m)))
  v <- setNames(res$value, res$name)
  expect_equal(unname(v[["T_xyz"]]), 2)        # sci notation OK
  expect_true(is.na(v[["T_xy"]]))              # comma not parsed; surfaced as NA, not a wrong number
})

# P-11 ------------------------------------------------------------------------
test_that("P-11: CRLF line endings and a UTF-8 BOM parse identically to LF", {
  r_lf   <- suppressMessages(read_legofit_output(write_legofit(valid_legofit, eol = "\n")))
  r_crlf <- suppressMessages(read_legofit_output(write_legofit(valid_legofit, eol = "\r\n", bom = TRUE)))
  expect_identical(r_lf$name, r_crlf$name)
  expect_equal(r_lf$value, r_crlf$value)
})

# P-12 ------------------------------------------------------------------------
test_that("P-12: bootci with trailing comment and blank lines recovers the data", {
  b <- c("# bootci.py run", "# confidence: 0.950",
         "       par             est             low            high",
         "      T_xy      0.50000000      0.46250000      0.53750000",
         "     T_xyz      2.00000000      1.90625000      2.09375000",
         "", "# end of pipeline", "")
  res <- parse_bootci_output(b)
  expect_setequal(res$parameter, c("T_xy", "T_xyz"))
  expect_equal(nrow(res), 2L)
  expect_setequal(names(res), c("parameter", "est", "low", "high"))
})
