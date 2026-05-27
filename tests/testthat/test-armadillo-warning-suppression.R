# Regression test for #20: low-level RcppArmadillo "solve(): system is
# singular; attempting approx solution" stderr warnings flood batch fits at
# higher K. Those warnings come from arma::arma_cerr (C++ stderr) and bypass
# R's warning condition system, so withCallingHandlers can't catch them; the
# only reliable check is capture.output(type = "message").
#
# The fix sets ARMA_WARN_LEVEL 1 in cpp_qpadm.cpp, cpp_qpgraph.cpp, and
# cpp_fstats.cpp -- errors still emit, warnings are dropped. The R-level
# near-singular diagnostic (warning + f4_var_rcond field) is unaffected and
# remains the canonical user-facing signal.
#
# Acceptance criteria from the issue:
#   1. A near-singular fit produces zero "system is singular" stderr lines.
#   2. The R-level diagnostic still fires for genuinely-near-singular fits.
# Both are checked below.

test_that("near-singular qpadm fit emits zero RcppArmadillo stderr warnings", {
  # Clone two populations to force a near-singular f4 variance matrix.
  # .clone_pop_in_f2_blocks() lives in helper-fixtures.R.
  data("example_f2_blocks", package = "admixtools", envir = environment())
  ef  = get("example_f2_blocks", envir = environment())
  ef2 = .clone_pop_in_f2_blocks(ef,  "Russia_Ust_Ishim.DG", "Ushim_clone1")
  ef3 = .clone_pop_in_f2_blocks(ef2, "Mbuti.DG",            "Mbuti_clone")

  # Capture stderr explicitly: capture.output(type = "message") is the only
  # way to see the C++-level arma_cerr writes. We also swallow R warnings to
  # avoid testthat third-edition's strict warning-is-error mode tripping on
  # the legitimate R-level near-singular warning that this test explicitly
  # asserts fires.
  stderr_lines = capture.output(
    suppressWarnings(
      res <- qpadm(
        ef3,
        left   = c("Russia_Ust_Ishim.DG", "Ushim_clone1", "Vindija.DG"),
        right  = c("Chimp.REF", "Altai_Neanderthal.DG", "Mbuti.DG", "Mbuti_clone"),
        target = "Switzerland_Bichon.SG",
        fudge  = 1e-14,
        verbose = FALSE)),
    type = "message")

  # The fix: zero "system is singular" lines from Armadillo.
  expect_false(any(grepl("system is singular", stderr_lines, fixed = TRUE)),
               info = paste0("Expected zero Armadillo singular-solve warnings, ",
                             "found ", sum(grepl("system is singular", stderr_lines,
                                                 fixed = TRUE)),
                             ". Lines: ",
                             paste(head(stderr_lines, 3), collapse = " | ")))
  expect_false(any(grepl("solve():", stderr_lines, fixed = TRUE)))
  expect_false(any(grepl("attempting approx solution", stderr_lines, fixed = TRUE)))

  # Fit must still complete successfully -- silencing the warning must not
  # break the fit itself.
  expect_true(is.list(res))
  expect_true("f4_var_rcond" %in% names(res))
})

test_that("near-singular fit still surfaces the R-level diagnostic", {
  # The whole point of the fix is to remove low-level noise without removing
  # the user-facing signal. The R-level signal is two-pronged: the
  # f4_var_rcond field (always surfaced; the programmatic signal) and a
  # human-facing `warning()` (verbose = TRUE only; the interactive signal).
  # Both are unaffected by the C++-level warning suppression.
  data("example_f2_blocks", package = "admixtools", envir = environment())
  ef  = get("example_f2_blocks", envir = environment())
  ef2 = .clone_pop_in_f2_blocks(ef,  "Russia_Ust_Ishim.DG", "Ushim_clone1")
  ef3 = .clone_pop_in_f2_blocks(ef2, "Mbuti.DG",            "Mbuti_clone")

  warns = character(0)
  # verbose = TRUE because the human-facing warning is verbose-gated in
  # qpadm.R (see the `&& verbose` clause around the rcond_concern branch).
  # The cli::cli_inform progress messages would otherwise also fire and
  # need swallowing; suppressMessages handles those.
  res = suppressMessages(withCallingHandlers(
    qpadm(
      ef3,
      left   = c("Russia_Ust_Ishim.DG", "Ushim_clone1", "Vindija.DG"),
      right  = c("Chimp.REF", "Altai_Neanderthal.DG", "Mbuti.DG", "Mbuti_clone"),
      target = "Switzerland_Bichon.SG",
      fudge  = 1e-14,
      verbose = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }))

  # R-level near-singular warning fires.
  expect_true(any(grepl("near-singular", warns, fixed = TRUE)),
              info = paste("warnings captured:",
                           paste(warns, collapse = " | ")))
  # f4_var_rcond is computed and surfaces a tiny number (the programmatic
  # signal, available regardless of verbose).
  expect_true(is.finite(res$f4_var_rcond))
  expect_lt(res$f4_var_rcond, 1e-8)
  # Loadings table is populated for downstream "which right pop is the
  # offender" introspection. Assert the specific structure AND content: the
  # right-set's clone pair (Mbuti.DG / Mbuti_clone) is the source of the
  # singularity, so those two pops must dominate the top two rows of the
  # loadings table. A wrong-axis bug in the SVD attribution (e.g. left vs
  # right pops swapped) would still leave the tibble non-NULL but would
  # rank the unrelated outgroup pops as the offenders instead.
  ld = res$f4_var_singular_loadings
  expect_s3_class(ld, "data.frame")
  expect_named(ld, c("right", "loading"), ignore.order = TRUE)
  top2 = ld$right[order(abs(ld$loading), decreasing = TRUE)][1:2]
  expect_setequal(top2, c("Mbuti.DG", "Mbuti_clone"))
})

test_that("clean fit emits zero stderr noise (sanity)", {
  # A well-conditioned fit on example_f2_blocks should produce no stderr
  # output at all. This catches a regression where ARMA_WARN_LEVEL gets
  # raised back to 2 (Armadillo's default) and ANY Armadillo warning creeps
  # back in -- not just `solve(): system is singular` but also `inv()`,
  # `chol()`, `svd()`, etc. We assert zero captured lines outright rather
  # than substring-match a fixed set of known warning strings.
  data("example_f2_blocks", package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())

  stderr_lines = capture.output(
    suppressMessages(suppressWarnings(
      res <- qpadm(
        ef,
        left   = c("Altai_Neanderthal.DG", "Vindija.DG"),
        right  = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG",
                   "Switzerland_Bichon.SG"),
        target = "Denisova.DG",
        verbose = FALSE))),
    type = "message")

  expect_identical(stderr_lines, character(0),
                   info = paste0("Expected empty stderr on a clean fit, got: ",
                                 paste(head(stderr_lines, 5), collapse = " | ")))
  expect_true(is.list(res))
})
