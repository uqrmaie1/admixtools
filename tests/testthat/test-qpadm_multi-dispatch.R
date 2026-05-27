# Tests for `.qpadm_multi_dispatch_dots()` — the private kwarg validation +
# splitting helper introduced to fix the round-4 finding that `qpadm_sweep`'s
# vignette-promised forwarding of qpadm-only kwargs (e.g., `singular_threshold`,
# `fudge`, `cm_file`) silently crashed before any qpadm() ran.
#
# The helper is package-private (leading dot). Tests access it via `:::`.

# Pin future::plan for the same reason as test-qpadm_sweep.R — these tests
# don't exercise furrr's parallel path, but pinning ensures the rest of the
# test_dir() invocation is sequential too.
if(requireNamespace("future", quietly = TRUE)) future::plan("sequential")


# ── unit tests for .qpadm_multi_dispatch_dots ──────────────────────────────

test_that(".qpadm_multi_dispatch_dots routes geno-only kwargs to $geno", {
  # cm_file is an f4blockdat_from_geno formal, NOT a qpadm formal.
  # The dispatch must put it in $geno only.
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(cm_file = "x.cm"), fit_fun = admixtools::qpadm)
  expect_equal(out$geno, list(cm_file = "x.cm"))
  # Empty fit bucket. Use expect_length rather than expect_equal because
  # R's `[` filtering of a named list with an all-FALSE logical preserves
  # an empty `character(0)` names attribute, while a fresh `list()` is
  # unnamed; expect_equal distinguishes them.
  expect_length(out$fit, 0)
})


test_that(".qpadm_multi_dispatch_dots routes fit-only kwargs to $fit", {
  # singular_threshold, fudge_twice, getcov are qpadm formals, NOT
  # f4blockdat_from_geno formals.
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(singular_threshold = 1e-12, fudge_twice = TRUE, getcov = FALSE),
    fit_fun = admixtools::qpadm)
  expect_equal(out$fit$singular_threshold, 1e-12)
  expect_equal(out$fit$fudge_twice, TRUE)
  expect_equal(out$fit$getcov, FALSE)
  expect_length(out$geno, 0)
})


test_that(".qpadm_multi_dispatch_dots routes shared kwargs to both buckets", {
  # auto_only, blgsize, poly_only are in both f4blockdat_from_geno and qpadm
  # formal lists. The helper does not deduplicate — each callee receives
  # what it needs, and downstream functions are responsible for using the
  # value consistently. (Same kwarg, same value → no semantic conflict.)
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(auto_only = FALSE, blgsize = 0.1, poly_only = TRUE),
    fit_fun = admixtools::qpadm)
  expect_equal(out$geno$auto_only, FALSE)
  expect_equal(out$geno$blgsize, 0.1)
  expect_equal(out$geno$poly_only, TRUE)
  expect_equal(out$fit$auto_only, FALSE)
  expect_equal(out$fit$blgsize, 0.1)
  expect_equal(out$fit$poly_only, TRUE)
})


test_that(".qpadm_multi_dispatch_dots routes correctly for the qpadm_p path", {
  # qpadm_p has a smaller formal set than qpadm. `singular_threshold` is on
  # qpadm but NOT on qpadm_p — passing it with fit_fun=qpadm_p must error
  # with the same Unknown-argument message as any other typo.
  qpadm_p = get("qpadm_p", envir = asNamespace("admixtools"))
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(singular_threshold = 1e-12), fit_fun = qpadm_p),
    "Unknown argument")
  # A kwarg that IS in qpadm_p's formals (`fudge`) routes correctly.
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(fudge = 0.01), fit_fun = qpadm_p)
  expect_equal(out$fit$fudge, 0.01)
})


test_that(".qpadm_multi_dispatch_dots raises a clear error for unknown kwargs", {
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(foo = 1, bar = 2), fit_fun = admixtools::qpadm),
    "Unknown argument\\(s\\) to qpadm_multi\\(\\)")
  # Error message names the offending arg AND points at the valid surface.
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(typo_argname = 1), fit_fun = admixtools::qpadm),
    "typo_argname")
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(typo_argname = 1), fit_fun = admixtools::qpadm),
    "f4blockdat_from_geno")
})


test_that(".qpadm_multi_dispatch_dots tolerates empty dots", {
  # Used heavily in practice (most qpadm_multi calls supply no `...`).
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(), fit_fun = admixtools::qpadm)
  expect_length(out$geno, 0)
  expect_length(out$fit, 0)
})


test_that(".qpadm_multi_dispatch_dots ignores positional (empty-name) kwargs", {
  # Positional `...` is rare in R and qpadm_multi has no semantic for them.
  # The helper passes them through silently (named-only filtering).
  named_dots = list(fudge = 0.01)
  # Simulate positional by stripping the name on one element.
  mixed = c(named_dots, list(42))
  names(mixed)[2] = ""
  out = admixtools:::.qpadm_multi_dispatch_dots(
    mixed, fit_fun = admixtools::qpadm)
  expect_equal(out$fit$fudge, 0.01)
  # Validation skipped for empty-name elements (the 42), no error.
})


# ── integration tests: qpadm_multi end-to-end on f2 path ───────────────────
#
# These cover the rebuild path where the round-3 fix was incomplete: passing
# qpadm-only kwargs through qpadm_multi. Pre-round-4, even with f2 data,
# qpadm_multi(data, ..., fudge_twice=TRUE) would forward fudge_twice via the
# inner future_map's `...`, then qpadm() with f4blocks=.x would route it to
# its `...`, hit the validation guard at qpadm.R:192, and stop().

test_that("qpadm_multi forwards qpadm-only kwargs to qpadm without crash", {
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  models = tibble::tibble(
    left   = list(c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
    right  = list(c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
    target = "Switzerland_Bichon.SG")

  # singular_threshold is a qpadm-only formal (not f4blockdat_from_geno).
  # Set it BELOW the actual rcond (~ 0.005 for this fixture, well above
  # machine eps) so the threshold doesn't trip and the call succeeds.
  # Pre-round-4, the same call would crash with "The following arguments
  # are not used: singular_threshold" because qpadm_multi forwarded `...`
  # into qpadm()'s `...` rather than its named formals, and qpadm()'s
  # validation guard at qpadm.R:192-198 rejects any non-formal kwarg when
  # f4blocks is supplied. This test pins the regression fix.
  res = suppressMessages(suppressWarnings(
    qpadm_multi(ef, models, verbose = FALSE,
                singular_threshold = 1e-12)))   # 1e-12 << 0.005 → no trip
  expect_length(res, 1)
  expect_true(is.list(res[[1]]))
  expect_true("f4_var_rcond" %in% names(res[[1]]))
  expect_true(is.finite(res[[1]]$f4_var_rcond))
})


test_that("qpadm_multi raises informative per-fit error when singular_threshold trips", {
  # When singular_threshold trips inside a per-fit qpadm() call, the
  # tryCatch wrapper in qpadm_multi re-wraps the error with the failing
  # model label so a sweep-style invocation surfaces *which* combo failed
  # rather than bubbling a context-less error from inside furrr.
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  models = tibble::tibble(
    left   = list(c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
    right  = list(c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
    target = "Switzerland_Bichon.SG")
  expect_error(
    suppressMessages(suppressWarnings(
      qpadm_multi(ef, models, verbose = FALSE,
                  singular_threshold = 1.0))),
    "qpadm fit failed for model")
})


test_that("qpadm_multi raises clear error for unknown kwarg via qpadm_sweep path", {
  # qpadm_sweep forwards `...` to qpadm_multi. An unknown kwarg surfaces at
  # qpadm_multi entry with a clear single error, not a cryptic downstream
  # error from f4blockdat_from_geno or qpadm() validation.
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  expect_error(
    suppressWarnings(qpadm_sweep(
      ef,
      targets     = "Switzerland_Bichon.SG",
      source_sets = list(s1 = c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
      right_sets  = list(r1 = c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
      verbose     = FALSE,
      typo_unknown_kwarg = 42)),
    "Unknown argument")
})
