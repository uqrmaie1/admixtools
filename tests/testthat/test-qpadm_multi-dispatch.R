# Tests for `.qpadm_multi_dispatch_dots()` — the private kwarg validation +
# splitting helper that fixes the round-4/5 findings on kwarg-forwarding:
#
#   * unknown names → single clear error at qpadm_multi entry
#   * reserved names (internal/positional) → rejected before reaching do.call
#   * geno-only names on the f2 path → rejected (silent-drop prevention)
#   * everything else → routed to $geno / $fit buckets for do.call
#
# The helper is package-private (leading dot). Tests access it via `:::`.

# Pin future::plan globally for the test_dir() invocation. See the matching
# comment in test-qpadm_sweep.R for the rationale (CI BLAS bit-differences
# under multisession can flake expect_identical on f4_var_rcond).
if(requireNamespace("future", quietly = TRUE)) future::plan("sequential")


# ── unit tests: dispatch routing ────────────────────────────────────────────

test_that(".qpadm_multi_dispatch_dots routes geno-only kwargs to $geno", {
  # cm_file is an f4blockdat_from_geno formal (geno-only). On the geno path
  # the dispatch must put it in $geno only.
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(cm_file = "x.cm"), fit_fun = admixtools::qpadm, on_f2_path = FALSE)
  expect_equal(out$geno, list(cm_file = "x.cm"))
  expect_length(out$fit, 0)
})


test_that(".qpadm_multi_dispatch_dots routes fit-only kwargs to $fit", {
  # singular_threshold, fudge_twice, getcov are qpadm formals, not
  # f4blockdat_from_geno formals.
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(singular_threshold = 1e-12, fudge_twice = TRUE, getcov = FALSE),
    fit_fun = admixtools::qpadm, on_f2_path = FALSE)
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
    fit_fun = admixtools::qpadm, on_f2_path = FALSE)
  expect_equal(out$geno$auto_only, FALSE)
  expect_equal(out$geno$blgsize, 0.1)
  expect_equal(out$geno$poly_only, TRUE)
  expect_equal(out$fit$auto_only, FALSE)
  expect_equal(out$fit$blgsize, 0.1)
  expect_equal(out$fit$poly_only, TRUE)
})


test_that(".qpadm_multi_dispatch_dots routes qpadm_p path correctly", {
  # qpadm_p has a smaller formal set than qpadm. singular_threshold is on
  # qpadm but not on qpadm_p — passing it with fit_fun = qpadm_p must error
  # with the same Unknown-argument message as any other typo.
  qpadm_p = get("qpadm_p", envir = asNamespace("admixtools"))
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(singular_threshold = 1e-12), fit_fun = qpadm_p, on_f2_path = FALSE),
    "Unknown argument")
  # A kwarg that IS in qpadm_p's formals (`fudge`) routes correctly.
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(fudge = 0.01), fit_fun = qpadm_p, on_f2_path = FALSE)
  expect_equal(out$fit$fudge, 0.01)
})


test_that(".qpadm_multi_dispatch_dots tolerates empty dots", {
  out = admixtools:::.qpadm_multi_dispatch_dots(
    list(), fit_fun = admixtools::qpadm, on_f2_path = FALSE)
  expect_length(out$geno, 0)
  expect_length(out$fit, 0)
})


test_that(".qpadm_multi_dispatch_dots ignores positional (empty-name) kwargs", {
  # Positional `...` is rare in R and qpadm_multi has no semantic for them.
  named_dots = list(fudge = 0.01)
  mixed = c(named_dots, list(42))
  names(mixed)[2] = ""
  out = admixtools:::.qpadm_multi_dispatch_dots(
    mixed, fit_fun = admixtools::qpadm, on_f2_path = FALSE)
  expect_equal(out$fit$fudge, 0.01)
})


# ── unit tests: error layers (the round-5 hardening) ────────────────────────

test_that(".qpadm_multi_dispatch_dots rejects unknown kwargs", {
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(foo = 1, bar = 2), fit_fun = admixtools::qpadm, on_f2_path = FALSE),
    "Unknown argument\\(s\\) to qpadm_multi\\(\\)")
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(typo_argname = 1), fit_fun = admixtools::qpadm, on_f2_path = FALSE),
    "typo_argname")
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(typo_argname = 1), fit_fun = admixtools::qpadm, on_f2_path = FALSE),
    "f4blockdat_from_geno")
})


test_that(".qpadm_multi_dispatch_dots rejects reserved names with a pointed message", {
  # The reserved set (data/pref/f2_data/left/right/target/f4blocks/popcombs/
  # verbose) is filled internally by qpadm_multi. Round-4 dispatch routed
  # them silently → do.call duplicate-name errors or silent positional
  # NULL-displacement. Round-5 rejects them at entry with a clear message
  # naming the offending kwarg and pointing at the correct argument.
  for(reserved in admixtools:::.QPADM_MULTI_RESERVED) {
    expect_error(
      admixtools:::.qpadm_multi_dispatch_dots(
        setNames(list(NULL), reserved),
        fit_fun = admixtools::qpadm, on_f2_path = FALSE),
      "Reserved argument",
      info = paste("reserved kwarg:", reserved))
  }
  # Multiple reserved names: error message lists all of them.
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(target = "X", verbose = TRUE),
      fit_fun = admixtools::qpadm, on_f2_path = FALSE),
    "Reserved argument\\(s\\).*(target.*verbose|verbose.*target)")
})


test_that(".qpadm_multi_dispatch_dots rejects geno-only kwargs on the f2 path", {
  # On the f2-data path, qpadm_multi never calls f4blockdat_from_geno.
  # Round-4 dispatch routed geno-only kwargs into $geno regardless — that
  # bucket then went unused, silently dropping the user's kwarg. Round-5
  # rejects with a clear error pointing at the path mismatch.
  expect_error(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(cm_file = "x.cm"), fit_fun = admixtools::qpadm, on_f2_path = TRUE),
    "cm_file.*apply only when `data` is a genotype file prefix")
  # On the geno path, the same kwarg is accepted (sanity check).
  expect_silent(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(cm_file = "x.cm"), fit_fun = admixtools::qpadm, on_f2_path = FALSE))
  # Shared kwargs (auto_only, blgsize, poly_only) are valid on the f2 path
  # because they're also fit-fun formals; rejection should only fire for
  # geno-only names.
  expect_silent(
    admixtools:::.qpadm_multi_dispatch_dots(
      list(auto_only = FALSE, blgsize = 0.1),
      fit_fun = admixtools::qpadm, on_f2_path = TRUE))
})


# ── integration tests: end-to-end through qpadm_multi / qpadm_sweep ────────

test_that("qpadm_multi forwards qpadm-only kwargs to qpadm without crash", {
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  models = tibble::tibble(
    left   = list(c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
    right  = list(c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
    target = "Switzerland_Bichon.SG")
  # Pre-round-4, this call would crash with "The following arguments are
  # not used: singular_threshold" because the wrapper chain routed `...`
  # through qpadm()'s `...` with f4blocks set, tripping qpadm's own
  # validation guard. Pre-round-5, attempts to fix it left the geno-path
  # the same. Round-5 fixes the f2-path too.
  res = suppressMessages(suppressWarnings(
    qpadm_multi(ef, models, verbose = FALSE,
                singular_threshold = 1e-12)))   # 1e-12 << rcond → no trip
  expect_length(res, 1)
  expect_true("f4_var_rcond" %in% names(res[[1]]))
  expect_true(is.finite(res[[1]]$f4_var_rcond))
})


test_that("qpadm_multi rejects reserved kwarg with a clear error at entry", {
  # User error: passing `target` via `...` (would have collided with the
  # base_args `target = TRUE` at do.call time and raised a confusing
  # "matched by multiple actual arguments"). Now caught at entry.
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  models = tibble::tibble(
    left   = list(c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
    right  = list(c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
    target = "Switzerland_Bichon.SG")
  expect_error(
    qpadm_multi(ef, models, verbose = FALSE, target = "other_pop"),
    "Reserved argument")
})


test_that("qpadm_multi rejects geno-only kwarg on f2 path", {
  # Pre-round-5, cm_file (an f4blockdat_from_geno formal) was silently
  # dropped when passed alongside f2 data. Now rejected at entry.
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  models = tibble::tibble(
    left   = list(c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
    right  = list(c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
    target = "Switzerland_Bichon.SG")
  expect_error(
    qpadm_multi(ef, models, verbose = FALSE, cm_file = "x.cm"),
    "apply only when `data` is a genotype file prefix")
})


test_that("qpadm_multi raises informative per-fit error when singular_threshold trips", {
  # When singular_threshold trips inside a per-fit qpadm() call, the
  # tryCatch wrapper re-wraps the error via rlang::abort(parent = e) so
  # the failing combo's model label is in the message AND the original
  # condition chain (class, call) is preserved for class-specific
  # downstream handlers.
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


test_that("qpadm_multi preserves classed errors via rlang parent chain", {
  # Round-4 used stop() which flattened any classed condition (rlang or
  # custom) to simpleError. Round-5 uses rlang::abort(parent = e) which
  # keeps the original condition reachable via $parent, so a caller's
  # class-specific tryCatch can still dispatch.
  data(example_f2_blocks, package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  models = tibble::tibble(
    left   = list(c("Mbuti.DG", "Russia_Ust_Ishim.DG")),
    right  = list(c("Chimp.REF", "Altai_Neanderthal.DG", "Vindija.DG")),
    target = "Switzerland_Bichon.SG")
  caught_parent = NULL
  tryCatch(
    suppressMessages(suppressWarnings(
      qpadm_multi(ef, models, verbose = FALSE, singular_threshold = 1.0))),
    error = function(e) {
      caught_parent <<- conditionCall(rlang::cnd_entrace(e)$parent) %||%
                        rlang::cnd_message(e$parent)
    })
  # We don't pin the exact class (qpadm's internal stop() raises a plain
  # simpleError at the time of writing), but we DO pin that `$parent`
  # exists on the rewrapped error — proof the rewrap preserves the chain.
  expect_false(is.null(caught_parent))
})


test_that("qpadm_sweep entry-level rejects unknown kwargs through dispatch chain", {
  # qpadm_sweep forwards `...` to qpadm_multi; an unknown kwarg surfaces
  # at qpadm_multi's entry with a clear single error, not a cryptic
  # downstream error from f4blockdat_from_geno or qpadm() validation.
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
