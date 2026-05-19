# Tests for the rcond + SVD-based collinearity diagnostic on qpadm()
# (closes pipeline issue #8). The diagnostic surfaces two fields on the
# qpadm return list:
#
#   f4_var_rcond            -- always present; reciprocal condition number
#                              of the (fudged) f4 variance matrix.
#   f4_var_singular_loadings -- populated only when f4_var_rcond < 1e-8;
#                              tibble with one row per right pop ranked by
#                              L2 loading on the smallest singular vector.
#
# Plus a new singular_threshold = NA_real_ argument that, when set,
# converts the silent pseudo-inverse fallback into a fail-loud error.
#
# Singular-case fixture: clone one of example_f2_blocks's populations
# into a new dimname with bit-identical f2 values. The cloned pair is
# perfectly collinear with respect to every other population, which
# makes the unfudged f4_var rank-deficient. qpadm's default `fudge =
# 1e-4` regularizes that deficiency away (rcond ~ 1e-4), so to exercise
# the auto-warn-at-1e-8 path we pass `fudge = 1e-12`; for the fail-loud
# path we keep default fudge and set `singular_threshold` above the
# fudge-induced rcond.

.clone_pop_in_f2_blocks = function(f2, src = "Mbuti.DG", dst = "Mbuti_clone") {
  nam = dimnames(f2)[[1]]
  if(!(src %in% nam)) stop("source pop not in f2_blocks: ", src)
  if(dst %in% nam)    stop("destination pop already exists: ", dst)
  n   = length(nam)
  nb  = dim(f2)[3]
  new_nam = c(nam, dst)
  out = array(NA_real_, dim = c(n + 1, n + 1, nb),
              dimnames = list(new_nam, new_nam, dimnames(f2)[[3]]))
  out[1:n, 1:n, ] = f2
  src_idx = match(src, nam)
  # Make column (n+1) and row (n+1) byte-identical to src's column/row.
  out[1:n,   n+1, ] = f2[1:n, src_idx, ]
  out[n+1,   1:n, ] = f2[src_idx, 1:n, ]
  out[n+1,   n+1, ] = 0
  out
}

# Standard 3-source / 4-right qpadm fit on example_f2_blocks.
.example_setup = function() {
  data("example_f2_blocks", package = "admixtools", envir = environment())
  list(
    f2     = get("example_f2_blocks", envir = environment()),
    target = "Denisova.DG",
    left   = c("Altai_Neanderthal.DG", "Vindija.DG"),
    right  = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG",
               "Switzerland_Bichon.SG")
  )
}

# ── clean case ──────────────────────────────────────────────────────────────

test_that("clean qpadm fit reports finite f4_var_rcond above the concern bar with no loadings table", {
  s = .example_setup()
  res = suppressMessages(suppressWarnings(
    qpadm(s$f2, left = s$left, right = s$right, target = s$target,
          verbose = FALSE)
  ))
  expect_true("f4_var_rcond" %in% names(res))
  expect_true(is.finite(res$f4_var_rcond))
  expect_gt(res$f4_var_rcond, 1e-8)
  expect_null(res$f4_var_singular_loadings)
})

# ── singular case, warn-only (auto-1e-8 bar) ───────────────────────────────

test_that("clone-pair with tiny fudge triggers the loadings table and a warning at the auto-1e-8 bar", {
  # Construct a clone of Mbuti.DG and put both copies in right[-1]. The two
  # generate identical f4 statistics for every (left, right[1], right[-1][i])
  # tuple, so f4_var has identical rows -> rank-deficient. With default
  # fudge = 1e-4 the regularization lifts rcond above the 1e-8 bar; with
  # fudge = 1e-12 the deficiency survives and the auto-bar trips. The
  # smallest singular vector is +/- 1/sqrt(2*nl) on the two clones, giving
  # row norms of sqrt(1/2) ~ 0.707.
  s = .example_setup()
  f2_clone = .clone_pop_in_f2_blocks(s$f2, src = "Mbuti.DG", dst = "Mbuti_clone")
  right = c("Chimp.REF", "Mbuti.DG", "Mbuti_clone", "Russia_Ust_Ishim.DG",
            "Switzerland_Bichon.SG")

  warn_msg = NULL
  res = withCallingHandlers(
    suppressMessages(qpadm(f2_clone, left = s$left, right = right,
                            target = s$target, fudge = 1e-12,
                            verbose = TRUE)),
    warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") })

  # rcond near machine epsilon, loadings populated, warning fired.
  expect_lt(res$f4_var_rcond, 1e-8)
  expect_s3_class(res$f4_var_singular_loadings, "tbl_df")
  expect_match(warn_msg, "near-singular", fixed = TRUE)

  # The loadings tibble has one row per right[-1] pop (length(right) - 1 = 4).
  # The two clones should both have loading ~ sqrt(1/2); the others ~ 0.
  expect_equal(nrow(res$f4_var_singular_loadings), length(right) - 1L)
  expect_true(all(c("Mbuti.DG", "Mbuti_clone") %in% res$f4_var_singular_loadings$right))
  clone_loadings = res$f4_var_singular_loadings$loading[
    res$f4_var_singular_loadings$right %in% c("Mbuti.DG", "Mbuti_clone")]
  expect_equal(clone_loadings, c(sqrt(0.5), sqrt(0.5)), tolerance = 1e-4)
})

# ── fail-loud case (singular_threshold set above the post-fudge rcond) ─────

test_that("singular_threshold converts the warning into a fail-loud error with the loadings table", {
  # Same clone-pair fixture, default fudge (rcond ~ 1e-4 after regularization),
  # but the caller sets singular_threshold = 1e-3, which is above the post-
  # fudge rcond. The diagnostic must fire at the threshold rather than at
  # the auto-1e-8 bar.
  s = .example_setup()
  f2_clone = .clone_pop_in_f2_blocks(s$f2, src = "Mbuti.DG", dst = "Mbuti_clone")
  right = c("Chimp.REF", "Mbuti.DG", "Mbuti_clone", "Russia_Ust_Ishim.DG",
            "Switzerland_Bichon.SG")

  expect_error(
    suppressMessages(suppressWarnings(
      qpadm(f2_clone, left = s$left, right = right, target = s$target,
            singular_threshold = 1e-3, verbose = FALSE)
    )),
    "near-singular"
  )

  # Error message includes the diagnostic loadings table inline.
  err = tryCatch(
    suppressMessages(suppressWarnings(
      qpadm(f2_clone, left = s$left, right = right, target = s$target,
            singular_threshold = 1e-3, verbose = FALSE)
    )),
    error = function(e) conditionMessage(e))
  expect_match(err, "Mbuti.DG",    fixed = TRUE)
  expect_match(err, "Mbuti_clone", fixed = TRUE)
  expect_match(err, "Drop one and refit", fixed = TRUE)
})

# ── qpwave plumbing ─────────────────────────────────────────────────────────

test_that("qpwave passes singular_threshold through to qpadm and errors the same way", {
  # qpwave is qpadm(target = NULL), so the singular_threshold gate must
  # plumb through identically.
  s = .example_setup()
  f2_clone = .clone_pop_in_f2_blocks(s$f2, src = "Mbuti.DG", dst = "Mbuti_clone")
  right = c("Chimp.REF", "Mbuti.DG", "Mbuti_clone", "Russia_Ust_Ishim.DG",
            "Switzerland_Bichon.SG")

  expect_error(
    suppressMessages(suppressWarnings(
      qpwave(f2_clone, left = s$left, right = right,
             singular_threshold = 1e-3, verbose = FALSE)
    )),
    "near-singular"
  )

  # And without the threshold (and with default fudge): rcond is regularized
  # above the auto-bar, so loadings is NULL even though the underlying f4_var
  # is rank-deficient. This is the documented behavior and matches the
  # design intent (fudge is supposed to swallow mild collinearity silently).
  res = suppressMessages(suppressWarnings(
    qpwave(f2_clone, left = s$left, right = right, verbose = FALSE)
  ))
  expect_true(is.finite(res$f4_var_rcond))
  expect_null(res$f4_var_singular_loadings)
})
