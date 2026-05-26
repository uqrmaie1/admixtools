# Tests for qpadm_sweep() -- all cases can run against example_f2_blocks
# (7 populations: Altai_Neanderthal.DG, Chimp.REF, Denisova.DG, Mbuti.DG,
#  Russia_Ust_Ishim.DG, Switzerland_Bichon.SG, Vindija.DG).
#
# Pop-set design:  each (target, source_set, right_set) combo must have no
# duplicate population names.  The fixtures below are pre-verified to work
# without triggering qpadm's "Duplicated pops" guard.
#
#   FEASIBLE model:   target=Switzerland_Bichon.SG,
#                     left  =c("Russia_Ust_Ishim.DG","Mbuti.DG"),
#                     right =c("Altai_Neanderthal.DG","Chimp.REF","Denisova.DG")
#                     → weights ≈ 0.90, 0.10 (all in [0,1])
#
#   INFEASIBLE model: target=Vindija.DG,
#                     left  =c("Chimp.REF","Mbuti.DG"),
#                     right =c("Altai_Neanderthal.DG","Russia_Ust_Ishim.DG","Switzerland_Bichon.SG")
#                     → weights ≈ 2.35, -1.35 (out of [0,1])

# ── fixtures ────────────────────────────────────────────────────────────────

.f2 = function() {
  data("example_f2_blocks", package = "admixtools", envir = environment())
  get("example_f2_blocks", envir = environment())
}

# Feasible 1×1×1 setup
.feas_args = function() list(
  data        = .f2(),
  targets     = "Switzerland_Bichon.SG",
  source_sets = list(modern = c("Russia_Ust_Ishim.DG", "Mbuti.DG")),
  right_sets  = list(refs   = c("Altai_Neanderthal.DG", "Chimp.REF", "Denisova.DG")),
  verbose     = FALSE
)

# Infeasible 1×1×1 setup
.infeas_args = function() list(
  data        = .f2(),
  targets     = "Vindija.DG",
  source_sets = list(archaic = c("Chimp.REF", "Mbuti.DG")),
  right_sets  = list(refs    = c("Altai_Neanderthal.DG", "Russia_Ust_Ishim.DG",
                                 "Switzerland_Bichon.SG")),
  verbose     = FALSE
)

# ── test 1: basic call returns the right dimensions and columns ──────────────

test_that("basic call returns (n_t * n_s * n_r)-row tibble with documented columns", {
  res = do.call(qpadm_sweep, .feas_args())

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1L)  # 1 target × 1 source_set × 1 right_set

  flat_cols = c("target", "source_set", "right_set",
                "left", "right", "f4rank", "p", "chisq", "dof", "feasible")
  for(col in flat_cols) expect_true(col %in% names(res), info = paste("missing column:", col))

  # full_results = TRUE (default) adds list-columns
  expect_true("weights"  %in% names(res))
  expect_true("rankdrop" %in% names(res))
})

test_that("row count scales as n_t * n_s * n_r", {
  f2 = .f2()
  res = qpadm_sweep(
    f2,
    targets     = "Switzerland_Bichon.SG",
    source_sets = list(modern = c("Russia_Ust_Ishim.DG", "Mbuti.DG")),
    right_sets  = list(
      r1 = c("Altai_Neanderthal.DG", "Chimp.REF", "Denisova.DG"),
      r2 = c("Altai_Neanderthal.DG", "Chimp.REF", "Vindija.DG")
    ),
    verbose = FALSE
  )
  expect_equal(nrow(res), 1L * 1L * 2L)
})

# ── test 2: auto-naming ──────────────────────────────────────────────────────

test_that("unnamed source_sets and right_sets get S1/S2... and R1/R2... labels", {
  f2 = .f2()
  res = qpadm_sweep(
    f2,
    targets     = "Switzerland_Bichon.SG",
    source_sets = list(c("Russia_Ust_Ishim.DG", "Mbuti.DG")),  # unnamed
    right_sets  = list(c("Altai_Neanderthal.DG", "Chimp.REF", "Denisova.DG")),  # unnamed
    verbose     = FALSE
  )
  expect_equal(res$source_set, "S1")
  expect_equal(res$right_set,  "R1")
})

# ── test 3: empty-name error ─────────────────────────────────────────────────

test_that("empty source_sets name raises a clear error", {
  f2 = .f2()
  bad_src = list(c("Russia_Ust_Ishim.DG", "Mbuti.DG"))
  names(bad_src) = ""
  expect_error(
    qpadm_sweep(f2, "Switzerland_Bichon.SG", bad_src,
                list(refs = c("Altai_Neanderthal.DG", "Chimp.REF")), verbose = FALSE),
    "named"
  )
})

test_that("empty right_sets name raises a clear error", {
  f2 = .f2()
  bad_rgt = list(c("Altai_Neanderthal.DG", "Chimp.REF"))
  names(bad_rgt) = ""
  expect_error(
    qpadm_sweep(f2, "Switzerland_Bichon.SG",
                list(modern = c("Russia_Ust_Ishim.DG", "Mbuti.DG")),
                bad_rgt, verbose = FALSE),
    "named"
  )
})

# ── test 4: duplicate-name error ─────────────────────────────────────────────

test_that("duplicate source_sets names raise an error", {
  f2 = .f2()
  dup = list(
    canonical = c("Russia_Ust_Ishim.DG"),
    canonical = c("Mbuti.DG")
  )
  expect_error(
    qpadm_sweep(f2, "Switzerland_Bichon.SG", dup,
                list(refs = c("Altai_Neanderthal.DG", "Chimp.REF")), verbose = FALSE),
    "unique"
  )
})

test_that("duplicate right_sets names raise an error", {
  f2 = .f2()
  dup = list(refs = c("Altai_Neanderthal.DG"), refs = c("Chimp.REF"))
  expect_error(
    qpadm_sweep(f2, "Switzerland_Bichon.SG",
                list(modern = c("Russia_Ust_Ishim.DG", "Mbuti.DG")),
                dup, verbose = FALSE),
    "unique"
  )
})

# ── test 5: non-character-vector error ───────────────────────────────────────

test_that("integer source_sets entry raises a clear error", {
  f2 = .f2()
  expect_error(
    qpadm_sweep(f2, "Switzerland_Bichon.SG",
                list(bad = 1:3),
                list(refs = c("Altai_Neanderthal.DG", "Chimp.REF")), verbose = FALSE),
    "character"
  )
})

test_that("integer right_sets entry raises a clear error", {
  f2 = .f2()
  expect_error(
    qpadm_sweep(f2, "Switzerland_Bichon.SG",
                list(modern = c("Russia_Ust_Ishim.DG", "Mbuti.DG")),
                list(bad = 1:3), verbose = FALSE),
    "character"
  )
})

# ── test 6: full_results = FALSE drops list-columns ──────────────────────────

test_that("full_results = FALSE drops weights and rankdrop list-columns", {
  res = do.call(qpadm_sweep, c(.feas_args(), list(full_results = FALSE)))

  expect_false("weights"  %in% names(res))
  expect_false("rankdrop" %in% names(res))

  flat_cols = c("target", "source_set", "right_set",
                "left", "right", "f4rank", "p", "chisq", "dof", "feasible")
  for(col in flat_cols) expect_true(col %in% names(res), info = paste("missing column:", col))
})

# ── test 7: feasible column ──────────────────────────────────────────────────

test_that("feasible is TRUE when all weights are in [0, 1]", {
  res = do.call(qpadm_sweep, .feas_args())
  expect_true(isTRUE(res$feasible))
})

test_that("feasible is FALSE when any weight is outside [0, 1]", {
  res = do.call(qpadm_sweep, .infeas_args())
  expect_true(isFALSE(res$feasible))
})

# ── test 8: f4_var_rcond + f4_var_singular_loadings (issue #16) ─────────────
#
# Acceptance criteria from fork issue #16: qpadm_sweep(..., full_results = TRUE)
# must surface f4_var_rcond as a numeric column and f4_var_singular_loadings as
# a list-column. Both values must match what a direct qpadm() call with the
# same inputs would return — proving the flatten preserves rather than
# recomputes the underlying qpadm() values.

test_that("qpadm_sweep surfaces f4_var_rcond + f4_var_singular_loadings (issue #16)", {
  res = suppressWarnings(do.call(qpadm_sweep, .feas_args()))

  # New columns are present and well-typed.
  expect_true("f4_var_rcond" %in% names(res))
  expect_true("f4_var_singular_loadings" %in% names(res))
  expect_type(res$f4_var_rcond, "double")
  expect_type(res$f4_var_singular_loadings, "list")
  expect_length(res$f4_var_rcond, nrow(res))
  expect_length(res$f4_var_singular_loadings, nrow(res))

  # Values match a direct qpadm() call with the swept row's inputs — proves
  # the flatten preserves rather than recomputes (a regression that recomputed
  # would silently de-couple the diagnostic from the actual fit).
  args = .feas_args()
  direct = suppressWarnings(qpadm(
    data   = args$data,
    target = args$targets[1],
    left   = args$source_sets[[1]],
    right  = args$right_sets[[1]],
    verbose = FALSE))

  expect_equal(res$f4_var_rcond[1], direct$f4_var_rcond, tolerance = 1e-12)
  expect_identical(res$f4_var_singular_loadings[[1]],
                   direct$f4_var_singular_loadings)
})

test_that("qpadm_sweep returns f4_var_rcond even when full_results = FALSE", {
  # Design choice: f4_var_rcond is a scalar diagnostic comparable to p / chisq —
  # surfaced unconditionally so a pruner using full_results = FALSE for memory
  # efficiency can still gate on rank deficiency. The list-column
  # f4_var_singular_loadings is gated.
  res = suppressWarnings(do.call(qpadm_sweep,
                                 c(.feas_args(), list(full_results = FALSE))))
  expect_true("f4_var_rcond" %in% names(res))
  expect_type(res$f4_var_rcond, "double")
  expect_false("f4_var_singular_loadings" %in% names(res))
})

# ── test 8: n=1 edge case ────────────────────────────────────────────────────

test_that("single target x single source-set x single right-set works without error", {
  res = do.call(qpadm_sweep, .feas_args())

  expect_equal(nrow(res), 1L)
  expect_true(is.numeric(res$p))
  expect_true(!is.na(res$p))
})
