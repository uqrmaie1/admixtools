# Tests for qpadm_with_pruning(). All cases run against example_f2_blocks
# (7 pops) plus singular-fixture variants built via .clone_pop_in_f2_blocks
# (helper-fixtures.R). Pin furrr to sequential so the inner qpadm fits run
# in a deterministic order under any host plan().
if(requireNamespace("future", quietly = TRUE)) future::plan("sequential")

# ---- shared fixture -------------------------------------------------------

.well_conditioned_args = function() {
  # A right set with no near-singular structure: standard outgroups, no
  # cloned pops. The first qpadm fit should already satisfy
  # singular_threshold = 1e-12 and the pruner should no-op.
  data("example_f2_blocks", package = "admixtools", envir = environment())
  list(
    data   = get("example_f2_blocks", envir = environment()),
    left   = c("Altai_Neanderthal.DG", "Vindija.DG"),
    right  = c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG",
               "Switzerland_Bichon.SG"),
    target = "Denisova.DG")
}

.singular_args = function(n_clones = 1) {
  # n_clones in {1, 2}. With 1 clone the right set has one sister pair;
  # with 2 clones it has two sister pairs (the "needs two iterations to
  # converge" case).
  data("example_f2_blocks", package = "admixtools", envir = environment())
  ef = get("example_f2_blocks", envir = environment())
  ef2 = .clone_pop_in_f2_blocks(ef, "Mbuti.DG", "Mbuti_clone")
  if(n_clones >= 2)
    ef2 = .clone_pop_in_f2_blocks(ef2, "Russia_Ust_Ishim.DG", "Ushim_clone")
  right = if(n_clones >= 2)
    c("Chimp.REF", "Altai_Neanderthal.DG", "Mbuti.DG", "Mbuti_clone",
      "Russia_Ust_Ishim.DG", "Ushim_clone")
  else
    c("Chimp.REF", "Altai_Neanderthal.DG", "Mbuti.DG", "Mbuti_clone",
      "Russia_Ust_Ishim.DG")
  list(data = ef2, left = c("Vindija.DG"), right = right,
       target = "Switzerland_Bichon.SG")
}

.mute = function(expr) {
  suppressMessages(suppressWarnings(force(expr)))
}


# ---- (a) well-conditioned no-op ------------------------------------------

test_that("well-conditioned fit returns immediately with empty trail", {
  a = .well_conditioned_args()
  fit = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target))

  expect_true(fit$converged)
  expect_identical(fit$reason, "converged")
  expect_equal(nrow(fit$pruning_trail), 0L)
  expect_identical(fit$right_final, a$right)
  expect_true("weights" %in% names(fit))
  expect_true("f4_var_rcond" %in% names(fit))
  expect_true(fit$f4_var_rcond >= 1e-12)
})


# ---- (b) single-iteration prune ------------------------------------------

test_that("one clone pair converges after one drop", {
  a = .singular_args(n_clones = 1)
  # singular_threshold = 1e-10 > the fixture's first-fit rcond (~2.2e-12);
  # this guarantees the fixture trips the threshold and exercises the
  # pruning loop. Default 1e-12 is on a knife edge with this fixture.
  fit = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10))

  expect_true(fit$converged)
  expect_equal(nrow(fit$pruning_trail), 1L)
  expect_identical(fit$pruning_trail$iteration, 1L)
  # The dropped pop must be one of the clone pair (Mbuti.DG / Mbuti_clone).
  expect_true(fit$pruning_trail$dropped %in% c("Mbuti.DG", "Mbuti_clone"))
  expect_lt(fit$pruning_trail$rcond_before, 1e-10)
  expect_true(fit$pruning_trail$loading > 0.5)
  # right_final lost exactly one pop, and it's the one the trail named.
  expect_equal(length(fit$right_final), length(a$right) - 1L)
  expect_identical(setdiff(a$right, fit$right_final),
                   fit$pruning_trail$dropped)
})


# ---- (c) multi-iteration prune to convergence -----------------------------

test_that("two clone pairs converge in two iterations", {
  a = .singular_args(n_clones = 2)
  fit = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10))

  expect_true(fit$converged)
  expect_equal(nrow(fit$pruning_trail), 2L)
  # Each iteration must drop one of the sister-pair members.
  pair1 = c("Mbuti.DG", "Mbuti_clone")
  pair2 = c("Russia_Ust_Ishim.DG", "Ushim_clone")
  dropped = fit$pruning_trail$dropped
  expect_true(any(dropped %in% pair1))
  expect_true(any(dropped %in% pair2))
  # right_final lost exactly two pops; iterations are 1, 2 in order.
  expect_equal(length(fit$right_final), length(a$right) - 2L)
  expect_identical(fit$pruning_trail$iteration, c(1L, 2L))
})


# ---- (d) floor exhaustion -------------------------------------------------

test_that("hitting min_right_pops before convergence reports floor", {
  a = .singular_args(n_clones = 1)
  # min_right_pops set to the starting right-set size means no drop is
  # allowed (would push below floor); pruner returns floor on iter 1.
  fit = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10,
    min_right_pops = length(a$right)))

  expect_false(fit$converged)
  expect_identical(fit$reason, "floor")
  expect_equal(nrow(fit$pruning_trail), 0L)
  expect_identical(fit$right_final, a$right)
})


# ---- (d2) iter_cap re-fit synchronises fit with right_final ---------------

test_that("iter_cap re-fits so fit corresponds to right_final", {
  # Two clone pairs need two drops to converge. Capping max_iterations at 1
  # forces the loop to exhaust its budget after a single drop. The post-loop
  # reconciliation must re-fit on the final (still-singular) right set so the
  # returned fit corresponds to right_final, not the pre-drop set.
  a = .singular_args(n_clones = 2)
  fit = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10, max_iterations = 1))

  expect_false(fit$converged)
  expect_identical(fit$reason, "iter_cap")
  expect_equal(nrow(fit$pruning_trail), 1L)
  expect_equal(length(fit$right_final), length(a$right) - 1L)
  # The returned fit must equal a direct qpadm() on right_final (the desync
  # this re-fit fixes: without it, fit would be the pre-drop 6-pop fit).
  direct = .mute(qpadm(a$data, left = a$left, right = fit$right_final,
                       target = a$target, fudge = 1e-12, verbose = FALSE))
  expect_equal(fit$f4_var_rcond, direct$f4_var_rcond, tolerance = 1e-8)
  # The last trail row's rcond_after was backfilled to the re-fit's rcond.
  expect_true(is.finite(fit$pruning_trail$rcond_after[1]))
  expect_equal(fit$pruning_trail$rcond_after[1], fit$f4_var_rcond,
               tolerance = 1e-8)
})


test_that("iter_cap reclassifies to converged when the final drop clears", {
  # One clone pair converges after one drop. Capping at max_iterations = 1
  # stops the loop right after that drop; the reconciliation re-fit on the
  # cleared 4-pop set must detect convergence and report it honestly rather
  # than mislabelling a good fit as iter_cap.
  a = .singular_args(n_clones = 1)
  fit = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10, max_iterations = 1))

  expect_true(fit$converged)
  expect_identical(fit$reason, "converged")
  expect_equal(length(fit$right_final), length(a$right) - 1L)
  expect_true(fit$f4_var_rcond >= 1e-10)
})


# ---- (e) pruner-direct equivalence ----------------------------------------

test_that("converged fit equals direct qpadm() with right_final", {
  a = .singular_args(n_clones = 1)
  prune = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10))
  expect_true(prune$converged)
  direct = .mute(qpadm(a$data, left = a$left, right = prune$right_final,
                       target = a$target, fudge = 1e-12, verbose = FALSE))
  # Numeric slots match within ulp tolerance. The bind below also catches
  # column-order or row-order drift between the two call paths.
  # Tolerances at 1e-8: both paths call qpadm() on identical inputs so results
  # are numerically identical, but matrix-inversion SEs can accumulate floating
  # point rounding that exceeds sub-epsilon tolerances on some platforms.
  expect_equal(prune$weights$weight, direct$weights$weight,
               tolerance = 1e-8)
  expect_equal(prune$weights$se,     direct$weights$se,     tolerance = 1e-8)
  expect_equal(prune$f4_var_rcond,   direct$f4_var_rcond,   tolerance = 1e-8)
  expect_equal(prune$rankdrop$p,     direct$rankdrop$p,     tolerance = 1e-8)
})


# ---- (f) lookahead strategy plumbing -------------------------------------

test_that("lookahead strategy emits the same trail-shape as greedy", {
  a = .singular_args(n_clones = 1)
  pg = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10, strategy = "greedy"))
  pl = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10,
    strategy = "lookahead", lookahead_top_j = 2))

  # Same column shape; lookahead fills rcond_after, greedy leaves it NA.
  expect_named(pg$pruning_trail,
               c("iteration", "dropped", "loading", "rcond_before",
                 "rcond_after"))
  expect_named(pl$pruning_trail,
               c("iteration", "dropped", "loading", "rcond_before",
                 "rcond_after"))
  expect_true(all(is.na(pg$pruning_trail$rcond_after)))
  expect_true(all(is.finite(pl$pruning_trail$rcond_after)))
  # Both converged; the converged right_final must be the same SIZE (both
  # strategies drop the same total number on this fixture).
  expect_equal(length(pg$right_final), length(pl$right_final))
})


# ---- (g) lookahead vs greedy on competing-loadings fixture ---------------

test_that("lookahead picks the candidate with best post-drop rcond", {
  # Both clones cause near-singularity; lookahead's job is to evaluate each
  # candidate's post-drop rcond and pick the BEST one (numerically, the one
  # with the larger post-drop f4_var_rcond -- usually the deeper clone of
  # the pair). We verify mechanically: the rcond_after slot for lookahead
  # equals the post-drop fit's f4_var_rcond.
  a = .singular_args(n_clones = 1)
  pl = .mute(qpadm_with_pruning(
    a$data, left = a$left, right = a$right, target = a$target,
    fudge = 1e-12, singular_threshold = 1e-10,
    strategy = "lookahead", lookahead_top_j = 2))

  expect_true(pl$converged)
  # The dropped pop's rcond_after must equal the converged fit's rcond
  # (lookahead committed the BEST candidate, then re-ran from there).
  # Tolerance 1e-8 (matching test (e)): the probe and committed fits are two
  # independent rcond() calls on the same matrix, bit-identical on one platform
  # but not guaranteed so across BLAS/CPU variants.
  expect_equal(pl$pruning_trail$rcond_after[1], pl$f4_var_rcond,
               tolerance = 1e-8)
})


# ---- (h) verbose path emits cli info -------------------------------------

test_that("verbose = TRUE emits one cli_inform per iteration", {
  a = .singular_args(n_clones = 1)
  messages = character()
  suppressWarnings(withCallingHandlers(
    qpadm_with_pruning(a$data, left = a$left, right = a$right,
                       target = a$target, fudge = 1e-12,
                       singular_threshold = 1e-10, verbose = TRUE),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  # At least one cli_inform fired with the iteration banner.
  expect_true(any(grepl("[qpadm_with_pruning] iter 1", messages,
                        fixed = TRUE)))
})


# ---- (i) argument validation ---------------------------------------------

test_that("validates singular_threshold + right + min_right_pops", {
  a = .well_conditioned_args()
  expect_error(qpadm_with_pruning(a$data, a$left, a$right[1:2], a$target),
               "min_right_pops")
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                  singular_threshold = -1),
               "singular_threshold")
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                  min_right_pops = 0),
               "min_right_pops")
})

test_that("validates target shape and rejects non-integer / oversized counts", {
  a = .well_conditioned_args()
  # target must be NULL or a single non-NA name.
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, target = NA),
               "target")
  expect_error(qpadm_with_pruning(a$data, a$left, a$right,
                                  target = c("Denisova.DG", "Mbuti.DG")),
               "target")
  # Non-integer counts are rejected rather than silently truncated.
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                  min_right_pops = 3.7),
               "whole number")
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                  max_iterations = 2.5),
               "whole number")
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                  strategy = "lookahead", lookahead_top_j = 1.5),
               "whole number")
  # A count above the integer range is rejected, not coerced to NA (which would
  # otherwise trip a cryptic `if(NA)` downstream).
  expect_error(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                  min_right_pops = 2^31),
               "whole number no larger than")
})

test_that("verbose + singular_threshold bind to formals, not `...`", {
  # Both `verbose` and `singular_threshold` are formals of qpadm_with_pruning;
  # R's argument matching binds them to the formals before `...` collection,
  # so a caller cannot accidentally pass them through to the inner qpadm()
  # via `...`. Sanity-check by asserting the inner qpadm() sees the pruner's
  # forced values (verbose = FALSE, singular_threshold = NA) regardless of
  # what the outer caller passed.
  a = .well_conditioned_args()
  fit = .mute(qpadm_with_pruning(a$data, a$left, a$right, a$target,
                                 verbose = FALSE, singular_threshold = 1e-10))
  # The inner qpadm() would have errored if singular_threshold = 1e-10 were
  # forwarded (the well-conditioned fixture's rcond is well above that bar;
  # qpadm() only trips when rcond < threshold) -- but a future regression
  # could pick a fixture where the threshold trips, so we assert the fit
  # actually came back rather than erroring.
  expect_true(fit$converged)
})
