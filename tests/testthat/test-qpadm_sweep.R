# Tests for qpadm_sweep() -- all cases can run against example_f2_blocks
# (7 populations: Altai_Neanderthal.DG, Chimp.REF, Denisova.DG, Mbuti.DG,
#  Russia_Ust_Ishim.DG, Switzerland_Bichon.SG, Vindija.DG).
#
# Future plan pin (file-scoped): qpadm_multi dispatches qpadm() to
# furrr::future_map. If a developer's .Rprofile or the CI runner has set
# future::plan('multisession') or 'multicore', the worker process's LAPACK
# rcond computation can return bit-different doubles from the main process
# (BLAS thread count, FMA, OMP variations); the expect_identical assertions
# on f4_var_rcond below would then flake.
#
# Setting the plan directly at file top is the simplest mechanism that
# actually works. (Round-3 used withr::local_options with .local_envir =
# parent.frame() inside a local() block — but parent.frame() inside local()
# resolves to the local()'s own frame, which exits immediately, restoring
# the prior plan before any test runs. Empirically: "After local block:
# plan = multisession" if multisession was active before.) The plan
# survives between this file's tests and subsequent test files in the same
# test_dir() invocation, which is the testthat-recommended scope.
if(requireNamespace("future", quietly = TRUE)) future::plan("sequential")

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
                "left", "right", "f4rank", "p", "chisq", "dof", "feasible",
                "f4_var_rcond")
  for(col in flat_cols) expect_true(col %in% names(res), info = paste("missing column:", col))

  # full_results = TRUE (default) adds list-columns
  expect_true("weights"  %in% names(res))
  expect_true("rankdrop" %in% names(res))
  expect_true("f4_var_singular_loadings" %in% names(res))
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

test_that("full_results = FALSE drops weights / rankdrop / loadings; keeps rcond", {
  res = do.call(qpadm_sweep, c(.feas_args(), list(full_results = FALSE)))

  # Gated list-columns are dropped.
  expect_false("weights"  %in% names(res))
  expect_false("rankdrop" %in% names(res))
  expect_false("f4_var_singular_loadings" %in% names(res))

  # Always-on flat surface (including the scalar rcond diagnostic) survives.
  flat_cols = c("target", "source_set", "right_set",
                "left", "right", "f4rank", "p", "chisq", "dof", "feasible",
                "f4_var_rcond")
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

# ── test 8: f4_var_rcond + f4_var_singular_loadings (fork issue #16) ─────────
#
# Coverage: shape (8a), multi-row alignment (8b), singular-positive round-trip
# (8c), defensive not-all-NA (8d), full_results = FALSE surface (8e).
#
# Implementation contract anchored by these tests:
#   - f4_var_rcond is the byte-for-byte value qpadm() returns (no recompute,
#     no fudge-of-fudge). expect_identical (no tolerance) catches drift from
#     a single-source-of-truth refactor.
#   - f4_var_singular_loadings is the actual tibble qpadm() emits, intact
#     through the lapply (no row reorder, no column rename, no recompute).
#   - Ordering: row i of qpadm_sweep's output corresponds to combo i of the
#     Cartesian product. The 2x2x2 sweep in 8b probes this directly.
#
# Narrow warning muffler: silence only the two predictable-noise patterns
# (the near-singular warning from qpadm() itself on the 8c fixture and
# furrr's UNRELIABLE VALUE warning). Anything else surfaces — keeping the
# regression-detection signal alive.

.mute_qpadm_sweep_noise = function(w) {
  msg = conditionMessage(w)
  # `fixed = TRUE` so noise patterns are interpreted as literal substrings.
  # A future maintainer adding a pattern containing regex metacharacters
  # (e.g. "[clone]", "(rcond)") gets literal match rather than silent
  # group / character-class semantics that would either misfire or error.
  noise = c("UNRELIABLE VALUE",   # furrr seed-handling on parallel plans
            "near-singular")      # auto-bar warning (may surface if qpadm_multi's
                                  # hard-coded verbose=FALSE is ever relaxed)
  if(any(vapply(noise, grepl, logical(1), msg, fixed = TRUE)))
    invokeRestart("muffleWarning")
}

test_that("8a: clean-case round-trip is byte-identical to direct qpadm()", {
  args = .feas_args()
  res  = withCallingHandlers(do.call(qpadm_sweep, args),
                             warning = .mute_qpadm_sweep_noise)

  # Shape + types.
  expect_true("f4_var_rcond" %in% names(res))
  expect_true("f4_var_singular_loadings" %in% names(res))
  expect_type(res$f4_var_rcond, "double")
  expect_type(res$f4_var_singular_loadings, "list")

  # Direct qpadm() with the swept row's exact inputs.
  direct = withCallingHandlers(
    qpadm(data = args$data, target = args$targets[1],
          left = args$source_sets[[1]], right = args$right_sets[[1]],
          verbose = FALSE),
    warning = .mute_qpadm_sweep_noise)

  # Byte-identical scalar — no tolerance, anchoring the no-recompute contract.
  expect_identical(res$f4_var_rcond[1], direct$f4_var_rcond)
  # Loadings cell matches the direct route. Either both NULL (the typical
  # clean-data case for this fixture) or both the same tibble.
  expect_identical(res$f4_var_singular_loadings[[1]],
                   direct$f4_var_singular_loadings)
})

test_that("8b: multi-row sweep preserves per-row alignment with direct qpadm()", {
  # Pop overlap rules out a true 2x2x2 in example_f2_blocks (7 pops, 2
  # disjoint targets eat 2, sources need ≥2 each, rights need ≥2 each,
  # and pairwise disjointness across all (target, source_set, right_set)
  # combos quickly exhausts the pool). 2 targets x 1 source_set x 2
  # right_sets = 4 rows is the largest valid grid here, and still exercises
  # the alignment invariant across two grid dimensions (target and right
  # vary; ordering bugs in either dimension surface).
  #
  # full_results = FALSE here (not the default TRUE) so this test exercises
  # the gated-column drop AND the always-on rcond surface simultaneously.
  # Combined with the all-finite assertion below it doubles as the defensive
  # guard against the "outer full_results=FALSE silently routes through
  # qpadm_p (no f4_var_rcond) -> all-NA column" regression scenario.
  f2 = .f2()
  targets = c("Switzerland_Bichon.SG", "Vindija.DG")
  sources = list(s1 = c("Russia_Ust_Ishim.DG", "Mbuti.DG"))
  rights  = list(r1 = c("Chimp.REF", "Altai_Neanderthal.DG"),
                 r2 = c("Chimp.REF", "Altai_Neanderthal.DG", "Denisova.DG"))

  res = withCallingHandlers(
    qpadm_sweep(f2, targets = targets, source_sets = sources,
                right_sets = rights, full_results = FALSE, verbose = FALSE),
    warning = .mute_qpadm_sweep_noise)

  expect_equal(nrow(res), 4L)
  expect_length(res$f4_var_rcond, 4L)

  # Cartesian-product completeness: every (target, source_set, right_set)
  # triple appears exactly once. Tested as a SET to be agnostic to the
  # underlying iteration order (base expand.grid vs tidyr::expand_grid
  # vary, and either is a valid implementation choice). A row-shuffle bug
  # that consistently shuffles inputs AND outputs together would otherwise
  # pass the per-row qpadm() round-trip below (because direct qpadm uses
  # the shuffled inputs); this set-equality check would still pass on a
  # consistent shuffle but fail on missing combos / duplicates / wrong
  # combos. Combined with the per-row qpadm() round-trip below, the two
  # together pin both completeness and per-row correctness.
  expected_triples = expand.grid(target = targets,
                                 source_set = names(sources),
                                 right_set = names(rights),
                                 stringsAsFactors = FALSE,
                                 KEEP.OUT.ATTRS = FALSE)
  actual_triples = data.frame(target = res$target,
                              source_set = res$source_set,
                              right_set = res$right_set,
                              stringsAsFactors = FALSE)
  # Set-equality via sorted-paste fingerprints (avoids dplyr).
  fp = function(d) sort(paste(d$target, d$source_set, d$right_set, sep = "|"))
  expect_identical(fp(actual_triples), fp(expected_triples))

  # Cartesian-product ORDERING (round-5 hardening). The round-trip loop below
  # checks "row i's inputs and row i's stored rcond are self-consistent",
  # which a row-shuffle that consistently shuffles inputs AND outputs
  # together would silently pass. This block pins the actual row order
  # against the documented expand.grid order — base expand.grid cycles the
  # first column fastest, so target varies inside source_set inside right_set.
  expect_identical(res$target,     expected_triples$target)
  expect_identical(res$source_set, expected_triples$source_set)
  expect_identical(res$right_set,  expected_triples$right_set)

  for(i in seq_len(nrow(res))) {
    # res$left and res$right are list-columns that round-trip the sweep's
    # row-to-combo mapping. Driving direct qpadm() from those columns
    # (rather than from the original Cartesian inputs) makes any
    # row-shuffling bug visible: the diagnostic on row i must match the
    # qpadm() call with row i's left and right.
    label = sprintf("row %d: target=%s source=%s right=%s",
                    i, res$target[i], res$source_set[i], res$right_set[i])
    direct = withCallingHandlers(
      qpadm(data = f2, target = res$target[i],
            left = res$left[[i]], right = res$right[[i]],
            verbose = FALSE),
      warning = .mute_qpadm_sweep_noise)
    expect_identical(res$f4_var_rcond[i], direct$f4_var_rcond,
                     label = label)
  }

  # Defensive: on the clean 4-row sweep every row should produce a finite
  # rcond. A future refactor that silently routes through qpadm_p (which
  # doesn't return f4_var_rcond) would produce all-NA here.
  expect_true(all(is.finite(res$f4_var_rcond)))
})

test_that("8c: singular fixture preserves the loadings tibble through the sweep", {
  # Closes the "NULL == NULL is vacuous" gap from earlier drafts of this
  # test. Drive rcond below the 1e-8 auto-bar by including a cloned pop
  # in the right set (perfectly collinear with its source) and using
  # fudge = 1e-12 so the unfudged-matrix singularity isn't regularized
  # away. Then verify the non-NULL loadings tibble survives the lapply
  # byte-for-byte.
  f2_clone = .clone_pop_in_f2_blocks(.f2(), src = "Mbuti.DG", dst = "Mbuti_clone")

  args = list(
    data        = f2_clone,
    targets     = "Switzerland_Bichon.SG",
    source_sets = list(s1 = c("Russia_Ust_Ishim.DG", "Vindija.DG")),
    right_sets  = list(r1 = c("Chimp.REF", "Altai_Neanderthal.DG",
                              "Mbuti.DG", "Mbuti_clone")),
    fudge       = 1e-12,
    verbose     = FALSE)

  res = withCallingHandlers(do.call(qpadm_sweep, args),
                            warning = .mute_qpadm_sweep_noise)
  direct = withCallingHandlers(
    qpadm(data = f2_clone, target = args$targets[1],
          left = args$source_sets[[1]], right = args$right_sets[[1]],
          fudge = 1e-12, verbose = FALSE),
    warning = .mute_qpadm_sweep_noise)

  # The fixture actually pushed rcond below the auto-bar. Source the
  # threshold from the package-private constant so a future re-tuning of
  # the auto-bar updates this test automatically (rather than silently
  # decoupling the test's literal from the implementation default).
  expect_lt(res$f4_var_rcond[1], admixtools:::.rcond_concern)

  # The loadings tibble is populated (not NULL) on both routes.
  expect_false(is.null(res$f4_var_singular_loadings[[1]]))
  expect_s3_class(res$f4_var_singular_loadings[[1]], "tbl_df")

  # Schema pin: column names + types. Catches a future qpadm refactor that
  # renames 'loading' to 'norm' (or adds columns) — the byte-identical
  # round-trip would still pass because both routes go through the same
  # code path, but downstream consumer code (vignette example arranges by
  # `loading`) would silently break.
  expect_named(res$f4_var_singular_loadings[[1]], c("right", "loading"))
  expect_type(res$f4_var_singular_loadings[[1]]$right,   "character")
  expect_type(res$f4_var_singular_loadings[[1]]$loading, "double")

  # Byte-for-byte preservation: tibble columns, row order, and values
  # are identical across the two routes. A recompute or reshape regression
  # surfaces here.
  expect_identical(res$f4_var_singular_loadings[[1]],
                   direct$f4_var_singular_loadings)
  expect_identical(res$f4_var_rcond[1], direct$f4_var_rcond)
})

test_that("8d: full_results = FALSE keeps f4_var_rcond (finite), drops loadings", {
  # Design choice pinned (and defensive all-NA guard): f4_var_rcond is a
  # scalar diagnostic comparable to p / chisq / dof — surfaced
  # unconditionally so a pruner using full_results = FALSE for smaller
  # outputs can still gate on rank deficiency. The is.finite() assertion
  # catches the regression scenario where a future qpadm_sweep refactor
  # propagates outer full_results=FALSE into qpadm_multi (routing through
  # qpadm_p, which doesn't return f4_var_rcond): the resulting all-NA
  # column would still pass an "exists + double type" check but fail this.
  # The list-column f4_var_singular_loadings is gated on full_results
  # because of variable per-row payload.
  res = withCallingHandlers(
    do.call(qpadm_sweep, c(.feas_args(), list(full_results = FALSE))),
    warning = .mute_qpadm_sweep_noise)
  expect_true("f4_var_rcond" %in% names(res))
  expect_type(res$f4_var_rcond, "double")
  expect_true(all(is.finite(res$f4_var_rcond)))
  expect_false("f4_var_singular_loadings" %in% names(res))
})

test_that("8e: full_results = TRUE adds popdrop list-column", {
  # popdrop is qpadm()'s leave-one-out-per-source table, populated whenever
  # target is non-null. qpadm_sweep always supplies target, so popdrop is
  # always populated. Surfaced as a list-column under full_results = TRUE
  # for method symmetry with weights / rankdrop / f4_var_singular_loadings.
  res = withCallingHandlers(do.call(qpadm_sweep, .feas_args()),
                            warning = .mute_qpadm_sweep_noise)
  expect_true("popdrop" %in% names(res))
  expect_type(res$popdrop, "list")
  expect_false(is.null(res$popdrop[[1]]))
  expect_s3_class(res$popdrop[[1]], "tbl_df")
  # Schema pin: popdrop's column names + types are the public surface that
  # downstream callers reach into. Test 8c pins the analogous contract for
  # f4_var_singular_loadings; same pattern here so a future qpadm refactor
  # that renames or drops a column surfaces as a clear test failure here
  # rather than only at downstream user code. Columns: `pat` (which pop was
  # dropped), `wt` (resulting weight), `p`, `f4rank`, and `feasible`.
  expect_true(all(c("pat", "wt", "p", "f4rank", "feasible") %in%
                  names(res$popdrop[[1]])))

  # And it's dropped under full_results = FALSE alongside the other
  # gated list-columns.
  res_flat = withCallingHandlers(
    do.call(qpadm_sweep, c(.feas_args(), list(full_results = FALSE))),
    warning = .mute_qpadm_sweep_noise)
  expect_false("popdrop" %in% names(res_flat))
})

# ── test 9: n=1 edge case ────────────────────────────────────────────────────

test_that("single target x single source-set x single right-set works without error", {
  res = do.call(qpadm_sweep, .feas_args())

  expect_equal(nrow(res), 1L)
  expect_true(is.numeric(res$p))
  expect_true(!is.na(res$p))
})
