# Regression tests for issue #89 / PR fix:
# qpadm(prefix, left, right, target = NULL, return_f4 = TRUE) used to crash
# with "Column `pop1` doesn't exist" because expand_grid(pop1 = NULL, ...)
# silently drops the pop1 column. The fix replaces `pop1 = target` with
# `pop1 = left[1]`, which is equivalent when target != NULL (the line above
# puts target at left[1]) and well-defined when target == NULL (the canonical
# anchor for the qpwave-style f4 parameterization).
#
# Strongest correctness argument: within the geno-prefix branch, the
# `return_f4 = TRUE` path (the one fixed by this PR) and the
# `return_f4 = FALSE` path (the always-working sibling) MUST produce
# bit-identical rankdrop output. They consume the same genotypes, use the
# same `f4blockdat_from_geno` block reader, apply the same SNP filters,
# and produce the same per-block f4 statistics -- the only difference
# is that the return_f4 = TRUE branch computes the full f4 matrix and
# filters it down to right[1]-anchored pairs (so it can expose `out$f4`
# for the caller), while return_f4 = FALSE builds the right[1]-anchored
# popcombs directly. The downstream f4blocks_to_f4stats / drop_ranks
# pipeline is identical.

.individual_pop_fixture = function(dir) {
  # build_pfile_fixture writes FID = "popA" / "popB" with 3 inds each. For
  # qpwave-style tests we need >=4 distinct populations, so we rewrite the
  # .fam to put each individual in its own population (FID = IID).
  fix = build_pfile_fixture(dir, with_fid = TRUE)
  pref = fix$bed_pref
  if(is.na(pref)) return(NA_character_)

  fam_path = paste0(pref, ".fam")
  fam = read.table(fam_path, stringsAsFactors = FALSE)
  fam$V1 = fam$V2
  write.table(fam, fam_path, sep = " ", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  pref
}

test_that("qpadm(prefix, target = NULL, return_f4 = TRUE) no longer crashes (#89)", {
  # Pre-fix this errored with `Column 'pop1' doesn't exist` because
  # expand_grid(pop1 = NULL, ...) silently dropped the pop1 column from
  # popcombs, and downstream summarize() / select() couldn't find it.
  withr::with_tempdir({
    pref = .individual_pop_fixture(getwd())
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    expect_no_error(
      suppressMessages(suppressWarnings(
        qpadm(pref,
              left   = c("popA_1", "popA_2"),
              right  = c("popB_1", "popB_2"),
              target = NULL,
              return_f4 = TRUE,
              blgsize = 100,
              verbose = FALSE)
      ))
    )
  })
})

test_that("qpadm(target = NULL, return_f4 = TRUE) returns the qpwave-shaped result", {
  # When target is NULL, qpadm's downstream branching skips the
  # target-dependent computations (weights, fitted_f4, popdrop). The
  # returned list should contain only `rankdrop` and `f4` -- no `weights`
  # and no `popdrop`. This pins the documented qpwave-via-qpadm contract.
  withr::with_tempdir({
    pref = .individual_pop_fixture(getwd())
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    res = suppressMessages(suppressWarnings(
      qpadm(pref,
            left   = c("popA_1", "popA_2"),
            right  = c("popB_1", "popB_2"),
            target = NULL,
            return_f4 = TRUE,
            blgsize = 100,
            verbose = FALSE)
    ))

    expect_true("rankdrop" %in% names(res))
    expect_true("f4"       %in% names(res))
    expect_false("weights" %in% names(res))
    expect_false("popdrop" %in% names(res))

    # rankdrop has the documented columns and finite values
    expect_true(all(c("f4rank", "dof", "chisq", "p") %in% names(res$rankdrop)))
    expect_true(all(is.finite(res$rankdrop$chisq)))
    expect_true(all(is.finite(res$rankdrop$p)))
  })
})

test_that("qpadm(target = NULL): return_f4 = TRUE matches return_f4 = FALSE bit-exactly (#89 correctness)", {
  # This is the strongest correctness check for the fix. The fixed branch
  # (return_f4 = TRUE) and the always-working sibling (return_f4 = FALSE)
  # consume the same genotypes through the same f4blockdat_from_geno reader.
  # They differ only in how `popcombs` is constructed:
  #
  #   return_f4 = TRUE:
  #     popcombs = expand_grid(pop1 = left[1], pop2 = left[-1],
  #                            pop3 = right,   pop4 = right) %>%
  #                filter(pop3 != pop4)
  #     -> compute all pairs
  #     -> filter to pop3 == right[1] (the right[1]-anchored f4 set)
  #
  #   return_f4 = FALSE:
  #     popcombs = expand_grid(pop1 = left[1], pop2 = left[-1],
  #                            pop3 = right[1], pop4 = right[-1])
  #     -> compute the right[1]-anchored f4 set directly
  #
  # Both reach the same set of (pop1, pop2, pop3, pop4) tuples in
  # f4blockdat, the same per-block f4 statistics, the same f4blocks 3D
  # array, the same f4blocks_to_f4stats output, and the same drop_ranks
  # output. So rankdrop must be BIT-IDENTICAL between the two calls.
  #
  # If this test ever fails, something downstream of `popcombs`
  # construction has gained a target-vs-NULL-or-return_f4-vs-FALSE
  # branch that's diverged the two paths.
  withr::with_tempdir({
    pref = .individual_pop_fixture(getwd())
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    left  = c("popA_1", "popA_2")
    right = c("popB_1", "popB_2")

    res_t = suppressMessages(suppressWarnings(
      qpadm(pref, left = left, right = right, target = NULL,
            return_f4 = TRUE, blgsize = 100, verbose = FALSE)
    ))
    res_f = suppressMessages(suppressWarnings(
      qpadm(pref, left = left, right = right, target = NULL,
            return_f4 = FALSE, blgsize = 100, verbose = FALSE)
    ))

    expect_identical(res_t$rankdrop$f4rank, res_f$rankdrop$f4rank)
    expect_identical(res_t$rankdrop$dof,    res_f$rankdrop$dof)
    expect_identical(res_t$rankdrop$chisq,  res_f$rankdrop$chisq)
    expect_identical(res_t$rankdrop$p,      res_f$rankdrop$p)
  })
})

test_that("qpadm(target = NULL, return_f4 = TRUE, allsnps = TRUE) runs on a geno prefix (#89)", {
  # The reporter's original call shape (and #69's documented use case):
  # qpwave-with-allsnps via qpadm(target = NULL, return_f4 = TRUE,
  # allsnps = TRUE). Pre-fix this errored with `Column 'pop1' doesn't exist`
  # before any allsnps-specific logic ran.
  withr::with_tempdir({
    pref = .individual_pop_fixture(getwd())
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    res = suppressMessages(suppressWarnings(
      qpadm(pref,
            left   = c("popA_1", "popA_2"),
            right  = c("popB_1", "popB_2"),
            target = NULL,
            return_f4 = TRUE,
            allsnps = TRUE,
            blgsize = 100,
            verbose = FALSE)
    ))
    expect_true("rankdrop" %in% names(res))
    expect_true(all(is.finite(res$rankdrop$chisq)))
  })
})

test_that("qpadm(target != NULL): the fix doesn't change behavior on the working path (regression)", {
  # The fix changes `pop1 = target` to `pop1 = left[1]` in qpadm.R's
  # return_f4 = TRUE branch. When target != NULL, the line right before
  # that does `left = c(target, setdiff(left, target))`, putting target
  # at left[1]. So `pop1 = target` and `pop1 = left[1]` were equivalent
  # for the target != NULL case before the fix and must REMAIN
  # equivalent. Pin that contract with an explicit before-and-after
  # comparison against return_f4 = FALSE on the same target != NULL call.
  withr::with_tempdir({
    pref = .individual_pop_fixture(getwd())
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    target_p = "popA_1"
    left   = c("popA_2", "popA_3")
    right  = c("popB_1", "popB_2", "popB_3")

    res_t = suppressMessages(suppressWarnings(
      qpadm(pref, left = left, right = right, target = target_p,
            return_f4 = TRUE, blgsize = 100, verbose = FALSE)
    ))
    res_f = suppressMessages(suppressWarnings(
      qpadm(pref, left = left, right = right, target = target_p,
            return_f4 = FALSE, blgsize = 100, verbose = FALSE)
    ))

    # The two return_f4 modes must agree on rankdrop for target != NULL too.
    expect_identical(res_t$rankdrop$chisq, res_f$rankdrop$chisq)
    expect_identical(res_t$rankdrop$p,     res_f$rankdrop$p)
    expect_identical(res_t$rankdrop$dof,   res_f$rankdrop$dof)
  })
})
