# Regression tests for NA handling in the block-jackknife leave-one-out path.
#
# tidyr::replace_na() dispatches row-wise via vctrs, so it no-ops on matrices /
# 3D arrays that carry partial NAs. The fixes:
#   - jack_mat_stats / jack_mat_stats2: base-R element-wise zeroing of the
#     centered deviation matrix xtau (0 = "no contribution", correct for a
#     pairwise-complete covariance).
#   - est_to_loo_nafix: NA-aware exclusion of missing blocks (delegates to
#     est_to_loo) instead of substituting 0 for a raw f4 block estimate.
#   - qpadm: weight SEs are computed over the surviving (non-missing) leave-one-
#     out blocks, the same set f4_var uses, instead of erroring / returning NaN.

test_that("est_to_loo_nafix matches a plain leave-one-out on complete data", {
  data('example_f2_blocks', package = 'admixtools')
  blk <- example_f2_blocks[1:3, 1:3, 1:25]
  expect_false(anyNA(blk))
  expect_silent(out <- admixtools:::est_to_loo_nafix(blk))
  expect_equal(out, admixtools:::est_to_loo(blk))
  expect_false(anyNA(out))
})

test_that("est_to_loo_nafix excludes missing blocks instead of substituting 0", {
  data('example_f2_blocks', package = 'admixtools')
  blk <- example_f2_blocks[1:3, 1:3, 1:25]
  blk[1, 2, 5] <- NA_real_
  blk[2, 1, 5] <- NA_real_

  # honest warning (the old text falsely claimed "Replacing N NAs with 0!")
  expect_warning(out <- admixtools:::est_to_loo_nafix(blk),
                 regexp = 'missing block entries excluded')

  # NA-aware exclusion == est_to_loo, NOT 0-substitution
  expect_equal(out, admixtools:::est_to_loo(blk))

  # the leave-one-out slice that removes the missing block stays NA in that cell
  # so downstream jackknife code drops it, rather than the 0-substituted value
  expect_true(is.na(out[1, 2, 5]))

  # 0-substitution would have produced a finite, nonzero (wrong) value there
  bl  <- readr::parse_number(dimnames(blk)[[3]])
  tot <- weighted.mean(blk[1, 2, ], bl, na.rm = TRUE)
  wrong_0sub <- tot / (1 - bl[5] / sum(bl))
  expect_true(is.finite(wrong_0sub) && wrong_0sub != 0)
})

test_that("jack_mat_stats2 returns a finite covariance on partial-NA input", {
  set.seed(1)
  loo <- matrix(rnorm(5 * 8, sd = 1e-3), nrow = 5)   # 5 popcombos x 8 blocks
  bl  <- rep(100, 8)

  v_full <- admixtools:::jack_mat_stats2(loo, bl)$var
  expect_true(all(is.finite(v_full)))

  # one block missing for popcomb 2 (some-but-not-all -> a partial-NA row).
  # Pre-fix, replace_na() left the NA and tcrossprod spread it across row/col 2.
  loo_na <- loo; loo_na[2, 3] <- NA_real_
  v_na <- admixtools:::jack_mat_stats2(loo_na, bl)$var
  expect_false(anyNA(v_na))
  expect_true(all(is.finite(v_na)))

  # covariance entries not touching the missing pairing are unchanged
  expect_equal(v_na[1, 4], v_full[1, 4])
})

test_that("qpadm returns finite weight SEs when an f2 block is partly missing", {
  data('example_f2_blocks', package = 'admixtools')
  f2 <- example_f2_blocks
  target <- "Denisova.DG"
  left   <- c("Altai_Neanderthal.DG", "Vindija.DG")
  right  <- c("Chimp.REF", "Mbuti.DG", "Russia_Ust_Ishim.DG", "Switzerland_Bichon.SG")

  # knock out one population pair in one block (f2 is symmetric)
  f2[target, right[1], 5] <- NA_real_
  f2[right[1], target, 5] <- NA_real_

  # Pre-fix this errored (R svd on an NA slice) or returned NaN SEs (cpp). Post-
  # fix the weight covariance uses the surviving-block set, matching f4_var.
  res <- suppressMessages(suppressWarnings(
    qpadm(f2, left = left, right = right, target = target, verbose = FALSE)
  ))
  expect_true(all(is.finite(res$weights$se)))
  expect_true(all(is.finite(res$weights$weight)))
  expect_true(is.finite(res$f4_var_rcond))
})
