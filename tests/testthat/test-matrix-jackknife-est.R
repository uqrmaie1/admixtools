# Tests for matrix_jackknife_est_full / matrix_jackknife_est (PR #108).
#
# Headline contract: the matrix-vectorized helpers produce per-popcomb
# bias-corrected jackknife estimates that are mathematically equivalent
# to the long-format dplyr chain in f4blockdat_to_f4out:
#
#   f4blockdat %>% group_by(pop1, pop2, pop3, pop4) %>%
#     est_to_loo_dat() %>% jack_dat_stats() %>% pull(est)
#
# The matrix helpers replace this chain in qpfstats to (a) avoid the
# super-linear dplyr group_by overhead at large npopcomb and (b) avoid
# the (npopcomb * nblocks)-row long-format build entirely. Tests
# verify the math equivalence on synthetic input.

# Build a synthetic (numer, cnt) pair plus the matching long-format
# tibble that the dplyr chain expects. Seeded so the test is
# deterministic across runs.
.build_synth_inputs = function(nblocks = 50L, npopcomb = 30L, na_frac = 0.05,
                               seed = 20260520L) {
  withr::with_seed(seed, {
    numer = matrix(rnorm(nblocks * npopcomb, sd = 0.01),
                   nrow = nblocks, ncol = npopcomb)
    cnt   = matrix(sample(50:200, nblocks * npopcomb, replace = TRUE),
                   nrow = nblocks, ncol = npopcomb)
    # Sprinkle NAs in numer to exercise the !is.finite(est) path in
    # est_to_loo_dat and jack_dat_stats.
    numer[sample.int(length(numer),
                     size = round(na_frac * length(numer)))] = NA_real_

    # Long-format tibble matching the shape f4blockdat_from_geno would
    # emit. Pop labels are synthesized to be unique per popcomb so the
    # dplyr group_by partitions exactly the columns of numer/cnt.
    popcomb = tibble::tibble(
      pop1 = sprintf("P%dA", seq_len(npopcomb)),
      pop2 = sprintf("P%dB", seq_len(npopcomb)),
      pop3 = sprintf("P%dC", seq_len(npopcomb)),
      pop4 = sprintf("P%dD", seq_len(npopcomb))
    )
    f4bd = popcomb %>%
      tidyr::expand_grid(block = seq_len(nblocks)) %>%
      dplyr::mutate(
        est = c(numer), n = c(cnt),
        length = sample(1000:5000, dplyr::n(), replace = TRUE),
        est = dplyr::if_else(is.finite(est) & n > 0, est, NA_real_))
  })
  list(numer = numer, cnt = cnt, f4bd = f4bd, popcomb = popcomb,
       nblocks = nblocks, npopcomb = npopcomb)
}

test_that("matrix_jackknife_est_full matches dplyr f4blockdat_to_f4out on synthetic input", {
  fx = .build_synth_inputs()

  ref = admixtools:::f4blockdat_to_f4out(fx$f4bd, FALSE)$est
  new = admixtools:::matrix_jackknife_est_full(fx$numer, fx$cnt)

  # Values should agree within FP summation-order noise. Observed max
  # abs diff on this fixture is ~1 ULP (2.2e-16); tolerance 1e-13
  # catches regressions that would introduce >100 ULPs while still
  # tolerating natural cross-platform FP variation.
  expect_equal(unname(ref), unname(new), tolerance = 1e-13)
})

test_that("matrix_jackknife_est_full handles high-NA-fraction input (30% NA)", {
  # Production f-stat data on ancient-DNA panels can have 30-60% NA per
  # popcomb after the per-block usesnps mask. Stress-test the matrix
  # path's NA propagation at that scale by inflating na_frac. The
  # equivalence vs dplyr should hold to the same FP-rounding tolerance
  # as at low-NA scale -- if the matrix path mishandles NA propagation,
  # it would show up here as a substantial diff.
  fx = .build_synth_inputs(nblocks = 100L, npopcomb = 40L, na_frac = 0.30)

  ref = admixtools:::f4blockdat_to_f4out(fx$f4bd, FALSE)$est
  new = admixtools:::matrix_jackknife_est_full(fx$numer, fx$cnt)

  # Same tolerance as the low-NA test -- high NA fraction shouldn't
  # amplify rounding beyond 1 ULP at this fixture size.
  expect_equal(unname(ref), unname(new), tolerance = 1e-13)

  # Sanity: with 30% NA per cell across 100 blocks, some popcombs will
  # have <2 valid blocks -- those return NA in both paths. NA pattern
  # must match between matrix and dplyr.
  expect_identical(is.na(ref), is.na(new))
})

test_that("matrix_jackknife_est (chunked) matches matrix_jackknife_est_full (unchunked) bitwise", {
  # Chunking is column-independent: every operation in
  # matrix_jackknife_est_full is a per-column sweep or colSums.
  # Chunked output must therefore be byte-identical to unchunked,
  # not just numerically equivalent.
  fx = .build_synth_inputs(nblocks = 40L, npopcomb = 25L)

  unchunked = admixtools:::matrix_jackknife_est_full(fx$numer, fx$cnt)
  chunked   = admixtools:::matrix_jackknife_est(fx$numer, fx$cnt, chunk_size = 7L)

  expect_identical(unchunked, chunked)
})

test_that("matrix_jackknife_est: NULL chunk_size auto-sizes and matches full path", {
  # The wrapper with chunk_size = NULL picks a chunk size that keeps
  # peak memory bounded. For small inputs the auto-sized chunk is
  # >= npopcomb and the wrapper bypasses the chunking machinery
  # entirely, returning matrix_jackknife_est_full's output directly.
  fx = .build_synth_inputs(nblocks = 30L, npopcomb = 20L)
  auto = admixtools:::matrix_jackknife_est(fx$numer, fx$cnt, chunk_size = NULL)
  full = admixtools:::matrix_jackknife_est_full(fx$numer, fx$cnt)
  expect_identical(auto, full)
})

test_that("all-NA popcomb returns NA in matrix_jackknife_est_full", {
  # Construct an input where one popcomb has every block's est = NA.
  # That row's rel_bl denominator is sum(n) > 0 but its tot is NaN
  # (0/0), and the final est should be NA (not crash, not Inf).
  fx = .build_synth_inputs(nblocks = 20L, npopcomb = 5L)
  fx$numer[, 3L] = NA_real_   # popcomb 3 has all-NA blocks

  new = admixtools:::matrix_jackknife_est_full(fx$numer, fx$cnt)
  expect_true(is.na(new[3]))
  # Other popcombs still have finite output.
  expect_true(all(!is.na(new[-3])))
})

test_that("single-valid-block popcomb matches dplyr (both return the surviving block's est)", {
  # Popcomb 1 has only block 1 valid (blocks 2-3 are NA in numer). The
  # math reduces to: tot = weighted.mean(est, n) = est[1]; after the
  # is.finite(loo) filter only block 1 survives; n_finite = 1; the
  # final formula yields est[1]. Verified against dplyr in the same
  # test (NOT a hard-coded NA assumption).
  numer = matrix(c(1.0, NA_real_, NA_real_,    # popcomb 1: 1 valid block
                   2.0, 3.0,      4.0),         # popcomb 2: 3 valid blocks
                 nrow = 3, ncol = 2)
  cnt   = matrix(rep(100, 6), nrow = 3, ncol = 2)

  new = admixtools:::matrix_jackknife_est_full(numer, cnt)

  # Build the matching long-format tibble for the dplyr reference.
  popcomb = tibble::tibble(pop1 = c("A1", "A2"), pop2 = c("B1", "B2"),
                           pop3 = c("C1", "C2"), pop4 = c("D1", "D2"))
  f4bd = popcomb %>%
    tidyr::expand_grid(block = 1:3) %>%
    dplyr::mutate(est = c(numer), n = c(cnt), length = 100L,
                  est = dplyr::if_else(is.finite(est) & n > 0, est, NA_real_))
  ref = admixtools:::f4blockdat_to_f4out(f4bd, FALSE)$est

  expect_equal(unname(ref), unname(new), tolerance = 1e-12)
  # Popcomb 1 should equal the surviving block's est (1.0).
  expect_equal(new[1], 1.0, tolerance = 1e-12)
})

test_that(".matrix_jackknife_chunk_size returns a positive integer scaling inversely with nblocks", {
  # ~5 GB target / 8 bytes per double = 625e6 doubles per intermediate;
  # divided by nblocks gives chunk size in popcombs.
  expect_gte(admixtools:::.matrix_jackknife_chunk_size(100L),
             admixtools:::.matrix_jackknife_chunk_size(1000L))
  expect_true(is.integer(admixtools:::.matrix_jackknife_chunk_size(100L)))
  expect_gte(admixtools:::.matrix_jackknife_chunk_size(1000L), 1L)
  # At pathologically large nblocks the chunk size must still be >= 1
  # (the helper clamps via max(1L, ...)).
  expect_gte(admixtools:::.matrix_jackknife_chunk_size(1e9), 1L)
})
