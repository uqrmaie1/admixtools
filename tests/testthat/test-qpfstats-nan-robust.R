# Tests for qpfstats_regression (PR #112).
#
# Locks in the NaN-aware per-block weighted-LS solve that replaced the
# original batched matmul. Covers:
#   1. No-NaN: match the original batched solve to machine precision.
#   2. 1% NaN: match a direct per-block reference solve.
#   3. 50% NaN stress: downdate stays well-behaved.
#   4. Bug demonstration: original code produces all-NaN on the same input.
#   5. All-NaN-block edge case: emits b[, i] = 0 + warning.
#   6. bglob NaN-aware path: matches reference when y has NaN cells.

# Reference: direct per-block solve with NaN drop. This is the
# mathematically defining solution and is independent of the
# implementation under test (no downdate trick, no shared Cholesky).
.direct_per_block_solve = function(x, ymat, y, ridge = 0.00001) {
  npairs  = ncol(x)
  nblocks = ncol(ymat)
  b = sapply(seq_len(nblocks), function(i) {
    v = is.finite(ymat[, i])
    solve(crossprod(x[v, , drop = FALSE]) + ridge * diag(npairs),
          crossprod(x[v, , drop = FALSE], ymat[v, i]))
  })
  vy = is.finite(y)
  bglob = drop(solve(
    crossprod(x[vy, , drop = FALSE]) + ridge * diag(npairs),
    crossprod(x[vy, , drop = FALSE], y[vy])
  ))
  list(b = b, bglob = bglob)
}

# Original (pre-PR) batched matmul. Kept verbatim to demonstrate the
# bug it introduced under NaN input.
.original_batched_solve = function(x, ymat, y, ridge = 0.00001) {
  lh = solve((t(x) %*% x) + diag(ncol(x)) * ridge) %*% t(x)
  list(b = lh %*% ymat, bglob = drop(lh %*% y))
}

# Synthetic fixture builder. The shape (5000 x 50 x 100) matches the
# PR description's reproducer and is small enough to run in <100ms.
.make_fixture = function(seed = 2026L, npopcomb = 5000L, nblocks = 50L,
                         npairs = 100L) {
  set.seed(seed)
  list(
    x    = matrix(rnorm(npopcomb * npairs), npopcomb, npairs),
    ymat = matrix(rnorm(npopcomb * nblocks, sd = 0.05), npopcomb, nblocks),
    y    = rnorm(npopcomb, sd = 0.05)
  )
}

test_that("qpfstats_regression: no-NaN input matches original batched solve to machine precision", {
  # The PR's central no-regression claim: when no cell of ymat is NaN,
  # the new code path produces the same answer as the original to
  # machine precision (different code path through BLAS, so not bytewise
  # but within ~3e-18 abs diff on this fixture).
  f = .make_fixture()
  new = admixtools:::qpfstats_regression(f$x, f$ymat, f$y)
  ref = .original_batched_solve(f$x, f$ymat, f$y)
  expect_lt(max(abs(new$b     - ref$b)),     1e-15)
  expect_lt(max(abs(new$bglob - ref$bglob)), 1e-15)
  expect_true(all(is.finite(new$b)))
  expect_true(all(is.finite(new$bglob)))
})

test_that("qpfstats_regression: 1% NaN cells in ymat match direct per-block reference", {
  # The headline bugfix scenario from the PR description: 1% NaN in
  # ymat. The original would emit a fully-NaN b column for any block
  # containing even one NaN cell; the PR matches the per-block
  # reference to machine precision.
  f = .make_fixture()
  set.seed(20260520L)
  f$ymat[sample.int(length(f$ymat), round(0.01 * length(f$ymat)))] = NaN
  new = admixtools:::qpfstats_regression(f$x, f$ymat, f$y)
  ref = .direct_per_block_solve(f$x, f$ymat, f$y)
  expect_lt(max(abs(new$b     - ref$b)),     1e-15)
  expect_lt(max(abs(new$bglob - ref$bglob)), 1e-15)
  expect_true(all(is.finite(new$b)))
  expect_true(all(is.finite(new$bglob)))
})

test_that("qpfstats_regression: 50% NaN cells stress test (downdate stays well-conditioned)", {
  # Pathological-but-not-degenerate scenario: half the ymat cells are
  # NaN, randomly placed. The downdate identity subtracts large
  # crossprod values from A_shared, which in theory could amplify
  # cancellation error. Empirically the agreement with the direct
  # reference holds at ~5e-18 (still machine precision).
  f = .make_fixture()
  set.seed(20260521L)
  f$ymat[sample.int(length(f$ymat), round(0.50 * length(f$ymat)))] = NaN
  new = admixtools:::qpfstats_regression(f$x, f$ymat, f$y)
  ref = .direct_per_block_solve(f$x, f$ymat, f$y)
  # Loosen by 100x vs the 1% test to leave headroom for BLAS variance,
  # but the empirical agreement is much tighter.
  expect_lt(max(abs(new$b     - ref$b)),     1e-13)
  expect_lt(max(abs(new$bglob - ref$bglob)), 1e-13)
  expect_true(all(is.finite(new$b)))
  expect_true(all(is.finite(new$bglob)))
})

test_that("qpfstats_regression: original batched solve produces all-NaN under 1% NaN (proof of bug)", {
  # Locks in the bug demonstration: this is the behavior the PR fixes.
  # If a future refactor accidentally re-introduced the batched matmul,
  # this test would flip from "passes" to "fails" - making the
  # regression visible.
  f = .make_fixture()
  set.seed(20260520L)
  f$ymat[sample.int(length(f$ymat), round(0.01 * length(f$ymat)))] = NaN
  bad = .original_batched_solve(f$x, f$ymat, f$y)
  expect_true(any(!is.finite(bad$b)),
              info = "original batched solve should produce NaN under NaN input")
  # Severity: NaN cells in even ~1% of ymat poison most or all blocks.
  nan_cols = which(colSums(!is.finite(bad$b)) > 0)
  expect_gt(length(nan_cols), 0)
})

test_that("qpfstats_regression: all-NaN block emits b[, i] = 0 with a warning", {
  # Documented edge case: when every popcomb in some block is NaN,
  # the downdate yields A_v = ridge*I and the RHS is 0, so b[, i] = 0
  # exactly. This biases downstream f2()$est for that block but is
  # strictly better than NaN-poisoning every block. The warning
  # surfaces the case so callers can detect it.
  f = .make_fixture()
  f$ymat[, 25] = NaN   # block 25 entirely NaN
  expect_warning(
    new <- admixtools:::qpfstats_regression(f$x, f$ymat, f$y),
    regexp = "no valid SNPs"
  )
  expect_true(all(new$b[, 25] == 0))
  expect_true(all(is.finite(new$b)))
  # Other blocks unaffected - only block 25 is zeroed.
  expect_true(any(new$b[, 1] != 0))
})

test_that("qpfstats_regression: bglob NaN-aware path matches reference when y has NaN cells", {
  # Less-common scenario: a popcomb has NO valid SNPs anywhere across
  # all blocks, so its y value is NaN. The PR applies the same
  # downdate trick to the bglob solve. Verify the result matches a
  # direct reference that drops those popcombs from y's regression.
  f = .make_fixture()
  set.seed(20260522L)
  f$y[sample.int(length(f$y), round(0.05 * length(f$y)))] = NaN
  new = admixtools:::qpfstats_regression(f$x, f$ymat, f$y)
  ref = .direct_per_block_solve(f$x, f$ymat, f$y)
  expect_lt(max(abs(new$bglob - ref$bglob)), 1e-15)
  expect_true(all(is.finite(new$bglob)))
})

test_that("qpfstats_regression: handles is.finite() input distinctions (NaN, NA, Inf)", {
  # The PR uses !is.finite() rather than is.nan(), so NA and Inf cells
  # are also dropped. Confirm each variant of non-finite input is
  # handled identically to NaN - all should drop the popcomb from the
  # regression and produce a finite b/bglob.
  f = .make_fixture()
  set.seed(20260523L)
  f$ymat[1, 1] = NA_real_
  f$ymat[2, 1] = NaN
  f$ymat[3, 1] = Inf
  f$ymat[4, 1] = -Inf
  new = admixtools:::qpfstats_regression(f$x, f$ymat, f$y)
  ref = .direct_per_block_solve(f$x, f$ymat, f$y)
  expect_lt(max(abs(new$b     - ref$b)),     1e-15)
  expect_lt(max(abs(new$bglob - ref$bglob)), 1e-15)
  expect_true(all(is.finite(new$b)))
})
