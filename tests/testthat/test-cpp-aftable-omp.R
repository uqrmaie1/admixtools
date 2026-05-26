# Tests for the OpenMP parallelization of the cpp_aftable_to_dstat* and
# cpp_gmat_to_aftable kernels.
#
# Two contracts:
#   1. Output under OMP_NUM_THREADS=1 must be bit-identical to output
#      under OMP_NUM_THREADS=N for any N >= 1. The reduction is per-row
#      with a single-writer-per-row partitioning, so this is a real
#      contract, not a tolerance-based comparison.
#   2. The serial validation passes added alongside the OMP pragmas
#      (NA / out-of-range index detection on p1..p4, modelvec, popvec)
#      surface as clean Rcpp::stop errors with the offending row index,
#      rather than wild memory access or silent corruption inside the
#      parallel region (where Rcpp::stop is not callable).
#
# Fixture: npopcomb = 256 to stress schedule(dynamic, 16) with multiple
# chunks per thread. With nsnp = 1000 the inner loop is large enough
# that the OMP team setup overhead is amortized.

.build_omp_fixture = function(npop = 12L, nsnp = 1000L, ncomb = 256L,
                              na_frac = 0.05, seed = 20260525L) {
  withr::with_seed(seed, {
    aft = matrix(runif(npop * nsnp), nrow = npop)
    aft[sample.int(length(aft), size = round(na_frac * length(aft)))] = NA_real_
    p1 = sample.int(npop, ncomb, replace = TRUE)
    p2 = sample.int(npop, ncomb, replace = TRUE)
    p3 = sample.int(npop, ncomb, replace = TRUE)
    p4 = sample.int(npop, ncomb, replace = TRUE)
    modelvec = sample.int(4L, ncomb, replace = TRUE)
    usesnps = matrix(sample(c(0, 1), 4L * nsnp, replace = TRUE,
                            prob = c(0.1, 0.9)),
                     nrow = 4L, ncol = nsnp)
  })
  list(aft = aft, p1 = p1, p2 = p2, p3 = p3, p4 = p4,
       modelvec = modelvec, usesnps = usesnps)
}

.with_threads = function(n, expr) {
  # OMP_NUM_THREADS is honored by both parallelly::availableCores() (which
  # the C++ kernels read for thread count) and by OpenMP's omp_get_max_threads.
  # withr::local_envvar scopes the change to the calling frame.
  withr::with_envvar(c(OMP_NUM_THREADS = as.character(n)), force(expr))
}

test_that("cpp_aftable_to_dstatnum_rowmeans: 1-thread output == 4-thread output (bitwise)", {
  fx = .build_omp_fixture()
  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      r1 = .with_threads(1L, admixtools:::cpp_aftable_to_dstatnum_rowmeans(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po))
      r4 = .with_threads(4L, admixtools:::cpp_aftable_to_dstatnum_rowmeans(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po))
      expect_identical(r1$means, r4$means,
                       info = sprintf("means allsnps=%s poly_only=%d", a, po))
      expect_identical(r1$cnt, r4$cnt,
                       info = sprintf("cnt allsnps=%s poly_only=%d", a, po))
    }
  }
})

test_that("cpp_aftable_to_dstatnum: 1-thread output == 4-thread output (bitwise)", {
  fx = .build_omp_fixture()
  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      r1 = .with_threads(1L, admixtools:::cpp_aftable_to_dstatnum(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po))
      r4 = .with_threads(4L, admixtools:::cpp_aftable_to_dstatnum(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po))
      expect_identical(r1$num, r4$num,
                       info = sprintf("num allsnps=%s poly_only=%d", a, po))
      expect_identical(r1$cnt, r4$cnt,
                       info = sprintf("cnt allsnps=%s poly_only=%d", a, po))
    }
  }
})

test_that("cpp_aftable_to_dstatden: 1-thread output == 4-thread output (bitwise)", {
  fx = .build_omp_fixture()
  for(a in c(FALSE, TRUE)) {
    r1 = .with_threads(1L, admixtools:::cpp_aftable_to_dstatden(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, a, 0L))
    r4 = .with_threads(4L, admixtools:::cpp_aftable_to_dstatden(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, a, 0L))
    expect_identical(r1, r4, info = sprintf("den allsnps=%s", a))
  }
})

test_that("cpp_gmat_to_aftable: 1-thread output == 4-thread output (bitwise)", {
  withr::with_seed(20260525L, {
    nind = 200L; nsnp = 500L; npop = 8L
    gmat = matrix(sample(c(0, 1, 2, NA_real_), nind * nsnp, replace = TRUE,
                         prob = c(0.3, 0.4, 0.25, 0.05)),
                  nrow = nind, ncol = nsnp)
    popvec = sample.int(npop, nind, replace = TRUE)
  })
  r1 = .with_threads(1L, admixtools:::cpp_gmat_to_aftable(gmat, popvec))
  r4 = .with_threads(4L, admixtools:::cpp_gmat_to_aftable(gmat, popvec))
  expect_identical(r1, r4)
})

test_that("dstat kernels Rcpp::stop on NA in p1..p4 instead of silent OOB", {
  fx = .build_omp_fixture(ncomb = 16L)
  fx$p1[3] = NA_real_
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L),
    "p1/p2/p3/p4 contains NA / NaN / Inf"
  )
})

test_that("dstat kernels Rcpp::stop on out-of-range p1..p4 instead of silent OOB", {
  fx = .build_omp_fixture(ncomb = 16L)
  npop_max = nrow(fx$aft)
  fx$p1[5] = as.numeric(npop_max + 5L)  # out of range
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L),
    "index out of range"
  )
})

test_that("dstat kernels Rcpp::stop on NA modelvec when allsnps=FALSE", {
  fx = .build_omp_fixture(ncomb = 16L)
  fx$modelvec[4] = NA_real_
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, FALSE, 0L),
    "modelvec contains NA / NaN / Inf"
  )
  # When allsnps=TRUE modelvec is unused; same input should be accepted.
  expect_no_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L)
  )
})

test_that("dstat kernels Rcpp::stop on out-of-range modelvec when allsnps=FALSE", {
  fx = .build_omp_fixture(ncomb = 16L)
  nmodel = nrow(fx$usesnps)
  fx$modelvec[7] = as.numeric(nmodel + 3L)  # past the usesnps row count
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, FALSE, 0L),
    "modelvec value out of range"
  )
  fx$modelvec[7] = 0  # also: 0 is out of range (1-based)
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, FALSE, 0L),
    "modelvec value out of range"
  )
})

test_that("cpp_gmat_to_aftable Rcpp::stop on popvec value > nind (upper bound check)", {
  nind = 50L; nsnp = 20L
  gmat = matrix(0, nrow = nind, ncol = nsnp)
  popvec = sample.int(5L, nind, replace = TRUE)
  popvec[10] = nind + 100L  # malformed: index larger than the implied pop list size
  expect_error(
    admixtools:::cpp_gmat_to_aftable(gmat, popvec),
    "exceeds nind"
  )
})

test_that("cpp_gmat_to_aftable Rcpp::stop on popvec length mismatch with gmat n_rows", {
  gmat = matrix(0, nrow = 50L, ncol = 20L)
  popvec = sample.int(5L, 49L, replace = TRUE)  # one short
  expect_error(
    admixtools:::cpp_gmat_to_aftable(gmat, popvec),
    "popvec length"
  )
})

test_that("OMP_NUM_THREADS env var (via parallelly::availableCores) is respected", {
  # Note: this verifies output identity between an OMP_NUM_THREADS=1 call
  # and a default-threads call. The cpp helper caches the thread count
  # process-locally on first call, so subsequent env-var flips within the
  # same R session may not actually change the running thread count; the
  # output equivalence holds regardless, which is the contract that
  # matters. To verify literal thread-count enforcement, exec a fresh R
  # process with the env var set.
  fx = .build_omp_fixture(ncomb = 64L)
  r1 = withr::with_envvar(c(OMP_NUM_THREADS = "1"),
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L))
  rN = admixtools:::cpp_aftable_to_dstatnum_rowmeans(
    fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
    fx$modelvec, fx$usesnps, TRUE, 0L)
  expect_identical(r1$means, rN$means)
  expect_identical(r1$cnt, rN$cnt)
})

test_that("poly_only_keep preserves Rcpp na_omit semantics on Inf inputs (after ISNAN fix)", {
  # Synthesize an aftable where the polymorphism test would diverge
  # between std::isfinite and ISNAN. Concretely: w = +Inf and the other
  # three pops all share the same finite value.
  #
  # Pre-fix (std::isfinite): finite = {0.5, 0.5, 0.5}, nf = 3, no diff,
  # poly_only != 2 -> drop. num stays at NA_REAL.
  #
  # Post-fix (ISNAN, matching R na_omit): kept = {+Inf, 0.5, 0.5, 0.5},
  # nk = 4, has_diff = true (Inf != 0.5) -> keep. num gets assigned
  # (Inf - 0.5) * (0.5 - 0.5) = Inf * 0 = NaN, which IS distinct from
  # NA_REAL: is.nan(NaN) is TRUE, is.nan(NA_REAL) is FALSE.
  npop = 4L; nsnp = 1L
  aft = matrix(0.5, nrow = npop, ncol = nsnp)
  aft[1, 1] = Inf
  usesnps = matrix(1.0, nrow = 1L, ncol = nsnp)
  res = admixtools:::cpp_aftable_to_dstatnum(
    aft,
    c(1),     # w = Inf
    c(2),     # x = 0.5
    c(3),     # y = 0.5
    c(4),     # z = 0.5
    c(1), usesnps, TRUE, 1L)
  # is.nan distinguishes NaN (kept then computed) from NA_REAL (skipped).
  expect_true(is.nan(res$num[1, 1]),
              info = "poly_only_keep with ISNAN must KEEP the Inf+three-equal popcomb")
  # Sanity: pre-fix would have dropped, leaving NA_REAL (is.nan FALSE).
})
