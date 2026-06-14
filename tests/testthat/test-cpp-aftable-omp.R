# Tests for the OpenMP parallelization of the cpp_aftable_to_dstat* and
# cpp_gmat_to_aftable kernels.
#
# The kernels take an explicit `nthreads` argument, resolved on the R side in
# the main process (see f4blockdat_from_geno / omp_thread_budget). The thread
# count is therefore set directly here, not via OMP_NUM_THREADS.
#
# Two contracts:
#   1. Output for nthreads = 1 (serial) must be bit-identical to nthreads = N.
#      Each popcomb row is reduced by a single thread with no cross-thread
#      accumulation, so equality is exact (not tolerance-based) and holds by
#      construction; the test guards that the serial and parallel code paths
#      really produce the same bytes. N = .NT = 2 keeps the test within CRAN's
#      2-core check limit; on a build without OpenMP both calls run serially
#      and the comparison still holds.
#   2. The serial validation pass surfaces bad p1..p4 / modelvec / popvec /
#      usesnps inputs as clean Rcpp::stop errors BEFORE the parallel region
#      (where Rcpp::stop is not callable).
#
# The fixture lives in helper-dstat-fixture.R (shared with
# test-cpp-aftable-rowmeans.R). Defaults give npopcomb = 256 so
# schedule(guided) produces several chunks per thread.

.NT = 2L  # parallel thread count for the equivalence tests (CRAN 2-core safe)

# Assert a kernel's output is bit-identical at 1 thread and .NT threads.
# `kernel` is a function of nthreads; `accessors` names the list elements to
# compare (NULL compares the whole return value, e.g. the dstatden matrix).
expect_threadcount_invariant = function(kernel, accessors = NULL, info = NULL) {
  r1 = kernel(1L)
  rN = kernel(.NT)
  if (is.null(accessors)) {
    expect_identical(r1, rN, info = info)
  } else {
    for (acc in accessors) {
      expect_identical(r1[[acc]], rN[[acc]], info = paste(info, acc))
    }
  }
}

test_that("cpp_aftable_to_dstatnum_rowmeans: serial == parallel (bitwise)", {
  fx = build_dstat_fixture()
  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      expect_threadcount_invariant(
        function(nt) admixtools:::cpp_aftable_to_dstatnum_rowmeans(
          fx$aft, fx$p1, fx$p2, fx$p3, fx$p4, fx$modelvec, fx$usesnps, a, po, nt),
        accessors = c("means", "cnt"),
        info = sprintf("allsnps=%s poly_only=%d", a, po))
    }
  }
})

test_that("cpp_aftable_to_dstatnum: serial == parallel (bitwise)", {
  fx = build_dstat_fixture()
  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      expect_threadcount_invariant(
        function(nt) admixtools:::cpp_aftable_to_dstatnum(
          fx$aft, fx$p1, fx$p2, fx$p3, fx$p4, fx$modelvec, fx$usesnps, a, po, nt),
        accessors = c("num", "cnt"),
        info = sprintf("allsnps=%s poly_only=%d", a, po))
    }
  }
})

test_that("cpp_aftable_to_dstatden: serial == parallel (bitwise)", {
  fx = build_dstat_fixture()
  for(a in c(FALSE, TRUE)) {
    expect_threadcount_invariant(
      function(nt) admixtools:::cpp_aftable_to_dstatden(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4, fx$modelvec, fx$usesnps, a, 0L, nt),
      info = sprintf("allsnps=%s", a))
  }
})

test_that("cpp_gmat_to_aftable: serial == parallel (bitwise)", {
  withr::with_seed(20260525L, {
    nind = 200L; nsnp = 500L; npop = 8L
    # Integer dosages in {0, 1, 2, NA_INTEGER}: this is what the readers hand
    # the kernel (cpp_read_plink -> IntegerMatrix, pgenlibr::ReadIntList ->
    # integer), and what cpp_gmat_to_aftable's arma::imat signature takes
    # natively. Sampling from a double vector would instead exercise R's
    # double->integer coercion fallback at the boundary, not the native path.
    gmat = matrix(sample(c(0L, 1L, 2L, NA_integer_), nind * nsnp, replace = TRUE,
                         prob = c(0.3, 0.4, 0.25, 0.05)),
                  nrow = nind, ncol = nsnp)
    popvec = sample.int(npop, nind, replace = TRUE)
  })
  expect_threadcount_invariant(
    function(nt) admixtools:::cpp_gmat_to_aftable(gmat, popvec, nt))
})

test_that("dstat kernels Rcpp::stop on NA in p1..p4 instead of silent OOB", {
  fx = build_dstat_fixture(ncomb = 16L)
  fx$p1[3] = NA_real_
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L),
    "p1/p2/p3/p4 contains NA / NaN / Inf"
  )
})

test_that("dstat kernels Rcpp::stop on out-of-range p1..p4 instead of silent OOB", {
  fx = build_dstat_fixture(ncomb = 16L)
  npop_max = nrow(fx$aft)
  fx$p1[5] = as.numeric(npop_max + 5L)  # out of range
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L),
    "index out of range"
  )
})

test_that("dstat kernels Rcpp::stop on finite-but-huge p index (range-checked as double)", {
  # A value above INT_MAX is finite (passes R_FINITE) but (int) casting it is
  # undefined behavior. The range check operates on the double BEFORE the cast,
  # so it is caught rather than slipping a garbage index into the OMP region.
  fx = build_dstat_fixture(ncomb = 16L)
  fx$p1[2] = 1e18
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, fx$usesnps, TRUE, 0L),
    "index out of range"
  )
})

test_that("dstat kernels Rcpp::stop on NA modelvec when allsnps=FALSE", {
  fx = build_dstat_fixture(ncomb = 16L)
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
  fx = build_dstat_fixture(ncomb = 16L)
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

test_that("dstat kernels Rcpp::stop on usesnps with fewer columns than aftable (allsnps=FALSE)", {
  # usesnps(m, i) is read for every i in [0, nsnp) inside the parallel region.
  # If usesnps has fewer columns than aftable, that read is out of bounds and
  # would throw across the OMP boundary (std::terminate). The validation pass
  # catches it serially first.
  fx = build_dstat_fixture(ncomb = 16L)
  bad_usesnps = fx$usesnps[, 1:10, drop = FALSE]  # 10 cols vs nsnp = 1000
  expect_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, bad_usesnps, FALSE, 0L),
    "usesnps has .* columns"
  )
  # allsnps=TRUE never reads usesnps, so a too-small usesnps is accepted.
  expect_no_error(
    admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
      fx$modelvec, bad_usesnps, TRUE, 0L)
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

test_that("poly_only_keep keeps an Inf allele frequency, matching R's na_omit(unique())", {
  # poly_only_keep uses ISNAN (drops NA / NaN, keeps +/-Inf), matching R's
  # na_omit(unique(c(w, x, y, z))). With w = +Inf and the other three equal,
  # kept = {Inf, 0.5, 0.5, 0.5} has > 1 distinct value, so the SNP is KEPT and
  # num = (Inf - 0.5) * (0.5 - 0.5) = NaN. is.nan(NaN) is TRUE, which
  # distinguishes a kept-then-computed cell from a dropped one (NA_REAL, where
  # is.nan is FALSE).
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
  expect_true(is.nan(res$num[1, 1]),
              info = "poly_only_keep with ISNAN must KEEP the Inf+three-equal popcomb")
})
