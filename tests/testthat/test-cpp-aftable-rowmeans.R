# Tests for cpp_aftable_to_dstatnum_rowmeans (PR #107).
#
# Headline contract: the streaming kernel produces output bytewise-
# equivalent to the materialized variant followed by rowMeans(., na.rm=TRUE):
#
#   ref = cpp_aftable_to_dstatnum(...)
#   new = cpp_aftable_to_dstatnum_rowmeans(...)
#   all.equal(rowMeans(ref$num, na.rm = TRUE), new$means)  # TRUE, tol = 0
#   identical(ref$cnt, new$cnt)                            # TRUE
#
# Both variants accept the same (aftable, p1..p4, modelvec, usesnps,
# allsnps, poly_only) inputs. The full equivalence space is the cross of
# (allsnps in {FALSE, TRUE}) x (poly_only in {0, 1, 2}), plus the all-NA
# row edge case where 0/0 must propagate to NaN to match rowMeans.

# The synthetic aftable + popcomb fixture is shared via helper-dstat-fixture.R
# (build_dstat_fixture). This file's original 8-pop / 200-snp / 60-comb /
# 3-model fixture is reproduced byte-for-byte by fixing nmodels = 3 and the
# original seed and forwarding the size args, so the call sites below are
# unchanged and the equivalence assertions see identical data.
.build_fixture = function(npop = 8L, nsnp = 200L, ncomb = 60L,
                          na_frac = 0.05, seed = 20260520L) {
  build_dstat_fixture(npop = npop, nsnp = nsnp, ncomb = ncomb,
                      nmodels = 3L, na_frac = na_frac, seed = seed)
}

test_that("streaming variant matches materialized + rowMeans across (allsnps x poly_only) grid", {
  # Cross of {FALSE, TRUE} x {0, 1, 2} = 6 cells. Each one exercises a
  # different predicate path in the inner loop.
  fx = .build_fixture()
  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      ref = admixtools:::cpp_aftable_to_dstatnum(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po)
      new = admixtools:::cpp_aftable_to_dstatnum_rowmeans(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po)

      # Counts must match exactly.
      expect_identical(as.numeric(ref$cnt), as.numeric(new$cnt),
                       info = sprintf("allsnps=%s poly_only=%d (cnt)", a, po))

      # Means must match bitwise. rowMeans(., na.rm=TRUE) on an all-NA
      # row returns NaN; the streaming variant computes 0/0 = NaN at the
      # same cells. expect_equal with tolerance = 0 enforces bitwise
      # equality but treats NaN == NaN correctly (unlike == in R).
      ref_means = unname(rowMeans(ref$num, na.rm = TRUE))
      expect_equal(ref_means, as.numeric(new$means), tolerance = 0,
                   info = sprintf("allsnps=%s poly_only=%d (means)", a, po))
    }
  }
})

test_that("all-NA rows yield NaN in both variants (matches rowMeans 0/0)", {
  # Force at least one popcomb to be all-NA after the usesnps filter by
  # using indices that point at NA-only rows of aftable. The streaming
  # variant must yield NaN (sum=0, cnt=0, 0/0 = NaN), matching
  # rowMeans(., na.rm=TRUE) on a row that's entirely NA.
  npop = 4L; nsnp = 10L
  aft = matrix(1.0, nrow = npop, ncol = nsnp)
  aft[1, ] = NA_real_   # pop 1 is all NA
  # Popcomb that's all-NA: (1, 1, 1, 1) -> w=x=y=z=NA -> (w-x)*(y-z)=NaN at every i
  p1 = c(1L, 2L); p2 = c(1L, 2L); p3 = c(1L, 2L); p4 = c(1L, 2L)
  modelvec = c(1L, 1L)
  usesnps = matrix(1L, nrow = 1L, ncol = nsnp)

  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      new = admixtools:::cpp_aftable_to_dstatnum_rowmeans(
        aft, p1, p2, p3, p4, modelvec, usesnps, a, po)
      expect_true(is.nan(new$means[1]),
                  info = sprintf("all-NA row, allsnps=%s poly_only=%d", a, po))
      expect_equal(new$cnt[1], 0,
                   info = sprintf("all-NA row cnt, allsnps=%s poly_only=%d", a, po))
      # Pop 2's row is fully populated with the value 1.0. Under
      # poly_only=0 it contributes finite cells (mean = 0 since
      # (w-x)*(y-z) = 0 when all pops are pop 2). Under poly_only>=1
      # with allsnps=TRUE, the polymorphism predicate excludes all SNPs
      # because uni = {1} has length 1 -- which is correct behavior, so
      # the streaming variant correctly returns NaN there too. Only
      # check the non-NaN expectation in the cells where the filter
      # doesn't fire.
      if(po == 0 || !a) {
        expect_false(is.nan(new$means[2]),
                     info = sprintf("non-NA row, allsnps=%s poly_only=%d", a, po))
      }
    }
  }
})

test_that("counts (cnt) are integer-valued and match the materialized variant", {
  # cnt is the count of finite cells contributing to each row's mean.
  # The materialized variant returns cnt as an arma::vec; the streaming
  # variant also returns arma::vec. Both should have integer-valued
  # entries (we never increment by anything other than 1.0).
  fx = .build_fixture(npop = 6L, nsnp = 50L, ncomb = 20L)
  for(a in c(FALSE, TRUE)) {
    for(po in 0:2) {
      ref = admixtools:::cpp_aftable_to_dstatnum(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po)
      new = admixtools:::cpp_aftable_to_dstatnum_rowmeans(
        fx$aft, fx$p1, fx$p2, fx$p3, fx$p4,
        fx$modelvec, fx$usesnps, a, po)
      expect_true(all(as.numeric(new$cnt) == floor(as.numeric(new$cnt))),
                  info = sprintf("cnt integer-valued allsnps=%s poly_only=%d", a, po))
      expect_identical(as.numeric(ref$cnt), as.numeric(new$cnt))
    }
  }
})
