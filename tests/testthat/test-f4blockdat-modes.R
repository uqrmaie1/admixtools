# Tests for f4blockdat_from_geno's branches that the streaming kernel
# doesn't directly touch but that need regression coverage:
#
#   1. snpwt non-NULL path: uses the materialized cpp_aftable_to_dstatnum
#      kernel + R-side `t(t(num) * snpwt)` column scaling + rowMeans, NOT
#      the streaming variant. We want to confirm this branch still works
#      and that snpwt = rep(1, nsnp) (the identity weighting) produces
#      output bytewise-identical to snpwt = NULL (which fires the
#      streaming variant). If those diverge, either the snpwt path is
#      broken or the streaming/materialized equivalence is broken.
#
#   2. f4mode = FALSE: emits the dstatden denominator alongside the f4
#      numerator. The streaming variant only replaces the numerator
#      path; the denominator still uses cpp_aftable_to_dstatden. We
#      smoke-test that the !f4mode branch completes and returns the
#      expected (est, n, den) column shape.

.modes_popcombs = function() {
  tibble::tibble(pop1 = "popA", pop2 = "popB",
                 pop3 = "popA", pop4 = "popB")
}

test_that("f4blockdat_from_geno: snpwt = rep(1, nsnp) (materialized identity) matches snpwt = NULL (streaming)", {
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    # snpwt is per-SNP weights over ALL SNPs in the geno file (not the
    # autosomal subset). For our fixture, that's just fix$nsnp values.
    # All ones means the materialized variant's column scaling
    # t(t(num) * 1) is a no-op, so the resulting rowMeans should match
    # exactly what the streaming variant produces directly.
    snpwt_ones = rep(1.0, fix$nsnp)

    res_stream = suppressMessages(suppressWarnings(
      admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref, popcombs = .modes_popcombs(),
        auto_only = FALSE, blgsize = 500L, f4mode = TRUE,
        allsnps = FALSE, snpwt = NULL, verbose = FALSE)))
    res_matwt = suppressMessages(suppressWarnings(
      admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref, popcombs = .modes_popcombs(),
        auto_only = FALSE, blgsize = 500L, f4mode = TRUE,
        allsnps = FALSE, snpwt = snpwt_ones, verbose = FALSE)))

    # Both paths run end-to-end and produce the same shape.
    expect_identical(dim(res_stream), dim(res_matwt))
    expect_identical(names(res_stream), names(res_matwt))

    # est and n columns should be bytewise-identical when snpwt = 1 acts
    # as a no-op. expect_equal with tolerance = 0 because est values are
    # produced by identical floating-point operations on both paths.
    expect_equal(res_stream$est, res_matwt$est, tolerance = 0)
    expect_equal(res_stream$n,   res_matwt$n,   tolerance = 0)
  })
})

test_that("f4blockdat_from_geno: snpwt path scales correctly with snpwt = 2 * 1s", {
  # Per-SNP weighting of 2.0 across all SNPs: the materialized variant
  # computes t(t(num) * 2) which doubles every cell. rowMeans then
  # returns 2x the unweighted mean (since the weighting affects the
  # numerator but not the divisor count). Direct check: snpwt = 2 path's
  # est should equal 2 * (snpwt = NULL path's est) bytewise.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    res_stream = suppressMessages(suppressWarnings(
      admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref, popcombs = .modes_popcombs(),
        auto_only = FALSE, blgsize = 500L, f4mode = TRUE,
        allsnps = FALSE, snpwt = NULL, verbose = FALSE)))
    res_twos = suppressMessages(suppressWarnings(
      admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref, popcombs = .modes_popcombs(),
        auto_only = FALSE, blgsize = 500L, f4mode = TRUE,
        allsnps = FALSE, snpwt = rep(2.0, fix$nsnp), verbose = FALSE)))

    # n (count) is unchanged by scaling (it's the count of finite cells).
    expect_equal(res_stream$n, res_twos$n, tolerance = 0)
    # est is exactly doubled (rowMeans of 2*x = 2*rowMeans(x) when only
    # finite cells contribute).
    expect_equal(2 * res_stream$est, res_twos$est, tolerance = 0)
  })
})

test_that("f4blockdat_from_geno: f4mode = FALSE emits the den (dstatden) column", {
  # With f4mode = FALSE, the function computes both the f4 numerator
  # AND the f3-style denominator (cpp_aftable_to_dstatden). The
  # streaming variant only replaces the numerator path; the denominator
  # is unchanged. Smoke-test that the !f4mode branch completes and
  # returns the documented (est, n, den) column set.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    res_f4 = suppressMessages(suppressWarnings(
      admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref, popcombs = .modes_popcombs(),
        auto_only = FALSE, blgsize = 500L, f4mode = TRUE,
        allsnps = FALSE, verbose = FALSE)))
    res_f3 = suppressMessages(suppressWarnings(
      admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref, popcombs = .modes_popcombs(),
        auto_only = FALSE, blgsize = 500L, f4mode = FALSE,
        allsnps = FALSE, verbose = FALSE)))

    # f4mode=TRUE result has no den column; f4mode=FALSE adds it.
    expect_false("den" %in% names(res_f4))
    expect_true("den" %in% names(res_f3))

    # The streaming path still computes est the same way regardless of
    # f4mode (it's the numerator). est should be identical between the
    # two modes for the same popcomb.
    expect_equal(res_f4$est, res_f3$est, tolerance = 0)
    expect_equal(res_f4$n,   res_f3$n,   tolerance = 0)

    # den should be a real number (or NaN for empty blocks); not NA all
    # the way. Tighter check: at least one block has a finite den value.
    expect_true(any(is.finite(res_f3$den)),
                info = "expected at least one block with a finite dstatden")
  })
})
