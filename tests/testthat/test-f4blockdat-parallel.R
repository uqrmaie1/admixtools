# Tests for the parallel-block path in f4blockdat_from_geno (PR #106).
#
# Equivalence contract: f4blockdat_from_geno must return identical output
# under plan(sequential) (the default) and plan(multisession, workers = N).
# Sequential is the no-op baseline; multisession exercises the actual
# parallel path through furrr::future_map and (for PFILE inputs) the
# per-worker lazy-open of pgenlibr handles in .make_block_reader.

# Common scaffolding: build a fixture, pick blgsize that gives multiple
# blocks (~4 on this fixture's 20 SNPs spanning bp positions 100..2000),
# and a 4-pop popcombs reusing the fixture's 2 pops in redundant slots
# (f4(A, B; A, B) = 0 trivially, but the code path is fully exercised).
.parallel_popcombs = function() {
  tibble::tibble(pop1 = "popA", pop2 = "popB",
                 pop3 = "popA", pop4 = "popB")
}

.run_under_plan = function(pref, plan_setup) {
  on.exit(future::plan(future::sequential), add = TRUE)
  do.call(future::plan, plan_setup)
  suppressMessages(suppressWarnings(
    admixtools:::f4blockdat_from_geno(
      pref = pref,
      popcombs = .parallel_popcombs(),
      auto_only = FALSE,        # fixture is single-chromosome at chr 1
      blgsize = 500L,           # bp-mode (>= 100): ~4 blocks across 20 SNPs
      f4mode = TRUE,
      allsnps = FALSE,
      verbose = FALSE
    )
  ))
}

test_that("f4blockdat_from_geno (BED) returns identical output under sequential and multisession", {
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    res_seq = .run_under_plan(fix$bed_pref, list(future::sequential))
    res_par = .run_under_plan(fix$bed_pref, list(future::multisession, workers = 2L))

    # Same shape, same column types, same values bit-for-bit.
    expect_identical(dim(res_seq),      dim(res_par))
    expect_identical(names(res_seq),    names(res_par))
    expect_identical(res_seq,           res_par)
  })
})

test_that("f4blockdat_from_geno (PFILE) returns identical output under sequential and multisession", {
  # This is the regression test for the PFILE-multisession lazy-open fix in
  # .make_block_reader. Pre-fix, pgenlibr external pointers serialized as
  # dead pointers on workers and f4blockdat_from_geno crashed at the first
  # block under plan(multisession). The fix lazy-opens handles per process
  # (main + each worker) using a Sys.getpid() check.
  testthat::skip_if_not_installed("pgenlibr")
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$pfile_pref), "fixture build skipped (plink2 not on PATH)")

    res_seq = .run_under_plan(fix$pfile_pref, list(future::sequential))
    res_par = .run_under_plan(fix$pfile_pref, list(future::multisession, workers = 2L))

    expect_identical(dim(res_seq),      dim(res_par))
    expect_identical(names(res_seq),    names(res_par))
    expect_identical(res_seq,           res_par)
  })
})

test_that("f4blockdat_from_geno emits no 'UNRELIABLE VALUE' warnings under plan(multisession)", {
  # PR #106's polish sets .options = furrr_options(seed = NULL) on the
  # future_map call, which asserts the mapped function is RNG-free.
  # Without it, furrr emits a "UNRELIABLE VALUE" warning every time it
  # can't statically prove RNG isn't used.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    future::plan(future::multisession, workers = 2L)
    on.exit(future::plan(future::sequential), add = TRUE)

    warns = testthat::capture_warnings(
      suppressMessages(admixtools:::f4blockdat_from_geno(
        pref = fix$bed_pref,
        popcombs = .parallel_popcombs(),
        auto_only = FALSE, blgsize = 500L,
        f4mode = TRUE, allsnps = FALSE,
        verbose = FALSE
      )))
    expect_false(any(grepl("UNRELIABLE VALUE", warns, fixed = TRUE)))
  })
})
