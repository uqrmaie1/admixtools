# omp_thread_budget() must honor OMP_NUM_THREADS (and the higher-precedence
# options(admixtools.omp_num_threads=)) as a hard cap on the kernel team size.
#
# These assert the RESOLVED thread count, not output identity. An earlier
# iteration shipped an "OMP_NUM_THREADS is respected" test that only checked
# output was identical with vs without the env var -- but the kernels are
# bit-identical at any thread count, so that test passed while OMP_NUM_THREADS
# in fact did nothing. Assert the count itself so the lever can't silently break.

test_that("omp_threads_cap reads OMP_NUM_THREADS and the option (option wins)", {
  withr::with_options(list(admixtools.omp_num_threads = NULL), {
    withr::with_envvar(c(OMP_NUM_THREADS = NA), {            # unset
      expect_true(is.na(admixtools:::omp_threads_cap()))
    })
    withr::with_envvar(c(OMP_NUM_THREADS = "3"), {
      expect_identical(admixtools:::omp_threads_cap(), 3L)
    })
    withr::with_envvar(c(OMP_NUM_THREADS = "4,2,1"), {       # nested-level form
      expect_identical(admixtools:::omp_threads_cap(), 4L)
    })
    withr::with_envvar(c(OMP_NUM_THREADS = "garbage"), {     # unparseable -> ignored
      expect_true(is.na(admixtools:::omp_threads_cap()))
    })
  })
  # The R option takes precedence over the env var.
  withr::with_envvar(c(OMP_NUM_THREADS = "8"), {
    withr::with_options(list(admixtools.omp_num_threads = 2L), {
      expect_identical(admixtools:::omp_threads_cap(), 2L)
    })
  })
})

test_that("omp_thread_budget caps the team size at OMP_NUM_THREADS", {
  withr::with_options(list(admixtools.omp_num_threads = NULL), {
    withr::with_envvar(c(OMP_NUM_THREADS = "1"), {
      expect_identical(admixtools:::omp_thread_budget(), 1L)         # forced serial
    })
    withr::with_envvar(c(OMP_NUM_THREADS = "2"), {
      expect_lte(admixtools:::omp_thread_budget(), 2L)               # capped at <= 2
    })
    # A cap above the auto budget never raises it (no oversubscription).
    base = withr::with_envvar(c(OMP_NUM_THREADS = NA),
                              admixtools:::omp_thread_budget())
    high = withr::with_envvar(c(OMP_NUM_THREADS = "1000"),
                              admixtools:::omp_thread_budget())
    expect_identical(high, base)
  })
  # Option path: 1 forces serial even with a permissive env var.
  withr::with_envvar(c(OMP_NUM_THREADS = "8"), {
    withr::with_options(list(admixtools.omp_num_threads = 1L), {
      expect_identical(admixtools:::omp_thread_budget(), 1L)
    })
  })
})
