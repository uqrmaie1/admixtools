# Tests for the parallel-readRDS path in read_f2() (PR #114).
#
# Equivalence contract: the result of read_f2(test_dir) must be identical
# under plan(sequential) (the default) and plan(multisession, workers = N).
# Sequential is the no-op baseline; multisession exercises the actual
# parallel path through furrr::future_map.

# Build a small on-disk f2 cache from example_f2_blocks. read_f2 needs:
#   - one .rds per (sorted) population pair under f2_dir/pop1/pop2_f2.rds
#     (written via write_f2)
#   - block_lengths_f2.rds in f2_dir (written manually, parsed from the
#     dimname-encoded block lengths "l<N>")
.build_test_cache = function(dir) {
  data("example_f2_blocks", package = "admixtools", envir = environment())
  f2 = get("example_f2_blocks", envir = environment())

  # Use only the first 4 populations to keep the cache tiny (6 pairs).
  pops4 = dimnames(f2)[[1]][1:4]
  f2_small = f2[pops4, pops4, , drop = FALSE]
  # Synthesize a counts array; values don't matter for read_f2's parallel
  # path (we only assert pairwise equality, not correctness of counts).
  counts = array(1L, dim = dim(f2_small), dimnames = dimnames(f2_small))

  write_f2(f2_small, counts, outdir = dir, id = 'f2', overwrite = TRUE)

  # block_lengths_f2.rds: read_f2 expects this file in the dir.
  block_lengths = readr::parse_number(dimnames(f2_small)[[3]])
  saveRDS(block_lengths, file.path(dir, 'block_lengths_f2.rds'))

  pops4
}

test_that("read_f2 returns identical output under plan(sequential) and plan(multisession)", {
  withr::with_tempdir({
    cache_dir = file.path(getwd(), "cache")
    dir.create(cache_dir)
    pops = .build_test_cache(cache_dir)

    # Sequential reference (the default plan).
    future::plan(future::sequential)
    res_seq = suppressMessages(suppressWarnings(
      read_f2(cache_dir, pops = pops, verbose = FALSE)))

    # Parallel under plan(multisession, workers = 2). Reset to sequential
    # on exit so we don't leave the future plan dirty for other tests.
    future::plan(future::multisession, workers = 2L)
    on.exit(future::plan(future::sequential), add = TRUE)
    res_par = suppressMessages(suppressWarnings(
      read_f2(cache_dir, pops = pops, verbose = FALSE)))

    # The 3D arrays must be byte-identical across plans. Same dim, same
    # dimnames, same values.
    expect_identical(dim(res_seq),       dim(res_par))
    expect_identical(dimnames(res_seq),  dimnames(res_par))
    expect_identical(res_seq,            res_par)
  })
})

test_that("read_f2 does not emit furrr 'UNRELIABLE VALUE' warnings under plan(multisession)", {
  # PR #114's fix sets .options = furrr_options(seed = NULL) on the
  # future_map call, which asserts the mapped function is RNG-free.
  # Without that, furrr emits a "UNRELIABLE VALUE" warning every time
  # it can't statically prove RNG isn't used.
  withr::with_tempdir({
    cache_dir = file.path(getwd(), "cache")
    dir.create(cache_dir)
    pops = .build_test_cache(cache_dir)

    future::plan(future::multisession, workers = 2L)
    on.exit(future::plan(future::sequential), add = TRUE)

    warns = testthat::capture_warnings(
      suppressMessages(read_f2(cache_dir, pops = pops, verbose = FALSE)))
    expect_false(any(grepl("UNRELIABLE VALUE", warns, fixed = TRUE)))
  })
})

test_that("read_f2 errors clearly when a pair file is missing", {
  # The path-resolution pass must surface the documented error message
  # for missing pair files (preserving the pre-PR behavior). The check
  # is now hoisted out of the inner loop, but the message is unchanged.
  withr::with_tempdir({
    cache_dir = file.path(getwd(), "cache")
    dir.create(cache_dir)
    pops = .build_test_cache(cache_dir)

    # Delete one pair file to trigger the not-found path. write_f2's
    # naming convention sorts pop names; pops[1] = "Altai_Neanderthal.DG",
    # pops[2] = "Chimp.REF" -- the file is .../Altai_Neanderthal.DG/Chimp.REF_f2.rds.
    sorted = sort(pops[1:2])
    target = file.path(cache_dir, sorted[1], paste0(sorted[2], "_f2.rds"))
    expect_true(file.exists(target))
    file.remove(target)

    expect_error(
      suppressMessages(suppressWarnings(read_f2(cache_dir, pops = pops, verbose = FALSE))),
      "not found"
    )
  })
})
