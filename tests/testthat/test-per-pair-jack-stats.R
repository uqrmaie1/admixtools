# Tests for per_pair_jack_stats (PR #115).
#
# Headline contract: the helper produces a tibble byte-for-byte
# identical to the pre-PR dplyr expression
#
#   out %>% group_by(pop1, pop2) %>%
#     summarize(f2dat = list(f2_blocks[pop1, pop2, ])) %>% ungroup %>%
#     mutate(sts = map(f2dat, ~statfun(., block_lengths)),
#            est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
#     mutate(se = sqrt(var)) %>%
#     select(pop1, pop2, est, se)
#
# on the same (out, f2_blocks, block_lengths, statfun) the PR branch
# builds upstream. Same kernel, same input, only the aggregation
# differs.

# Inline reconstruction of the pre-PR dplyr expression, kept verbatim so
# regressions show up clearly in the diff if it ever drifts.
.dplyr_reference = function(out, f2_blocks, block_lengths, statfun) {
  out %>% dplyr::group_by(pop1, pop2) %>%
    dplyr::summarize(f2dat = list(f2_blocks[pop1, pop2, ]),
                     .groups = "drop") %>%
    dplyr::mutate(sts = purrr::map(f2dat, ~statfun(., block_lengths)),
                  est = purrr::map_dbl(sts, "est"),
                  var = purrr::map_dbl(sts, "var")) %>%
    dplyr::mutate(se = sqrt(var)) %>%
    dplyr::select(pop1, pop2, est, se)
}

test_that("per_pair_jack_stats is bit-identical to the pre-PR dplyr expression (jackknife)", {
  data("example_f2_blocks", package = "admixtools")
  f2_blocks = admixtools:::est_to_loo(example_f2_blocks)
  block_lengths = readr::parse_number(dimnames(f2_blocks)[[3]])
  statfun = admixtools:::cpp_jack_vec_stats

  out = admixtools:::fstat_get_popcombs(
    example_f2_blocks, pop1 = NULL, pop2 = NULL,
    sure = FALSE, unique_only = TRUE, fnum = 2)

  cand = admixtools:::per_pair_jack_stats(out, f2_blocks, block_lengths, statfun)
  ref  = .dplyr_reference(out, f2_blocks, block_lengths, statfun)

  expect_identical(cand, ref)
})

test_that("per_pair_jack_stats is bit-identical to the pre-PR dplyr expression (bootstrap)", {
  # The boot path uses cpp_boot_vec_stats and the boot-resampled blocks
  # produced by est_to_boo. est_to_boo uses sample() internally, so the
  # f2_blocks input must be built once and shared between candidate and
  # reference -- otherwise the two draws differ and the comparison
  # becomes pointless. withr::with_seed pins the resample so the test is
  # reproducible across runs.
  data("example_f2_blocks", package = "admixtools")
  f2_blocks = withr::with_seed(20260519L,
    admixtools:::est_to_boo(example_f2_blocks, nboot = 30L))
  block_lengths = readr::parse_number(dimnames(f2_blocks)[[3]])
  statfun = admixtools:::cpp_boot_vec_stats

  out = admixtools:::fstat_get_popcombs(
    example_f2_blocks, pop1 = NULL, pop2 = NULL,
    sure = FALSE, unique_only = TRUE, fnum = 2)

  cand = admixtools:::per_pair_jack_stats(out, f2_blocks, block_lengths, statfun)
  ref  = .dplyr_reference(out, f2_blocks, block_lengths, statfun)

  expect_identical(cand, ref)
})

test_that("f2() returns the documented shape, types, and sort order", {
  data("example_f2_blocks", package = "admixtools")
  out = f2(example_f2_blocks)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("pop1", "pop2", "est", "se"))
  expect_type(out$pop1, "character")
  expect_type(out$pop2, "character")
  expect_type(out$est, "double")
  expect_type(out$se, "double")

  # Sort order must be locale-independent radix order over (pop1, pop2).
  # If the rows are already in radix order, order() returns 1..n.
  expect_identical(
    order(out$pop1, out$pop2, method = "radix"),
    seq_len(nrow(out)))
})

test_that("fst() returns the documented shape, types, and sort order", {
  data("example_f2_blocks", package = "admixtools")
  out = suppressMessages(suppressWarnings(fst(example_f2_blocks)))

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("pop1", "pop2", "est", "se"))
  expect_type(out$est, "double")
  expect_type(out$se, "double")

  expect_identical(
    order(out$pop1, out$pop2, method = "radix"),
    seq_len(nrow(out)))
})

test_that("f2() output is byte-identical to dplyr group_by() under mixed-case pop labels", {
  # Falsification test for the `method = 'radix'` fix in
  # per_pair_jack_stats. Under any UTF-8 locale, default order()
  # sorts case-insensitively while radix uses byte order:
  #
  #   order(c('anc1', 'ANC2', 'anc3'))                  # C.UTF-8: anc1, ANC2, anc3
  #   order(c('anc1', 'ANC2', 'anc3'), method='radix')  # ANC2, anc1, anc3 (byte: 'A'=0x41 < 'a'=0x61)
  #
  # dplyr::group_by() uses radix sort, so it matches the right-hand
  # side. Without `method = 'radix'` in per_pair_jack_stats, the
  # helper would use the LEFT-hand sort and f2()'s row order would
  # disagree with dplyr's -- expect_identical(cand, ref) would fail
  # on locale != C, even though the per-pair (est, se) values would
  # still match. This test exercises that divergent path.
  #
  # We force LC_COLLATE = en_US.UTF-8 via withr::local_locale so the
  # test is reproducible across systems. If that locale isn't
  # available (e.g. minimal Docker images), skip.
  if(is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8")))) {
    skip("en_US.UTF-8 locale not available")
  }
  # Undo the probe; withr::local_locale will set + restore properly.
  Sys.setlocale("LC_COLLATE", "")
  withr::local_locale(c(LC_COLLATE = "en_US.UTF-8"))

  # Sanity: confirm the locale actually diverges from radix on our
  # test labels. If a future R or libc collation change makes these
  # agree, the test stops being a falsification and we should pick
  # new labels.
  testcase = c("anc1", "ANC2", "anc3")
  stopifnot(!identical(order(testcase), order(testcase, method = "radix")))

  data("example_f2_blocks", package = "admixtools")
  blk = example_f2_blocks
  pops = dimnames(blk)[[1]]
  pops[1:3] = testcase
  dimnames(blk)[[1]] = pops
  dimnames(blk)[[2]] = pops

  # Candidate: f2() output (uses per_pair_jack_stats with radix sort).
  cand = f2(blk)

  # Reference: rebuild the same upstream inputs f2() builds, then run
  # the dplyr group_by chain on them. Same kernel, same input slice,
  # only the aggregation strategy differs. dplyr's group_by uses
  # radix sort, so this sorts byte-order.
  f2_blocks_loo = admixtools:::est_to_loo(blk)
  block_lengths = readr::parse_number(dimnames(f2_blocks_loo)[[3]])
  statfun = admixtools:::cpp_jack_vec_stats
  out = admixtools:::fstat_get_popcombs(
    blk, pop1 = NULL, pop2 = NULL,
    sure = FALSE, unique_only = TRUE, fnum = 2)
  ref = .dplyr_reference(out, f2_blocks_loo, block_lengths, statfun)

  expect_identical(cand, ref)
})
