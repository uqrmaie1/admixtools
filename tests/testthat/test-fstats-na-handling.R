test_that("xmats_to_pairarrs zeroes missing genotypes element-wise (replace_na no-ops on a matrix)", {

  # SNP x population dosage matrix. fix_ploidy reads every column as diploid
  # because each contains 0, 1 and 2, so no pseudohaploid rescaling occurs.
  base <- matrix(c(0, 1, 2,
                   2, 0, 1,
                   1, 2, 0,
                   2, 1, 0), nrow = 4, byrow = TRUE,
                 dimnames = list(NULL, c("A", "B", "C")))

  # Two extra SNPs missing in B only (NA in column B, present in A and C).
  # These are exactly the partial-NA rows that survive an extract with
  # maxmiss > 0 and that replace_na() fails to zero on a matrix.
  add <- matrix(c(1, NA, 2,
                  0, NA, 1), nrow = 2, byrow = TRUE,
                dimnames = list(NULL, c("A", "B", "C")))
  ext <- rbind(base, add)

  ratio_pair <- function(m, p, q) {
    arrs <- admixtools:::xmats_to_pairarrs(m, m)
    bl <- nrow(m)
    (admixtools:::block_arr_mean(arrs$aa, bl) /
       admixtools:::block_arr_mean(arrs$nn, bl))[p, q, 1]
  }

  # 1) Regression: no NA leaks into the allele-product / count arrays.
  arrs_ext <- admixtools:::xmats_to_pairarrs(ext, ext)
  expect_true(all(is.finite(arrs_ext$aa)))
  expect_true(all(is.finite(arrs_ext$nn)))

  # 2) Correctness: a SNP missing in B contributes 0 to both aa and nn for
  #    B-pairs, so adding such SNPs must not change the B-pair allele-product
  #    ratio. Under the bug the dropped-from-aa-but-not-nn mismatch inflated it.
  expect_equal(ratio_pair(ext, "A", "B"), ratio_pair(base, "A", "B"))
  expect_equal(ratio_pair(ext, "B", "C"), ratio_pair(base, "B", "C"))

  # Sanity: the added SNPs are present for the A-C pair, so that ratio DOES
  # change. This confirms the invariance above is specific to the missing pop
  # and the test is exercising real arithmetic.
  expect_false(isTRUE(all.equal(ratio_pair(ext, "A", "C"),
                                ratio_pair(base, "A", "C"))))
})
