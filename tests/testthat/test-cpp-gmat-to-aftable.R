# Tests for cpp_gmat_to_aftable (PR #113).
#
# Mathematical equivalence to the R one-liner is the headline contract:
#   rowsum(gmat, popvec, na.rm = TRUE) / rowsum((!is.na(gmat))+0, popvec) / 2
# The C++ kernel takes an integer dosage matrix in {0, 1, 2, NA_INTEGER}
# (the BED and PFILE readers produce one natively) and checks g != NA_INTEGER,
# the integer analogue of R's !is.na.

# Reference R implementation, the one-liner that the C++ replaces.
.r_gmat_to_aftable = function(gmat, popvec) {
  rowsum(gmat, popvec, na.rm = TRUE) / rowsum((!is.na(gmat))+0, popvec) / 2
}

# Build a realistic synthetic genotype matrix: integer values in {0, 1, 2, NA}.
# Integer storage matches what the readers now hand the kernel (cpp_read_plink
# -> IntegerMatrix, pgenlibr::ReadIntList -> integer), so the equivalence tests
# exercise the native integer path rather than R's double->integer coercion.
.synthetic_gmat = function(nind, nsnp, na_frac = 0.1, seed = 12345L) {
  withr::with_seed(seed, {
    g = sample(0:2, nind * nsnp, replace = TRUE)
    g[sample.int(length(g), size = round(length(g) * na_frac))] = NA_integer_
    matrix(as.integer(g), nrow = nind, ncol = nsnp)
  })
}

test_that("cpp_gmat_to_aftable matches the R one-liner on realistic genotype data", {
  # 60 pops x 600 inds (10 inds/pop) x 200 SNPs, 10% NA.
  nind = 60L * 10L
  nsnp = 200L
  gmat = .synthetic_gmat(nind, nsnp, na_frac = 0.1)
  popvec = rep(seq_len(60L), each = 10L)

  r_out   = .r_gmat_to_aftable(gmat, popvec)
  cpp_out = admixtools:::cpp_gmat_to_aftable(gmat, as.integer(popvec))

  # rowsum sets rownames on its result ("1", "2", ...); the C++ version
  # returns a plain matrix. Strip rownames before comparing so the
  # comparison is on values only (downstream callers don't read these
  # rownames anyway).
  r_out = unname(r_out)
  expect_equal(dim(cpp_out), dim(r_out))
  expect_equal(cpp_out, r_out, tolerance = 0)
})

test_that("cpp_gmat_to_aftable returns NaN rows for all-NA pops", {
  # 3 pops, but pop 3 has every genotype NA -> result row for pop 3 is
  # all NaN (0/0 from rowsum + na.rm = TRUE).
  gmat = matrix(c(
    0, 1, 2, 0,    # pop 1, ind 1
    1, 0, 1, 2,    # pop 1, ind 2
    0, 2, NA, 1,   # pop 2, ind 3
    2, 1, 2, 0,    # pop 2, ind 4
    NA, NA, NA, NA, # pop 3, ind 5
    NA, NA, NA, NA  # pop 3, ind 6
  ), nrow = 6, byrow = TRUE)
  popvec = c(1L, 1L, 2L, 2L, 3L, 3L)

  r_out   = .r_gmat_to_aftable(gmat, popvec)
  cpp_out = admixtools:::cpp_gmat_to_aftable(gmat, popvec)
  r_out = unname(r_out)
  expect_equal(cpp_out, r_out, tolerance = 0)
  # Pop 3 row is all NaN
  expect_true(all(is.nan(cpp_out[3, ])))
})

test_that("cpp_gmat_to_aftable handles single-pop / single-SNP / single-ind boundaries", {
  # Single pop, single SNP, single individual
  gmat = matrix(1, nrow = 1, ncol = 1)
  popvec = 1L
  cpp_out = admixtools:::cpp_gmat_to_aftable(gmat, popvec)
  expect_equal(dim(cpp_out), c(1L, 1L))
  expect_equal(cpp_out[1, 1], 0.5)  # 1/1/2 = 0.5

  # Single pop, multiple SNPs, multiple inds (no NAs)
  gmat = matrix(c(0, 2, 1, 2,
                  1, 1, 2, 0,
                  2, 0, 1, 1), nrow = 3, byrow = TRUE)
  popvec = c(1L, 1L, 1L)
  cpp_out = admixtools:::cpp_gmat_to_aftable(gmat, popvec)
  r_out = .r_gmat_to_aftable(gmat, popvec)
  r_out = unname(r_out)
  expect_equal(cpp_out, r_out, tolerance = 0)
})

test_that("cpp_gmat_to_aftable treats NA_INTEGER as missing (matches R's !is.na contract)", {
  # gmat arrives as an integer matrix (the BED and PFILE readers produce
  # IntegerMatrix natively; a double matrix from an older caller gets
  # auto-coerced by R, with Inf -> NA per as.integer's contract). The C++
  # side checks `g != NA_INTEGER` rather than ISNAN(g); both are the
  # equivalent R `!is.na()` semantic on their respective types.
  gmat = matrix(c(
    0L, 1L,
    NA_integer_, 2L,
    1L, NA_integer_
  ), nrow = 3, byrow = TRUE)
  popvec = c(1L, 1L, 2L)

  r_out   = .r_gmat_to_aftable(gmat, popvec)
  cpp_out = admixtools:::cpp_gmat_to_aftable(gmat, popvec)
  r_out = unname(r_out)
  # Pop 1 col 1: 0 / 1 / 2 = 0      (NA at row 2 dropped)
  # Pop 1 col 2: (1 + 2) / 2 / 2 = 0.75
  # Pop 2 col 1: 1 / 1 / 2 = 0.5
  # Pop 2 col 2: NA -> 0/0 -> NaN
  expect_equal(cpp_out, r_out, tolerance = 0)
  expect_equal(cpp_out[1, 1], 0)
  expect_equal(cpp_out[1, 2], 0.75)
  expect_true(is.nan(cpp_out[2, 2]))
})

test_that("cpp_gmat_to_aftable errors on popvec containing 0 / negative / NA", {
  gmat = matrix(c(0, 1, 2, 0), nrow = 4)
  expect_error(
    admixtools:::cpp_gmat_to_aftable(gmat, c(1L, 0L, 1L, 1L)),
    "popvec contains values < 1"
  )
  expect_error(
    admixtools:::cpp_gmat_to_aftable(gmat, c(1L, -2L, 1L, 1L)),
    "popvec contains values < 1"
  )
  expect_error(
    admixtools:::cpp_gmat_to_aftable(gmat, c(1L, NA_integer_, 1L, 1L)),
    "popvec contains values < 1"
  )
})
