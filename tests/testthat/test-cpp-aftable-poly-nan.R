# Equivalence guard for the poly_only_keep refactor in src/cpp_fstats.cpp.
#
# poly_only_keep replaced the Rcpp idiom
#     na_omit(unique(NumericVector::create(w, x, y, z)))
# with a pure C++ gate so the popcomb loop could be parallelised. The
# replacement turns on the exact NA / NaN / Inf semantics of na_omit, but the
# shared kernel fixture (helper-dstat-fixture.R) only seeds NA, never a bare
# NaN or an infinity, so the serial vs parallel tests never exercise the input
# class the helper was actually written for.
#
# These tests pin the kernels against an independent R oracle (ref_keep below)
# derived from the documented pre refactor behaviour, on a hand built aftable
# that contains NaN, +Inf, all 0 and all 1 columns. na_omit drops NA and NaN
# (is.na(NaN) is TRUE) and keeps +/-Inf, so the gate keeps a SNP whose four
# values include an infinity but drops one that is monomorphic once NaN is
# removed. allsnps = TRUE throughout, which is the path poly_only_keep gates.

# na_omit(unique(c(w, x, y, z))) polymorphism gate, written in R and
# independent of the C++ helper it is meant to pin.
ref_keep = function(quad, poly_only) {
  if (poly_only == 0L) return(TRUE)
  vals = quad[!is.na(quad)]                 # drop NA and NaN, keep +/-Inf
  uni = unique(vals)
  if (length(uni) > 1L) return(TRUE)        # length(uni) > 1 branch
  if (poly_only == 2L && length(uni) >= 1L &&
      max(uni) > 0.0001 && max(uni) < 0.9999) return(TRUE)
  FALSE
}

# Per cell numerator for the single popcomb (rows 1..4 are w, x, y, z),
# matching the kernel: NA_real_ where the SNP is dropped, otherwise
# (w - x) * (y - z) under IEEE rules. NA is promoted to NaN first because C++
# arithmetic on a kept value never reaches this oracle (the fixture has no kept
# pure NA column), and promoting keeps the helper total.
ref_num_row = function(aft, poly_only) {
  vapply(seq_len(ncol(aft)), function(i) {
    quad = aft[, i]
    if (!ref_keep(quad, poly_only)) return(NA_real_)
    ieee = quad; ieee[is.na(ieee)] = NaN
    (ieee[1] - ieee[2]) * (ieee[3] - ieee[4])
  }, numeric(1))
}

# Classify each cell as NA / NaN / Inf / -Inf / finite, so the structure is
# asserted exactly without depending on which NaN payload the FPU produces.
cell_class = function(v) {
  out = rep("finite", length(v))
  out[is.na(v) & !is.nan(v)] = "NA"
  out[is.nan(v)] = "NaN"
  inf = is.infinite(v)
  out[inf] = ifelse(v[inf] > 0, "Inf", "-Inf")
  out
}

# One popcomb over four populations; each column is a (w, x, y, z) quad chosen
# to land on a distinct branch of the gate.
poly_nan_aftable = function() {
  cols = list(
    c(0.10, 0.20, 0.30, 0.40),    # distinct finite       keep all;        v = 0.01
    c(0.50, 0.50, 0.50, 0.50),    # monomorphic 0.5        keep poly2 only; v = 0
    c(0.00, 0.00, 0.00, 0.00),    # all 0                  drop poly1 & poly2
    c(1.00, 1.00, 1.00, 1.00),    # all 1                  drop poly1 & poly2
    c(NaN,  0.50, 0.50, 0.50),    # NaN dropped, rest 0.5  keep poly2;      v = NaN
    c(NaN,  0.00, 0.00, 0.00),    # NaN dropped, rest 0    drop poly1 & poly2
    c(Inf,  0.50, 0.70, 0.20),    # Inf kept, polymorphic  keep all;        v = Inf
    c(Inf,  Inf,  0.50, 0.50),    # Inf kept, uni{Inf,0.5} keep all;        v = NaN
    c(Inf,  Inf,  Inf,  Inf)      # all Inf, uni{Inf}      drop poly1 & poly2
  )
  matrix(unlist(cols), nrow = 4L)
}

p1 = 1; p2 = 2; p3 = 3; p4 = 4              # popcomb over rows 1..4
modelvec = 1; usesnps = matrix(0)           # unused on the allsnps = TRUE path

test_that("cpp_aftable_to_dstatnum: poly_only gate matches the na_omit oracle on NaN / Inf", {
  aft = poly_nan_aftable()
  for (po in 0:2) {
    res = admixtools:::cpp_aftable_to_dstatnum(
      aft, p1, p2, p3, p4, modelvec, usesnps, TRUE, po, 1L)
    got = as.numeric(res$num[1, ])
    ref = ref_num_row(aft, po)
    # NA vs NaN vs Inf structure, exact.
    expect_identical(cell_class(got), cell_class(ref),
                     info = sprintf("poly_only=%d cell structure", po))
    # finite values agree.
    fin = is.finite(ref)
    expect_equal(got[fin], ref[fin],
                 info = sprintf("poly_only=%d finite values", po))
    # cnt counts exactly the finite cells.
    expect_equal(as.numeric(res$cnt), sum(fin),
                 info = sprintf("poly_only=%d cnt", po))
  }
})

test_that("cpp_aftable_to_dstatnum: literal poly_only=2 values cross-check the oracle", {
  # Hand verified expected values for poly_only = 2 on poly_nan_aftable(),
  # independent of ref_keep / ref_num_row, so an oracle that shares a
  # misconception with the kernel cannot pass silently.
  #   col 1 distinct {.1,.2,.3,.4}   keep   v = (.1-.2)*(.3-.4) = 0.01
  #   col 2 mono 0.5                 keep   v = 0
  #   col 3 all 0                    drop   (max 0 not > 0.0001)
  #   col 4 all 1                    drop   (max 1 not < 0.9999)
  #   col 5 NaN dropped, rest 0.5    keep   v = (NaN-.5)*0 = NaN
  #   col 6 NaN dropped, rest 0      drop
  #   col 7 Inf, polymorphic         keep   v = (Inf-.5)*(.7-.2) = Inf
  #   col 8 Inf, uni{Inf, 0.5}       keep   v = (Inf-Inf)*0 = NaN
  #   col 9 all Inf, uni{Inf}        drop   (max Inf not < 0.9999)
  aft = poly_nan_aftable()
  res = admixtools:::cpp_aftable_to_dstatnum(
    aft, p1, p2, p3, p4, modelvec, usesnps, TRUE, 2L, 1L)
  got = as.numeric(res$num[1, ])
  expect_identical(
    cell_class(got),
    c("finite", "finite", "NA", "NA", "NaN", "NA", "Inf", "NaN", "NA"))
  expect_equal(got[c(1, 2)], c(0.01, 0))
  expect_equal(as.numeric(res$cnt), 2)
})

test_that("cpp_aftable_to_dstatnum_rowmeans: poly_only gate matches the oracle on NaN / Inf", {
  aft = poly_nan_aftable()
  for (po in 0:2) {
    res = admixtools:::cpp_aftable_to_dstatnum_rowmeans(
      aft, p1, p2, p3, p4, modelvec, usesnps, TRUE, po, 1L)
    ref = ref_num_row(aft, po)
    fin = is.finite(ref)
    expect_equal(as.numeric(res$cnt), sum(fin),
                 info = sprintf("poly_only=%d cnt", po))
    expect_equal(as.numeric(res$means)[1],
                 if (any(fin)) sum(ref[fin]) / sum(fin) else NaN,
                 info = sprintf("poly_only=%d mean", po))
  }
})

test_that("poly_only gate treats NA and NaN identically (ISNAN semantics)", {
  # Column 1 is monomorphic once the top value is dropped, so it is dropped by
  # the gate under poly_only 1 and 2 whether the top value is NA or NaN, and no
  # arithmetic ever touches it. Column 2 is distinct and kept. The two fixtures
  # must therefore produce byte identical kernel output, which is the
  # observable consequence of ISNAN treating NA and NaN the same.
  mk = function(top) matrix(c(top, 0, 0, 0,  0.1, 0.2, 0.3, 0.4), nrow = 4L)
  aft_na  = mk(NA_real_)
  aft_nan = mk(NaN)
  for (po in 1:2) {
    r_na  = admixtools:::cpp_aftable_to_dstatnum(
      aft_na,  p1, p2, p3, p4, modelvec, usesnps, TRUE, po, 1L)
    r_nan = admixtools:::cpp_aftable_to_dstatnum(
      aft_nan, p1, p2, p3, p4, modelvec, usesnps, TRUE, po, 1L)
    expect_identical(r_na$num, r_nan$num, info = sprintf("poly_only=%d num", po))
    expect_identical(r_na$cnt, r_nan$cnt, info = sprintf("poly_only=%d cnt", po))
  }
})
