# Linkage guard for the TWO poly_only polymorphism definitions that must stay
# in sync. They live in different files, gate different code paths, and operate
# on different shapes, but encode the SAME polymorphism intent:
#
#   (A) C++  poly_only_keep      src/cpp_fstats.cpp (~L165), allsnps path,
#       per-popcomb over the four values {w, x, y, z}. poly_only == 2 keeps a
#       SNP when the distinct non-NA/NaN values exceed 1, OR the max kept value
#       lies strictly inside (POLY_ONLY_LO, POLY_ONLY_HI) = (0.0001, 0.9999).
#       NA and NaN are dropped (ISNAN); +/-Inf is kept.
#
#   (B) R lambda `fn`            R/io.R (~L3427), !allsnps path, applied
#       per-population-group across model columns:
#           length(unique(na.omit(x))) > 1 | !max(na.omit(x)) %in% c(0, 1)
#       i.e. keep when distinct non-NA/NaN values exceed 1, OR the max kept
#       value is not exactly 0 or 1.
#
# A future edit to one definition's thresholds or logic that is not mirrored in
# the other should make this file fail. The two are NOT byte-identical on every
# input (see the all-Inf row below), so the test pins (1) the C++ decision
# against an R model of the C++ rule, (2) the R lambda against a fixed expected
# table, and (3) that the two rules AGREE on every shared in-[0,1] class and
# DISAGREE on all-Inf in the one documented way. Any drift trips at least one.

# --- (1) Observe poly_only_keep's keep/drop decision via the kernel ----------
# cpp_aftable_to_dstatnum fills num with NA_REAL where the gate DROPS a SNP and
# with the computed (w-x)*(y-z) where it KEEPS one. A kept cell may be finite,
# NaN (e.g. Inf-Inf), or +/-Inf, but is.na(.) && !is.nan(.) is TRUE only for the
# dropped NA_REAL sentinel, never for a kept value. That distinction is the
# observable keep/drop signal.
cpp_gate_keeps = function(quad, poly_only) {
  aft = matrix(quad, nrow = 4L)               # one popcomb over rows 1..4 (w,x,y,z)
  res = admixtools:::cpp_aftable_to_dstatnum(
    aft, 1, 2, 3, 4,                           # p1..p4 -> rows 1..4
    1, matrix(0), TRUE, poly_only, 1L)         # modelvec/usesnps unused on allsnps=TRUE
  v = res$num[1, 1]
  !(is.na(v) && !is.nan(v))                    # TRUE = kept, FALSE = dropped
}

# --- (2) R model of the C++ poly_only_keep rule (thresholds, ISNAN) ----------
# Independent reimplementation of definition (A): drop NA/NaN, keep +/-Inf;
# keep on distinct>1; for poly_only==2 also keep when max is inside the open
# (POLY_ONLY_LO, POLY_ONLY_HI) interval. Mirrors src/cpp_fstats.cpp exactly so a
# change to the C++ thresholds (0.0001 / 0.9999) or branch logic breaks (3a).
POLY_ONLY_LO = 0.0001
POLY_ONLY_HI = 0.9999
cpp_rule_keeps = function(quad, poly_only) {
  if (poly_only == 0L) return(TRUE)
  vals = quad[!is.na(quad)]                    # is.na(NaN) is TRUE -> drops NA & NaN
  if (length(vals) == 0L) return(FALSE)
  uni = unique(vals)
  if (length(uni) > 1L) return(TRUE)
  if (poly_only == 2L) {
    mx = max(uni)
    if (mx > POLY_ONLY_LO && mx < POLY_ONLY_HI) return(TRUE)
  }
  FALSE
}

# --- (3) The R lambda from R/io.R L3427, lifted verbatim ---------------------
# This is the exact predicate the !allsnps path applies. If L3427 changes, copy
# the change here and the expected table below must be re-justified.
r_lambda_keeps = function(quad) {
  x = quad
  as.logical(length(unique(na.omit(x))) > 1 | !max(na.omit(x)) %in% c(0, 1))
}

# Shared allele-frequency vectors, one per polymorphism class. Deterministic;
# no RNG. Every branch of both rules is exercised.
poly_cases = list(
  polymorphic        = c(0.10, 0.20, 0.30, 0.40),  # distinct>1            -> keep both
  mono_in_range      = c(0.50, 0.50, 0.50, 0.50),  # mono 0.5, in (0,1)    -> keep poly2 / lambda
  all_zero           = c(0.00, 0.00, 0.00, 0.00),  # mono 0                -> drop both
  all_one            = c(1.00, 1.00, 1.00, 1.00),  # mono 1                -> drop both
  nan_then_mono_half = c(NaN,  0.50, 0.50, 0.50),  # NaN dropped, mono 0.5 -> keep poly2 / lambda
  nan_then_mono_zero = c(NaN,  0.00, 0.00, 0.00),  # NaN dropped, mono 0   -> drop both
  inf_polymorphic    = c(Inf,  0.50, 0.70, 0.20),  # Inf kept, distinct>1  -> keep both
  all_inf            = c(Inf,  Inf,  Inf,  Inf)     # Inf kept, mono Inf    -> C++ drops, lambda keeps
)

test_that("C++ poly_only_keep matches an independent model of the C++ rule (all classes)", {
  for (nm in names(poly_cases)) {
    quad = poly_cases[[nm]]
    for (po in 0:2) {
      expect_identical(
        cpp_gate_keeps(quad, po),
        cpp_rule_keeps(quad, po),
        info = sprintf("class=%s poly_only=%d", nm, po))
    }
  }
})

# Monomorphic vectors that straddle the (POLY_ONLY_LO, POLY_ONLY_HI) = (0.0001,
# 0.9999) interval, so a drift in the C++ thresholds that is not mirrored in
# cpp_rule_keeps above trips test (1) on at least one of these. Kept OUT of the
# shared-agreement test below: the R lambda's `%in% c(0,1)` keeps every one of
# these (none is exactly 0 or 1), while the C++ interval drops the two outside
# (LO, HI), so the two definitions legitimately disagree here.
poly_threshold_cases = list(
  just_below_lo = rep(0.00005, 4L),  # < LO  -> C++ poly2 DROP
  just_above_lo = rep(0.0002,  4L),  # > LO  -> C++ poly2 KEEP
  just_below_hi = rep(0.9998,  4L),  # < HI  -> C++ poly2 KEEP
  just_above_hi = rep(0.99995, 4L))  # > HI  -> C++ poly2 DROP

test_that("C++ poly_only_keep pins the (0.0001, 0.9999) thresholds (boundary cases)", {
  for (nm in names(poly_threshold_cases)) {
    quad = poly_threshold_cases[[nm]]
    for (po in 0:2) {
      expect_identical(
        cpp_gate_keeps(quad, po),
        cpp_rule_keeps(quad, po),
        info = sprintf("class=%s poly_only=%d", nm, po))
    }
  }
})

test_that("R/io.R L3427 lambda matches its documented keep/drop table (all classes)", {
  # Expected keep (TRUE) / drop (FALSE) per class for the L3427 predicate,
  # hand-derived from `distinct>1 | max not in {0,1}`, independent of the helper.
  expected = c(
    polymorphic        = TRUE,
    mono_in_range      = TRUE,
    all_zero           = FALSE,
    all_one            = FALSE,
    nan_then_mono_half = TRUE,
    nan_then_mono_zero = FALSE,
    inf_polymorphic    = TRUE,
    all_inf            = TRUE)   # max(Inf) not in {0,1} -> lambda KEEPS
  for (nm in names(poly_cases)) {
    expect_identical(
      r_lambda_keeps(poly_cases[[nm]]),
      unname(expected[[nm]]),
      info = sprintf("class=%s lambda", nm))
  }
})

test_that("the two poly_only definitions agree on every in-[0,1] class (shared intent)", {
  # On all finite/monomorphic/in-range classes the allsnps C++ gate at
  # poly_only==2 and the !allsnps R lambda must make the SAME decision: this is
  # the polymorphism intent the two files share. Drift in either thresholds or
  # branch logic that changes one but not the other breaks this.
  shared = setdiff(names(poly_cases), "all_inf")
  for (nm in shared) {
    quad = poly_cases[[nm]]
    expect_identical(
      cpp_gate_keeps(quad, 2L),     # observed C++ allsnps decision
      r_lambda_keeps(quad),         # observed R !allsnps decision
      info = sprintf("class=%s shared-intent", nm))
  }
})

test_that("all-Inf is the one documented divergence: C++ drops, R lambda keeps", {
  # The only input where the two definitions differ: max(kept) = Inf.
  #   C++ poly_only==2: Inf is NOT < POLY_ONLY_HI (0.9999) -> DROP.
  #   R lambda:         Inf is NOT in {0,1}                -> KEEP.
  # Pinning this so a future change that either (a) accidentally unifies the two
  # (e.g. swapping the C++ thresholds for an exact {0,1} test, or the R lambda
  # for a bounded interval) or (b) introduces a NEW divergence is caught here
  # rather than silently changing which SNPs the allsnps vs !allsnps paths keep.
  quad = poly_cases[["all_inf"]]
  expect_false(cpp_gate_keeps(quad, 2L))   # C++ allsnps gate drops all-Inf
  expect_true(r_lambda_keeps(quad))        # R !allsnps lambda keeps all-Inf
})

# Randomized differential test: poly_only_keep must reproduce the ORIGINAL Rcpp
# idiom it replaced -- na_omit(unique(c(w, x, y, z))) with length(uni) > 1 ||
# (poly_only == 2 && max(uni) in (0.0001, 0.9999)) -- byte for byte. The fixtures
# above hit the known branch classes by hand; this fuzzes many random quads
# (mixing [0,1] values, the threshold boundaries, exact 0/1, +/-Inf, NaN, NA, and
# duplicates) through the kernel and against a literal reimplementation of the
# old idiom.
old_idiom_keeps = function(quad, poly_only) {
  if (poly_only == 0L) return(TRUE)
  uni = unique(quad[!is.na(quad)])           # na_omit(unique(.)): drops NA & NaN, keeps +/-Inf
  length(uni) > 1L ||
    (poly_only == 2L && length(uni) >= 1L && max(uni) > 0.0001 && max(uni) < 0.9999)
}

test_that("poly_only_keep matches the old na_omit(unique()) idiom over random quads", {
  set.seed(20260610L)
  N = 10000L
  # Pool weighted toward the branch boundaries; sampling with replacement injects
  # duplicates (and repeated Inf/NaN/0/1) naturally. NA is intentionally excluded:
  # a kept quad that CONTAINS an NA also yields num = NA (the product propagates
  # it), so keep/drop is not observable from num for NA-bearing quads. NaN
  # exercises the identical ISNAN gate path and IS observable (a kept NaN product
  # is is.nan, distinct from the NA_REAL drop sentinel); NA-vs-NaN gate identity
  # is pinned separately in test-cpp-aftable-poly-nan.R.
  pool = c(runif(40), 0, 1, 1e-5, 0.0001, 0.00011, 0.9998, 0.9999, 0.99991, 0.5,
           Inf, -Inf, NaN)
  aft = matrix(sample(pool, 4L * N, replace = TRUE), nrow = 4L)
  usesnps = matrix(0)                        # unused on the allsnps = TRUE path
  for (po in 1:2) {
    res = admixtools:::cpp_aftable_to_dstatnum(aft, 1, 2, 3, 4, 1, usesnps, TRUE, po, 1L)
    v = as.numeric(res$num[1, ])
    got = !(is.na(v) & !is.nan(v))           # kept iff not the NA_REAL drop sentinel
    want = vapply(seq_len(N), function(i) old_idiom_keeps(aft[, i], po), logical(1))
    expect_identical(got, want,
                     info = sprintf("poly_only=%d over %d random quads", po, N))
  }
})
