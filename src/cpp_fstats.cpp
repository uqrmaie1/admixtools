// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// Defensive uniformity with cpp_qpadm.cpp / cpp_qpgraph.cpp: this TU contains
// no arma::solve/inv/pinv/chol/svd calls today, but the macro is set here so
// that any future helper added to this file (e.g. a covariance-matrix
// regularization path) inherits the package-wide warning policy without an
// easy-to-miss audit. ARMA_WARN_LEVEL 1 keeps errors, drops level-2 stderr
// warnings. See src/cpp_qpadm.cpp for the full rationale.
#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// Resolve the OpenMP team size for the parallel kernels in this file
// by delegating to parallelly::availableCores(). This composes with the
// package's future::plan() convention:
//
//   * plan(sequential): full core count, OMP uses all cores.
//   * plan(multisession, workers = W): inside each future worker,
//     parallelly returns 1 (future sets the relevant fallback), so
//     each worker runs serial OMP, blocks dispatch via furrr, total
//     active threads = W, no oversubscription.
//   * Honors OMP_NUM_THREADS, MC_CORES, _R_CHECK_LIMIT_CORES_, and the
//     parallelly.availableCores.* options as explicit overrides.
//
// Call this before entering a #pragma omp parallel region (the R API
// is not thread-safe). The integer result is cached process-locally;
// parallelly::availableCores() is ~300 us on aarch64 macOS and ~1 ms
// on x86_64 Linux (probes /proc and cgroups), and a dense qpdstat
// run invokes the kernels ~1400 times, so caching matters. The cache
// is fixed for the life of the current R session (or future worker);
// users who flip OMP_NUM_THREADS mid-session should restart R.
static int admixtools_omp_threads() {
  static int cached_n = -1;
  if (cached_n < 0) {
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("parallelly");
    Rcpp::Function availableCores = pkg["availableCores"];
    int n = Rcpp::as<int>(availableCores());
    cached_n = (n > 0) ? n : 1;
  }
  return cached_n;
}

// Serial pre-pass shared by the three cpp_aftable_to_dstat* functions.
// Catches NA / NaN / Inf and out-of-range integer values in p1..p4
// and modelvec before the parallel region. Errors surface as a clean
// Rcpp::stop with the offending row index, instead of a wild memory
// access or silent corruption inside the OMP region (Rcpp::stop is
// not callable from parallel code). O(npopcomb), negligible vs the
// inner i loop.
static inline void validate_dstat_inputs(
    const arma::mat& aftable,
    const arma::vec& p1, const arma::vec& p2,
    const arma::vec& p3, const arma::vec& p4,
    const arma::vec& modelvec, const arma::mat& usesnps,
    bool allsnps, int npopcomb) {
  const int npop   = (int)aftable.n_rows;
  const int nmodel = (int)usesnps.n_rows;
  for (int j = 0; j < npopcomb; ++j) {
    if (!R_FINITE(p1(j)) || !R_FINITE(p2(j)) ||
        !R_FINITE(p3(j)) || !R_FINITE(p4(j))) {
      Rcpp::stop("p1/p2/p3/p4 contains NA / NaN / Inf at index %d", j + 1);
    }
    const int i1 = (int)p1(j) - 1, i2 = (int)p2(j) - 1,
              i3 = (int)p3(j) - 1, i4 = (int)p4(j) - 1;
    if (i1 < 0 || i1 >= npop || i2 < 0 || i2 >= npop ||
        i3 < 0 || i3 >= npop || i4 < 0 || i4 >= npop) {
      Rcpp::stop("p1/p2/p3/p4 index out of range [1..%d] at j=%d", npop, j + 1);
    }
    if (!allsnps) {
      if (!R_FINITE(modelvec(j))) {
        Rcpp::stop("modelvec contains NA / NaN / Inf at index %d", j + 1);
      }
      const int m = (int)modelvec(j) - 1;
      if (m < 0 || m >= nmodel) {
        Rcpp::stop("modelvec value out of range [1..%d] at j=%d", nmodel, j + 1);
      }
    }
  }
}

// Pure C++ replacement for the Rcpp idiom
//   uni = na_omit(unique(NumericVector::create(w, x, y, z)));
//   valid = uni.length() > 1 || (poly_only == 2 && max(uni) in (0, 1));
// Thread-safe (no Rcpp R-API allocations) so the j loop can be
// parallelized via OpenMP. Uses ISNAN to match na_omit / R's !is.na
// (dropping NA and NaN but keeping +/-Inf), preserving byte-for-byte
// equivalence with the upstream idiom on all inputs.
static inline bool poly_only_keep(double w, double x, double y, double z, int poly_only) {
  if (poly_only == 0) return true;
  double kept[4];
  int nk = 0;
  if (!ISNAN(w)) kept[nk++] = w;
  if (!ISNAN(x)) kept[nk++] = x;
  if (!ISNAN(y)) kept[nk++] = y;
  if (!ISNAN(z)) kept[nk++] = z;
  if (nk == 0) return false;
  // length(uni) > 1 -- any two kept values differ
  bool has_diff = false;
  for (int k = 1; k < nk; ++k) {
    if (kept[k] != kept[0]) { has_diff = true; break; }
  }
  if (has_diff) return true;
  if (poly_only == 2) {
    double maxv = kept[0];
    for (int k = 1; k < nk; ++k) if (kept[k] > maxv) maxv = kept[k];
    if (maxv > 0.0001 && maxv < 0.9999) return true;
  }
  return false;
}

// [[Rcpp::export]]
List cpp_aftable_to_dstatnum_old(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4) {

  // aftable is npop x nsnp
  // num is npopcomb x nsnp
  mat num(p1.n_elem, aftable.n_cols);
  vec cnt = zeros(p1.n_elem);
  double w, x, y, z;
  int i1, i2, i3, i4;
  for(int j = 0; j < p1.n_elem; j++) {
    i1 = p1(j)-1;
    i2 = p2(j)-1;
    i3 = p3(j)-1;
    i4 = p4(j)-1;
    for(int i = 0; i < aftable.n_cols; i++) {
      w = aftable(i1, i);
      x = aftable(i2, i);
      y = aftable(i3, i);
      z = aftable(i4, i);
      num(j, i) = (w - x) * (y - z);
      if(std::isfinite(num(j, i))) cnt(j) += 1;
    }
  }
  //return num;
  return Rcpp::List::create(_["num"] = num, _["cnt"] = cnt);
}

// [[Rcpp::export]]
List cpp_aftable_to_dstatnum(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4,
                             arma::vec& modelvec, arma::mat& usesnps, bool allsnps, int poly_only) {

  // aftable is npop x nsnp
  // num is npopcomb x nsnp
  mat num(p1.n_elem, aftable.n_cols);
  num.fill(NA_REAL);
  vec cnt = zeros(p1.n_elem);

  const int npopcomb = (int)p1.n_elem;
  const int nsnp     = (int)aftable.n_cols;
  validate_dstat_inputs(aftable, p1, p2, p3, p4, modelvec, usesnps, allsnps, npopcomb);

  // Parallel over the popcomb dimension. Each thread writes a distinct
  // row of num and a distinct slot of cnt, so no reduction is needed.
  // poly_only_keep is pure-C++ (no Rcpp allocations) so it is safe
  // inside the parallel region.
  const int nthreads = admixtools_omp_threads();
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for (int j = 0; j < npopcomb; ++j) {
    const int i1 = (int)p1(j) - 1;
    const int i2 = (int)p2(j) - 1;
    const int i3 = (int)p3(j) - 1;
    const int i4 = (int)p4(j) - 1;
    const int m  = allsnps ? 0 : (int)modelvec(j) - 1;
    double row_cnt = 0.0;
    for (int i = 0; i < nsnp; ++i) {
      if (!allsnps && !usesnps(m, i)) continue;
      const double w = aftable(i1, i);
      const double x = aftable(i2, i);
      const double y = aftable(i3, i);
      const double z = aftable(i4, i);
      bool valid = true;
      if (poly_only && allsnps) valid = poly_only_keep(w, x, y, z, poly_only);
      if (valid) {
        const double v = (w - x) * (y - z);
        num(j, i) = v;
        if (std::isfinite(v)) row_cnt += 1.0;
      }
    }
    cnt(j) = row_cnt;
  }
  return Rcpp::List::create(_["num"] = num, _["cnt"] = cnt);
}


// Streaming variant of cpp_aftable_to_dstatnum.
//
// Mathematically equivalent to:
//   tmp = cpp_aftable_to_dstatnum(...)
//   list(means = rowMeans(tmp$num, na.rm = TRUE), cnt = tmp$cnt)
//
// but accumulates per-row sum and count in scalars instead of materializing
// the full (p1.n_elem x aftable.n_cols) matrix. This serves callers that
// only need the row means (currently the f4mode path in
// f4blockdat_from_geno when no per-SNP weights are applied).
//
// Motivation: on dense f-statistic runs (npopcomb in the millions, blocks
// of a few thousand SNPs) the materialized matrix in cpp_aftable_to_dstatnum
// is a many-GB transient that exceeds production servers' RAM. At
// 1.57M popcombs x 3608 SNPs/block, the matrix is 1.57e6 * 3608 * 8 bytes
// = ~45 GB per call. R either errors with "cannot allocate vector of size
// X" or the OS OOM-kills the R process during the fill. Streaming brings
// per-block peak from O(p1.n_elem * aftable.n_cols) down to O(p1.n_elem)
// -- ~12 MB at the same scale -- so f4blockdat_from_geno completes on
// machines where the materialized variant doesn't fit.
//
// (Earlier framings of this PR mentioned 32-bit integer overflow in arma's
// indexing. On 64-bit builds RcppArmadillo auto-enables ARMA_64BIT_WORD,
// so arma::uword is u64 and the matrix indexing has no 32-bit ceiling at
// realistic admixtools scales. The bottleneck is memory, not arithmetic.)
//
// The poly_only / usesnps / allsnps semantics, NA propagation, and `cnt`
// definition match cpp_aftable_to_dstatnum exactly; the only change is
// that `means` (= sum / cnt, which yields NaN when cnt == 0, matching
// rowMeans on an all-NA row) is returned in place of `num`.
//
// Precision: the per-row sum is accumulated in `long double`, matching
// R's rowMeans which uses LDOUBLE internally (see src/main/array.c's
// do_colsum / do_rowmeans). On x86_64 that's 80-bit x87 extended
// precision; on arm64 / Apple Silicon long double == double (64-bit).
// Using long double preserves byte-identical output vs the
// pre-streaming-PR materialized + rowMeans path on every platform.
// (Switching from double accumulation to long double costs ~5-10% in
// the hot sum loop on x86_64 because x87 long-double ops aren't SSE-
// vectorizable, but this is well below the win from skipping the
// O(npopcomb * nsnp) matrix allocation + Rcpp::wrap + R-side rowMeans
// pass that the streaming path replaces.)
//
// [[Rcpp::export]]
List cpp_aftable_to_dstatnum_rowmeans(arma::mat& aftable,
                                      arma::vec& p1, arma::vec& p2,
                                      arma::vec& p3, arma::vec& p4,
                                      arma::vec& modelvec,
                                      arma::mat& usesnps,
                                      bool allsnps, int poly_only) {

  // sum_v in long double matches R's rowMeans precision; cnt is an
  // integer count where 64-bit double is more than enough.
  const int npopcomb = (int)p1.n_elem;
  const int nsnp     = (int)aftable.n_cols;
  validate_dstat_inputs(aftable, p1, p2, p3, p4, modelvec, usesnps, allsnps, npopcomb);

  std::vector<long double> sum_v(npopcomb, 0.0L);
  vec cnt = zeros<vec>(npopcomb);

  // Parallel over the popcomb dimension. Each thread writes a distinct
  // sum_v[j] and cnt(j); aftable / usesnps are read-only shared. The
  // poly_only path uses poly_only_keep (pure C++, no Rcpp), so the
  // parallel region is safe regardless of poly_only && allsnps.
  const int nthreads = admixtools_omp_threads();
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for (int j = 0; j < npopcomb; ++j) {
    const int i1 = (int)p1(j) - 1;
    const int i2 = (int)p2(j) - 1;
    const int i3 = (int)p3(j) - 1;
    const int i4 = (int)p4(j) - 1;
    const int m  = allsnps ? 0 : (int)modelvec(j) - 1;
    long double s = 0.0L;
    double      c = 0.0;
    for (int i = 0; i < nsnp; ++i) {
      if (!allsnps && !usesnps(m, i)) continue;
      const double w = aftable(i1, i);
      const double x = aftable(i2, i);
      const double y = aftable(i3, i);
      const double z = aftable(i4, i);
      bool valid = true;
      if (poly_only && allsnps) valid = poly_only_keep(w, x, y, z, poly_only);
      if (!valid) continue;
      const double val = (w - x) * (y - z);
      if (std::isfinite(val)) {
        // double -> long double promotion in the add, matching R's
        // LDOUBLE accumulator behavior.
        s += val;
        c += 1.0;
      }
    }
    sum_v[j] = s;
    cnt(j)   = c;
  }

  vec means(npopcomb);
  for (int j = 0; j < npopcomb; ++j) {
    // 0/0 -> NaN matches rowMeans(., na.rm=TRUE) on an all-NA row.
    // Cast back to double for the return value, matching R's rowMeans
    // `mean = sum / cnt` (LDOUBLE -> double on assignment).
    means(j) = (double)(sum_v[j] / (long double)cnt(j));
  }

  return Rcpp::List::create(_["means"] = means, _["cnt"] = cnt);
}


// [[Rcpp::export]]
arma::mat cpp_aftable_to_dstatden(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4,
                                  arma::vec& modelvec, arma::mat& usesnps, bool allsnps, int poly_only) {

  const int npopcomb = (int)p1.n_elem;
  const int nsnp     = (int)aftable.n_cols;
  validate_dstat_inputs(aftable, p1, p2, p3, p4, modelvec, usesnps, allsnps, npopcomb);

  mat den(npopcomb, nsnp);
  den.fill(NA_REAL);

  // Parallel over popcomb dimension. Each thread writes a distinct row
  // of den; aftable / usesnps are read-only shared.
  const int nthreads = admixtools_omp_threads();
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for (int j = 0; j < npopcomb; ++j) {
    const int i1 = (int)p1(j) - 1;
    const int i2 = (int)p2(j) - 1;
    const int i3 = (int)p3(j) - 1;
    const int i4 = (int)p4(j) - 1;
    const int m  = allsnps ? 0 : (int)modelvec(j) - 1;
    for (int i = 0; i < nsnp; ++i) {
      if (!allsnps && !usesnps(m, i)) continue;
      const double w = aftable(i1, i);
      const double x = aftable(i2, i);
      const double y = aftable(i3, i);
      const double z = aftable(i4, i);
      den(j, i) = (w + x - 2*w*x) * (y + z - 2*y*z);
    }
  }
  return den;
}



// [[Rcpp::export]]
arma::cube cpp_outer_array_mul(arma::mat& m1, arma::mat& m2) {

  int d1 = m1.n_cols;
  int d2 = m2.n_cols;
  int d3 = m1.n_rows;
  cube out(d1, d2, d3);
  for(int i = 0; i < d1; i++) {
    for(int j = 0; j < d2; j++) {
      for(int k = 0; k < d3; k++) {
        out(i, j, k) = m1(k, i) * m2(k, j);
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
arma::cube cpp_outer_array_plus(arma::mat& m1, arma::mat& m2) {

  int d1 = m1.n_cols;
  int d2 = m2.n_cols;
  int d3 = m1.n_rows;
  cube out(d1, d2, d3);
  for(int i = 0; i < d1; i++) {
    for(int j = 0; j < d2; j++) {
      for(int k = 0; k < d3; k++) {
        out(i, j, k) = m1(k, i) + m2(k, j);
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::cube cpp_outer_array_minus(arma::mat& m1, arma::mat& m2) {

  int d1 = m1.n_cols;
  int d2 = m2.n_cols;
  int d3 = m1.n_rows;
  cube out(d1, d2, d3);
  for(int i = 0; i < d1; i++) {
    for(int j = 0; j < d2; j++) {
      for(int k = 0; k < d3; k++) {
        out(i, j, k) = m1(k, i) - m2(k, j);
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::cube cpp_mats_to_f2_arr(arma::mat& afmat1, arma::mat& afmat2,
                              arma::mat& countmat1, arma::mat& countmat2,
                              bool apply_corr) {

  int nsnp = afmat1.n_rows;
  int nc1 = afmat1.n_cols;
  int nc2 = afmat2.n_cols;
  cube out = pow(cpp_outer_array_minus(afmat1, afmat2), 2);
  if(apply_corr) {
    mat denom1 = ones<mat>(nsnp, nc1);
    mat denom2 = ones<mat>(nsnp, nc2);
    denom1 = arma::max(denom1, countmat1 - 1);
    denom2 = arma::max(denom2, countmat2 - 1);
    // mat denom1 = countmat1-1;
    // mat denom2 = countmat2-1;
    mat corr1 = afmat1 % (1 - afmat1)/denom1;
    mat corr2 = afmat2 % (1 - afmat2)/denom2;
    out = out - cpp_outer_array_plus(corr1, corr2);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector row_prods(NumericMatrix x) {
  // Stolen from Rfast
  const int n=x.nrow();

  NumericVector f(n);
  mat X = mat(x.begin(), n, x.ncol(), false);
  colvec ff(f.begin(),n,false);
  ff = prod(X, 1);
  return f;
}


// Mathematically equivalent to the R one-liner
//   rowsum(gmat, popvec, na.rm = TRUE) / rowsum((!is.na(gmat))+0, popvec) / 2
// in one pass over gmat, one (npop x nsnp) allocation, no temporary
// (!is.na + 0) intermediate. Called per-block from f4blockdat_from_geno
// and qpdstat_geno; for a 60-pop / 1357-block / ~600-sample run this is
// 1357 calls each accumulating ~960K cells. The original R version
// allocates two (nind x nsnp) doubles per call (the !is.na cast and the
// rowsum output of it); this version allocates zero block-scale
// temporaries.
//
// popvec is the 1-based per-individual pop assignment vector that
// match() against the file's pop list produces upstream (length nind,
// values in 1..npop). NA / NaN genotypes in gmat are skipped via R's
// ISNAN macro (matches R's `!is.na(x)`, which returns FALSE for both
// NA and NaN but TRUE for +/-Inf). Using std::isfinite or arma::is_finite
// here would diverge from the R version on Inf inputs -- harmless on
// real genotype data (values in {0,1,2,NA}) but a real semantic
// difference, so we use ISNAN for strict equivalence.
//
// Returns an (npop x nsnp) matrix of reference-allele frequencies in
// [0, 1], with NaN where a (pop, snp) cell had zero valid genotypes
// (matches the original R version's 0/0 behavior).
//
// popvec is validated up front: values must be >= 1. NA_integer_
// (INT_MIN), zero, and negatives all trigger an error rather than
// undefined behavior at the sum_g(p, s) index step.
//
// [[Rcpp::export]]
arma::mat cpp_gmat_to_aftable(arma::mat& gmat, arma::ivec& popvec) {
  const int nind = (int)gmat.n_rows;
  const int nsnp = (int)gmat.n_cols;
  if ((int)popvec.n_elem != nind) {
    Rcpp::stop("popvec length (%d) != gmat n_rows (%d)",
               (int)popvec.n_elem, nind);
  }
  int npop = 0;
  for (int i = 0; i < nind; ++i) {
    if (popvec(i) < 1) Rcpp::stop("popvec contains values < 1 (or NA); expected 1-based pop indices");
    // Upper bound: popvec is a 1-based index from R-side match() against
    // a pop list of length <= nind. Larger values indicate a malformed
    // caller and would otherwise cause an enormous npop*nsnp allocation.
    if (popvec(i) > nind) {
      Rcpp::stop("popvec(%d) = %d exceeds nind (%d); expected 1-based index into a pop list with length <= nind",
                 i + 1, (int)popvec(i), nind);
    }
    if (popvec(i) > npop) npop = popvec(i);
  }

  mat sum_g = zeros<mat>(npop, nsnp);
  mat cnt   = zeros<mat>(npop, nsnp);

  // Parallel over SNPs. Each thread writes a distinct column of sum_g
  // and cnt (gmat and sum_g/cnt are arma::mat, column-major), so column
  // s is a contiguous slice owned by one thread; no cross-thread writes
  // collide. Rcpp::stop is not callable inside the parallel region;
  // the popvec validation above runs serially before we enter.
  const int nthreads = admixtools_omp_threads();
  #pragma omp parallel for num_threads(nthreads) schedule(static)
  for (int s = 0; s < nsnp; ++s) {
    for (int i = 0; i < nind; ++i) {
      const double g = gmat(i, s);
      if (!ISNAN(g)) {
        const int p = popvec(i) - 1;   // 1-based -> 0-based
        sum_g(p, s) += g;
        cnt(p, s)   += 1.0;
      }
    }
  }

  // Element-wise sum/cnt/2. 0/0 -> NaN, matching rowsum's behavior in R.
  return sum_g / cnt / 2.0;
}

