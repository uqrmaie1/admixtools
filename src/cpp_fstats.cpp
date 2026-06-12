// [[Rcpp::depends(RcppArmadillo)]]
//
// The C++ standard (CXX_STD = CXX17) and the OpenMP flags
// ($(SHLIB_OPENMP_CXXFLAGS)) are supplied by src/Makevars and
// src/Makevars.win, which are the single source of truth for a package
// (R CMD INSTALL) build. The Rcpp::plugins(cpp11)/plugins(openmp) attributes
// only take effect under Rcpp::sourceCpp, so they are intentionally omitted
// here: keeping them would imply the package build depends on them, and the
// stale cpp11 plugin also contradicted CXX17.
//
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

// Shared OpenMP parallel-for boilerplate for the four kernels. Expands to the
// parallel-for pragma (team size clamped to >= 1) under OpenMP, or to a no-op
// that suppresses the otherwise-unused `nthreads` without it. `sched` is a bare
// schedule token (guided/static). Expansion is identical to the hand-written
// pragmas, so kernel output is unchanged.
#define ADMIXTOOLS_STRINGIFY_(x) #x
#define ADMIXTOOLS_STRINGIFY(x) ADMIXTOOLS_STRINGIFY_(x)
#ifdef _OPENMP
  #define ADMIXTOOLS_OMP_PARALLEL_FOR(nthreads, sched)                     \
      _Pragma(ADMIXTOOLS_STRINGIFY(omp parallel for                        \
              num_threads((nthreads) > 0 ? (nthreads) : 1) schedule(sched)))
#else
  #define ADMIXTOOLS_OMP_PARALLEL_FOR(nthreads, sched) (void)(nthreads);
#endif

// OpenMP + fork() hazard mitigation. libgomp's thread pool, once created in a
// process, is inherited across fork() but its worker threads are not, so a
// forked child that enters an OpenMP region deadlocks on pool locks held by
// threads that no longer exist in the child. This bites a caller who runs an
// f-stat under plan(sequential) -- which executes these kernels in the main
// process and creates the pool -- and then switches to plan(multicore), which
// forks. Draining the pool in a pthread_atfork prepare handler makes every
// fork start clean; the next parallel region in parent or child rebuilds it.
// No effect on plan(sequential) or plan(multisession) (whose workers are fresh
// spawned processes, not forks), or when OpenMP is unavailable. The drain uses
// omp_pause_resource_all, an OpenMP 5.0 routine. The guard enables it when the
// compiler advertises OpenMP 5.0 (_OPENMP >= 201811, e.g. recent clang / libomp,
// including macOS) or for real gcc >= 9, whose libgomp ships the symbol even
// though gcc pins _OPENMP to 4.5. It deliberately does NOT key on
// __clang_major__: Apple clang's version does not track the linked libomp, so a
// new-clang / old-libomp pairing could wrongly claim the symbol and fail to
// compile. On runtimes outside both branches the handler is omitted and callers
// should prefer plan(multisession) over plan(multicore).
//
// Handler lifecycle. pthread_atfork handlers cannot be unregistered, so the
// image holding this handler must stay mapped: admixtools_omp_atfork_active_()
// tells .onUnload() (R/admixtools.R) not to unload the .so while a handler is
// live. That also prevents re-registration, since R then reuses the loaded DLL
// on a reload instead of mapping a second image, so the lambda below runs once
// per process and the handler is never stacked or left dangling.

// True while a live at-fork handler exists in this image. Outside the OpenMP
// guard so the accessor compiles on every build; set true only in the gated
// lambda, so it stays false when built without OpenMP.
static bool admixtools_omp_atfork_active = false;

#if defined(_OPENMP) && \
    (_OPENMP >= 201811 || \
     (defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 9))
#include <pthread.h>
static void admixtools_omp_drain_pool_before_fork() {
  omp_pause_resource_all(omp_pause_hard);
}
static const int admixtools_omp_atfork_registered = []() {
  // Register the drain handler once, when this image loads (see the lifecycle
  // note above for why once is enough).
  pthread_atfork(admixtools_omp_drain_pool_before_fork, nullptr, nullptr);
  admixtools_omp_atfork_active = true;   // a handler is live in this image
  return 0;
}();
#endif

// Liveness flag for .onUnload(): when true, the .so must stay mapped or the
// registered at-fork handler would dangle and crash the next fork(). Always
// defined (false on non-OpenMP builds).
// [[Rcpp::export]]
bool admixtools_omp_atfork_active_() {
  return admixtools_omp_atfork_active;
}

// The OpenMP team size is no longer resolved here. The caller passes an
// explicit `nthreads` (see the kernels below), resolved once on the R side
// in the main process where parallelly::availableCores() sees the true
// machine / cgroup / scheduler core count and future::nbrOfWorkers() reflects
// the active plan() (f4blockdat_from_geno / f3blockdat_from_geno in R/io.R).
// Resolving it there -- rather than calling availableCores() from inside a
// kernel -- means a dispatched or forked future worker never re-resolves the
// count (which would either oversubscribe via an inherited cache or collapse
// to 1), and the thread budget is split across the future workers so total
// active threads stay at or below the core count under every plan.

// Polymorphism thresholds for the poly_only == 2 rule. NOTE: a sibling
// polymorphism filter lives in R, in the `usesnps` construction inside
// f4blockdat_from_geno / f3blockdat_from_geno (R/io.R). That R filter runs on
// the !allsnps path across a model's population groups, whereas
// poly_only_keep() below runs on the allsnps path across the four popcomb
// members, so the two are deliberately separate code paths encoding the same
// biological intent. If you change the polymorphism definition in one place,
// review the other.
static const double POLY_ONLY_LO = 0.0001;
static const double POLY_ONLY_HI = 0.9999;

// Per-popcomb 0-based population indices, plus the usesnps model row. Shared
// by the three dstat kernels.
struct dstat_popcomb { int i1, i2, i3, i4, m; };

static inline dstat_popcomb decode_popcomb(
    const arma::vec& p1, const arma::vec& p2, const arma::vec& p3,
    const arma::vec& p4, const arma::vec& modelvec, bool allsnps, int j) {
  dstat_popcomb ix;
  ix.i1 = (int)p1(j) - 1;   // validate_dstat_inputs guarantees [1, npop],
  ix.i2 = (int)p2(j) - 1;   // so these (int) casts are in range.
  ix.i3 = (int)p3(j) - 1;
  ix.i4 = (int)p4(j) - 1;
  ix.m  = allsnps ? 0 : (int)modelvec(j) - 1;   // validated in [1, nmodel]
  return ix;
}

// Load the four allele frequencies for SNP i of one popcomb, returning false
// (skip this SNP) when the !allsnps usesnps mask excludes it. Shared by the
// three dstat kernels so the usesnps gate and the four reads live in one
// place. No Rcpp R-API allocations, so it is safe inside the parallel region.
static inline bool load_quad(
    const arma::mat& aftable, const arma::mat& usesnps,
    const dstat_popcomb& ix, int i, bool allsnps,
    double& w, double& x, double& y, double& z) {
  if (!allsnps && usesnps(ix.m, i) == 0.0) return false;
  w = aftable(ix.i1, i);
  x = aftable(ix.i2, i);
  y = aftable(ix.i3, i);
  z = aftable(ix.i4, i);
  return true;
}

// Serial pre-pass shared by the three cpp_aftable_to_dstat* functions, run
// before the parallel region because Rcpp::stop is not callable from OpenMP
// code: a bad index would otherwise become a wild read or an Armadillo bounds
// exception thrown across the omp boundary (std::terminate). Checks, with the
// offending row index:
//   * p1..p4 (and modelvec when !allsnps) are long enough to index, so this
//     pass and the kernels never read past the vectors;
//   * p1..p4 are finite 1-based indices in [1, npop] -- the range test is on
//     the double value, BEFORE the (int) cast, so a finite-but-out-of-int
//     value cannot slip through an undefined float->int conversion into the
//     region;
//   * when !allsnps, modelvec is a finite 1-based index in [1, nmodel] and
//     usesnps has at least nsnp columns, so usesnps(m, i) is in bounds for
//     every i. O(npopcomb), negligible vs the inner i loop.
//
// The length/column guards always run; they bound what the parallel region
// dereferences. The per-element scan is block-invariant in the hot paths, so
// callers run it once (validate = true on block 1) and skip it after with
// validate = false, which then trusts the caller's indices.
static inline void validate_dstat_inputs(
    const arma::mat& aftable,
    const arma::vec& p1, const arma::vec& p2,
    const arma::vec& p3, const arma::vec& p4,
    const arma::vec& modelvec, const arma::mat& usesnps,
    bool allsnps, int npopcomb, int nsnp, bool validate = true) {
  const int npop   = (int)aftable.n_rows;
  const int nmodel = (int)usesnps.n_rows;
  // npopcomb is p1.n_elem by construction (set by the caller); only p2..p4
  // (and modelvec, below) can be shorter, which would index past their ends.
  if ((int)p2.n_elem < npopcomb || (int)p3.n_elem < npopcomb ||
      (int)p4.n_elem < npopcomb) {
    Rcpp::stop("p2/p3/p4 shorter than p1 (%d popcombs)", npopcomb);
  }
  if (!allsnps) {
    if ((int)modelvec.n_elem < npopcomb) {
      Rcpp::stop("modelvec shorter than the number of popcombs (%d)", npopcomb);
    }
    if ((int)usesnps.n_cols < nsnp) {
      Rcpp::stop("usesnps has %d columns but aftable has %d SNPs", (int)usesnps.n_cols, nsnp);
    }
  }
  // The cheap bounds guards above always run; the per-element scan below is the
  // redundant-across-blocks part the hot paths skip with validate = false.
  if (!validate) return;
  for (int j = 0; j < npopcomb; ++j) {
    if (!R_FINITE(p1(j)) || !R_FINITE(p2(j)) ||
        !R_FINITE(p3(j)) || !R_FINITE(p4(j))) {
      Rcpp::stop("p1/p2/p3/p4 contains NA / NaN / Inf at index %d", j + 1);
    }
    // Range-check the double value before the (int) cast in decode_popcomb.
    if (p1(j) < 1.0 || p1(j) > (double)npop || p2(j) < 1.0 || p2(j) > (double)npop ||
        p3(j) < 1.0 || p3(j) > (double)npop || p4(j) < 1.0 || p4(j) > (double)npop) {
      Rcpp::stop("p1/p2/p3/p4 index out of range [1..%d] at j=%d", npop, j + 1);
    }
    if (!allsnps) {
      if (!R_FINITE(modelvec(j))) {
        Rcpp::stop("modelvec contains NA / NaN / Inf at index %d", j + 1);
      }
      if (modelvec(j) < 1.0 || modelvec(j) > (double)nmodel) {
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
    if (maxv > POLY_ONLY_LO && maxv < POLY_ONLY_HI) return true;
  }
  return false;
}

// [[Rcpp::export]]
List cpp_aftable_to_dstatnum(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4,
                             arma::vec& modelvec, arma::mat& usesnps, bool allsnps, int poly_only,
                             int nthreads = 1, bool validate = true) {

  // aftable is npop x nsnp
  // num is npopcomb x nsnp
  mat num(p1.n_elem, aftable.n_cols);
  num.fill(NA_REAL);
  vec cnt = zeros(p1.n_elem);

  const int npopcomb = (int)p1.n_elem;
  const int nsnp     = (int)aftable.n_cols;
  validate_dstat_inputs(aftable, p1, p2, p3, p4, modelvec, usesnps, allsnps, npopcomb, nsnp, validate);

  // Parallel over the popcomb dimension. Each thread writes a distinct row of
  // num and a distinct slot of cnt, so no reduction is needed; the helpers
  // are pure C++ (no Rcpp allocations) so they are safe inside the region.
  // schedule(guided) keeps large chunks for the dense (many-popcomb) case
  // while still spreading small/medium npopcomb across threads. The inner SNP
  // loop stays serial per row, which is what makes the result bit-identical
  // regardless of nthreads, so a single statistic (npopcomb == 1) is
  // intentionally not split further (an OpenMP reduction would reorder the
  // sum and break that contract).
  ADMIXTOOLS_OMP_PARALLEL_FOR(nthreads, guided)
  for (int j = 0; j < npopcomb; ++j) {
    const dstat_popcomb ix = decode_popcomb(p1, p2, p3, p4, modelvec, allsnps, j);
    double row_cnt = 0.0;
    double w, x, y, z;
    for (int i = 0; i < nsnp; ++i) {
      if (!load_quad(aftable, usesnps, ix, i, allsnps, w, x, y, z)) continue;
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
// indexing. On 64-bit builds this does not arise: RcppArmadillo's
// compiler_setup.hpp auto-defines ARMA_64BIT_WORD unless ARMA_32BIT_WORD is
// set (this package never sets it) or the platform is 32-bit, so arma::uword
// is u64 and a matrix has no 2^31-element ceiling at realistic admixtools
// scales. The int loop counters here do assume npopcomb and nsnp stay below
// 2^31, which holds comfortably (millions of popcombs, blocks of a few
// thousand SNPs). The binding constraint is memory, not index width.)
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
                                      bool allsnps, int poly_only,
                                      int nthreads = 1, bool validate = true) {

  // sum_v in long double matches R's rowMeans precision; cnt is an
  // integer count where 64-bit double is more than enough.
  const int npopcomb = (int)p1.n_elem;
  const int nsnp     = (int)aftable.n_cols;
  validate_dstat_inputs(aftable, p1, p2, p3, p4, modelvec, usesnps, allsnps, npopcomb, nsnp, validate);

  std::vector<long double> sum_v(npopcomb, 0.0L);
  vec cnt = zeros<vec>(npopcomb);

  // Parallel over the popcomb dimension. Each thread writes a distinct
  // sum_v[j] and cnt(j); aftable / usesnps are read-only shared. The inner
  // SNP sum stays serial per row, so the per-row long-double accumulation
  // order -- and thus the output -- is identical for any nthreads. See the
  // schedule(guided) / single-statistic note on cpp_aftable_to_dstatnum.
  ADMIXTOOLS_OMP_PARALLEL_FOR(nthreads, guided)
  for (int j = 0; j < npopcomb; ++j) {
    const dstat_popcomb ix = decode_popcomb(p1, p2, p3, p4, modelvec, allsnps, j);
    long double s = 0.0L;
    double      c = 0.0;
    double w, x, y, z;
    for (int i = 0; i < nsnp; ++i) {
      if (!load_quad(aftable, usesnps, ix, i, allsnps, w, x, y, z)) continue;
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
                                  arma::vec& modelvec, arma::mat& usesnps, bool allsnps, int poly_only,
                                  int nthreads = 1, bool validate = true) {

  const int npopcomb = (int)p1.n_elem;
  const int nsnp     = (int)aftable.n_cols;
  validate_dstat_inputs(aftable, p1, p2, p3, p4, modelvec, usesnps, allsnps, npopcomb, nsnp, validate);

  mat den(npopcomb, nsnp);
  den.fill(NA_REAL);

  // Parallel over popcomb dimension. Each thread writes a distinct row of
  // den; aftable / usesnps are read-only shared. schedule(guided) as in the
  // numerator kernels.
  ADMIXTOOLS_OMP_PARALLEL_FOR(nthreads, guided)
  for (int j = 0; j < npopcomb; ++j) {
    const dstat_popcomb ix = decode_popcomb(p1, p2, p3, p4, modelvec, allsnps, j);
    double w, x, y, z;
    for (int i = 0; i < nsnp; ++i) {
      if (!load_quad(aftable, usesnps, ix, i, allsnps, w, x, y, z)) continue;
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
// values in 1..npop). gmat is an integer dosage matrix; missing cells
// hold NA_INTEGER (INT_MIN) and are skipped via `g != NA_INTEGER`, the
// integer analogue of R's `!is.na(x)`. A double matrix passed by an
// older caller is coerced to integer at the R/C++ boundary (Inf/NaN ->
// NA_INTEGER with R's "NAs introduced by coercion" warning); real
// genotype data is {0,1,2,NA} so this is value-preserving.
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
arma::mat cpp_gmat_to_aftable(arma::imat& gmat, arma::ivec& popvec, int nthreads = 1) {
  // gmat is an integer dosage matrix in {0, 1, 2, NA_INTEGER}. arma::imat
  // (vs the previous arma::mat) lets every reader's output pass through
  // without R coercing it to a double matrix per block: the BED/
  // PACKEDANCESTRYMAP/EIGENSTRAT readers return IntegerMatrix and PFILE's
  // pgenlibr::ReadIntList is already integer. Missing is checked with
  // == NA_INTEGER (INT_MIN), the integer analogue of R's !is.na.
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
  // and cnt (sum_g/cnt are arma::mat, column-major), so column s is a
  // contiguous slice owned by one thread; no cross-thread writes collide.
  // Rcpp::stop is not callable inside the parallel region; the popvec
  // validation above runs serially before we enter. nthreads is the
  // caller-supplied budget (see f4blockdat_from_geno in R/io.R).
  // schedule(static): per-column work is uniform (one pass over nind).
  ADMIXTOOLS_OMP_PARALLEL_FOR(nthreads, static)
  for (int s = 0; s < nsnp; ++s) {
    for (int i = 0; i < nind; ++i) {
      const int g = gmat(i, s);
      if (g != NA_INTEGER) {
        const int p = popvec(i) - 1;   // 1-based -> 0-based
        sum_g(p, s) += (double)g;
        cnt(p, s)   += 1.0;
      }
    }
  }

  // Element-wise sum/cnt/2. 0/0 -> NaN, matching rowsum's behavior in R.
  return sum_g / cnt / 2.0;
}

