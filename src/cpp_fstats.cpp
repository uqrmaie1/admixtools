// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
      if(is_finite(num(j, i))) cnt(j) += 1;
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
  NumericVector uni;

  double w, x, y, z;
  int i1, i2, i3, i4, m, valid = 0;
  for(int j = 0; j < p1.n_elem; j++) {
    i1 = p1(j)-1;
    i2 = p2(j)-1;
    i3 = p3(j)-1;
    i4 = p4(j)-1;
    if(!allsnps) m = modelvec(j)-1;
    for(int i = 0; i < aftable.n_cols; i++) {
      if(allsnps || usesnps(m, i)) {
        w = aftable(i1, i);
        x = aftable(i2, i);
        y = aftable(i3, i);
        z = aftable(i4, i);
        if(!(poly_only && allsnps)) {
          valid = 1;
        } else {
          // original admixtools only excludes polymorphic SNPs if they are all 0 or all 1
          uni = na_omit(unique(NumericVector::create(w, x, y, z)));
          if(poly_only == 0 || uni.length() > 1 || (poly_only == 2 && (max(uni) > 0.0001 && max(uni) < 0.9999))) {
            valid = 1;
          }
        }
        if(valid) {
          num(j, i) = (w - x) * (y - z);
          valid = 0;
        }
        if(is_finite(num(j, i))) cnt(j) += 1;
      }
    }
  }
  //return num;
  return Rcpp::List::create(_["num"] = num, _["cnt"] = cnt);
}


// [[Rcpp::export]]
arma::mat cpp_aftable_to_dstatden(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4,
                                  arma::vec& modelvec, arma::mat& usesnps, bool allsnps, int poly_only) {

  mat den(p1.n_elem, aftable.n_cols);
  den.fill(NA_REAL);

  double w, x, y, z;
  int i1, i2, i3, i4, m;
  for(int j = 0; j < p1.n_elem; j++) {
    i1 = p1(j)-1;
    i2 = p2(j)-1;
    i3 = p3(j)-1;
    i4 = p4(j)-1;
    if(!allsnps) m = modelvec(j)-1;
    for(int i = 0; i < aftable.n_cols; i++) {
      if(allsnps || usesnps(m, i)) {
        w = aftable(i1, i);
        x = aftable(i2, i);
        y = aftable(i3, i);
        z = aftable(i4, i);
        den(j, i) = (w + x - 2*w*x) * (y + z - 2*y*z);
      }
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
// but in one pass over gmat, one (npop x nsnp) allocation, no temporary
// (!is.na + 0) intermediate. Called per-block from f4blockdat_from_geno
// and qpdstat_geno; for a 60-pop / 1357-block / ~600-sample run this is
// 1357 calls each accumulating ~960K cells. The original R version
// allocates two (nind x nsnp) doubles per call (the !is.na cast and the
// rowsum output of it); this version allocates zero block-scale
// temporaries.
//
// popvec is the 1-based per-individual pop assignment vector that
// match() against the file's pop list produces upstream (length nind,
// values in 1..npop). NA / missing genotypes in gmat are skipped.
// Returns an (npop x nsnp) matrix of reference-allele frequencies in
// [0, 1], with NaN where a (pop, snp) cell had zero valid genotypes
// (matches the original R version's 0/0 behavior).
//
// [[Rcpp::export]]
arma::mat cpp_gmat_to_aftable(arma::mat& gmat, arma::ivec& popvec) {
  int nind = gmat.n_rows;
  int nsnp = gmat.n_cols;
  int npop = 0;
  for(int i = 0; i < nind; i++) if(popvec(i) > npop) npop = popvec(i);

  mat sum_g = zeros<mat>(npop, nsnp);
  mat cnt   = zeros<mat>(npop, nsnp);

  for(int s = 0; s < nsnp; s++) {
    for(int i = 0; i < nind; i++) {
      double g = gmat(i, s);
      if(std::isfinite(g)) {
        int p = popvec(i) - 1;   // 1-based -> 0-based
        sum_g(p, s) += g;
        cnt(p, s)   += 1.0;
      }
    }
  }

  // Element-wise sum/cnt/2. 0/0 -> NaN, matching rowsum's behavior in R.
  return sum_g / cnt / 2.0;
}

