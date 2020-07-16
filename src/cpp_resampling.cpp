// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
IntegerVector cpp_get_block_lengths(IntegerVector chr, DoubleVector pos, double dist = 0.05) {

  double fpos = -1e20;
  int chrom;
  int lchrom = -1;
  int xsize = 0;
  int n = 0;
  int nsnps = chr.length();
  double gpos, dis;
  IntegerVector bsize(nsnps);

  for(int i = 0; i < nsnps; i++) {
    chrom = chr(i);
    gpos = pos(i);
    dis = gpos - fpos;
    if((chrom != lchrom) || (dis >= dist)) {
      if(xsize > 0) {
        bsize(n) = xsize;
        n++;
      }
      lchrom = chrom;
      fpos = gpos;
      xsize = 0;
    }
    xsize++;
  }
  if (xsize > 0) {
    bsize(n) = xsize;
    n++;
  }
  return bsize[Rcpp::Range(0, n-1)];
}


// [[Rcpp::export]]
List cpp_jack_vec_stats(NumericVector loo_vec, NumericVector block_lengths) {
  // input is a vector of leave-one-out estimates
  // output is list with jackknife mean and covariance
  // should give same results as 'jack_arr_stats' and 'jack_mat_stats'

  block_lengths = block_lengths[!is_na(loo_vec)];
  loo_vec = loo_vec[!is_na(loo_vec)];

  NumericVector w = 1-block_lengths/sum(block_lengths);
  double tot = sum(loo_vec * w)/sum(w);
  double est = mean(loo_vec);
  NumericVector y = sum(block_lengths)/block_lengths;
  NumericVector xtau = (tot * y - loo_vec * (y-1) - est) / sqrt(y-1);
  double var = mean(pow(xtau, 2.0));

  //return List::create(est, var);
  return Rcpp::List::create(_["est"] = est, _["var"] = var);
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
arma::cube cpp_mats_to_f2_arr(arma::mat& afmat1, arma::mat& afmat2, arma::mat& countmat1, arma::mat& countmat2) {

  int nsnp = afmat1.n_rows;
  int nc1 = afmat1.n_cols;
  int nc2 = afmat2.n_cols;
  mat denom1 = ones<mat>(nsnp, nc1);
  mat denom2 = ones<mat>(nsnp, nc2);
  denom1 = arma::max(denom1, countmat1 - 1);
  denom2 = arma::max(denom2, countmat2 - 1);
  mat pq1 = afmat1 % (1 - afmat1)/denom1;
  mat pq2 = afmat2 % (1 - afmat2)/denom2;
  cube out = pow(cpp_outer_array_minus(afmat1, afmat2), 2) - cpp_outer_array_plus(pq1, pq2);
  return out;
}











