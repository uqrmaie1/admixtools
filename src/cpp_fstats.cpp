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
                             arma::vec& modelvec, arma::mat& usesnps, bool allsnps) {

  // aftable is npop x nsnp
  // num is npopcomb x nsnp
  mat num(p1.n_elem, aftable.n_cols);
  num.fill(NA_REAL);
  vec cnt = zeros(p1.n_elem);

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
        num(j, i) = (w - x) * (y - z);
        if(is_finite(num(j, i))) cnt(j) += 1;
      }
    }
  }
  //return num;
  return Rcpp::List::create(_["num"] = num, _["cnt"] = cnt);
}


// [[Rcpp::export]]
arma::mat cpp_aftable_to_dstatden(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4,
                                  arma::vec& modelvec, arma::mat& usesnps, bool allsnps) {

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

