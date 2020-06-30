// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat cpp_opt_A(const arma::mat& B, const arma::mat& xvec, const arma::mat& qinv, int nr, double fudge) {
  // return A which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  // A: nr * rnk
  // B: rnk * nc
  // xvec: nr * nc
  // tdim: rnk * max(nr, nc)
  // coeffs: tdim * tdim
  // rhs: rnk * nc
  // qinv: nr*nc * nr*nc
  mat B2 = kron(eye(nr, nr), B);
  mat coeffs = B2 * qinv * B2.t();
  vec rhs = B2 * qinv * xvec;
  coeffs.diag() += fudge * sum(coeffs.diag());
  vec A2 = solve(coeffs, rhs, solve_opts::fast);
  return mat(&A2[0], B.n_rows, nr, false).t();
}

// [[Rcpp::export]]
arma::mat cpp_opt_B(const arma::mat& A, const arma::vec& xvec, const arma::mat& qinv, int nc, double fudge) {
  // return B which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  // A: nr * rnk
  // B: rnk * nc
  // xvec: nr * nc
  // tdim: rnk * max(nr, nc)
  // coeffs: tdim * tdim
  // rhs: rnk * nc
  // qinv: nr*nc * nr*nc
  mat A2 = kron(A, eye(nc, nc));
  mat coeffs = A2.t() * qinv * A2;
  vec rhs = A2.t() * qinv * xvec;
  coeffs.diag() += fudge * sum(coeffs.diag());
  vec B2 = solve(coeffs, rhs, solve_opts::fast);
  return mat(&B2[0], nc, A.n_cols, false).t();
}


// [[Rcpp::export]]
List cpp_qpadm_weights(const arma::mat& xmat, const arma::mat& qinv,
                       int rnk, double fudge = 0.0001, int iterations = 20,
                       bool constrained = false, Function qpsolve = R_NilValue) {

  if(rnk == 0) return Rcpp::List::create(_["weights"] = 1.0);
  mat U, V, x, y, rhs, lhs;
  vec s;
  svd(U, s, V, xmat);
  mat B = V.cols(0, rnk-1).t();
  mat A = xmat * B.t();
  vec xvec = vectorise(xmat.t());
  int nr = xmat.n_rows;
  int nc = xmat.n_cols;
  for(int i=0; i<iterations; i++) {
    A = cpp_opt_A(B, xvec, qinv, nr, fudge);
    B = cpp_opt_B(A, xvec, qinv, nc, fudge);
  }
  //vec w = solve(join_rows(A, ones(A.n_rows)).t(), join_cols(zeros(rnk), ones(1)));
  x = join_rows(A, ones(A.n_rows)).t();
  y = join_cols(zeros(rnk), ones(1));
  rhs = x.t() * x;
  lhs = x.t() * y;
  vec w;
  if(constrained) w = -as<vec>(qpsolve(rhs, lhs, -mat(nr, nr, fill::eye), zeros(nr)));
  else w = solve(rhs, lhs);
  vec weights = w / sum(w);

  return Rcpp::List::create(_["weights"] = weights, _["A"] = A, _["B"] = B);
}

// duplicated to have fast version that doesn't return List (used in covariance function)
vec cpp_qpadm_weights2(const arma::mat& xmat, const arma::mat& qinv,
                       int rnk, double fudge = 0.0001, bool constrained = false,
                       Function qpsolve = R_NilValue, int iterations = 20) {

  if(rnk == 0) return ones(1);
  mat U, V, x, y, rhs, lhs;
  vec s;
  svd(U, s, V, xmat);
  mat B = V.cols(0, rnk-1).t();
  mat A = xmat * B.t();
  vec xvec = vectorise(xmat.t());
  int nr = xmat.n_rows;
  int nc = xmat.n_cols;
  for(int i=0; i<iterations; i++) {
    A = cpp_opt_A(B, xvec, qinv, nr, fudge);
    B = cpp_opt_B(A, xvec, qinv, nc, fudge);
  }
  //vec w = solve(join_rows(A, ones(A.n_rows)).t(), join_cols(zeros(rnk), ones(1)));
  x = join_rows(A, ones(A.n_rows)).t();
  y = join_cols(zeros(rnk), ones(1));
  rhs = x.t() * x;
  lhs = x.t() * y;
  vec w;
  if(constrained) w = -as<vec>(qpsolve(rhs, lhs, -mat(nr, nr, fill::eye), zeros(nr)));
  else w = solve(rhs, lhs);
  return w / sum(w);
}


// [[Rcpp::export]]
arma::mat cpp_get_weights_covariance(arma::cube f4_lo, arma::mat qinv, arma::vec block_lengths,
                                     double fudge=0.0001, int boot = 0,
                                     bool constrained = false, Function qpsolve = R_NilValue) {
  int rnk = f4_lo.n_rows-1;
  int numreps = f4_lo.n_slices;
  mat wmat(numreps, f4_lo.n_rows);
  for(int i=0; i<numreps; i++) {
      wmat.row(i) = cpp_qpadm_weights2(f4_lo.slice(i), qinv, rnk, fudge,
                                       constrained, qpsolve).t();
  }
  if(!boot) wmat *= sqrt(numreps-1);
  return(cov(wmat));
}


// [[Rcpp::export]]
arma::mat cpp_aftable_to_dstatnum(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4) {

  mat num(p1.n_elem, aftable.n_cols);
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
    }
  }
  return num;
}


// [[Rcpp::export]]
arma::mat cpp_aftable_to_dstatden(arma::mat& aftable, arma::vec& p1, arma::vec& p2, arma::vec& p3, arma::vec& p4) {

  mat den(p1.n_elem, aftable.n_cols);
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
      den(j, i) = (w + x - 2*w*x) * (y + z - 2*y*z);
    }
  }
  return den;
}

