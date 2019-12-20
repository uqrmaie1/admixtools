// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat cpp_opt_A(const arma::mat& B, const arma::mat& xmat, const arma::mat& qinv, double fudge) {
// return A which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
// A: nr * rnk
// B: rnk * nc
// xmat: nr * nc
// tdim: rnk * max(nr, nc)
// coeffs: tdim * tdim
// rhs: rnk * nc
// qinv: nr*nc * nr*nc
  mat B2 = kron(eye(xmat.n_rows, xmat.n_rows), B);
  mat coeffs = B2 * qinv * B2.t();
  vec rhs = B2 * qinv * vectorise(xmat.t());
  coeffs.diag() += fudge;
  vec A2 = solve(coeffs, rhs, solve_opts::fast);
  return mat(&A2[0], B.n_rows, xmat.n_rows, false).t();
}

// [[Rcpp::export]]
arma::mat cpp_opt_B(const arma::mat& A, const arma::mat& xmat, const arma::mat& qinv, double fudge) {
  // return B which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  // A: nr * rnk
  // B: rnk * nc
  // xmat: nr * nc
  // tdim: rnk * max(nr, nc)
  // coeffs: tdim * tdim
  // rhs: rnk * nc
  // qinv: nr*nc * nr*nc
  mat A2 = kron(A, eye(xmat.n_cols, xmat.n_cols));
  mat coeffs = A2.t() * qinv * A2;
  vec rhs = A2.t() * qinv * vectorise(xmat.t());
  coeffs.diag() += fudge;
  vec B2 = solve(coeffs, rhs, solve_opts::fast);
  return mat(&B2[0], xmat.n_cols, A.n_cols, false).t();
}


// [[Rcpp::export]]
arma::vec cpp_qpadm_weights(const arma::mat& xmat, const arma::mat& qinv, int rnk, double fudge=0.0001, int iterations=20) {

  mat U;
  vec s;
  mat V;
  svd(U, s, V, xmat);
  mat B = V.cols(0, rnk-1).t();
  mat A = xmat * B.t();
  for(int i=0; i<iterations; i++) {
    A = cpp_opt_A(B, xmat, qinv, fudge);
    B = cpp_opt_B(A, xmat, qinv, fudge);
  }
  vec w = solve(join_rows(A, ones(A.n_rows)).t(), join_cols(zeros(rnk), ones(1)));
  return w / sum(w);
}


// [[Rcpp::export]]
arma::mat cpp_get_weights_covariance(arma::cube f4_blocks, arma::mat qinv, double fudge=0.0001) {

  int rnk = f4_blocks.n_rows-1;
  int numblocks = f4_blocks.n_slices;
  mat wmat(numblocks, f4_blocks.n_rows);
  for(int i=0; i<numblocks; i++) {
      wmat.row(i) = cpp_qpadm_weights(f4_blocks.slice(i), qinv, rnk, fudge).t();
  }
  vec jackmeans = mean(wmat, 0).t();
  mat wmat2 = -wmat.t();
  mat mnc = wmat2.each_col() + jackmeans;
  return (numblocks-1) / (double)numblocks * (mnc * mnc.t());
}







