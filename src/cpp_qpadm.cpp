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




// [[Rcpp::export]]
arma::vec cpp_is_polymorphic(arma::mat geno) {

  int nr = geno.n_rows;
  int nc = geno.n_cols;
  double first, val;
  bool isfirst = true;
  vec out = zeros<vec>(nr);
  for(int i = 0; i < nr; i++) {
    isfirst = true;
    for(int j = 0; j < nc; j++) {
      val = geno(i,j);
      if(is_finite(val)) {
        if(isfirst) {
          first = val;
          isfirst = false;
        } else {
          if(first != val) {
            out(i) = 1;
            break;
          }
        }
      }
    }
  }
  return out;
}

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



