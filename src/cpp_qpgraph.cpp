// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <RcppEigen.h>

#include <iostream>
#include <sstream>
#include <string>

//#include <Eigen/Dense>
//#include "eiquadprog.h"

using namespace Rcpp;
using namespace arma;
//using namespace Eigen;

// #include "eigen-qp.hpp"
//
// namespace EigenQP {
// template<typename Scalar, int NVars, int NIneq>
// void quadprog(Eigen::Matrix<Scalar,NVars,NVars> &Q, Eigen::Matrix<Scalar,NVars,1> &c,
//              Eigen::Matrix<Scalar,NIneq,NVars> &A, Eigen::Matrix<Scalar,NIneq,1> &b,
//              Eigen::Matrix<Scalar,NVars,1> &x);
//
//
// template<typename Scalar, int NVars, int NIneq>
// void quadprog2(Eigen::Matrix<Scalar,NVars,NVars> &Q);
// }


struct QPException : public std::exception
{
  std::string s;
  QPException(const std::string &ss) : s(ss) {}
  ~QPException() throw () {} // Updated
  const char* what() const throw() { return s.c_str(); }
};



// [[Rcpp::export]]
arma::vec cpp_opt_edge_lengths(const arma::mat& ppwts_2d, const arma::mat& ppinv,
                               const arma::vec& f3_est, Function qpsolve,
                               const arma::vec& lower, const arma::vec& upper, double fudge) {
  mat pppp = ppwts_2d.t() * ppinv;
  mat cc = pppp * ppwts_2d;
  int nc = cc.n_cols;
  double trace = 0.0;
  for(int i = 0; i < nc; i++) {
    for(int j = i; j < nc; j++) {
      cc(i,j) = cc(j,i);
      if(i == j) trace += cc(j,i);
    }
  }
  for(int i = 0; i < nc; i++) {
    //cc(i,i) += fudge;
    cc(i,i) += fudge * trace;
  }
  vec q1 = pppp * f3_est;
  mat CI(nc, nc, fill::eye);
  return as<vec>(qpsolve(cc, q1, join_horiz(CI, -CI), join_vert(lower, -upper)));
}



// [[Rcpp::export]]
arma::mat cpp_fill_pwts(arma::mat& pwts, const arma::vec& weights,
                        const arma::mat& path_edge_table,
                        const arma::mat& path_admixedge_table, int numpaths) {
  // puts weights onto pwts, using index matrix and vectors
  // returns pwts to be consistent with R function. However, it edits pwts in place. Used to return NULL.

  if(weights.n_elem == 0) return pwts;

  vec pathweights(numpaths, fill::ones);
  vec pwtsfill(pwts.n_elem, fill::zeros);
  vec pwtsfillind(pwts.n_elem, fill::zeros);

  int path, admixedge, row, column, elem;
  for(int i=0; i<path_admixedge_table.n_rows; i++) {
    path = path_admixedge_table(i, 0)-1;
    admixedge = path_admixedge_table(i, 1)-1;
    if(admixedge % 2 == 0) pathweights(path) *= weights(admixedge/2);
    else pathweights(path) *= 1-weights(admixedge/2);
  }
  for(int i=0; i<path_edge_table.n_rows; i++) {
    row = path_edge_table(i, 2)-1;
    column = path_edge_table(i, 4)-1;
    path = path_edge_table(i, 0)-1;
    elem = row + pwts.n_rows*column;
    pwtsfill(elem) += pathweights(path);
    pwtsfillind(elem) = 1;
  }
  for(int i=0; i<pwts.n_elem; i++) {
    if(pwtsfillind(i) > 0) pwts(i) = pwtsfill(i);
  }
  return pwts;
}



// [[Rcpp::export]]
double cpp_optimweightsfun(arma::vec weights, List args) {
  mat pwts = args[0];
  mat ppinv = args[1];
  vec f3_est = args[2];
  mat path_edge_table = args[3];
  mat path_admixedge_table = args[4];
  int numpaths = args[5];
  // args[6] is cmb matrix with column combinations; not needed here
  Function qpsolve = args[7];
  vec lower = args[8];
  vec upper = args[9];
  double fudge = args[10];
  double nr = pwts.n_rows;
  double nc = pwts.n_cols;
  cpp_fill_pwts(pwts, weights, path_edge_table, path_admixedge_table, numpaths);
  mat ppwts_2d(nc*(nc+1)/2, nr);
  for(int i=0; i<nr; i++) {
    int c=0;
    for(int j=0; j<nc; j++) {
      for(int k=j; k<nc; k++) {
        ppwts_2d(c, i) = pwts(i, j) * pwts(i, k);
        c++;
      }
    }
  }
  vec q2 = cpp_opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper, fudge);
  vec w2 = (ppwts_2d * q2) - f3_est;
  vec lik = w2.t() * ppinv * w2;
  return lik(0);
}

int choose2(int k) {
  return k*(k-1)/2;
}

// [[Rcpp::export]]
NumericVector cpp_get_pairindex(const NumericVector perm) {

  int numpop = perm.size();
  int numpair = choose2(numpop);
  int c, c1, c2, num, i1, i2, v1, v2;
  NumericVector new_order(numpair);
  NumericVector perm2(numpop-1);
  c = 0;
  for(int i=0; i<numpop; i++) {
    num = perm(i);
    if(num != 1) {
      perm2(c) = perm(i);
      c++;
    }
  }
  c1 = 0;
  c2 = 0;
  for(int i=0; i<numpair; i++) {
    if(c2 == numpop-1) {
      c1++;
      c2 = c1;
    }
    i1 = perm2(c1);
    i2 = perm2(c2);
    v1 = std::min(i1, i2);
    v2 = std::max(i1, i2);
    new_order(i) = choose2(numpop) - choose2(numpop-v1+1) + v2 - numpop;
    c2++;
  }
  return new_order;
}

// Rcpp optimization using Roptim (cpp_lbfgsb / Optimfunctor in c++) is not actually faster than using R optim.

/*
// [[Rcpp::depends(roptim)]]
#include <roptim.h>
using namespace roptim;

class Optimfunctor : public Functor {

public:
  Optimfunctor(const arma::mat &pwts, const arma::mat &ppinv, const arma::vec &f3_est, const arma::mat &path_edge_table,
               const arma::mat &path_admixedge_table, Function qpsolve, int numpaths) : qpsolve2(qpsolve) {
    thispwts = pwts;
    thisppinv = ppinv;
    thisf3_est = f3_est;
    thispath_edge_table = path_edge_table;
    thispath_admixedge_table = path_admixedge_table;
    thisnumpaths = numpaths;
  }

  double operator()(const arma::vec &weights) override {

    double nr = thispwts.n_rows;
    double nc = thispwts.n_cols;

    cpp_fill_pwts(thispwts, weights, thispath_edge_table, thispath_admixedge_table, thisnumpaths);

    arma::mat ppwts_2d(nc*(nc+1)/2, nr);

    for(int i=0; i<nr; i++) {
      int c=0;
      for(int j=0; j<nc; j++) {
        for(int k=j; k<nc; k++) {
          ppwts_2d(c, i) = thispwts(i, j) * thispwts(i, k);
          c++;
        }
      }
    }

    arma::vec q2 = cpp_opt_edge_lengths(ppwts_2d, thisppinv, thisf3_est, qpsolve2);
    arma::vec w2 = (ppwts_2d * q2) - thisf3_est;
    arma::vec lik = w2.t() * thisppinv * w2;
    return lik(0);
  }

  private:
    arma::mat thispwts;
    arma::mat thisppinv;
    arma::vec thisf3_est;
    arma::mat thispath_edge_table;
    arma::mat thispath_admixedge_table;
    Function qpsolve2;
    int thisnumpaths;
};


// [[Rcpp::export]]
List cpp_lbfgsb(arma::vec weights, const arma::mat &pwts, const arma::mat &ppinv,
                const arma::vec &f3_est, const arma::mat &path_edge_table,
                const arma::mat &path_admixedge_table, int numpaths, Function qpsolve) {

  Optimfunctor optfu(pwts, ppinv, f3_est, path_edge_table,
                     path_admixedge_table, qpsolve, numpaths);
  Roptim<Optimfunctor> opt("L-BFGS-B");

  arma::vec lower(weights.n_elem, arma::fill::zeros);
  arma::vec upper(weights.n_elem, arma::fill::ones);
  arma::vec w2 = weights;

  opt.control.fnscale = 1;
  opt.set_lower(lower);
  opt.set_upper(upper);
  opt.minimize(optfu, weights);

  NumericVector opar = Rcpp::wrap(opt.par());
  opar.attr("dim") = R_NilValue;

  return Rcpp::List::create(_["value"] = opt.value(), _["par"] = opar);
}

*/



