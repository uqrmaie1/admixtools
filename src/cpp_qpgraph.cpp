// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec cpp_opt_edge_lengths(const arma::mat& ppwts_2d, const arma::mat& ppinv, const arma::vec& f3_jest, Function qpsolve) {
  arma::mat pppp = ppwts_2d.t() * ppinv;
  arma::mat cc = pppp * ppwts_2d;
  int nc = cc.n_cols;
  for(int i = 0; i < nc; i++) {
    for(int j = i+1; j < nc; j++) {
      cc(i,j) = cc(j,i);
    }
  }
  for(int i = 0; i < nc; i++) {
    cc(i,i) += 0.0001;
  }
  arma::vec q1 = -(pppp * f3_jest);
  arma::mat CI(nc, nc, arma::fill::eye);
  arma::vec ci0 = arma::zeros<arma::vec>(nc);
  return -as<arma::vec>(qpsolve(cc, q1, -CI, ci0));
}


// [[Rcpp::export]]
void cpp_fill_pwts(arma::mat& pwts, const arma::vec& weights, const arma::mat& path_edge_table, const arma::mat& path_admixedge_table, int numpaths) {
  // puts weights onto pwts, using index matrix and vectors

  if(weights.n_elem == 0) return;

  arma::vec pathweights(numpaths, arma::fill::ones);
  arma::vec pwtsfill(pwts.n_elem, arma::fill::zeros);
  arma::vec pwtsfillind(pwts.n_elem, arma::fill::zeros);

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
}



// [[Rcpp::export]]
double cpp_optimweightsfun(arma::vec weights, List args) {

  arma::mat pwts = args[0];
  arma::mat ppinv = args[1];
  arma::vec f3_jest = args[2];
  arma::mat path_edge_table = args[3];
  arma::mat path_admixedge_table = args[4];
  int numpaths = args[5];
  // args[6] is cmb matrix with column combinations; not needed here
  Function qpsolve = args[7];

  double nr = pwts.n_rows;
  double nc = pwts.n_cols;
  cpp_fill_pwts(pwts, weights, path_edge_table, path_admixedge_table, numpaths);
  arma::mat ppwts_2d(nc*(nc+1)/2, nr);
  for(int i=0; i<nr; i++) {
    int c=0;
    for(int j=0; j<nc; j++) {
      for(int k=j; k<nc; k++) {
        ppwts_2d(c, i) = pwts(i, j) * pwts(i, k);
        c++;
      }
    }
  }
  arma::vec q2 = cpp_opt_edge_lengths(ppwts_2d, ppinv, f3_jest, qpsolve);
  arma::vec w2 = (ppwts_2d * q2) - f3_jest;
  arma::vec lik = w2.t() * ppinv * w2;
  return lik(0);
}


// Rcpp optimization using Roptim (cpp_lbfgsb / Optimfunctor in c++) is not actually faster than using R optim.

/*
// [[Rcpp::depends(roptim)]]
#include <roptim.h>
using namespace roptim;

class Optimfunctor : public Functor {

public:
  Optimfunctor(const arma::mat &pwts, const arma::mat &ppinv, const arma::vec &f3_jest, const arma::mat &path_edge_table,
               const arma::mat &path_admixedge_table, Function qpsolve, int numpaths) : qpsolve2(qpsolve) {
    thispwts = pwts;
    thisppinv = ppinv;
    thisf3_jest = f3_jest;
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

    arma::vec q2 = cpp_opt_edge_lengths(ppwts_2d, thisppinv, thisf3_jest, qpsolve2);
    arma::vec w2 = (ppwts_2d * q2) - thisf3_jest;
    arma::vec lik = w2.t() * thisppinv * w2;
    return lik(0);
  }

  private:
    arma::mat thispwts;
    arma::mat thisppinv;
    arma::vec thisf3_jest;
    arma::mat thispath_edge_table;
    arma::mat thispath_admixedge_table;
    Function qpsolve2;
    int thisnumpaths;
};


// [[Rcpp::export]]
List cpp_lbfgsb(arma::vec weights, const arma::mat &pwts, const arma::mat &ppinv,
                const arma::vec &f3_jest, const arma::mat &path_edge_table,
                const arma::mat &path_admixedge_table, int numpaths, Function qpsolve) {

  Optimfunctor optfu(pwts, ppinv, f3_jest, path_edge_table,
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






