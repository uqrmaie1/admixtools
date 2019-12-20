//  Copyright (C) 2018 Yi Pan <ypan1988@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  https://www.R-project.org/Licenses/

#ifndef ROPTIM_H_
#define ROPTIM_H_

#include <cassert>
#include <cmath>
#include <cstddef>

#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <RcppArmadillo.h>

#include "applic.h"
#include "functor.h"
#include "samin.h"

namespace roptim {

template <typename Derived>
class Roptim {
 public:
  std::string method_;
  arma::vec lower_, upper_;
  bool hessian_flag_ = false;
  arma::mat hessian_;

  arma::vec lower() const { return lower_; }
  arma::vec upper() const { return upper_; }

  arma::vec par() const { return par_; }
  double value() const { return val_; }
  int fncount() const { return fncount_; }
  int grcount() const { return grcount_; }
  int convergence() const { return fail_; }
  std::string message() const { return message_; }
  arma::mat hessian() const { return hessian_; }

  void print() const;

 private:
  arma::vec par_;
  double val_ = 0.0;
  int fncount_ = 0;
  int grcount_ = 0;
  int fail_ = 0;
  std::string message_ = "NULL";

 public:
  struct RoptimControl {
    int trace = 0;
    double fnscale = 1.0;
    arma::vec parscale;
    arma::vec ndeps;
    std::size_t maxit = 100;
    double abstol = R_NegInf;
    double reltol = sqrt(2.220446e-16);
    double alpha = 1.0;
    double beta = 0.5;
    double gamma = 2.0;
    int REPORT = 10;
    bool warn_1d_NelderMead = true;
    int type = 1;
    int lmm = 5;
    double factr = 1e7;
    double pgtol = 0.0;
    double temp = 10.0;
    int tmax = 10;
  } control;

  Roptim(const std::string method = "Nelder-Mead") : method_(method) {
    if (method_ != "Nelder-Mead" && method_ != "BFGS" && method_ != "CG" &&
        method_ != "L-BFGS-B" && method_ != "SANN")
      Rcpp::stop("Roptim::Roptim(): unknown 'method'");

    // Sets default value for maxit & REPORT (which depend on method)
    if (method_ == "Nelder-Mead") {
      control.maxit = 500;
    } else if (method_ == "SANN") {
      control.maxit = 10000;
      control.REPORT = 100;
    }
  }

  void set_method(const std::string &method) {
    if (method != "Nelder-Mead" && method != "BFGS" && method != "CG" &&
        method != "L-BFGS-B" && method != "SANN")
      Rcpp::stop("Roptim::set_method(): unknown 'method'");
    else
      method_ = method;

    // Sets default value for maxit & REPORT (which depend on method)
    if (method_ == "Nelder-Mead") {
      control.maxit = 500;
      control.REPORT = 10;
    } else if (method_ == "SANN") {
      control.maxit = 10000;
      control.REPORT = 100;
    } else {
      control.maxit = 100;
      control.REPORT = 10;
    }
  }

  void set_lower(const arma::vec &lower) {
    if (method_ != "L-BFGS-B")
      Rcpp::warning(
          "Roptim::set_lower(): bounds can only be used with method L-BFGS-B");
    method_ = "L-BFGS-B";
    lower_ = lower;
  }

  void set_upper(const arma::vec &upper) {
    if (method_ != "L-BFGS-B")
      Rcpp::warning(
          "Roptim::set_upper(): bounds can only be used with method L-BFGS-B");
    method_ = "L-BFGS-B";
    upper_ = upper;
  }

  void set_hessian(bool flag) { hessian_flag_ = flag; }

  void minimize(Derived &func, arma::vec &par);
};

template <typename Derived>
inline void Roptim<Derived>::print() const {
  par_.t().print(".par()");
  Rcpp::Rcout << "\n.value()\n" << val_ << std::endl;
  Rcpp::Rcout << "\n.fncount()\n" << fncount_ << std::endl;

  if (method_ == "Nelder-Mead" || method_ == "SANN")
    Rcpp::Rcout << "\n.grcount()\nNA" << std::endl;
  else
    Rcpp::Rcout << "\n.grcount()\n" << grcount_ << std::endl;

  Rcpp::Rcout << "\n.convergence()\n" << fail_ << std::endl;
  Rcpp::Rcout << "\n.message()\n" << message_ << std::endl;
  if (hessian_flag_) hessian_.print("\n.hessian()");
  Rcpp::Rcout << std::endl;
}

template <typename Derived>
inline void Roptim<Derived>::minimize(Derived &func, arma::vec &par) {
  int debug = 0;

  // PART 1: optim()

  // Checks if lower and upper bounds is used
  if ((!lower_.is_empty() || !upper_.is_empty()) && method_ != "L-BFGS-B") {
    Rcpp::warning("bounds can only be used with method L-BFGS-B");
    method_ = "L-BFGS-B";
  }

  // Sets the parameter size
  std::size_t npar = par.size();

  // Sets default value for parscale & ndeps (which depend on npar)
  if (control.parscale.is_empty())
    control.parscale = arma::ones<arma::vec>(npar);
  if (control.ndeps.is_empty())
    control.ndeps = arma::ones<arma::vec>(npar) * 1e-3;

  // Checks control variable trace
  if (control.trace < 0)
    Rcpp::warning("read the documentation for 'trace' more carefully");
  else if (method_ == "SANN" && control.trace && control.REPORT == 0)
    Rcpp::stop("'trace != 0' needs 'REPORT >= 1'");

  // Note that "method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol'
  // and 'abstol'". There is no simple way to detect whether users set new
  // values for 'reltol' and 'abstol'.

  // Gives warning of 1-dim optimization by Nelder-Mead
  if (npar == 1 && method_ == "Nelder-Mead" && control.warn_1d_NelderMead)
    Rcpp::warning("one-dimensional optimization by Nelder-Mead is unreliable");

  // Sets default value for lower_
  if (method_ == "L-BFGS-B" && lower_.is_empty()) {
    lower_ = arma::zeros<arma::vec>(npar);
    lower_.for_each([](arma::mat::elem_type &val) { val = R_NegInf; });
  }
  // Sets default value for upper_
  if (method_ == "L-BFGS-B" && upper_.is_empty()) {
    upper_ = arma::zeros<arma::vec>(npar);
    upper_.for_each([](arma::mat::elem_type &val) { val = R_PosInf; });
  }

  // PART 2: C_optim()

  func.os.usebounds_ = 0;
  func.os.fnscale_ = control.fnscale;
  func.os.parscale_ = control.parscale;

  if (control.ndeps.size() != npar)
    Rcpp::stop("'ndeps' is of the wrong length");
  else
    func.os.ndeps_ = control.ndeps;

  arma::vec dpar = arma::zeros<arma::vec>(npar);
  arma::vec opar = arma::zeros<arma::vec>(npar);

  dpar = par / control.parscale;

  if (method_ == "Nelder-Mead") {
    if (debug) std::cout << "Nelder-Mead:" << std::endl;

    nmmin(npar, dpar.memptr(), opar.memptr(), &val_, fminfn, &fail_,
          control.abstol, control.reltol, &func, control.alpha, control.beta,
          control.gamma, control.trace, &fncount_, control.maxit);

    par = opar % control.parscale;
    grcount_ = 0;

  } else if (method_ == "SANN") {
    if (debug) std::cout << "SANN:" << std::endl;

    int trace = control.trace;
    if (trace) trace = control.REPORT;

    if (control.tmax == NA_INTEGER || control.tmax < 1)
      Rcpp::stop("'tmax' is not a positive integer");

    internal::samin(npar, dpar.memptr(), &val_, fminfn, control.maxit,
                    control.tmax, control.temp, trace, &func);
    par = dpar % control.parscale;
    fncount_ = npar > 0 ? control.maxit : 1;
    grcount_ = 0;

  } else if (method_ == "BFGS") {
    if (debug) std::cout << "BFGS:" << std::endl;

    arma::ivec mask = arma::ones<arma::ivec>(npar);
    vmmin(npar, dpar.memptr(), &val_, fminfn, fmingr, control.maxit,
          control.trace, mask.memptr(), control.abstol, control.reltol,
          control.REPORT, &func, &fncount_, &grcount_, &fail_);

    par = dpar % control.parscale;
  } else if (method_ == "CG") {
    cgmin(npar, dpar.memptr(), opar.memptr(), &val_, fminfn, fmingr, &fail_,
          control.abstol, control.reltol, &func, control.type, control.trace,
          &fncount_, &grcount_, control.maxit);

    par = opar % control.parscale;
  } else if (method_ == "L-BFGS-B") {
    arma::vec lower(npar);
    arma::vec upper(npar);
    arma::ivec nbd = arma::zeros<arma::ivec>(npar);
    char msg[60];

    for (std::size_t i = 0; i != npar; ++i) {
      lower(i) = lower_(i) / func.os.parscale_(i);
      upper(i) = upper_(i) / func.os.parscale_(i);
      if (!std::isfinite(lower(i))) {
        if (!std::isfinite(upper(i)))
          nbd(i) = 0;
        else
          nbd(i) = 3;
      } else {
        if (!std::isfinite(upper(i)))
          nbd(i) = 1;
        else
          nbd(i) = 2;
      }
    }

    func.os.usebounds_ = 1;
    func.os.lower_ = lower;
    func.os.upper_ = upper;

    lbfgsb(npar, control.lmm, dpar.memptr(), lower.memptr(), upper.memptr(),
           nbd.memptr(), &val_, fminfn, fmingr, &fail_, &func, control.factr,
           control.pgtol, &fncount_, &grcount_, control.maxit, msg,
           control.trace, control.REPORT);

    par = dpar % control.parscale;
    message_ = msg;
  } else
    Rcpp::stop("Roptim::minimize(): unknown 'method'");

  par_ = par;
  val_ *= func.os.fnscale_;

  if (hessian_flag_) func.ApproximateHessian(par_, hessian_);
}

}  // namespace roptim

#endif  // ROPTIM_H_
