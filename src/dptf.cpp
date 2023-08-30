#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
extern "C" {
#include "tf_dp.h"
}
#include "dptf.h"

// [[Rcpp::export]]
Rcpp::NumericVector dptf(Rcpp::NumericVector y, double lam) {
  int n = y.size();
  Rcpp::NumericVector beta(n);
  tf_dp(n, y.begin(), lam, beta.begin());  // get pointers to the vectors
  return beta;
}

/**
 * Filter mean (i.e., exponential of natural parameter) trend for 0th-order
 * weighted trend filtering
 * @param y signal length
 * @param lam hyperparameter
 * @param w signal weights
 * @return filtered mean (exponential of natural parameter) trend
 */
// [[Rcpp::export]]
Rcpp::NumericVector weight_dptf(Rcpp::NumericVector y,
                                double lam,
                                Rcpp::NumericVector w) {
  int n = y.size();
  Rcpp::NumericVector beta(n);
  lam = n * lam;
  tf_dp_weight(n, y.begin(), w.begin(), lam, beta.begin());
  return beta;
}
