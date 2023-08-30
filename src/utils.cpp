#include <Eigen/Sparse>
#include <RcppEigen.h>
#include <dspline.h>
#include "utils.h"

using namespace Rcpp;
using Eigen::VectorXd;

NumericVector evec_to_nvec(VectorXd evec) {
  NumericVector nvec(wrap(evec));
  return nvec;
}

VectorXd nvec_to_evec(NumericVector nvec) {
  VectorXd evec = as<Eigen::Map<VectorXd> >(nvec);
  return evec;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_Dtil(int k, NumericVector xd) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, false, Rcpp::seq(0, n - k - 1), true);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_D(int k, NumericVector xd) {
  k++;
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, true, Rcpp::seq(0, n - k - 1), true);
}

void create_lambda(NumericVector& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol) {
  if (all(lambda == 0).is_false()) {
    lambdamin = min(lambda);
    lambdamax = max(lambda);
    nsol = lambda.size();
  } else {
    lambdamin = (lambdamin < 0) ? lambda_min_ratio * lambdamax : lambdamin;
    double lmpad = 1e-20;
    double ns = static_cast<double>(nsol);
    double p;
    if (lambdamin > lmpad) {
      p = pow(lambdamin / lambdamax, 1 / (ns - 1));
      lambda(0) = lambdamax;
      for (int i = 1; i < nsol; i++) lambda[i] = lambda[i - 1] * p;
    } else {
      p = pow(lmpad / lambdamax, 1 / (ns - 2));
      lambda(0) = lambdamax;
      for (int i = 1; i < nsol - 1; i++) lambda[i] = lambda[i - 1] * p;
      lambda(nsol - 1) = lambdamin;
    }
  }
}

// [[Rcpp::export]]
NumericVector create_lambda_test(NumericVector lambda,
                                 double lambdamin,
                                 double lambdamax,
                                 double lambda_min_ratio,
                                 int nsol) {
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);
  return lambda;
}

// [[Rcpp::export]]
NumericVector doDv(NumericVector v, int k, NumericVector xd) {
  return dspline::rcpp_d_mat_mult(v, k, xd, false, false);
}

// [[Rcpp::export]]
NumericVector doDtv(NumericVector v, int k, NumericVector xd) {
  return dspline::rcpp_d_mat_mult(v, k, xd, false, true);
}
