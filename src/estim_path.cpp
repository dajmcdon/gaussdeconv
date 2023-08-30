#include <Eigen/Sparse>
#include <RcppEigen.h>
#include "admm.h"
#include "utils.h"
#include "dptf.h"

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qr;

using namespace Rcpp;

// [[Rcpp::export]]
List estim_path(NumericVector y,
                NumericVector x,
                NumericVector w,
                int korder,
                NumericVector lambda,
                double lambdamax = -1,
                double lambdamin = -1,
                int nsol = 100,
                double rho = -1,
                int maxiter = 1e5,
                double tolerance = 1e-3,
                double lambda_min_ratio = 1e-4,
                int verbose = 0) {
  int n = y.size();

  // Placeholders for solutions
  NumericMatrix theta(n, nsol);
  NumericVector niter(nsol);
  NumericVector dof(nsol);

  // Build D matrices as needed
  Eigen::SparseMatrix<double> D;
  Eigen::SparseMatrix<double> Dk;
  Eigen::SparseMatrix<double> DkDk;
  D = get_D(korder, x);
  qr.compute(D.transpose());
  int m = n;
  if (korder > 0) {
    Dk = get_Dtil(korder, x);
    DkDk = Dk.transpose() * Dk;
    m = Dk.rows();
  }

  // Generate lambda sequence if necessary
  if (abs(lambda[nsol - 1]) < tolerance / 100 && lambdamax <= 0) {
    VectorXd b(n - korder);
    VectorXd wy = nvec_to_evec(w * y);
    b = qr.solve(wy);
    NumericVector bp = evec_to_nvec(b);
    lambdamax = max(abs(bp));
  }
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);

  // ADMM parameters
  double _rho;

  // ADMM variables
  NumericVector beta(n);
  NumericVector alpha(m);
  NumericVector u(m);
  int iters = 0;
  int nsols = nsol;

  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {
    if (verbose > 0) Rcout << ".";
    Rcpp::checkUserInterrupt();

    if (korder == 0) {
      beta = weight_dptf(y, lambda[i], w);
      niter[i] = 0;
    } else {
      _rho = (rho < 0) ? lambda[i] : rho;
      admm_gauss(maxiter, n, korder, y, x, w,
                 beta, alpha, u, _rho, lambda[i] / _rho,
                 DkDk, tolerance, iters);
      niter[i] = iters;
      maxiter -= iters + 1;
      if (maxiter < 0) nsols = i + 1;
    }

    // Store solution
    theta(_, i) = beta;
    if (korder == int(0)) {
      dof[i] = sum(abs(diff(beta)) > tolerance);
    } else {
      dof[i] = sum(abs(diff(alpha)) > tolerance);
    }

    // Verbose handlers
    if (verbose > 1) Rcout << niter(i);
    if (verbose > 2) Rcout << "(" << lambda(i) << ")";
    if (verbose > 0) Rcout << std::endl;
    if (maxiter < 0) break;
  }

  // Return
  List out = List::create(
    Named("ct") = theta(_, Range(0, nsols - 1)),
    Named("lambda") = lambda[Range(0, nsols - 1)],
    Named("korder") = korder,
    Named("dof") = dof[Range(0, nsols - 1)],
    Named("niter") = niter[Range(0, nsols - 1)]
  );

  return out;
}
