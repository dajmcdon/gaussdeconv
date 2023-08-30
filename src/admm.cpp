#include <RcppEigen.h>
#include <Eigen/Sparse>
#include "utils.h"
#include "dptf.h"
#include "admm.h"

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qradmm;

using namespace Rcpp;

/**
 * ADMM for Gaussian trend filtering
 * @param M maximum iteration of the algos
 * @param n signal length
 * @param korder degree of trend filtering
 * @param y observed signals
 * @param x signal locations
 * @param w signal weights
 * @param theta primal variable of length `n`
 * @param z auxiliary variable of length `n-korder`
 * @param u dual variable of length `n-korder`
 * @param rho Lagrangian parameter of ADMM
 * @param lam_z hyperparameter of the auxiliary step of ADMM
 * @param DD D^T * D
 * @param tol tolerance of stopping criteria
 * @param iter interation index
 */
void admm_gauss(int M,
                int n,
                int korder,
                Rcpp::NumericVector const& y,
                Rcpp::NumericVector const& x,
                Rcpp::NumericVector const& w,
                Rcpp::NumericVector& theta,
                Rcpp::NumericVector& z,
                Rcpp::NumericVector& u,
                double rho,
                double lam_z,
                Eigen::SparseMatrix<double> const& DD,
                double tol,
                int& iter) {
  double r_norm = 0.0;
  double s_norm = 0.0;
  NumericVector z_old = clone(z);
  NumericVector tmp_n(n);
  NumericVector tmp_m(z.size());
  NumericVector Dth(z.size());
  VectorXd tmp_theta(n);
  SparseMatrix<double> cDD = DD * n * rho;  // a copy that doesn't change
  VectorXd eW = nvec_to_evec(w);
  for (int i = 0; i < n; i++) {
    cDD.diagonal()(i) += eW(i);
  }
  qradmm.compute(cDD);

  for (iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0) Rcpp::checkUserInterrupt();
    // solve for primal variable - theta:
    tmp_n = doDtv(z - u, korder, x) * n * rho;
    tmp_n += w * y;
    tmp_theta = nvec_to_evec(tmp_n);
    tmp_theta = qradmm.solve(tmp_theta);
    theta = evec_to_nvec(tmp_theta);
    // solve for alternating variable - z:
    Dth = doDv(theta, korder, x);
    tmp_m = Dth + u;
    z = dptf(tmp_m, lam_z);
    // update dual variable - u:
    u += Dth - z;

    // primal residuals:
    r_norm = sqrt(mean(pow(Dth - z, 2)));
    // dual residuals:
    tmp_n = doDtv(z - z_old, korder, x);
    s_norm = rho * sqrt(mean(pow(tmp_n, 2)));
    // stopping criteria check:
    if (r_norm < tol && s_norm < tol) break;

    // auxiliary variables update:
    z_old = z;
  }
}
