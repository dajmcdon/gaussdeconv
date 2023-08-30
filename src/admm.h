#ifndef __ADMM_H
#define __ADMM_H

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
                int& iter);

#endif
