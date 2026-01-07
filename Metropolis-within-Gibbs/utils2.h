#ifndef UTILS2_H
#define UTILS2_H

#include <RcppArmadillo.h>
#include <random>
#include <chrono>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

// Aux functions
void reportProgress(int iter, int nBurnin, int nSims, int nReport, int save_idx = 0,
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now());
arma::vec rig(const int N, const arma::vec& mu, const double lambda);
arma::vec RandomMultiNormalC(const arma::mat& Q, const arma::vec& b);
arma::mat inv_covariance_matrix(
    const double precision, 
    const double decay, 
    const double varsigma, 
    const double varphi, 
    const double cprec,
    const arma::mat& dmat,
    const arma::vec& dvec,
    const arma::mat& dmatc
);
arma::mat inv_conv_covariance_matrix(
    const double precision, 
    const double decay, 
    const double varsigma, 
    const double varphi, 
    const double cprec, 
    const arma::mat& dmat,
    const arma::mat& dmatcoast,
    const arma::mat& dr,
    const double lencoast
);

#endif