#ifndef UTILS_MOD_H
#define UTILS_MOD_H

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
arma::vec rig(const int N, const arma::vec& mu, const double lambda, 
  const bool parallel = false, const int nThreads = 0);
arma::vec RandomMultiNormalC(const arma::mat& Q, const arma::vec& b);
arma:: mat inv_covariance_matrix()(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const arma::mat& dmat,
    const arma::vec& dvec,
    const arma::mat& dmatc
);

#endif