#include "utils.mod.h"

void reportProgress(
    int iter,
    int nBurnin,
    int nSims,
    int nReport,
    int save_idx,
    std::chrono::steady_clock::time_point start_time) {
  
  int iter_from_start = iter - (1 - nBurnin);
  int total_iters = nBurnin + nSims;
  
  if (nReport > 0 && (iter_from_start % nReport) == 0) {
    
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count() / 1000.0;
    
    double pct = 100.0 * double(iter_from_start) / double(std::max(1, total_iters));
    
    double est_total = elapsed / std::max(0.001, pct / 100.0); // total estimado
    double remaining = std::max(0.0, est_total - elapsed);
    
    if (iter <= 0) {
      int burned = iter + nBurnin;
      Rcpp::Rcout << "Burn-in: iteration " << burned << " / " << nBurnin
                  << " (" << std::floor(pct) << "%)"
                  << " | elapsed: " << std::round(elapsed) << "s"
                  << " | ETA: " << std::round(remaining) << "s"
                  << "\n";
    } else {
      Rcpp::Rcout << "Sampling: iteration " << iter << " / " << nSims
                  << ", saved: " << save_idx
                  << " (" << std::floor(pct) << "%)"
                  << " | elapsed: " << std::round(elapsed) << "s"
                  << " | ETA: " << std::round(remaining) << "s"
                  << "\n";
    }
  }
}


arma::vec RandomMultiNormalC(
    const arma::mat& Q, 
    const arma::vec& b) {
  
  int n = Q.n_rows;
  arma::vec x(n);
  
  arma::mat L = arma::chol(Q, "lower");
  arma::vec y = arma::solve(arma::trimatl(L), b);
  arma::vec mu = arma::solve(arma::trimatu(L.t()), y);
  arma::vec z = arma::randn(n);
  arma::vec u = arma::solve(arma::trimatu(L.t()), z);
  x = mu + u;
  
  return x;
}

arma::vec rig(
    const int N,
    const arma::vec& mu,
    const double lambda) {
  
  arma::vec X(N);
  arma::vec V = arma::square(arma::randn(N)); // chi^2(1)
  arma::vec U = arma::randu(N);               // uniform
    
  arma::vec C = mu / (2.0 * lambda);
  double mu_i, X_i, W, P1;
  
  for (int i = 0; i < N; ++i) {
    mu_i = mu(i);
    W = mu_i * V(i);
    X_i = mu_i + C(i) * (W - std::sqrt(W * (4.0 * lambda + W)));
    P1 = mu_i / (mu_i + X_i);
    if (U(i) > P1)
      X_i = mu_i * mu_i / X_i;
    X(i) = X_i;
  }
    
    return X;
}



arma:: mat inv_covariance_matrix()(
  const double precision, //hp(0, m)
  const double decay, //hp(1, m)
  const double varsigma, //hp(2, m)
  const double varphi, //hp(3, m)
  const arma::mat& dmat,
  const arma::vec& dvec,
  const arma::mat& dmatc
){
  // new definition of the covariance matrices
  arma::vec expdc = arma::exp(- varsigma* dvec); //nx1
  arma::mat Mcoast = arma::exp(- varphi * dmatc); //nxn
  // column multiplication
  Mcoast.each_row() %= expdc.t();
  // row multiplication
  Mcoast.each_col() %= expdc;
  
  arma::mat R = arma::inv_sympd(1.0 / precision * exp(- decay * dmat) + Mcoast);
  
  return R;
}