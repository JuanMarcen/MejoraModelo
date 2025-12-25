#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Log-densidad normal multivariante
// [[Rcpp::export]]
double log_dmvnorm(vec x, vec mu, mat Sigma) {
  int n = x.n_elem;
  mat iS = inv_sympd(Sigma);
  double sign, logdet;
  log_det(logdet, sign, Sigma);
  vec diff = x - mu;
  return -0.5 * (logdet + as_scalar(diff.t() * iS * diff) + n * log(2 * M_PI));
}

// Exponencial kernel
// [[Rcpp::export]]
mat exp_cov(mat D, double sigma2, double phi) {
  return sigma2 * exp(-phi * D);
}
