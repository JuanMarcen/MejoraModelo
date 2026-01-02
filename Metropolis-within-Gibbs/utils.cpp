#include "utils.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

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

// @description Computes the distance matrix.
// @param coords \eqn{n \times 2} coordinates matrix.
// [[Rcpp::export]]
arma::mat dist1(const arma::mat& coords) {
  double dx, dy;
  int n = coords.n_rows;
  arma::mat D(n, n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      dx = coords(i,0) - coords(j,0);
      dy = coords(i,1) - coords(j,1);
      D(i,j) = std::sqrt(dx*dx + dy*dy);
      D(j,i) = D(i,j);
    }
  }
  return D;
}

// @description Computes the squared distance matrix.
// @param coords \eqn{n \times m} coordinates matrix.
// [[Rcpp::export]]
arma::mat dist2(const arma::mat& coords) {
  double dx, dy;
  int n = coords.n_rows;
  arma::mat D(n, n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      dx = coords(i,0) - coords(j,0);
      dy = coords(i,1) - coords(j,1);
      D(i,j) = dx*dx + dy*dy;
      D(j,i) = D(i,j);
    }
  }
  return D;
}

// @description Computes the distance between two sets of coordinates.
// @param A,B \eqn{n_{A} \times 2} and \eqn{n_{B} \times 2} coordinate matrices.
// [[Rcpp::export]]
arma::mat dist_mat(const arma::mat& A, const arma::mat& B) {
  int n1 = A.n_rows;
  int n2 = B.n_rows;
  arma::mat D(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      D(i,j) = arma::norm(A.row(i) - B.row(j));
    }
  }
  return D;
}

// [[Rcpp::export]]
double dtnorm(double x, double mu, double sigma, double a, double b) { // log
  return R::dnorm(x, mu, sigma, 1) - std::log(R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0));
}

// [[Rcpp::export]]
double rtnorm(double mu, double sigma, double a, double b) {
  double alpha = (a - mu) / sigma;
  double beta  = (b - mu) / sigma;
  double pAlpha = R::pnorm(alpha, 0, 1, 1, 0);
  double pBeta  = R::pnorm(beta, 0, 1, 1, 0);
  return R::qnorm(pAlpha + R::runif(0, 1) * (pBeta - pAlpha), 0, 1, 1, 0) * sigma + mu;
}

// [[Rcpp::export]]
double psi1(double x, double alpha, double lambda) {
  Rcpp::NumericVector coshx = { x };
  coshx = Rcpp::cosh(coshx);
  return - alpha * (coshx(0) - 1) - lambda * (exp(x) - x - 1);
}
// [[Rcpp::export]]
double psi2(double x, double alpha, double lambda) {
  Rcpp::NumericVector sinhx = { x };
  sinhx = Rcpp::sinh(sinhx);
  return - alpha * sinhx(0) - lambda * (std::exp(x) - 1);
}
// [[Rcpp::export]]
arma::vec rgig(
    const int N, 
    const double a, 
    const arma::vec b, 
    const double nu) {
  
  arma::vec X(N);
  static const double cosh1 = 1.543081;
  
  double lambda = nu;
  arma::vec omega = arma::sqrt(a * b);
  
  double alpha;
  double t, s;
  double tp, sp, q;
  double eta, zeta, theta, xi;
  double p, r;
  double chi;
  double U, V, W;
  
  double aux;
  
  Rcpp::NumericVector v(2);
  
  for (int i = 0; i < N; ++i) {
    alpha = std::sqrt(omega(i) * omega(i) + lambda * lambda) - lambda;
    aux = - psi1(1, alpha, lambda);
    if (aux > 2) {
      t = std::sqrt(2 / (alpha + lambda));
    } else if (aux < 0.5) {
      t = std::log(4 / (alpha + 2 * lambda));
    } else {
      t = 1;
    }
    aux = - psi1(-1, alpha, lambda);
    if (aux > 2) {
      s = std::sqrt(4 / (alpha * cosh1 + lambda));
    } else if (aux < 0.5) {
      aux = 1 / alpha;
      v = { 1 / lambda, log(1 + aux + std::sqrt(aux * aux + 2 * aux)) };
      s = Rcpp::min(v);
    } else {
      s = 1;
    }
    
    eta   = - psi1( t, alpha, lambda);
    zeta  = - psi2( t, alpha, lambda);
    theta = - psi1(-s, alpha, lambda);
    xi    =   psi2(-s, alpha, lambda);
    
    p = 1 / xi;
    r = 1 / zeta;
    
    tp = t - r * eta;
    sp = s - p * theta;
    q  = tp + sp;
    
    do {
      U = R::runif(0, 1);
      V = R::runif(0, 1);
      W = R::runif(0, 1);
      if (U < (q / (p + q + r))) {
        X(i) = - sp + q * V;
      } else if (U < ((q + r) / (p + q + r))) {
        X(i) = tp + r * std::log(1 / V);
      } else {
        X(i) = - sp - p * std::log(1 / V);
      }
      if (X(i) > tp) {
        chi = exp(- eta - zeta * (X(i) - t));
      } else if (X(i) < (- sp)) {
        chi = exp(- theta + xi * (X(i) + s));
      } else {
        chi = 1;
      }
    } while ((W * chi) > std::exp(psi1(X(i), alpha, lambda)));
    
    aux = lambda / omega(i);
    X(i) = (aux + std::sqrt(1 + aux * aux)) * std::exp(X(i));
  }
  
  return X % arma::sqrt(b / a);
}

// [[Rcpp::export]]
arma::vec rig(
    const int N,
    const arma::vec& mu,
    const double lambda,
    const bool parallel,
    int nThreads) {
  
  arma::vec X(N);
  arma::vec V = arma::square(arma::randn(N)); // chi^2(1)
  arma::vec U = arma::randu(N);               // uniform
  
  if (!parallel) {
    
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
  
#ifdef _OPENMP
  if (nThreads <= 0)
    nThreads = omp_get_max_threads();
  omp_set_num_threads(nThreads);
#endif
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < N; ++i) {
    double mu_i = mu(i);
    double W = mu_i * V(i);
    double C = mu_i / (2.0 * lambda);
    double X_i = mu_i + C * (W - std::sqrt(W * (4.0 * lambda + W)));
    double P1 = mu_i / (mu_i + X_i);
    if (U(i) > P1)
      X_i = mu_i * mu_i / X_i;
    X(i) = X_i;
  }

  return X;
}

// [[Rcpp::export]]
arma::vec ralRcpp(
    const arma::vec sigma,
    const double tau) {
  
  const int N = sigma.n_elem;
  const double aux = (tau * (1 - tau));
  const double theta = (1 - 2 * tau) / aux;
  const double w2 = 2 / aux;
  
  arma::vec U = - log(arma::randu(N));
  U = (theta * U + sqrt(w2 * U) % arma::randn(N)) % sigma;
  
  return U;
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
arma::mat ffbs(
    const arma::mat& Y,       // n x T : observations
    const double rho,         // AR parameter
    const arma::mat& Sigma_w, // state conditional covariance
    const Rcpp::List& D_list  // list of length T with vec of diag matrices
) {
  const int n = Y.n_rows;
  const int T = Y.n_cols;
  
  // Output: states w_{1:T}
  arma::mat W(n, T);
  
  // Storage for filtering
  std::vector<arma::vec> m_tt(T+1);   // m_{t|t}
  std::vector<arma::mat> P_tt(T+1);   // P_{t|t}
  std::vector<arma::vec> m_ttm1(T+1); // m_{t|t-1}
  std::vector<arma::mat> P_ttm1(T+1); // P_{t|t-1}
  
  // ---- INITIALIZATION (t = 0) ----
  m_tt[0] = arma::zeros(n);
  P_tt[0] = Sigma_w / (1.0 - rho*rho);   // stationary marginal covariance
  
  // aux
  arma::vec v(n), D(n), mu(n);
  arma::mat S(n,n), K(n,n);
  arma::mat J(n,n), Sigma(n,n);
  
  // ---- FORWARD FILTERING ----
  for (int t = 1; t <= T; t++) {
    
    // Predict
    m_ttm1[t] = rho * m_tt[t-1];
    P_ttm1[t] = rho*rho * P_tt[t-1] + Sigma_w;
    
    // Innovation
    v = Y.col(t-1) - m_ttm1[t];
    
    // Observation noise
    D = Rcpp::as<arma::vec>(D_list[t-1]);
    S = P_ttm1[t]; 
    S.diag() += D;
    
    // Compute Kalman gain: K = P_ttm1 * S^{-1}
    // K = P_ttm1[t] * arma::inv_sympd(S);
    K = arma::solve(S, P_ttm1[t], arma::solve_opts::fast).t();
    
    // Update
    m_tt[t] = m_ttm1[t] + K * v;
    P_tt[t] = P_ttm1[t] - K * S * K.t();
  }
  
  // ---- BACKWARD SAMPLING ----
  // Sample w_T ~ N(m_{T|T}, P_{T|T})
  W.col(T-1) = m_tt[T] + arma::chol(P_tt[T], "lower") * arma::randn(n);
  
  // Backwards recursion
  for (int t = T-1; t >= 1; t--) {
    
    // J_t = rho P_{t|t} * P_{t+1|t}^{-1}
    // J = rho * P_tt[t] * arma::inv_sympd(P_ttm1[t+1]);
    J = arma::solve(P_ttm1[t+1], rho * P_tt[t], arma::solve_opts::fast).t();
    
    // Conditional mean/covariance
    mu = m_tt[t] + J * (W.col(t) - m_ttm1[t+1]);
    Sigma = P_tt[t] - J * P_ttm1[t+1] * J.t();
    
    // Sample
    W.col(t-1) = mu + arma::chol(Sigma, "lower") * arma::randn(n);
  }
  
  return W;
}


//// [[Rcpp::export]]
//arma::vec priorSpatialGP(
//  arma::vec v, 
//  double mean, double prec, double decay,
//  double Rsum, arma::mat Rinv
//  ) {
//  
//  const double na = 0;
//  const double nb = 0.000001;
//  
//  double ZtRZ = arma::as_scalar(Z.t() * Rinv * Z);
//  
//  // mean
//  V    = 1 / (prec * Rsum + nb); 
//  mean = V * (prec * arma::as_scalar(onen.t() * Rinv * v) + nb * na);
//  mean = R::rnorm(mean, sqrt(V));
//  
//  // prec
//  Z = v - mean;
//  prec = R::rgamma(n / 2 + ga, 1 / (ZtRZ / 2 + gb));
//  
//  // decay
//  decay_aux   = rtnorm(decay, sd, Umin, Umax);
//  Rinv_aux    = arma::inv_sympd(exp(- decay_aux * dist));
//  Rlogdet_aux = arma::log_det_sympd(Rinv_aux);
//  ZtRZ_aux    = arma::as_scalar(Z.t() * Rinv_aux * Z);  
//  A = (Rlogdet_aux - prec * ZtRZ_aux) / 2 - 
//    (Rlogdet - prec * ZtRZ) / 2 + 
//    dtnorm(decay, decay_aux, sd, Umin, Umax) - 
//   dtnorm(decay_aux, decay, sd, Umin, Umax);
//  if (log(R::runif(0, 1)) <= A) {
//    decay = decay_aux;
//    Rinv = Rinv_aux;
//    Rlogdet = Rlogdet_aux;
//    Rsum = arma::accu(Rinv);
//    ZtRZ = ZtRZ_aux;
//  }
//  
//  return { mean, prec, decay };
//}
