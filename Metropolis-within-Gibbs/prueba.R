rm(list = ls())
library(Rcpp)
library(RcppArmadillo)
sourceCpp("Metropolis-within-Gibbs/mcmc_utils.cpp")

# simulated data
set.seed(123)

n  <- 10
T  <- 5
N  <- n * T
p  <- 2

coords <- matrix(runif(2*n), n, 2)
Dland  <- as.matrix(dist(coords))

X <- array(rnorm(N * p), dim = c(N, p))
Y <- rnorm(N)

Z <- cbind(1, coords)   # intercept + spatial covariates

set.seed(123)

n  <- 10
T  <- 5
N  <- n * T
p  <- 2

coords <- matrix(runif(2*n), n, 2)
Dland  <- as.matrix(dist(coords))

X <- array(rnorm(N * p), dim = c(N, p))
Y <- rnorm(N)

Z <- cbind(1, coords)   # intercept + spatial covariates


#algorithm
#initialization 
beta  <- rnorm(n)
gamma <- rnorm(ncol(Z))
xi    <- rexp(N)
sigma <- 1

sigma2_k <- 1
phi_k   <- 1
vars_k  <- 1
varphi_k <- 1


niter <- 1000
store_beta <- matrix(NA, niter, n)

for (it in 1:niter) {
  
  ## ---- 1. xi | ...
  R <- Y - X[,1] * beta[rep(1:n, each = T)]
  for (i in 1:N) {
    # aquí usarías rgig; lo dejamos aproximado
    xi[i] <- 1 / rgamma(1, 1, 1)
  }
  
  ## ---- 2. sigma | ...
  shape <- 1.5 * N + 1
  rate  <- sum(xi) + sum((R - xi)^2 / xi)
  sigma <- 1 / rgamma(1, shape, rate)
  
  ## ---- 3. gamma | beta
  Ek <- exp_cov(Dland, sigma2_k, phi_k)
  Vk <- solve(t(Z) %*% solve(Ek) %*% Z + diag(ncol(Z)))
  mk <- Vk %*% t(Z) %*% solve(Ek) %*% beta
  gamma <- as.vector(mvtnorm::rmvnorm(1, mk, Vk))
  
  ## ---- 4. beta | rest
  Vk_beta <- solve(diag(N) + solve(Ek))
  mk_beta <- Vk_beta %*% (X[,1] * Y + solve(Ek) %*% Z %*% gamma)
  beta <- as.vector(mvtnorm::rmvnorm(1, mk_beta, Vk_beta))
  
  ## ---- 5. MH: log(phi_k)
  logphi_prop <- log(phi_k) + rnorm(1, 0, 0.1)
  phi_prop <- exp(logphi_prop)
  
  Ek_prop <- exp_cov(Dland, sigma2_k, phi_prop)
  
  logpost_prop <- log_dmvnorm(beta, Z %*% gamma, Ek_prop) +
    dgamma(phi_prop, 2, 2, log = TRUE)
  
  logpost_curr <- log_dmvnorm(beta, Z %*% gamma, Ek) +
    dgamma(phi_k, 2, 2, log = TRUE)
  
  if (log(runif(1)) < logpost_prop - logpost_curr) {
    phi_k <- phi_prop
    Ek <- Ek_prop
  }
  
  store_beta[it, ] <- beta
}
