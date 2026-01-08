# documento de pruebas para mis funciones modificadas
rm(list = ls())
tau <- 0.5

library(qs)
df <- qread('C:/Users/PC/Desktop/Juan/df_jun_ag.qs')
stations <- readRDS('stations.rds')

Y <- df$Y
X <- matrix(0, nrow = length(Y), ncol = 0)
V <- cbind(1, df$s.1, df$c.1, df$g300, df$g500, df$g700)

# X_alpha
aux.vars <- c('elev', 'dist')
X_alpha <- list()
for (i in 1:ncol(V)){
  X_alpha[[i]] <- cbind(1, unique(df$elev), unique(df$dist))
}

dist <- readRDS('maps coast/dist.matrix.rds')

#priors
M <- matrix(0, nrow = 0, ncol = 1)
P <- matrix(0, nrow = 0, ncol = 0)

M_beta_alpha <- list()
P_beta_alpha <- list()
for (i in 1:length(X_alpha)){
  M_beta_alpha[[i]] <- rep(0, 3)
  P_beta_alpha[[i]] <- 0.0001 * diag(3)
}

da <- 38
db <- 7400
ga <- 2
gb <- 1
ra <- 83
rb <- 24600
na <- 0.1
nb <- 0.1

# initial
beta <- matrix(0, nrow = 0, ncol = 1)
alpha <- matrix(0.1, nrow = nrow(stations), ncol = ncol(V))
prec <- 1

hp <- matrix(0, nrow = 5, ncol = ncol(V))
hp[1, ] <- 1 #precision
hp[2, ] <- 3/600 #decay
hp[3, ] <- 1 #varsigma
hp[4, ] <- 3/900 # varphi
hp[5, ] <- 1

beta_alpha <- list()
for (i in 1:length(X_alpha)){
  beta_alpha[[i]] <- rep(0.1, times = ncol(X_alpha[[i]]))
}

# more constants 
N <- nrow(df)
n <- nrow(stations)
p <- 0
r <- ncol(V)

p_alpha <- unlist(lapply(X_alpha, ncol))
s <- rep(0:39, each = 5888)

nSims <- 100000
nThin <- 100
nBurnin <- 100000
nReport <- 1000

#more distances
dist_coast <- readRDS('maps coast/dist.vec.rds')
dist_coast_points <- readRDS('maps coast/dist.coast.points2.rds')

dmatcoast_conv <- readRDS('maps coast/phimat.rds')
drmat_conv <- readRDS('maps coast/dr.rds')
lencoast_conv <- drmat_conv[1, 200]

Rcpp::sourceCpp("Metropolis-within-Gibbs/mcmc3.cpp")

#try my mcmc function function
basura <- inv_covariance_matrix(hp[1,1], hp[2,1], hp[3,1], hp[4,1], hp[5,1], dist, dist_coast, dist_coast_points)
class(basura)
dim(basura)

basura <- inv_conv_covariance_matrix(hp[1,1], hp[2,1], hp[3,1], hp[4,1], hp[5,1], 
                                     dist, dmatcoast_conv, drmat_conv, lencoast_conv)
class(basura)
dim(basura)
repeat {
  basura <- try(spQuantileRcpp(
    tau = tau,
    Y = Y,
    X = X,
    V = V,
    X_alpha = X_alpha,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha,
    P_beta_alpha = P_beta_alpha,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha,
    prec = prec,
    hp = hp,
    beta_alpha = beta_alpha,
    N = N,
    n = n,
    p = p,
    r = r,
    p_alpha = p_alpha,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(basura, "try-error")) {
    return(basura)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}

basura <-spQuantileRcpp(
  tau = tau,
  Y = Y,
  X = X,
  V = V,
  X_alpha = X_alpha,
  dist = dist,
  dist_coast = dist_coast,
  dist_coast_point = dist_coast_points,
  dmatcoast_conv = dmatcoast_conv,
  drmat_conv = drmat_conv,
  lencoast_conv = lencoast_conv,
  M = M,
  P = P,
  M_beta_alpha = M_beta_alpha,
  P_beta_alpha = P_beta_alpha,
  da = da,
  db = db,
  ga = ga,
  gb = gb,
  ra = ra,
  rb = rb,
  na = na, 
  nb = nb,
  beta = beta,
  alpha = alpha,
  prec = prec,
  hp = hp,
  beta_alpha = beta_alpha,
  N = N,
  n = n,
  p = p,
  r = r,
  p_alpha = p_alpha,
  nSims = nSims,
  nThin = nThin,
  nBurnin = nBurnin,
  nReport = nReport,
  s = s)

# traceplots pf GP

names <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700')
for (j in 0:5){
  filename <- paste0('Metropolis-within-Gibbs/trplots/trplots_', names[j +1], '.pdf')
  pdf(filename, width = 12, height = 15)
  par(mfrow = c(8, 5))
  for (i in 1:40){
    plot(basura$process[, i + j * 48], type = 'l', main = paste0(names[j + 1], '(', stations$NAME2[i], ')'))
  }
  dev.off()
}

# tarceplots of parameters in GP

for (i in 0:5){
  filename <- paste0('Metropolis-within-Gibbs/trplots/hp_trplots_', names[i +1], '.pdf')
  pdf(filename, width = 7, height = 10)
  par(mfrow = c(4,2))
  plot(basura$process[, 41 + i * 48], type = 'l', main = paste0('mu(intercept)_', names[i + 1]))
  plot(basura$process[, 42 + i * 48], type = 'l', main = paste0('mu(elev)_', names[i + 1]))
  plot(basura$process[, 43 + i * 48], type = 'l', main = paste0('mu(dist)_', names[i + 1]))
  plot(basura$process[, 44 + i * 48], type = 'l', main = paste0('1/sigma_k^2_', names[i + 1])) #1/simga_k^2
  plot(basura$process[, 45 + i * 48], type = 'l', main = paste0('decay_', names[i + 1])) # decay
  plot(basura$process[, 46 + i * 48], type = 'l', main = paste0('varsgima_', names[i + 1])) # varsigma
  plot(basura$process[, 47 + i * 48], type = 'l', main = paste0('varphi_', names[i + 1])) #varphi
  plot(basura$process[, 48 + i * 48], type = 'l', main = paste0('prec.coast_', names[i + 1])) # prec.coast
  dev.off()
}



plot(basura$params[, 1], type = 'l')

plot(basura$process[, 50], type = 'l')
