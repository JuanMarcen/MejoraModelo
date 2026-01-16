# fitting of convolution models M1 M2 and M3
# both quantiles
tau <- 0.5

library(qs)
df <- qread('df_jun_ag.qs')
df <- qread('C:/Users/PC/Desktop/Juan/df_jun_ag.qs')
stations <- readRDS('stations.rds')

# data preparation
library(dplyr)
library(lubridate)
df <- df %>%
  select(Date, station, Y, l, t, s.1, c.1, elev, dist, g300, g500, g700) %>%
  mutate(
    month = month(Date),
    `t:month6` = ifelse(month == 6, t, 0),
    `t:month7` = ifelse(month == 7, t, 0),
    `t:month8` = ifelse(month == 8, t, 0)
  )


Y <- df$Y
X <- matrix(0, nrow = length(Y), ncol = 0)

# definition of arguments for MCMC function
V.M1 <- cbind(1, df$s.1, df$c.1, df$`t:month6`, df$`t:month7`, df$`t:month8`)
V.M2 <- cbind(1, df$s.1, df$c.1, df$g300, df$g500, df$g700)
V.M3 <- cbind(1, df$s.1, df$c.1, df$g300, df$g500, df$g700, df$`t:month6`, df$`t:month7`, df$`t:month8`)

# X_alpha
aux.vars <- c('elev', 'dist')
X_alpha.M1 <- list()
X_alpha.M2 <- list()
X_alpha.M3 <- list()
#intercept
X_alpha.M1[[1]] <- cbind(1, unique(df$elev), unique(df$dist))
X_alpha.M2[[1]] <- cbind(1, unique(df$elev), unique(df$dist))
X_alpha.M3[[1]] <- cbind(1, unique(df$elev), unique(df$dist))
for (i in 2:ncol(V.M1)){
  X_alpha.M1[[i]] <- matrix(1, ncol = 1, nrow = nrow(stations))
}
for (i in 2:ncol(V.M2)){
  X_alpha.M2[[i]] <- matrix(1, ncol = 1, nrow = nrow(stations))
}
for (i in 2:ncol(V.M3)){
  X_alpha.M3[[i]] <- matrix(1, ncol = 1, nrow = nrow(stations))
}

dist <- readRDS('maps coast/dist.matrix.rds')

#priors
M <- matrix(0, nrow = 0, ncol = 1)
P <- matrix(0, nrow = 0, ncol = 0)

M_beta_alpha.M1 <- list()
P_beta_alpha.M1 <- list()
M_beta_alpha.M1[[1]] <- rep(0, 3)
P_beta_alpha.M1[[1]] <- 0.0001 * diag(3)
for (i in 2:length(X_alpha.M1)){
  M_beta_alpha.M1[[i]] <- rep(0, 1)
  P_beta_alpha.M1[[i]] <- 0.0001 * diag(1)
}

M_beta_alpha.M2 <- list()
P_beta_alpha.M2 <- list()
M_beta_alpha.M2[[1]] <- rep(0, 3)
P_beta_alpha.M2[[1]] <- 0.0001 * diag(3)
for (i in 2:length(X_alpha.M2)){
  M_beta_alpha.M2[[i]] <- rep(0, 1)
  P_beta_alpha.M2[[i]] <- 0.0001 * diag(1)
}

M_beta_alpha.M3 <- list()
P_beta_alpha.M3 <- list()
M_beta_alpha.M3[[1]] <- rep(0, 3)
P_beta_alpha.M3[[1]] <- 0.0001 * diag(3)
for (i in 2:length(X_alpha.M3)){
  M_beta_alpha.M3[[i]] <- rep(0, 1)
  P_beta_alpha.M3[[i]] <- 0.0001 * diag(1)
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
alpha.M1 <- matrix(0.1, nrow = nrow(stations), ncol = ncol(V.M1))
alpha.M2 <- matrix(0.1, nrow = nrow(stations), ncol = ncol(V.M2))
alpha.M3 <- matrix(0.1, nrow = nrow(stations), ncol = ncol(V.M3))
prec <- 1

hp.M1 <- matrix(0, nrow = 5, ncol = ncol(V.M1))
hp.M1[1, ] <- 1 #precision
hp.M1[2, ] <- 3/600 #decay
hp.M1[3, ] <- 1 #varsigma
hp.M1[4, ] <- 3/900 # varphi
hp.M1[5, ] <- 1

hp.M2 <- matrix(0, nrow = 5, ncol = ncol(V.M2))
hp.M2[1, ] <- 1 #precision
hp.M2[2, ] <- 3/600 #decay
hp.M2[3, ] <- 1 #varsigma
hp.M2[4, ] <- 3/900 # varphi
hp.M2[5, ] <- 1

hp.M3 <- matrix(0, nrow = 5, ncol = ncol(V.M3))
hp.M3[1, ] <- 1 #precision
hp.M3[2, ] <- 3/600 #decay
hp.M3[3, ] <- 1 #varsigma
hp.M3[4, ] <- 3/900 # varphi
hp.M3[5, ] <- 1

beta_alpha.M1 <- list()
for (i in 1:length(X_alpha.M1)){
  beta_alpha.M1[[i]] <- rep(0.1, times = ncol(X_alpha.M1[[i]]))
}

beta_alpha.M2 <- list()
for (i in 1:length(X_alpha.M2)){
  beta_alpha.M2[[i]] <- rep(0.1, times = ncol(X_alpha.M2[[i]]))
}

beta_alpha.M3 <- list()
for (i in 1:length(X_alpha.M3)){
  beta_alpha.M3[[i]] <- rep(0.1, times = ncol(X_alpha.M3[[i]]))
}


# more constants 
N <- nrow(df)
n <- nrow(stations)
p <- 0

r.M1 <- ncol(V.M1)
r.M2 <- ncol(V.M2)
r.M3 <- ncol(V.M3)

p_alpha.M1 <- unlist(lapply(X_alpha.M1, ncol))
p_alpha.M2 <- unlist(lapply(X_alpha.M2, ncol))
p_alpha.M3 <- unlist(lapply(X_alpha.M3, ncol))

s <- rep(0:39, each = 5888)

nSims <- 10
nThin <- 1
nBurnin <- 10
nReport <- 100

#more distances
dist_coast <- readRDS('maps coast/dist.vec.rds')
dist_coast_points <- readRDS('maps coast/dist.coast.points2.rds')

dmatcoast_conv <- readRDS('maps coast/phimat.rds')
drmat_conv <- readRDS('maps coast/dr.rds')
lencoast_conv <- drmat_conv[1, 200]

# --- MCMC ---
Rcpp::sourceCpp("Metropolis-within-Gibbs/mcmc3.cpp")
# M1
repeat {
  conv.M1.q0.50 <- try(spQuantileRcpp(
    tau = 0.50,
    Y = Y,
    X = X,
    V = V.M1,
    X_alpha = X_alpha.M1,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha.M1,
    P_beta_alpha = P_beta_alpha.M1,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha.M1,
    prec = prec,
    hp = hp.M1,
    beta_alpha = beta_alpha.M1,
    N = N,
    n = n,
    p = p,
    r = r.M1,
    p_alpha = p_alpha.M1,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = 2
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv.M1.q0.50, "try-error")) {
    return(conv.M1.q0.50)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}
saveRDS(conv.M1.q0.50, 'conv.models/conv.M1.q0.50.rds')
repeat {
  conv.M1.q0.95 <- try(spQuantileRcpp(
    tau = 0.95,
    Y = Y,
    X = X,
    V = V.M1,
    X_alpha = X_alpha.M1,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha.M1,
    P_beta_alpha = P_beta_alpha.M1,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha.M1,
    prec = prec,
    hp = hp.M1,
    beta_alpha = beta_alpha.M1,
    N = N,
    n = n,
    p = p,
    r = r.M1,
    p_alpha = p_alpha.M1,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = 2
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv.M1.q0.95, "try-error")) {
    return(conv.M1.q0.95)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}
saveRDS(conv.M1.q0.95, 'conv.models/conv.M1.q0.95.rds')

# M2
repeat {
  conv.M2.q0.50 <- try(spQuantileRcpp(
    tau = 0.50,
    Y = Y,
    X = X,
    V = V.M2,
    X_alpha = X_alpha.M2,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha.M2,
    P_beta_alpha = P_beta_alpha.M2,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha.M2,
    prec = prec,
    hp = hp.M2,
    beta_alpha = beta_alpha.M2,
    N = N,
    n = n,
    p = p,
    r = r.M2,
    p_alpha = p_alpha.M2,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = 2
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv.M2.q0.50, "try-error")) {
    return(conv.M2.q0.50)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}
saveRDS(conv.M2.q0.50, 'conv.models/conv.M2.q0.50.rds')
repeat {
  conv.M2.q0.95 <- try(spQuantileRcpp(
    tau = 0.95,
    Y = Y,
    X = X,
    V = V.M2,
    X_alpha = X_alpha.M2,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha.M2,
    P_beta_alpha = P_beta_alpha.M2,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha.M2,
    prec = prec,
    hp = hp.M2,
    beta_alpha = beta_alpha.M2,
    N = N,
    n = n,
    p = p,
    r = r.M2,
    p_alpha = p_alpha.M2,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = 2
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv.M2.q0.95, "try-error")) {
    return(conv.M2.q0.95)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}
saveRDS(conv.M2.q0.95, 'conv.models/conv.M2.q0.95.rds')

# M3
repeat {
  conv.M3.q0.50 <- try(spQuantileRcpp(
    tau = 0.50,
    Y = Y,
    X = X,
    V = V.M3,
    X_alpha = X_alpha.M3,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha.M3,
    P_beta_alpha = P_beta_alpha.M3,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha.M3,
    prec = prec,
    hp = hp.M3,
    beta_alpha = beta_alpha.M3,
    N = N,
    n = n,
    p = p,
    r = r.M3,
    p_alpha = p_alpha.M3,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = 2
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv.M3.q0.50, "try-error")) {
    return(conv.M3.q0.50)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}
saveRDS(conv.M3.q0.50, 'conv.models/conv.M3.q0.50.rds')
repeat {
  conv.M3.q0.95 <- try(spQuantileRcpp(
    tau = 0.95,
    Y = Y,
    X = X,
    V = V.M3,
    X_alpha = X_alpha.M3,
    dist = dist,
    dist_coast = dist_coast,
    dist_coast_point = dist_coast_points,
    dmatcoast_conv = dmatcoast_conv,
    drmat_conv = drmat_conv,
    lencoast_conv = lencoast_conv,
    M = M,
    P = P,
    M_beta_alpha = M_beta_alpha.M3,
    P_beta_alpha = P_beta_alpha.M3,
    da = da,
    db = db,
    ga = ga,
    gb = gb,
    ra = ra,
    rb = rb,
    na = na, 
    nb = nb,
    beta = beta,
    alpha = alpha.M3,
    prec = prec,
    hp = hp.M3,
    beta_alpha = beta_alpha.M3,
    N = N,
    n = n,
    p = p,
    r = r.M3,
    p_alpha = p_alpha.M3,
    nSims = nSims,
    nThin = nThin,
    nBurnin = nBurnin,
    nReport = nReport,
    s = s,
    model = 2
  ),
  silent = TRUE
  )
  
  
  # Si NO dio error → salimos del bucle
  if (!inherits(conv.M3.q0.95, "try-error")) {
    return(conv.M3.q0.95)
  }
  
  message("⚠️  Error en spTm() → reintentando...")
  
  # Opcional: pequeña pausa para no saturar la CPU
  Sys.sleep(0.5)
}
saveRDS(conv.M3.q0.95, 'conv.models/conv.M3.q0.95.rds')

#---- traceplots pf GP (do in big computer) ----
names <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700', 
           't_month6', 't_month7', 't_month8')
for (j in 0:8){
  filename <- paste0('Metropolis-within-Gibbs/trplots/trplots_', names[j +1], '.pdf')
  pdf(filename, width = 12, height = 15)
  par(mfrow = c(8, 5))
  for (i in 1:40){
    plot(basura$process[, i + j * 48], type = 'l', main = paste0(names[j + 1], '(', stations$NAME2[i], ')'))
  }
  dev.off()
}

# traceplots of parameters in GP
for (i in 0:8){
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

############################################################
# for launching in big computer
# KRIGING AND MAP

# distances
dmatcoast_conv2 <- readRDS('maps coast/phimat.grid.rds')
newcoords <- readRDS('grid_km.rds')
coords <- readRDS('coords.stations.rds')

dist_coast_points.grid <- readRDS('maps coast/dist.coast.points.grid.rds')
dist_coast.grid <- readRDS('maps coast/dist.vec.grid.rds')

dist_coast_points.comb <- readRDS('maps coast/r.stations.grid.rds')

# function for kriging all parameters
kriging <- function(chain, X_alpha, vars){
  
  out <- list()
  
  # intercept different from the rest
  w <- chain$process[, 1:40]
  Z <- X_alpha[[1]]
  w.def <- w - t(Z %*% t(chain$process[, 41:43]))
  out[[vars[1]]] <- krigeBayesRcpp(
    w = w.def,
    hp = chain$process[ , 44:48, drop = FALSE],
    coords = coords,
    newcoords = newcoords,
    dr = drmat_conv,
    dcoast = dmatcoast_conv,
    newdcoast = dmatcoast_conv2,
    lencoast = lencoast_conv,
    dvec = dist_coast,
    newdvec = dist_coast.grid,
    dmatc = dist_coast_points,
    newdmatc = dist_coast_points.grid,
    combdmatc = dist_coast_points.comb,
    model = 2
  )
  
  for (i in 0:(length(vars) - 2)){
    # substract mean of the GP
    w <- chain$process[, 49:88 + i * 46]
    Z <- X_alpha[[i + 2]]
    w.def <- w - t(Z %*% t(chain$process[, 89 + i * 46]))
    
    out[[vars[i + 2]]] <- krigeBayesRcpp(
      w = w.def,
      hp = chain$process[ , 90:94 + i * 46, drop = FALSE],
      coords = coords,
      newcoords = newcoords,
      dr = drmat_conv,
      dcoast = dmatcoast_conv,
      newdcoast = dmatcoast_conv2,
      lencoast = lencoast_conv,
      dvec = dist_coast,
      newdvec = dist_coast.grid,
      dmatc = dist_coast_points,
      newdmatc = dist_coast_points.grid,
      combdmatc = dist_coast_points.comb,
      model = 2
    )
  }
  
  return(out)
}

vars.M1 <- c('intercept', 's.1', 'c.1', 
             't_month6', 't_month7', 't_month8')
vars.M2 <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700')
vars.M3 <- c('intercept', 's.1', 'c.1', 'g300', 'g500', 'g700', 
          't_month6', 't_month7', 't_month8')

kr.conv.M1.q0.50 <- kriging(conv.M1.q0.50, X_alpha.M1, vars.M1)
saveRDS(kr.conv.M1.q0.50, 'conv.models/kr.conv.M1.q0.50.rds')
kr.conv.M1.q0.95 <- kriging(conv.M1.q0.95, X_alpha.M1, vars.M1)
saveRDS(kr.conv.M1.q0.95, 'conv.models/kr.conv.M1.q0.95.rds')

kr.conv.M2.q0.50 <- kriging(conv.M2.q0.50, X_alpha.M2, vars.M2)
saveRDS(kr.conv.M2.q0.50, 'conv.models/kr.conv.M2.q0.50.rds')
kr.conv.M2.q0.95 <- kriging(conv.M2.q0.95, X_alpha.M2, vars.M2)
saveRDS(kr.conv.M2.q0.95, 'conv.models/kr.conv.M2.q0.95.rds')

kr.conv.M3.q0.50 <- kriging(conv.M3.q0.50, X_alpha.M3, vars.M3)
saveRDS(kr.conv.M3.q0.50, 'conv.models/kr.conv.M3.q0.50.rds')
kr.conv.M3.q0.95 <- kriging(conv.M3.q0.95, X_alpha.M3, vars.M3)
saveRDS(kr.conv.M3.q0.95, 'conv.models/kr.conv.M3.q0.95.rds')

save(conv.M1.q0.50, conv.M1.q0.95, kr.conv.M1.q0.50, kr.conv.M1.q0.95,
     conv.M2.q0.50, conv.M2.q0.95, kr.conv.M2.q0.50, kr.conv.M2.q0.95,
     conv.M3.q0.50, conv.M3.q0.95, kr.conv.M3.q0.50, kr.conv.M3.q0.95,
     file = 'M1.M2.M3.conv.RData')

# map plots (in big computer)
library(viridis)
library(ggplot2)
library(sf)
library(sp)
limits <- readRDS('limits.rds')
background <- readRDS('background.rds')
grid <- readRDS('grid.rds')
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")],
      data = stations,
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    ),
    'sf'
  ),
  2062
)

plot.gp <- function(kriging, grid, stations, limits, background, var, tau){
  mu <- colMeans(kriging[[var]], na.rm = T)
  
  grid_coords <- cbind(st_coordinates(grid), mu = mu)
  grid_coords <- na.omit(grid_coords)
  grid.new <- st_sf(mu = round(as.vector(mu), 3), geometry = grid)
  
  ggplot(data = background) + 
    geom_sf(fill = "antiquewhite") + 
    xlab("Longitud (º)") + ylab("Latitud (º)") + ggtitle(bquote( .(var) * .(' (') * tau *.(' = ') *.(0.5) *.(')'))) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6, angle = 90),
          axis.title=element_text(size = 10, face = "bold")) + 
    geom_tile(data = grid_coords, aes(X, Y, fill = mu)) +
    geom_tile(data = grid.new, ggplot2::aes(x = st_coordinates(grid.new)[, 1], y = st_coordinates(grid.new)[, 2], fill = mu)) +
    geom_sf(data = stations, aes(color = mu), 
            shape = 21,        # forma con relleno y borde
            color = "black",   # contorno
            size = 3,
            stroke = 1 )+
    scale_fill_gradient2( low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                          space = "Lab", midpoint = 0, limits = c(-5, 5), name = "Distance (km)") +
    scale_color_gradient2( low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                           space = "Lab", midpoint = 0, limits = c(-5, 5), name = "Distance (km)") +
    
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])
}

# for all vars and save images
for (var in vars){
  filename <- paste0('Metropolis-within-Gibbs/mapsGP/conv/GP.', var, '.pdf')
  g <- plot.gp(kr.conv, grid, stations, limits, background, var, 0.50)
  ggsave(filename, g, height = 8, width = 10)
}

for (var in vars){
  filename <- paste0('Metropolis-within-Gibbs/mapsGP/coastal/GP.', var, '.pdf')
  g <- plot.gp(kr.coastal, grid, stations, limits, background, var, 0.50)
  ggsave(filename, g, height = 8, width = 10)
}


##############################
plot(basura$process[, 1 + 4 * 48], type = 'l', main = paste0(names[4 + 1], '(', stations$NAME2[1], ')'),
     xlab = 'Sample', ylab = 'Value')
