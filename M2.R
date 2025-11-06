# M1: g300, g500 and g700

# only need for the moment the datframe 
# with elev, dist, l, t, harmonics, and month dummy

rm(list = ls())

library(qs)

stations <- readRDS('stations.rds')
df_jun_ag <- qread('df_jun_ag.qs')


# subset of the data frame
library(dplyr)
library(lubridate)
df.M2 <- df_jun_ag %>%
  select(Date, station, Y, l, t, g300, g500, g700, elev, dist) 

#----BAYESIAN MODELS----
library(spTReg)
library(quantreg)
library(sf)
library(sp)
library(coda)
stations_dist <- readRDS('stations_dist.rds')
#coordinates in km of stations
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")], 
      data = stations[c("STAID", "STANAME", "LON", "LAT", "HGHT","color",'NAME1','NAME2')],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

coords_km <- st_coordinates(stations) / 1000

# formula of model
formula <- as.formula('Y ~ g300 + g500 + g700 + elev + dist' )
vars <- c('g300', 'g500', 'g700')

# BAYESIAN MODELS (Run at other computer)
source('metrics_bay.R')

# PARALLEL COMPUTATION  
n.fixed.eff <- length(vars) + 3 # intercept, dist and elev
n.random.eff <- length(vars) + 1 # random effects only for the covariates and intercept
# for the initial values of the random effects we need one strating value per variable and station
# if desired it can be fized to only one value. In this case we just put the fixed value
n.stations <- nrow(stations)
# different starting points
set.seed(2025)

start_beta_1 <- rnorm(n.fixed.eff, mean = 0, sd = 0.01)   # close to 0
start_beta_2 <- rnorm(n.fixed.eff, mean = 0.05, sd = 0.02) # little positive bias
start_beta_3 <- rnorm(n.fixed.eff, mean = -0.05, sd = 0.02) # little negative bias

inic_proc_1 <- matrix(rnorm(n.random.eff * n.stations, mean = 0, sd = 0.01),
                      nrow = n.stations, ncol = n.random.eff)
inic_proc_2 <- matrix(rnorm(n.random.eff, mean = 0.02, sd = 0.02),
                      nrow = n.stations, ncol = n.random.eff)
inic_proc_3 <- matrix(rnorm(n.random.eff, mean = -0.02, sd = 0.02),
                      nrow = n.stations, ncol = n.random.eff)

inits_list <- list(
  list(start_beta = start_beta_1, inic_proc = inic_proc_1),
  list(start_beta = start_beta_2, inic_proc = inic_proc_2),
  list(start_beta = start_beta_3, inic_proc = inic_proc_3)
)

library(parallel)

# number of nodes
cl <- makeCluster(3)

# export objects and functions needed
clusterExport(cl, c("mod_bay", "formula", "df.M2", "vars", "coords_km"))

# load needed library
clusterEvalQ(cl, library(spTReg))

# Run 3 parallel fitting of models
resultados <- parLapply(
  cl,
  X = inits_list,
  fun = function(inits){
    mod_bay(
      formula = formula, 
      data = df.M2, 
      tau = 0.95,                      # Aquí tu tau fijo
      vars = vars, 
      coords = coords_km, 
      start_beta = inits$start_beta,   # Cambia beta
      inic_procesos = inits$inic_proc, # Cambia procesos
      n.samples = 100, 
      n.burnin = 100, 
      n.thin = 1, 
      n.report = 1
    )
  }
)

stopCluster(cl)

mod1 <- resultados[[1]]
mod2 <- resultados[[2]]
mod3 <- resultados[[3]]

plot(as.numeric(mod3$p.params.samples[, 1]), type = 'l')
lines(as.numeric(mod1$p.params.samples[, 1]), col = 'tomato')
lines(as.numeric(mod2$p.params.samples[, 1]), col = 'grey')

#final chain of parameters
#join the 3 chains estimated
final.chain <-rbind(
  mod1$p.params.samples,
  mod2$p.params.samples,
  mod3$p.params.samples
                    )
final.chain <- as.data.frame(final.chain)
plot(as.numeric(final.chain[,1]), type = 'l')


#----METRICS TO ASSES GOODNESS OF FIT----
#mean values of the chains of parameters
betas_q0.95 <- betas(vars, final.chain)


#predictions
elev_sc <- scale(stations_dist$HGHT)
dist_sc <- scale(stations_dist$DIST)

pred_q0.95 <- predictions(vars, betas_q0.95, df.M2, cuantil = 0.95)
pred_q0.95 <- cbind(df.M2[,c('Date', 'station', 'Y')], pred_q0.95)

# R1
R1_bay_q0.95 <- R1_bay(pred_q0.95, 0.95, df.M2)


# rho
rho_q0.95 <- rho_bay(pred_q0.95, 0.95)


# 
# save(df.M1,
#      vars,
#      stations,
#      stations_dist,
#      elev_sc, dist_sc,
#      local.models.q0.50,
#      local.models.q0.95,
#      formula,
#      mod_q0.50_bay,
#      mod_q0.95_bay,
#      pred_q0.95,
#      pred_q0.50,
#      R1_bay_q0.50,
#      R1_bay_q0.95,
#      rho_q0.50,
#      rho_q0.95,
#      file = 'dataM1.RData')
