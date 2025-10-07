# M1: Only harmonics, elev, dist and t:month

# only need for the moment the datframe 
# with elev, dist, l, t, harmonics, and month dummy

rm(list = ls())

stations <- readRDS('stations.rds')
df <- readRDS('df_jun_ag.rds') ## MODIFICAR EN UN FUTURO

# subset of the data frame
library(dplyr)
library(lubridate)
df.M1 <- df %>%
  select(Date, station, Y, l, t, s.1, c.1, elev, dist) %>%
  mutate(
    month = month(Date),
    `t:month6` = ifelse(month == 6, t, 0),
    `t:month7` = ifelse(month == 7, t, 0),
    `t:month8` = ifelse(month == 8, t, 0)
  )

# trying things
ind <- which(df$station == stations$STAID[1])
formula <- as.formula('Y ~ s.1 + c.1 + `t:month6` + `t:month7` + `t:month8` ' )

library(quantreg)
mod_null <- rq(Y ~ 1, tau = 0.50, data = df.M1, subset = ind)
basura <- rq(formula, tau = 0.50, data = df.M1, subset = ind)
1 - basura$rho/mod_null$rho
summary(basura)
step(basura)

mod_null <- rq(Y ~ 1, tau = 0.95, data = df.M1, subset = ind)
basura2 <- rq(formula, tau = 0.95, data = df.M1, subset = ind)
summary(basura2)
1 - basura2$rho/mod_null$rho

source('functions.R')

df_dia <- rho_day(basura, basura2, df.M1)
df_year <- rho_year(basura, basura2, df.M1)
par(mfrow = c(1,2))
plot(1:92, df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (armónicos)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:92, df_dia$rho_l_q0.95)
abline(h = 0.95, col = 'red')

plot(1:64, df_year$rho_t_q0.5, type='l', 
     main = 'Madrid (Retiro) (años) (armónicos)',
     ylab = expression(rho[t](tau)), xlab = 't', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:64, df_year$rho_t_q0.95)
abline(h = 0.95, col = 'red')


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

# local models for starting points NO GUARDA BIEN LAS DUMMYS
vars <- c('s.1', 'c.1', '`t:month6`', '`t:month7`', '`t:month8`')

local.models <- function(vars, tau, df, which.dummy = NULL){
  vars <- ifelse(grepl("^I\\(.*\\^2\\)$", vars),
                 paste0("`", vars, "`"),
                 vars)
  
  var.dummy <- grep(paste0(':', which.dummy), vars)
  var.dummy <- vars[var.dummy]
  vars <- c(setdiff(vars, var.dummy), var.dummy)
  formula <- as.formula(paste('Y ~', paste(vars, collapse = '+')))
  
  cat(deparse(formula), '\n')
  
  
  vars <- gsub('`', '', vars)
  
  
  if (!is.null(which.dummy)){
    aux <- length(unique(df[, which.dummy])) - 1
  }else{
    aux <- 0
  }
  
  names.vars <- setdiff(vars, var.dummy)
  if (!is.null(which.dummy)){
    for (i in levels(df[, which.dummy])){
      names.vars <- c(names.vars, paste0(var.dummy, i))
    }
  }else{
    names.vars <- vars
  }
  
  
  modelos_finales_df <- matrix(NA, nrow=dim(stations)[1], ncol = length(vars) + aux + 4)
  modelos_finales_df <- as.data.frame(modelos_finales_df)
  colnames(modelos_finales_df) <- c('stations', 'NAME2', 'intercept', names.vars, 'R1')
  modelos_finales_df$stations <- stations$STAID
  modelos_finales_df$NAME2 <- stations$NAME2
  
  for(i in 1:dim(stations)[1]){
    cat(stations$NAME2[i], '\n')
    
    ind <- which(df$station == stations$STAID[i])
    cat(ind[c(1, 2)], '\n')
    
    mod <- rq(formula = formula, tau = tau, data = df, subset = ind)
    mod_nulo <- rq(formula = Y ~ 1, tau = tau, data = df, subset = ind)
    
    cat('Intercepto:', coef(mod)[1], '\n')
    
    R1 <- 1 - mod$rho / mod_nulo$rho
    cat('R1: ', R1, '\n')
    
    modelos_finales_df[i, 3:(2 + length(vars) + 1 + aux)]<-c(coef(mod))
    modelos_finales_df[i, 'R1'] <- R1
    
  }
  
  return(list(
    modelos = modelos_finales_df,
    formula = formula))
}

local.models.q0.95 <- local.models(vars, 0.95, df.M1)[["modelos"]]
formula <- local.models(vars, 0.95, df.M1)[["formula"]]
local.models.q0.50 <- local.models(vars, 0.50, df.M1)[["modelos"]]

#formula update (add elev and dist)
formula <- update(formula, .~. + elev + dist)

# STARTING POINTS
starting.values <- function(local.models){
  start_beta <- apply(
    local.models[,3:(ncol(local.models) - 1)],
    MARGIN = 2, FUN = mean)
  
  elev_sc <- scale(stations_dist$HGHT)
  dist_sc <- scale(stations_dist$DIST)
  mod <- lm(local.models$intercept ~ elev_sc + dist_sc)
  
  elev_inic <- coef(mod)[2]
  dist_inic <- coef(mod)[3]
  
  inic_procesos <- as.matrix(
    cbind(
      mod$residuals,
      sweep(local.models[,4:(ncol(local.models) - 1)],
            2, start_beta[-1])))
  inic_procesos <- inic_procesos + rnorm(40 * (ncol(local.models)-3), 0, 0.1)
  
  start_beta <- c(start_beta, elev_inic, dist_inic)
  
  return(list(
    start_beta = start_beta,
    inic_procesos = inic_procesos
  ))
}

# OJO QUE HAY DUMMY VARIABLES

start_beta_q0.95 <- starting.values(local.models.q0.95)[['start_beta']]

inic_proc_q0.95 <- starting.values(local.models.q0.95)[['inic_procesos']]

start_beta_q0.50 <- starting.values(local.models.q0.50)[['start_beta']]
inic_proc_q0.50 <- starting.values(local.models.q0.95)[['inic_procesos']]


# BAYESIAN MODELS (Run at other computer)
mod_bay <- function(formula, data, tau, vars, 
                    coords, start_beta, inic_procesos,
                    n.samples, n.burnin, n.thin, n.report){
  
  vars <- gsub('`', '', vars)
  ind <- length(vars) + 3
  
  mod<-spTm(formula,
            data = data,
            method = 'q',
            quantile = tau,
            coords = coords,
            v = as.matrix(cbind(1,data[, vars])),
            priors = list(
              "beta" = list(M = rep(0, ind), P = 0.0001 * 
                              diag(ind)),
              "sigma" = c(0.1, 0.1),
              "phi" = c(38, 7400),
              "mu" = c(0, 0.0001)),
            starting = list(
              "beta" = start_beta,
              "sigma" = 1,
              "alpha" = inic_procesos,
              "hp" = c("mu" = 0, "sigma" = 1, "phi" = 3 / 600)),
            n.samples = n.samples,
            n.burnin = n.burnin,
            n.thin = n.thin,
            n.report = n.report
  )
  
  return(mod)
}

# do the try-catch as in JorgeScript.R
mod_q0.95_bay <- mod_bay(
  formula = formula, 
  data = df.M1, 
  tau = 0.95, 
  vars = vars, 
  coords = coords_km, 
  start_beta = start_beta_q0.95, 
  inic_procesos = inic_proc_q0.95, 
  n.samples = 100, 
  n.burnin = 100, 
  n.thin = 1, 
  n.report = 10)

mod_q0.50_bay <- mod_bay(
  formula = formula, 
  data = df.M1, 
  tau = 0.50, 
  vars = vars, 
  coords = coords_km, 
  start_beta = start_beta_q0.50, 
  inic_procesos = inic_proc_q0.50, 
  n.samples = 100, 
  n.burnin = 100, 
  n.thin = 1, 
  n.report = 10)

#----METRICS TO ASSES GOODNESS OF FIT----
#mean values of the chains of parameters
traducir_nombres_coef <- function(nombres_coef) { 
  traducidos <- character(length(nombres_coef))
  
  for (i in seq_along(nombres_coef)) {
    nombre <- nombres_coef[i]
    
    if (grepl("^poly\\((.+), 2\\)1$", nombre)) {
      base <- sub("^poly\\((.+), 2\\)1$", "\\1", nombre)
      traducidos[i] <- base
    } else if (grepl("^poly\\((.+), 2\\)2$", nombre)) {
      base <- sub("^poly\\((.+), 2\\)2$", "\\1", nombre)
      traducidos[i] <- paste0("I(", base, "^2)")
    } else {
      traducidos[i] <- nombre
    }
  }
  
  return(traducidos)
}

betas <- function(vars, mod, cuantil){
  vars <- gsub('`', '', vars)
  ind <- length(vars) + 1
  params <- as.data.frame(mod$p.params.samples)
  
  tr <- traducir_nombres_coef(colnames(params)[2:ind])
  colnames(params)[2:ind] <- gsub('`', '', tr)
  #intercepto
  int <- mean(params[['(Intercept)']])
  int <- rep(int, length=dim(stations)[1])
  
  #beta_fija
  betas_fijas <- apply(params[, vars], 2, mean)
  betas_fijas <- matrix(rep(betas_fijas, each = dim(stations)[1]), nrow = dim(stations)[1])
  
  elev <- mean(params[,'elev'])
  elev <- matrix(rep(elev, each = dim(stations)[1]), nrow = dim(stations)[1])
  dist <- mean(params[,'dist'])
  dist <- matrix(rep(dist, each = dim(stations)[1]), nrow = dim(stations)[1])
  
  #betas espaciales (beta1,...,beta16)
  cols <- grep('beta', names(params), value=T)
  mu <- apply(params[, cols], 2, mean)
  betas_esp <- matrix(mu, nrow = dim(stations)[1])
  
  #juntar en data frame
  betas <- cbind(int, elev, dist, betas_fijas, betas_esp)
  betas <- as.data.frame(betas, row.names = stations$NAME2)
  colnames(betas) <- c('intercept', 'elev', 'dist',
                       vars, paste0('beta',1:(length(vars) + 1)))
  
  
  return(betas)
}

betas_q0.95 <- betas(vars, mod_q0.95_bay)
betas_q0.50 <- betas(vars, mod_q0.50_bay)

#predictions
predictions <- function(vars, betas, df, cuantil){
  vars <- gsub('`', '', vars)
  pred <- numeric(nrow(df))
  for (i in 1:dim(stations)[1]){
    ind <- which(df$station == stations$STAID[i])
    for (j in ind){
      #inicializar en interceptos
      comp_esp <- betas[i, 'beta1'] #intercepto espacial para la estacion i
      
      comp_fija <- betas[i, 'intercept'] #intercepto fijo para la estacion i
      
      for (k in 1:(length(vars))){ #beta 1 es la componente espacial del intercepto
        comp_esp <- comp_esp + betas[i, paste0('beta', k + 1)] * df[j, vars[k]]
        comp_fija <- comp_fija + betas[i, vars[k]] * df[j, vars[k]]
      }
      
      # if (cuantil==0.5){
      #   pred[j]<-comp_esp+comp_fija + betas[i,'elev']*elev_sc[i] + betas[i,'dist']*dist_sc[i]
      # }
      
      if (cuantil == 0.95 || cuantil == 0.90 || cuantil == 0.75){
        pred[j] <- comp_esp + comp_fija + betas[i, 'elev'] * elev_sc[i] + betas[i, 'dist'] * dist_sc[i]
      }
      
    }
  }
  
  return(pred)
}

elev_sc <- scale(stations_dist$HGHT)
dist_sc <- scale(stations_dist$DIST)

pred_q0.95 <- predictions(vars, betas_q0.95, df.M1, cuantil = 0.95)
pred_q0.95 <- cbind(df.M1[,c('Date', 'station', 'Y')], pred_q0.95)

pred_q0.50 <- predictions(vars, betas_q0.50, df.M1, cuantil = 0.50)
pred_q0.50 <- cbind(df.M1[,c('Date', 'station', 'Y')], pred_q0.50)


# R1
check <- function(u, tau) {
  return(u * (tau - (u < 0)))  # Implements the quantile loss function
}

R1_bay <- function(pred, tau, data){
  pred_clean <- na.omit(pred)
  
  #dataframe para global
  df <- matrix(NA, nrow=1, ncol=1)
  df <- as.data.frame(df)
  colnames(df) <- c('R1_bay_global')
  
  #dataframe para locales
  df_local <- matrix(NA, nrow=dim(stations)[1], ncol=1)
  df_local <- as.data.frame(df_local, row.names = stations$NAME2)
  colnames(df_local) <- c('R1_bay_local')
  
  #modelos nulos, son para todas variables igual
  mod_nulo_f <- rq(Y ~ as.factor(station), data = data, tau = tau)
  
  rho_estacion <- rep(NA, dim(stations)[1])
  R1_nulo_est <- rep(NA, dim(stations)[1])
  for (j in 1:length(rho_estacion)){
    ind <- which(pred_clean$station == stations$STAID[j])
    rho_estacion[j] <- sum(check(pred_clean$Y[ind] - pred_clean[ind,paste0('pred_q',format(tau, nsmall = 2))],tau = tau))
    R1_nulo_est[j] <- sum(check(mod_nulo_f$residuals[ind], tau = tau))
  }
  
  df['R1_bay_global'] <- 1 - sum(rho_estacion) / mod_nulo_f$rho
  df_local['R1_bay_local'] <- 1 - rho_estacion / R1_nulo_est
  
  return(list(R1_globales = df, R1_locales = df_local))
}

R1_bay_q0.95 <- R1_bay(pred_q0.95, 0.95, df.M1)
R1_bay_q0.50 <- R1_bay(pred_q0.50, 0.50, df.M1)

# rho
rho_bay <- function(predicciones, tau){
  
  #global
  #dataframe para global
  df <- matrix(NA, nrow=1, ncol=1)
  df <- as.data.frame(df)
  colnames(df) <- c('rho_bay_global')
  pred <- predicciones[[paste0('pred_q',format(tau, nsmall = 2))]]
  dif<- predicciones$Y - pred
  df['rho_bay_global'] <- sum(dif < 0, na.rm = T) / ( 40 * 64 * 92)
  
  #estaciones
  df_est <- matrix(NA, nrow=dim(stations)[1], ncol=1)
  df_est <- as.data.frame(df_est, row.names = stations$NAME2)
  colnames(df_est) <- c('rho_bay_est')
  for (i in 1:dim(stations)[1]){
    ind <- which(predicciones$station == stations$STAID[i])
    dif <- predicciones$Y[ind]-pred[ind]
    df_est[i,] <- sum(dif < 0, na.rm = T) / (64 * 92)
  }
  
  #dias
  df_dia_list <- list() #lista para dias por estacion
  
  day_month <- unique(format(predicciones$Date, "%d-%m"))
  df_dia <- matrix(NA, nrow = length(day_month), ncol=1)
  df_dia <- as.data.frame(df_dia, row.names = day_month)
  colnames(df_dia) <- c('rho_bay_dia')
  
  for (i in 1:length(day_month)){
    ind <- which(format(predicciones$Date, "%d-%m") == day_month[i])
    dif <- predicciones$Y[ind] - pred[ind]
    df_dia[i,] <- sum(dif < 0, na.rm = T) / (64 * 40)
    
    #por estaciones
    for (j in 1:dim(stations)[1]){
      nombre <- stations$NAME2[j]
      
      # Si la estación aún no está en la lista, inicialízala
      if (!(nombre %in% names(df_dia_list))) {
        df_dia_list[[nombre]] <- data.frame(rho_bay_dia = rep(NA, length(day_month)), row.names = day_month)
      }
      df_temp <- df_dia_list[[nombre]]
      
      ind_2 <- which(predicciones$station == stations$STAID[j])
      ind_2 <- ind_2[which(ind_2 %in% ind)]
      dif <- predicciones$Y[ind_2] - pred[ind_2]
      
      #guardado
      df_temp[i, 1] <- sum(dif < 0, na.rm = T) / 64
      df_dia_list[[nombre]] <- df_temp
      
      
    }
    
    
  }
  
  #años
  df_year_list <- list()#lista para años por estacion
  
  year <- unique(year(predicciones$Date))
  df_year <- matrix(NA, nrow = length(year), ncol = 1)
  df_year <- as.data.frame(df_year, row.names = year)
  colnames(df_year) <- c('rho_bay_year')
  for (i in 1:length(year)){
    ind <- which(year(predicciones$Date) == year[i])
    dif <- predicciones$Y[ind] - pred[ind]
    df_year[i,] <- sum(dif < 0, na.rm = T) / (40 * 92)
    
    for (j in 1:dim(stations)[1]){
      nombre <- stations$NAME2[j]
      
      # Si la estación aún no está en la lista, inicialízala
      if (!(nombre %in% names(df_year_list))) {
        df_year_list[[nombre]] <- data.frame(rho_bay_year = rep(NA, length(year)), row.names = year)
      }
      df_temp <- df_year_list[[nombre]]
      
      ind_2 <- which(predicciones$station == stations$STAID[j])
      ind_2 <- ind_2[which(ind_2%in%ind)]
      dif <- predicciones$Y[ind_2] - pred[ind_2]
      
      #guardado
      df_temp[i, 1] <- sum(dif < 0, na.rm = T) / 92
      df_year_list[[nombre]] <- df_temp
    }
    
  }
  
  return(list(rho_globales=df,
              rho_estaciones=df_est,
              rho_años=df_year,
              rho_dias=df_dia,
              rho_dias_est=df_dia_list,
              rho_años_est=df_year_list))
}

rho_q0.95 <- rho_bay(pred_q0.95, 0.95)
rho_q0.50 <- rho_bay(pred_q0.50, 0.50)


save(df.M1,
     vars,
     stations,
     stations_dist,
     elev_sc, dist_sc,
     local.models.q0.50,
     local.models.q0.95,
     formula,
     mod_q0.50_bay,
     mod_q0.95_bay,
     pred_q0.95,
     pred_q0.50,
     R1_bay_q0.50,
     R1_bay_q0.95,
     rho_q0.50,
     rho_q0.95,
     file = 'dataM1.RData')
