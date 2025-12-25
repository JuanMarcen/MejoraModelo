#fitting of 3 tyoes of model
rm(list = ls())

#libraries and data importation
library(qs)
library(dplyr)
library(lubridate)
library(spTReg)
library(quantreg)
library(sf)
library(sp)
library(coda)
stations <- readRDS('stations.rds')
df <- qread('df_jun_ag.qs')
stations_dist <- readRDS('stations_dist.rds')


#dataframe for each model
df.M1 <- df %>%
  select(Date, station, Y, l, t, s.1, c.1, elev, dist) %>%
  mutate(
    month = month(Date),
    `t:month6` = ifelse(month == 6, t, 0),
    `t:month7` = ifelse(month == 7, t, 0),
    `t:month8` = ifelse(month == 8, t, 0)
  )

df.M2 <- df %>%
  select(Date, station, Y, l, t, g300, g500, g700, elev, dist)

df.M3 <- df %>%
  select(Date, station, Y, l, t, s.1, c.1, elev, dist, g300, g500, g700) %>%
  mutate(
    month = month(Date),
    `t:month6` = ifelse(month == 6, t, 0),
    `t:month7` = ifelse(month == 7, t, 0),
    `t:month8` = ifelse(month == 8, t, 0)
  )

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


# variables and formula for each model
vars.M1 <- c('s.1', 'c.1', '`t:month6`', '`t:month7`', '`t:month8`')
formula.M1 <- as.formula(paste('Y ~', paste(vars.M1, collapse = '+'), '+ elev + dist'))

vars.M2 <- c('g300', 'g500', 'g700')
formula.M2 <- as.formula(paste('Y ~', paste(vars.M2, collapse = '+'), '+ elev + dist'))

vars.M3 <- c(vars.M1, vars.M2)
formula.M3 <- as.formula(paste('Y ~', paste(vars.M3, collapse = '+'), '+ elev + dist'))

#----FITTING OF MODELS (PARALELLIZED TO OBTAIN 3 CHAINS OF PARAMETERS)----
source('metrics_bay.R')

# STARTING POINTS
starting_points <- function(vars, stations, seed = 05052002){
  set.seed(seed)
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
  
  return(inits_list)
}

init.list.M1 <- starting_points(vars.M1, stations)
init.list.M2 <- starting_points(vars.M2, stations)
init.list.M3 <- starting_points(vars.M3, stations)

# paralellized function
library(parallel)
mod.bay.parallel <- function(n.mod, formula, data, vars, coords, init.list, 
                             tau, n.samples, n.burnin, n.thin, n.report){
  
  # number of nodes
  cl <- makeCluster(n.mod)
  
  # Export objects and functions needed
  clusterExport(cl, varlist = c("formula", "data", "vars", "coords", "init.list",
                                "tau", "n.samples", "n.burnin", "n.thin", "n.report",
                                "mod_bay"),
                envir = environment())
  
  # load needed libraries
  clusterEvalQ(cl, library(spTReg))
  
  # Fitting of 3 models in parallel
  results <- parLapply(
    cl,
    X = init.list,
    fun = function(inits){
      mod_bay(
        formula = formula, 
        data = data, 
        tau = tau,                      # Aquí tu tau fijo
        vars = vars, 
        coords = coords, 
        start_beta = inits$start_beta,   # Cambia beta
        inic_procesos = inits$inic_proc, # Cambia procesos
        n.samples = n.samples, 
        n.burnin = n.burnin, 
        n.thin = n.thin, 
        n.report = n.report
      )
    }
  )
  
  stopCluster(cl)
  return(results)
}

final.chain <- function(models.list){
  chains <- lapply(models.list, function(mod) mod$p.params.samples)

  final.chain <- do.call(rbind, chains)
  
  final.chain <- as.data.frame(final.chain)
  
  return(final.chain)
}

elev_sc <- scale(stations_dist$HGHT)
dist_sc <- scale(stations_dist$DIST)
#----q 0.50----
models.q0.50.M1 <- mod.bay.parallel(n.mod = 3,
                              formula = formula.M1,
                              data = df.M1,
                              vars = vars.M1,
                              coords = coords_km,
                              init.list = init.list.M1,
                              tau = 0.50,
                              n.samples = 10,
                              n.burnin = 10,
                              n.thin = 1,
                              n.report = 1)

final.chain.q0.50.M1 <- final.chain(models.q0.50.M1)

models.q0.50.M2 <- mod.bay.parallel(n.mod = 3,
                                    formula = formula.M2,
                                    data = df.M2,
                                    vars = vars.M2,
                                    coords = coords_km,
                                    init.list = init.list.M2,
                                    tau = 0.50,
                                    n.samples = 10,
                                    n.burnin = 10,
                                    n.thin = 1,
                                    n.report = 1)

final.chain.q0.50.M2 <- final.chain(models.q0.50.M2)

models.q0.50.M3 <- mod.bay.parallel(n.mod = 3,
                                    formula = formula.M3,
                                    data = df.M3,
                                    vars = vars.M3,
                                    coords = coords_km,
                                    init.list = init.list.M3,
                                    tau = 0.50,
                                    n.samples = 10,
                                    n.burnin = 10,
                                    n.thin = 1,
                                    n.report = 1)

final.chain.q0.50.M3 <- final.chain(models.q0.50.M3)

#betas, predictions, metrics...
betas.q0.50.M1 <- betas(vars.M1, final.chain.q0.50.M1)
betas.q0.50.M2 <- betas(vars.M2, final.chain.q0.50.M2)
betas.q0.50.M3 <- betas(vars.M3, final.chain.q0.50.M3)

pred.q0.50.M1 <- predictions(vars.M1, betas.q0.50.M1, df.M1, cuantil = 0.50)
pred.q0.50.M1 <- cbind(df.M1[,c('Date', 'station', 'Y')], pred.q0.50.M1)
colnames(pred.q0.50.M1)[4] <- 'pred_q0.50'

pred.q0.50.M2 <- predictions(vars.M2, betas.q0.50.M2, df.M2, cuantil = 0.50)
pred.q0.50.M2 <- cbind(df.M2[,c('Date', 'station', 'Y')], pred.q0.50.M2)
colnames(pred.q0.50.M2)[4] <- 'pred_q0.50'

pred.q0.50.M3 <- predictions(vars.M3, betas.q0.50.M3, df.M3, cuantil = 0.50)
pred.q0.50.M3 <- cbind(df.M3[,c('Date', 'station', 'Y')], pred.q0.50.M3)
colnames(pred.q0.50.M3)[4] <- 'pred_q0.50'
                       
R1.bay.q0.50.M1 <- R1_bay(pred.q0.50.M1, 0.50, df.M1)
R1.bay.q0.50.M2 <- R1_bay(pred.q0.50.M2, 0.50, df.M2)
R1.bay.q0.50.M3 <- R1_bay(pred.q0.50.M3, 0.50, df.M3)

rho.q0.50.M1 <- rho_bay(pred.q0.50.M1, 0.50)
rho.q0.50.M2 <- rho_bay(pred.q0.50.M2, 0.50)
rho.q0.50.M3 <- rho_bay(pred.q0.50.M3, 0.50)

#----q 0.95----
models.q0.95.M1 <- mod.bay.parallel(n.mod = 3,
                                    formula = formula.M1,
                                    data = df.M1,
                                    vars = vars.M1,
                                    coords = coords_km,
                                    init.list = init.list.M1,
                                    tau = 0.95,
                                    n.samples = 100,
                                    n.burnin = 100,
                                    n.thin = 1,
                                    n.report = 1)

final.chain.q0.95.M1 <- final.chain(models.q0.95.M1)

models.q0.95.M2 <- mod.bay.parallel(n.mod = 3,
                                    formula = formula.M2,
                                    data = df.M2,
                                    vars = vars.M2,
                                    coords = coords_km,
                                    init.list = init.list.M2,
                                    tau = 0.95,
                                    n.samples = 100,
                                    n.burnin = 100,
                                    n.thin = 1,
                                    n.report = 1)

final.chain.q0.95.M2 <- final.chain(models.q0.95.M2)

models.q0.95.M3 <- mod.bay.parallel(n.mod = 3,
                                    formula = formula.M3,
                                    data = df.M3,
                                    vars = vars.M3,
                                    coords = coords_km,
                                    init.list = init.list.M3,
                                    tau = 0.95,
                                    n.samples = 100,
                                    n.burnin = 100,
                                    n.thin = 1,
                                    n.report = 1)

final.chain.q0.95.M3 <- final.chain(models.q0.95.M3)

#betas, predictions, metrics...
betas.q0.95.M1 <- betas(vars.M1, final.chain.q0.95.M1)
betas.q0.95.M2 <- betas(vars.M2, final.chain.q0.95.M2)
betas.q0.95.M3 <- betas(vars.M3, final.chain.q0.95.M3)

pred.q0.95.M1 <- predictions(vars.M1, betas.q0.95.M1, df.M1, cuantil = 0.95)
pred.q0.95.M1 <- cbind(df.M1[,c('Date', 'station', 'Y')], pred.q0.95.M1)

pred.q0.95.M2 <- predictions(vars.M2, betas.q0.95.M2, df.M2, cuantil = 0.95)
pred.q0.95.M2 <- cbind(df.M2[,c('Date', 'station', 'Y')], pred.q0.95.M2)

pred.q0.95.M3 <- predictions(vars.M3, betas.q0.95.M3, df.M3, cuantil = 0.95)
pred.q0.95.M3 <- cbind(df.M3[,c('Date', 'station', 'Y')], pred.q0.95.M3)

R1.bay.q0.95.M1 <- R1_bay(pred.q0.95.M1, 0.95, df.M1)
R1.bay.q0.95.M2 <- R1_bay(pred.q0.95.M2, 0.95, df.M2)
R1.bay.q0.95.M3 <- R1_bay(pred.q0.95.M3, 0.95, df.M3)

rho.q0.95.M1 <- rho_bay(pred.q0.95.M1, 0.95)
rho.q0.95.M2 <- rho_bay(pred.q0.95.M2, 0.95)
rho.q0.95.M3 <- rho_bay(pred.q0.95.M3, 0.95)


save(models.q0.50.M1, models.q0.50.M2, models.q0.50.M3,
     betas.q0.50.M1, betas.q0.50.M2, betas.q0.50.M3,
     final.chain.q0.50.M1, final.chain.q0.50.M2, final.chain.q0.50.M3,
     pred.q0.50.M1, pred.q0.50.M2, pred.q0.50.M3,
     R1.bay.q0.50.M1, R1.bay.q0.50.M2, R1.bay.q0.50.M3,
     rho.q0.50.M1, rho.q0.50.M2, rho.q0.50.M3,
     models.q0.95.M1, models.q0.95.M2, models.q0.95.M3,
     betas.q0.95.M1, betas.q0.95.M2, betas.q0.95.M3,
     final.chain.q0.95.M1, final.chain.q0.95.M2, final.chain.q0.95.M3,
     pred.q0.95.M1, pred.q0.95.M2, pred.q0.95.M3,
     R1.bay.q0.95.M1, R1.bay.q0.95.M2, R1.bay.q0.95.M3,
     rho.q0.95.M1, rho.q0.95.M2, rho.q0.95.M3,
     file = 'fullmodels.RData')

#----CONVERGENCE----
#shrink factor and traceplots
# modify to see if all are less that 1.1??
load('fullmodels.RData')

png('graphs/convergence/conv.q0.50.M3.int.png',height = 2000,width = 4000,res=150)
par(mfrow=c(5,8))
shrink.f.tr.plot(models.q0.50.M2, type = 'coef', vars = vars.M2)
dev.off()

# ESS
ESS(final.chain.q0.50.M3, type = 'coef', vars = vars.M3)

# uncertainty
uncertainty(final.chain.q0.50.M1, 'coef', vars.M1)



# extra 
plot(1:64, unlist(rho.q0.50.M3$rho_años), type = 'l',
     ylim = c(0,1), lwd = 2)
abline(h = 0.5, col = 'red')
lines(1:64, unlist(rho.q0.50.M1$rho_años), col = 'blue')
lines(1:64, unlist(rho.q0.50.M2$rho_años), col = 'forestgreen')

plot(1:92, unlist(rho.q0.50.M3$rho_dias), type = 'l',
     ylim = c(0,1), lwd = 2)
abline(h = 0.5, col = 'red')
lines(1:92, unlist(rho.q0.50.M1$rho_dias), col = 'blue')
lines(1:92, unlist(rho.q0.50.M2$rho_dias), col = 'forestgreen')

#----Writing tables----
orden <- readRDS('orden.rds')
R1.df <- data.frame(
  M1.q0.50 = round(R1.bay.q0.50.M1$R1_locales, 3),
  M2.q0.50 = round(R1.bay.q0.50.M2$R1_locales, 3),
  M3.q0.50 = round(R1.bay.q0.50.M3$R1_locales, 3),
  M1.q0.95 = round(R1.bay.q0.95.M1$R1_locales, 3),
  M2.q0.95 = round(R1.bay.q0.95.M2$R1_locales, 3),
  M3.q0.95 = round(R1.bay.q0.95.M3$R1_locales, 3)
)
rownames(R1.df) <- stations$NAME2
R1.df <- R1.df[orden, ]

esc_tabla_negrita<-function(tabla,colq0.5,colq0.95,negrita=T){
  
  fila <- rownames(tabla)
  columnas <- ncol(tabla)
  
  media_cant<-round(apply(tabla[1:6,],FUN = mean,MARGIN = 2),3)
  media_med<-round(apply(tabla[7:18,],FUN = mean,MARGIN = 2),3)
  media_centro<-round(apply(tabla[19:40,],FUN = mean,MARGIN = 2),3)
  media_todo<-round(apply(tabla,FUN = mean,MARGIN = 2),3)
  
  for (i in 1:dim(tabla)[1]){
    
    maxq0.5<-max(tabla[i,colq0.5])
    maxq0.95<-max(tabla[i,colq0.95])
    
    cat(fila[i],'& ')
    
    
    for (j in 1:columnas){
      
      x <- tabla[i,j]
      
      if (j %in% colq0.5 && x == maxq0.5 && negrita==T){
        cat('$\\mathbf{',format(x,nsmall=3),'}$')
      }else if(j %in% colq0.95 && x == maxq0.95 && negrita==T){
        cat('$\\mathbf{',format(x,nsmall=3),'}$')
      }else{
        cat('$',format(x,nsmall=3),'$')
      }
      
      
      
      #separo valor y salto linea
      if (j < columnas){
        cat(' & ')
      }else{
        cat(' \\\\\n')
      }
      
    }
    
    #separacion zonas
    if(i %in% c(6,18)){
      cat('[2pt]','\\hline', '\\\\ [-10pt]')
    }
    if(i==40){
      cat('[2pt]','\\hline \\\\ [-10pt]')
    }
    
    
    if (i==6){
      cat(' Media & ')
      for (j in 1:columnas){
        m<-media_cant[j]
        cat('$',format(m,nsmall=3),'$')
        if (j < columnas){
          cat(' & ')
        }else{
          cat(' \\\\ [2pt] \\hline \\\\[-10pt] \n')
        }
      }
    }
    
    if (i==18){
      cat(' Media & ')
      for (j in 1:columnas){
        m<-media_med[j]
        cat('$',format(m,nsmall=3),'$')
        if (j < columnas){
          cat(' & ')
        }else{
          cat(' \\\\ [2pt] \\hline \\\\[-10pt] \n')
        }
      }
    }
    
    if (i==40){
      cat(' Media & ')
      for (j in 1:columnas){
        m<-media_centro[j]
        cat('$',format(m,nsmall=3),'$')
        if (j < columnas){
          cat(' & ')
        }else{
          cat(' \\\\ [2pt] \\hline \\\\[-10pt] \n')
        }
      }
    }
  }
  
  cat(' Media total & ')
  for (j in 1:columnas){
    m<-media_todo[j]
    cat('$',format(m,nsmall=3),'$')
    if (j<columnas){
      cat(' & ')
    }else{
      cat(' \\\\ [2pt] \\hline')
    }
  }
  
}

esc_tabla_negrita(R1.df, colq0.5 = c(1, 2, 3), colq0.95 = c(4, 5, 6), negrita = T)

#----plots rho----
pdf('rho_3est.pdf', width = 12, height = 8)
rho.q0.95.M1.madrid <- rho.q0.95.M1$rho_dias_est$`Madrid (Barajas)`
rho.q0.95.M2.madrid <- rho.q0.95.M2$rho_dias_est$`Madrid (Barajas)`
rho.q0.95.M3.madrid <- rho.q0.95.M3$rho_dias_est$`Madrid (Barajas)`

par(mfrow = c(2,3))

plot(1:92, rho.q0.95.M1.madrid$rho_bay_dia, type = 'l',
     ylim = c(0.75,1), col = 'forestgreen',
     ylab = expression(rho[l](0.95 * '; ' * s)),
     xlab = 'l',
     main = 'Madrid (Barajas)')
lines(1:92, rho.q0.95.M2.madrid$rho_bay_dia, col = 'blue')
lines(1:92, rho.q0.95.M3.madrid$rho_bay_dia, lw = 2)
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2))

rho.q0.95.M1.madrid <- rho.q0.95.M1$rho_dias_est$Gijón
rho.q0.95.M2.madrid <- rho.q0.95.M2$rho_dias_est$Gijón
rho.q0.95.M3.madrid <- rho.q0.95.M3$rho_dias_est$Gijón

plot(1:92, rho.q0.95.M1.madrid$rho_bay_dia, type = 'l',
     ylim = c(0.75,1), col = 'forestgreen',
     ylab = expression(rho[l](0.95 * '; ' * s)),
     xlab = 'l',
     main = 'Gijón')
lines(1:92, rho.q0.95.M2.madrid$rho_bay_dia, col = 'blue')
lines(1:92, rho.q0.95.M3.madrid$rho_bay_dia, lw = 2)
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2))

rho.q0.95.M1.madrid <- rho.q0.95.M1$rho_dias_est$`Barcelona (Aeropuerto)`
rho.q0.95.M2.madrid <- rho.q0.95.M2$rho_dias_est$`Barcelona (Aeropuerto)`
rho.q0.95.M3.madrid <- rho.q0.95.M3$rho_dias_est$`Barcelona (Aeropuerto)`

plot(1:92, rho.q0.95.M1.madrid$rho_bay_dia, type = 'l',
     ylim = c(0.75,1), col = 'forestgreen',
     ylab = expression(rho[l](0.95 * '; ' * s)),
     xlab = 'l',
     main = 'Barcelona (Airport)')
lines(1:92, rho.q0.95.M2.madrid$rho_bay_dia, col = 'blue')
lines(1:92, rho.q0.95.M3.madrid$rho_bay_dia, lw = 2)
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2))

rho.q0.95.M1.madrid <- rho.q0.95.M1$rho_años_est$`Madrid (Barajas)`
rho.q0.95.M2.madrid <- rho.q0.95.M2$rho_años_est$`Madrid (Barajas)`
rho.q0.95.M3.madrid <- rho.q0.95.M3$rho_años_est$`Madrid (Barajas)`

plot(1:64, rho.q0.95.M1.madrid$rho_bay_year, type = 'l',
     ylim = c(0.75,1), col = 'forestgreen',
     ylab = expression(rho[t](0.95 * '; ' * s)),
     xlab = 't')
lines(1:64, rho.q0.95.M2.madrid$rho_bay_year, col = 'blue')
lines(1:64, rho.q0.95.M3.madrid$rho_bay_year, lw = 2)
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2))

rho.q0.95.M1.madrid <- rho.q0.95.M1$rho_años_est$Gijón
rho.q0.95.M2.madrid <- rho.q0.95.M2$rho_años_est$Gijón
rho.q0.95.M3.madrid <- rho.q0.95.M3$rho_años_est$Gijón

plot(1:64, rho.q0.95.M1.madrid$rho_bay_year, type = 'l',
     ylim = c(0.75,1), col = 'forestgreen',
     ylab = expression(rho[t](0.95 * '; ' * s)),
     xlab = 't')
lines(1:64, rho.q0.95.M2.madrid$rho_bay_year, col = 'blue')
lines(1:64, rho.q0.95.M3.madrid$rho_bay_year, lw = 2)
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2))

rho.q0.95.M1.madrid <- rho.q0.95.M1$rho_años_est$`Barcelona (Aeropuerto)`
rho.q0.95.M2.madrid <- rho.q0.95.M2$rho_años_est$`Barcelona (Aeropuerto)`
rho.q0.95.M3.madrid <- rho.q0.95.M3$rho_años_est$`Barcelona (Aeropuerto)`

plot(1:64, rho.q0.95.M1.madrid$rho_bay_year, type = 'l',
     ylim = c(0.75,1), col = 'forestgreen',
     ylab = expression(rho[t](0.95 * '; ' * s)),
     xlab = 't')
lines(1:64, rho.q0.95.M2.madrid$rho_bay_year, col = 'blue')
lines(1:64, rho.q0.95.M3.madrid$rho_bay_year, lw = 2)
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2))
dev.off()

pdf('rho_global_q0.50.pdf', width = 12, height = 8*2/3)
a <- rho.q0.50.M1$rho_dias
b <- rho.q0.50.M2$rho_dias
c <- rho.q0.50.M3$rho_dias

par(mfrow = c(1,2))

plot(1:92, a$rho_bay_dia, type = 'l',
     ylim = c(0.25,0.75), col = 'forestgreen',
     ylab = expression(rho[l](0.50)),
     xlab = 'l')
lines(1:92, b$rho_bay_dia, col = 'blue')
lines(1:92, c$rho_bay_dia, lw = 2)
abline(h = 0.50, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2),
       bg = 'white')

a <- rho.q0.50.M1$rho_años
b <- rho.q0.50.M2$rho_años
c <- rho.q0.50.M3$rho_años

plot(1:64, a$rho_bay_year, type = 'l',
     ylim = c(0.25,0.75), col = 'forestgreen',
     ylab = expression(rho[t](0.50)),
     xlab = 't')
lines(1:64, b$rho_bay_year, col = 'blue')
lines(1:64, c$rho_bay_year, lw = 2)
abline(h = 0.50, col = 'red')
legend('bottom', legend = c('M1', 'M2', 'M3'),
       ncol = 3,
       col = c('forestgreen', 'blue', 'black'), lty = 1, lwd = c(1,1,2),
       bg = 'white')

dev.off()
#----CI betas----
final.chain <- final.chain.q0.50.M1
sta <- c('Madrid (Barajas)', 'Gijón', 'Barcelona (Aeropuerto)')
ind.s <- which(stations$NAME2 %in% sta)

df <- data.frame(matrix(NA,nrow = nrow(final.chain)))
#M1
for (i in ind.s){
  for (j in 1:length(vars.M1)){
    var <- vars.M1[j]
    df[[paste0(var, '.M1.', stations$NAME2[i])]] <- final.chain[[var]] + 
      final.chain[[paste0('beta',j,'(s',i,')')]]
  }
}

#M2
final.chain <- final.chain.q0.50.M2
for (i in ind.s){
  for (j in 1:length(vars.M2)){
    var <- vars.M2[j]
    df[[paste0(var, '.M2.', stations$NAME2[i])]] <- final.chain[[var]] + 
      final.chain[[paste0('beta',j,'(s',i,')')]]
  }
}

#M3
final.chain <- final.chain.q0.50.M3
for (i in ind.s){
  for (j in 1:length(vars.M3)){
    var <- vars.M3[j]
    df[[paste0(var, '.M3.', stations$NAME2[i])]] <- final.chain[[var]] + 
      final.chain[[paste0('beta',j,'(s',i,')')]]
  }
}

df <- df[, -1]

CI.coef <- apply(df, 2, quantile, probs = c(0.025, 0.0975))
