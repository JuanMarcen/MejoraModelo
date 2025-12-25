rm(list = setdiff(ls(), c('vars.M1', 'final.chain.q0.50.M1', 'df.M1',
                          'vars.M2', 'final.chain.q0.50.M2', 'df.M2',
                          'vars.M3', 'final.chain.q0.50.M3', 'df.M3',
                          'final.chain.q0.95.M1', 
                          'final.chain.q0.95.M2', 
                          'final.chain.q0.95.M3', 
                          'coords_km', 'stations', 'stations_dist'
                          )))

load('spatial.interplation.RData')
load('fullmodels.RData')
# mantengo 
# vars, final.chain, coords, stations_dist, stations, df,
library(stringr)
library(coda)
library(spTReg)
library(dplyr)
library(lubridate)
grid <- readRDS('grid.rds')
grid_km <- readRDS('grid_km.rds')
grid_elev <- readRDS('grid_elev.rds')
grid_dist <- readRDS('grid_dist.rds')

source('metrics_bay.R')

#----BAYESIAN KRIGING----
bay.kriging <- function(final.chain, vars, coords, grid, stations_dist){
  #number of spatial coeff
  ind <- length(vars) + 1
  vars <- gsub('`', '', vars)
  # bayesian kriging of rest of grid
  final.chain.mcmc <- as.mcmc(final.chain)
  for (i in 1:ind) {
    v<-krigeBayes(final.chain.mcmc[, paste0('beta', i, '(s', 1:40, ')')],
                  hp = final.chain.mcmc[, paste0(c('mu', 'sigma', 'phi'), i)],
                  coords = coords, newcoords = grid)
    assign(paste0("b", i), v)
  }
  
  #Objective: b_i = b_i + b_i(s)
  
  #beta1 to betai are ordered by vars
  params <- as.data.frame(final.chain)
  # in case of polys
  params_c <- params
  colnames(params_c)[2:ind] <- gsub('`', '', colnames(params_c)[2:ind])
  tr <- traducir_nombres_coef(colnames(params)[2:ind])
  tr <- gsub('`', '', tr) #quito comas
  colnames(params)[2:ind] <- tr
  betas.fixed <- params[, c('(Intercept)', vars)]
  
  #dataframe con nombres poly para saber localizarlos
  names <- colnames(params_c)[match(vars, tr) + 1]
  
  #All ordered now. betas.fixed[,1] corresponds to predicted b1 (spatial)
  
  # NOT NECESSARY NOW, I STILL DON'T HAVE POLYS
  #ahora tegno que ver, que variables estan con poly
  variables_poly <- str_match(names, "poly\\(([^,]+),")[,2]
  variables_poly <- variables_poly[!is.na(variables_poly)]
  variables_poly <- sort(unique(variables_poly))
  # y sacar sus coeficientes del polinomio ortogonal (a,b,c,d,e)
  coef_orth <- data.frame(matrix(
    NA, ncol = length(variables_poly), nrow = 5),
    row.names = c('a','b','c','d','e'))
  colnames(coef_orth) <- variables_poly
  for (var in variables_poly){
    p <- poly.orth(df[,var],2)
    coef_orth['a',var] <- coef(p[[2]])[1]
    coef_orth['b',var] <- coef(p[[2]])[2]
    coef_orth['c',var] <- coef(p[[3]])[1]
    coef_orth['d',var] <- coef(p[[3]])[2]
    coef_orth['e',var] <- coef(p[[3]])[3]
  }
  
  coef_orth <- cbind(coef_orth, coef_orth) #duplicado para luego bucle entero mejor
  if (ncol(coef_orth) > 0){
    colnames(coef_orth)[(length(variables_poly) + 1):ncol(coef_orth)] <- paste0(
      'I(', colnames(coef_orth)[1:length(variables_poly)], '^2)')
  }
  
  
  
  # NOT NECESSARY. I WORK WITH ANOMALIES, NOT SCALED COVARIATES
  #medias y desviaciones típicas ORDENADAS NO SE INCLUYE EL INTERCEPTO. 
  # Tambien se trata el caso en que solo se de una variable
  # if (length(vars)==1){
  #   medias<-mean(df[,vars])
  #   desv<-sd(df[,vars])
  # }else{
  #   # medias<-colMeans(df[,vars])
  #   # desv<-apply(df[,vars],2,sd)
  #   medias<-colMeans(df_medias)
  #   desv<-apply(df_medias,2,sd)
  # }
  
  # calculate how the elevation and distance were scaled, for the new elevations and distances
  elev_mu <- mean(stations_dist$HGHT)
  dist_mu <- mean(stations_dist$DIST)
  elev_sd <- sd(stations_dist$HGHT)
  dist_sd <- sd(stations_dist$DIST)
  
  #spatial betas. When haivng polys, look at the TFM code
  for (i in 1:ind){
    
    if (i == 1){ #intercept
      print('INTERCEPTO')
      res <- betas.fixed[, 1] + get("b1") #intercept
      # add elev and dist
      res <- res - final.chain[, 'elev'] * elev_mu / elev_sd - final.chain[, 'dist'] * dist_mu / dist_sd
      cont <- 1
      for (j in 2:ind) {
        if (vars[j - 1] %in% colnames(coef_orth)){ #poly
          if (cont <= length(variables_poly)){ #polys lineals
            print('poly')
            res <- res + (betas.fixed[, j] + get(paste0("b", j))) * (coef_orth['a',vars[j-1]])
            cont <-cont + 1
          }else{
            print('poly_sq')
            res <- res + (betas.fixed[, j] + get(paste0("b", j))) * (coef_orth['c', vars[j - 1]])
            cont <- cont +1
          }
        }else{ #normales
          print('normal')
          res <- res
        }
      }
      
      assign(paste0('b',1),res)
      
      
      
    }else{ #rest of covariates
      print('REST OF COVARIATES')
      if(vars[i-1] %in% colnames(coef_orth)){#poly
        # distincion en si quiero coeficinete de sgundo grado o lineal
        if(grepl('I',vars[i-1])==F){#lineal
          print('poly lineal')
          var_aux<-paste0('I(',vars[i-1],'^2)')
          index<-which(var_aux==vars)
          beta<-(betas.fixed[,i]+get(paste0('b',i)))*coef_orth['b',vars[i-1]]/desv[i-1] + (betas.fixed[,index+1]+get(paste0('b',index+1)))*coef_orth['d',var_aux]/desv[index]
          assign(paste0('b',i), beta )
        }else{#cuadratico
          print('poly sq')
          beta<-(betas.fixed[,i]+get(paste0('b',i)))*coef_orth['e',vars[i-1]]/desv[i-1]
          assign(paste0('b',i), beta )
        }
      }else{
        #print('elementos no poly')
        assign(paste0('b', i), (betas.fixed[, i]+get(paste0('b', i))))
      }
      
    }
    
    
    
  }
  
  return(
    setNames(
      mget(paste0("b", 1:ind)), 
      paste0("beta", 1:ind)
    )
  )
}

bay.kriging.q0.50.M1 <- bay.kriging(final.chain.q0.50.M1, vars.M1, 
                                    coords_km, grid_km, stations_dist)
bay.kriging.q0.50.M2 <- bay.kriging(final.chain.q0.50.M2, vars.M2, 
                                    coords_km, grid_km, stations_dist)
bay.kriging.q0.50.M3 <- bay.kriging(final.chain.q0.50.M3, vars.M3, 
                                    coords_km, grid_km, stations_dist)

bay.kriging.q0.95.M1 <- bay.kriging(final.chain.q0.95.M1, vars.M1, 
                                    coords_km, grid_km, stations_dist)
bay.kriging.q0.95.M2 <- bay.kriging(final.chain.q0.95.M2, vars.M2, 
                                    coords_km, grid_km, stations_dist)
bay.kriging.q0.95.M3 <- bay.kriging(final.chain.q0.95.M3, vars.M3, 
                                    coords_km, grid_km, stations_dist)

#----QUANTILE PREDICTION----
library(qs)
X.grid <- qread('X_grid.qs')


#quantile prediction using means
quantile.pred <- function(X.grid, data.model,
                          vars, final.chain, betas.bay.kriging){
  
  pred.df <- matrix(NA, nrow = length(unique(data.model$Date)), ncol = 790)
  pred.df <- as.data.frame(pred.df)
  
  vars <- gsub('`', '', vars)
  
  # betas elevation and distance
  elev <- mean(final.chain[, 'elev'])
  dist <- mean(final.chain[, 'dist'])
 
  elev_sd <- sd(stations_dist$HGHT)
  dist_sd <- sd(stations_dist$DIST)
  
  for (i in 1:nrow(pred.df)){
    date <- unique(data.model$Date)[i]
    ind <- which(X.grid$Date == date)
    vars.grid <- X.grid[ind, vars]
    # Extraer las variables cuadráticas (entre I(...) )
    # cuad_vars <- gsub("^I\\((.*)\\^2\\)$", "\\1", vars[grepl("^I\\(.*\\^2\\)$", vars)])
    # # Crear las columnas cuadradas a partir de G
    # cuad_df <- as.data.frame(sapply(cuad_vars, function(var) G[[var]]^2))
    # names(cuad_df) <- paste0("I(", cuad_vars, "^2)")
    # # Extraer las variables normales
    norm_vars <- vars[!grepl("^I\\(.*\\^2\\)$", vars)]
    norm_df <- vars.grid[, norm_vars, drop = FALSE]
    # # Combinar en orden según vars
    # df_final <- cbind(norm_df, cuad_df)[, vars]
    
    df_final <- norm_df[, vars]
    
    pred <- matrix(colMeans(betas.bay.kriging[['beta1']]), nrow=1)
    
    for (j in 1:length(vars)){
      pred <- pred + matrix(colMeans(betas.bay.kriging[[paste0('beta', j + 1)]], na.rm = T),nrow=1) * t(df_final[, vars[j]])
    }
    
    #add distance and elevation
    pred.df[i, ] <- pred + elev * grid_elev / elev_sd + dist * grid_dist/dist_sd

    print(paste(date, pred.df[i, 1]))
  }
  
  rownames(pred.df) <- unique(data.model$Date)
  
  return(pred.df)
}


quantile.q0.50.M1 <- quantile.pred(X.grid = X.grid, data.model = df.M1, vars = vars.M1,
                                   final.chain = final.chain.q0.50.M1,
                                   betas.bay.kriging = bay.kriging.q0.50.M1)
quantile.q0.50.M2 <- quantile.pred(X.grid = X.grid, data.model = df.M2, vars = vars.M2,
                                   final.chain = final.chain.q0.50.M2,
                                   betas.bay.kriging = bay.kriging.q0.50.M2)
quantile.q0.50.M3 <- quantile.pred(X.grid = X.grid, data.model = df.M3, vars = vars.M3,
                                   final.chain = final.chain.q0.50.M3,
                                   betas.bay.kriging = bay.kriging.q0.50.M3)

quantile.q0.95.M1 <- quantile.pred(X.grid = X.grid, data.model = df.M1, vars = vars.M1,
                                   final.chain = final.chain.q0.95.M1,
                                   betas.bay.kriging = bay.kriging.q0.95.M1)
quantile.q0.95.M2 <- quantile.pred(X.grid = X.grid, data.model = df.M2, vars = vars.M2,
                                   final.chain = final.chain.q0.95.M2,
                                   betas.bay.kriging = bay.kriging.q0.95.M2)
quantile.q0.95.M3 <- quantile.pred(X.grid = X.grid, data.model = df.M3, vars = vars.M3,
                                   final.chain = final.chain.q0.95.M3,
                                   betas.bay.kriging = bay.kriging.q0.95.M3)

#full quantile prediction (using whole chain)
quantile.pred.total <- function(X.grid, data.model, vars, 
                                final.chain, betas.bay.kriging, 
                                year1, year2, month) {
  

  elev <- final.chain[, 'elev']
  dist <- final.chain[, 'dist']
  
  elev_sd <- sd(stations_dist$HGHT)
  dist_sd <- sd(stations_dist$DIST)
  
  vars <- gsub('`', '', vars)
  
  dates <- unique(data.model$Date[year(data.model$Date) >= year1 &
                            year(data.model$Date) <= year2 &
                            month(data.model$Date) == month])
  
  # list for predictions storage
  pred_list <- vector("list", length(dates))
  
  # cuad_vars <- gsub("^I\\((.*)\\^2\\)$", "\\1", vars[grepl("^I\\(.*\\^2\\)$", vars)])
  norm_vars <- vars[!grepl("^I\\(.*\\^2\\)$", vars)]
  
  for (i in seq_along(dates)) {
    date <- dates[i]
    ind <- which(X.grid$Date == date)
    vars.grid <- X.grid[ind, vars, drop = FALSE]
    
    # Crea data frame de variables cuadradas
    # cuad_df <- as.data.frame(lapply(cuad_vars, function(var) vars.grid[[var]]^2))
    # names(cuad_df) <- paste0("I(", cuad_vars, "^2)")
    
    # Combina y ordena columnas
    norm_df <- vars.grid[, norm_vars, drop = FALSE]
    # df_final <- cbind(norm_df, cuad_df)[, vars, drop = FALSE]
    df_final <- cbind(norm_df)[, vars, drop = FALSE]
    
    # Initialize pred with intercept
    pred <- betas.bay.kriging[['beta1']]
    
    # Add effect of variables
    for (j in seq_along(vars)) {
      pred <- pred + sweep(betas.bay.kriging[[paste0('beta', j + 1)]], 2, t(df_final[[vars[j]]]), FUN = '*')
    }
    
    # add terms of distance and elevation
    elev.matrix <- outer(elev, grid_elev, "*") / elev_sd
    dist.matrix <- outer(dist, grid_dist, "*") / dist_sd
    
    pred <- pred + elev.matrix + dist.matrix
    
    pred_list[[i]] <- pred
    if (i %% 10 == 0 || i == 1) cat("Progreso:", i, "/", length(dates), "\n")
  }
  
  # join all to calculate uncertainty
  pred.all <- do.call(rbind, pred_list)
  return(pred.all)
}

# quantile 0.50
# first decade 
#M1
quantile.q0.50.M1.1dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.50.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M1,
                                                  year1 = 1960, year2 = 1969, month = '6')
quantile.q0.50.M1.1dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.50.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M1,
                                                  year1 = 1960, year2 = 1969, month = '7')
quantile.q0.50.M1.1dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.50.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M1,
                                                  year1 = 1960, year2 = 1969, month = '8')
#M2
quantile.q0.50.M2.1dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.50.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M2,
                                                  year1 = 1960, year2 = 1969, month = '6')
quantile.q0.50.M2.1dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.50.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M2,
                                                  year1 = 1960, year2 = 1969, month = '7')
quantile.q0.50.M2.1dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.50.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M2,
                                                  year1 = 1960, year2 = 1969, month = '8')

#M3
quantile.q0.50.M3.1dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                vars = vars.M3,
                                                final.chain = final.chain.q0.50.M3, 
                                                betas.bay.kriging = bay.kriging.q0.50.M3,
                                                year1 = 1960, year2 = 1969, month = '6')
quantile.q0.50.M3.1dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.50.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M3,
                                                  year1 = 1960, year2 = 1969, month = '7')
quantile.q0.50.M3.1dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.50.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M3,
                                                  year1 = 1960, year2 = 1969, month = '8')

# SECOND decade 
#M1
quantile.q0.50.M1.2dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.50.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M1,
                                                  year1 = 2014, year2 = 2023, month = '6')
quantile.q0.50.M1.2dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.50.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M1,
                                                  year1 = 2014, year2 = 2023, month = '7')
quantile.q0.50.M1.2dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.50.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M1,
                                                  year1 = 2014, year2 = 2023, month = '8')
#M2
quantile.q0.50.M2.2dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.50.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M2,
                                                  year1 = 2014, year2 = 2023, month = '6')
quantile.q0.50.M2.2dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.50.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M2,
                                                  year1 = 2014, year2 = 2023, month = '7')
quantile.q0.50.M2.2dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.50.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M2,
                                                  year1 = 2014, year2 = 2023, month = '8')

#M3
quantile.q0.50.M3.2dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.50.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M3,
                                                  year1 = 2014, year2 = 2023, month = '6')
quantile.q0.50.M3.2dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.50.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M3,
                                                  year1 = 2014, year2 = 2023, month = '7')
quantile.q0.50.M3.2dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.50.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.50.M3,
                                                  year1 = 2014, year2 = 2023, month = '8')




# quantile 0.95 
# first decade 
#M1
quantile.q0.95.M1.1dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.95.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M1,
                                                  year1 = 1960, year2 = 1969, month = '6')
quantile.q0.95.M1.1dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.95.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M1,
                                                  year1 = 1960, year2 = 1969, month = '7')
quantile.q0.95.M1.1dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.95.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M1,
                                                  year1 = 1960, year2 = 1969, month = '8')
#M2
quantile.q0.95.M2.1dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.95.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M2,
                                                  year1 = 1960, year2 = 1969, month = '6')
quantile.q0.95.M2.1dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.95.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M2,
                                                  year1 = 1960, year2 = 1969, month = '7')
quantile.q0.95.M2.1dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.95.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M2,
                                                  year1 = 1960, year2 = 1969, month = '8')

#M3
quantile.q0.95.M3.1dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.95.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M3,
                                                  year1 = 1960, year2 = 1969, month = '6')
quantile.q0.95.M3.1dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.95.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M3,
                                                  year1 = 1960, year2 = 1969, month = '7')
quantile.q0.95.M3.1dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.95.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M3,
                                                  year1 = 1960, year2 = 1969, month = '8')

# SECOND decade 
#M1
quantile.q0.95.M1.2dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.95.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M1,
                                                  year1 = 2014, year2 = 2023, month = '6')
quantile.q0.95.M1.2dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.95.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M1,
                                                  year1 = 2014, year2 = 2023, month = '7')
quantile.q0.95.M1.2dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M1, 
                                                  vars = vars.M1,
                                                  final.chain = final.chain.q0.95.M1, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M1,
                                                  year1 = 2014, year2 = 2023, month = '8')
#M2
quantile.q0.95.M2.2dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.95.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M2,
                                                  year1 = 2014, year2 = 2023, month = '6')
quantile.q0.95.M2.2dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.95.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M2,
                                                  year1 = 2014, year2 = 2023, month = '7')
quantile.q0.95.M2.2dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M2, 
                                                  vars = vars.M2,
                                                  final.chain = final.chain.q0.95.M2, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M2,
                                                  year1 = 2014, year2 = 2023, month = '8')

#M3
quantile.q0.95.M3.2dec.jun <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.95.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M3,
                                                  year1 = 2014, year2 = 2023, month = '6')
quantile.q0.95.M3.2dec.jul <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.95.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M3,
                                                  year1 = 2014, year2 = 2023, month = '7')
quantile.q0.95.M3.2dec.aug <- quantile.pred.total(X.grid = X.grid, data.model = df.M3, 
                                                  vars = vars.M3,
                                                  final.chain = final.chain.q0.95.M3, 
                                                  betas.bay.kriging = bay.kriging.q0.95.M3,
                                                  year1 = 2014, year2 = 2023, month = '8')


#----DIFFERENCES AND SEE IF 0 IN CREDIBILITY INTERVAL----
dif.q0.50.M1.jun <- quantile.q0.50.M1.2dec.jun - quantile.q0.50.M1.1dec.jun
dif.q0.50.M1.jul <- quantile.q0.50.M1.2dec.jul - quantile.q0.50.M1.1dec.jul
dif.q0.50.M1.aug <- quantile.q0.50.M1.2dec.aug - quantile.q0.50.M1.1dec.aug
rm('quantile.q0.50.M1.2dec.jun', 'quantile.q0.50.M1.1dec.jun',
   'quantile.q0.50.M1.2dec.jul', 'quantile.q0.50.M1.1dec.jul',
   'quantile.q0.50.M1.2dec.aug', 'quantile.q0.50.M1.1dec.aug')
gc()

dif.q0.50.M2.jun <- quantile.q0.50.M2.2dec.jun - quantile.q0.50.M2.1dec.jun
dif.q0.50.M2.jul <- quantile.q0.50.M2.2dec.jul - quantile.q0.50.M2.1dec.jul
dif.q0.50.M2.aug <- quantile.q0.50.M2.2dec.aug - quantile.q0.50.M2.1dec.aug
rm('quantile.q0.50.M2.2dec.jun', 'quantile.q0.50.M2.1dec.jun',
   'quantile.q0.50.M2.2dec.jul', 'quantile.q0.50.M2.1dec.jul',
   'quantile.q0.50.M2.2dec.aug', 'quantile.q0.50.M2.1dec.aug')
gc()


dif.q0.50.M3.jun <- quantile.q0.50.M3.2dec.jun - quantile.q0.50.M3.1dec.jun
dif.q0.50.M3.jul <- quantile.q0.50.M3.2dec.jul - quantile.q0.50.M3.1dec.jul
dif.q0.50.M3.aug <- quantile.q0.50.M3.2dec.aug - quantile.q0.50.M3.1dec.aug
rm('quantile.q0.50.M3.2dec.jun', 'quantile.q0.50.M3.1dec.jun',
   'quantile.q0.50.M3.2dec.jul', 'quantile.q0.50.M3.1dec.jul',
   'quantile.q0.50.M3.2dec.aug', 'quantile.q0.50.M3.1dec.aug')
gc()


dif.q0.95.M1.jun <- quantile.q0.95.M1.2dec.jun - quantile.q0.95.M1.1dec.jun
dif.q0.95.M1.jul <- quantile.q0.95.M1.2dec.jul - quantile.q0.95.M1.1dec.jul
dif.q0.95.M1.aug <- quantile.q0.95.M1.2dec.aug - quantile.q0.95.M1.1dec.aug
rm('quantile.q0.95.M1.2dec.jun', 'quantile.q0.95.M1.1dec.jun',
   'quantile.q0.95.M1.2dec.jul', 'quantile.q0.95.M1.1dec.jul',
   'quantile.q0.95.M1.2dec.aug', 'quantile.q0.95.M1.1dec.aug')
gc()

dif.q0.95.M2.jun <- quantile.q0.95.M2.2dec.jun - quantile.q0.95.M2.1dec.jun
dif.q0.95.M2.jul <- quantile.q0.95.M2.2dec.jul - quantile.q0.95.M2.1dec.jul
dif.q0.95.M2.aug <- quantile.q0.95.M2.2dec.aug - quantile.q0.95.M2.1dec.aug
rm('quantile.q0.95.M2.2dec.jun', 'quantile.q0.95.M2.1dec.jun',
   'quantile.q0.95.M2.2dec.jul', 'quantile.q0.95.M2.1dec.jul',
   'quantile.q0.95.M2.2dec.aug', 'quantile.q0.95.M2.1dec.aug')
gc()

dif.q0.95.M3.jun <- quantile.q0.95.M3.2dec.jun - quantile.q0.95.M3.1dec.jun
dif.q0.95.M3.jul <- quantile.q0.95.M3.2dec.jul - quantile.q0.95.M3.1dec.jul
dif.q0.95.M3.aug <- quantile.q0.95.M3.2dec.aug - quantile.q0.95.M3.1dec.aug
rm('quantile.q0.95.M3.2dec.jun', 'quantile.q0.95.M3.1dec.jun',
   'quantile.q0.95.M3.2dec.jul', 'quantile.q0.95.M3.1dec.jul',
   'quantile.q0.95.M3.2dec.aug', 'quantile.q0.95.M3.1dec.aug')
gc()


is.0.dif <- function(dif){
  index <- rep(1:1500, times = nrow(dif) / 1500)
  sum <- rowsum(dif, group = index)
  mean <- sum /(nrow(dif) / 1500)
  
  IC <- apply(mean, 2, FUN = quantile, probs = c(0.025, 0.975))
  is0 <- IC[1, ] <= 0 & IC[2, ] >= 0
  
  return(is0)
}

is.0.q0.50.M1.jun <- is.0.dif(dif.q0.50.M1.jun)
is.0.q0.50.M1.jul <- is.0.dif(dif.q0.50.M1.jul)
is.0.q0.50.M1.aug <- is.0.dif(dif.q0.50.M1.aug)

is.0.q0.50.M2.jun <- is.0.dif(dif.q0.50.M2.jun)
is.0.q0.50.M2.jul <- is.0.dif(dif.q0.50.M2.jul)
is.0.q0.50.M2.aug <- is.0.dif(dif.q0.50.M2.aug)

is.0.q0.50.M3.jun <- is.0.dif(dif.q0.50.M3.jun)
is.0.q0.50.M3.jul <- is.0.dif(dif.q0.50.M3.jul)
is.0.q0.50.M3.aug <- is.0.dif(dif.q0.50.M3.aug)

is.0.q0.95.M1.jun <- is.0.dif(dif.q0.95.M1.jun)
is.0.q0.95.M1.jul <- is.0.dif(dif.q0.95.M1.jul)
is.0.q0.95.M1.aug <- is.0.dif(dif.q0.95.M1.aug)

is.0.q0.95.M2.jun <- is.0.dif(dif.q0.95.M2.jun)
is.0.q0.95.M2.jul <- is.0.dif(dif.q0.95.M2.jul)
is.0.q0.95.M2.aug <- is.0.dif(dif.q0.95.M2.aug)

is.0.q0.95.M3.jun <- is.0.dif(dif.q0.95.M3.jun)
is.0.q0.95.M3.jul <- is.0.dif(dif.q0.95.M3.jul)
is.0.q0.95.M3.aug <- is.0.dif(dif.q0.95.M3.aug)


save(bay.kriging.q0.50.M1, quantile.q0.50.M1, 
     is.0.q0.50.M1.jun, is.0.q0.50.M1.jul, is.0.q0.50.M1.aug,
     bay.kriging.q0.50.M2, quantile.q0.50.M2, 
     is.0.q0.50.M2.jun, is.0.q0.50.M2.jul, is.0.q0.50.M2.aug,
     bay.kriging.q0.50.M3, quantile.q0.50.M3, 
     is.0.q0.50.M3.jun, is.0.q0.50.M3.jul, is.0.q0.50.M3.aug,
     bay.kriging.q0.95.M1, quantile.q0.95.M1, 
     is.0.q0.95.M1.jun, is.0.q0.95.M1.jul, is.0.q0.95.M1.aug,
     bay.kriging.q0.95.M2, quantile.q0.95.M2, 
     is.0.q0.95.M2.jun, is.0.q0.95.M2.jul, is.0.q0.95.M2.aug,
     bay.kriging.q0.95.M3, quantile.q0.95.M3, 
     is.0.q0.95.M3.jun, is.0.q0.95.M3.jul, is.0.q0.95.M3.aug,
     file = 'spatial.interpolation.RData')

