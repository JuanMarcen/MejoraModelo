rm(list = ls())
# mantengo 
# vars, final.chain, coords, stations_dist, stations, df,
library(stringr)
grid <- readRDS('grid.rds')
grid_km <- readRDS('grid_km.rds')
grid_elev <- readRDS('grid_elev.rds')
grid_dist <- readRDS('grid_dist.rds')

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

basura <- bay.kriging(final.chain.q0.50.M2, vars.M2, coords_km, grid_km, stations_dist)
