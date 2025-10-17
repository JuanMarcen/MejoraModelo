# CL-BIC
library(quantreg)
library(sf)
library(sp)

exp_weights <- function(dist_matrix, h, scale = TRUE){
  
  dist_matrix <- dist_matrix / h
  
  dist_kernel <- exp(-dist_matrix)
  
  w <- numeric()
  
  for (i in 1:dim(dist_matrix)[1]){
    w[i] <- 1 / sum(dist_kernel[i, ])
  }
  
  if (scale == TRUE){
    n <- dim(dist_matrix)[1]
    w <- n * w / sum(w) 
  }
  
  return(w)
}

log_lik_rq <- function(model){
  
  loss <- model$rho
  tau <- model$tau
  n <- length(model$fitted.values)
  
  loglik <- n * log(tau * (1 - tau)) - n * log( (1 / n) * loss) - n
  
  return(loglik)
}

logCL_rq <- function(weights = 1, logliks){
  
  logCL <- sum(weights * logliks)
  
  return(logCL)
}

eff_number_param <- function(models, weights = 1, stations_df){
  
  # Computation of sensitivity matrix H
  # tau and n are the same for every station
  tau <- models[[1]]$tau
  n <- length(models[[1]]$fitted.values)
  num_coef <- length(models[[1]]$coefficients)
  H <- matrix(0, ncol = num_coef, nrow = num_coef)
  # sigmas for each model
  sigmas <- numeric()
  for (i in 1:dim(stations_df)[1]){
    #sigma = 1/n * loss(beta)
    model <- models[[i]]
    n <- length(model$fitted.values)
    sigmas[i] <- model$rho / n
  }
  # x_tl(s_i) %*% t(x_tl(s_i)) the same for each station also
  X_H <- list()
  for (i in 1:dim(stations_df)[1]){
    model <- models[[i]]
    
    x <- matrix(0, ncol = num_coef, nrow = num_coef)
    for (tl in 1:n){
      x <- x +  (tau * (1 - tau) / sigmas[i]^2) * model$x[tl, ] %*% t(model$x[tl, ])
    }
    
    X_H[[as.character(stations_df$STAID[i])]] <- x
  }
  
  # I need the weights (argument)
  if (length(weights) == 1 && weights == 1) {
    weights <- rep(1, times = dim(stations_df)[1])
  }
  # cat('weights: ', weights, '\n')
  for (i in 1:dim(stations_df)[1]){
    H <- H + weights[i] * X_H[[i]]
  }
  
  # Computation of variability matrix J
  psi_tau <- function(u, tau){
    return(tau - as.numeric(u < 0))
  }
  
  # i have the weights and sigmas
  X_J <- list()
  for (i in 1:dim(stations_df)[1]){
    model <- models[[i]]
    check_sq <- psi_tau(model$y - model$x %*% model$coefficients, tau)^2
    
    x <- matrix(0, ncol = num_coef, nrow = num_coef)
    for (tl in 1:n){
      x <- x + (1 / sigmas[i]^2) * check_sq[tl] * model$x[tl, ] %*% t(model$x[tl, ])
    }
    
    X_J[[as.character(stations_df$STAID[i])]] <- x
  }
  
  J <- matrix(0, ncol = num_coef, nrow = num_coef)
  for (i in 1:dim(stations_df)[1]){
    J <- J + weights[i]^2 * X_J[[i]]
  }
  
  # Computation of Godobame informatin matrix G
  G <- H %*% solve(J) %*% H
  
  eff_param <- sum(diag(H %*% solve(G))) - 1 #substract the intercept
  
  return(eff_param)
}

CLBIC <- function(models, weights = 1, eff_param = FALSE, gamma = 0, p = 100,
                  stations_df){
  # models as a list
  
  loglik <- numeric(length(models))
  logCL <- numeric()
  
  for (i in 1:length(models)){
    loglik[i] <- log_lik_rq(models[[i]])
    #cat('modelo',i,'\t',loglik[i],'\n')
  }
  #cat("weights:", weights, "\n")
  #cat("loglik:", loglik, "\n")
  logCL <- logCL_rq(weights, loglik)
  
  #cat('logCL', logCL ,'\n')
  # supposedly the models have the same amount of parameters and observations
  # n <- length(models[[1]]$fitted.values)
  # n is equal to the number of iid groups (40)
  n <- dim(stations_df)[1]
  if (eff_param == TRUE){
    if (length(weights) == 1 && weights == 1) {
      weights <- rep(1, times = dim(stations_df)[1])
    }
    k <- eff_number_param(models, weights, stations_df = stations_df)
    #cat('Número de parámetros efectivo', k, '\n')
  }else{
    k <- length(coef(models[[1]])) - 1
    #cat('Número de parámetros', k, '\n')
  }
  
  #cat(n,k,'\n')
  
  CLBIC <- suppressWarnings(
    -2 * logCL + k * log(n) + 2 * gamma * log(choose(p,k))
  )
  
  
  return(CLBIC)
}

R1_global <- function(formula, stations_df, data, tau){
  # R1 
  data.aux <- data[is.element(data$station, stations_df$STAID), ]
  mod_null <- rq(Y ~ as.factor(station), data = data.aux, tau = tau)
  rho_station <- rep(NA, dim(stations_df)[1])
  
  for (i in 1:length(rho_station)){
    ind <- which(data$station == stations_df$STAID[i])
    mod.aux <- rq(formula, data = data[ind, ], tau = tau)
    rho_station[i] <- mod.aux$rho
  }
  
  return(1 - sum(rho_station)/mod_null$rho)
}

step_rq_CLBIC<-function(initial_models,
                        null_models,
                        data, 
                        stations_df,
                        scope, 
                        weights = 1,
                        eff_param = FALSE,
                        gamma = 0, 
                        p = 100,
                        trace = TRUE,
                        #harmonics = FALSE,
                        replacements = list(
                          c('g300','g500','g300_g500'),
                          c('g300','g700','g300_g700')
                          )){
  
  # in order to substract the scale part
  strip_scale <- function(x) sub("^scale\\((.*)\\)$", "\\1", x)
  
  if (eff_param == TRUE){
    cat('Se va a utilizar el número efectivo de parámetros\n')
  }else{
    cat('Se va a utlizar el número de parámetros\n')
  }
  
  # size of covariates set
  vars <- labels(terms(scope)) #formula terms
  #print(vars)
  p <- length(vars) # at the moment not necessary because gamma = 0
  response <- as.character(scope[[2]]) #response variable
  #print(response)
  tau <- initial_models[[1]]$tau
  #print(tau)
  #data <- get(as.character(initial_models$call$data))
  #print(dim(data))
  
  # if (harmonics == TRUE){
  #   n <- length(vars) / 2
  #   vars <- sapply(1:n, function(i) paste0('c.',i,' + s.',i))
  # }
  
  # Initial model
  formula_current <- initial_models[[1]]$formula
  cat('Initial model: \n')
  print(formula_current)
  cat('\n')
  model_current <- initial_models[[1]]
  
  best_CLBIC <- CLBIC(initial_models, weights, eff_param, gamma, p, stations_df = stations_df)
  
  selected_vars <- attr(terms(formula(model_current)), "term.labels")
  #print(selected_vars)
  remaining_vars <- setdiff(vars, selected_vars)
  #print(remaining_vars)
  
  improved <- TRUE
  combos_probados <- character(0)
  
  models_solution <- list(models = initial_models,
                          formula = formula_current,
                          CLBIC = best_CLBIC)
  
  while (improved && length(remaining_vars) > 0){
    improved <- FALSE
    CLBICs <- c()
    models <- list()
    error_occurred <- FALSE
    
    for (var in remaining_vars){
      error_occurred <- FALSE
      models_stations <- list()
      formula_try <- as.formula(
        paste(response, '~', paste(c(selected_vars, var), collapse = '+'))
      )
      
      for (i in 1:dim(stations_df)[1]){
        ind <- which(data$station == stations_df$STAID[i])
        model_try <- tryCatch(
          rq(formula_try, data = data, subset = ind, tau = tau),
          error = function(e) return(NULL)
        )
        
        if (is.null(model_try)) {
          if (trace) cat("Saltando variable (error en ajuste):", var, '\n')
          error_occurred <- TRUE
          break
        }
        
        model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
        models_stations[[as.character(stations_df$STAID[i])]] <- model_try
      }
      
      if (error_occurred) {
        CLBICs <- c(CLBICs, Inf)
        models[[var]] <- list(models = NULL, formula = formula_try, CLBIC = Inf)
        next  # ⚠️ SALTA A LA SIGUIENTE VARIABLE
      }
      
      # Si todo fue bien:
      CLBIC_val <- CLBIC(models_stations, weights, eff_param, gamma, p, stations_df = stations_df)
      CLBICs <- c(CLBICs, CLBIC_val)
      models[[var]] <- list(models = models_stations, 
                            formula = formula_try, 
                            CLBIC = CLBIC_val)
    }
    
    #print(CLBICs)
    min_CLBIC <- min(CLBICs)
    
    #print(min_CLBIC)
    if (is.na(min_CLBIC) || is.na(min_CLBIC)){
      min_CLBIC <- Inf
    }
    if(min_CLBIC < best_CLBIC){
      best_var <- remaining_vars[which.min(CLBICs)]
      selected_vars <- c(selected_vars, best_var)
      remaining_vars <- setdiff(remaining_vars, best_var)
      
      model_current <- models[[best_var]]$models[[1]][1]
      formula_current <- models[[best_var]]$formula
      best_CLBIC <- min_CLBIC
      
      improved <- TRUE
      
      if (trace){
        cat("Added:", strip_scale(best_var), "| CLBIC =", round(best_CLBIC, 2), "\n")
      }
      
      models_solution <- models[[best_var]]
    }
    
    if (improved && length(replacements) > 0){
      for (combo in replacements){
        v1 <- as.character(combo[1])
        v2 <- as.character(combo[2])
        v1_v2 <- as.character(combo[3])
        
        # check if v1 and v2 in selected vars
        if (all(c(v1, v2) %in% strip_scale(selected_vars)) && 
            v1_v2 %in% strip_scale(vars)){
          
          #cat('1', '\n')
          
          combo_id <- v1_v2
          if (combo_id %in% combos_probados){
            next
          }
          
          # new_selected <- setdiff(selected_vars, paste0('scale(', c(v1,v2), ')'))
          # new_selected <- c(new_selected, paste0('scale(', v1_v2, ')'))
          
          new_selected <- setdiff(selected_vars, c(v1,v2))
          new_selected <- c(new_selected, v1_v2)
          
          models_stations <- list()
          
          formula_try <- as.formula(
            paste(response, '~', paste(new_selected, collapse = '+'))
          )
          
          error_occurred <- FALSE
          
          for (i in 1:dim(stations_df)[1]){
            ind <- which(data$station == stations_df$STAID[i])
            model_try <- tryCatch(
              rq(formula_try, data = data, subset = ind, tau = tau),
              error = function(e) return(NULL)
            )
            
            if (is.null(model_try)) {
              if (trace) cat("Saltando variable (error en ajuste):", var, '\n')
              error_occurred <- TRUE
              break
            }
            
            model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
            models_stations[[as.character(stations_df$STAID[i])]] <- model_try
          }
          
          if (error_occurred != TRUE){
            CLBIC_val <- CLBIC(models_stations, weights, eff_param, gamma, p, stations_df = stations_df)
            
            combos_probados <- c(combos_probados, combo_id)
            
            if(CLBIC_val < best_CLBIC){
              if (trace){
                cat("Reemplazando",
                    paste0(strip_scale(v1), "y", strip_scale(v2)),
                    "→", strip_scale(v1_v2),
                    "| CLBIC mejora a", round(CLBIC_val, 2), "\n")
              }
              
              # actualizar estado
              selected_vars  <- new_selected
              remaining_vars <- setdiff(vars, selected_vars)
              models[[best_var]] <- list(models = models_stations, 
                                         formula = formula_try, 
                                         CLBIC = CLBIC_val)
              model_current <- models[[best_var]]$models[[1]][1]
              formula_current <- models[[best_var]]$formula
              best_CLBIC <- CLBIC_val
              
              # Forzar que el bucle vuelva a empezar con
              # el nuevo conjunto (por si hay más mejoras)
              improved <- TRUE
              models_solution <- models[[best_var]]
              break 
              
            }else{
              cat('Reemplazamiento por ', v1_v2, ' no mejora', '\n')
            }
            
          }else{
            cat('Reemplazamiento por ', v1_v2, ' hay error', '\n')
          }
        }
        
        if (all(c(v1_v2, v1) %in% strip_scale(selected_vars)) && 
            v2 %in% strip_scale(vars)){
          
          #cat('2', '\n')
          
          combo_id <- v1_v2
          if (combo_id %in% combos_probados){
            next
          }
          
          # new_selected <- setdiff(selected_vars, paste0('scale(', c(v1,v2), ')'))
          # new_selected <- c(new_selected, paste0('scale(', v1_v2, ')'))
          
          new_selected <- setdiff(selected_vars, c(v1_v2))
          new_selected <- c(new_selected, v2)
          #cat(new_selected, '\n')
          
          models_stations <- list()
          
          formula_try <- as.formula(
            paste(response, '~', paste(new_selected, collapse = '+'))
          )
          
          error_occurred <- FALSE
          for (i in 1:dim(stations_df)[1]){
            ind <- which(data$station == stations_df$STAID[i])
            model_try <- tryCatch(
              rq(formula_try, data = data, subset = ind, tau = tau),
              error = function(e) return(NULL)
            )
            
            if (is.null(model_try)) {
              if (trace) cat("Saltando variable (error en ajuste):", var, '\n')
              error_occurred <- TRUE
              break
            }
            
            model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
            models_stations[[as.character(stations_df$STAID[i])]] <- model_try
          }
          
          if (error_occurred != TRUE){
            CLBIC_val <- CLBIC(models_stations, weights, eff_param, gamma, p, stations_df = stations_df)
            
            combos_probados <- c(combos_probados, combo_id)
            
            if (trace){
              cat("Reemplazando",
                  paste0(strip_scale(v1_v2)),
                  "→", strip_scale(v2),
                  "| CLBIC cambia a", round(CLBIC_val, 2), "\n")
            }
            
            # actualizar estado
            selected_vars  <- new_selected
            remaining_vars <- setdiff(vars, selected_vars)
            models[[best_var]] <- list(models = models_stations, 
                                       formula = formula_try, 
                                       CLBIC = CLBIC_val)
            model_current <- models[[best_var]]$models[[1]][1]
            formula_current <- models[[best_var]]$formula
            best_CLBIC <- CLBIC_val
            
            # Forzar que el bucle vuelva a empezar con
            # el nuevo conjunto (por si hay más mejoras)
            improved <- TRUE
            models_solution <- models[[best_var]]
            break 
            
          }else{
            cat('Reemplazamiento por ', v2, ' hay error', '\n')
        
        }
        }
        
        if (all(c(v1_v2, v2) %in% strip_scale(selected_vars)) && 
            v1 %in% strip_scale(vars)){
          
          #cat('3', '\n')
          
          combo_id <- v1_v2
          if (combo_id %in% combos_probados){
            next
          }
          
          # new_selected <- setdiff(selected_vars, paste0('scale(', c(v1,v2), ')'))
          # new_selected <- c(new_selected, paste0('scale(', v1_v2, ')'))
          
          new_selected <- setdiff(selected_vars, c(v1_v2))
          new_selected <- c(new_selected, v1)
          #cat(new_selected, '\n')
          
          models_stations <- list()
          
          formula_try <- as.formula(
            paste(response, '~', paste(new_selected, collapse = '+'))
          )
          
          error_occurred <- FALSE
          for (i in 1:dim(stations_df)[1]){
            ind <- which(data$station == stations_df$STAID[i])
            model_try <- tryCatch(
              rq(formula_try, data = data, subset = ind, tau = tau),
              error = function(e) return(NULL)
            )
            
            if (is.null(model_try)) {
              if (trace) cat("Saltando variable (error en ajuste):", var, '\n')
              error_occurred <- TRUE
              break
            }
            
            model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
            models_stations[[as.character(stations_df$STAID[i])]] <- model_try
          }
          
          if (error_occurred != TRUE){
            CLBIC_val <- CLBIC(models_stations, weights, eff_param, gamma, p, stations_df = stations_df)
            
            combos_probados <- c(combos_probados, combo_id)
            
            if (trace){
              cat("Reemplazando",
                  paste0(strip_scale(v1_v2)),
                  "→", strip_scale(v1),
                  "| CLBIC cambia a", round(CLBIC_val, 2), "\n")
            }
            
            # actualizar estado
            selected_vars  <- new_selected
            remaining_vars <- setdiff(vars, selected_vars)
            models[[best_var]] <- list(models = models_stations, 
                                  formula = formula_try, 
                                  CLBIC = CLBIC_val)
            model_current <- models[[best_var]]$models[[1]][1]
            formula_current <- models[[best_var]]$formula
            best_CLBIC <- CLBIC_val
            
            # Forzar que el bucle vuelva a empezar con
            # el nuevo conjunto (por si hay más mejoras)
            improved <- TRUE
            
            models_solution <- models[[best_var]]
            break 
            
          }else{
            cat('Reemplazamiento por ', v1, ' hay error', '\n')
            
          }
        }
      
    }
    
    }
    
  }
  
  
  # eBIC
  model_current$CLBIC <- best_CLBIC
  
  
  models_solution$R1 <- suppressWarnings(
    R1_global(formula_current, stations_df, data, tau)
  )
  # R1
  # model_null <- suppressWarnings(
  #   rq(paste(response, '~ 1'), data = data, tau = tau)
  # )
  # model_current$R1 <- 1 - model_current$rho / model_null$rho
  
  # if (trace && harmonics == TRUE) {
  #   cat("\nFinal model:\n")
  #   print(formula_current)
  #   cat('R1_final: ', model_current$R1, '\n')
  #   cat('eBIC final: ', model_current$eBIC, '\n')
  # }
  
  if (trace) {
    cat("\nFinal model:\n")
    print(formula_current)
    if (eff_param ==TRUE) cat('Número de parámetros efectivo: ', 
                              eff_number_param(models_solution$models, weights, stations_df = stations_df),
                              '\n')
    cat('CLBIC initial: ', CLBIC(initial_models, weights, eff_param, gamma, p, stations_df = stations_df),
        '| CLBIC final: ', model_current$CLBIC, '\n')
    
    cat('Global R1: ', models_solution$R1)
  }
  
  return(models_solution)
}




# #pruebas para elección de EFFECTIVE NUMBER OF PARAMETERS
# modelos_prueba <- list()
# formula <- as.formula('Y ~ s.1 + c.1 + g300 + g500 + g700 +
#   g300_g300_lag + g500_g500_lag + g700_g700_lag')
# for (i in 1:dim(stations)[1]){
#   ind <- which(df_jun_ag$station == stations$STAID[i])
#   mod <- rq(formula, tau = 0.95, data = df_jun_ag, subset = ind)
#   # mod$R1 <- 1 - mod$rho / models_null[[as.character(stations$STAID[i])]]$rho
#   modelos_prueba[[as.character(stations$STAID[i])]] <- mod
# }

