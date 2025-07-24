# CL-BIC
library(quantreg)
library(sf)
library(sp)

exp_weights <- function(dist_matrix, h, scale = TRUE){
  
  dist_matrix <- dist_matrix / h
  
  dist_kernel <- exp(-dist_matrix)
  
  w <- numeric()
  
  for (i in 1:dim(stations)[1]){
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

CLBIC <- function(models, weights = 1, gamma = 0, p = 100){
  
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
  n <- length(models[[1]]$fitted.values)
  k <- length(coef(models[[1]])) - 1
  #cat(n,k,'\n')
  
  
  CLBIC <- -2 * logCL + k * log(n) + 2 * gamma * log(choose(p,k))
  
  
  return(CLBIC)
}

step_rq_CLBIC<-function(initial_models,
                        null_models,
                        data, 
                        stations,
                        scope, 
                        weights = 1,
                        gamma = 0, 
                        trace = TRUE,
                        #harmonics = FALSE,
                        replacements = list(
                          c('g300','g500','g300_g500'),
                          c('g300','g700','g300_g700')
                          )){
  
  # in order to substract the scale part
  strip_scale <- function(x) sub("^scale\\((.*)\\)$", "\\1", x)
  
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
  best_CLBIC <- CLBIC(initial_models, weights, gamma, p)
  
  selected_vars <- attr(terms(formula(model_current)), "term.labels")
  #print(selected_vars)
  remaining_vars <- setdiff(vars, selected_vars)
  #print(remaining_vars)
  
  improved <- TRUE
  combos_probados <- character(0)
  
  while (improved && length(remaining_vars) > 0){
    improved <- FALSE
    CLBICs <- c()
    models <- list()
    error_occurred <- FALSE
    
    for (var in remaining_vars){
      models_stations <- list()
      formula_try <- as.formula(
        paste(response, '~', paste(c(selected_vars, var), collapse = '+'))
      )
      
      for (i in 1:dim(stations)[1]){
        ind <- which(data$station == stations$STAID[i])
        model_try <- tryCatch(
          rq(formula_try, data = data, subset = ind, tau = tau),
          error = function(e) return(NULL)
        )
        
        if (is.null(model_try)) {
          if (trace) cat("Saltando variable (error en ajuste):", var, '\n',
                         'Cometido en la estación', stations$STAID[i], "\n")
          error_occurred <- TRUE
          break
        }
        
        model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
        models_stations[[as.character(stations$STAID[i])]] <- model_try
      }
      
      if (error_occurred) {
        CLBICs <- c(CLBICs, Inf)
        models[[var]] <- list(models = NULL, formula = formula_try, CLBIC = Inf)
        next  # ⚠️ SALTA A LA SIGUIENTE VARIABLE
      }
      
      # Si todo fue bien:
      CLBIC_val <- CLBIC(models_stations, weights, gamma, p)
      CLBICs <- c(CLBICs, CLBIC_val)
      models[[var]] <- list(models = models_stations, 
                            formula = formula_try, 
                            CLBIC = CLBIC_val)
    }
    
    min_CLBIC <- min(CLBICs)
    
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
          
          cat('1', '\n')
          
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
          
          for (i in 1:dim(stations)[1]){
            ind <- which(data$station == stations$STAID[i])
            model_try <- tryCatch(
              rq(formula_try, data = data, subset = ind, tau = tau),
              error = function(e) return(NULL)
            )
            
            if (is.null(model_try)) {
              if (trace) cat("Saltando variable (error en ajuste):", var, '\n',
                             'Cometido en la estación', stations$STAID[i], "\n")
              error_occurred <- TRUE
              break
            }
            
            model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
            models_stations[[as.character(stations$STAID[i])]] <- model_try
          }
          
          if (error_occurred != TRUE){
            CLBIC_val <- CLBIC(models_stations, weights, gamma, p)
            
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
          
          cat('2', '\n')
          
          combo_id <- v1_v2
          if (combo_id %in% combos_probados){
            next
          }
          
          # new_selected <- setdiff(selected_vars, paste0('scale(', c(v1,v2), ')'))
          # new_selected <- c(new_selected, paste0('scale(', v1_v2, ')'))
          
          new_selected <- setdiff(selected_vars, c(v1_v2))
          new_selected <- c(new_selected, v2)
          cat(new_selected, '\n')
          
          models_stations <- list()
          
          formula_try <- as.formula(
            paste(response, '~', paste(new_selected, collapse = '+'))
          )
          
          for (i in 1:dim(stations)[1]){
            ind <- which(data$station == stations$STAID[i])
            model_try <- tryCatch(
              rq(formula_try, data = data, subset = ind, tau = tau),
              error = function(e) return(NULL)
            )
            
            if (is.null(model_try)) {
              if (trace) cat("Saltando variable (error en ajuste):", var, '\n',
                             'Cometido en la estación', stations$STAID[i], "\n")
              error_occurred <- TRUE
              break
            }
            
            model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
            models_stations[[as.character(stations$STAID[i])]] <- model_try
          }
          
          if (error_occurred != TRUE){
            CLBIC_val <- CLBIC(models_stations, weights, gamma, p)
            
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
          
          cat('3', '\n')
          
          combo_id <- v1_v2
          if (combo_id %in% combos_probados){
            next
          }
          
          # new_selected <- setdiff(selected_vars, paste0('scale(', c(v1,v2), ')'))
          # new_selected <- c(new_selected, paste0('scale(', v1_v2, ')'))
          
          new_selected <- setdiff(selected_vars, c(v1_v2))
          new_selected <- c(new_selected, v1)
          cat(new_selected, '\n')
          
          models_stations <- list()
          
          formula_try <- as.formula(
            paste(response, '~', paste(new_selected, collapse = '+'))
          )
          
          for (i in 1:dim(stations)[1]){
            ind <- which(data$station == stations$STAID[i])
            model_try <- tryCatch(
              rq(formula_try, data = data, subset = ind, tau = tau),
              error = function(e) return(NULL)
            )
            
            if (is.null(model_try)) {
              if (trace) cat("Saltando variable (error en ajuste):", var, '\n',
                             'Cometido en la estación', stations$STAID[i], "\n")
              error_occurred <- TRUE
              break
            }
            
            model_try$R1 <- 1 - model_try$rho / null_models[[i]]$rho
            models_stations[[as.character(stations$STAID[i])]] <- model_try
          }
          
          if (error_occurred != TRUE){
            CLBIC_val <- CLBIC(models_stations, weights, gamma, p)
            
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
    cat('CLBIC initial: ', CLBIC(initial_models), '| CLBIC final: ', model_current$CLBIC, '\n')
  }
  
  return(models_solution)
}

