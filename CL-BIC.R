# CL-BIC

log_lik_rq <- function(model){
  
  loss <- model$rho
  tau <- model$tau
  n <- length(model$fitted.values)
  
  loglik <- n * log(tau * (1 - tau)) - n * log( (1 / n) * loss) - n
  
  return(loglik)
}

logCL_rq <- function(weights, logliks){
  logCL <- sum(weights * logliks)
  return(logCL)
}

CLBIC <- function(models, gamma = 0, p = 100){
  
  loglik <- numeric(length(models))
  logCL <- numeric()
  
  for (i in 1:length(models)){
    loglik[i] <- log_lik_rq(models[[i]])
    #cat('modelo',i,'\t',loglik[i],'\n')
  }
  
  logCL <- logCL_rq(1, loglik)
  cat('logCL', logCL ,'\n')
  # supposedly the models have the same amount of parameters and observations
  n <- length(models[[1]]$fitted.values)
  k <- length(coef(models[[1]])) - 1
  cat(n,k,'\n')
  
  
  CLBIC <- -2 * logCL + k * log(n) + 2 * gamma * log(choose(p,k))
  
  
  return(CLBIC)
}

step_rq_CLBIC<-function(initial_models, data, scope, 
                       gamma = 0, 
                       trace = TRUE, harmonics = FALSE,
                       replacements = list(
                         c('g300','g500','g300_g500'),
                         c('g300','g700','g300_g700')
                       )){
  
  # in order to subustract the scale part
  strip_scale <- function(x) sub("^scale\\((.*)\\)$", "\\1", x)
  
  # size of covariates set
  vars <- labels(terms(scope)) #formula terms
  #print(vars)
  p <- length(vars) 
  response <- as.character(scope[[2]]) #response variable
  #print(response)
  tau <- initial_model$tau
  #print(tau)
  #data <- get(as.character(initial_model$call$data))
  #print(dim(data))
  
  if (harmonics == TRUE){
    n <- length(vars) / 2
    vars <- sapply(1:n, function(i) paste0('c.',i,' + s.',i))
  }
  
  # Initial model
  formula_current <- initial_model$formula
  model_current <- initial_model
  best_eBIC <- my_eBIC(model_current, gamma, p)
  
  selected_vars <- attr(terms(formula(model_current)), "term.labels")
  #print(selected_vars)
  remaining_vars <- setdiff(vars, selected_vars)
  #print(remaining_vars)
  
  improved <- TRUE
  combos_probados <- character(0)
  
  while (improved && length(remaining_vars) > 0){
    improved <- FALSE
    eBICs <- c()
    models <- list()
    
    for (var in remaining_vars){
      formula_try <- as.formula(
        paste(response, '~', paste(c(selected_vars, var), collapse = '+'))
      )
      # Manejo del error en rq()
      model_try <- tryCatch(
        rq(formula_try, data = data, tau = tau),
        error = function(e) return(NULL)
      )
      
      # Si el modelo falló, pasa al siguiente
      if (is.null(model_try)) {
        if (trace) cat("Saltando variable (error en ajuste):", var, "\n")
        eBIC_val <- Inf
        eBICs <- c(eBICs, eBIC_val)
        models[[var]] <- list(model = NULL, formula = formula_try, eBIC = eBIC_val)
        next
      }
      
      eBIC_val <- my_eBIC(model_try, gamma, p)
      eBICs <- c(eBICs,eBIC_val)
      models[[var]] <- list(model = model_try, formula = formula_try, eBIC = eBIC_val)
    }
    
    min_eBIC <- min(eBICs)
    
    if(min_eBIC < best_eBIC){
      best_var <- remaining_vars[which.min(eBICs)]
      selected_vars <- c(selected_vars, best_var)
      remaining_vars <- setdiff(remaining_vars, best_var)
      
      model_current <- models[[best_var]]$model
      formula_current <- models[[best_var]]$formula
      best_eBIC <- min_eBIC
      
      #steps[[length(steps) + 1]] <- list(formula = formula_current, eBIC = best_eBIC)
      
      improved <- TRUE
      
      if (trace){
        cat("Added:", strip_scale(best_var), "| eBIC =", round(best_eBIC, 2), "\n")
      }
    }
    
    if (improved && length(replacements) > 0){
      for (combo in replacements){
        v1 <- as.character(combo[1])
        v2 <- as.character(combo[2])
        v1_v2 <- as.character(combo[3])
        
        # check if v1 and v2 in selected vars
        if (all(c(v1, v2) %in% strip_scale(selected_vars)) && 
            v1_v2 %in% strip_scale(vars)){
          
          combo_id <- v1_v2
          if (combo_id %in% combos_probados){
            next
          }
          
          new_selected <- setdiff(selected_vars, paste0('scale(', c(v1,v2), ')'))
          new_selected <- c(new_selected, paste0('scale(', v1_v2, ')'))
          
          formula_try <- as.formula(
            paste(response, '~', paste(new_selected, collapse = '+'))
          )
          
          model_try <- tryCatch(
            rq(formula_try, data = data, tau = tau),
            error = function(e) NULL
          )
          
          if (!is.null(model_try)){
            eBIC_val <- my_eBIC(model_try, gamma, p)
            
            combos_probados <- c(combos_probados, combo_id)
            
            if(eBIC_val < best_eBIC){
              if (trace){
                cat("Reemplazando",
                    paste0(strip_scale(v1), "y", strip_scale(v2)),
                    "→", strip_scale(v1_v2),
                    "| eBIC mejora a", round(eBIC_val, 2), "\n")
              }
              
              # actualizar estado
              selected_vars  <- new_selected
              remaining_vars <- setdiff(vars, selected_vars)
              model_current  <- model_try
              formula_current<- formula_try
              best_eBIC      <- eBIC_val
              
              # Forzar que el bucle vuelva a empezar con
              # el nuevo conjunto (por si hay más mejoras)
              improved <- TRUE
              break 
              
            }else{
              cat('Reemplazamiento por ', v1_v2, ' no mejora ', '\n')
            }
            
          }
        }
        
      }
      
    }
    
  }
  
  
  # eBIC
  model_current$eBIC <- best_eBIC
  
  # R1
  model_null <- suppressWarnings(
    rq(paste(response, '~ 1'), data = data, tau = tau)
  )
  model_current$R1 <- 1 - model_current$rho / model_null$rho
  
  if (trace && harmonics == TRUE) {
    cat("\nFinal model:\n")
    print(formula_current)
    cat('R1_final: ', model_current$R1, '\n')
    cat('eBIC final: ', model_current$eBIC, '\n')
  }
  
  if (trace && harmonics == FALSE) {
    cat("\nFinal model:\n")
    print(formula_current)
    cat('R1 initial: ', initial_model$R1, '| R1_final: ', model_current$R1, '\n')
    cat('eBIC initial: ', initial_model$eBIC, '| eBIC final: ', model_current$eBIC, '\n')
  }
  
  return(model_current)
}

