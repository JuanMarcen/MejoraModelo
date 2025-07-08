# Extendeded BIC (eBIC) for quantile regression models

# obtained with 'rq'

# generalize to bayesian models
library(quantreg)

my_eBIC<-function(model, gamma, p){
  
  loss <- model$rho
  
  n <- length(model$fitted.values)
  k <- length(coef(model)) - 1

  
  eBIC <- 2 * n * log((1 / n) * loss) + k * log(n) + 2 * gamma * log(choose(p,k))
  #AIC <- - 2 * n * log (0.5*0.5) +  2*n*log((1/n)*loss) + 2*n + k*2
  
  return(eBIC)
  
}


step_rq_eBIC<-function(initial_model, data, scope, 
                       gamma = 0, 
                       trace = TRUE, harmonics = FALSE){
  
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
  steps <- list()
  steps[[1]] <- list(formula = formula_current, eBIC = best_eBIC)
  
  improved <- TRUE
  
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
        cat("Added:", best_var, "| eBIC =", round(best_eBIC, 2), "\n")
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


