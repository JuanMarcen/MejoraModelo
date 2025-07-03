# Extendeded BIC (eBIC) for quantile regression models

# obtained with 'rq'

# generalize to bayesian models

my_eBIC<-function(model, gamma, p){
  
  loss <- model$rho
  
  n <- length(model$fitted.values)
  k <- length(coef(model)) 
  
  eBIC <- 2 * n * log((1 / n) * loss) + k * log(n) + 2 * gamma * log(choose(p,k))
  #AIC <- - 2 * n * log (0.5*0.5) +  2*n*log((1/n)*loss) + 2*n + k*2
  
  return(eBIC)
  
}


step_rq_eBIC<-function(data, response, tau = 0.5, gamma = 0.5, trace = TRUE){
  
  # size of covariates set
  vars <- paste0('`',setdiff(names(data),response),'`')
  p <- length(vars) + 1
  
  # Null model
  formula_current <- as.formula(paste(response, '~ 1'))
  model_current <- rq(formula = formula_current, data = data, tau = tau)
  best_eBIC <- my_eBIC(model_current, gamma, p)
  
  selected_vars <- c()
  remaining_vars <- vars
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
      model_try <- rq(formula_try, data = data, tau = tau)
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
  
  if (trace) {
    cat("\nFinal model:\n")
    print(formula_current)
  }
  
  return(list(model = model_current, eBIC = best_eBIC))
}


