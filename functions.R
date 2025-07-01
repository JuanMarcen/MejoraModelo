# Functions

# Rho by days and years
rho_day <- function(model_q1, model_q2, df){
  
  day_month <- unique(df$l)
  df_dia <- matrix(NA, nrow=length(day_month), ncol=2)
  df_dia <- as.data.frame(df_dia, row.names = day_month)
  colnames(df_dia) <- c('rho_l_q0.5', 'rho_l_q0.95')
  
  year <- unique(df$t)
  
  for (i in 1:length(day_month)){
    
    day <- which(df$l == day_month[i])
    df_dia[i,'rho_l_q0.5']<-sum(model_q1$residuals[day] < 0, na.rm=T) / length(year)
    df_dia[i,'rho_l_q0.95']<-sum(model_q2$residuals[day] < 0, na.rm=T) / length(year)
    
  }
 
  return(df_dia) 
}

rho_year <- function(model_q1, model_q2, df){
  
  year <- unique(df$t)
  df_year <- matrix(NA, nrow=length(year), ncol=2)
  df_year <- as.data.frame(df_year, row.names = year)
  colnames(df_year) <- c('rho_t_q0.5','rho_t_q0.95')
  
  day_month <- unique(df$l)
  
  for (i in 1:length(year)){
    
    j <- which(df$t == year[i])
    df_year[i,'rho_t_q0.5'] <- sum(model_q1$residuals[j] < 0, na.rm=T) / length(day_month)
    df_year[i,'rho_t_q0.95'] <- sum(model_q2$residuals[j] < 0, na.rm=T) / length(day_month)
    
  }
  
  return(df_year)
}