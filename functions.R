# Functions

# Rho by days and years
rho_day <- function(model_q1, model_q2){
  
  day_month <- unique(format(Y$Date, "%d-%m"))
  df_dia <- matrix(NA, nrow=length(day_month), ncol=2)
  df_dia <- as.data.frame(df_dia, row.names = day_month)
  colnames(df_dia) <- c('rho_l_q0.5', 'rho_l_q0.95')
  
  for (i in 1:length(day_month)){
    
    day <- which(format(Y$Date[ind], "%d-%m") == day_month[i])
    df_dia[i,'rho_l_q0.5']<-sum(model_q1$residuals[day] < 0, na.rm=T) / 64
    df_dia[i,'rho_l_q0.95']<-sum(model_q2$residuals[day] < 0, na.rm=T) / 64
    
  }
 
  return(df_dia) 
}

rho_year <- function(model_q1, model_q2){
  
  year <- unique(year(Y$Date[ind]))
  df_year <- matrix(NA, nrow=length(year), ncol=2)
  df_year <- as.data.frame(df_year, row.names = year)
  colnames(df_year) <- c('rho_t_q0.5','rho_t_q0.95')
  
  for (i in 1:length(year)){
    
    j <- which(year(Y$Date[ind]) == year[i])
    df_year[i,'rho_t_q0.5'] <- sum(model_q1$residuals[j] < 0, na.rm=T) / 92
    df_year[i,'rho_t_q0.95'] <- sum(model_q2$residuals[j] < 0, na.rm=T) / 92
    
  }
  
  return(df_year)
}