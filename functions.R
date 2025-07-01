# Functions

# Rho by days and years
rho_day <- function(model_q1, model_q2, df, extra = F){
  
  if (extra == F){
    
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
    
  }else{
    
    day_month <- unique(df$l)
    df_dia <- matrix(NA, nrow=length(day_month), ncol=2*3)
    df_dia <- as.data.frame(df_dia, row.names = day_month)
    colnames(df_dia) <- c('rho_l_q0.5_p1', 'rho_l_q0.95_p1',
                          'rho_l_q0.5_p2', 'rho_l_q0.95_p2',
                          'rho_l_q0.5_p3', 'rho_l_q0.95_p3')
    
    p1 <- 1:21
    p2 <- 22:42
    p3 <- 43:64
    
    for (i in 1:length(day_month)){
      
      day <- which(df$l == day_month[i] & df$t %in% p1)
      df_dia[i,'rho_l_q0.5_p1']<-sum(model_q1$residuals[day] < 0, na.rm=T) / length(p1)
      df_dia[i,'rho_l_q0.95_p1']<-sum(model_q2$residuals[day] < 0, na.rm=T) / length(p1)
      
      day <- which(df$l == day_month[i] & df$t %in% p2)
      df_dia[i,'rho_l_q0.5_p2']<-sum(model_q1$residuals[day] < 0, na.rm=T) / length(p2)
      df_dia[i,'rho_l_q0.95_p2']<-sum(model_q2$residuals[day] < 0, na.rm=T) / length(p2)
      
      day <- which(df$l == day_month[i] & df$t %in% p3)
      df_dia[i,'rho_l_q0.5_p3']<-sum(model_q1$residuals[day] < 0, na.rm=T) / length(p3)
      df_dia[i,'rho_l_q0.95_p3']<-sum(model_q2$residuals[day] < 0, na.rm=T) / length(p3)
      
    }
    
  }
  
  return(df_dia) 
}

rho_year <- function(model_q1, model_q2, df, extra = F){
  
  if (extra == F){
    
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
    
  }else{
    
    year <- unique(df$t)
    df_year <- matrix(NA, nrow=length(year), ncol=2*3)
    df_year <- as.data.frame(df_year, row.names = year)
    colnames(df_year) <- c('rho_t_q0.5_jun','rho_t_q0.95_jun',
                           'rho_t_q0.5_jul','rho_t_q0.95_jul',
                           'rho_t_q0.5_ag','rho_t_q0.95_ag')
    
    jun <- 1:30
    jul <- 31:61
    ag <- 62:92
    
    for (i in 1:length(year)){
      
      j <- which(df$t == year[i] & df$l %in% jun)
      df_year[i,'rho_t_q0.5_jun'] <- sum(model_q1$residuals[j] < 0, na.rm=T) / length(jun)
      df_year[i,'rho_t_q0.95_jun'] <- sum(model_q2$residuals[j] < 0, na.rm=T) / length(jun)
      
      j <- which(df$t == year[i] & df$l %in% jul)
      df_year[i,'rho_t_q0.5_jul'] <- sum(model_q1$residuals[j] < 0, na.rm=T) / length(jul)
      df_year[i,'rho_t_q0.95_jul'] <- sum(model_q2$residuals[j] < 0, na.rm=T) / length(jul)
      
      j <- which(df$t == year[i] & df$l %in% ag)
      df_year[i,'rho_t_q0.5_ag'] <- sum(model_q1$residuals[j] < 0, na.rm=T) / length(ag)
      df_year[i,'rho_t_q0.95_ag'] <- sum(model_q2$residuals[j] < 0, na.rm=T) / length(ag)
      
  }
  
  }
  
  return(df_year)
}