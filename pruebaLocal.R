# Improving existing models
# example in one station

rm(list=ls())

# Data uploading
df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
vars_q0.5 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/vars_q0.5.rds")
vars_q0.95 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/vars_q0.95.rds")
stations <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/stations.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

# Subset Madrid (Retiro)
ind <- which(Y$station==230)

# Pruebas
vars_q0.5 <- ifelse(grepl("^I\\(.*\\)$", vars_q0.5),
                    paste0("`", vars_q0.5, "`"),
                    vars_q0.5)
vars_q0.95 <- ifelse(grepl("^I\\(.*\\)$", vars_q0.95),
                     paste0("`", vars_q0.95, "`"),
                     vars_q0.95)

formula_q0.5 <- as.formula(
  paste('Y~', paste(
    vars_q0.5,collapse = '+'
  ),'+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))')
)
formula_q0.95 <- as.formula(
  paste('Y~', paste(
    vars_q0.95,collapse = '+'
  ),'+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))')
)


# Modelos
library(quantreg)
<<<<<<< HEAD
mod_q0.5 <- rq(formula_q0.5,data=df_conj_filled_sc,subset=ind,tau=0.50)
mod_q0.95 <- rq(formula_q0.95,data=df_conj_filled_sc,subset=ind,tau=0.95)

# R1
mod_nulo_q0.5 <- rq(Y~1,data=df_conj_filled_sc,subset=ind,tau=0.5)  
mod_nulo_q0.95 <- rq(Y~1,data=df_conj_filled_sc,subset=ind,tau=0.95) 

R1_q0.5 <- 1 - mod_q0.5$rho / mod_nulo_q0.5$rho
R1_q0.95 <- 1 - mod_q0.95$rho / mod_nulo_q0.95$rho

# rho

rho_q0.5 <- sum(mod_q0.5$residuals<0) / length(ind)
rho_q0.95 <- sum(mod_q0.95$residuals<0) / length(ind)

# Día 
library(lubridate)

day_month <- unique(format(Y$Date, "%d-%m"))
df_dia <- matrix(NA, nrow=length(day_month), ncol=2)
df_dia <- as.data.frame(df_dia, row.names = day_month)
colnames(df_dia) <- c('rho_l_q0.5', 'rho_l_q0.95')

for (i in 1:length(day_month)){
  
  day <- which(format(Y$Date[ind], "%d-%m") == day_month[i])
  df_dia[i,'rho_l_q0.5']<-sum(mod_q0.5$residuals[day] < 0, na.rm=T) / 64
  df_dia[i,'rho_l_q0.95']<-sum(mod_q0.95$residuals[day] < 0, na.rm=T) / 64
  
}

# Años 
year <- unique(year(Y$Date[ind]))
df_year <- matrix(NA, nrow=length(year), ncol=2)
df_year <- as.data.frame(df_year, row.names = year)
colnames(df_year) <- c('rho_t_q0.5','rho_t_q0.95')

for (i in 1:length(year)){
  
  j <- which(year(Y$Date[ind]) == year[i])
  df_year[i,'rho_t_q0.5'] <- sum(mod_q0.5$residuals[j] < 0, na.rm=T) / 92
  df_year[i,'rho_t_q0.95'] <- sum(mod_q0.95$residuals[j] < 0, na.rm=T) / 92
  
}


# Gráficos
setwd('C:/Users/jumar/OneDrive/Escritorio/Github/MejoraModelo')
png("solo_armonicos.png", width = 1400, height = 600, res = 150)
par(mfrow = c(1,2))
plot(1:92, df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (armónicos)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:92, df_dia$rho_l_q0.95)
abline(h = 0.95, col = 'red')

plot(1:64, df_year$rho_t_q0.5, type='l', 
     main = 'Madrid (Retiro) (años) (armónicos)',
     ylab = expression(rho[t](tau)), xlab = 't', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:64, df_year$rho_t_q0.95)
abline(h = 0.95, col = 'red')
dev.off()

