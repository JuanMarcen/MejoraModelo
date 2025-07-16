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
    vars_q0.5,collapse = '+')
    , '+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))'
    #, '+', paste0('I(sin(2*pi*l/365)):',vars_q0.5, collapse = '+')
    #, '+', paste0('I(cos(2*pi*l/365)):',vars_q0.5, collapse = '+')
    #, '+', paste0('t:',vars_q0.5, collapse = '+')
    , '+ t'
  )
)
formula_q0.95 <- as.formula(
  paste('Y~', paste(
    vars_q0.95,collapse = '+')
    , '+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))'
    #, '+', paste0('I(sin(2*pi*l/365)):',vars_q0.95, collapse = '+')
    #, '+', paste0('I(cos(2*pi*l/365)):',vars_q0.95, collapse = '+')
    #, '+', paste0('t:',vars_q0.95, collapse = '+')
    , '+ t'
  )
)


# Modelos
library(quantreg)
# sin interacciones con armonicos o t
mod_q0.5 <- rq(formula_q0.5,data=df_conj_filled_sc,subset=ind,tau=0.50)
mod_q0.95 <- rq(formula_q0.95,data=df_conj_filled_sc,subset=ind,tau=0.95)

# modelos nulos
mod_nulo_q0.5 <- rq(Y~1,data=df_conj_filled_sc,subset=ind,tau=0.5)  
mod_nulo_q0.95 <- rq(Y~1,data=df_conj_filled_sc,subset=ind,tau=0.95) 

# seleccion con interacciones
mod_q0.5<-step(mod_nulo_q0.5, scope = formula_q0.5, direction = 'forward',
               pen = log(length(ind)), trace = F)
mod_q0.95<-step(mod_nulo_q0.95, scope = formula_q0.95, direction = 'forward',
               pen = log(length(ind)), trace = F)

# R1
R1_q0.5 <- 1 - mod_q0.5$rho / mod_nulo_q0.5$rho
R1_q0.95 <- 1 - mod_q0.95$rho / mod_nulo_q0.95$rho

# rho
rho_q0.5 <- sum(mod_q0.5$residuals<0) / length(ind)
rho_q0.95 <- sum(mod_q0.95$residuals<0) / length(ind)

# Día 
library(lubridate)
source('functions.R')

df_dia <- rho_day(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind,])
df_year <- rho_year(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind,])
df_geop <- rho_geop(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind, ], 'g300', 50)

df_dia_p <- rho_day(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind,], extra = T)
df_year_m <- rho_year(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind,], extra = T)



# Gráficos
setwd('C:/Users/jumar/OneDrive/Escritorio/Github/MejoraModelo')
png("armonicos_t.png", width = 1400, height = 600, res = 150)
par(mfrow = c(1,2))
plot(1:92, df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (armónicos)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:92, df_dia_p$rho_l_q0.5_p1, col = 'blue')
lines(1:92, df_dia_p$rho_l_q0.5_p2, col = 'forestgreen')
lines(1:92, df_dia_p$rho_l_q0.5_p3, col = 'purple')
lines(1:92, df_dia$rho_l_q0.95)
lines(1:92, df_dia_p$rho_l_q0.95_p1, col = 'blue')
lines(1:92, df_dia_p$rho_l_q0.95_p2, col = 'forestgreen')
lines(1:92, df_dia_p$rho_l_q0.95_p3, col = 'purple')
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('Todo', '1960-1980', '1981-2001', '2002-2023'),
       col = c('black', 'blue', 'forestgreen', 'purple'), lty = 1, horiz = T,
       cex = 0.5)

plot(1:64, df_year$rho_t_q0.5, type='l', 
     main = 'Madrid (Retiro) (años) (armónicos)',
     ylab = expression(rho[t](tau)), xlab = 't', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:64, df_year_m$rho_t_q0.5_jun, col = 'blue')
lines(1:64, df_year_m$rho_t_q0.5_jul, col = 'forestgreen')
lines(1:64, df_year_m$rho_t_q0.5_ag, col = 'purple')
lines(1:64, df_year$rho_t_q0.95)
lines(1:64, df_year_m$rho_t_q0.95_jun, col = 'blue')
lines(1:64, df_year_m$rho_t_q0.95_jul, col = 'forestgreen')
lines(1:64, df_year_m$rho_t_q0.95_ag, col = 'purple')
abline(h = 0.95, col = 'red')
legend('bottom', legend = c('Todo', 'Junio', 'Julio', 'Agosto'),
       col = c('black', 'blue', 'forestgreen', 'purple'), lty = 1, horiz = T,
       cex = 0.5)

dev.off()

plot(1:50, df_geop$rho_g_q0.5, type='l', 
     main = 'Madrid (Retiro) (geop) (pol gr3)',
     ylab = expression(rho[g](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:50, df_geop$rho_g_q0.95)
abline(h = 0.95, col = 'red')

