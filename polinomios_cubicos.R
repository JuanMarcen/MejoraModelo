rm(list=ls())

# Data uploading
df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
vars_q0.5 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/vars_q0.5.rds")
vars_q0.95 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/vars_q0.95.rds")
stations <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/stations.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

# Subset Madrid (Retiro)
ind <- which(Y$station==230)

# Adding 3rd degree polynomials
library(quantreg)

formula <- as.formula(
  paste('Y ~', paste(colnames(df_conj_filled_sc)[2:16], collapse = '+'), '+',
  paste0('`', colnames(df_conj_filled_sc)[17:46], '`', collapse = '+'))
  )

# Armónicos
formula <- as.formula(
  paste('Y ~', paste(colnames(df_conj_filled_sc)[2:16], collapse = '+'), '+',
        paste0('`', colnames(df_conj_filled_sc)[17:46], '`', collapse = '+'),
        '+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))'
        , '+ t'
        )
)

# With interactions
# vars <- paste0('`', colnames(df_conj_filled_sc)[2:46], '`')
# 
# terminos_principales <- paste(vars, collapse = ' + ')
# 
# interacciones_vars <- combn(vars, 2, FUN = function(x) paste(x, collapse = ' : '))
# terminos_interaccion_vars <- paste(interacciones_vars, collapse = ' + ')
# 
# seno <- 'I(sin(2*pi*l/365))'
# coseno <- 'I(cos(2*pi*l/365))'
# 
# interacciones_estacionales <- c(
#   paste(vars, seno, sep = ' : '),
#   paste(vars, coseno, sep = ' : ')
# )
# terminos_interaccion_estacionales <- paste(interacciones_estacionales, collapse = ' + ')
# 
# formula_completa <- as.formula(
#   paste('Y ~',
#         terminos_principales, '+',
#         terminos_interaccion_vars, '+',
#         seno, '+', coseno, '+',
#         terminos_interaccion_estacionales
#   )
# )



# Models
mod_nulo_q0.5 <- rq(Y ~ 1, data = df_conj_filled_sc, subset = ind, tau = 0.5)
mod_q0.5<-step(mod_nulo_q0.5, scope = formula, direction = 'forward',
               pen = log(length(ind)), trace = F)

mod_nulo_q0.95 <- rq(Y ~ 1, data = df_conj_filled_sc, subset = ind, tau = 0.95)
mod_q0.95<-step(mod_nulo_q0.95, scope = formula, direction = 'forward',
               pen = log(length(ind)), trace = F)


# R1
R1_q0.5 <- 1 - mod_q0.5$rho / mod_nulo_q0.5$rho
R1_q0.95 <- 1 - mod_q0.95$rho / mod_nulo_q0.95$rho


# Rho
rho_q0.5 <- sum(mod_q0.5$residuals < 0) / length(ind)
rho_q0.95 <- sum(mod_q0.95$residuals < 0) / length(ind)

library(lubridate)
source('functions.R')

df_dia <- rho_day(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind, ])
df_year <- rho_year(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind, ])
df_geop <- rho_geop(mod_q0.5, mod_q0.95, df_conj_filled_sc[ind, ], 'g300', 50)

# Gráficos
setwd('C:/Users/jumar/OneDrive/Escritorio/Github/MejoraModelo')
png("pol_gr_3_armonicos_t.png", width = 1400, height = 600, res = 150)
par(mfrow = c(1,2))
plot(1:92, df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (pol gr3)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:92, df_dia$rho_l_q0.95)
abline(h = 0.95, col = 'red')

plot(1:64, df_year$rho_t_q0.5, type='l', 
     main = 'Madrid (Retiro) (años) (pol gr3)',
     ylab = expression(rho[t](tau)), xlab = 't', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:64, df_year$rho_t_q0.95)
abline(h = 0.95, col = 'red')
dev.off()

plot(1:50, df_geop$rho_g_q0.5, type='l', 
     main = 'Madrid (Retiro) (geop) (pol gr3)',
     ylab = expression(rho[g](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:50, df_geop$rho_g_q0.95)
abline(h = 0.95, col = 'red')

