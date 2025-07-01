rm(list=ls())

df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

ind <- which(Y$station==230)

library(lubridate)
df_madrid <- df_conj_filled_sc[ind,]


# Dataframe with lags
library(dplyr)

df_madrid <- df_madrid %>%
  group_by(t) %>%
  mutate(across(all_of(names(df_madrid)[2:46]), ~lag(.), .names = "{.col}_lag")) %>%
  ungroup() %>%
  as.data.frame() %>%
  na.omit()

# Model

formula <- as.formula(paste(
  'Y ~',
  paste(paste0('`',colnames(df_madrid)[c(2:46,51:95)],'`'),collapse='+'),
  '+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))'
  , '+ t'
)
)

library(quantreg)
mod_nulo_q0.5 <- rq(Y ~ 1, data = df_madrid, tau = 0.5)
mod_q0.5<-step(mod_nulo_q0.5, scope = formula, direction = 'forward',
               pen = log(dim(df_madrid)[1]))

mod_nulo_q0.95 <- rq(Y ~ 1, data = df_madrid, tau = 0.95)
mod_q0.95<-step(mod_nulo_q0.95, scope = formula, direction = 'forward',
                pen = log(length(ind)))


# R1
R1_q0.5 <- 1 - mod_q0.5$rho / mod_nulo_q0.5$rho
R1_q0.95 <- 1 - mod_q0.95$rho / mod_nulo_q0.95$rho


# Rho
rho_q0.5 <- sum(mod_q0.5$residuals < 0) / (dim(df_madrid)[1])
rho_q0.95 <- sum(mod_q0.95$residuals < 0) / (dim(df_madrid)[1])

source('functions.R')

df_dia <- rho_day(mod_q0.5, mod_q0.95, df_madrid)
df_year <- rho_year(mod_q0.5, mod_q0.95, df_madrid)

# Gráficos
setwd('C:/Users/jumar/OneDrive/Escritorio/Github/MejoraModelo')
png("lags.png", width = 1400, height = 600, res = 150)
par(mfrow = c(1,2))
plot(unique(df_madrid$l), df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (pol gr3)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(unique(df_madrid$l), df_dia$rho_l_q0.95)
abline(h = 0.95, col = 'red')

plot(unique(df_madrid$t), df_year$rho_t_q0.5, type='l', 
     main = 'Madrid (Retiro) (años) (pol gr3)',
     ylab = expression(rho[t](tau)), xlab = 't', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(unique(df_madrid$t), df_year$rho_t_q0.95)
abline(h = 0.95, col = 'red')
dev.off()
