rm(list=ls())

#df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc_lag.rds")
#Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y_lag.rds")

ind <- which(Y$station==230)

df_madrid <- df_conj_filled_sc[ind,]

source('harmonics.R')
df_harm <- as.data.frame(cbind(df_madrid, cs(df_madrid$l, 1:10)))
colnames(df_harm)[1] <- 'Y'

source('eBIC.R')
mod_nulo_q0.5 <- rq(Y~1, tau= 0.5, data = df_harm)
mod_nulo_q0.95 <- rq(Y~1, tau= 0.95, data = df_harm)

mod_harm <- step_rq_eBIC(
  mod_nulo_q0.5, data = df_harm,
  scope = as.formula(paste('Y ~',paste0('c.', 1:5, ' + s.', 1:5,collapse = '+'))),
  gamma = 1,
  harmonics = TRUE
) 

# my_eBIC(mod_nulo_q0.95,gamma = 1, 10)
# for (i in 1:5){
#   mod <- rq(as.formula(paste('Y ~',paste0('c.', 1:i, ' + s.', 1:i,collapse = '+'))),
#             data = df_harm, tau = 0.95)
#   cat('armonico',i,'eBIC:',my_eBIC(mod,gamma = 1,p = 10),
#       'R1:',1-mod$rho/mod_nulo_q0.95$rho , '\n' )
# }

# After this, we know we need harmonics of order 1
# Next step: Create the df with all the needed data as Jorge said

jun_ag <- which(df_madrid$l >= 5 & df_madrid$l <= 96) # motnhs of june and august
# this subset is the one that may have all the possibe data available

library(dplyr)
df <- df_harm %>%
  select(Y,l,t,c.1,s.1) %>%
  as.data.frame()

# g300 residuals (substract obs of 28th may)
not_28 <- which(df$l >=2)
mod <- lm(df_madrid$g300[not_28] ~ df$s.1[not_28] + df$s.1[not_28])
g300_res <- scale(mod$residuals)

#g500 res
mod <- lm(df_madrid$g500[not_28] ~ df$s.1[not_28] + df$s.1[not_28])
g500_res <- scale(mod$residuals)

#g700 res
mod <- lm(df_madrid$g700[not_28] ~ df$s.1[not_28] + df$s.1[not_28])
g700_res <- scale(mod$residuals)













## EXTRA 

mod_q0.5 <- rq(as.formula(paste('Y ~',paste0('c.', 1:2, ' + s.', 1:2,collapse = '+'))),
          data = df_harm, tau = 0.5)
mod_q0.95 <- rq(as.formula(paste('Y ~',paste0('c.', 1:2, ' + s.', 1:2,collapse = '+'))),
          data = df_harm, tau = 0.95)
library(lubridate)
source('functions.R')

df_dia <- rho_day(mod_q0.5, mod_q0.95, df_harm)
df_dia_2 <- rho_day(mod_q0.5, mod_q0.95, df_harm)
df_year <- rho_year(mod_q0.5, mod_q0.95, df_harm)
plot(1:92, df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (armónicos)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:92, df_dia_2$rho_l_q0.5,col = 'blue')
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