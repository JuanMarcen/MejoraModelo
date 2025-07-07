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

#----Harmonics----
mod_harm <- step_rq_eBIC(
  mod_nulo_q0.5, data = df_harm,
  scope = as.formula(paste('Y ~',paste0('c.', 1:5, ' + s.', 1:5,collapse = '+'))),
  gamma = 0.5,
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

#----Next step: Create the df with all the needed data as Jorge said----
library(dplyr)

# Period of reference 1991-2020 (t: 32-61) nad without 28th may
ref_not_28 <- which(df_madrid$t >= 32 & df_madrid$t <= 61 & df_madrid$l >=2)

# g300 residuals (substract obs of 28th may)
mod <- lm(df_madrid$g300[ref_not_28] ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g300 <- scale(mod$residuals)

#g500 res
mod <- lm(df_madrid$g500[ref_not_28] ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g500 <- scale(mod$residuals)

#g700 res
mod <- lm(df_madrid$g700[ref_not_28] ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g700 <- scale(mod$residuals)

#g300-g500 res
mod <- lm(df_madrid$g300[ref_not_28]-df_madrid$g500[ref_not_28] ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g300_g500 <- scale(mod$residuals)

#g300-g700 res
mod <- lm(df_madrid$g300[ref_not_28]-df_madrid$g700[ref_not_28] ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g300_g700 <- scale(mod$residuals)


#g300-g300_lag
aux <- df_madrid[,c('g300','g500','g700','t','l')] %>% group_by(t) %>% 
  mutate(across(c(g300,g500,g700),
                .fns = ~lag(.),
                .names = '{.col}_lag')) %>%
  as.data.frame() %>% na.omit()

mod <- lm(aux$g300-aux$g300_lag ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g300_g300_lag <- scale(mod$residuals)

#g500-g500_lag
mod <- lm(aux$g500-aux$g500_lag ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g500_g500_lag <- scale(mod$residuals)

#g700-g700_lag
mod <- lm(aux$g700-aux$g700_lag ~ df_harm$s.1[ref_not_28] + df_harm$c.1[ref_not_28])
g700_g700_lag <- scale(mod$residuals)


# DATA FRAME FINAL: OBS JUNE - AUGUST
jun_ag <- which(df_madrid$l >= 5 & df_madrid$l <= 96) # months of june and august
# this subset is the one that may have all the possibe data available
ind_aux <- match(jun_ag, ref_not_28) # selection of of the values of the residuals in june - august
df <- df_harm[jun_ag,] %>%
  select(Y,l,t,c.1,s.1) %>%
  mutate(
    g300 = g300[ind_aux],
    g500 = g500[ind_aux],
    g700 = g700[ind_aux],
    g300_g500 = g300_g500[ind_aux],
    g300_g700 = g300_g700[ind_aux],
    g300_g300_lag = g300_g300_lag[ind_aux],
    g500_g500_lag = g500_g500_lag[ind_aux],
    g700_g700_lag = g700_g700_lag[ind_aux]
  ) %>%
  as.data.frame()


# recalculation of harmonics
df$l <- df$l - min(df$l) + 1 # 1 == 1 jun, 92 == 31 aug
df$c.1 <- cs(df$l, 1)[, 1]
df$s.1 <- cs(df$l, 1)[, 2]

# After having the data frame we can begin with the construction of a model

#----Step 1: We start with the harmonic model chosen----
formula <- as.formula(paste('Y ~', paste(names(df)[4:ncol(df)], collapse = '+')))

mod_step1 <- step_rq_eBIC(
  initial_model = mod_harm,
  data = df,
  scope = formula
)

#----Step 2: Add the lags of the chosen variables (not harmonics)----

# lags of order 1, 2 and 3 of the residuals
lag_res <- data.frame(
  t = df_madrid$t[ref_not_28],
  g300, g500, g700,
  g300_g500, g300_g700,
  g300_g300_lag, g500_g500_lag, g700_g700_lag
) %>%
  group_by(t) %>%
  mutate(across(
    everything(),
    .fns = list(
      lag1 = ~lag(.x, 1),
      lag2 = ~lag(.x, 2),
      lag3 = ~lag(.x, 3)
    ),
    .names = "{.col}_{fn}"
  )) %>% 
  na.omit() %>%
  as.data.frame()

# final data frame for all the following steps
lags_only <- lag_res %>% select(matches("_lag[123]$"))
df_final <- cbind(df, lags_only)

# chosen variables after step 1
vars <- names(mod_step1$coefficients)[2:length(names(mod_step1$coefficients))]
vars_g_lag1 <- paste0(vars[grepl("^g", vars)], "_lag1")

formula <- as.formula(
  paste('Y ~',
        paste(c(vars, vars_g_lag1), collapse = '+'))
)

mod_step2 <- step_rq_eBIC(
  initial_model = mod_step1,
  data = df_final,
  scope = formula
)

#----Step 3: Add lag order 2 of those selected with lag order 1----
vars <- names(mod_step2$coefficients)[2:length(names(mod_step2$coefficients))]
vars_lag1 <- vars[grepl("lag1$", vars)]
vars_lag2 <- sub("lag1$", "lag2", vars_lag1)

formula <- as.formula(
  paste('Y ~',
        paste(c(vars, vars_lag2), collapse = '+'))
)

mod_step3 <- step_rq_eBIC(
  initial_model = mod_step2,
  data = df_final,
  scope = formula
)

#----Step 4: Add lag order 3 of those selected with lag order 2----
vars <- names(mod_step3$coefficients)[2:length(names(mod_step3$coefficients))]
vars_lag3 <- sub("lag2$", "lag3", vars_lag2)

formula <- as.formula(
  paste('Y ~',
        paste(c(vars, vars_lag3), collapse = '+'))
)

mod_step4 <- step_rq_eBIC(
  initial_model = mod_step3,
  data = df_final,
  scope = formula
)


#----Step 5: interactions between all variables chosen----
vars <- names(mod_step4$coefficients)[2:length(names(mod_step4$coefficients))]

interactions <- combn(vars,2)
interactions <- apply(interactions, 2, function(x) paste(x[1], x[2], sep=":"))
interactions <- interactions[2: length(interactions)]

formula <- as.formula(
  paste('Y ~',
        paste(c(vars, interactions), collapse = '+'))
)

mod_step5 <- step_rq_eBIC(
  initial_model = mod_step4,
  data = df_final,
  scope = formula
)

#----Step 6: Polynomials of order 2 and 3 of all the chosen variables----
vars <- names(mod_step5$coefficients)[2:length(names(mod_step5$coefficients))]

polynomials <- c(paste0(
  'I(', vars, '^2)'
), paste(
  'I(', vars, '^3)'
))

formula <- as.formula(
  paste('Y ~',
        paste(c(vars, polynomials), collapse = '+'))
)

mod_step6 <- step_rq_eBIC(
  initial_model = mod_step5,
  data = df_final,
  scope = formula
)

##----EXTRA---- 

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