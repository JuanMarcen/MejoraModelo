rm(list=ls())


#df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc_lag.rds")

#Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y_lag.rds")

ind <- which(Y$station==230)

df_madrid <- df_conj_filled_sc[ind,] # original data frame

source('harmonics.R')
df_harm <- as.data.frame(cbind(df_madrid, cs(df_madrid$l, 1:2)))
colnames(df_harm)[3] <- 'Y'

source('eBIC.R')
mod_nulo_q0.5 <- rq(Y~1, tau= 0.5, data = df_harm)
mod_nulo_q0.95 <- rq(Y~1, tau= 0.95, data = df_harm)

#----Harmonics----
mod_harm <- step_rq_eBIC(
  mod_nulo_q0.5, data = df_harm,
  scope = as.formula(paste('Y ~',paste0('c.', 1:2, ' + s.', 1:2,collapse = '+'))),
  gamma = 0,
  harmonics = TRUE
) 

# After this, we know we need harmonics of order 1

#----Next step: Create the df with all the needed data as Jorge said----
library(dplyr)

# Anomalies of all grid cells of the data frame
ind <- which(df_harm$station == 230)
ind_jja <- which(df_harm$t[ind] >= 22 & df_harm$t[ind] <= 51 & df_harm$l[ind] >= 152)
  
for (j in 4:18){
  var <- names(df_harm)[j]
  formula <- as.formula(paste(var, "~ s.1 + c.1"))
  mod <- lm(formula, data = df_harm[ind,], subset = ind_jja)
  preds <- predict(mod, newdata = data.frame(
    c.1 = df_harm$c.1[ind],
    s.1 = df_harm$s.1[ind]
  ))
  
  res <- df_harm[ind, var] - preds
  print(sum(preds[ind_jja] - mod$fitted.values <= 1e-10))
  
  df_harm[ind,var] <- res 
}

not_28 <- which(df_madrid$l >= 149)

# lags of the anomalies
aux <- data.frame(
  t = df_madrid$t,
  l = df_madrid$l,
  g300 = df_harm$g300,
  g500 = df_harm$g500,
  g700 = df_harm$g700
)

aux <- aux %>%
  group_by(t) %>%
  mutate(across(c(g300,g500,g700),
                .fns = ~lag(.),
                .names = '{.col}_lag')) %>%
  as.data.frame() %>% na.omit()

g300_g300_lag <- aux$g300 - aux$g300_lag

#g500-g500_lag
g500_g500_lag <- aux$g500 - aux$g500_lag

#g700-g700_lag
g700_g700_lag <- aux$g700 - aux$g700_lag

#substract observation of day 28th may
g300 <- df_harm$g300[not_28]
g500 <- df_harm$g500[not_28]
g700 <- df_harm$g700[not_28]
g300_g500 <- df_harm$g300[not_28] - df_harm$g500[not_28]
g300_g700 <- df_harm$g300[not_28] - df_harm$g700[not_28]

# DATA FRAME FINAL: OBS JUNE - AUGUST
jun_ag <- which(df_madrid$l >= 152 & df_madrid$l <= 243) # months of june and august
# this subset is the one that may have all the possibe data available
ind_aux <- match(jun_ag, not_28) # selection of of the values of the residuals in june - august
df <- df_harm[jun_ag,] %>%
  select(Date, Y, l, t, c.1, s.1) %>%
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

# After having the data frame we can begin with the construction of a model

#----Step 1: We start with the harmonic model chosen----
formula <- as.formula(paste('Y ~ scale(c.1) + scale(s.1) +', paste(paste0('scale(',names(df)[7:ncol(df)],')'), collapse = '+')))

mod_step1 <- step_rq_eBIC(
  initial_model = mod_harm,
  data = df,
  scope = formula
)

#----Step 2: Add the lags of the chosen variables (not harmonics)----

# lags of order 1, 2 and 3 of the residuals
lag_res <- data.frame(
  t = df_madrid$t[not_28],
  l = df_madrid$l[not_28],
  g300, g500, g700,
  g300_g500, g300_g700,
  g300_g300_lag, g500_g500_lag, g700_g700_lag
) %>%
  group_by(t) %>%
  mutate(across(
    -c(1),
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
strip_scale <- function(x) sub("^scale\\((.*)\\)$", "\\1", x)
vars <- strip_scale(
  names(mod_step1$coefficients)[2:length(names(mod_step1$coefficients))]
)
vars_g_lag1 <- paste0(vars[grepl("^g", vars)], "_lag1")

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], vars_g_lag1), ')')), collapse = '+'))
)

mod_step2 <- step_rq_eBIC(
  initial_model = mod_step1,
  data = df_final,
  scope = formula
)

#----Step 3: Add lag order 2 of those selected with lag order 1----
vars <- strip_scale(
  names(mod_step2$coefficients)[2:length(names(mod_step2$coefficients))]
)
vars_lag1 <- vars[grepl("lag1$", vars)]
vars_lag2 <- sub("lag1$", "lag2", vars_lag1)

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], vars_lag2), ')')), collapse = '+'))
)

mod_step3 <- step_rq_eBIC(
  initial_model = mod_step2,
  data = df_final,
  scope = formula
)

#----Step 4: Add lag order 3 of those selected with lag order 2----
vars <- strip_scale(
  names(mod_step3$coefficients)[2:length(names(mod_step3$coefficients))]
)
vars_lag2 <- vars[grepl("lag2$", vars)]
vars_lag3 <- sub("lag2$", "lag3", vars_lag2)

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], vars_lag3), ')')), collapse = '+'))
)

mod_step4 <- step_rq_eBIC(
  initial_model = mod_step3,
  data = df_final,
  scope = formula
)


#----Step 5: interactions between all variables chosen----
vars <- strip_scale(
  names(mod_step4$coefficients)[2:length(names(mod_step4$coefficients))]
)
interactions <- combn(vars,2)
interactions <- apply(interactions, 2, function(x) paste(x[1], x[2], sep=":"))
interactions <- interactions[2: length(interactions)]

for (inter in interactions) {
  # Separar las dos variables que componen la interacción
  aux <- unlist(strsplit(inter, ":"))
  var1 <- aux[1]
  var2 <- aux[2]
  
  # Crear un nombre para la nueva variable
  new_var <- paste0(var1, "_x_", var2)
  
  # Calcular la interacción y escalarla
  df_final[[new_var]] <- df_final[[var1]] * df_final[[var2]]
}

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], gsub(':','_x_',interactions)), ')')), collapse = '+'))
)

mod_step5 <- step_rq_eBIC(
  initial_model = mod_step4,
  data = df_final,
  scope = formula
)

#----Step 6: Polynomials of order 2 ----
vars <- strip_scale(
  names(mod_step5$coefficients)[2:length(names(mod_step5$coefficients))]
)
polynomials <- c(paste0(
  'I(', vars[3:length(vars)], '^2)'
))

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], polynomials), ')')), collapse = '+'))
)


mod_step6 <- step_rq_eBIC(
  initial_model = mod_step5,
  data = df_final,
  scope = formula
)

#----Step 7: Polynomials of order 3 (those chosen for order 2) ----
vars <- strip_scale(
  names(mod_step6$coefficients)[2:length(names(mod_step6$coefficients))]
)
polynomials <- grep("^I\\(.*\\^2\\)$", vars, value = TRUE)
polynomials <- gsub("\\^2\\)", "^3)", polynomials) 

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], polynomials), ')')), collapse = '+'))
)

mod_step7 <- step_rq_eBIC(
  initial_model = mod_step6,
  data = df_final,
  scope = formula
)


vars_air_column <- strip_scale(
  names(mod_step7$coefficients)[2:length(names(mod_step7$coefficients))]
)


#----SAME PROCESS BUT FOR THE CORNERS---- 
corners_df <- data.frame(matrix(nrow = length(not_28), ncol = 0))

# lags of the anomalies
aux <- data.frame(
  t = df_madrid$t,
  l = df_madrid$l,
  g300_45_.10 = df_harm$g300_45_.10,
  g300_45_5 = df_harm$g300_45_5,
  g300_35_.10 = df_harm$g300_35_.10,
  g300_35_5 = df_harm$g300_35_5,
  g500_45_.10 = df_harm$g500_45_.10,
  g500_45_5 = df_harm$g500_45_5,
  g500_35_.10 = df_harm$g500_35_.10,
  g500_35_5 = df_harm$g500_35_5,
  g700_45_.10 = df_harm$g700_45_.10,
  g700_45_5 = df_harm$g700_45_5,
  g700_35_.10 = df_harm$g700_35_.10,
  g700_35_5 = df_harm$g700_35_5
)

aux <- aux %>%
  group_by(t) %>%
  mutate(across(-c(1),
                .fns = ~lag(.),
                .names = '{.col}_lag')) %>%
  as.data.frame() %>% na.omit()


corners <- c('_45_.10', '_45_5', '_35_.10', '_35_5')

# normal and lags
for (g in c('g300','g500', 'g700')){
  for (c in corners){
    corners_df[[paste0(g, c)]] <- df_harm[not_28, paste0(g, c)]
    corners_df[[paste0(g, c , '_', g, c, '_lag')]] <- aux[[paste0(g, c)]] - aux[[paste0(g, c, '_lag')]]
  }
}

# substractions
for (c in corners){
  for (g in c('g500', 'g700')){
    corners_df[[paste0('g300', c, '_', g, c)]] <- df_harm[not_28, paste0('g300', c)] - df_harm[not_28, paste0(g, c)]
  }
}

corners_df_final <- df_harm[jun_ag,] %>%
  select(Date, Y, l, t, c.1, s.1) %>%
  bind_cols(corners_df[ind_aux,]) %>%
  as.data.frame()
df <- cbind(corners_df_final, df[,7:ncol(df)])

# Step 8
formula <- as.formula(
  paste('Y ~', 
        paste(vars_air_column, collapse = '+'), ' + ', 
        paste(paste0('scale(',names(corners_df_final)[7:ncol(corners_df_final)],')'), collapse = '+')))
mod_step8 <- step_rq_eBIC(
  initial_model = mod_step7,
  data = df,
  scope = formula,
  replacements = list(
    c('g300_45_.10','g500_45_.10','g300_45_.10_g500_45_.10'),
    c('g300_45_.10','g700_45_.10','g300_45_.10_g700_45_.10'),
    c('g300_45_5','g500_45_5','g300_45_5_g500_45_5'),
    c('g300_45_5','g700_45_5','g300_45_5_g700_45_5'),
    c('g300_35_.10','g500_35_.10','g300_35_.10_g500_35_.10'),
    c('g300_35_.10','g700_35_.10','g300_35_.10_g700_35_.10'),
    c('g300_35_5','g500_35_5','g300_35_5_g500_35_5'),
    c('g300_35_5','g700_35_5','g300_35_5_g700_35_5')
  ))


# Step 2

# lags of order 1, 2 and 3 of the residuals
lag_res <- cbind(df_madrid$t[not_28], df_madrid$l[not_28], corners_df)
colnames(lag_res)[1:2] <- c('t', 'l')
lag_res <- lag_res %>%
  group_by(t) %>%
  mutate(across(
    -c(1),
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
df_final <- cbind(df_final, lags_only)

# chosen variables after step 1
strip_scale <- function(x) sub("^scale\\((.*)\\)$", "\\1", x)
vars <- strip_scale(
  names(mod_step8$coefficients)[2:length(names(mod_step8$coefficients))]
)
vars <- setdiff(vars, vars_air_column)
vars_g_lag1 <- paste0(vars[grepl("^g", vars)], "_lag1")

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], vars_g_lag1), ')')), collapse = '+'))
)

mod_step2 <- step_rq_eBIC(
  initial_model = mod_step1,
  data = corners_df_final,
  scope = formula
)

#----Step 3: Add lag order 2 of those selected with lag order 1----
vars <- strip_scale(
  names(mod_step2$coefficients)[2:length(names(mod_step2$coefficients))]
)
vars_lag1 <- vars[grepl("lag1$", vars)]
vars_lag2 <- sub("lag1$", "lag2", vars_lag1)

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], vars_lag2), ')')), collapse = '+'))
)

mod_step3 <- step_rq_eBIC(
  initial_model = mod_step2,
  data = corners_df_final,
  scope = formula
)

#----Step 4: Add lag order 3 of those selected with lag order 2----
vars <- strip_scale(
  names(mod_step3$coefficients)[2:length(names(mod_step3$coefficients))]
)
vars_lag2 <- vars[grepl("lag2$", vars)]
vars_lag3 <- sub("lag2$", "lag3", vars_lag2)

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], vars_lag3), ')')), collapse = '+'))
)

mod_step4 <- step_rq_eBIC(
  initial_model = mod_step3,
  data = corners_df_final,
  scope = formula
)


#----Step 5: interactions between all variables chosen----
vars <- strip_scale(
  names(mod_step4$coefficients)[2:length(names(mod_step4$coefficients))]
)
interactions <- combn(vars,2)
interactions <- apply(interactions, 2, function(x) paste(x[1], x[2], sep=":"))
interactions <- interactions[2: length(interactions)]

for (inter in interactions) {
  # Separar las dos variables que componen la interacción
  aux <- unlist(strsplit(inter, ":"))
  var1 <- aux[1]
  var2 <- aux[2]
  
  # Crear un nombre para la nueva variable
  new_var <- paste0(var1, "_x_", var2)
  
  # Calcular la interacción y escalarla
  corners_df_final[[new_var]] <- corners_df_final[[var1]] * corners_df_final[[var2]]
}

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], gsub(':','_x_',interactions)), ')')), collapse = '+'))
)

mod_step5 <- step_rq_eBIC(
  initial_model = mod_step4,
  data = corners_df_final,
  scope = formula
)

#----Step 6: Polynomials of order 2 ----
vars <- strip_scale(
  names(mod_step5$coefficients)[2:length(names(mod_step5$coefficients))]
)
polynomials <- c(paste0(
  'I(', vars[3:length(vars)], '^2)'
))

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], polynomials), ')')), collapse = '+'))
)


mod_step6 <- step_rq_eBIC(
  initial_model = mod_step5,
  data = corners_df_final,
  scope = formula
)

#----Step 7: Polynomials of order 3 (those chosen for order 2) ----
vars <- strip_scale(
  names(mod_step6$coefficients)[2:length(names(mod_step6$coefficients))]
)
polynomials <- grep("^I\\(.*\\^2\\)$", vars, value = TRUE)
polynomials <- gsub("\\^2\\)", "^3)", polynomials) 

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', c(vars[3:length(vars)], polynomials), ')')), collapse = '+'))
)

mod_step7 <- step_rq_eBIC(
  initial_model = mod_step6,
  data = corners_df_final,
  scope = formula
)


vars_corners<- strip_scale(
  names(mod_step7$coefficients)[2:length(names(mod_step7$coefficients))]
)




#----FINAL MODEL----
DF_FINAL <- cbind(df_final, corners_df_final)

vars_final <- unique(c(vars_air_column,vars_corners))

formula <- as.formula(
  paste('Y ~',
        paste(c(vars[1:2], paste0('scale(', vars_final[3:length(vars_final)], ')')), collapse = '+'))
)

# 1. model with 80 covariates
mod_final_1 <- rq(formula, data = DF_FINAL, tau = 0.50)
mod_nulo_q0.5 <- rq(Y ~ 1, data = DF_FINAL, tau = 0.50)

R1_final_1 <- 1 - mod_final_1$rho / mod_nulo_q0.5$rho

# 2. step between the 80 covariates
mod_final_2 <- step_rq_eBIC(
  initial_model = mod_harm,
  data = DF_FINAL,
  scope = formula,
  replacements = list(
    c('g300','g500','g300_g500'),
    c('g300','g700','g300_g700'),
    c('g300_45_.10','g500_45_.10','g300_45_.10_g500_45_.10'),
    c('g300_45_.10','g700_45_.10','g300_45_.10_g700_45_.10'),
    c('g300_45_5','g500_45_5','g300_45_5_g500_45_5'),
    c('g300_45_5','g700_45_5','g300_45_5_g700_45_5'),
    c('g300_35_.10','g500_35_.10','g300_35_.10_g500_35_.10'),
    c('g300_35_.10','g700_35_.10','g300_35_.10_g700_35_.10'),
    c('g300_35_5','g500_35_5','g300_35_5_g500_35_5'),
    c('g300_35_5','g700_35_5','g300_35_5_g700_35_5')
  ))

# final model with 53 covariates

R1_final_2 <- mod_final_2$R1

# saving of important things
save(
  DF_FINAL,
  df_harm,
  mod_final_1,
  mod_final_2,
  mod_nulo_q0.5,
  vars_air_column,
  vars_corners,
  vars_final,
  file = 'data_q0.5_madrid.RData'
)


##
formula <- as.formula(paste('Y ~', paste(vars_air_column, collapse = '+')))
mod <- rq(formula, data = DF_FINAL)
1 -mod$rho/mod_nulo_q0.5$rho

#----EXTRA----
library(lubridate)
source('functions.R')

df_dia <- rho_day(mod_final_1, mod_final_2, df_madrid[jun_ag,])
df_year <- rho_year(mod_final_1, mod_final_2, df_madrid[jun_ag,])

plot(1:92, df_dia$rho_l_q0.5, type='l', 
     main = 'Madrid (Retiro) (días) (armónicos)',
     ylab = expression(rho[l](tau)), xlab = 'l', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:92, df_dia$rho_l_q0.95, col = 'blue')

plot(1:64, df_year$rho_t_q0.5, type='l', 
     main = 'Madrid (Retiro) (años) (armónicos)',
     ylab = expression(rho[t](tau)), xlab = 't', ylim = c(0, 1))
abline(h = 0.5, col = 'red')
lines(1:64, df_year$rho_t_q0.95, col = 'blue')


mod <- lm(scale(Y) ~ 1, data = df_madrid)
plot(unique(df_madrid$t), tapply(mod$residuals, df_madrid$t, mean))

mod <- lm(scale(Y) ~ s.1+c.1+poly(g300,3)*poly(g500,3)*poly(g700,3) +
            poly(g300_lag,3)*poly(g500_lag,3)*poly(g700_lag,3) 
          + df_madrid$g300_45_.10[not_28], data = aux)
plot(unique(aux$t), tapply(mod$residuals, aux$t, mean))
summary(mod)
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