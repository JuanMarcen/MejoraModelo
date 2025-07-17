rm(list=ls())

df <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc_lag.rds")
stations <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/stations.rds")

library(quantreg)
library(lubridate)
library(dplyr)

#----Dataframe----
df <- df[,c(1:18, (ncol(df) - 3):ncol(df))]

# Harmonics
source('harmonics.R')
df_harm <- cbind(df, cs(df$l,1))

# Anomalies. Ref 1981-2010
for (i in 1:dim(stations)[1]){
  ind <- which(df_harm$station == stations$STAID[i])
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
    #print(sum(preds[ind_jja] - mod$fitted.values <= 1e-10))
    
    df_harm[ind,var] <- res 
  }
  
}

# diferencias de anomalias
not_28 <- not_28 <- which(df_harm$l >= 149)
df_final <- df_harm[not_28, ]

corners <- c('','_45_.10', '_45_5', '_35_.10', '_35_5')

for (c in corners){
  for (g in c('g500', 'g700')){
    df_final[[paste0('g300', c, '_', g, c)]] <- df_harm[not_28, paste0('g300', c)] - df_harm[not_28, paste0(g, c)]
  }
}

# lags
aux <- data.frame(
  station = df_harm$station,
  t = df_harm$t,
  l = df_harm$l,
  g300 = df_harm$g300,
  g300_45_.10 = df_harm$g300_45_.10,
  g300_45_5 = df_harm$g300_45_5,
  g300_35_.10 = df_harm$g300_35_.10,
  g300_35_5 = df_harm$g300_35_5,
  g500 = df_harm$g500,
  g500_45_.10 = df_harm$g500_45_.10,
  g500_45_5 = df_harm$g500_45_5,
  g500_35_.10 = df_harm$g500_35_.10,
  g500_35_5 = df_harm$g500_35_5,
  g700 = df_harm$g700,
  g700_45_.10 = df_harm$g700_45_.10,
  g700_45_5 = df_harm$g700_45_5,
  g700_35_.10 = df_harm$g700_35_.10,
  g700_35_5 = df_harm$g700_35_5
)

aux <- aux %>%
  group_by(t, station) %>%
  mutate(across(-c(1),
                .fns = ~lag(.),
                .names = '{.col}_lag')) %>%
  as.data.frame() %>% na.omit()

for (g in c('g300','g500', 'g700')){
  for (c in corners){
    df_final[[paste0(g, c , '_', g, c, '_lag')]] <- aux[[paste0(g, c)]] - aux[[paste0(g, c, '_lag')]]
  }
}

# solo JJA
jun_ag <- which(df_final$l >= 152)
df_jun_ag <- df_final[jun_ag, ] #this dataframe is the one we use to fit the models

#----Step 0: null models and initial models (only harmonics)----
models_null <- list()
for (i in 1:dim(stations)[1]){
  ind <- which(df_jun_ag$station == stations$STAID[i])
  mod_null <- suppressWarnings(
    rq(Y ~ 1, tau = 0.50, data = df_jun_ag, subset = ind)
  )
  models_null[[as.character(stations$STAID[i])]] <- mod_null
}

harmonics <- c('c.1', 's.1')

formula <- as.formula(paste('Y ~', paste(harmonics, collapse = '+')))

models_harmonics <- list()

for (i in 1:dim(stations)[1]){
  ind <- which(df_jun_ag$station == stations$STAID[i])
  mod <- rq(formula, tau = 0.50, data = df_jun_ag, subset = ind)
  mod$R1 <- 1 - mod$rho / models_null[[as.character(stations$STAID[i])]]$rho
  models_harmonics[[as.character(stations$STAID[i])]] <- mod
}

#----Step 1: variables or the air column in each station (8)----
source('CL-BIC.R')
vars_air_column <- c(
  'g300', 'g500', 'g700',
  'g300_g500', 'g300_g700',
  'g300_g300_lag', 'g500_g500_lag', 'g700_g700_lag'
)

formula <- as.formula(
  paste('Y ~', 
        paste(c(harmonics, vars_air_column), collapse = '+')))

mod_step1 <- step_rq_CLBIC(
  initial_models = models_harmonics,
  null_models = models_null,
  data = df_jun_ag,
  stations = stations,
  scope = formula
)
