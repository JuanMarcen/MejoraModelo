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
not_28 <- which(df_harm$l >= 149)
df_final <- df_harm[not_28, ]

corners <- c('','_45_.10', '_45_5', '_35_.10', '_35_5')

for (c in corners){
  for (g in c('g500', 'g700')){
    df_final[[paste0('g300', c, '_', g, c)]] <- df_harm[not_28, paste0('g300', c)] - df_harm[not_28, paste0(g, c)]
    df_final[[paste0('g300', c, '_', g, c, '.lag')]] <- lag(df_final[[paste0('g300', c, '_', g, c)]])
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
                .names = '{.col}.lag')) %>%
  as.data.frame() %>% na.omit()

for (g in c('g300','g500', 'g700')){
  for (c in corners){
    df_final[[paste0(g, c, '.lag')]] <- aux[[paste0(g, c, '.lag')]]
    df_final[[paste0(g, c , '_', g, c, '_lag')]] <- aux[[paste0(g, c)]] - aux[[paste0(g, c, '.lag')]]
    df_final[[paste0(g, c , '_', g, c, '_lag.lag')]] <- lag(df_final[[paste0(g, c , '_', g, c, '_lag')]])
  }
}

# solo JJA
jun_ag <- which(df_final$l >= 152)
df_jun_ag <- df_final[jun_ag, ] #this dataframe is the one we use to fit the models

rm(list = setdiff(ls(), c('df_jun_ag', 'stations')))

#----Step 0: null models and initial models (only harmonics)----
models_null <- list()
for (i in 1:dim(stations)[1]){
  ind <- which(df_jun_ag$station == stations$STAID[i])
  mod_null <- suppressWarnings(
    rq(Y ~ 1, tau = 0.95, data = df_jun_ag, subset = ind)
  )
  models_null[[as.character(stations$STAID[i])]] <- mod_null
}

harmonics <- c('c.1', 's.1')

formula <- as.formula(paste('Y ~', paste(harmonics, collapse = '+')))

models_harmonics <- list()

for (i in 1:dim(stations)[1]){
  ind <- which(df_jun_ag$station == stations$STAID[i])
  mod <- rq(formula, tau = 0.95, data = df_jun_ag, subset = ind)
  mod$R1 <- 1 - mod$rho / models_null[[as.character(stations$STAID[i])]]$rho
  models_harmonics[[as.character(stations$STAID[i])]] <- mod
}

rm(list = setdiff(ls(), c('df_jun_ag', 'stations', 'models_harmonics', 'models_null', 'harmonics')))
#----Step 1: variables or the air column in each station (8)----
source('CL-BIC.R')


# mod_step1 <- step_rq_CLBIC(
#   initial_models = models_harmonics,
#   null_models = models_null,
#   data = df_jun_ag,
#   stations = stations,
#   scope = formula,
#   weights = 1,
#   eff_param = F
# )

# with different weights
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")], 
      data = stations[c("STAID", "STANAME", "LON", "LAT", "HGHT","color",'NAME1','NAME2')],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

coords_km <- st_coordinates(stations) / 1000
dist <- as.matrix(dist(coords_km))
w <- exp_weights(dist, h = 100, scale = TRUE) # h will be the CV parameter
names(w) <- stations$NAME2

#with effective number of parameters takes a while (for each CLBIC comp. it has to calculate them)
vars_air_column <- c(
  'g300', 'g500', 'g700',
  'g300_g500', 'g300_g700',
  'g300_g300_lag', 'g500_g500_lag', 'g700_g700_lag'
)

formula <- as.formula(
  paste('Y ~', 
        paste(c(harmonics, vars_air_column), collapse = '+')))
log_file <- 'CLCAIC_sel_q0.95.it.w.txt'
header <- "=== SELECCIÓN MODELO PASO A PASO CLCAIC ===\n"
writeLines(header, log_file)
cat('\n--- Step 1 ---\n\n', file = log_file, append = TRUE)
sink(log_file, append = TRUE)

mod_step1 <- step_rq_CLBIC(
  initial_models = models_harmonics,
  null_models = models_null,
  data = df_jun_ag,
  stations = stations,
  scope = formula,
  weights = w,
  eff_param = TRUE,
  it.weights = T,
  tol = 0.01,
  max.iter = 5,
  replacements = list(
    c('g300','g500','g300_g500'),
    c('g300','g700','g300_g700')
    )
)

sink()

# prep for step 2
# choose variable lag to introduce alone and replacement
aux <- labels(terms(mod_step1$formula))[grepl('g', labels(terms(mod_step1$formula)))]
# aux <- aux[grepl('_lag', aux)]
# aux <- sub("^(.+?)_\\1.*$", "\\1", aux)
vars2 <- paste0(aux, '.lag')

formula <- update(mod_step1$formula, 
                  as.formula(paste(". ~ . +", paste(vars2, collapse = " + ")))
                  )

cat('\n\n--- Step 2 ---\n\n', file = log_file, append = TRUE)
sink(log_file, append = TRUE)
mod_step2 <- step_rq_CLBIC(
  initial_models = mod_step1$models,
  null_models = models_null,
  data = df_jun_ag,
  stations = stations,
  scope = formula,
  weights = w,
  eff_param = T,
  replacements = list(
    c('g500', 'g500.lag', 'g500_g500_lag'),
    c('g300', 'g300.lag', 'g300_g300_lag'),
    c('g700', 'g700.lag', 'g700_g700_lag')
  )
)
sink()

# step 3
aux <- labels(terms(mod_step2$formula))[grepl('g', labels(terms(mod_step2$formula)))]
vars3 <- paste0('I(', aux, '^2)')

formula <- update(mod_step2$formula, 
                  as.formula(paste(". ~ . +", paste(vars3, collapse = " + ")))
)

cat('\n\n--- STEP 3 ---\n\n', file = log_file, append = TRUE)
sink(log_file, append = TRUE)
mod_step3 <- step_rq_CLBIC(
  initial_models = mod_step2$models,
  null_models = models_null,
  data = df_jun_ag,
  stations = stations,
  scope = formula,
  weights = w,
  eff_param = T,
  replacements = list()
)

sink()

# step 4
aux <- labels(terms(mod_step3$formula))[grepl('I', labels(terms(mod_step3$formula)))]
aux <- gsub("I\\((.*)\\^2\\)", "\\1", aux)
vars4 <- paste0('I(', aux, '^3)')

formula <- update(mod_step3$formula, 
                  as.formula(paste(". ~ . +", paste(vars4, collapse = " + ")))
)

cat('\n\n--- STEP 4 ---\n\n', file = log_file, append = TRUE)
sink(log_file, append = TRUE)
mod_step4 <- step_rq_CLBIC(
  initial_models = mod_step3$models,
  null_models = models_null,
  data = df_jun_ag,
  stations = stations,
  scope = formula,
  weights = w,
  eff_param = T,
  replacements = list()
)

sink()

# SAME STEPS BUT WITH ALL THE CORNERS
vars_air_column <- labels(terms(mod_step4$formula))

vars_corners <- c()
corners <- c('_45_.10', '_45_5', '_35_.10', '_35_5')
g <- c('g300', 'g500', 'g700')

for (k in g){
  vars_corners <- c(vars_corners, paste0(k, corners))
}


for (c in corners){
  for (g in c('g500', 'g700')){
    vars_corners <- c(vars_corners, paste0('g300', c, '_', g, c))
  }
}

for (g in c('g300','g500', 'g700')){
  for (c in corners){
    vars_corners <- c(vars_corners, paste0(g, c , '_', g, c, '_lag'))
  }
}

formula <- update(mod_step4$formula, 
                  as.formula(paste(". ~ . +", paste(vars_corners, collapse = " + ")))
)
#replacements (total of 8)
replacement.corners <- list()
i <- 1
for (c in corners){
  for (g in c('g500', 'g700')){
    replacement.corners[[i]] <- c(paste0('g300', c), paste0(g, c),
                                  paste0('g300', c, '_', g, c))
    i <- i + 1
  }
}


cat('\n\n--- STEP 5 (corners) ---\n\n', file = log_file, append = TRUE)
sink(log_file, append = TRUE)
mod_step5 <- step_rq_CLBIC(
  initial_models = mod_step4$models,
  null_models = models_null,
  data = df_jun_ag,
  stations = stations,
  scope = formula,
  weights = w,
  eff_param = T,
  replacements = replacement.corners
)

sink()



# comparison of R1 
# mod1 <- numeric()
# mod2 <- numeric()
# for (i in as.character(stations$STAID)){
#   mod1[i] <- mod_step1[["models"]][[i]][["R1"]]
#   mod2[i] <- mod_step1_2[["models"]][[i]][["R1"]]
# }
# 
# comp_R1 <- cbind(mod1, mod2, mod1 - mod2)



#----simple example of CV for h?----
source('CL-BIC.R')
vars_air_column <- c(
  'g300'
)

formula <- as.formula(
  paste('Y ~', 
        paste(c(harmonics, vars_air_column), collapse = '+')))
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")], 
      data = stations[c("STAID", "STANAME", "LON", "LAT", "HGHT","color",'NAME1','NAME2')],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

coords_km <- st_coordinates(stations) / 1000
dist <- as.matrix(dist(coords_km))
w <- exp_weights(dist, h = mean(dist), scale = TRUE) # h will be the CV parameter
names(w) <- stations$NAME2
# 1. dividir en bloques stations
# ojear el tema de dividir en bloques cercanos de las estaciones y no aleatorios?
create_station_folds <- function(stations, k = 10, seed = 1234) {
  if (!is.vector(stations) && !is.factor(stations)) stop("stations debe ser un vector o factor de ids")
  n <- length(stations)
  if (k < 2) stop("k debe ser al menos 2")
  if (k > n) stop("k no puede ser mayor que el número de estaciones")
  
  
  set.seed(seed)
  # barajar índices y luego cortar en k bloques lo más parejos posible
  idx <- sample(seq_len(n))
  sizes <- rep(floor(n / k), k)
  rem <- n - sum(sizes)
  if (rem > 0) sizes[1:rem] <- sizes[1:rem] + 1
  
  
  folds <- vector("list", k)
  start <- 1
  for (i in seq_len(k)) {
    end <- start + sizes[i] - 1
    folds[[i]] <- stations[idx[start:end]]
    start <- end + 1
  }
  names(folds) <- paste0("Fold", seq_len(k))
  return(folds)
}
get_train_test_from_folds <- function(folds, held_out_fold) {
  k <- length(folds)
  if (held_out_fold < 1 || held_out_fold > k) stop("held_out_fold fuera de rango")
  test <- folds[[held_out_fold]]
  train <- unlist(folds[-held_out_fold], use.names = FALSE)
  return(list(train = train, test = test))
}

K <- 10
folds <- create_station_folds(stations$STAID, k = K)

R1 <- c()
for (fold in seq_along(folds)){
  cat('Fold ', fold, '\n')
  train_test <- get_train_test_from_folds(folds, fold)
  train <- train_test$train
  test <- train_test$test
  
  stations.train <- stations[is.element(stations$STAID, train), ]
  stations.test <- stations[is.element(stations$STAID, test), ]
  coords_km.train <- st_coordinates(stations.train) / 1000
  dist.train <- as.matrix(dist(coords_km.train))
  h <- mean(dist.train) # CV PARAMETER
  w.train <- exp_weights(dist.train, h = h, scale = TRUE) # h will be the CV parameter
  names(w.train) <- stations.train$NAME2
  
  #initial models in those stations too
  models_null.train <- list()
  for (i in 1:dim(stations.train)[1]){
    ind <- which(df_jun_ag$station == stations.train$STAID[i])
    mod_null <- suppressWarnings(
      rq(Y ~ 1, tau = 0.50, data = df_jun_ag, subset = ind)
    )
    models_null.train[[as.character(stations.train$STAID[i])]] <- mod_null
  }
  
  harmonics <- c('c.1', 's.1')
  formula.harm <- as.formula(paste('Y ~', paste(harmonics, collapse = '+')))
  models_harmonics.train <- list()
  for (i in 1:dim(stations.train)[1]){
    ind <- which(df_jun_ag$station == stations.train$STAID[i])
    mod <- rq(formula.harm, tau = 0.50, data = df_jun_ag, subset = ind)
    mod$R1 <- 1 - mod$rho / models_null[[as.character(stations.train$STAID[i])]]$rho
    models_harmonics.train[[as.character(stations.train$STAID[i])]] <- mod
  }
  

  mod_step1.train <- step_rq_CLBIC(
    initial_models = models_harmonics.train,
    null_models = models_null.train,
    data = df_jun_ag,
    stations_df = stations.train,
    scope = formula,
    weights = w.train,
    eff_param = T,
    replacements = list(
      c('g300','g500','g300_g500'),
      c('g300','g700','g300_g700')
    ),
    trace = F
  )
  
  #ajuste de ese modelo en conunto de datos 
  formula.new <- mod_step1.train$formula
  R1.test <- R1_global(formula.new, stations.test, df_jun_ag, 0.5)
  
  R1 <- c(R1, R1.test)
}
