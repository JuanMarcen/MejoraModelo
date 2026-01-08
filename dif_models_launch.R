rm(list = ls())

# data
library(qs)
library(spTReg)
library(quantreg)
df <- qread('df_jun_ag.qs')
stations <- readRDS('stations.rds')

source('metrics_bay.R')

# ---- 1. IIDM MODEL ----
formula <- as.formula('Y ~ s.1 + c.1 + g300 + g500 + g700')

mod1.iidm <- iidm(
  formula = formula,
  data = df,
  method = 'quantile',
  quantile = 0.5,
  priors = list(
    beta = c(0, 0.0001),
    sigma = c(0.1, 0.1)
  ),
  n.samples = 1000,
  n.thin = 1,
  n.burnin = 1000,
  n.report = 10
)

plot(mod1.iidm$p.params.samples)

params.df <- as.data.frame(mod1.iidm$p.params.samples)
betas <- apply(params.df[-ncol(params.df)], 2, mean)

pred <- numeric(nrow(df))
for (i in 1:length(pred)){
  pred[i] <- betas[1] + betas[2] * df[i, 's.1'] + betas[3] * df[i, 'c.1'] + 
    betas[4] * df[i, 'g300'] + betas[5] * df[i, 'g500'] + betas[6] * df[i, 'g700']
}
pred<- cbind(df[,c('Date', 'station', 'Y')], pred)
colnames(pred)[4] <- 'pred_q0.50'

R1.mod1.iidm <- R1_bay(pred, 0.50, df)
R1.mod1.iidm$R1_locales
R1.mod1.iidm$R1_globales

mod_nulo_f <- rq(Y ~ as.factor(station), tau = 0.50, data = df)
mod_nulo <- rq(Y ~ 1, tau = 0.50, data = df)

mod_global <- rq(formula, data = df, tau = 0.50)

quantile(df$Y, probs = 0.50)
1 - mod_global$rho/mod_nulo$rho

betas
mod_global$coefficients

# ---- 2. MODEL WITH GP (NO COAST EFFECT) (M2) ----
load('fullmodels.RData')
load('M2.RData')

par(mfrow = c(4,5))
for (i in 1:40){
  plot(final.chain.q0.50.M2[, paste0('beta6(s', i, ')')] + 
         final.chain.q0.50.M2[, 'g700'], type = 'l',
       main = stations$NAME2[i])
}

# ---- 3. MODEL WITH GP (COASTAL EFFECT) ----

# ---- 4. MODEL WITH GP (COASTAL CONVOLUTION) ----
conv.q0.50.100k <- readRDS('conv.q0.50.100k.rds')
