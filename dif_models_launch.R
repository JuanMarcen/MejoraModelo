rm(list = ls())

# data
library(qs)
library(spTReg)
library(quantreg)
df <- qread('df_jun_ag.qs')
stations <- readRDS('stations.rds')

source('metrics_bay.R')

# ---- 1. IIDM MODEL ----
formula <- as.formula('Y ~ s.1 + c.1 + g300 * (elev + dist) + 
                      g500 * (elev + dist) + g700 * (elev + dist)')

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
  n.report = 100
)

plot(mod1.iidm$p.params.samples)

params.df <- as.data.frame(mod1.iidm$p.params.samples)
betas <- apply(params.df[-ncol(params.df)], 2, mean)

pred <- numeric(nrow(df))
for (i in 1:length(pred)){
  pred[i] <- betas[1] + betas[2] * df[i, 's.1'] + betas[3] * df[i, 'c.1'] + 
    betas[4] * df[i, 'g300'] + betas[7] * df[i, 'g500'] + betas[8] * df[i, 'g700'] + 
    betas[5] * df[i, 'elev'] + betas[6] * df[i, 'dist'] + 
    betas[9] * df[i, 'elev'] * df[i, 'g300'] + betas[10] * df[i, 'dist'] * df[i, 'g300']+ 
    betas[11] * df[i, 'elev'] * df[i, 'g500'] + betas[12] * df[i, 'dist'] * df[i, 'g500']+ 
    betas[13] * df[i, 'elev'] * df[i, 'g700'] + betas[14] * df[i, 'dist'] * df[i, 'g700']
}
pred<- cbind(df[,c('Date', 'station', 'Y')], pred)
colnames(pred)[4] <- 'pred_q0.50'

R1.mod1.iidm <- R1_bay(pred, 0.50, df)
R1.mod1.iidm$R1_locales
R1.mod1.iidm$R1_globales

# mod_nulo_f <- rq(Y ~ as.factor(station), tau = 0.50, data = df)
# mod_nulo <- rq(Y ~ 1, tau = 0.50, data = df)
# 
# mod_global <- rq(formula, data = df, tau = 0.50)
# 
# quantile(df$Y, probs = 0.50)
# 1 - mod_global$rho/mod_nulo$rho
# 
# betas
# mod_global$coefficients
# 
# #chequeo 
# pred$pred_q0.50 <- mod_global$fitted.values
# R1.mod1.iidm <- R1_bay(pred, 0.50, df)
# R1.mod1.iidm$R1_locales
# R1.mod1.iidm$R1_globales
# ---- 2. MODEL WITH GP (NO COAST EFFECT) (M2) ----
load('fullmodels.RData')
load('M2.RData')

cbind(R1.mod1.iidm$R1_locales, R1.bay.q0.50.M2$R1_locales)

par(mfrow = c(4,5))
for (i in 1:40){
  plot(final.chain.q0.50.M2[, paste0('beta6(s', i, ')')] + 
         final.chain.q0.50.M2[, 'g700'], type = 'l',
       main = stations$NAME2[i])
}

# ---- 3. MODEL WITH GP (COASTAL EFFECT) ----
coastal.q0.50.100k <- readRDS('coastal.q0.50.100k.rds')

betas <- matrix(apply(coastal.q0.50.100k$process, 2, mean))

pred.coastal <- numeric(nrow(df))
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  X <- as.matrix(cbind(1, df[ind, c('s.1', 'c.1', 'g300', 'g500', 'g700')]))
  betas.aux <- betas[c(i, i + 48, i + 48 * 2, i + 48 * 3, i + 48 * 4, i + 48 * 5), ,drop = FALSE]
  pred.coastal[ind] <- X %*% betas.aux
}

pred.coastal <- cbind(df[,c('Date', 'station', 'Y')], pred.coastal)
colnames(pred.coastal)[4] <- 'pred_q0.50'

R1.coastal <- R1_bay(pred.coastal, 0.50, df)
R1.coastal$R1_locales
R1.coastal$R1_globales

# ---- 4. MODEL WITH GP (COASTAL CONVOLUTION) ----
conv.q0.50.100k <- readRDS('conv.nofactor.q0.50.100k.rds')
betas <- matrix(apply(conv.q0.50.100k$process, 2, mean))

pred.conv <- numeric(nrow(df))
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  X <- as.matrix(cbind(1, df[ind, c('s.1', 'c.1', 'g300', 'g500', 'g700')]))
  betas.aux <- betas[c(i, i + 48, i + 48 * 2, i + 48 * 3, i + 48 * 4, i + 48 * 5), ,drop = FALSE]
  pred.conv[ind] <- X %*% betas.aux
}

pred.conv <- cbind(df[,c('Date', 'station', 'Y')], pred.conv)
colnames(pred.conv)[4] <- 'pred_q0.50'

R1.conv<- R1_bay(pred.conv, 0.50, df)
R1.conv$R1_locales
R1.conv$R1_globales

# ---- fast comparison ----
R1 <- data.frame(og = R1.bay.q0.50.M2$R1_locales, coastal = R1.coastal$R1_locales, conv = R1.conv$R1_locales)
names(R1) <- c('no coast', 'coastal', 'conv')
R1$station <- rownames(R1)

library(dplyr)
library(tidyr)
R1 <- R1 %>%
  pivot_longer(cols = -station, names_to = 'model', values_to = 'R1')

# ---- DIC ----

loglik_values <- apply(conv.q0.50.100k, 1, logLik, y = df$Y)

D_theta <- -2 * loglik_values
D_bar <- mean(D_theta)

theta_bar <- colMeans(theta)
D_theta_bar <- -2 * loglik(theta_bar, y)

DIC <- 2 * D_bar - D_theta_bar
DIC
