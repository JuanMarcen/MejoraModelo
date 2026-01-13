rm(list = ls())

# data
library(qs)
library(spTReg)
library(quantreg)
df <- qread('df_jun_ag.qs')
library(dplyr)
library(lubridate)
df <- df %>%
  select(Date, station, Y, l, t, s.1, c.1, elev, dist, g300, g500, g700) %>%
  mutate(
    month = month(Date),
    `t:month6` = ifelse(month == 6, t, 0),
    `t:month7` = ifelse(month == 7, t, 0),
    `t:month8` = ifelse(month == 8, t, 0)
  )
stations <- readRDS('stations.rds')

source('metrics_bay.R')

# ---- 1. IIDM MODEL ----
formula <- as.formula('Y ~ s.1 + c.1 + g300 * (elev + dist) + 
                      g500 * (elev + dist) + g700 * (elev + dist) + 
                      `t:month6`*(elev + dist) + `t:month7`*(elev + dist) +
                      `t:month8`*(elev + dist)')

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
    betas[12] * df[i, 'elev'] * df[i, 'g300'] + betas[13] * df[i, 'dist'] * df[i, 'g300']+ 
    betas[14] * df[i, 'elev'] * df[i, 'g500'] + betas[15] * df[i, 'dist'] * df[i, 'g500']+ 
    betas[16] * df[i, 'elev'] * df[i, 'g700'] + betas[17] * df[i, 'dist'] * df[i, 'g700'] +
    
    betas[9] * df[i, 't:month6'] + betas[10] * df[i, 't:month7'] + betas[11] * df[i, 't:month8'] + 
    betas[18] * df[i, 'elev'] * df[i, 't:month6'] + betas[19] * df[i, 'dist'] * df[i, 't:month6']+ 
    betas[20] * df[i, 'elev'] * df[i, 't:month7'] + betas[21] * df[i, 'dist'] * df[i, 't:month7']+ 
    betas[22] * df[i, 'elev'] * df[i, 't:month8'] + betas[23] * df[i, 'dist'] * df[i, 't:month8'] 
}
pred<- cbind(df[,c('Date', 'station', 'Y')], pred)
colnames(pred)[4] <- 'pred_q0.50'

R1.mod1.iidm <- R1_bay(pred, 0.50, df)
R1.mod1.iidm$R1_locales
R1.mod1.iidm$R1_globales


rho.mod1.iidm <- rho_bay(pred, 0.50)
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
coastal.q0.50.100k <- readRDS('coastal.trend.q0.50.100k.rds')

betas <- matrix(apply(coastal.q0.50.100k$process, 2, mean))

pred.coastal <- numeric(nrow(df))
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  X <- as.matrix(cbind(1, df[ind, c('s.1', 'c.1', 'g300', 'g500', 'g700')]))
  betas.aux <- betas[c(i, i + 48, i + 48 * 2, i + 48 * 3, i + 48 * 4, i + 48 * 5), ,drop = FALSE]
  pred.coastal[ind] <- X %*% betas.aux
}

pred.coastal <- cbind(df[,c('Date', 'station', 'Y')], coastal.q0.50.100k$fitted)
colnames(pred.coastal)[4] <- 'pred_q0.50'

R1.coastal <- R1_bay(pred.coastal, 0.50, df)
R1.coastal$R1_locales
R1.coastal$R1_globales

rho.coastal <- rho_bay(pred.coastal, 0.50)

# ---- 4. MODEL WITH GP (COASTAL CONVOLUTION) ----
conv.q0.50.100k <- readRDS('conv.nofactor.q0.50.100k.rds')
conv.q0.50.100k <- readRDS('conv.trend.q0.50.100k.rds')
betas <- matrix(apply(conv.q0.50.100k$process, 2, mean))

pred.conv <- numeric(nrow(df))
# i could use the fitted component of the list
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  X <- as.matrix(cbind(1, df[ind, c('s.1', 'c.1', 'g300', 'g500', 'g700', 
                                    't:month6', 't:month7', 't:month8')]))
  betas.aux <- betas[c(i, i + 48, i + 48 * 2, i + 48 * 3, i + 48 * 4, i + 48 * 5,
                       i + 48 * 6, i + 48 * 7, i + 48 * 8), ,drop = FALSE]
  pred.conv[ind] <- X %*% betas.aux
}

pred.conv <- cbind(df[,c('Date', 'station', 'Y')], pred.conv)
colnames(pred.conv)[4] <- 'pred_q0.50'

R1.conv<- R1_bay(pred.conv, 0.50, df)
R1.conv$R1_locales
R1.conv$R1_globales

rho.conv <- rho_bay(pred.conv, 0.50)



  

# ---- fast comparison ----
R1 <- data.frame(iidm = R1.mod1.iidm$R1_locales,
                 og = R1.bay.q0.50.M3$R1_locales, 
                 coastal = R1.coastal$R1_locales, 
                 conv = R1.conv$R1_locales)
names(R1) <- c('iidm','no coast', 'coastal', 'conv')
R1$station <- rownames(R1)


library(ggplot2)
library(tidyr)
aux <- cbind(rho.conv$rho_dias, rho.coastal$rho_dias, rho.q0.50.M3$rho_dias, 
             rho.mod1.iidm$rho_dias)
aux$year <- 1:92
names(aux)[1:4] <- c('conv', 'coastal', 'og', 'iidm')

aux.long <- aux %>%
  pivot_longer(cols = - year, names_to = 'Model', values_to = 'rho')

ggplot(data = aux.long, aes(x = year, y = rho, color = Model)) +
  geom_line() +
  labs(
    x = 't',
    y = 'rho',
    title = 'rho'
  ) +
  ylim(0,1) +
  geom_hline(yintercept = 0.5, color = 'red') +
  theme_minimal(base_size = 12)


plot(1:92, rho.mod1.iidm$rho_dias_est$`Madrid (Barajas)`[, 1], type = 'l')
abline(h = 0.5, col= 'red')
# ---- DIC ----
# DIC of the original one (i need the fcking XI)
N <- nrow(df)

# ---- inference
