rm(list=ls())

df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

ind <- which(Y$station==230)

df_madrid<-df_conj_filled_sc[ind,]

library(corrplot)

corrplot(cor(df_madrid[2:46]),
         order = 'hclust',
         method = 'color')

# Annual trend of a geopotencial
library(dplyr)

for (i in 2:46){
  var <- colnames(df_madrid)[i]
  
  media <- df_madrid %>%
    group_by(t) %>%
    summarise(media = mean(.data[[var]], na.rm = T)) %>%
    as.data.frame()
  
  plot(media[, 't'], media[, 'media'], type = 'b', pch = 19, main = var)
}

# Anual trend of temperature
media <- df_madrid %>%
  group_by(t) %>%
  summarise(media = mean(Y, na.rm = T)) %>%
  as.data.frame()

plot(media[, 't'], media[, 'media'], type = 'b', pch = 19, main = 'Y')


# Empirical Quantile
quantile <- df_madrid %>%
  group_by(t) %>%
  summarise(q0.5 = quantile(Y, na.rm = T, probs = 0.5), q0.95 = quantile(Y, na.rm = T, probs = 0.95)) %>%
  ungroup() %>%
  as.data.frame()

plot(quantile[, 't'], quantile[, 'q0.5'], type = 'b', 
     pch = 19, main = 'Quantiles', ylim = c(25,40))
lines(quantile[, 't'], quantile[, 'q0.95'], type = 'b', col = 'red', pch = 19)


# Quantile model with only t
library(quantreg)

mod_q0.5 <- rq(Y ~ t, tau = 0.5, data = df_madrid)
mod_q0.95 <- rq(Y ~ t, tau = 0.95, data = df_madrid)

plot(df_madrid$t, df_madrid$Y, pch = 19)
lines(df_madrid$t, predict(mod_q0.5, df_madrid), col='blue', lwd= 2)
lines(df_madrid$t, predict(mod_q0.95, df_madrid), col='red', lwd= 2)


# Linear models of the empirical quantiles with t
mod1 <- lm(q0.5 ~ t, data = quantile)
mod2 <- lm(q0.95 ~ t, data = quantile)

lines(quantile$t, predict(mod1, quantile), col='green', lwd= 2)
lines(quantile$t, predict(mod2, quantile), col='purple', lwd= 2)

# CV 
library(boot)
miNIP <- 815607
maxdegree<-20
eep_q0.5 <- rep(NA,maxdegree)
eep_q0.95 <- rep(NA,maxdegree)
for (grado in 1:maxdegree) {
  # q0.5
  set.seed(miNIP)
  eep_q0.5[grado] <- cv.glm(
    data = quantile,
    glm(q0.5 ~ poly(t, grado), data = quantile),
    K = 10
  )$delta[2]
  
  # q0.95
  set.seed(miNIP)
  eep_q0.95[grado] <- cv.glm(
    data = quantile,
    glm(q0.95 ~ poly(t, grado), data = quantile),
    K = 10
  )$delta[2]
}

plot(eep_q0.5, type = 'b')
gr.q0.5 <- which.min(eep_q0.5)
gr.eep.q0.5 <- min(eep_q0.5)
points(gr.q0.5, gr.eep.q0.5, col='red', pch=19)

plot(eep_q0.95, type = 'b')
gr.q0.95 <- which.min(eep_q0.95)
gr.eep.q0.95 <- min(eep_q0.95)
points(gr.q0.95, gr.eep.q0.95, col='red', pch=19)


mod3 <- lm(q0.5 ~ poly(t, 7), data = quantile)
mod4 <- lm(q0.95 ~ poly(t, 7), data = quantile)

lines(quantile$t, predict(mod3, quantile), col='darkgreen', lwd= 2)
lines(quantile$t, predict(mod4, quantile), col='magenta', lwd= 2)

