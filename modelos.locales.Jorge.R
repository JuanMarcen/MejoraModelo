rm(list = ls())

library(qs)
df <- qread('df_jun_ag.qs')
stations <- readRDS('stations.rds')
orden <- readRDS('orden.rds')

# modelos nulos
library(quantreg)
null.models.q0.50 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  null.models.q0.50[[i]] <- rq(Y ~ 1, data = df, subset = ind, tau = 0.50)
}

null.models.q0.95 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  null.models.q0.95[[i]] <- rq(Y ~ 1, data = df, subset = ind, tau = 0.95)
}

#M1
formula <- as.formula('Y ~  s.1 + c.1 + g300 + g500 + g700')
M1.q0.50 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  M1.q0.50[[i]] <- rq(formula, data = df, subset = ind, tau = 0.50)
}

M1.q0.95 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  M1.q0.95[[i]] <- rq(formula, data = df, subset = ind, tau = 0.95)
}

R1.M1.q0.50 <- c()
R1.M1.q0.95 <- c()
for (i in 1:40){
  R1.M1.q0.50 <- c(R1.M1.q0.50, 1 - M1.q0.50[[i]]$rho / null.models.q0.50[[i]]$rho)
  R1.M1.q0.95 <- c(R1.M1.q0.95, 1 - M1.q0.95[[i]]$rho / null.models.q0.95[[i]]$rho)
}

#M2
formula <- as.formula('Y ~ s.1 + c.1 + g300 + g300_45_.10 + g300_45_5 + g300_35_.10 + g300_35_5 + g500 + g500_45_.10 + g500_45_5 + g500_35_.10 + g500_35_5+ g700 + g700_45_.10 + g700_45_5 + g700_35_.10 + g700_35_5')
M2.q0.50 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  M2.q0.50[[i]] <- rq(formula, data = df, subset = ind, tau = 0.50)
}

M2.q0.95 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  M2.q0.95[[i]] <- rq(formula, data = df, subset = ind, tau = 0.95)
}

R1.M2.q0.50 <- c()
R1.M2.q0.95 <- c()
for (i in 1:40){
  R1.M2.q0.50 <- c(R1.M2.q0.50, 1 - M2.q0.50[[i]]$rho / null.models.q0.50[[i]]$rho)
  R1.M2.q0.95 <- c(R1.M2.q0.95, 1 - M2.q0.95[[i]]$rho / null.models.q0.95[[i]]$rho)
}

#M3
formula_q0.95 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/formula_q0.95.rds")
formula_q0.5 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/formula_q0.5.rds")

formula_q0.5 <- update(formula_q0.5, .~. + s.1 + c.1)
formula_q0.95 <- update(formula_q0.95, .~. + s.1 + c.1)

M3.q0.50 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  M3.q0.50[[i]] <- rq(formula_q0.5, data = df, subset = ind, tau = 0.50)
}

M3.q0.95 <- list()
for (i in 1:40){
  ind <- which(df$station == stations$STAID[i])
  
  M3.q0.95[[i]] <- rq(formula_q0.95, data = df, subset = ind, tau = 0.95)
}

R1.M3.q0.50 <- c()
R1.M3.q0.95 <- c()
for (i in 1:40){
  R1.M3.q0.50 <- c(R1.M3.q0.50, 1 - M3.q0.50[[i]]$rho / null.models.q0.50[[i]]$rho)
  R1.M3.q0.95 <- c(R1.M3.q0.95, 1 - M3.q0.95[[i]]$rho / null.models.q0.95[[i]]$rho)
}



# R1
R1.df <- data.frame(
  R1.M1.q0.50 = round(R1.M1.q0.50, 3),
  R1.M2.q0.50 = round(R1.M2.q0.50, 3),
  R1.M3.q0.50 = round(R1.M3.q0.50, 3),
  R1.M1.q0.95 = round(R1.M1.q0.95, 3),
  R1.M2.q0.95 = round(R1.M2.q0.95, 3),
  R1.M3.q0.95 = round(R1.M3.q0.95, 3)
)
rownames(R1.df) <- stations$NAME2

R1.df <- R1.df[orden, ]

esc_tabla_negrita(R1.df, colq0.5 = c(1,2,3), colq0.95 = c(4,5,6), negrita = T)
