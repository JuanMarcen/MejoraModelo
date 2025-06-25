# Improving existing models
# example in one station

rm(list=ls())

# Data uploading
df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
vars_q0.5 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/vars_q0.5.rds")
vars_q0.95 <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/vars_q0.95.rds")
stations <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/stations.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

# Subset Madrid (Retiro)
ind <- which(Y$station==230)

# Pruebas
vars_q0.5 <- ifelse(grepl("^I\\(.*\\)$", vars_q0.5),
                    paste0("`", vars_q0.5, "`"),
                    vars_q0.5)
vars_q0.95 <- ifelse(grepl("^I\\(.*\\)$", vars_q0.95),
                     paste0("`", vars_q0.95, "`"),
                     vars_q0.95)

formula_q0.5 <- as.formula(
  paste('Y~', paste(
    vars_q0.5,collapse = '+'
  ),'+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))')
)
formula_q0.95 <- as.formula(
  paste('Y~', paste(
    vars_q0.95,collapse = '+'
  ),'+ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))')
)


# Modelos
library(quantreg)
mod_q0.5<-rq(formula_q0.5,data=df_conj_filled_sc,subset=ind,tau=0.50)
mod_q0.95<-rq(formula_q0.95,data=df_conj_filled_sc,subset=ind,tau=0.95)