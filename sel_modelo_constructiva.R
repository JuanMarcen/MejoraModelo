rm(list=ls())

df_conj_filled_sc <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/df_conj_filled_sc.rds")
Y <- readRDS("C:/Users/jumar/OneDrive/Escritorio/TFM/Datos/Y.rds")

ind <- which(Y$station==230)

df_madrid <- df_conj_filled_sc[ind,]

source('harmonics.R')
df_harm <- as.data.frame(cbind(df_madrid, cs(df_madrid$t, 1:5)))
colnames(df_harm)[1] <- 'Y'

source('eBIC.R')
mod_nulo_q0.5 <- rq(Y~1, tau= 0.5, data = df_harm)
mod_harm <- step_rq_eBIC(
  mod_nulo_q0.5, data = df_harm,
  scope = as.formula(paste('Y ~',paste0('c.', 1:5, ' + s.', 1:5,collapse = '+')))
) 


mod2<- step_rq_eBIC(
  mod_harm, data = df_harm,
  scope = as.formula(paste('Y~c.1+s.1+',paste(names(df_harm[2:16]),collapse = '+')))
)
