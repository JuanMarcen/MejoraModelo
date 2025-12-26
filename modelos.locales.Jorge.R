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


#######################################
# NEW LOCAL MODELS
# M3
formula <- as.formula('Y ~ s.1 + c.1 + g300 + g500 + g700')
vars <- c('s.1', 'c.1', 'g300', 'g500', 'g700')

M3.models.q0.50 <- list()
M3.models.q0.95 <- list()
df.ic.q0.50 <- data.frame(
  station = NA,
  term = NA,
  est = NA,
  lower = NA,
  upper = NA
)
df.ic.q0.95 <- data.frame(
  station = NA,
  term = NA,
  est = NA,
  lower = NA,
  upper = NA
)
for (station in stations$NAME2){
  ind <- which(df$station == stations$STAID[stations$NAME2 == station])
  M3.models.q0.50[[station]] <- rq(formula, tau = 0.50,
                                   data = df,
                                   subset = ind)
  M3.models.q0.95[[station]] <- rq(formula, tau = 0.95,
                                   data = df,
                                   subset = ind)
  
  #confidence intervals
  boot <- boot.rq(x = cbind(1, df[ind, vars]),
                    y = df[ind, 'Y'], 
                    tau = 0.95)
  IC <- t(apply(boot$B, 2, quantile, probs = c(0.025, 0.975)))
  
  aux <- data.frame(
    station = rep(station, 6),
    term  = names(M3.models.q0.95[[station]]$coefficients),
    est = as.numeric(M3.models.q0.95[[station]]$coefficients),
    lower = IC[, 1],
    upper = IC[, 2],
    stringsAsFactors = FALSE
  )
  
  df.ic.q0.95 <- rbind(df.ic.q0.95, aux)
}

saveRDS(df.ic.q0.95, 'df.ic.q0.95.rds')

library(dplyr)
library(ggplot2)
dist <- readRDS('dist.full.coast.rds')

#plots CI
for (i in 1:length(vars)){
  basura <- df.ic.q0.95 %>%
    filter(term == vars[i]) %>% 
    mutate(dist = dist) %>%
    arrange(dist) %>%
    mutate(station = factor(station, levels = station))
  
  g <- ggplot(basura, aes(x = station, y = est)) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.2
    ) +
    labs(
      title = paste(vars[i], 'coefficients q0.95'),
      x = "Station (ordered by distance to coast)",
      y = "Estimated value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 1
      )
    )
  filename <- paste0('graphs/local models/', vars[i], '.q0.95.pdf')
  ggsave(filename, g, height = 8, width = 12)
  
}
