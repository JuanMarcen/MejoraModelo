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
library(quantreg)
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
  # boot <- boot.rq(x = cbind(1, df[ind, vars]),
  #                   y = df[ind, 'Y'], 
  #                   tau = 0.95)
  # IC <- t(apply(boot$B, 2, quantile, probs = c(0.025, 0.975)))
  # 
  # aux <- data.frame(
  #   station = rep(station, 6),
  #   term  = names(M3.models.q0.95[[station]]$coefficients),
  #   est = as.numeric(M3.models.q0.95[[station]]$coefficients),
  #   lower = IC[, 1],
  #   upper = IC[, 2],
  #   stringsAsFactors = FALSE
  # )
  # 
  # df.ic.q0.95 <- rbind(df.ic.q0.95, aux)
}

saveRDS(df.ic.q0.95, 'df.ic.q0.95.rds')

library(dplyr)
library(ggplot2)
library(viridis)
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

# maps of parameters
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")],
      data = stations,
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    ),
    'sf'
  ),
  2062
)

limits <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-10.2, 5.2), Y = c(34.8, 44)),
      data = data.frame(X = c(-10.2, 5.2), Y = c(34.8, 44)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    ),
    'sf'
  ),
  2062
)

world_map <- ne_countries(scale = "large", returnclass = 'sf')
european_union <- c("Algeria", "Andorra", "France", "Gibraltar", "Morocco", 
                    "Portugal", "Spain")
european_union_map <- world_map %>% filter(name %in% european_union)
background <- st_transform(european_union_map, 2062)

# paramater values as df
map.coef <- function(var){
  aux.df <- stations
  aux.val.q0.50 <- c()
  aux.val.q0.95 <- c()
  for (station in aux.df$NAME2){
    aux.val.q0.50 <- c(aux.val.q0.50, M3.models.q0.50[[station]]$coefficients[var])
    aux.val.q0.95 <- c(aux.val.q0.95, M3.models.q0.95[[station]]$coefficients[var])
  }
  aux.df$varq0.50 <- aux.val.q0.50
  aux.df$varq0.95 <- aux.val.q0.95
  
  #m3 q0.50
  g1 <- ggplot(data = background) +
    geom_sf(fill = "antiquewhite") +
    xlab("Longitude") + 
    ylab("Latitude") + 
    ggtitle(paste0(var, ' coefficient (q0.50)')) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6, angle = 90),
          axis.title=element_text(size = 10, face = "bold")) + 
    #values
    geom_sf(data = aux.df, 
            aes(color = varq0.50, size = varq0.50)) + 
    scale_color_viridis(name = var,
                        option = 'viridis',
                        limits = c(min(aux.df$varq0.50),max(aux.df$varq0.50)),
                        direction = 1) +
    scale_size_continuous(
      range = c(1, 3),   # tamaño mínimo y máximo
      guide = "none"     # quita la leyenda del tamaño
    ) +
    # limits
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) 
  
  # m3 q0.95
  g2 <- ggplot(data = background) +
    geom_sf(fill = "antiquewhite") +
    xlab("Longitude") + 
    ylab("Latitude") + 
    ggtitle(paste0(var, ' coefficient (q0.95)')) + 
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6, angle = 90),
          axis.title=element_text(size = 10, face = "bold")) + 
    #values
    geom_sf(data = aux.df, 
            aes(color = varq0.95, size = varq0.95)) + 
    scale_color_viridis(name = var,
                        option = 'viridis',
                        limits = c(min(aux.df$varq0.95),max(aux.df$varq0.95)),
                        direction = 1) +
    scale_size_continuous(
      range = c(1, 3),   # tamaño mínimo y máximo
      guide = "none"     # quita la leyenda del tamaño
    ) +
    # limits
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) 
  
  g <- ggpubr::ggarrange(g1, g2, ncol = 2)
  return(g)
}

#saving 
for(var in vars){
  g <- map.coef(var)
  filename <- paste0('graphs/local models/maps/', var, '.pdf')
  ggsave(filename,
         g, height = 7*0.75, width = 16*0.75)
}

