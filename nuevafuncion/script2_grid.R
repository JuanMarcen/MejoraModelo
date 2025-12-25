# install.packages(c("sf", "rnaturalearth", "rnaturalearthdata"))

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

esp_prov <- ne_states(country = "Spain", returnclass = "sf")

provincias_aragon <- esp_prov[esp_prov$name %in% c("Navarra", "Huesca", "Zaragoza", "Teruel"), ]

aragon <- st_union(provincias_aragon)

aragon <- st_sf(geometry = aragon)

#aragon_utm <- st_transform(aragon, 25830)

cellsize <- .1  # ~ 10 km

grid_aragon <- st_make_grid(
  aragon, #_utm
  cellsize = cellsize,
  square = TRUE
)

grid_aragon <- st_sf(geometry = grid_aragon)

grid_aragon_clip <- st_intersection(grid_aragon, aragon) # _utm

plot(st_geometry(grid_aragon_clip), border = "grey")

centroides <- st_centroid(grid_aragon_clip)
centroides$X <- st_coordinates(centroides)[,1]
centroides$Y <- st_coordinates(centroides)[,2]

plot_grid_with_significance <- function(
    grid_sf, values, lower, upper,
    title = "Map with significance",
    value_range = NULL,
    nonsig_method = c("hatch", "shade")
) {
  
  library(sf)
  library(ggplot2)
  library(dplyr)
  
  nonsig_method <- match.arg(nonsig_method)
  
  if (!inherits(grid_sf, "sf"))
    stop("grid_sf must be an sf object")
  if (any(lengths(list(values, lower, upper)) != nrow(grid_sf)))
    stop("values, lower, and upper must match number of grid cells")
  
  nonsig <- (lower <= 0 & upper >= 0)
  
  grid_sf$estimate <- values
  grid_sf$nonsig  <- nonsig
  
  p <- ggplot(grid_sf) +
    geom_sf(aes(fill = estimate), color = NA)
  
  if (is.null(value_range)) {
    p <- p + scale_fill_viridis_c(option = "plasma")
  } else {
    p <- p + scale_fill_viridis_c(option = "plasma",
                                  limits = value_range,
                                  oob = scales::squish)
  }
  
  if (nonsig_method == "hatch") {
    p <- p +
      geom_sf(data = subset(grid_sf, nonsig),
              fill = NA, color = "white", size = 0.2)
  }
  
  if (nonsig_method == "shade") {
    p <- p +
      geom_sf(data = subset(grid_sf, nonsig),
              fill = "grey70", alpha = 0.5, color = NA)
  }
  
  p +
    labs(fill = "Estimate", title = title) +
    theme_minimal()
}

### fit models for quantiles 05, 50, 95
install.packages("spTReg")
if (!require("remotes")) install.packages("remotes")
remotes::install_github("JorgeCastilloMateo/spTReg")
library(lubridate)
library(spTReg)

modelp05 <- spTReg::spTm(
  Tmax ~ -1 +
    scale(sin(2*pi*yday/365)) + scale(cos(2*pi*yday/365)), 
  u = matrix(0, 2, 2), v = cbind(1, scale(df0$year)),
  x_alpha = list(
    x1 = cbind(1, scale(coords$altitude)),
    x2 = matrix(1, nrow = 15, ncol = 1)),
  coords = as.matrix(coords[,c("lat", "lon")]),
  dates = as.Date("2025-12-01"),
  priors = list("beta" = c(0, 1 / 10^16), 
                "sigma" = c(0.1, 0.1),
                "phi" = c(2, 10),
                "mu" = c(0, 1 / 10^16)),
  starting = list(beta = 0.01, sigma = 1, gamma = matrix(0, 2, 2), alpha = 0, hp = c(0, 1, 3/100)),
  data = df0, method = "quantile", quantile = 0.05, na.action = na.pass,
  n.samples = 1000, n.thin = 1, n.burnin = 1000, verbose = TRUE,
  n.report = 10)

modelp50 <- spTReg::spTm(
  Tmax ~ -1 +
    scale(sin(2*pi*yday/365)) + scale(cos(2*pi*yday/365)), 
  u = matrix(0, 2, 2), v = cbind(1, scale(df0$year)),
  x_alpha = list(
    x1 = cbind(1, scale(coords$altitude)),
    x2 = matrix(1, nrow = 15, ncol = 1)),
  coords = as.matrix(coords[,c("lat", "lon")]),
  dates = as.Date("2025-12-01"),
  priors = list("beta" = c(0, 1 / 10^16), 
                "sigma" = c(0.1, 0.1),
                "phi" = c(2, 10),
                "mu" = c(0, 1 / 10^16)),
  starting = list(beta = 0.01, sigma = 1, gamma = matrix(0, 2, 2), alpha = 0, hp = c(0, 1, 3/100)),
  data = df0, method = "quantile", quantile = 0.50, na.action = na.pass,
  n.samples = 1000, n.thin = 1, n.burnin = 1000, verbose = TRUE,
  n.report = 100)

modelp95 <- spTReg::spTm(
  Tmax ~ -1 +
    scale(sin(2*pi*yday/365)) + scale(cos(2*pi*yday/365)), 
  u = matrix(0, 2, 2), v = cbind(1, scale(df0$year)),
  x_alpha = list(
    x1 = cbind(1, scale(coords$altitude)),
    x2 = matrix(1, nrow = 15, ncol = 1)),
  coords = as.matrix(coords[,c("lat", "lon")]),
  dates = as.Date("2025-12-01"),
  priors = list("beta" = c(0, 1 / 10^16), 
                "sigma" = c(0.1, 0.1),
                "phi" = c(2, 10),
                "mu" = c(0, 1 / 10^16)),
  starting = list(beta = 0.01, sigma = 1, gamma = matrix(0, 2, 2), alpha = 0, hp = c(0, 1, 3/100)),
  data = df0, method = "quantile", quantile = 0.95, na.action = na.pass,
  n.samples = 1000, n.thin = 1, n.burnin = 1000, verbose = TRUE,
  n.report = 100)

### obtain kriged values
krigep05.1 <- spTReg::krigeBayes(modelp05$p.params.samples[,1:15+21],
                                 modelp05$p.params.samples[,36+c(1,2,3)],
                                 coords = as.matrix(coords[,c("lon", "lat")]),
                                 newcoords = st_coordinates(centroides),
                                 parallel = TRUE)

krigep50.1 <- spTReg::krigeBayes(modelp50$p.params.samples[,1:15+21],
                                 modelp50$p.params.samples[,36+c(1,2,3)],
                                 coords = as.matrix(coords[,c("lon", "lat")]),
                                 newcoords = st_coordinates(centroides),
                                 parallel = TRUE)

krigep95.1 <- spTReg::krigeBayes(modelp95$p.params.samples[,1:15+21],
                                 modelp95$p.params.samples[,36+c(1,2,3)],
                                 coords = as.matrix(coords[,c("lon", "lat")]),
                                 newcoords = st_coordinates(centroides),
                                 parallel = TRUE)

### plots
plot_grid_with_significance(
  grid_aragon_clip, colMeans(krigep05.1) / sd(df0$year),
  lower = apply(krigep05.1, 2, quantile, prob = c(0.05)) / sd(df0$year),
  upper = apply(krigep05.1, 2, quantile, prob = c(0.95)) / sd(df0$year),
  value_range = c(-0.2,1.2), nonsig_method = 'shade') + geom_point(aes(x=lon, y=lat), data = coords[,c("lon", "lat")])
plot_grid_with_significance(
  grid_aragon_clip, colMeans(krigep50.1) / sd(df0$year),
  lower = apply(krigep50.1, 2, quantile, prob = c(0.05)) / sd(df0$year),
  upper = apply(krigep50.1, 2, quantile, prob = c(0.95)) / sd(df0$year),
  value_range = c(-0.2,1.2), nonsig_method = 'shade') + geom_point(aes(x=lon, y=lat), data = coords[,c("lon", "lat")])
plot_grid_with_significance(
  grid_aragon_clip, colMeans(krigep95.1) / sd(df0$year),
  lower = apply(krigep95.1, 2, quantile, prob = c(0.05)) / sd(df0$year),
  upper = apply(krigep95.1, 2, quantile, prob = c(0.95)) / sd(df0$year),
  value_range = c(-0.2,1.2),  nonsig_method = 'shade') + geom_point(aes(x=lon, y=lat), data = coords[,c("lon", "lat")])
