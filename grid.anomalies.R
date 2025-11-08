# grid of anomalies
rm(list = ls())
#mantengo coords_km
library(qs)
library(lubridate)
library(dplyr)

grid_km <- readRDS('grid_km.rds')
coord <- readRDS('coord.rds') #coordinate of ERA5 geop grid points
coords_km <- readRDS('coords_km.rds')
g300_jja <- readRDS("g300_jja.rds")
g500_jja <- readRDS("g500_jja.rds")
g700_jja <- readRDS("g700_jja.rds")
scale.info <- readRDS('scale.info.rds')
grid_elev <- readRDS('grid_elev.rds')
grid_dist <- readRDS('grid_dist.rds')

## nearest neighbor
library(RANN)
nn <- nn2(data = coords_km, query = grid_km, k = 1)
indices_mas_cercanos <- nn$nn.idx[, 1]
estacion_mas_cercana <- coord$station[indices_mas_cercanos]

grid_cercanos <- as.data.frame(cbind(1:790, grid_km, estacion_mas_cercana))
colnames(grid_cercanos) <- c('station_grid','lon','lat','station_2')

#construction of data frame that contain all info (6144*790 rows)
#later filter in our months
# corners not needed
#
order <- match(grid_cercanos$station_2, colnames(g300_jja)[-1]) #a estos indices hay que sumar uno por el de fecha que ignaramos
g300_jja_ord <- g300_jja[, c(1, order + 1)]
Xg300_vector <- as.vector(as.matrix(g300_jja_ord[,-1]))
head(Xg300_vector)
Xg300 <- data.frame(
  Date = rep(g300_jja_ord$Date, times = ncol(g300_jja_ord[,-1])),  # Repetir fechas
  Station = rep(colnames(g300_jja_ord[,-1]), each = nrow(g300_jja_ord[,-1])),  # Repetir estaciones
  Value = Xg300_vector
)
head(Xg300)

#g500 cercano
order<-match(grid_cercanos$station_2, colnames(g500_jja)[-1]) #a estos indices hay que sumar uno por el de fecha que ignaramos
g500_jja_ord <- g500_jja[, c(1, order + 1)]

Xg500_vector <- as.vector(as.matrix(g500_jja_ord[,-1]))
head(Xg500_vector)
Xg500 <- data.frame(
  Date = rep(g500_jja_ord$Date, times = ncol(g500_jja_ord[,-1])),  # Repetir fechas
  Station = rep(colnames(g500_jja_ord[,-1]), each = nrow(g500_jja_ord[,-1])),  # Repetir estaciones
  Value = Xg500_vector
)
head(Xg500)

#g700 cercano
order<-match(grid_cercanos$station_2, colnames(g700_jja)[-1]) #a estos indices hay que sumar uno por el de fecha que ignaramos
g700_jja_ord <- g700_jja[, c(1, order + 1)]

Xg700_vector <- as.vector(as.matrix(g700_jja_ord[,-1]))
head(Xg700_vector)
Xg700 <- data.frame(
  Date = rep(g700_jja_ord$Date, times = ncol(g700_jja_ord[,-1])),  # Repetir fechas
  Station = rep(colnames(g700_jja_ord[,-1]), each = nrow(g700_jja_ord[,-1])),  # Repetir estaciones
  Value = Xg700_vector
)
head(Xg700)



X_grid <- data.frame(
  Date=Xg700$Date,
  station_x = Xg700$Station,
  g300 = Xg300$Value,
  g500 = Xg500$Value,
  g700 = Xg700$Value
)
head(X_grid)

#add station, t, l
l <- 148:243
l <- rep(l, times = 790 * 64)
t <- year(X_grid$Date) - 1960 + 1

X_grid$l <- l
X_grid$t <- t

X_grid$station<-NA

X_grid <- X_grid[, c('Date', 't', 'l', 'station', 'station_x', 
                     'g300', 'g500', 'g700')]

X_grid$station<-rep(1:790,each = 96*64)

for (i in c('g300', 'g500', 'g700')){
  X_grid[[i]]<-X_grid[[i]]/1000
}
head(X_grid)

# CALCULATION OF ANOMALIES USING THE ORIGINAL ECALATION INFO
# scale
X_grid[, c('g300', 'g500', 'g700')] <-  sweep(X_grid[, c('g300', 'g500', 'g700')], 
                                           2,
                                           scale.info$means, "-")
X_grid[, c('g300', 'g500', 'g700')] <-  sweep(X_grid[, c('g300', 'g500', 'g700')], 
                                              2,
                                              scale.info$sd, "/")

#anomalies
source('harmonics.R')
X_grid <- cbind(X_grid, cs(X_grid$l,1))

# Anomalies. Ref 1981-2010
for (i in 1:790){
  ind <- which(X_grid$station == i)
  ind_jja <- which(X_grid$t[ind] >= 22 & X_grid$t[ind] <= 51 & X_grid$l[ind] >= 152)
  
  for (j in c('g300', 'g500', 'g700')){
    var <- j
    formula <- as.formula(paste(var, "~ s.1 + c.1"))
    mod <- lm(formula, data = X_grid[ind,], subset = ind_jja)
    preds <- predict(mod, newdata = data.frame(
      c.1 = X_grid$c.1[ind],
      s.1 = X_grid$s.1[ind]
    ))
    
    res <- X_grid[ind, var] - preds
    #print(sum(preds[ind_jja] - mod$fitted.values <= 1e-10))
    
    X_grid[ind, var] <- res 
  }
  
}


#añadido de elev y distancia a costa
X_grid <- cbind(X_grid, rep(grid_elev, each = 96 * 64), rep(grid_dist, each = 96 * 64))
colnames(X_grid)[c(ncol(X_grid) - 1, ncol(X_grid))]<-c('elev','dist')

X_grid <- X_grid %>%
  filter(l >= 152)

qsave(X_grid,'X_grid.qs')

