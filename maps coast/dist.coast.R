rm(list=ls())

# Paquetes
library(sf)
library(sp)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)

# --- Import stations ---
stations <- readRDS('stations.rds')

# stations to sf (WGS84 -> 2062)
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

# --- peninsular limits and grid ---
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

# Spain
spain <- ne_countries(scale = "large", country = "Spain", returnclass = "sf")
spain <- st_transform(spain, 2062)

# Extraction of peninsula
spain_coords <- Polygons(
  list(Polygon(st_coordinates(spain)[st_coordinates(spain)[,"L2"] == 3, 1:2])),
  ID = "spain")
spain_coords <- SpatialPolygons(list(spain_coords))
spain_coords <- as(spain_coords, "sf")
st_crs(spain_coords) <- st_crs(spain)

# Grid (centers 10x10 km) and intersection with peninsula
grid <- st_make_grid(spain, cellsize = 10000, what = "centers")
grid <- st_intersection(grid, spain_coords)
grid <- st_sf(geometry = grid)  # sf

# ---------------------------
# MEDITERRANEAN COAST
# ---------------------------

# 1) Download gloab coastline
coastline_world <- ne_coastline(scale = "large", returnclass = "sf")

# 2) Define a mask that has more or less the mediterranean in CRS 4326
# Adjust values for different results
med_bbox <- st_bbox(c(xmin = -5.5, ymin = 35.0, xmax = 5.0, ymax = 43), crs = st_crs(4326))
med_poly_geo <- st_as_sfc(med_bbox)

# 3) Transfrom to working CRS
med_poly_2062 <- st_transform(med_poly_geo, 2062)
coastline_2062 <- st_transform(coastline_world, 2062)

# 4) Trim coastline to mask defined
med_coastline_2062 <- st_intersection(coastline_2062, med_poly_2062)

# MULTILINESTRING/GEOMETRYCOLLECTION --> LINESTRING
med_coastline_2062 <- st_cast(med_coastline_2062, "LINESTRING", warn = FALSE)

# cHECK EXISTENCE
if(nrow(med_coastline_2062) == 0) {
  stop("La intersección devolvió 0 líneas de costa mediterránea. Ajusta el med_bbox.")
}

# ---------------------------
# DISTANCES (km)
# ---------------------------

# Grid distance
dist_grid_to_med <- st_distance(grid, med_coastline_2062)
# Minimum distance
min_dist_grid_m <- apply(as.matrix(dist_grid_to_med), 1, min)
grid$dist_med_km <- round(min_dist_grid_m / 1000, 3)

# stations dist
dist_st_to_med <- st_distance(stations$geometry, med_coastline_2062)
min_dist_st_m <- apply(as.matrix(dist_st_to_med), 1, min)
stations$dist_med <- round(min_dist_st_m / 1000, 3)

# ---------------------------
# PLOT 
# ---------------------------
library(RColorBrewer)
# Background
world_map <- ne_countries(scale = "large", returnclass = 'sf')
european_union <- c("Algeria", "Andorra", "France", "Gibraltar", "Morocco", 
                    "Portugal", "Spain")
european_union_map <- world_map %>% filter(name %in% european_union)
background <- st_transform(european_union_map, 2062)

ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  ggtitle("Distance to Mediterranean coast (km)") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold")) +
  geom_tile(data = grid, 
            aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], 
                fill = dist_med_km)) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(9, "Reds"),
    limits = c(0, max(grid$dist_med_km)),
    name = "Dist (km)"
  ) +
  #coastline
  geom_sf(data = med_coastline_2062, inherit.aes = FALSE, 
          size = 1, color = 'blue') +
  #stations
  geom_sf(data = stations, inherit.aes = FALSE, pch = 19, size = 1.5) +
  #limits
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

saveRDS(stations$dist_med, 'dist.med.rds')

### SAME PROCESS BUT FOR CANTABRIAN COAST

# ---------------------------
# cANTANBRIAN COAST
# ---------------------------
# 2) Define a mask that has more or less the mediterranean in CRS 4326
# Adjust values for different results
cant_bbox <- st_bbox(c(xmin = -8.6, ymin = 43, xmax = 0, ymax = 44), crs = st_crs(4326))
cant_poly_geo <- st_as_sfc(cant_bbox)

# 3) Transfrom to working CRS
cant_poly_2062 <- st_transform(cant_poly_geo, 2062)

# 4) Trim coastline to mask defined
cant_coastline_2062 <- st_intersection(coastline_2062, cant_poly_2062)

# MULTILINESTRING/GEOMETRYCOLLECTION --> LINESTRING
cant_coastline_2062 <- st_cast(cant_coastline_2062, "LINESTRING", warn = FALSE)

# cHECK EXISTENCE
if(nrow(med_coastline_2062) == 0) {
  stop("La intersección devolvió 0 líneas de costa mediterránea. Ajusta el med_bbox.")
}

# ---------------------------
# DISTANCES (km)
# ---------------------------

# Grid distance
dist_grid_to_cant <- st_distance(grid, cant_coastline_2062)
# Minimum distance
min_dist_grid_m <- apply(as.matrix(dist_grid_to_cant), 1, min)
grid$dist_cant_km <- round(min_dist_grid_m / 1000, 3)

# stations dist
dist_st_to_cant <- st_distance(stations$geometry, cant_coastline_2062)
min_dist_st_m <- apply(as.matrix(dist_st_to_cant), 1, min)
stations$dist_cant <- round(min_dist_st_m / 1000, 3)

# ---------------------------
# PLOT 
# ---------------------------
ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  ggtitle("Distance to Mediterranean coast (km)") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold")) +
  geom_tile(data = grid, 
            aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], 
                fill = dist_cant_km)) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(9, "Reds"),
    limits = c(0, max(grid$dist_cant_km)),
    name = "Dist (km)"
  ) +
  #coastline
  geom_sf(data = cant_coastline_2062, inherit.aes = FALSE, 
          size = 1, color = 'blue') +
  #stations
  geom_sf(data = stations, inherit.aes = FALSE, pch = 19, size = 1.5) +
  #limits
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

saveRDS(stations$dist_cant, 'dist.cant.rds')
