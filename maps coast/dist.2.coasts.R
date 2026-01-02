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

# DEFINITION 2
med_poly <- st_sfc(
  st_polygon(list(matrix(
    c(
      -5.585162405175785, 35.93925648881908,
      0.7694599948747527, 37.05530764678005,
      1.1649678073747527, 40.04456280698777,
      4.548756869874753, 41.91875057870956,
     3.2466502058274216, 42.46183953068916,
      1.906168827119541, 42.59447052504181,
     -0.4760705045471014, 41.914922391934574,
      -4.519039254547101, 39.245413417434705,
     -5.585162405175785, 35.93925648881908
    ),
    ncol = 2,
    byrow = TRUE
  ))),
  crs = 4326
)

# Full world coastline
coastline_world <- ne_coastline(scale = "large", returnclass = "sf")

# Transform to working CRS
coastline_2062 <- st_transform(coastline_world, 2062)
med_coastline_2062 <- st_transform(med_poly, 2062)

med_coastline_2062 <- st_intersection(coastline_2062, med_coastline_2062)

# Cast geometry to LINESTRING
med_coastline_2062 <- st_cast(med_coastline_2062, "LINESTRING", warn = FALSE)

# Check
if (nrow(med_coastline_2062) == 0) {
  stop("The intersection returned 0 coastline segments.")
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

g1 <- ggplot(data = background) +
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

ggsave('graphs/maps/dist.med.coast.pdf',
       g1, height = 8, width = 8)

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

# DEFINITION 2 (ANTLANTIC COAST)
cant_poly <- st_sfc(
  st_polygon(list(matrix(
    c(
      -5.9587733947874355, 39.0034389187482,
        -1.1281003938232703, 43.24800793912893,
        -1.7708005891357703, 43.39387635311883,
        -2.5437567373064507, 44.02353016467899,
        -9.79473329980645, 44.46426631165177,
        -10.159808655475588, 36.98942985231237,
      -5.585162405175785, 35.93925648881908,
      -5.9587733947874355, 39.0034389187482
    ),
    ncol = 2,
    byrow = TRUE
  ))),
  crs = 4326
)

# Full world coastline
coastline_world <- ne_coastline(scale = "large", returnclass = "sf")

# Transform to working CRS
coastline_2062 <- st_transform(coastline_world, 2062)
cant_coastline_2062 <- st_transform(cant_poly, 2062)

cant_coastline_2062 <- st_intersection(coastline_2062, cant_coastline_2062)

# Cast geometry to LINESTRING
cant_coastline_2062 <- st_cast(cant_coastline_2062, "LINESTRING", warn = FALSE)

# Check
if (nrow(cant_coastline_2062) == 0) {
  stop("The intersection returned 0 coastline segments.")
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
g2 <- ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  ggtitle("Distance to Atlantic coast (km)") +
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

ggsave('graphs/maps/dist.cant.coast.pdf',
       g2, height = 8, width = 8)

saveRDS(stations$dist_cant, 'dist.cant.rds')


# -----------------------------------
# NEAREST COAST POINTS PER STATION to both coasts
# -----------------------------------

med_coast_union <- st_union(med_coastline_2062)
cant_coast_union <- st_union(cant_coastline_2062)

med_nearest_lines <- st_nearest_points(
  stations,
  med_coast_union
)
cant_nearest_lines <- st_nearest_points(
  stations,
  cant_coast_union
)

med_nearest_lines_sf <- st_sf(
  id = seq_along(med_nearest_lines),
  geometry = med_nearest_lines
)
cant_nearest_lines_sf <- st_sf(
  id = seq_along(cant_nearest_lines),
  geometry = cant_nearest_lines
)

med_nearest_points_coast <- st_cast(
  med_nearest_lines_sf,
  "POINT"
) %>%
  group_by(id) %>%
  slice(2) %>%
  ungroup()
cant_nearest_points_coast <- st_cast(
  cant_nearest_lines_sf,
  "POINT"
) %>%
  group_by(id) %>%
  slice(2) %>%
  ungroup()

g3 <- ggplot(background) +
  geom_sf(fill = "antiquewhite") +
  #mediterranean coast
  geom_sf(
    data = med_nearest_lines_sf,
    color = "#bdd7e7",
    size = 0.5
  ) +
  geom_sf(
    data = med_coastline_2062,
    color = "#6baed6",
    size = 1
  ) +
  geom_sf(
    data = med_nearest_points_coast,
    color = "#08519c",
    size = 2,
    shape = 15
  ) +
  #cantabrian coast
  geom_sf(
    data = cant_nearest_lines_sf,
    color = "#bae4b3",
    size = 0.5
  ) +
  geom_sf(
    data = cant_coastline_2062,
    color = "#74c476",
    size = 1
  ) +
  geom_sf(
    data = cant_nearest_points_coast,
    color = "#006d2c",
    size = 2,
    shape = 17
  ) +
  geom_sf(
    data = stations,
    color = "black",
    size = 1.8
  ) +
  
  coord_sf(
    xlim = st_coordinates(limits)[, 1],
    ylim = st_coordinates(limits)[, 2]
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggtitle("Nearest coastal point for each station")

ggsave('graphs/maps/points.med.cant.coast.pdf',
       g3, height = 8, width = 8)

