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

dist.matrix <- units::drop_units(round(st_distance(stations)/1000, 3))
saveRDS(dist.matrix, 'maps coast/dist.matrix.rds')
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

# Portugal
portugal <- ne_countries(
  scale = "large",
  country = "Portugal",
  returnclass = "sf"
)
portugal <- st_transform(portugal, 2062)

# Unión España peninsular + Portugal
iberian_peninsula <- st_union(spain_coords, portugal)
iberian_peninsula <- st_make_valid(iberian_peninsula)
iberian_boundary <- st_boundary(iberian_peninsula)

# spain_peninsula <- st_make_valid(spain_coords)
# spain_boundary <- st_boundary(spain_peninsula)

# Grid (centers 10x10 km) and intersection with peninsula
grid <- st_make_grid(spain, cellsize = 10000, what = "centers")
grid <- st_intersection(grid, spain_coords)
grid <- st_sf(geometry = grid)  # sf

# ---------------------------
# FULL COAST
# ---------------------------

# 1) Download gloab coastline
coastline_world <- ne_coastline(scale = "large", returnclass = "sf")

# 2) Define a mask that has more or less the mediterranean in CRS 4326
# Adjust values for different results
bbox <- st_bbox(c(xmin = -10, ymin = 35.0, xmax = 5.0, ymax = 44), crs = st_crs(4326))
poly_geo <- st_as_sfc(bbox)

# 3) Transfrom to working CRS
poly_2062 <- st_transform(poly_geo, 2062)
coastline_2062 <- st_transform(coastline_world, 2062)

# 4) Trim coastline to mask defined
full_coastline_2062 <- st_intersection(coastline_2062, poly_2062)

# MULTILINESTRING/GEOMETRYCOLLECTION --> LINESTRING
full_coastline_2062 <- st_cast(full_coastline_2062, "LINESTRING", warn = FALSE)
# cHECK EXISTENCE
if(nrow(full_coastline_2062) == 0) {
  stop("La intersección devolvió 0 líneas de costa mediterránea. Ajusta el med_bbox.")
}

# DEFINITION 2 (CUSTOM POLYGON) GOOD ONE
# Iberian Peninsula polygon (EPSG:4326)
iberian_poly_4326 <- st_sfc(
  st_polygon(list(matrix(
    c(
      -1.7973703292203114, 43.420658518267985,
      -10.685085614545189, 44.277700140621846,
      -9.806179364545189, 35.24679268225082,
      -5.398430544071187, 35.95919861144038,
      0.32249797600901964, 36.981823151967866,
      1.0695682885090196, 40.17583174219498,
      4.23363078850902, 41.60425160613797,
      3.301653557815749, 42.45836262739984,
      -1.7973703292203114, 43.420658518267985
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
iberian_poly_2062 <- st_transform(iberian_poly_4326, 2062)

iberian_coast_2062 <- st_intersection(coastline_2062, iberian_poly_2062)
iberian_coast_2062 <- iberian_coast_2062[1, ]
# Cast geometry to LINESTRING
iberian_coast_2062 <- st_cast(iberian_coast_2062, "LINESTRING", warn = FALSE)

# Check
if (nrow(iberian_coast_2062) == 0) {
  stop("The intersection returned 0 coastline segments.")
}

# ---------------------------
# DISTANCES (km)
# ---------------------------

full_coastline_2062 <- iberian_coast_2062

# Grid distance
dist_grid_to_coast <- st_distance(grid, full_coastline_2062)
# Minimum distance
min_dist_grid_m <- apply(as.matrix(dist_grid_to_coast), 1, min)
grid$dist <- round(min_dist_grid_m / 1000, 3)

# stations dist
dist_st_to_coast <- st_distance(stations$geometry, full_coastline_2062)
min_dist_st_m <- apply(as.matrix(dist_st_to_coast), 1, min)
stations$dist <- round(min_dist_st_m / 1000, 3)

saveRDS(stations$dist, 'maps coast/dist.vec.rds')

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

g <- ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  ggtitle("Distance to coast (km)") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold")) +
  geom_tile(data = grid, 
            aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], 
                fill = dist)) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(9, "Reds"),
    limits = c(0, max(grid$dist)),
    name = "Dist (km)"
  ) +
  #coastline
  geom_sf(data = full_coastline_2062, inherit.aes = FALSE, 
          size = 1, color = 'blue') +
  #stations
  geom_sf(data = stations, inherit.aes = FALSE, pch = 19, size = 1.5) +
  #limits
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

ggsave('graphs/maps/dist.coast.pdf',
       g, height = 8, width = 8)



# -----------------------------------
# NEAREST COAST POINTS PER STATION
# -----------------------------------

coast_union <- st_union(full_coastline_2062)

nearest_lines <- st_nearest_points(
  stations,
  coast_union
)

nearest_lines_sf <- st_sf(
  id = seq_along(nearest_lines),
  geometry = nearest_lines
)

nearest_points_coast <- st_cast(
  nearest_lines_sf,
  "POINT"
) %>%
  group_by(id) %>%
  slice(2) %>%
  ungroup()


dist.coast.points <- units::drop_units(st_distance(nearest_points_coast))
saveRDS(dist.coast.points, 'maps coast/dist.coast.points.rds')

g2 <- ggplot(background) +
  geom_sf(fill = "antiquewhite") +
  geom_sf(
    data = nearest_lines_sf,
    color = "gray40",
    size = 0.5
  ) +
  geom_sf(
    data = nearest_points_coast,
    color = "red",
    size = 2
  ) +
  geom_sf(
    data = stations,
    color = "black",
    size = 1.8
  ) +
  geom_sf(
    data = full_coastline_2062,
    color = "blue",
    size = 0.5
  ) +
  coord_sf(
    xlim = st_coordinates(limits)[, 1],
    ylim = st_coordinates(limits)[, 2]
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue")
  ) +
  ggtitle("Nearest coastal point for each station")

ggsave('graphs/maps/points.coast.pdf',
       g2, height = 8, width = 8)

# -----------------------------------
# PARAMETRIZATION r(s)
# -----------------------------------

#coast_line <- st_line_merge(coast_union)
coast_coords <- st_coordinates(coast_union)[, 1:2]
dr <- sqrt(
  diff(coast_coords[,1])^2 +
    diff(coast_coords[,2])^2
)

r_coast <- c(0, cumsum(dr))
coast_length_km <- round(r_coast[length(r_coast)]/1000, 3)
coast_proj <- st_coordinates(nearest_points_coast)

nearest_vertex <- apply(coast_proj, 1, function(p) {
  which.min(
    (coast_coords[,1] - p[1])^2 +
      (coast_coords[,2] - p[2])^2
  )
})


stations$r <- r_coast[nearest_vertex] / 1000  # km

#matix of distances
dist.coast.points2 <- abs(outer(stations$r, stations$r, '-'))
saveRDS(dist.coast.points2, 'maps coast/dist.coast.points2.rds')

ggplot(background) +
  geom_sf(fill = "antiquewhite") +
  geom_sf(data = full_coastline_2062, color = "blue") +
  
  # estaciones
  geom_sf(data = stations, aes(color = r), size = 2) +
  
  # puntos de la costa seleccionados
  geom_sf(
    data = nearest_points_coast,
    color = "black",
    shape = 17,
    size = 5
  ) + 
  geom_point(aes(x = X, y = Y), data = coast_coords[nearest_vertex, ],
             color = 'red',
             size = 5,
             shape = 17) +
  
  scale_color_viridis_c(name = "r(s) [km]") +
  coord_sf(
    xlim = st_coordinates(limits)[,1],
    ylim = st_coordinates(limits)[,2]
  )



# PARAM in 200 NODES
plot(coast_union)
points(st_coordinates(coast_union), pch = 19, cex = 0.25)
coast_line <- coast_union
s <- seq(0, 1, length.out = 200)
coast_nodes <- st_line_sample(coast_line, sample = s)
points(st_coordinates(coast_nodes), col = 'red', pch = 19)

coastlen <- st_length(coast_line)
r <- coastlen * s

# matrices for MCMC
coast_points <- st_cast(coast_nodes, "POINT")
phimat <- round(units::drop_units(st_distance(stations$geometry, coast_points)/1000), 3)
saveRDS(phimat, 'maps coast/phimat.rds')

#ri - rj
dr <- units::drop_units(abs(outer(r/1000, r/1000, '-')))
saveRDS(dr, 'maps coast/dr.rds')

# -----------------------------------
# GRID AND DISTNACES TO 200 nodes
# -----------------------------------
grid <- st_make_grid(spain, cellsize = 25000, what = "centers")
grid <- st_intersection(grid, spain_coords)
grid <- st_sf(geometry = grid)  # sf

phimat.grid <- round(units::drop_units(st_distance(grid$geometry, coast_points)/1000), 3)
saveRDS(phimat.grid, 'maps coast/phimat.grid.rds')

#quick check
grid$dist_coast_10 <- phimat.grid[, 10]
ggplot(background) +
  geom_sf(fill = "antiquewhite") + 
  
  # estaciones
  geom_sf(data = stations, size = 2) +
  geom_sf(data = coast_points[10], col = 'blue') + 
  geom_sf(data = grid, shape = 15, aes(color = dist_coast_10)) +
  scale_color_viridis_c() +
  coord_sf(
    xlim = st_coordinates(limits)[,1],
    ylim = st_coordinates(limits)[,2]
  )


# checks for future covariances matrices
stations.basura <- rbind(grid[, 0], stations[, 0])
phimat.basura <- round(units::drop_units(st_distance(stations.basura$geometry, coast_points)/1000), 3)

K.coast <- exp(-0.003*dr)

final <- phimat.basura %*% K.coast %*% t(phimat.basura)
R21 <- final[791:830, 1:790]

final2 <- phimat %*% K.coast %*% t(phimat.grid)

all.equal(R21, final2)
