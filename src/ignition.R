library(terra)
library(sf)

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

grd <- read_sf("./data/rdat/grd.gpkg")

lc <- rast("./data/landcover/BMRlandcover.tif")

roads <- lc == 31

# transform road raster mask to road polygons
roads_vec <- terra:: as.polygons(roads)

roads_sf <- sf::st_as_sf(roads_vec) %>% 
  filter(BMRlandcover == 1) %>%
  st_transform(roads_sf, crs = st_crs(grd))

road_dist <- st_distance(grd, roads_sf)

grd$road_distance <- as.numeric(road_dist)

# Normalize the distance to create a probability surface
max_distance <- max(grd$road_distance, na.rm = TRUE)
grd$ignition_probability_surface <- 1 - (grd$road_distance / max_distance)

# Ensure all cells have a probability larger than zero
min_prob <- 0.01
grd$ignition_probability_surface <- grd$ignition_probability_surface + min_prob
grd$ignition_probability_surface <- grd$ignition_probability_surface / sum(grd$ignition_probability_surface, na.rm = TRUE)

# save the grid
write_sf(grd, "./data/rdat/grd.gpkg")

# plot using leaflet
library(leaflet)
mapview::mapview(grd, zcol = "ignition_probability_surface")


# draw 100000 cells and verify that the probability surface is respected
# selected_cells <- data.frame(selected_cells = sample(x = unique(grd$id), 
#                                                      size = 10000000, 
#                                                      replace = TRUE,
#                                                      prob = grd$ignition_probability_surface))
# 
# # count cell selections and compare to the probability surface
# cell_count <- selected_cells %>%
#   group_by(selected_cells) %>%
#   summarise(Freq = n()) %>%
#   left_join(data.frame(id = grd$id, prob = grd$ignition_probability_surface), by = c("selected_cells" = "id"))
# 
# plot(cell_count$prob, cell_count$Freq)
