# BCN fire simulation - Data preparation

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcn_fire/data/Re_ Omnsicape Algorihm")

library(terra)
library(sf)
library(tidyverse)
library(exactextractr)


bm_bau <- rast("bmr_ROSbau_raster.tif")
names(bm_bau) <- "biomass_bau"
cn_bau <- rast("bmr_FireConnHbau_raster.tif")
names(cn_bau) <- "connectivity_bau"

# load salvage logging bm and con
bm_sal <- rast("bmr_ROSsal_raster.tif")
names(bm_sal) <- "biomass_sal"
cn_sal <- rast("bmr_FireConnHsal_raster.tif")
names(cn_sal) <- "connectivity_sal"

# plot
plot(bm_bau)
plot(cn_bau)

# Get raster dimensions and extent
cols <- ncol(bm_bau)
rows <- nrow(bm_bau)
ext <- ext(bm_bau)

# Create the sf grid
grd <- st_make_grid(
  st_as_sfc(st_bbox(as.vector(ext))), # Extent as sf bbox
  n = c(cols, rows), # Number of cells in x and y
  what = "polygons", # Output as polygons
  crs = crs(bm_bau) # Maintain the original CRS
) %>% 
  st_sf() # Convert to sf object

# add id column
grd$id <- 1:nrow(grd)

# Add the biomass and connectivity values to the grid
grd$biomass_bau <- exact_extract(bm_bau, grd, fun = "mean")
grd$connectivity_bau <- exact_extract(cn_bau, grd, fun = "mean")
grd$biomass_sal <- exact_extract(bm_sal, grd, fun = "mean")
grd$connectivity_sal <- exact_extract(cn_sal, grd, fun = "mean")

# hist con
hist(grd$connectivity_bau)
hist(grd$connectivity_sal)

max_con <- max(log(1+grd$connectivity_bau), na.rm = T)


# stadardize logged con values to 0-1
grd$connectivity_bau <- log(1+grd$connectivity_bau) / max_con
grd$connectivity_sal <- log(1+grd$connectivity_sal) / max_con

# subset data to non missing biomass
grd <- grd[!is.na(grd$biomass_bau),]

# replace con NA with 0
grd$connectivity_bau[is.na(grd$connectivity_bau)] <- 0
grd$connectivity_sal[is.na(grd$connectivity_sal)] <- 0


# save grid
write_sf(grd, "grd_bau.gpkg")


# work with nb indices instead of ids
# Create a sparse adjacency matrix
adjacent <- st_intersects(grd, sparse = TRUE)

# exclude self-intersections
neighbor_idx <- list()

for (i in 1:length(adjacent)) {
  neighbor_idx[[i]] <- adjacent[[i]] %>% subset(!.==i)
  #neighbor_idx[[i]] <- grd$id[neighbor_index]
}

# save neighbors
save(neighbor_idx, file = "neighbor_idx.RData")
