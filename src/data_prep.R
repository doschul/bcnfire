# BCN fire simulation - Data preparation

##### Setup #####

rm(list=ls())

# working directory
setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

#setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire/data/Re_ Omnsicape Algorihm")

# libraries
library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(spdep) # for poly2nb

rcl_lc <- FALSE # whether or not to reclassify raw landcover raster


# data paths

path_topo <- "./data/topo/topo_bcn.tif"
path_wui <- "./data/wui/WUI_Catalonia.tif"

path_lc_raw <- "./data/landcover/BMRlandcover.tif"
path_lc_rcl <- "./data/landcover/BMRlandcover_rcl.tif"

path_roads_sf <- "./data/landcover/roads_sf.gpkg"


path_ros_bau_avg <- "./data/Re_ Omnsicape Algorihm/bmr_ROSbau_raster.tif"
path_con_bau_avg <- "./data/Re_ Omnsicape Algorihm/bmr_FireConnHbau_raster.tif"
path_ros_sal_avg <- "./data/Re_ Omnsicape Algorihm/bmr_ROSsal_raster.tif"
path_con_sal_avg <- "./data/Re_ Omnsicape Algorihm/bmr_FireConnHsal_raster.tif"


path_annual_mask_rast <- "./data/Re_ Omnsicape Algorihm/LOG_rasters"
path_annual_ros_rast <- "./data/Re_ Omnsicape Algorihm/ROS_rasters"
path_annual_con_rast <- "./data/Re_ Omnsicape Algorihm/CON_rasters"


##### Load and prepare data #####

# load average scenarios, BAU and SAL, Con and Ros
bm_bau <- rast(path_ros_bau_avg)
cn_bau <- rast(path_con_bau_avg)
bm_sal <- rast(path_ros_sal_avg)
cn_sal <- rast(path_con_sal_avg)

r.bm <- c(bm_bau, bm_sal)
r.cn <- c(cn_bau, cn_sal)

names(r.bm) <- c("biomass_bau_avg", "biomass_sal_avg")
names(r.cn) <- c("connectivity_bau_avg", "connectivity_sal_avg")


# load and c rasters 
annual_mask_rast <- list.files(path_annual_mask_rast, pattern = ".tif$", full.names = TRUE) |>
  lapply(terra::rast) |>
  do.call(c, args = _)

annual_ros_rast <- list.files(path_annual_ros_rast, pattern = ".tif$",full.names = TRUE) %>%
  lapply(terra::rast) %>%
  do.call(c, .) 

annual_con_rast <- list.files(path_annual_con_rast, pattern = ".tif$",full.names = TRUE) %>%
  lapply(terra::rast) %>%
  do.call(c, .) 

# Assign names to the SpatRaster object using base R's names() function
names(annual_mask_rast) <- gsub(".tif", "", list.files(path_annual_mask_rast, pattern = ".tif$"))
names(annual_ros_rast) <- gsub(".tif", "", list.files(path_annual_ros_rast, pattern = ".tif$"))
names(annual_con_rast) <- gsub(".tif", "", list.files(path_annual_con_rast, pattern = ".tif$"))

# load topography raster
topo <- rast(path_topo)
topo$aspect <- NULL


# relclassify landcover into aggregate categories, if needed.
if(rcl_lc){
  # extract landcover data
  lc <- rast(path_lc_raw)
  
  # 1. Agricultural Areas (Arees agrícoles)
  # Conreus herbacis (1): Herbaceous crops
  # Horta, vivers i conreus forçats (2): Orchards, nurseries, and forced crops (greenhouse cultivation)
  # Vinyes (3): Vineyards
  # Oliverars (4): Olive groves
  # Altres conreus llenyosos (5): Other woody crops
  # Conreus en transformació (6): Crops in transformation (areas under conversion)
  
  # 2. Forest and Natural Areas (Arees forestals i naturals)
  # Boscos densos d'aciculifolis (7): Dense coniferous forests
  # Boscos densos de caducifolis, planifolis (8): Dense deciduous broadleaf forests
  # Boscos densos d'esclerofil-les i laurifolis (9): Dense sclerophyllous and laurel forests
  # Matollar (10): Shrubland
  # Boscos clars d'aciculifolis (11): Open coniferous forests
  # Boscos clars de caducifolis, planifolis (12): Open deciduous broadleaf forests
  # Boscos clars d'esclerofil-les i laurifolis (13): Open sclerophyllous and laurel forests
  # Prats i herbassars (14): Meadows and grasslands
  # Bosc de ribera (15): Riparian forest
  # Sòl nu forestal (16): Bare forest soil
  # Zones cremades (17): Burned areas
  # Roquissars i congestes (18): Rocky areas and scree
  # Platges (19): Beaches
  # Zones humides (20): Wetlands
  
  # 3. Urbanized Areas (Arees urbanitzades)
  # Casc urbà (21): Urban core
  # Eixample (22): Urban expansion (often referring to planned extensions)
  # Zones Urbanes laxes (23): Sparse urban areas
  # Edificacions aïllades en l'espai rural (24): Isolated buildings in rural areas
  # Arees residencials aïllades (25): Isolated residential areas
  # Zones verdes (26): Green areas (parks)
  # Zones industrials, comercials i/o de serveis (27): Industrial, commercial, and/or service areas
  # Zones esportives i de lleure (28): Sports and leisure areas
  # Zones d'extracció minera i/o abocadors (29): Mining and/or landfill areas
  # Zones en transformació (30): Areas in transformation (under development)
  # Xarxa viària (31): Road network
  # Sòl nu urbà (32): Bare urban soil
  # Zones aeroportuàries (33): Airport areas
  # Xarxa ferroviària (34): Railway network
  # Zones portuàries (35): Port areas
  
  # 4. Water Bodies (Masses d'aigua)
  # Embassaments (36): Reservoirs
  # Llacs i llacunes (37): Lakes and lagoons
  # Cursos d'aigua (38): Watercourses (rivers, streams)
  # Basses (39): Ponds
  # Canals artificials (40): Artificial canals
  # Mar (41): Sea
  # Other
  #                  
  # Sense dades (0): No data
  
  # Land cover classification:
  # agriculture <- 1 to 6
  # wildland <- 7 to 17 and 20
  # rock <- 18 and 19
  # artificial <- 21 to 35 (roads <- 31)
  # water <- 36 to 41
  
  rcl_m <- matrix(c(1,   6, 1, # agriculture
                    7,  17, 2, # wildland
                    18, 19, 3, # rock
                    20, 20, 2, # wildland
                    21, 35, 4, # artificial
                    36, 41, 5),# water
                  ncol = 3, byrow = T)
  
  # reclassify landcover
  lc.rcl <- terra::classify(lc, rcl_m, right = NA)
  
  # export reclasified landcover
  writeRaster(lc.rcl, path_lc_rcl, overwrite = T)
  
  # extract roads and save as shape
  roads <- lc == 31
  
  # transform road raster mask to road polygons
  roads_vec <- terra:: as.polygons(roads)
  
  roads_sf <- sf::st_as_sf(roads_vec) %>% 
    filter(BMRlandcover == 1)
  
  # save as gpkg 
  write_sf(roads_sf, dsn = path_roads_sf)
  
} else {
  lc.rcl <- rast(path_lc_rcl)
  roads_sf <- read_sf(path_roads_sf)
}

##### Make grid #####

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

grd_ctrd <- st_centroid(grd)


##### Extract values #####

#### Road distance & Ignition probability ####
roads_sf <- st_transform(roads_sf, crs = st_crs(grd))
road_dist <- st_distance(grd, roads_sf)
grd$road_distance <- as.numeric(road_dist)

# Normalize the distance to create a probability surface
max_distance <- max(grd$road_distance, na.rm = TRUE)
grd$ignition_probability_surface <- 1 - (grd$road_distance / max_distance)

# Ensure all cells have a probability larger than zero
min_prob <- 0.01
grd$ignition_probability_surface <- grd$ignition_probability_surface + min_prob
grd$ignition_probability_surface <- grd$ignition_probability_surface / sum(grd$ignition_probability_surface, na.rm = TRUE)


#### Land cover and topography ####
lc_shares <- exact_extract(lc.rcl, grd, fun = "frac")
names(lc_shares) <- c("lc_agriculture", "lc_wildland", "lc_rock", "lc_artificial", "lc_water")

e.topo <- exact_extract(topo, grd, fun = c("mean", "max"))

grd <- grd %>%
  cbind(lc_shares) %>%
  cbind(e.topo)

#### Average scenarios ####

e.avg.bm <- terra::extract(r.bm, grd_ctrd, ID = F)
e.avg.cn <- terra::extract(r.cn, grd_ctrd, ID = F)

# annual values of masks, biomass and connectivity
e.annual.msk <- terra::extract(annual_mask_rast, grd_ctrd, ID = F)
e.annual.ros <- terra::extract(annual_ros_rast, grd_ctrd, ID = F)
e.annual.con <- terra::extract(annual_con_rast, grd_ctrd, ID = F)

# merge with grd data
grd <- grd %>% 
  cbind(. , e.avg.bm) %>%
  cbind(. , e.avg.cn) %>%
  cbind(. , e.annual.ros) %>%
  cbind(. , e.annual.con) %>%
  cbind(. , e.annual.msk)

# save full grid (unfiltered)
write_sf(grd, "./data/rdat/grd_full.gpkg")

apply(grd, 2, function(x){sum(is.na(x))})

# Filter data and save
#max_con <- max(grd$connectivity_bau, na.rm = T)

# stadardize logged con values to 0-1
#grd$connectivity_bau <- grd$connectivity_bau / max_con
#grd$connectivity_sal <- grd$connectivity_sal / max_con

# subset data to non missing biomass
grd_f <- grd %>%
  filter(!is.na(biomass_bau_avg),
         !is.na(mean.elev)) %>%
  # replace all NA with 0
  replace(is.na(.), 0)


row.names(grd) <- paste0("cell_", 1:nrow(grd)) # Set row names to "cell_1", "cell_2", etc.
row.names(grd_f) <- paste0("cell_", 1:nrow(grd_f)) # Set row names to "cell_1", "cell_2", etc.

# load wildland urba interface
wui <- rast(path_wui)
wui_e <- exactextractr::exact_extract(wui, grd_f, fun="frac")

names(wui_e) <- c("wui_low", "wui_medium", "wui_nohouse", "wui_verylow", "wui_intermix", "wui_interface")
grd_f <- bind_cols(grd_f, wui_e)

# save grid
write_sf(grd_f, "./data/rdat/grd_filt.gpkg")

##### Create neighborhood list #####

# 1. Read the spatial data frame
grd <- st_read("./data/rdat/grd_filt.gpkg")
row.names(grd) <- paste0("cell_", 1:nrow(grd)) # Set row names to "cell_1", "cell_2", etc.

# 2. Get all queen and all rook neighbors (these already return lists of numeric indices)
nb_queen_all_indices <- poly2nb(grd, queen = TRUE)
nb_rook_indices <- poly2nb(grd, queen = FALSE)

# 3. Directly create the differentiated neighbor list using lapply
# This builds the list where each element is a list(rook = ..., queen = ...)
# and the outer list is naturally indexed by cell position (1 to nrow(grd)).
neighbor_idx_differentiated <- lapply(seq_along(nb_queen_all_indices), function(i) {
  current_rook_indices <- nb_rook_indices[[i]]
  current_all_queen_indices <- nb_queen_all_indices[[i]]
  
  # Determine exclusive queen neighbors by removing rook neighbors from all queen neighbors
  exclusive_queen_indices <- setdiff(current_all_queen_indices, current_rook_indices)
  
  list(
    rook = current_rook_indices,
    queen = exclusive_queen_indices
  )
})

# 4. Assign original row names to the list elements for direct lookup later
# This ensures that neighbor_idx_differentiated[[ignition_idx]] works correctly
# when ignition_idx is a row name string.
names(neighbor_idx_differentiated) <- row.names(grd)

# The variable 'neighbor_idx_differentiated' now contains the pre-processed list
# in the exact format expected by get_burners_optimized.

# save the neighbor list for later use
save(neighbor_idx_differentiated, file = "./data/rdat/neighbor_idx_differentiated.RData")

