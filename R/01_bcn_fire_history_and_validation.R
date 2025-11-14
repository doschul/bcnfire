# Historic wildfire ignition pattern

# Setup
rm(list=ls())

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")


# ---------- Setup ----------
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(tidyverse)
  library(stringr)
  library(osmextract)
  library(fs)
  library(units)
  library(exactextractr)
  library(here)
  library(ggpubr)
  library(scales) 
  library(MASS) # Added for the fitdistr function
  library(fitdistrplus)
  library(randomForest)
  library(blockCV)
  library(pROC)
})



if(F){
  # ---------------------------------------
  # 1. Download Geofabrik extract (Catalonia)
  # ---------------------------------------
  # Default download folder (~/.osmextract by default)
  
  its_match = oe_match("catalonia", quiet = TRUE)
  
  oe_download(
    file_url = its_match$url,
    file_size = its_match$file_size,
    file_basename = "catalonia",
    provider = "geofabrik", 
    download_directory = "./data/osm"
  )
}


# Input files -------------------------------------------------
pbf_path      <- here("data", "osm", "geofabrik_catalonia.pbf")
topo_path     <- here("data", "topo", "topo_bcn.tif")
landcover_path <- here("data", "landcover", "BMRlandcover_rcl.tif")
grid_path     <- here("data", "rdat", "catalunya_grd.rds")   # assumes saved grid

# Output file -------------------------------------------------
covariate_raster_out <- here("data", "rdat", "covariates_cat.tif")


#### Load data, create Catalunya grid ####

# study region: Catalunya
# load grd shape
grd <- st_read("./data/rdat/grd_filt.gpkg")

comarcas <- st_read("./data/divisions-administratives-v2r1-20250101/divisions-administratives-v2r1-comarques-1000000-20250101.json") |> 
  st_transform(crs = st_crs(grd)) |> 
  st_sf() |>       
  st_make_valid()

aoi <- st_union(comarcas)


# make 1km grid across study region
catalunya_grd <- st_make_grid(aoi, 
                              cellsize = units::as_units(500, "meter")) |> 
  st_make_valid() |> 
  st_intersection(aoi)
catalunya_grd <- st_sf(id = 1:length(catalunya_grd), geometry = catalunya_grd)
saveRDS(catalunya_grd, grid_path)

wildfires <- list.files("./data/wildfires/wildfires_shp", pattern = ".shp$", 
                        recursive = TRUE, full.names = TRUE, include.dirs = TRUE) |>
  lapply(st_read) |>
  bind_rows()

wildfires <- wildfires |>
  filter(!is.na(CODI_FINAL)) |>
  mutate(
    year = case_when(nchar(DATA_INCEN) == 7 ~ as.integer(paste0("20", substr(DATA_INCEN, 6, 7))),
                     nchar(DATA_INCEN) == 8 ~ as.integer(paste0("20", substr(DATA_INCEN, 7, 8))),
                     nchar(DATA_INCEN) >= 10 ~ as.integer(substr(DATA_INCEN, 7, 10)),
                     TRUE ~ NA),
    month = case_when(nchar(DATA_INCEN) == 7 ~ as.integer(substr(DATA_INCEN, 3, 4)),
                      nchar(DATA_INCEN) >= 8 ~ as.integer(substr(DATA_INCEN, 4, 5)),
                      TRUE ~ NA),
    # day is anything before first slash
    day = as.integer(substring(DATA_INCEN, 1, regexpr("/", DATA_INCEN) - 1)),
    date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d"),
    area_ha = as.numeric(st_area(geometry)) / 10000,  # convert from m² to ha
    MUNICIPI = str_to_title(MUNICIPI)  # capitalize first letter of each word
  ) |>
  dplyr::select(CODI_FINAL, MUNICIPI, date, year, month, area_ha, geometry)


# transform crs
wildfires <- st_transform(wildfires, crs = st_crs(grd)) %>%
  st_make_valid() |>
  # valid geom only
  filter(st_is_valid(geometry)) |>
  mutate(log_area_ha = log(area_ha)) #|> filter(year >= 2015 & year <= 2025)


# export for use in Python
wf_observed <- wildfires %>%
  st_set_geometry(NULL)
write_csv(wf_observed, "./data/wildfires/wf_observed.csv")
wf_observed <- read_csv("./data/wildfires/wf_observed.csv")

wf_aoi <- wf_observed
# optional: Study region only
# wf_aoi <- wf_aoi[st_intersects(wf_aoi, ch, sparse = FALSE), ]

##### Wildfire descriptive analysis #####

# Create Descriptive plot 
wf_summary <- wf_aoi %>%
  #st_set_geometry(NULL) %>%
  # Group and summarise the observed data
  group_by(year) %>%
  summarise(
    n = n(),
    area_tot = sum(area_ha),
    .groups = 'drop' # Recommended to drop grouping immediately after summarise
  ) %>%
  # Complete the time series
  complete(
    year = seq(min(year), max(year)), # Create the sequence of all years in the observed range
    fill = list(n = 0, area_tot = 0)  # Fill missing rows with zeros
  )


# --- 2. CALCULATE SCALING FACTOR ---
# The line data (area_tot) must be scaled down to fit the bar data (n) scale.
# The scaling factor is the ratio of the maximum values of the two metrics.
scaling_factor <- max(wf_summary$area_tot) / max(wf_summary$n)


# --- 3. CREATE DUAL-AXIS PLOT (p1) ---
p1 <- wf_summary %>%
  ggplot(aes(x = year)) +
  
  # 1. Bar Geometry (Primary Axis: Frequency 'n')
  geom_col(
    aes(y = n),
    fill = "#1f78b4", # Blue color for bars
    alpha = 0.7 
  ) +
  
  # 2. Line Geometry (Secondary Axis: Area 'area_tot')
  # We divide 'area_tot' by the scaling_factor to plot it on the primary 'n' scale.
  geom_line(
    aes(y = area_tot / scaling_factor, group = 1),
    color = "#e31a1c", # Red color for line
    linewidth = 1.2
  ) +
  # Add points to the line for visual clarity
  geom_point(
    aes(y = area_tot / scaling_factor),
    color = "#e31a1c",
    size = 3
  ) +
  
  # 3. Define the Y-Axis Scales
  scale_y_continuous(
    # Primary Y-Axis (Left) - for the bars ('n')
    name = "Frequency (Count)",
    
    # Secondary Y-Axis (Right) - for the line ('area_tot')
    sec.axis = sec_axis(
      trans = ~ . * scaling_factor,
      name = "Total Area (ha)"
    )
  ) +
  
  # 4. Enhance Aesthetics
  scale_x_continuous(
    name = "Year",
    breaks = wf_summary$year, 
    guide = guide_axis(n.dodge = 2) 
  ) +
  labs(
    title = "A: Annual Event Frequency and Total Area",
    subtitle = paste0("Scaling Factor: ", round(scaling_factor, 2))
  ) +
  theme_minimal() +
  theme(
    # Custom colors to match the geoms
    axis.title.y.left = element_text(color = "#1f78b4", margin = margin(r = 10)),
    axis.text.y.left = element_text(color = "#1f78b4"),
    axis.title.y.right = element_text(color = "#e31a1c", margin = margin(l = 10)),
    axis.text.y.right = element_text(color = "#e31a1c"),
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    plot.margin = margin(5, 5, 5, 5) 
  )


# --- 4. CREATE FREQUENCY HISTOGRAM (p2) - WITH FITTED COUNT DISTRIBUTIONS ---

# Fit Poisson and Negative Binomial distributions to the annual event counts (n)
counts <- wf_summary$n
fit_pois <- MASS::fitdistr(counts, "poisson")
fit_nbinom <- MASS::fitdistr(counts, "negative binomial")

# Prepare data frame for plotting the fitted lines
max_count <- max(counts)
x_val <- min(counts):max_count
df_fit_counts <- data.frame(n = x_val)

# Calculate theoretical counts (PMF * Total Observations, which is 10 years)
df_fit_counts <- df_fit_counts %>%
  mutate(
    pois_expected = dpois(n, lambda = fit_pois$estimate["lambda"]) * length(counts),
    nbinom_expected = dnbinom(n, size = fit_nbinom$estimate["size"], mu = fit_nbinom$estimate["mu"]) * length(counts)
  )

# Extract parameters for annotation
pois_lambda <- format(fit_pois$estimate["lambda"], digits = 3)
nbinom_mu <- format(fit_nbinom$estimate["mu"], digits = 3)
nbinom_size <- format(fit_nbinom$estimate["size"], digits = 3)
p2_subtitle <- paste0(
  "Poisson [dashed] (\u03bb): ", pois_lambda, 
  " | NegBin [dotted] (\u03bc, size): ", nbinom_mu, ", ", nbinom_size
)

p2 <- wf_summary %>%
  ggplot(aes(x = n)) +
  # Note: binwidth is used here to group the 10 annual counts for visualization
  geom_histogram(bins = max_count, fill = "#a6cee3", color = "white") + 
  
  # Add fitted distributions (plotted as expected counts/frequencies)
  geom_line(data = df_fit_counts, aes(x = n, y = pois_expected), 
            color = "#1b9e77", linewidth = 1, linetype = "dashed") + # Green for Poisson
  geom_line(data = df_fit_counts, aes(x = n, y = nbinom_expected), 
            color = "#d95f02", linewidth = 1, linetype = "dotted") + # Orange for NegBin
  
  labs(
    title = "B: Annual Event Frequency Distribution",
    subtitle = p2_subtitle,
    x = "Events per Year (n)",
    y = "Count of Years (Frequency)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5) 
  )


# --- 5. CREATE AREA PER EVENT HISTOGRAM (p3) - WITH FITTED CONTINUOUS DISTRIBUTIONS (NORMAL, STUDENT T, and GEV on log-transformed data) ---

# 1. Log-transform the non-zero area data
log_area_data <- wf_aoi$log_area_ha[wf_aoi$area_ha > 0]

# 2a. Fit Normal distribution to the logged data (for Log-Normal fit)
fit_norm <- MASS::fitdistr(log_area_data, "normal")
meanlog_fit <- fit_norm$estimate["mean"]
sdlog_fit <- fit_norm$estimate["sd"]

# 2b. Fit Student-t distribution to the logged data
# We'll use custom MLE fitting through fitdistrplus::fitdist

fit_t <- MASS::fitdistr(log_area_data, "t")

t_mean <- fit_t$estimate["m"]
t_scale <- fit_t$estimate["s"]
t_df <- fit_t$estimate["df"]

# 2c. Fit Generalized Extreme Value (GEV) distribution to the logged data
fit_gev <- evd::fgev(log_area_data, std.err = FALSE)
gev_loc <- fit_gev$estimate["loc"]
gev_scale <- fit_gev$estimate["scale"]
gev_shape <- fit_gev$estimate["shape"]

# 3. Format parameters for subtitle
lognorm_meanlog <- format(meanlog_fit, digits = 2)
lognorm_sdlog <- format(sdlog_fit, digits = 2)
t_mean_f <- format(t_mean, digits = 2)
t_scale_f <- format(t_scale, digits = 2)
t_df_f <- format(t_df, digits = 2)
gev_loc_f <- format(gev_loc, digits = 2)
gev_scale_f <- format(gev_scale, digits = 2)
gev_shape_f <- format(gev_shape, digits = 2)

# Subtitle text
p3_subtitle <- paste0(
  "Norm [dashed] (μ, σ): ", lognorm_meanlog, ", ", lognorm_sdlog,
  " | Student-t [dotted] (m, s, df): ", t_mean_f, ", ", t_scale_f, ", ", t_df_f,
  " | GEV [solid] (μ, σ, ξ): ", gev_loc_f, ", ", gev_scale_f, ", ", gev_shape_f
)

# 4. Create the plot using the logged data on the X-axis
p3 <- wf_aoi %>%
  ggplot(aes(x = log_area_ha)) +
  
  # Histogram of log-transformed data (density scale)
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "#fb9a99",
    color = "white"
  ) +
  
  # 1. Normal Fit (Blue Dashed)
  stat_function(
    fun = dnorm,
    args = list(mean = meanlog_fit, sd = sdlog_fit),
    color = "#1f78b4",
    linewidth = 1,
    linetype = "dashed"
  ) +
  
  # 2. Student-t Fit (Red Dotted)
  stat_function(
    fun = function(x) dt((x - t_mean) / t_scale, df = t_df) / t_scale,
    color = "#e31a1c",
    linewidth = 1,
    linetype = "dotted"
  ) +
  
  # 3. GEV Fit (Purple Solid)
  stat_function(
    fun = evd::dgev,
    args = list(loc = gev_loc, scale = gev_scale, shape = gev_shape),
    color = "#7570b3",
    linewidth = 1,
    linetype = "solid"
  ) +
  
  labs(
    title = "C: Area Per Event Distribution (Log-Scale with Dual Fit)",
    subtitle = p3_subtitle,
    x = "Area per Event (ha)",
    y = "Density"
  ) +
  
  # Display log10 axis labels (convert from log-values)
  scale_x_continuous(
    labels = function(x) {
      parse(text = paste0("10^", round(x / log(10), 0)))
    },
    breaks = log(c(1, 10, 100, 1000, 10000, 100000))
  ) +
  
  coord_cartesian(xlim = range(log_area_data)) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 7, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )


# --- 6. COMBINE PLOTS using ggpubr ---
final_plot <- ggarrange(
  p1, p2, p3,
  ncol = 1, nrow = 3,
  heights = c(0.4, 0.3, 0.3) # 50% for p1, 50% for the stacked histograms
)
final_plot
ggsave(final_plot, filename = "./out/fig/wildfire_distr.png",
       width = 5, height = 6, bg = "white")


##### Predictive ignition grid based on covariates #####

##### Prapare covariate rasters #####


# ============================================================
# Build Catalonia Environmental Covariate Raster Stack
# Input: OSM PBF, existing 200m grid catalunya_grd (sf polygon grid)
# Output: One multi-band raster: data/rdat/covariates_cat.tif
# Bands: road densities, trail densities, campsite distance, topo, landcover fractions
# ============================================================



# ---------- Load grid ----------
catalunya_grd <- readRDS(grid_path)   # sf polygon grid with id column
crs_proj <- st_crs(catalunya_grd)

# ---------- Load OSM data ----------
osm_lines  <- oe_read(pbf_path, layer = "lines")
osm_points <- oe_read(pbf_path, layer = "points")

# Filter OSM features ----------------------------------------
major_highways <- c("primary", "secondary", "tertiary", "unclassified")
trail_types    <- c("footway", "path", "pedestrian", "residential", "track")

osm_major <- osm_lines %>% filter(highway %in% major_highways)
osm_cycle <- osm_lines %>% filter(highway == "cycleway")
osm_trail <- osm_lines %>% filter(highway %in% trail_types)

camp_sites <- osm_points %>% 
  filter(grepl("tourism", other_tags) & grepl("camp_site", other_tags))

# Transform CRS to grid CRS -----------------------------------
osm_major <- st_transform(osm_major, crs_proj)
osm_cycle <- st_transform(osm_cycle, crs_proj)
osm_trail <- st_transform(osm_trail, crs_proj)
camp_sites <- st_transform(camp_sites, crs_proj)

# ---------- Intersection with grid + length densities ----------
compute_density <- function(line_sf, grid, id_col = "id") {
  if (nrow(line_sf) == 0) {
    return(grid %>% mutate(density = 0))
  }
  inter <- st_intersection(line_sf, grid)
  inter$len_km <- drop_units(set_units(st_length(inter), "km"))
  dens <- inter %>%
    st_drop_geometry() %>%
    group_by(!!sym(id_col)) %>%
    summarise(total_km = sum(len_km, na.rm=TRUE))
  grid %>%
    left_join(dens, by = id_col) %>%
    mutate(total_km = replace_na(total_km, 0),
           area_km2 = drop_units(set_units(st_area(.), "km^2")),
           density = total_km / area_km2)
}

major_test <- compute_density(osm_major, catalunya_grd) %>% 
  rename(major_road_density = density) %>%
  st_set_geometry(NULL) %>%
  dplyr::select(all_of(c("id", "major_road_density")))

catalunya_grd <- catalunya_grd %>% 
  mutate(id = as.character(id))

catalunya_grd <- catalunya_grd %>%
  left_join(compute_density(osm_major, catalunya_grd) %>% 
              rename(major_road_density = density) %>%
              st_set_geometry(NULL) %>%
              dplyr::select(all_of(c("id", "major_road_density"))),
            by = "id") %>%
  left_join(compute_density(osm_cycle, catalunya_grd) %>% 
              rename(cycleway_density = density) %>%
              st_set_geometry(NULL) %>%
              dplyr::select(all_of(c("id", "cycleway_density"))),
            by = "id") %>%
  left_join(compute_density(osm_trail, catalunya_grd) %>% 
              rename(trail_density = density) %>%
              st_set_geometry(NULL) %>%
              dplyr::select(all_of(c("id", "trail_density"))),
            by = "id")

catalunya_grd <- catalunya_grd %>%
  mutate(across(contains("density"), ~replace_na(.x, 0)))

# ---------- Distance to nearest campsite ----------
centroids <- st_centroid(catalunya_grd)
dist_mat <- st_distance(centroids, camp_sites)  # meters
catalunya_grd$dist_camp_km <- apply(dist_mat, 1, min) / 1000   # km

# ---------- Extract topography & landcover to grid ----------
rast_topo <- rast(topo_path)
rast_lc   <- rast(landcover_path)

topo_vals <- exact_extract(rast_topo, catalunya_grd, "mean")
lc_vals   <- exact_extract(rast_lc, catalunya_grd, "frac")

catalunya_grd <- bind_cols(catalunya_grd, topo_vals, lc_vals)

catalunya_grd <- catalunya_grd %>%
  mutate(across(starts_with("mean."), ~replace_na(.x, 0)),
         across(starts_with("frac_"), ~replace_na(.x, 0)))

# ---------- Rasterize all covariates ----------
grid_vect <- vect(catalunya_grd)
r_template <- rast(ext(grid_vect), resolution = 500, crs = crs(grid_vect))

covars <- c(
  rasterize(grid_vect, r_template, field = "major_road_density"),
  rasterize(grid_vect, r_template, field = "cycleway_density"),
  rasterize(grid_vect, r_template, field = "trail_density"),
  rasterize(grid_vect, r_template, field = "dist_camp_km"),
  rasterize(grid_vect, r_template, field = "mean.elev"),
  rasterize(grid_vect, r_template, field = "mean.slope"),
  rasterize(grid_vect, r_template, field = "mean.aspect"),
  rasterize(grid_vect, r_template, field = "frac_1"),
  rasterize(grid_vect, r_template, field = "frac_2"),
  rasterize(grid_vect, r_template, field = "frac_3"),
  rasterize(grid_vect, r_template, field = "frac_4"),
  rasterize(grid_vect, r_template, field = "frac_5")
)

names(covars) <- c(
  "major_road_density",
  "cycleway_density",
  "trail_density",
  "campsite_distance_km",
  "elevation_mean",
  "slope_mean",
  "aspect_mean",
  "lc_frac_1",
  "lc_frac_2",
  "lc_frac_3",
  "lc_frac_4",
  "lc_frac_5"
)

writeRaster(covars, covariate_raster_out, overwrite = TRUE)

##### Predict fire ignition using covariates #####

# ============================================================
# Model wildfire ignition probability using covariate raster
# Method: Random Forest presence-background model
# Output: ignition_prob_rf.tif + metrics + variable importance
# ============================================================


# ---------- Load data ----------
covs <- rast(covariate_raster_out)
ign_sf <- st_centroid(wildfires[wildfires$area_ha<1000,])

# Ensure same CRS
ign_sf <- st_transform(ign_sf, st_crs(covs))

# ---------- Extract covariates at ignition points (presence = 1) ----------
pres_xy <- st_coordinates(ign_sf)
pres_vals <- terra::extract(covs, pres_xy)
pres_df <- cbind(presence = 1, pres_vals)

# ---------- Sample background points (pseudo-absence) ----------
set.seed(123)
n_bg <- max(nrow(pres_df) * 2, 2000)  # 2x presence points or minimum 2000

bg_xy <- spatSample(x = covs, size = n_bg, method = "random", as.points = TRUE, na.rm = TRUE)
bg_vals <- terra::extract(covs, bg_xy, ID = F)
bg_df <- cbind(presence = 0, bg_vals)

# ---------- Combine & clean ----------
dat <- rbind(pres_df, bg_df)

coords <- rbind(
  cbind(pres_vals, pres_xy) %>% as.data.frame(),
  cbind(bg_vals, st_coordinates(st_as_sf(bg_xy))) %>% as.data.frame()
)

dat$X <- coords$X
dat$Y <- coords$Y

dat <- dat %>% filter(complete.cases(.))

# Optional: remove highly correlated features (if needed)
dat <- dat %>% dplyr::select(-aspect_mean)

# ---------- Train/Test Split ----------
set.seed(42)
idx <- sample(seq_len(nrow(dat)), 0.7 * nrow(dat))
train <- dat[idx, !names(dat) %in% c("X", "Y")]
test  <- dat[-idx, !names(dat) %in% c("X", "Y")]

# ---------- Train Random Forest ----------
rf_model <- randomForest(
  factor(presence) ~ ., 
  data = train,
  ntree = 2000,
  mtry = floor(sqrt(ncol(train) - 1)),
  importance = TRUE
)

print(rf_model)

# ---------- Evaluate Model ----------
# Predictions
test$pred <- predict(rf_model, test, type = "prob")[,2]

# AUC
auc_val <- roc(test$presence, test$pred)$auc
cat("AUC:", auc_val, "\n")

# Confusion Matrix (using best threshold / Youden J)
opt_th <- as.numeric(coords(roc(test$presence, test$pred), "best")["threshold"])
test$pred_class <- ifelse(test$pred >= opt_th, 1, 0)
table(True = test$presence, Pred = test$pred_class)

# Variable importance
varImpPlot(rf_model, main="Variable Importance - RF Ignition Model")

# ---------- Predict over region ----------
rf_pred <- predict(covs, rf_model, type = "prob")$X1

# Scale to [0,1]
rf_pred <- clamp(rf_pred, lower = 0, upper = 1)

names(rf_pred) <- "ignition_prob_rf"

# Save raster
writeRaster(rf_pred, out_rast_path, overwrite = TRUE)


par(mfrow=c(1,1))
plot(rf_pred)

# Boosted regression trees
library(ranger)
library(xgboost)


# Make sf
dat_sf <- st_as_sf(dat, coords = c("X", "Y"), crs = st_crs(ign_sf))

# Create spatial blocks (e.g. 10 km)
set.seed(123)
cv_blocks <- spatialBlock(
  speciesData = dat_sf,
  species = "presence",
  theRange = 10000,   # 10 km blocks
  k = 5,              # 5-fold CV
  selection = "random",
  iteration = 100,
  showBlocks = FALSE
)


cv_auc <- c()

for (fold in 1:cv_blocks$k) {
  
  train_idx <- which(cv_blocks$foldID != fold)
  test_idx  <- which(cv_blocks$foldID == fold)
  
  train <- dat[train_idx,]
  test  <- dat[test_idx,]
  
  rf <- ranger(
    presence ~ ., 
    data = train,
    importance = "impurity",
    probability = TRUE,
    num.trees = 700,
    mtry = floor(sqrt(ncol(train)-1))
  )
  
  test$pred <- predict(rf, test)$predictions[,2]
  cv_auc[fold] <- roc(test$presence, test$pred)$auc
}

mean(cv_auc)







# 
# X <- as.matrix(dat[, setdiff(names(dat), c("presence","X", "Y" ))])
# y <- dat$presence
# # store feature names
# feat_names <- colnames(X)
# 
# # Convert to DMatrix
# dtrain <- xgb.DMatrix(data = X, label = y)
# 
# params <- list(
#   booster = "gbtree",
#   objective = "binary:logistic",
#   eval_metric = "auc",
#   eta = 0.05,
#   max_depth = 5,
#   subsample = 0.8,
#   colsample_bytree = 0.8
# )
# 
# set.seed(123)
# xgb_model <- xgb.train(
#   params = params,
#   data = dtrain,
#   nrounds = 500,
#   watchlist = list(train = dtrain),
#   verbose = 0
# )
# 
# # In-sample check
# pred_xgb <- predict(xgb_model, X)
# auc_xgb <- roc(y, pred_xgb)$auc
# auc_xgb
# 
# # 1. Get the feature importance table
# importance_matrix <- xgb.importance(feature_names = feat_names, model = xgb_model)
# 
# # 2. (Optional) Plot the top N features
# xgb.plot.importance(importance_matrix, top_n = 10)
# 
# 
# 
# cv_auc_xgb <- c()
# 
# for (fold in 1:cv_blocks$k) {
#   train <- dat[cv_blocks$foldID != fold,]
#   test  <- dat[cv_blocks$foldID == fold,]
#   
#   Xtr <- as.matrix(train[, setdiff(names(train), "presence")])
#   ytr <- train$presence
#   Xte <- as.matrix(test[, setdiff(names(test), "presence")])
#   yte <- test$presence
#   
#   dtrain <- xgb.DMatrix(Xtr, label = ytr)
#   
#   model <- xgb.train(params = params, data = dtrain, nrounds = 500, verbose = 0)
#   pred <- predict(model, Xte)
#   
#   cv_auc_xgb[fold] <- roc(yte, pred)$auc
# }
# 
# mean(cv_auc_xgb)
# 
# 
# # convert raster to matrix
# covs_df <- as.data.frame(covs) %>%
#   filter(complete.cases(.))
# 
# # enforce same names and column order
# covs_df <- covs_df[, feat_names, drop = FALSE]
# 
# dmat <- xgb.DMatrix(data = as.matrix(covs_df))
# pred_vec <- predict(xgb_model, dmat)
# 
# # pred_vec must be your XGBoost prediction vector (numeric)
# pred_r <- rast(covs[[1]])
# non_na_indices <- which(!is.na(values(covs[[1]])))
# 
# # Initialize the prediction raster with NAs everywhere
# values(pred_r) <- NA
# 
# # Assign the predicted values to the non-NA cells
# # The length of pred_vec must match the length of non_na_indices
# values(pred_r)[non_na_indices] <- pred_vec
# 
# # Rename and plot
# names(pred_r) <- "ignition_prob_xgb"
# 
# plot(pred_r)
# 
# writeRaster(pred_r, here("data", "rdat", "ignition_prob_xgb.tif"), overwrite = TRUE)




library(mgcv)

gam_model <- gam(presence ~ s(elevation_mean) + s(slope_mean) + s(campsite_distance_km) +
                   s(major_road_density) + s(trail_density) + 
                   lc_frac_1 + lc_frac_2 + lc_frac_3 + lc_frac_4 + lc_frac_5,
                 data = train, family = binomial)

summary(gam_model)
# ---------- Predict over region ----------
gam_pred <- predict(covs, gam_model, type = "response")

# Scale to [0,1]
#gam_pred <- clamp(gam_pred, lower = 0, upper = 1)

names(gam_pred) <- "ignition_prob_gam"


plot(gam_pred)
plot(rf_pred)





##### Seconda dataset (comarca, not spatial explicit) #####




mapview::mapview(catalunya_grd)

# load wildfire data, calculate centroids, 
wildfires <- list.files("./data/wildfires/wildfires_shp", pattern = ".shp$", 
                        recursive = TRUE, full.names = TRUE, include.dirs = TRUE) |>
  lapply(st_read) |>
  bind_rows()

wildfires <- wildfires |>
  filter(!is.na(CODI_FINAL)) |>
  mutate(
    year = case_when(nchar(DATA_INCEN) == 7 ~ as.integer(paste0("20", substr(DATA_INCEN, 6, 7))),
                     nchar(DATA_INCEN) == 8 ~ as.integer(paste0("20", substr(DATA_INCEN, 7, 8))),
                     nchar(DATA_INCEN) >= 10 ~ as.integer(substr(DATA_INCEN, 7, 10)),
                     TRUE ~ NA),
    month = case_when(nchar(DATA_INCEN) == 7 ~ as.integer(substr(DATA_INCEN, 3, 4)),
                      nchar(DATA_INCEN) >= 8 ~ as.integer(substr(DATA_INCEN, 4, 5)),
                      TRUE ~ NA),
    # day is anything before first slash
    day = as.integer(substring(DATA_INCEN, 1, regexpr("/", DATA_INCEN) - 1)),
    date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d"),
    area_ha = as.numeric(st_area(geometry)) / 10000,  # convert from m² to ha
    MUNICIPI = str_to_title(MUNICIPI)  # capitalize first letter of each word
  ) |>
  select(CODI_FINAL, MUNICIPI, date, year, month, area_ha, geometry)


# transform crs
wildfires <- st_transform(wildfires, crs = st_crs(grd)) %>%
  st_make_valid()

wildfire_centroids <- st_centroid(wildfires)

mapview::mapview(wildfire_centroids)

# join with grid data
wf_dat <- st_join(catalunya_grd, wildfire_centroids, join = st_intersects) |>
  mutate(any_wildfire = ifelse(!is.na(CODI_FINAL), 1, 0))

table(wf_dat$any_wildfire)

# bivariate histogram
library(ggplot2)

# histogram of total burned area
wf_df %>%
  filter(year >=2015) %>%
  filter(COMARCA %in% com_aoi) %>%
  filter(HAFORESTAL>4) %>%
  ggplot(aes(x = HAFORESTAL)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of Total Burned Area (ha)",
       x = "Total Burned Area (ha)",
       y = "Frequency") +
  theme_minimal()

# histogram of simulated annual fire frequency
ggplot(sim_dt, aes(x = avg_fire_duration)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
  labs(title = "Distribution of Simulated Annual Fire Duration",
       x = "Avg. Fire Duration",
       y = "Frequency") +
  theme_minimal()

ggplot(sim_dt, aes(x = annual_fire_frequency)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
  labs(title = "Distribution of Simulated Annual Fire Frequency",
       x = "Annual Fire Frequency",
       y = "Frequency") +
  theme_minimal()


ggplot(sim_dt, aes(x = annual_fire_frequency, y = avg_fire_duration)) +
  geom_bin2d(bins = 300) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Bivariate Histogram of Fire Frequency and Duration",
       x = "Annual Fire Frequency",
       y = "Average Fire Duration (days)") +
  theme_minimal()
