# BCN raster exploration
# Distribution of size and cost of SAL interventions across years

library(tidyverse)
library(sf)
library(terra)

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")


##### Test targeting with average scenario ####

# load grd data with average extracts
grd <- st_read("./data/rdat/grd.gpkg")
grd_ctrd <- st_centroid(grd)

annual_mask_rast_dir <- "./data/Re_ Omnsicape Algorihm/LOG_rasters"
annual_ros_rast_dir <- "./data/Re_ Omnsicape Algorihm/ROS_rasters"

# load and c rasters 
annual_mask_rast <- list.files(annual_mask_rast_dir, pattern = ".tif$", full.names = TRUE) |>
  lapply(terra::rast) |>
  do.call(c, args = _)

annual_ros_rast <- list.files(annual_ros_rast_dir, pattern = ".tif$",full.names = TRUE) %>%
  lapply(terra::rast) %>%
  do.call(c, .) 

# Assign names to the SpatRaster object using base R's names() function
names(annual_mask_rast) <- gsub(".tif", "", list.files(annual_mask_rast_dir, pattern = ".tif$"))
names(annual_ros_rast) <- gsub(".tif", "", list.files(annual_ros_rast_dir, pattern = ".tif$"))

annual_maps <- c(annual_mask_rast, annual_ros_rast)

# extract mask and ros values for each year to grd_ctrd
e1 <- terra::extract(annual_maps, grd_ctrd)

# replace NA by 0

# combine with grd
target_df <- cbind(grd, e1)

# remove all but target_df
rm(list = setdiff(ls(), "target_df"))
gc()
source("./src/bcn_targeting.R")


# iterate over thresholds in 10% increase
years <- 2040:2050
thresholds <- seq(0.1, 1, by = 0.1)
strategies <- c("random", "road_distance", "lc_artificial", "biomass_bau", "connectivity_bau")
order_decreasing <- c(FALSE, FALSE, TRUE, TRUE, TRUE)


for (y in years) {
  cat(paste0("Year: ", y, "\n"))
  for (s in 1:length(strategies)) {
    cat(paste0("Target: ", strategies[s], " (decr.:", order_decreasing[s], ") \n"))
    for (i in thresholds) {
      target_df <- BCN_target_df(
        data = target_df,
        original_col = paste0("bau_ROS_", y),
        mask_col = paste0("Logging_", y),
        replacement_col = paste0("sal_ROS_", y),
        target_col = strategies[s], # Optional: if not provided, random assignment is used
        threshold = i,
        budget = NULL,
        relative = TRUE,
        year = y,
        target_decreasing = order_decreasing[s],
        verbose = FALSE
      )
    }
  }
}

saveRDS(target_df, file = "./data/rdat/target_df.rds")
# load target_df
target_df <- readRDS(file = "./data/rdat/target_df.rds")

# --- Reshape the dataframe from wide to long format ---
target_df_long = target_df %>%
  st_set_geometry(NULL) %>%
  pivot_longer(
    cols = matches("^(new_mask|modified_values)_[a-zA-Z0-9_]+_t\\d+_y\\d+$"), # Select columns matching the new pattern
    names_to = c(".value", "strategy", "threshold", "year"), # Use .value, and new 'strategy' and 'threshold' columns
    names_pattern = "^(new_mask|modified_values)_([a-zA-Z0-9_]+)_t(\\d+)_y(\\d+)$" # Regex to capture type, strategy, and threshold
  ) %>%
  mutate(
    threshold = as.numeric(threshold), # Convert the 'threshold' column to numeric
    cost_scaler = road_distance / max(road_distance),
    sal_costs_ha = 1750 + (400 * cost_scaler),
    sal_ha = new_mask * lc_wildland * 4, # Convert to hectares (assuming each cell is 4 ha)
    sal_costs_cell = sal_costs_ha  * sal_ha
  )

# aggregate by threshold (mean and sum)
target_df_agg <- target_df_long %>%
  #st_set_geometry(NULL) %>%
  group_by(strategy, threshold, year) %>%
  summarise(
    sum_SAL_px = sum(new_mask, na.rm = TRUE),
    mean_SAL_ha = mean(sal_ha, na.rm = TRUE),
    sum_SAL_ha = sum(sal_ha, na.rm = TRUE),
    mean_SAL_costs = mean(sal_costs_cell, na.rm = TRUE),
    sum_SAL_costs = sum(sal_costs_cell, na.rm = TRUE)
  ) |>
  # rename strategies to more readable
  mutate(strategy = case_when(
    strategy == "random" ~ "Random",
    strategy == "road_distance" ~ "Lowest road distance",
    strategy == "lc_artificial" ~ "Highest artificial",
    strategy == "biomass_bau" ~ "Highest biomass first",
    strategy == "connectivity_bau" ~ "Highest connectivity"
  ))


# maximum sum_SAL_costs per year
target_df_agg$max_cost <- ave(target_df_agg$sum_SAL_costs, target_df_agg$year, FUN = max)
target_df_agg$max_cost_scaled <- round(target_df_agg$max_cost / max(target_df_agg$max_cost), 4)
target_df_agg$max_sal_area_km2 <- ave(target_df_agg$sum_SAL_ha, target_df_agg$year, FUN = max) / 100 # convert ha to km2

target_df_agg$facet_lab <- paste0(target_df_agg$year, " (max. area: ", round(target_df_agg$max_sal_area_km2, 2), " km²)")

# plot mask share (Y) vs threshold (X)
p.cum_sal <- ggplot(target_df_agg, aes(x = threshold, y = sum_SAL_costs, colour = strategy)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = max_cost_scaled),
            inherit.aes = FALSE, alpha = 1) + # alpha = 1 ensures full fill
  scale_fill_gradient(low = "white", high = "grey60", guide = "none") + # Lightest to darkest grey
  geom_line(lwd = 1.2) +
  labs(
    x = "Percentage of area salvage-logged",
    y = "Costs of salvage logging (EUR)"
  ) +
  theme_minimal()+
  # Format y-axis labels to show values in millions
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  # facet grid by year
  facet_wrap(~ facet_lab, scales = "free_y", ncol = 4) + # Added ncol = 4
  theme(
    legend.position = c(0.975, 0.025), # Placed legend in the approximate 12th facet position (bottom-right)
    legend.justification = c("right", "bottom") # Adjust justification to place it correctly within the box
  )

p.cum_sal

ggsave(p.cum_sal, filename = "./out/fig/cum_sal_costs.png", 
       width = 12, height = 8, dpi = 300, bg = "white")

# alternative: express in relative terms cmopared to random scenario
target_df_rdm <- target_df_agg %>%
  filter(strategy == "Random") %>%
  # drop groups
  ungroup() %>%
  select(year, threshold, sum_SAL_costs) %>%
  rename(rdm_sum_SAL_costs = sum_SAL_costs)

target_df_rel <- target_df_agg %>%
  filter(strategy != "Random") %>%
  left_join(target_df_rdm, by = c("year", "threshold")) %>%
  mutate(
    rel_sum_SAL_costs = (sum_SAL_costs / rdm_sum_SAL_costs) - 1 
  ) %>%
  select(-rdm_sum_SAL_costs)

# make area plot
ggplot(target_df_rel, aes(x = threshold, y = rel_sum_SAL_costs, color = strategy)) +
  geom_line(lwd = 1.2) +
  # add horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Percentage of area salvage-logged",
    y = "Relative costs of salvage logging (compared to random scenario)"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") + # Use a color palette for better visibility
  facet_wrap(~ year, scales = "free_y", ncol = 4) + # Added ncol = 4
  theme(
    legend.position = c(0.975, 0.025), # Placed legend in the approximate 12th facet position (bottom-right)
    legend.justification = c("right", "bottom") # Adjust justification to place it correctly within the box
  )

# save relative plot
ggsave("./out/fig/cum_sal_costs_relative.png", width = 12, height = 8, dpi = 300, bg = "white")


# make area plot
p.rel.smooth <- ggplot(target_df_rel, aes(x = threshold, y = rel_sum_SAL_costs, color = strategy)) +
  geom_point(alpha = 0.2) +
  # smooth line aggregating years
  geom_smooth(method = "loess", se = TRUE, lwd = 1.2) +
  # add horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Percentage of area salvage-logged",
    y = "Relative costs of salvage logging (compared to random scenario)"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") + # Use a color palette for better visibility
  #facet_wrap(~ year, scales = "free_y", ncol = 4) + # Added ncol = 4
  theme(
    legend.position = c(0.975, 0.025), # Placed legend in the approximate 12th facet position (bottom-right)
    legend.justification = c("right", "bottom") # Adjust justification to place it correctly within the box
  )

# save smooth plot
ggsave(p.rel.smooth, filename = "./out/fig/cum_sal_costs_relative_smooth.png", 
       width = 6, height = 4, dpi = 300, bg = "white")


#### simulate burning patterns across scenarios and starting locations ####

##### Wildfire Simulation  #####
spread_fun = function (advancement_level, speed, attraction){
  advancement_level + (speed * attraction)
}
# burn function (now takes scenario name instead of suffix
get_burners_optimized <- function(time_horizon,
                                  ignition_cell,
                                  cell_size = 200, # in meters
                                  time_step = 15,    # in minutes
                                  scenario = "bau", # select scenario (bau, sal)
                                  full = FALSE, # return full data frame or last column
                                  neighbor_list_differentiated, # Pre-processed neighbor list
                                  grd, # Pass grd directly to avoid global variable dependency
                                  verbose = FALSE # Verbose output for debugging
) {
  
  # Define variable names based on scenario
  # speed_var <- paste0("biomass_", scenario)
  # connectivity_var <- paste0("connectivity_", scenario)
  speed_var <- scenario
  connectivity_var <- "connectivity_bau"
  
  grd_geo <- grd
  grd_attr <- grd %>% st_set_geometry(NULL) # Attributes only for faster access
  
  # Prepare empty grid for simulation results
  n_cells <- nrow(grd_attr) # Number of cells
  bdf <- matrix(0, nrow = n_cells, ncol = time_horizon + 1) # Initialize as a matrix
  rownames(bdf) <- row.names(grd_attr) # Set row names for easier lookup
  
  # Find the row index of the ignition cell
  ignition_idx <- which(row.names(grd_attr) == ignition_cell)
  bdf[ignition_idx, 1] <- 1  # Set ignition cell
  
  # Pre-calculate speed and attraction vectors
  speed_values <- grd_attr[[speed_var]]
  attraction_values <- grd_attr[[connectivity_var]]
  
  # Pre-calculate distance factors for all cells based on their type (rook or queen from origin)
  # This part is more complex as it depends on the origin cell, so we'll handle it within the loop,
  # but by pre-indexing.
  
  neighbor_idx_differentiated <- neighbor_list_differentiated
  # Convert neighbor_list_differentiated to use indices directly for faster access
  # neighbor_idx_differentiated <- lapply(neighbor_list_differentiated, function(cell_neighbors) {
  #   list(
  #     rook = which(row.names(grd_attr) %in% cell_neighbors$rook),
  #     queen = which(row.names(grd_attr) %in% cell_neighbors$queen)
  #   )
  # })
  #names(neighbor_idx_differentiated) <- row.names(grd_attr) # Maintain names
  
  for (t in 2:(time_horizon + 1)) {
    prev_t <- t - 1
    
    # Identify currently burning cells (indices)
    burning_idx <- which(bdf[, prev_t] >= 1)
    
    # Carry over burning status
    bdf[burning_idx, t] <- 1
    
    # Initialize a temporary vector to store burn pressure for potential burners in this time step
    # This will accumulate pressure from all burning origin cells
    current_time_step_burn_pressure <- rep(0, n_cells)
    
    if (length(burning_idx) > 0) { # Only proceed if there are burning cells
      for (origin_cell_idx in burning_idx) {
        # Get differentiated neighbor indices for the current origin cell
        current_cell_neighbors_idx <- neighbor_idx_differentiated[[origin_cell_idx]] # Access by index directly
        
        rook_neighbors_idx <- current_cell_neighbors_idx$rook
        queen_neighbors_idx <- current_cell_neighbors_idx$queen
        
        # Process Rook Neighbors
        if (length(rook_neighbors_idx) > 0) {
          # Filter out already burning cells
          potential_rook_burners_idx <- rook_neighbors_idx[!(rook_neighbors_idx %in% burning_idx)]
          
          if (length(potential_rook_burners_idx) > 0) {
            advancement_level <- bdf[potential_rook_burners_idx, prev_t]
            speed <- speed_values[potential_rook_burners_idx]
            attraction <- attraction_values[potential_rook_burners_idx]
            
            # Rook neighbors have a distance factor of 1
            speed_per_time_step_adjusted <- (speed / cell_size) * time_step
            burn_pressure <- spread_fun(advancement_level, speed_per_time_step_adjusted, attraction)
            
            # Accumulate burn pressure
            current_time_step_burn_pressure[potential_rook_burners_idx] <-
              pmax(current_time_step_burn_pressure[potential_rook_burners_idx], burn_pressure)
          }
        }
        
        # Process Queen (Exclusive) Neighbors
        if (length(queen_neighbors_idx) > 0) {
          # Filter out already burning cells
          potential_queen_burners_idx <- queen_neighbors_idx[!(queen_neighbors_idx %in% burning_idx)]
          
          if (length(potential_queen_burners_idx) > 0) {
            advancement_level <- bdf[potential_queen_burners_idx, prev_t]
            speed <- speed_values[potential_queen_burners_idx]
            attraction <- attraction_values[potential_queen_burners_idx]
            
            # Queen neighbors have a distance factor of sqrt(2)
            distance_factor_queen <- sqrt(2)
            speed_per_time_step_adjusted <- (speed / (cell_size * distance_factor_queen)) * time_step
            burn_pressure <- spread_fun(advancement_level, speed_per_time_step_adjusted, attraction)
            
            # Accumulate burn pressure
            current_time_step_burn_pressure[potential_queen_burners_idx] <-
              pmax(current_time_step_burn_pressure[potential_queen_burners_idx], burn_pressure)
          }
        }
      }
    }
    # Apply the accumulated burn pressure to the bdf for the current time step
    # Only update cells that received some pressure and are not yet fully burned
    cells_to_update <- which(current_time_step_burn_pressure > 0)
    bdf[cells_to_update, t] <- pmax(bdf[cells_to_update, t], current_time_step_burn_pressure[cells_to_update])
    
    # Ensure advancement level doesn't exceed 1 for the previous time step (already processed)
    bdf[, prev_t] <- pmin(1, bdf[, prev_t])
  }
  
  # Convert back to data frame for merging
  bdf <- as.data.frame(bdf)
  
  last_col_name <- tail(names(bdf), 1)
  
  if (full) {
    bdf <- bdf %>%
      filter(!!sym(last_col_name) > 0)
    return(bdf)
  } else {
    # merge with spatial data, filter where last column is larger 0
    bdf <- bind_cols(grd_geo, bdf) %>%
      filter(!!sym(last_col_name) > 0) %>%
      dplyr::select(all_of(c("id", last_col_name))) %>%
      rename(burning = !!sym(last_col_name))
    
    return(bdf)
  }
}

# load differentiated neighbor list
# load grd
#grd <- st_read("./data/rdat/grd.gpkg")
#row.names(grd) <- paste0("cell_", 1:nrow(grd))
#load("./data/rdat/nb_list_rq_short.RData")
target_df <- readRDS(file = "./data/rdat/target_df.rds")
row.names(target_df) <- paste0("cell_", 1:nrow(target_df))

# takes 5-10 minutes to run
neighbor_idx_differentiated <- lapply(neighbor_list, function(cell_neighbors) {
  list(
    rook = which(row.names(grd) %in% cell_neighbors$rook),
    queen = which(row.names(grd) %in% cell_neighbors$queen)
  )
})
names(neighbor_idx_differentiated) <- row.names(grd)

# save the neighbor list for later use
save(neighbor_idx_differentiated, file = "./data/rdat/neighbor_idx_differentiated.RData")
load("./data/rdat/neighbor_idx_differentiated.RData")



# run lapply in parallel using future apply and show progress bar
library(future)
library(future.apply)
library(progressr)

# calculate 850mb limit:
# 850*1024^2 = 891289600
options(future.globals.maxSize= 891289600, 
        future.seed = TRUE)

handlers(handler_progress(format="[:bar] :percent :eta :message"))

plan(multisession, workers=2)
rdm_ign_cells <- sample(x = row.names(target_df), 
                        size = 30, 
                        replace = FALSE,
                        prob = target_df$ignition_probability_surface)
mc_time_step <- 5
mc_time_horizon <- 12
# Define the scenarios
scenarios <- c("bau_ROS_2045",  
               #"sal_ROS_2045",
               #"modified_values_lc_artificial_t50_y2045", 
               #"modified_values_connectivity_bau_t50_y2045", 
               #"modified_values_random_t50_y2045",
               "modified_values_road_distance_t50_y2045") # Add all your scenarios here

# --- Efficiently loop over scenarios ---
all_scenario_results <- list() # To store results from each scenario

aggregate_cells <- FALSE

with_progress({
  p_scenarios <- progressor(along = scenarios)
  
  # Outer loop: iterate over scenarios
  for (s_name in scenarios) {
    p_scenarios(message = paste("Processing scenario:", s_name)) # Progress for scenarios
    
    # Inner parallel loop: iterate over ignition cells for the current scenario
    mc_res_current_scenario <- future_lapply(rdm_ign_cells, function(x) {
      # You can add a progressor here if you want a nested progress bar,
      # but it might get cluttered. The outer progress bar is usually sufficient.
      # p_cells() # Uncomment and define p_cells if you want inner progress
      
      get_burners_optimized(x,
                            time_horizon = mc_time_horizon,
                            time_step = mc_time_step,
                            full = TRUE,
                            scenario = s_name, # Use the current scenario name
                            grd = target_df,
                            neighbor_list_differentiated = neighbor_idx_differentiated)
    })
    
    
    
    if(aggregate_cells){
      # Aggregate results for the current scenario
      agg_current_scenario <- mc_res_current_scenario %>%
        lapply(., function(x) apply(x, 2, sum)) %>%
        bind_rows() %>% # Use bind_rows directly
        as.data.frame() %>%
        mutate(scenario = s_name) # Assign the correct scenario name
      
      # Store the aggregated results for the current scenario
      all_scenario_results[[s_name]] <- agg_current_scenario
    } else if(!aggregate_cells){
      all_scenario_results[[s_name]] <- mc_res_current_scenario
    }
  }
})

long_df_raw <- map_df(names(all_scenario_results), function(scenario_name) {
  
  # Get the list of dataframes for the current scenario
  scenario_dfs <- all_scenario_results[[scenario_name]]
  
  # Process each dataframe within the scenario
  map_df(scenario_dfs, function(df) {
    
    ignition_id <- row.names(df[df$V1 == 1,])  # Extract ignition_id from the first row
    
    df %>%
      tibble::rownames_to_column(var = "cell_id") %>% # Convert row names to a column
      mutate(scenario = scenario_name, 
             ignition_id = ignition_id) %>%                # Add the scenario name
      pivot_longer(-c(cell_id, scenario, ignition_id), names_to = "time", values_to = "burned") %>%
      mutate(time = as.numeric(gsub("V", "", time)) * mc_time_step)
  })
})

# save for later merging with land cover data
saveRDS(long_df_raw, file = "./data/rdat/mc_res_raw.rds")
long_df_raw <- readRDS(file = "./data/rdat/mc_res_raw.rds")

# merge with lc data
lc_dat <- target_df %>%
  st_set_geometry(NULL) %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  select(all_of(c("cell_id", "lc_agriculture", "lc_wildland","lc_artificial")))



long_df_agg <- long_df_raw %>%
  left_join(lc_dat, by = "cell_id") %>%
  group_by(scenario, ignition_id, time) %>%
  summarise(burned_land_ha = sum(burned * 4, na.rm = TRUE),
            burned_agri_ha = sum(burned * lc_agriculture * 4, na.rm = TRUE),
            burned_wild_ha = sum(burned * lc_wildland * 4, na.rm = TRUE),
            burned_urbn_ha = sum(burned * lc_artificial * 4, na.rm = TRUE), 
            .groups = 'drop') %>%
  ungroup()



saveRDS(long_df_agg, file = "./data/rdat/mc_res.rds")
long_df_agg <- readRDS(file = "./data/rdat/mc_res.rds")

# ridge plot using full distributions of long_df_agg
library(ggridges)
library(tidyverse)
library(ggplot2)
p.mc.ridge <- ggplot(long_df_agg |> 
                       filter(time %in% c(10, 30, 60, 180),
                              !grepl("mask", scenario)), 
                     aes(x = burned_land_ha, y = factor(scenario))) +
  geom_density_ridges(scale = 0.9,
                      quantile_lines = TRUE, alpha = 0.75,
                      quantiles = c(0.05, 0.5, 0.95)) +
  labs(x = "Burned hectares", y = "Time (minutes)") +
  theme_minimal() +
  theme(legend.position = "right") +
  # log-scale x axis
  scale_x_log10() +
  # facet by scenario
  facet_wrap(~ factor(time), scales = "free_x")


p.mc.ridge



# 2. Calculate time-scenario specific summaries


# same as above, but summarise all burned_land_ha, burned_agri_ha, burned_wild_ha, burned_urbn_ha
scenario_summaries <- long_df_agg %>%
  group_by(time, scenario) %>%
  summarise_at(
    vars(burned_land_ha, burned_agri_ha, burned_wild_ha, burned_urbn_ha),
    list(
      mean = ~ mean(.x, na.rm = TRUE),
      median = ~ median(.x, na.rm = TRUE),
      min = ~ min(.x, na.rm = TRUE),
      max = ~ max(.x, na.rm = TRUE),
      percentile95 = ~ quantile(.x, 0.95, na.rm = TRUE),
      percentile5 = ~ quantile(.x, 0.05, na.rm = TRUE)
    ),
    .groups = 'drop' # Important to drop grouping after summarise
  ) %>%
  # long format
  pivot_longer(
    cols = -c(time, scenario),
    names_to = c("metric", "type", "unit", "statistic"),
    names_sep = "_",
    values_to = "value"
  )



# 3. Calculate the difference to the 'biomass_bau' scenario

# First, isolate the 'biomass_bau' summaries
bau_summaries <- scenario_summaries %>%
  filter(scenario == "bau_ROS_2045") %>%
  rename(value_bau = "value")

scenario_differences <- scenario_summaries %>%
  # Join by time, land_cover, and statistic to match the correct BAU value for comparison
  left_join(bau_summaries, by = c("time", "scenario","metric", 
                                  "type", "unit", "statistic")) %>%
  # Filter out the BAU scenario itself, as we're calculating differences *from* it
  filter(scenario != "bau_ROS_2045") %>%
  # Calculate the difference
  mutate(
    difference = value - value_bau
  ) %>%
  # Select only the relevant columns for plotting differences
  select(time, scenario, type, statistic, difference)

# You can now use 'scenario_differences' for plotting.
# Each row represents a specific time, scenario, land_cover type, and statistic,
# with its calculated difference from the BAU scenario.

ggplot(scenario_differences, aes(x = time, y = difference, color = scenario)) +
  geom_line() +
  facet_grid(statistic ~ land_cover, scales = "free_y") + # Separate plots for each statistic and land cover
  labs(
    title = "Difference in Burned Land Cover Metrics vs. BAU",
    y = "Difference from BAU (ha)",
    x = "Time"
  ) +
  theme_minimal()

p.differences
