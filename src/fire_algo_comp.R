library(data.table)
library(tidyverse)
library(future)
library(future.apply)
library(progressr)


#### Functions ####
# spread_fun remains the same
spread_fun = function(advancement_level, speed, attraction){
  advancement_level + (speed * attraction)
}

# burn function (now takes scenario name instead of suffix
get_burners_optimized <- function(time_horizon,
                                  ignition_cell,
                                  cell_size = 200, # in meters
                                  time_step = 15,    # in minutes
                                  scenario = "bau", # select scenario (bau, sal)
                                  full = FALSE, # return full data frame or last column
                                  neighbor_list_differentiated, # Pre-processed neighbor list (should use indices or names consistent with grd)
                                  grd, # Pass grd directly (now a normal data frame)
                                  con_scaler = 2000, # rescale unitless connectivity
                                  verbose = FALSE # Verbose output for debugging
) {
  
  # Input Validation for grd: Ensure it's a data.frame
  if (!inherits(grd, "data.frame")) {
    stop("Input 'grd' must be a data.frame object.")
  }
  
  # Store original row names as a column, this will serve as the 'id'
  original_row_names <- row.names(grd)
  
  # Define variable names based on scenario
  speed_var <- paste0("ros_trg_", scenario)
  connectivity_var <- paste0("con_trg_", scenario)
  
  # grd_attr is now simply grd
  grd_attr <- grd
  
  # Prepare empty grid for simulation results
  n_cells <- nrow(grd_attr) # Number of cells
  bdf <- matrix(0, nrow = n_cells, ncol = time_horizon + 1) # Initialize as a matrix
  rownames(bdf) <- original_row_names # Set row names for easier lookup
  
  # Find the row index of the ignition cell
  # Ensure ignition_cell exists in the row names
  if (!ignition_cell %in% original_row_names) {
    stop(paste0("Ignition cell '", ignition_cell, "' not found in row names of 'grd'."))
  }
  ignition_idx <- which(original_row_names == ignition_cell)
  bdf[ignition_idx, 1] <- 1  # Set ignition cell
  
  # Pre-calculate speed and attraction vectors
  # Ensure scenario-dependent columns exist
  if (!speed_var %in% names(grd_attr)) {
    stop(paste0("Speed variable column '", speed_var, "' not found in 'grd'."))
  }
  if (!connectivity_var %in% names(grd_attr)) {
    stop(paste0("Connectivity variable column '", connectivity_var, "' not found in 'grd'."))
  }
  speed_values <- grd_attr[[speed_var]]
  attraction_values <- grd_attr[[connectivity_var]] / con_scaler
  
  # Pre-calculate distance factors for all cells based on their type (rook or queen from origin)
  # This part is more complex as it depends on the origin cell, so we'll handle it within the loop,
  # but by pre-indexing.
  
  # neighbor_idx_differentiated expects indices matching bdf rows.
  # If neighbor_list_differentiated was built using row.names, ensure it's still accessible by index.
  # The original code's commented out section for converting names to indices might be useful
  # if neighbor_list_differentiated uses names instead of numeric indices.
  # For now, assuming neighbor_list_differentiated uses 1-based numeric indices directly
  # or can be directly indexed by original_row_names.
  
  # Ensure neighbor_list_differentiated is indexed by row names if it comes from string names.
  # If it's already indexed by 1-based integers (which is more efficient), no change needed.
  # The previous version commented out `lapply` conversion, assuming it's pre-converted.
  # Let's add a check and conversion if it's not indexed by numeric indices.
  if (is.list(neighbor_list_differentiated) && !is.null(names(neighbor_list_differentiated)) && !all(sapply(names(neighbor_list_differentiated), is.numeric))) {
    if (verbose) {
      message("Converting neighbor_list_differentiated from names to indices for faster access.")
    }
    # Create a mapping from row names to their 1-based numeric indices
    name_to_idx_map <- setNames(seq_along(original_row_names), original_row_names)
    
    neighbor_idx_differentiated <- lapply(neighbor_list_differentiated, function(cell_neighbors) {
      list(
        rook = name_to_idx_map[cell_neighbors$rook],
        queen = name_to_idx_map[cell_neighbors$queen]
      )
    })
    # Ensure all converted indices are valid
    if(anyNA(unlist(neighbor_idx_differentiated))) {
      stop("Some neighbor IDs in neighbor_list_differentiated do not match row names in 'grd'.")
    }
  } else {
    # Assume neighbor_list_differentiated is already indexed by integers (efficient)
    neighbor_idx_differentiated <- neighbor_list_differentiated
  }
  
  for (t in 2:(time_horizon + 1)) {
    prev_t <- t - 1
    
    # Identify currently burning cells (indices)
    burning_idx <- which(bdf[, prev_t] >= 1)
    
    # Carry over burning status
    bdf[burning_idx, t] <- 1
    
    # Initialize a temporary vector to store burn pressure for potential burners in this time step
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
    cells_to_update <- which(current_time_step_burn_pressure > 0)
    bdf[cells_to_update, t] <- pmax(bdf[cells_to_update, t], current_time_step_burn_pressure[cells_to_update])
    
    # Ensure advancement level doesn't exceed 1 for the previous time step (already processed)
    bdf[, prev_t] <- pmin(1, bdf[, prev_t])
  }
  
  # Convert back to data frame
  bdf <- as.data.frame(bdf)
  
  last_col_name <- tail(names(bdf), 1)
  
  # Add the original row names as the 'id' column
  bdf$id <- original_row_names
  
  if (full) {
    # Filter where the last column is greater than 0
    bdf <- bdf %>%
      filter(!!sym(last_col_name) > 0)
    
    # Reorder columns to have 'id' first
    bdf <- bdf %>%
      dplyr::select(all_of(c("id", names(bdf)[!names(bdf) %in% "id"])))
    
    return(bdf)
  } else {
    # Filter where the last column is larger than 0, select 'id' and the last column
    bdf <- bdf %>%
      filter(!!sym(last_col_name) > 0) %>%
      dplyr::select(all_of(c("id", last_col_name))) %>%
      rename(burning = !!sym(last_col_name))
    
    return(bdf)
  }
}

get_burners_optimized_dt <- function(time_horizon,
                                     ignition_cell, # This should now be the 'id_col' value of the ignition cell
                                     cell_size = 200, # in meters
                                     time_step = 15,    # in minutes
                                     scenario = "bau", # select scenario (bau, sal)
                                     full = FALSE, # return full data frame or last column
                                     neighbor_list_differentiated_preprocessed, # NOW PRE-PROCESSED
                                     grd, # Pass grd directly (should be a data.table with 'id_col' and keyed)
                                     con_scaler = 2000, # rescale unitless connectivity
                                     verbose = FALSE # Verbose output for debugging
) {
  # Input Validation for grd: Ensure it's a data.table and has 'id_col'
  if (!inherits(grd, "data.table")) {
    stop("Input 'grd' must be a data.table object.")
  }
  if (!"id_col" %in% names(grd)) {
    stop("Input 'grd' must contain an 'id_col' column (integer unique identifier).")
  }
  if (!key(grd) == "id_col") {
    setkey(grd, id_col) # Ensure key is set for faster joins
  }
  
  # Define variable names based on scenario
  speed_var <- paste0("ros_trg_", scenario)
  connectivity_var <- paste0("con_trg_", scenario)
  
  # Check if scenario-dependent columns exist early
  if (!speed_var %in% names(grd)) {
    stop(paste0("Speed variable column '", speed_var, "' not found in 'grd'."))
  }
  if (!connectivity_var %in% names(grd)) {
    stop(paste0("Connectivity variable column '", connectivity_var, "' not found in 'grd'."))
  }
  
  # Prepare empty data.table for simulation results (bdf)
  n_cells <- nrow(grd)
  bdf <- data.table(id_col = grd$id_col)
  for (col_idx in 0:time_horizon) {
    set(bdf, j = paste0("V", col_idx), value = 0)
  }
  setkey(bdf, id_col) # Key bdf for efficient lookups
  
  # Ensure ignition_cell is a valid id_col
  if (!ignition_cell %in% grd$id_col) {
    stop(paste0("Ignition cell '", ignition_cell, "' not found in 'grd$id_col'."))
  }
  
  # Set ignition cell
  set(bdf, i = which(bdf$id_col == ignition_cell), j = "V0", value = 1)
  
  # Pre-calculate speed and attraction values as new columns in grd
  # This avoids repeated column indexing inside the loop
  grd[, `_speed_values_` := get(speed_var)]
  grd[, `_attraction_values_` := get(connectivity_var) / con_scaler]
  
  # The neighbor_list_differentiated_preprocessed is ALREADY in the correct format
  # and named by id_col values, so no conversion logic needed here.
  # Let's rename it for clarity in the function.
  neighbor_map <- neighbor_list_differentiated_preprocessed
  
  # Loop for time steps
  for (t_idx in 1:time_horizon) {
    prev_t_col_name <- paste0("V", t_idx - 1)
    curr_t_col_name <- paste0("V", t_idx)
    
    # Identify currently burning cells (id_col values)
    burning_ids <- bdf[get(prev_t_col_name) >= 1, id_col]
    
    # Carry over burning status for already burning cells.
    if (length(burning_ids) > 0) {
      set(bdf, i = which(bdf$id_col %in% burning_ids), j = curr_t_col_name, value = 1)
    }
    
    # Initialize a temporary data.table for burn pressure for this time step
    temp_burn_pressure_dt <- data.table(id_col = grd$id_col, pressure = 0)
    setkey(temp_burn_pressure_dt, id_col)
    
    if (length(burning_ids) > 0) {
      # Data.table version for processing all burning cells' neighbors at once
      # Collect all potential neighbors (rook and queen) from all burning cells
      all_potential_rook_neighbors <- unlist(lapply(as.character(burning_ids), function(id) neighbor_map[[id]]$rook))
      all_potential_queen_neighbors <- unlist(lapply(as.character(burning_ids), function(id) neighbor_map[[id]]$queen))
      
      # Filter out already burning cells from all potential neighbors
      potential_rook_burners_ids <- setdiff(all_potential_rook_neighbors, burning_ids)
      potential_queen_burners_ids <- setdiff(all_potential_queen_neighbors, burning_ids)
      
      # Process Rook Neighbors (vectorized)
      if (length(potential_rook_burners_ids) > 0) {
        # Get relevant data for potential rook burners
        dt_rook_burners <- grd[.(potential_rook_burners_ids),
                               .(`_speed_values_`, `_attraction_values_`),
                               on = "id_col"]
        dt_rook_burners[, advancement_level := bdf[.(potential_rook_burners_ids), get(prev_t_col_name), on = "id_col"]]
        
        speed_per_time_step_adjusted <- (dt_rook_burners$`_speed_values_` / cell_size) * time_step
        burn_pressure_rook <- spread_fun(dt_rook_burners$advancement_level, speed_per_time_step_adjusted, dt_rook_burners$`_attraction_values_`)
        
        # Update temp_burn_pressure_dt
        for (i in seq_along(potential_rook_burners_ids)) {
          current_id <- potential_rook_burners_ids[i]
          current_pressure <- burn_pressure_rook[i]
          # Use direct update on temp_burn_pressure_dt
          current_idx_in_temp <- which(temp_burn_pressure_dt$id_col == current_id)
          if (length(current_idx_in_temp) > 0) {
            set(temp_burn_pressure_dt, i = current_idx_in_temp, j = "pressure",
                value = pmax(temp_burn_pressure_dt[current_idx_in_temp, pressure], current_pressure))
          }
        }
      }
      
      # Process Queen (Exclusive) Neighbors (vectorized)
      if (length(potential_queen_burners_ids) > 0) {
        dt_queen_burners <- grd[.(potential_queen_burners_ids),
                                .(`_speed_values_`, `_attraction_values_`),
                                on = "id_col"]
        dt_queen_burners[, advancement_level := bdf[.(potential_queen_burners_ids), get(prev_t_col_name), on = "id_col"]]
        
        distance_factor_queen <- sqrt(2)
        speed_per_time_step_adjusted <- (dt_queen_burners$`_speed_values_` / (cell_size * distance_factor_queen)) * time_step
        burn_pressure_queen <- spread_fun(dt_queen_burners$advancement_level, speed_per_time_step_adjusted, dt_queen_burners$`_attraction_values_`)
        
        # Update temp_burn_pressure_dt
        for (i in seq_along(potential_queen_burners_ids)) {
          current_id <- potential_queen_burners_ids[i]
          current_pressure <- burn_pressure_queen[i]
          current_idx_in_temp <- which(temp_burn_pressure_dt$id_col == current_id)
          if (length(current_idx_in_temp) > 0) {
            set(temp_burn_pressure_dt, i = current_idx_in_temp, j = "pressure",
                value = pmax(temp_burn_pressure_dt[current_idx_in_temp, pressure], current_pressure))
          }
        }
      }
    } # End if burning_ids > 0
    
    # Apply the accumulated burn pressure to the bdf for the current time step
    cells_to_update_pressure_dt <- temp_burn_pressure_dt[pressure > 0]
    if (nrow(cells_to_update_pressure_dt) > 0) {
      bdf[cells_to_update_pressure_dt,
          (curr_t_col_name) := pmax(get(curr_t_col_name), i.pressure),
          on = "id_col"]
    }
    
    # Ensure advancement level doesn't exceed 1 for the previous time step
    set(bdf, j = prev_t_col_name, value = pmin(1, bdf[[prev_t_col_name]]))
  } # End time horizon loop
  
  # Final cleanup of temporary columns in grd
  grd[, `_speed_values_` := NULL]
  grd[, `_attraction_values_` := NULL]
  
  # Post-processing (using dplyr for convenience, can be converted to data.table if strictly needed)
  # Convert 'bdf' back to a data.frame if 'full = TRUE' or if it's the expected output type.
  # Ensure 'id' column from original 'grd' is handled correctly.
  # Assuming 'grd_original' (from pre-processing) has the cell_id column that maps to original row names.
  # Or, if 'grd' already has a 'cell_id' column, use that.
  # For now, let's assume 'grd' has an 'original_row_name' column or similar if you need that specific output.
  # Otherwise, 'id_col' is the unique identifier.
  
  # Merge with original row names if needed for the final output format.
  # This requires 'grd_original' (or similar) to be available or to store the mapping.
  # For simplicity, let's assume 'grd' carries the `cell_id` or similar original identifier.
  # For this example, let's just make 'id_col' the 'id' in the final output.
  final_bdf <- copy(bdf) # Create a copy to avoid modifying bdf by reference if it's used elsewhere.
  
  last_col_name <- paste0("V", time_horizon)
  
  if (full) {
    # If the original row names are important for output, you need to join them back.
    # Assuming grd_original from the preprocessing step could be used, or you store the mapping.
    # For simplicity, if 'grd' always has the mapping of 'id_col' to 'original_name_string':
    # final_bdf <- final_bdf[grd[, list(id_col, original_name_string)], on = "id_col"]
    # setnames(final_bdf, "original_name_string", "id")
    
    # If `id_col` is the only ID needed:
    final_bdf <- final_bdf[get(last_col_name) > 0]
    setnames(final_bdf, "id_col", "id")
    # No `dplyr::select` needed if all columns are relevant and id is first
    return(as.data.frame(final_bdf)) # Return as data.frame as per original function
  } else {
    final_bdf <- final_bdf[get(last_col_name) > 0, .(id_col, (!!sym(last_col_name)))]
    setnames(final_bdf, c("id_col", last_col_name), c("id", "burning"))
    return(as.data.frame(final_bdf)) # Return as data.frame as per original function
  }
}


get_burners_optimized_dt_final <- function(time_horizon,
                                           ignition_cell, # This should be the 'id_col' value
                                           cell_size = 200, # in meters
                                           time_step = 15,    # in minutes
                                           scenario = "bau", # select scenario (bau, sal)
                                           full = FALSE, # return full data frame or last column
                                           neighbor_connections_dt, # NEW: Pre-processed data.table of connections
                                           grd, # Pass grd directly (must be a data.table with 'id_col' and keyed)
                                           con_scaler = 2000, # rescale unitless connectivity
                                           verbose = FALSE # Verbose output for debugging
) {
  # --- Input Validation ---
  if (!inherits(grd, "data.table") || !"id_col" %in% names(grd) || key(grd)[1] != "id_col") {
    stop("Input 'grd' must be a data.table object with an 'id_col' key.")
  }
  if (!inherits(neighbor_connections_dt, "data.table") || !all(c("origin_id", "neighbor_id", "type", "distance_factor") %in% names(neighbor_connections_dt))) {
    stop("Input 'neighbor_connections_dt' must be a data.table with 'origin_id', 'neighbor_id', 'type', 'distance_factor' columns.")
  }
  
  speed_var <- paste0("ros_trg_", scenario)
  connectivity_var <- paste0("con_trg_", scenario)
  
  if (!speed_var %in% names(grd) || !connectivity_var %in% names(grd)) {
    stop(paste0("Scenario-dependent columns (", speed_var, ", ", connectivity_var, ") not found in 'grd'."))
  }
  
  # --- Prepare grd for calculation (add temporary columns) ---
  grd_local <- copy(grd) # Create a local copy to modify, avoid modifying original `grd` passed by reference
  
  # --- FIX: Use simple, syntactic names for temporary columns ---
  grd_local[, temp_speed := get(speed_var)]
  grd_local[, temp_attraction := get(connectivity_var) / con_scaler]
  
  # --- Initialize bdf (burn dataframe) ---
  n_cells <- nrow(grd_local)
  bdf <- data.table(id_col = grd_local$id_col)
  for (col_idx in 0:time_horizon) {
    set(bdf, j = paste0("V", col_idx), value = 0)
  }
  setkey(bdf, id_col)
  
  # Set ignition cell
  if (!ignition_cell %in% grd_local$id_col) {
    stop(paste0("Ignition cell '", ignition_cell, "' not found in 'grd$id_col'."))
  }
  set(bdf, i = which(bdf$id_col == ignition_cell), j = "V0", value = 1)
  
  # --- Main Simulation Loop ---
  for (t_idx in 1:time_horizon) {
    prev_t_col_name <- paste0("V", t_idx - 1)
    curr_t_col_name <- paste0("V", t_idx)
    
    # 1. Identify currently burning cells (ids)
    burning_ids <- bdf[get(prev_t_col_name) >= 1, id_col]
    
    # 2. Carry over burning status to current time step
    if (length(burning_ids) > 0) {
      set(bdf, i = which(bdf$id_col %in% burning_ids), j = curr_t_col_name, value = 1)
    }
    
    # If no cells are burning, nothing more to do this step
    if (length(burning_ids) == 0) {
      set(bdf, j = prev_t_col_name, value = pmin(1, bdf[[prev_t_col_name]])) # Ensure previous is capped
      next # Skip to next time step
    }
    
    # 3. Get all potential neighbor connections from currently burning cells
    active_connections <- neighbor_connections_dt[.(burning_ids), on = "origin_id", nomatch = 0]
    
    # 4. Filter out neighbors that are already burning
    potential_new_burners_connections <- active_connections[! (neighbor_id %in% burning_ids)]
    
    if (nrow(potential_new_burners_connections) > 0) {
      # 5. Get attributes for these potential new burners
      #    Perform joins sequentially to avoid ambiguity and correctly reference columns
      #    First, join `grd_local` with `potential_new_burners_connections`
      temp_data <- grd_local[potential_new_burners_connections,
                             on = .(id_col = neighbor_id),
                             nomatch = 0]
      
      #    Then, join the result with `bdf` to get advancement_level
      potential_new_burners_data <- bdf[temp_data,
                                        on = .(id_col)][
                                          , .(origin_id = i.origin_id, # `i.` refers to columns from `temp_data` (right side)
                                              neighbor_id = id_col, # `id_col` is the matched column from `bdf` (left side)
                                              type = i.type,
                                              distance_factor = i.distance_factor,
                                              advancement_level = get(prev_t_col_name),
                                              speed = i.temp_speed, # --- FIX: Use new clean name ---
                                              attraction = i.temp_attraction # --- FIX: Use new clean name ---
                                          )]
      
      # 6. Calculate burn pressure for all potential new burners
      potential_new_burners_data[, speed_per_time_step_adjusted := (speed / (cell_size * distance_factor)) * time_step]
      potential_new_burners_data[, burn_pressure_raw := spread_fun(advancement_level, speed_per_time_step_adjusted, attraction)]
      
      # 7. Aggregate burn pressure by neighbor_id
      aggregated_burn_pressure <- potential_new_burners_data[, max(burn_pressure_raw), by = neighbor_id]
      setnames(aggregated_burn_pressure, "V1", "max_burn_pressure")
      
      # 8. Apply the aggregated burn pressure to bdf for the current time step
      bdf[aggregated_burn_pressure,
          (curr_t_col_name) := pmax(get(curr_t_col_name), i.max_burn_pressure),
          on = .(id_col = neighbor_id)]
    }
    
    # Ensure advancement level doesn't exceed 1 for the previous time step
    set(bdf, j = prev_t_col_name, value = pmin(1, bdf[[prev_t_col_name]]))
  } # End time horizon loop
  
  # --- Cleanup Temporary Columns in grd_local ---
  grd_local[, temp_speed := NULL]
  grd_local[, temp_attraction := NULL]
  
  # --- Final Post-processing ---
  final_bdf <- copy(bdf)
  last_col_name <- paste0("V", time_horizon)
  
  if (full) {
    final_bdf <- final_bdf[get(last_col_name) > 0]
    setnames(final_bdf, "id_col", "id")
    return(as.data.frame(final_bdf))
  } else {
    final_bdf <- final_bdf[get(last_col_name) > 0, .(id_col, (!!sym(last_col_name)))]
    setnames(final_bdf, c("id_col", last_col_name), c("id", "burning"))
    return(as.data.frame(final_bdf))
  }
}


#### Loop ####

# Compare performance of the two functions
# load data

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")


# calculate 850mb limit:
# 850*1024^2 = 891289600
options(future.globals.maxSize= 891289600, 
        future.seed = TRUE)

handlers(handler_progress(format="[:bar] :percent :eta :message"))


library(future)
library(future.apply)
library(progressr)
library(dplyr)

# ----------------------------------
# Load data
# ----------------------------------
target_df <- readRDS("./data/rdat/target_df.rds")
load("./data/rdat/neighbor_idx_differentiated.RData")

# Add scenarios (as before)
for (y in 2040:2050) {
  target_df[paste0("ros_trg_baseline_y", y)] <- target_df[paste0("bau_ROS_", y)]
  target_df[paste0("con_trg_baseline_y", y)] <- target_df[paste0("bmr_FireConn", y, "b_raster")]
  target_df[paste0("ros_trg_full_y", y)] <- target_df[paste0("sal_ROS_", y)]
  target_df[paste0("con_trg_full_y", y)] <- target_df[paste0("bmr_FireConn", y, "s_raster")]
}

# ----------------------------------
# Parameters
# ----------------------------------
set.seed(123)
rdm_ign_cells <- sample(
  x = row.names(target_df),
  size = 1000,
  replace = FALSE,
  prob = target_df$ignition_probability_surface
)

batch_size <- 1
mc_time_step <- 5
mc_time_horizon <- 180

all_scenarios <- names(target_df) %>%
  subset(grepl("ros_trg_|con_trg_", .)) %>%
  gsub("ros_trg_|con_trg_", "", .) %>%
  unique()

aggregate_cells <- FALSE

# ----------------------------------
# Prepare batching + output
# ----------------------------------
out_dir <- "./data/rdat/mc_fire_batches/"
dir.create(out_dir, showWarnings = FALSE)

ignition_batches <- split(rdm_ign_cells, ceiling(seq_along(rdm_ign_cells) / batch_size))

# ----------------------------------
# Function for one batch
# ----------------------------------


library(future)
library(future.apply)
library(progressr)
library(dplyr)

# --- Setup cluster and export big globals once ---
cl <- makeClusterPSOCK(20)  # adjust worker count
clusterExport(cl, c("target_df", "neighbor_idx_differentiated",
                    "get_burners_optimized",
                    "mc_time_horizon", "mc_time_step", "aggregate_cells"))
plan(cluster, workers = cl)

# --- Define batch function (sequential inside each batch) ---
process_batch <- function(batch_id, ign_cells, scenarios) {
  res_batch <- list()
  
  for (s_name in scenarios) {
    res_scenario <- lapply(ign_cells, function(x) {
      get_burners_optimized(
        x,
        time_horizon = mc_time_horizon,
        time_step = mc_time_step,
        full = TRUE,
        scenario = s_name,
        grd = target_df,
        con_scaler = 2000,
        neighbor_list_differentiated = neighbor_idx_differentiated
      )
    })
    
    if (aggregate_cells) {
      res_scenario <- res_scenario %>%
        lapply(., function(x) apply(x, 2, sum)) %>%
        bind_rows() %>%
        as.data.frame() %>%
        mutate(scenario = s_name)
    }
    
    res_batch[[s_name]] <- res_scenario
  }
  
  saveRDS(res_batch, file = file.path(out_dir, paste0("batch_", batch_id, ".rds")))
  return(TRUE)
}

# --- Resume logic: find todo batches ---
done_batches <- list.files(out_dir, pattern = "batch_.*\\.rds") %>%
  gsub("batch_|\\.rds", "", .) %>%
  as.integer()
todo_batches <- setdiff(seq_along(ignition_batches), done_batches)

# --- Run batches in parallel with progress bar ---
handlers(global = TRUE)
with_progress({
  p <- progressor(along = todo_batches)
  
  future_lapply(todo_batches, function(b) {
    ign_cells <- ignition_batches[[b]]
    process_batch(b, ign_cells, all_scenarios)
    p(sprintf("Completed batch %d", b))
  })
})

# --- cleanup cluster ---
plan(sequential)
parallel::stopCluster(cl)


# --- Final combine step ---

library(dplyr)
library(tidyr)
library(purrr)

# --- Land cover data for merging ---
lc_dat <- target_df %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  select(all_of(c("cell_id", "lc_agriculture", "lc_wildland","lc_artificial")))

library(data.table)

# --- Land cover table (data.table, not tibble) ---
lc_dat <- as.data.table(target_df)
lc_dat[, cell_id := .I]   # row index as cell id
lc_dat <- lc_dat[, .(cell_id, lc_agriculture, lc_wildland, lc_artificial)]
lc_dat$cell_id <- paste0("cell_", lc_dat$cell_id) 

# --- Process each batch in a streaming way ---
batch_files <- list.files(out_dir, pattern = "batch_.*\\.rds", full.names = TRUE)


for (i in seq_along(batch_files)) {
  message("Processing ", batch_files[i])
  
  # Load batch
  batch_res <- readRDS(batch_files[i])
  
  print("loaded.")
  
  # Empty list to collect aggregated scenario results (only small objects)
  batch_agg_list <- vector("list", length(batch_res))
  
  # Loop over scenarios in this batch
  for (s in seq_along(batch_res)) {
    scenario_name <- names(batch_res)[s]
    scenario_dfs <- batch_res[[s]]
    
    print(scenario_name)
    
    # List to hold per-ignition aggregated data
    scen_agg_list <- vector("list", length(scenario_dfs))
    
    for (j in seq_along(scenario_dfs)) {
      df <- as.data.table(scenario_dfs[[j]])
      
      # Identify ignition id
      ignition_id <- df$id[df$V1 == 1]
      
      # Wide -> long
      df_long <- melt(df, id.vars = "id",
                      variable.name = "time",
                      value.name = "burned")
      
      # Clean time variable
      df_long[, time := as.numeric(gsub("V", "", time)) * mc_time_step]
      
      # Merge with land cover
      df_long <- merge(df_long, lc_dat, by.x = "id", by.y = "cell_id", all.x = TRUE)
      
      # Aggregate immediately
      df_agg <- df_long[, .(
        burned_land_ha = sum(burned * 4, na.rm = TRUE),
        burned_agri_ha = sum(burned * lc_agriculture * 4, na.rm = TRUE),
        burned_wild_ha = sum(burned * lc_wildland * 4, na.rm = TRUE),
        burned_urbn_ha = sum(burned * lc_artificial * 4, na.rm = TRUE)
      ), by = .(time)]
      
      # Add scenario + ignition info
      df_agg[, scenario := scenario_name]
      df_agg[, ignition_id := ignition_id]
      
      # Store aggregated result
      scen_agg_list[[j]] <- df_agg
    }
    
    # Combine all ignitions for this scenario
    batch_agg_list[[s]] <- rbindlist(scen_agg_list)
  }
  
  # Combine all scenarios in this batch
  batch_agg <- rbindlist(batch_agg_list)
  
  # Save aggregated batch result
  saveRDS(batch_agg, file = file.path(out_dir, paste0("batch_lc_agg_", i, ".rds")))
  rm(batch_res, batch_agg_list, batch_agg); gc()
}

# --- Final combine step ---
agg_files <- list.files(out_dir, pattern = "batch_lc_agg_.*\\.rds", full.names = TRUE)
all_agg_results <- rbindlist(lapply(agg_files, readRDS))
saveRDS(all_agg_results, "./data/rdat/mc_res_agg_all.rds")

message("✅ Finished: combined aggregated results saved to mc_res_agg_all.rds")




















# 23 seconds

### Now with data.table 

# 2. Convert to data.table and add 'id_col'
target_dt <- as.data.table(target_df)
target_dt[, id_col := .I] # CRUCIAL: Creates a unique integer ID for each row

# 3. Set the key for target_dt
setkey(target_dt, id_col) # CRUCIAL: Sets 'id_col' as the key, enabling fast lookups

# 4. Create the name_to_id_map
name_to_id_map <- setNames(target_dt$id_col, rownames(target_df))

# 5. Pre-process neighbor_list_differentiated using the map
raw_neighbor_list_differentiated <- neighbor_idx_differentiated

# 3. Pre-process neighbor_list_differentiated into a data.table format
#    This is the NEW and CRITICAL pre-processing step.
neighbor_connections_dt <- data.table(
  origin_id = integer(),
  neighbor_id = integer(),
  type = character(), # "rook" or "queen"
  distance_factor = numeric()
)

for (origin_cell_name in names(raw_neighbor_list_differentiated)) {
  origin_id_val <- name_to_id_map[[origin_cell_name]]
  neighbors <- raw_neighbor_list_differentiated[[origin_cell_name]]
  
  # Process rook neighbors
  if (length(neighbors$rook) > 0) {
    rook_neighbor_ids <- name_to_id_map[neighbors$rook]
    new_rows <- data.table(
      origin_id = origin_id_val,
      neighbor_id = rook_neighbor_ids,
      type = "rook",
      distance_factor = 1 # Rook distance factor
    )
    neighbor_connections_dt <- rbindlist(list(neighbor_connections_dt, new_rows))
  }
  
  # Process queen neighbors
  if (length(neighbors$queen) > 0) {
    queen_neighbor_ids <- name_to_id_map[neighbors$queen]
    new_rows <- data.table(
      origin_id = origin_id_val,
      neighbor_id = queen_neighbor_ids,
      type = "queen",
      distance_factor = sqrt(2) # Queen distance factor
    )
    neighbor_connections_dt <- rbindlist(list(neighbor_connections_dt, new_rows))
  }
}

# Ensure neighbor_connections_dt is keyed for efficient joins
setkey(neighbor_connections_dt, origin_id)

# IMPORTANT: Check for NAs after mapping to ensure all names were found
if (anyNA(neighbor_connections_dt$neighbor_id) || anyNA(neighbor_connections_dt$origin_id)) {
  stop("Pre-processing error: Some neighbor names could not be mapped to 'id_col' in 'grd'.")
}


all_scenario_results_dt <- list() # To store results from each scenario

raw_row_names <- as.integer(gsub("cell_", "", rdm_ign_cells))

# setup parallel processing
plan(multisession, workers=5)

# start time
start_time <- Sys.time()

# Use future_lapply for the outer loop over scenarios
all_scenario_results_dt <- future_lapply(scenarios, function(s_name) {
  p_scenarios(message = paste("Processing scenario:", s_name)) # Progress for scenarios
  
  # Inner parallel loop (or sequential, depending on your plan)
  # The inner lapply remains the same as it's already vectorized/parallelizable
  mc_res_current_scenario <- lapply(raw_row_names, function(x) {
    get_burners_optimized_dt_final(x,
                          time_horizon = mc_time_horizon,
                          time_step = mc_time_step,
                          full = TRUE,
                          scenario = s_name, # Use the current scenario name
                          grd = target_dt,
                          con_scaler = 2000,
                          neighbor_connections_dt = neighbor_connections_dt)
  })
  
  if (aggregate_cells) {
    # Aggregate results for the current scenario
    agg_current_scenario <- mc_res_current_scenario %>%
      lapply(., function(x) apply(x, 2, sum)) %>%
      bind_rows() %>%
      as.data.frame() %>%
      mutate(scenario = s_name) # Assign the correct scenario name
    
    return(agg_current_scenario)
  } else {
    return(mc_res_current_scenario)
  }
})

# If aggregate_cells is TRUE, all_scenario_results will be a list of data.frames.
# If aggregate_cells is FALSE, all_scenario_results will be a list of lists.
# You might want to name the elements of all_scenario_results based on scenario names
# after the future_lapply call, as future_lapply doesn't preserve names by default.
names(all_scenario_results_dt) <- scenarios

# End time
end_time <- Sys.time()
# Print the time taken for the entire operation
time_taken <- end_time - start_time
message(paste("Total time taken for all scenarios (DT):", time_taken))






t1 <- get_burners_optimized_dt_final(40843,
                                     time_horizon = mc_time_horizon,
                                     time_step = mc_time_step,
                                     full = TRUE,
                                     scenario = scenarios[1], # Use the current scenario name
                                     grd = target_dt,
                                     con_scaler = 2000,
                                     neighbor_connections_dt = neighbor_connections_dt)

