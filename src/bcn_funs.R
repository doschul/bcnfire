##### Functions #####


library(sf)
library(spdep)
library(dplyr)


# Make sure you have dplyr loaded if using `%>%`, `filter`, `select`, `rename`, `bind_cols`, `sym`, `all_of`.
# library(dplyr) # Uncomment if not already loaded

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




one_cell_scenario_path = function(ignition_cell, 
                                  scenarios = c("bau", "sal"),
                                  time_horizon = 6,
                                  time_step = 10
                                  ){
  
  # prepare viz_list with scenario results
  viz_list <- lapply(scenarios, function(x){
    get_burners(time_horizon = time_horizon,
                ignition_cell = ignition_cell,
                time_step = time_step,    # in minutes
                scenario = x, # select scenario (bau, sal)
                full = TRUE)})
    
    # add name column to each list df
    for (s in 1:length(scenarios)) {
      
      viz_list[[s]] <- apply(viz_list[[s]], 2, sum) %>% t() %>%
        as.data.frame(.) %>%
        mutate(scenario = scenarios[s]) %>%
        pivot_longer(-scenario, names_to = "time", values_to = "burned") %>%
        mutate(time = as.numeric(gsub("V", "", time)) * time_step)
    }
    
    viz_df <- do.call(bind_rows, viz_list)
    
    p <- ggplot(viz_df) +
      geom_line(aes(x = time, y = burned, color = scenario)) +
      labs(x = "Time (minutes)", y = "Burned cells") +
      theme_minimal()
    
    return(p)
}


BCN_target_df = function(
    data,
    original_cols, # Must be length 2: [ROS_col, Con_col]
    mask_col,
    replacement_cols, # Must be length 2: [ROS_repl_col, Con_repl_col]
    target_col = "random", # Optional: if not provided, random assignment is used
    threshold = NULL,
    budget = NULL,
    relative = TRUE,
    target_decreasing = FALSE,
    year = NULL, # This is the year parameter for the targeting scenario
    verbose = FALSE
) {
  # --- Input Validation ---
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data.frame object.")
  }
  
  if (!is.character(original_cols) || length(original_cols) != 2) {
    stop("Argument 'original_cols' must be a character vector of length 2 (ROS, Connectivity).")
  }
  if (!is.character(replacement_cols) || length(replacement_cols) != 2) {
    stop("Argument 'replacement_cols' must be a character vector of length 2 (ROS, Connectivity).")
  }
  # Lengths are checked above, so this specific check for equality of lengths is redundant but harmless
  if (length(original_cols) != length(replacement_cols)) {
    stop("The length of 'original_cols' must match the length of 'replacement_cols'.")
  }
  
  mandatory_single_col_args = c("mask_col")
  for (arg in mandatory_single_col_args) {
    if (is.null(get(arg))) {
      stop(paste("Argument '", arg, "' must be provided.", sep = ""))
    }
  }
  
  all_specified_cols = c(original_cols, mask_col, replacement_cols)
  if (target_col != "random") {
    all_specified_cols = c(all_specified_cols, target_col)
  }
  
  if (!all(all_specified_cols %in% names(data))) {
    missing_cols = setdiff(all_specified_cols, names(data))
    stop(paste("The following specified columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (is.null(threshold) && is.null(budget)) {
    stop("Either 'threshold' or 'budget' must be provided.")
  }
  
  if (!is.logical(relative)) {
    stop("'relative' must be a boolean value (TRUE or FALSE).")
  }
  
  if (!all(data[[mask_col]] %in% c(0, 1, NA))) {
    stop(paste("Column '", mask_col, "' must contain only 0s, 1s, or NAs.", sep = ""))
  }
  
  # --- Calculate Intervention Cells ---
  original_intervention_indices = which(data[[mask_col]] == 1 & !is.na(data[[mask_col]]))
  intervention_cells = length(original_intervention_indices)
  
  # --- Generate Column Name Parts for Suffixes ---
  param_suffix <- if (relative) {
    if (!is.null(threshold)) paste0("_t", as.character(threshold * 100), "_y", year) else ""
  } else {
    if (!is.null(budget)) paste0("_b", as.character(budget), "_y", year) else ""
  }
  
  strategy_name_part <- if (target_col != "random") paste0("_", target_col) else "_random"
  
  # New mask column name
  new_mask_col_name = paste0("new_mask", strategy_name_part, param_suffix)
  
  # --- If no cells are marked for intervention, return the original data with a message ---
  if (intervention_cells == 0) {
    if (verbose) {
      message("No cells found with 'intervention_mask' = 1. Returning original data.")
    }
    data[[new_mask_col_name]] = 0
    
    # Initialize new ROS and Connectivity columns
    ros_trg_col_name = paste0("ros_trg", strategy_name_part, param_suffix)
    con_trg_col_name = paste0("con_trg", strategy_name_part, param_suffix)
    
    data[[ros_trg_col_name]] = data[[original_cols[1]]] # First original_col is ROS
    data[[con_trg_col_name]] = data[[original_cols[2]]] # Second original_col is Connectivity
    
    return(data)
  }
  
  # --- Determine Target Cell Count ---
  if (relative) {
    if (is.null(threshold) || (threshold <= 0 || threshold > 1)) {
      stop("If 'relative' is TRUE, 'threshold' (between 0 and 1) must be provided.")
    }
    target_cells = round(intervention_cells * threshold)
  } else {
    if (is.null(budget)) {
      stop("If 'relative' is FALSE, 'budget' must be provided.")
    }
    target_cells = budget
  }
  
  if (target_cells > intervention_cells) {
    warning(paste("Target cells (", target_cells, ") exceed total intervention cells (", intervention_cells, "). Setting target_cells to intervention_cells."))
    target_cells = intervention_cells
  }
  if (target_cells < 0) {
    warning(paste("Target cells (", target_cells, ") is negative. Setting target_cells to 0."))
    target_cells = 0
  }
  
  # --- Initialize New Columns ---
  data[[new_mask_col_name]] = 0
  
  # Define the new column names directly
  ros_trg_col_name = paste0("ros_trg", strategy_name_part, param_suffix)
  con_trg_col_name = paste0("con_trg", strategy_name_part, param_suffix)
  
  # Initialize new ROS and Connectivity columns with their original values
  data[[ros_trg_col_name]] = data[[original_cols[1]]] # First original_col is ROS
  data[[con_trg_col_name]] = data[[original_cols[2]]] # Second original_col is Connectivity
  
  # --- Select Target Cells for Modification ---
  target_cells_global_indices = integer(0)
  
  has_valid_target_surface = (target_col != "random") &&
    any(!is.na(data[[target_col]][original_intervention_indices]))
  
  if (!has_valid_target_surface) {
    if (target_cells > 0) {
      target_cells_global_indices = sample(
        original_intervention_indices,
        size = target_cells,
        replace = FALSE
      )
    }
  } else {
    temp_selection_data = data.frame(
      original_index = original_intervention_indices,
      target_value = data[[target_col]][original_intervention_indices]
    )
    
    if (nrow(temp_selection_data) == 0) {
      if (verbose) {
        message("No valid 'target_col' values within intervention mask for ordering. Randomly selecting target cells.")
      }
      if (target_cells > 0) {
        target_cells_global_indices = sample(
          original_intervention_indices,
          size = target_cells,
          replace = FALSE
        )
      }
    } else {
      ordered_temp_data = temp_selection_data[order(temp_selection_data$target_value, decreasing = target_decreasing), ]
      
      if (target_cells > 0) {
        target_cells_global_indices = ordered_temp_data$original_index[1:min(target_cells, nrow(ordered_temp_data))]
      }
    }
  }
  
  # --- Apply Replacement Values and Update New Mask ---
  if (length(target_cells_global_indices) > 0) {
    # Apply ROS replacement
    data[[ros_trg_col_name]][target_cells_global_indices] = data[[replacement_cols[1]]][target_cells_global_indices]
    
    # Apply Connectivity replacement
    data[[con_trg_col_name]][target_cells_global_indices] = data[[replacement_cols[2]]][target_cells_global_indices]
    
    data[[new_mask_col_name]][target_cells_global_indices] = 1
  }
  
  # --- Verbose Output ---
  if (verbose) {
    message(paste("Number of target cells selected for modification:", length(target_cells_global_indices)))
    message(paste("Total number of intervention-eligible cells considered:", intervention_cells))
    message(paste("Calculated target cells threshold/budget:", target_cells))
  }
  
  return(data)
}
