##### Functions #####


library(sf)
library(spdep)
library(dplyr)

spread_fun = function (advancement_level, speed, attraction){
  advancement_level + (speed * attraction)
}

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
  speed_var <- paste0("biomass_", scenario)
  connectivity_var <- paste0("connectivity_", scenario)
  
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
  
  # Convert neighbor_list_differentiated to use indices directly for faster access
  neighbor_idx_differentiated <- lapply(neighbor_list_differentiated, function(cell_neighbors) {
    list(
      rook = which(row.names(grd_attr) %in% cell_neighbors$rook),
      queen = which(row.names(grd_attr) %in% cell_neighbors$queen)
    )
  })
  names(neighbor_idx_differentiated) <- row.names(grd_attr) # Maintain names
  
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


# get_burners <- function(time_horizon,
#                         ignition_cell,
#                         cell_size = 200, # in meters
#                         time_step = 15,    # in minutes
#                         scenario = "bau", # select scenario (bau, sal)
#                         full = FALSE # return full data frame or last column
# ) {
#   
#   # define variable names based on scenario
#   speed_var <- paste0("biomass_", scenario)
#   connectivity_var <- paste0("connectivity_", scenario)
#   
#   grd_geo <- grd
#   grd <- grd %>% st_set_geometry(NULL)
#   
#   # prepare empty grid for simulation results
#   n_cells <- length(unique(grd$id)) # Number of cells
#   bdf <- matrix(0, nrow = n_cells, ncol = time_horizon + 1) # Initialize as a matrix
#   
#   # Find the row index of the ignition cell
#   ignition_idx <- which(grd$id == ignition_cell)  # Assuming grd has the cell IDs
#   bdf[ignition_idx, 1] <- 1  # Set ignition cell
#   
#   #cat("Ignition started \n")
#   
#   for (t in 2:(time_horizon + 1)) {
#     #cat("Time step:", t - 1, "\n")
#     
#     prev_t <- t - 1
#     
#     burning_idx <- which(bdf[, prev_t] >= 1)
#     
#     bdf[burning_idx, t] <- 1
#     
#     for (origin_cell in burning_idx) {
#       #cat("Origin cell:", origin_cell, "\n")
#       
#       potential_burners <- neighbor_idx[[origin_cell]]
#       potential_burners <- potential_burners[!(potential_burners %in% burning_idx)]
#       
#       #cat("Potential burners:", potential_burners, "\n")
#       
#       for (target_cell in potential_burners) {
#         
#         #target_idx <- which(grd$id == target_cell) # Use grd to find index
#         advancement_level <- bdf[target_cell, prev_t]
#         speed <- grd[,speed_var][target_cell]
#         speed_per_time_step <- speed / cell_size * time_step
#         attraction <- grd[,connectivity_var][target_cell]
#         
#         burn_pressure <- spread_fun(advancement_level, 
#                                     speed_per_time_step, 
#                                     attraction)
#         
#         #cat("Burn pressure:", burn_pressure, "\n")
#         
#         bdf[target_cell, t] <- burn_pressure
#         bdf[, prev_t] <- pmin(1, bdf[, prev_t]) #  min
#       }
#     }
#   }
#   
#   #cat("Simulation finished \n")
#   
#   # Convert back to data frame for merging
#   bdf <- as.data.frame(bdf)
#   
#   last_col_name <- tail(names(bdf), 1)
#   
#   if (full) {
#     bdf <- bdf %>% 
#       filter(!!sym(last_col_name) > 0)
#     return(bdf)
#   } else {
#     
#     # merge with spatial data, filter where last column is larger 0
#     bdf <- bind_cols(grd_geo, bdf) %>% 
#       filter(!!sym(last_col_name) > 0) %>%
#       dplyr::select(all_of(c("id", last_col_name))) %>%
#       rename(burning = !!sym(last_col_name))
#     
#     return(bdf)
#   }
# }

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

