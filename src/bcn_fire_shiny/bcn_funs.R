##### Functions #####
spread_fun = function (advancement_level, speed, attraction){
  advancement_level + (speed * attraction)
}


get_burners <- function(time_horizon,
                        ignition_cell,
                        cell_size = 200, # in meters
                        time_step = 15,    # in minutes
                        scenario = "bau", # select scenario (bau, sal)
                        full = FALSE, # return full data frame or last column
                        grd = grd,
                        neighbor_idx = neighbor_idx
) {
  
  # define variable names based on scenario
  speed_var <- paste0("biomass_", scenario)
  connectivity_var <- paste0("connectivity_", scenario)
  
  grd_geo <- grd
  grd <- grd %>% st_set_geometry(NULL)
  
  # prepare empty grid for simulation results
  n_cells <- length(unique(grd$id)) # Number of cells
  bdf <- matrix(0, nrow = n_cells, ncol = time_horizon + 1) # Initialize as a matrix
  
  # Find the row index of the ignition cell
  ignition_idx <- which(grd$id == ignition_cell)  # Assuming grd has the cell IDs
  bdf[ignition_idx, 1] <- 1  # Set ignition cell
  
  #cat("Ignition started \n")
  
  for (t in 2:(time_horizon + 1)) {
    #cat("Time step:", t - 1, "\n")
    
    prev_t <- t - 1
    
    burning_idx <- which(bdf[, prev_t] >= 1)
    
    bdf[burning_idx, t] <- 1
    
    for (origin_cell in burning_idx) {
      #cat("Origin cell:", origin_cell, "\n")
      
      potential_burners <- neighbor_idx[[origin_cell]]
      potential_burners <- potential_burners[!(potential_burners %in% burning_idx)]
      
      #cat("Potential burners:", potential_burners, "\n")
      
      for (target_cell in potential_burners) {
        
        #target_idx <- which(grd$id == target_cell) # Use grd to find index
        advancement_level <- bdf[target_cell, prev_t]
        speed <- grd[,speed_var][target_cell]
        speed_per_time_step <- speed / cell_size * time_step
        attraction <- grd[,connectivity_var][target_cell]
        
        burn_pressure <- spread_fun(advancement_level, 
                                    speed_per_time_step, 
                                    attraction)
        
        #cat("Burn pressure:", burn_pressure, "\n")
        
        bdf[target_cell, t] <- burn_pressure
        bdf[, prev_t] <- pmin(1, bdf[, prev_t]) #  min
      }
    }
  }
  
  #cat("Simulation finished \n")
  
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

