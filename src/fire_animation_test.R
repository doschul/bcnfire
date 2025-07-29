library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(spdep) # for poly2nb

# test new burner funcction
# Sample grid creation (as before)
grd <- st_make_grid(st_bbox(c(xmin=0, ymin=0, xmax=100, ymax=100)), n=100, square=TRUE) %>%
  st_sf() %>%
  mutate(id = 1:nrow(.)) # Add an ID column

row.names(grd) <- paste0("cell_", 1:nrow(grd))
#row.names(grd) <- paste0(1:nrow(grd))

# Add dummy biomass and connectivity data
grd$biomass_bau <- runif(nrow(grd), 0.1, 1)
grd$connectivity_bau <- runif(nrow(grd), 0.1, 1)

# Generate differentiated neighbor list (as in previous steps)
nb_queen_all <- poly2nb(grd, queen = TRUE)
nb_rook <- poly2nb(grd, queen = FALSE)

neighbor_list_differentiated <- vector("list", length = nrow(grd))
names(neighbor_list_differentiated) <- row.names(grd)

for (i in seq_along(nb_queen_all)) {
  current_cell_name <- row.names(grd)[i]
  rook_neighbors_indices <- nb_rook[[i]]
  rook_neighbors_names <- if (length(rook_neighbors_indices) > 0) {
    row.names(grd)[rook_neighbors_indices]
  } else {
    character(0)
  }
  all_queen_neighbors_indices <- nb_queen_all[[i]]
  exclusive_queen_neighbors_indices <- setdiff(all_queen_neighbors_indices, rook_neighbors_indices)
  exclusive_queen_neighbors_names <- if (length(exclusive_queen_neighbors_indices) > 0) {
    row.names(grd)[exclusive_queen_neighbors_indices]
  } else {
    character(0)
  }
  neighbor_list_differentiated[[current_cell_name]] <- list(
    rook = rook_neighbors_names,
    queen = exclusive_queen_neighbors_names
  )
}

# Run the optimized function
system.time({
  result_optimized <- get_burners_optimized(
    time_horizon = 20,
    ignition_cell = "cell_250", # Example ignition cell
    cell_size = 200,
    time_step = 5,
    scenario = "bau",
    full = TRUE,
    neighbor_list_differentiated = neighbor_list_differentiated,
    grd = grd # Pass grd to the function
  )
})

# Convert to long format and merge with geometry
res_long <- result_optimized %>%
  mutate(id = gsub("cell_", "", row.names(.))) %>% # Use . to refer to the current dataframe
  pivot_longer(cols = starts_with("V"),
               names_to = "time_step",
               values_to = "burning") %>%
  mutate(time_step = as.numeric(gsub("V", "", time_step))) %>%
  dplyr::select(id, time_step, burning) %>% # Use dplyr::select to avoid conflicts
  mutate(id = as.integer(id)) %>%
  # Merge with grd geometry
  left_join(grd, by = "id") %>%
  # Filter out cells that never burn (burning == 0 for all time steps)
  # This makes the animation cleaner by only showing relevant cells
  group_by(id) %>%
  filter(any(burning > 0)) %>%
  ungroup() |>
  st_as_sf() # Convert to sf object for plotting

# --- 4. Create the animation using gganimate ---


###########
library(tmap)
library(sf)
library(dplyr)
library(tidyr)


# Set tmap to plot mode (interactive or plot)
tmap_mode("plot") # Use "view" for interactive mapview-like output in RStudio viewer

# Define ignition_cell and time_step from your simulation parameters for labs
# (Assuming these are available from the broader script scope or you can manually define them)
# For example:
ignition_cell <- "cell_250"
cell_size <- 200
time_step <- 15

res_long$time_factor <- factor(res_long$time_step, 
                                  levels = unique(res_long$time_step),
                                  labels = paste("Time Step:", unique(res_long$time_step)))


# Create the tmap animation object
# We use a static background for the entire grid (grd)
# Then layer the burning cells from res_long, colored by 'burning' value
# Create the tmap animation object

tmap_options(facet.max = 100)

fire_animation <- tm_shape(grd) +
  tm_fill("grey90", fill_alpha = 0.5) + # Use fill_alpha for tmap v4
  tm_borders("white", lwd = 0.5) +
  
  tm_shape(res_long) +
  tm_fill("burning",
          palette = c("grey90","yellow", "red"),
          breaks = c(0, 0.01, 0.5, 1),
          labels = c("Unburned", "Low", "Medium", "High"),
          # Migrate title and other legend arguments to fill.legend
          fill.legend = tm_legend(
            fill.title = "Burning Advancement",
            fill.colorNA.show = FALSE # Explicitly hide NA in legend
          ),
          # Explicitly set colorNA to transparent or a non-plotting value
          colorNA = "transparent" # Or "grey90", or "#FFFFFF00" for full transparency
  ) +
  tm_borders(col = "white", lwd = 0.5) +
  tm_title("Fire Spread Simulation") +
  
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right"
  ) +
  tm_facets(by = "time_factor", nrow = 1, ncol = 1,
            free.coords = FALSE,
            animation.frame = "Time Step: {.$time_step} (Total Time: {.$time_step * time_step} minutes)")


tmap_animation(fire_animation, filename = "fire_spread_tmap_animation.gif", 
               width=600, height = 600, fps = 10, outer.margins = 0)
} 