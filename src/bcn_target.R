# Barcelona Spatial targeting of salvage logging

##### Setup #####

rm(list=ls())

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

source("./src/bcn_funs.R")

# libraries
library(tidyverse)
library(sf)
library(terra)

# load data
grd <- st_read("./data/rdat/grd_filt.gpkg")

set.seed(123)
grd <- grd %>%
  mutate(
    inaccessibility = scale(slope_max, center = F) + scale(road_distance, center = F),
    cost_scaler = inaccessibility / max(inaccessibility),
    sal_costs_ha = 1750 + (400 * cost_scaler),
    sal_costs_cell = sal_costs_ha  * lc_wildland * 4,
    random_target = runif(nrow(grd), 0, 1),
    wui_risk_score = wui_nohouse*0.05 + wui_verylow*0.2 + wui_low*0.4 + 
      wui_medium*0.6 + wui_intermix*0.8 + wui_interface*1
  )


# iterate over thresholds in 10% increase
years <- 2040:2050
thresholds <- seq(0.1, 1, by = 0.1)
strategies <- c("random_target", "sal_costs_cell", "wui_risk_score", "biomass", "connectivity")
order_decreasing <- c(FALSE, FALSE, TRUE, TRUE, TRUE)

target_df <- grd %>% st_set_geometry(NULL)

# run targeting
for (y in years) {
  cat(paste0("Year: ", y, "\n"))
  for (s in 1:length(strategies)) {
    
    # optional: target by temporally explicit ROS/CON
    if(strategies[s] == "biomass"){
      strtgy <- paste0("bau_ROS_", y)
    } else if(strategies[s] == "connectivity"){
      strtgy <- paste0("bmr_FireConn", y,  "b_raster")
    } else {
      strtgy <- strategies[s]
    }
    cat(paste0("- Target: ", strtgy, " (decr.:", order_decreasing[s], ") \n"))
    
    
    for (i in thresholds) {
      
      #cat(paste0("Threshold: ", i))
      
      target_df <- BCN_target_df(
        data = target_df,
        original_col = c( paste0("bau_ROS_", y), paste0("bmr_FireConn", y, "b_raster") ),
        mask_col = paste0("Logging_", y),
        replacement_col = c( paste0("sal_ROS_", y), paste0("bmr_FireConn", y, "s_raster") ),
        target_col = strtgy, # Optional: if not provided, random assignment is used
        threshold = i,
        budget = NULL,
        relative = TRUE,
        year = y,
        target_decreasing = order_decreasing[s],
        verbose = F
      )
    }
  }
}


saveRDS(target_df, file = "./data/rdat/target_df.rds")
target_df <- readRDS(file = "./data/rdat/target_df.rds")


rm_cols <- names(grd)[!names(grd) %in% c(names(grd)[1:11], "sal_costs_cell", "geom")]

# --- Reshape the dataframe from wide to long format ---
target_df_long = target_df %>%
  select(-all_of(rm_cols)) %>%
  #st_set_geometry(NULL) %>%
  pivot_longer(
    cols = matches("^(new_mask|ros_trg|con_trg)_[a-zA-Z0-9_]+_t\\d+_y\\d+$"), # Select columns matching the new pattern
    names_to = c(".value", "strategy", "threshold", "year"), # Use .value, and new 'strategy' and 'threshold' columns
    names_pattern = "^(new_mask|ros_trg|con_trg)_([a-zA-Z0-9_]+)_t(\\d+)_y(\\d+)$" # Regex to capture type, strategy, and threshold
  ) %>%
  mutate(
    threshold = as.numeric(threshold), # Convert the 'threshold' column to numeric
    sal_costs_cell = new_mask * sal_costs_cell
  )

# aggregate by threshold (mean and sum)
target_df_agg <- target_df_long %>%
  #st_set_geometry(NULL) %>%
  group_by(strategy, threshold, year) %>%
  summarise(
    sum_SAL_ha = sum(new_mask * lc_wildland * 4, na.rm = TRUE),
    sum_SAL_costs = sum(sal_costs_cell * new_mask, na.rm = TRUE)
  ) |>
  # rename strategies to more readable
  mutate(strategy = case_when(
    strategy == "random_target" ~ "Random",
    strategy == "sal_costs_cell"  ~ "Low logging costs",
    strategy == "road_distance" ~ "Lowest road distance",
    strategy == "wui_risk_score" ~ "Highest artificial",
    grepl("biomass|bau_ROS", strategy)  ~ "Highest biomass",
    grepl("connectivity|FireConn", strategy)  ~ "Highest connectivity",
    TRUE ~ strategy
  ))


# maximum sum_SAL_costs per year
target_df_agg$max_cost <- ave(target_df_agg$sum_SAL_costs, target_df_agg$year, FUN = max)
target_df_agg$max_cost_scaled <- round(target_df_agg$max_cost / max(target_df_agg$max_cost), 4)
target_df_agg$max_sal_area_km2 <- ave(target_df_agg$sum_SAL_ha, target_df_agg$year, FUN = max) / 100 # convert ha to km2
target_df_agg$facet_lab <- paste0(target_df_agg$year, " (max. area: ", round(target_df_agg$max_sal_area_km2, 2), " km²)")

# save data
saveRDS(target_df_agg, file = "./data/rdat/target_df_agg.rds")
target_df_agg <- readRDS(file = "./data/rdat/target_df_agg.rds")


# aggregate years (only sum_SAL_ha and sum_SAL_costs)
target_df_agg_annual <- target_df_agg %>%
  filter(threshold==100, strategy == "Random")


# barplot total area per year
p.sal.area <- ggplot(target_df_agg_annual, aes(x = factor(year), y = sum_SAL_ha)) +
  geom_bar(stat = "identity", position = "dodge") +
  # format Y to show values in km2
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k"),
                     expand = expansion(mult = c(0, 0.2))) + # Convert ha to km2
  # add costs as label above each bar
  geom_text(aes(label = scales::label_number(scale = 1e-6, suffix = "M")(sum_SAL_costs)), 
            vjust = -0.5, size = 3.5) +
  labs(
    x = "Share of cells salvage-logged (%)",
    y = "Targeted area (ha)"
  ) +
  theme_minimal()
p.sal.area
ggsave(p.sal.area, filename = "./out/fig/sal_area_per_year.png", 
       width = 6, height = 3, dpi = 300, bg = "white")


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

# 1. Define the color palette and labels for the 'strategy' variable
strategy_colors <- c(
  "Highest biomass"  = "darkgreen", 
  "Highest connectivity" = "purple",
  "Highest artificial" = "red",
  "Random" = "orange",
  "Low logging costs" = "blue"
)

# make area plot
p.rel.smooth <- ggplot(target_df_rel, aes(x = threshold, y = rel_sum_SAL_costs*100, color = strategy)) +
  geom_point(alpha = 0.2) +
  # smooth line aggregating years
  geom_smooth(method = "loess", se = TRUE, lwd = 1.2) +
  # add horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Share of cells salvage-logged (%)",
    y = "Relative costs (%)"
  ) +
  theme_minimal() +
  # multiply y-axis by 100 to express in percentage
  scale_y_continuous(labels = scales::label_percent(scale = 1, accuracy = 1), 
                     breaks = seq(-70, 10, by = 10),
                     expand = expansion(mult = c(0, 0.1))) + # Expand y-axis to avoid clipping
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  scale_color_manual(values = strategy_colors) #+
  theme(
    legend.position = c(0.975, 0.025), # Placed legend in the approximate 12th facet position (bottom-right)
    legend.justification = c("right", "bottom") # Adjust justification to place it correctly within the box
  )
p.rel.smooth

# save smooth plot
ggsave(p.rel.smooth, filename = "./out/fig/cum_sal_costs_relative_smooth.png", 
       width = 6, height = 4, dpi = 300, bg = "white")



# plot absolute area on X, Costs on Y, color strategy, symbol as year and size as threshold
p.sal.area.costs <- ggplot(target_df_agg, aes(x = factor(threshold), y = sum_SAL_costs, color = strategy)) +
  geom_boxplot(position = "dodge") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"),
                     transform = "log") + # Convert EUR to million EUR
  # horizontal line at y = 1000000
  geom_hline(yintercept = 1000000, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 10000000, linetype = "dashed", color = "grey50") +
  # add text labels for hlines
  annotate("text", x = factor(70), y = 1500000, label = "1M EUR", color = "grey50", size = 3) +
  annotate("text", x = factor(30), y = 15000000, label = "10M EUR", color = "grey50", size = 3) +
  labs(
    x = "",
    y = "Total costs (€)",
    shape = "Year",
    size = "Threshold (%)"
  ) +
  theme_minimal()+
  scale_color_manual(values = strategy_colors)

p.sal.area.costs

ggsave(p.sal.area.costs, filename = "./out/fig/sal_area_costs.png", 
       width = 6, height = 4, dpi = 300, bg = "white")


# combine plots
library(ggpubr)

p.combined <- ggarrange(
  p.sal.area.costs,
  p.rel.smooth,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  align = "hv",
  common.legend = TRUE, legend = "bottom"
)
p.combined
ggsave(p.combined, filename = "./out/fig/sal_area_costs_combined.png", 
       width = 7, height = 6, dpi = 300, bg = "white")
