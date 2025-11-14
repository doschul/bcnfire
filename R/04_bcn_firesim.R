# ==============================================
# 04_bcn_firesim.R
# Description  : Run fire spread scenarios across ignition cells/years/strategies and aggregate landcover impacts per timestep.
# Inputs       : data/rdat/target_df.rds, data/rdat/neighbor_list.rds
# Outputs      : data/rdat/long_df_agg.rds, data/rdat/perc_data_long.rds, out/figures/cost_effectiveness_heatmap.png
# Author       : Dario Schulz
# Created      : 31.10.2025
# ==============================================

rm(list = ls())
setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

# ---------- Setup ----------
suppressPackageStartupMessages({
  library(stringr)
  library(here)
  library(tidyverse)
  library(future)
  library(future.apply)
  library(progressr)
  library(data.table)
  library(parallel)
})

set.seed(123)
options(stringsAsFactors = FALSE)
# calculate 850mb limit:
# 850*1024^2 = 891289600
options(future.globals.maxSize= 891289600, 
        future.seed = TRUE)

handlers(handler_progress(format="[:bar] :percent :eta :message"))

# Project root assumed to be Git repo top; use here::here()
data_raw_dir  <- here::here("data", "raw")
data_proc_dir <- here::here("data", "rdat")
out_fig_dir   <- here::here("out", "fig")
out_tab_dir   <- here::here("out", "tab")

# Source helper functions
funs_path <- here::here("src", "bcn_funs.R")
if (file.exists(funs_path)) source(funs_path)

# ---------- Analysis ----------
# ----------------------------------
# Load data
# ----------------------------------
target_df <- readRDS("./data/rdat/target_df.rds")
load("./data/rdat/neighbor_idx_differentiated.RData")

# Add scenarios
for (y in 2040:2050) {
  target_df[[paste0("ros_trg_baseline_y", y)]] <- target_df[[paste0("bau_ROS_", y)]]
  target_df[[paste0("con_trg_baseline_y", y)]] <- target_df[[paste0("bmr_FireConn", y, "b_raster")]]
  target_df[[paste0("ros_trg_full_y", y)]]     <- target_df[[paste0("sal_ROS_", y)]]
  target_df[[paste0("con_trg_full_y", y)]]     <- target_df[[paste0("bmr_FireConn", y, "s_raster")]]
}

# ----------------------------------
# Parameters
# ----------------------------------
set.seed(123)
rdm_ign_cells <- sample(
  x = row.names(target_df),
  size = 1000,
  replace = FALSE,
  prob = target_df$e_ign_prob
)

batch_size     <- 1   # process 20 ignitions at a time
mc_time_step   <- 5
mc_time_horizon<- 180

all_scenarios <- names(target_df) |>
  (\(x) x[grepl("^(ros_trg_|con_trg_)", x)])() |>
  (\(x) sub("^(ros_trg_|con_trg_)", "", x))() |>
  unique()

#all_scenarios <- all_scenarios[grepl("y2046", all_scenarios)]
# ----------------------------------
# Prepare batching + output
# ----------------------------------
out_dir <- "./data/rdat/mc_fire_batches/"
dir.create(out_dir, showWarnings = FALSE)
ignition_batches <- split(rdm_ign_cells, ceiling(seq_along(rdm_ign_cells) / batch_size))

# ----------------------------------
# Land cover data
# ----------------------------------
target_df$cell_id <- row.names(target_df)
lc_dat <- as.data.table(target_df)
#lc_dat$cell_id <- paste0("cell_", lc_dat$cell_id)
#lc_dat[, cell_id := as.character(.I)]
lc_dat <- lc_dat[, .(cell_id, lc_agriculture, lc_wildland, lc_artificial,
                     wui_nohouse, wui_verylow, wui_low, 
                     wui_medium, wui_intermix, wui_interface)]


# ----------------------------------
# Batch processor (aggregates directly)
# ----------------------------------
process_batch <- function(batch_id, ign_cells, scenarios) {
  batch_agg_list <- vector("list", length(scenarios))
  
  for (s in seq_along(scenarios)) {
    s_name <- scenarios[s]
    scen_agg_list <- vector("list", length(ign_cells))
    
    for (j in seq_along(ign_cells)) {
      df <- as.data.table(
        get_burners_optimized(
          ign_cells[j],
          time_horizon = mc_time_horizon,
          time_step    = mc_time_step,
          full         = TRUE,
          scenario     = s_name,
          grd          = target_df,
          con_scaler   = 2000,
          neighbor_list_differentiated = neighbor_idx_differentiated
        )
      )
      
      ignition_id <- df$id[df$V1 == 1]
      
      df_long <- melt(df, id.vars = "id",
                      variable.name = "time", value.name = "burned")
      df_long[, time := as.numeric(gsub("V", "", time)) * mc_time_step]
      df_long[, id := as.character(id)]
      
      df_long <- merge(df_long, lc_dat, by.x = "id", by.y = "cell_id", all.x = TRUE)
      
      df_agg <- df_long[, .(
        burned_land_ha = sum(burned * 4, na.rm = TRUE),
        burned_agri_ha = sum(burned * lc_agriculture * 4, na.rm = TRUE),
        burned_wild_ha = sum(burned * lc_wildland * 4, na.rm = TRUE),
        burned_urbn_ha = sum(burned * lc_artificial * 4, na.rm = TRUE),
        burned_wui_nohouse_ha   = sum(burned * wui_nohouse * 4, na.rm = TRUE),
        burned_wui_verylow_ha   = sum(burned * wui_verylow * 4, na.rm = TRUE),
        burned_wui_low_ha       = sum(burned * wui_low * 4, na.rm = TRUE),
        burned_wui_medium_ha    = sum(burned * wui_medium * 4, na.rm = TRUE),
        burned_wui_intermix_ha  = sum(burned * wui_intermix * 4, na.rm = TRUE),
        burned_wui_interface_ha = sum(burned * wui_interface * 4, na.rm = TRUE)
      ), by = .(time)]
      
      df_agg[, scenario := s_name]
      df_agg[, ignition_id := ignition_id]
      
      scen_agg_list[[j]] <- df_agg
    }
    
    batch_agg_list[[s]] <- rbindlist(scen_agg_list)
  }
  
  batch_agg <- rbindlist(batch_agg_list)
  saveRDS(batch_agg, file = file.path(out_dir, paste0("batch_lc_agg_", batch_id, ".rds")))
  invisible(TRUE)
}

# ----------------------------------
# Run in parallel with progress bar
# ----------------------------------
done_batches <- list.files(out_dir, pattern = "batch_lc_agg_.*\\.rds") |>
  gsub("batch_lc_agg_|\\.rds", "", x = _) |>
  as.integer()

todo_batches <- setdiff(seq_along(ignition_batches), done_batches)

cl <- makeClusterPSOCK(28)
clusterExport(cl, c("target_df","neighbor_idx_differentiated",
                    "get_burners_optimized","mc_time_horizon","mc_time_step","lc_dat"))
plan(cluster, workers = cl)

handlers(global = TRUE)
with_progress({
  p <- progressor(along = todo_batches)
  future_lapply(todo_batches, function(b) {
    ign_cells <- ignition_batches[[b]]
    process_batch(b, ign_cells, all_scenarios)
    p(sprintf("Completed batch %d", b))
  })
})

plan(sequential)
parallel::stopCluster(cl)

# ----------------------------------
# Final combine
# ----------------------------------
agg_files <- list.files(out_dir, pattern = "batch_lc_agg_.*\\.rds", full.names = TRUE)
all_agg_results <- rbindlist(lapply(agg_files, readRDS))



# extract variable information from scenarios
long_df_agg <- all_agg_results %>%
  mutate(
    # year = last four characters of scenario
    year = substr(scenario, nchar(scenario) - 3, nchar(scenario)),
    # logging intensity 0 for baseline, 100 for full
    # pattern between _t and _y
    sal_int = case_when(grepl("baseline", scenario) ~ "0",
                        grepl("full", scenario) ~ "100",
                        TRUE ~ sub(".*_t(\\d+).*", "\\1", scenario)),
    # strategy is string before _t pattern starts
    strategy = sub("_(t|y).*", "", scenario),
    strategy = case_when(
      grepl("bau_ROS", strategy) ~ "High Biomass",
      grepl("bmr_FireConn", strategy) ~ "High Connectivity",
      grepl("random", strategy) ~ "Random",
      grepl("sal_costs_cell", strategy) ~ "Low SAL cost",
      grepl("wui_risk_score", strategy) ~ "High Built-up",
      TRUE ~ strategy
    ),
    # year and sal_int as integers
    year = as.integer(year),
    sal_int = as.integer(sal_int),
    time = as.integer(time) # Ensure time is an integer
  )

#saveRDS(long_df_agg, "./data/rdat/mc_res_agg_WUI_TEST.rds")

saveRDS(long_df_agg, "./data/rdat/mc_res_agg.rds")

mc_res_agg_random_sal0 <- long_df_agg %>% 
  filter(sal_int == 0, 
         strategy == "baseline") %>%
  select(ignition_id, time, burned_land_ha, year)


write_csv(mc_res_agg_random_sal0, 
          "./data/rdat/mc_res_agg_random_sal0.csv")
mc_res_agg_random_sal0 <- read_csv("./data/rdat/mc_res_agg_random_sal0.csv")

t1 <- mc_res_agg_random_sal0

### Aggregate Percentiles with data.table



# ============================================================
# 0. Load & prepare
# ============================================================
long_df_agg <- readRDS("./data/rdat/mc_res_agg.rds")
# remove duplicate
long_df_agg <- long_df_agg %>%
  filter(!strategy == "full")

setDT(long_df_agg)

# ============================================================
# 1. Duplicate baseline/full rows for each non-(baseline|full) strategy
# ============================================================

unique_strategies <- unique(long_df_agg[!grepl("baseline|full", strategy), strategy])
baseline_rows <- long_df_agg[grepl("baseline", scenario)]
#full_rows     <- long_df_agg[grepl("full", scenario)]

# Helper table for join
sDT <- data.table(strategy = unique_strategies, dummy = 1L)

# Add dummy column to template rows
baseline_rows[, dummy := 1L]
#full_rows[, dummy := 1L]

# Cartesian merge separately
base_expanded <- sDT[baseline_rows, on = "dummy", allow.cartesian = TRUE][
  , dummy := NULL][]

#full_expanded <- sDT[full_rows, on = "dummy", allow.cartesian = TRUE][
#  , dummy := NULL][]

# Combine and append
#rep_rows <- rbindlist(list(base_expanded, full_expanded), use.names = TRUE, fill = TRUE)

# Append and drop original baseline/full
long_df_agg <- rbindlist(list(long_df_agg, base_expanded), use.names = TRUE, fill = TRUE)
cat("After append:", nrow(long_df_agg), "\n")

long_df_agg <- long_df_agg[!grepl("baseline|full", strategy)]
cat("After removing originals:", nrow(long_df_agg), "\n")

# Optional sanity check: intensity distribution
print(table(long_df_agg$sal_int))

# ============================================================
# 2. Compute avoided burned areas FIRST (before percentiles)
# ============================================================
# Identify measures and grouping vars
meas <- c("burned_land_ha", "burned_agri_ha", "burned_wild_ha", "burned_urbn_ha",
          "burned_wui_nohouse_ha", "burned_wui_verylow_ha", "burned_wui_low_ha",
          "burned_wui_medium_ha", "burned_wui_intermix_ha",  "burned_wui_interface_ha")

# Add baseline values (sal_int==0) within each (year, time, strategy, ignition ID)
id_col <- "ignition_id"

long_df_agg[
  , paste0("base_", meas) := lapply(.SD, function(x) x[sal_int == 0]),
  by = .(year, time, strategy, get(id_col)),
  .SDcols = meas
]

# Compute absolute avoided areas (BAU - scenario)
long_df_agg[
  , c("abs_avoid_tot_ha", "abs_avoid_agri_ha", "abs_avoid_wild_ha", "abs_avoid_urbn_ha",
      "abs_avoid_wui_nohouse_ha", "abs_avoid_wui_verylow_ha", "abs_avoid_wui_low_ha",
      "abs_avoid_wui_medium_ha", "abs_avoid_wui_intermix_ha",  "abs_avoid_wui_interface_ha") :=
    .(base_burned_land_ha - burned_land_ha,
      base_burned_agri_ha - burned_agri_ha,
      base_burned_wild_ha - burned_wild_ha,
      base_burned_urbn_ha - burned_urbn_ha,
      
      base_burned_wui_nohouse_ha - burned_wui_nohouse_ha,
      base_burned_wui_verylow_ha - burned_wui_verylow_ha,
      base_burned_wui_low_ha - burned_wui_low_ha,
      base_burned_wui_medium_ha - burned_wui_medium_ha,
      base_burned_wui_intermix_ha - burned_wui_intermix_ha,
      base_burned_wui_interface_ha - burned_wui_interface_ha)
]

# ============================================================
# 3. Sanity check: avoided values should be >= 0
# ============================================================
check_cols <- c("abs_avoid_tot_ha", "abs_avoid_agri_ha", "abs_avoid_wild_ha", "abs_avoid_urbn_ha",
                "abs_avoid_wui_nohouse_ha", "abs_avoid_wui_verylow_ha", "abs_avoid_wui_low_ha",
                "abs_avoid_wui_medium_ha", "abs_avoid_wui_intermix_ha",  "abs_avoid_wui_interface_ha")
for (v in check_cols) {
  neg_n <- long_df_agg[get(v) < -1e-6, .N]
  if (neg_n > 0) {
    warning(sprintf("Found %d negative avoided values in %s!", neg_n, v))
  } else {
    message(sprintf("Check passed: no negative values in %s.", v))
  }
}

# ============================================================
# 4. Compute percentiles over avoided values
# ============================================================
probs <- seq(0, 1, by = 0.01)
pnames <- paste0("p", sprintf("%02d", as.integer(round(probs * 100, 0))))
avoid_meas <- check_cols

# Group by higher level (year, strategy, sal_int, time)
perc_data_wide <- long_df_agg[
  , {
    out <- list()
    for (m in avoid_meas) {
      qs <- quantile(.SD[[m]], probs = probs, na.rm = TRUE)
      names(qs) <- paste0(m, "_", pnames)
      out <- c(out, as.list(qs))
    }
    out
  },
  by = .(year, strategy, sal_int, time),
  .SDcols = avoid_meas
]



# ============================================================
# 5. Pivot to long (optional for plotting)
# ============================================================
perc_data_long <- melt(
  perc_data_wide,
  id.vars = c("year", "strategy", "sal_int", "time"),
  measure.vars = patterns(
    "^abs_avoid_tot_ha_",
    "^abs_avoid_agri_ha_",
    "^abs_avoid_wild_ha_",
    "^abs_avoid_urbn_ha_",
    "^abs_avoid_wui_nohouse_ha",
    "^abs_avoid_wui_verylow_ha",
    "^abs_avoid_wui_low_ha",
    "^abs_avoid_wui_medium_ha",
    "^abs_avoid_wui_intermix_ha",
    "^abs_avoid_wui_interface_ha"
  ),
  variable.name = "p_idx",
  value.name   = c("abs_avoid_tot_ha", "abs_avoid_agri_ha", "abs_avoid_wild_ha", "abs_avoid_urbn_ha",
                   "abs_avoid_wui_nohouse_ha", "abs_avoid_wui_verylow_ha", "abs_avoid_wui_low_ha",
                   "abs_avoid_wui_medium_ha", "abs_avoid_wui_intermix_ha",  "abs_avoid_wui_interface_ha")
)[
  , `:=`(
    percentile      = probs[p_idx],
    percentile_name = paste0("p", sprintf("%02d", p_idx))
  )
][]



# ============================================================
# 6. Save
# ============================================================
saveRDS(perc_data_long, "./data/rdat/perc_data_long.rds")









# 2. Calculate time-scenario specific summaries

# --- 2. Prepare Data for Spaghetti Plot (Calculate Percentiles) ---
# For a spaghetti plot of percentiles, you need to calculate quantiles
# for 'hectares' for each unique combination of 'year', 'strategy', and 'sal_int'.

summary_data <- long_df_agg %>%
  group_by(year, strategy, sal_int, time) %>%
  summarise(
    # Calculate all desired percentiles
    across(.cols = c(burned_land_ha, burned_agri_ha,
                     burned_wild_ha, burned_urbn_ha),
           .fns = list(
             p50 = ~quantile(., probs = 0.50, na.rm = TRUE),
             p90 = ~quantile(., probs = 0.90, na.rm = TRUE),
             p95 = ~quantile(., probs = 0.95, na.rm = TRUE),
             p99 = ~quantile(., probs = 0.99, na.rm = TRUE)
           ),
           # FIX: Use "{.col}_{.fn}" to prepend original column name
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  # Now pivot longer in two steps:
  # 1. Pivot to separate the original burned_type from the percentile type
  pivot_longer(
    cols = ends_with(c("_p50", "_p90", "_p95", "_p99")), # Select all newly created percentile columns
    names_to = c("burned_type", "percentile_type"), # Create two new columns
    names_pattern = "(.*)_(p\\d+)", # Regex to capture the burned type and the percentile
    values_to = "hectares_value"
  ) %>%
  # Convert percentile_type to a factor for ordered legends
  mutate(percentile_type = factor(percentile_type,
                                  levels = paste0("p", c(50, 90, 95, 99)),
                                  labels = c("Median (50th Percentile)", "90th Percentile", "95th Percentile", "99th Percentile")),
         # Optional: Clean up burned_type for better labels if needed
         burned_type = gsub("_ha$", "", burned_type) # Remove "_ha" suffix if desired
  )



# --- 3. Create the Spaghetti Plot ---

ggplot(summary_data %>%
         filter(time %in% c(60, 120, 180)), aes(x = sal_int, y = hectares_value, color = strategy, linetype = percentile_type, group = interaction(strategy, percentile_type))) +
  geom_line(size = 0.8) + # Lines for trends
  geom_point(size = 2, aes(shape = percentile_type)) + # Points for specific budget steps
  facet_wrap(year ~ time, scales = "free_y") + # Facet by year, adjust ncol as needed
  # You might want to manually set colors for strategies if you have a specific palette
  # scale_color_manual(values = c("Strategy_A" = "blue", "Strategy_B" = "red", ...)) +
  scale_linetype_manual(values = c("Median (50th Percentile)" = "solid",
                                   "90th Percentile" = "dashed",
                                   "95th Percentile" = "dotted",
                                   "99th Percentile" = "dotdash")) +
  labs(
    title = "Wildfire Hectares by Targeting Budget, Scenario, and Year (Percentiles)",
    x = "Logging (%)",
    y = "Hectares Burned",
    color = "Scenario", # Legend title for colors
    linetype = "Percentile Type", # Legend title for linetypes
    shape = "Percentile Type" # Legend title for shapes
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"), # Make facet labels bold
    panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if they overlap
  )


# express relative to BAU scenario, i.e., sal_int = 0
# 1. Prepare Data for Difference Calculation
# First, we need to summarize the data for each scenario and time point


summary_relative <- summary_data %>%
  group_by(burned_type, strategy, year, time, percentile_type) %>%
  mutate(area_relative = hectares_value / hectares_value[sal_int == 0] - 1) %>%
  ungroup()

saveRDS(summary_relative, file = "./data/rdat/summary_relative.rds")
summary_relative <- read_rds("./data/rdat/summary_relative.rds")

# 2. Create the Spaghetti Plot for Relative Areas

perc_to_plot <- c("Median (50th Percentile)", "99th Percentile")
lc_to_plot <- c("burned_agri", "burned_wild", "burned_urbn")

# instead of year facet, make smooth line across years and aggregate
ggplot(summary_relative %>%
         filter(time %in% c(120, 180),
                percentile_type %in% perc_to_plot,
                burned_type %in% lc_to_plot), 
       aes(x = sal_int, y = area_relative, 
           color = strategy, linetype = percentile_type)) +
  geom_smooth(aes(weight = time),
              method = "loess", size = 0.8, se = F) + # Smooth lines for trends
  geom_point(size = 1, alpha = 0.2) +
  facet_wrap(year ~ burned_type, scales = "free_y") + # Facet by year, adjust ncol as needed
  scale_linetype_manual(values = c("Median (50th Percentile)" = "solid",
                                   "90th Percentile" = "dashed",
                                   "95th Percentile" = "dotted",
                                   "99th Percentile" = "dotdash")) +
  labs(
    title = "Relative Wildfire Hectares by Targeting Budget and Scenario (Percentiles)",
    x = "Logging (%)",
    y = "Relative Hectares Burned (relative to BAU)",
    color = "Scenario", # Legend title for colors
    linetype = "Percentile Type", # Legend title for linetypes
    shape = "Percentile Type" # Legend title for shapes
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"), # Make facet labels bold
    panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if they overlap
  )
