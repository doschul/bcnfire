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
funs_path <- here::here("src", "R", "bcn_funs.R")
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

saveRDS(long_df_agg, "./data/rdat/mc_res_agg.rds")

mc_res_agg_random_sal0 <- long_df_agg %>% 
  filter(sal_int == 0, 
         strategy == "baseline") %>%
  select(ignition_id, time, burned_land_ha, year)


write_csv(mc_res_agg_random_sal0, 
          "./data/rdat/mc_res_agg_random_sal0.csv")
mc_res_agg_random_sal0 <- read_csv("./data/rdat/mc_res_agg_random_sal0.csv")

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
# 4. Compute percentiles over avoided AND burned values
# ============================================================
probs <- seq(0, 1, by = 0.01)
pnames <- paste0("p", sprintf("%02d", as.integer(round(probs * 100, 0))))
avoid_meas <- check_cols
burned_meas <- c("burned_land_ha", "burned_agri_ha", "burned_wild_ha", "burned_urbn_ha")
all_meas <- c(avoid_meas, burned_meas)

# Group by higher level (year, strategy, sal_int, time)
perc_data_wide <- long_df_agg[
  , {
    out <- list()
    for (m in all_meas) {
      qs <- quantile(.SD[[m]], probs = probs, na.rm = TRUE)
      names(qs) <- paste0(m, "_", pnames)
      out <- c(out, as.list(qs))
    }
    out
  },
  by = .(year, strategy, sal_int, time),
  .SDcols = all_meas
]



# ============================================================
# 5. Pivot to long
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
    "^abs_avoid_wui_interface_ha",
    "^burned_land_ha_",
    "^burned_agri_ha_",
    "^burned_wild_ha_",
    "^burned_urbn_ha_"
  ),
  variable.name = "p_idx",
  value.name   = c("abs_avoid_tot_ha", "abs_avoid_agri_ha", "abs_avoid_wild_ha", "abs_avoid_urbn_ha",
                   "abs_avoid_wui_nohouse_ha", "abs_avoid_wui_verylow_ha", "abs_avoid_wui_low_ha",
                   "abs_avoid_wui_medium_ha", "abs_avoid_wui_intermix_ha",  "abs_avoid_wui_interface_ha",
                   "burned_land_ha", "burned_agri_ha", "burned_wild_ha", "burned_urbn_ha")
)[
  , `:=`(
    percentile      = probs[p_idx],
    percentile_name = paste0("p", sprintf("%02d", p_idx))
  )
][]



# ============================================================
# 6. Save unified percentile data
# ============================================================
saveRDS(perc_data_long, "./data/rdat/perc_data_long.rds")


# ============================================================
# ============================================================
# ENTRY POINT: everything below runs from saved .rds files
# ============================================================
# ============================================================

# -----------------------------
# 7. Load saved data
# -----------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(mgcv)
library(patchwork)
library(scales)
library(tibble)
library(gt)
library(data.table)
library(ggpubr)

options(scipen = 999)

perc_data_long <- readRDS("./data/rdat/perc_data_long.rds")
target_df_agg  <- readRDS(file = "./data/rdat/target_df_agg.rds")

# rename strategies to match
target_df_agg <- target_df_agg %>%
  mutate(strategy = case_when(
    strategy == "Highest artificial" ~ "High Built-up",
    strategy == "Highest biomass" ~ "High Biomass",
    strategy == "Highest connectivity" ~ "High Connectivity",
    strategy == "Low logging costs" ~ "Low SAL cost",
    TRUE ~ strategy
  ))

# -----------------------------
# 8. Derive summary_relative from perc_data_long
# -----------------------------
pct_map <- c(`0.5` = "Median (50th Percentile)",
             `0.9` = "90th Percentile",
             `0.95` = "95th Percentile",
             `0.99` = "99th Percentile")

summary_data <- perc_data_long %>%
  filter(percentile %in% c(0.50, 0.90, 0.95, 0.99)) %>%
  select(year, strategy, sal_int, time, percentile,
         burned_land_ha, burned_agri_ha, burned_wild_ha, burned_urbn_ha) %>%
  pivot_longer(
    cols = c(burned_land_ha, burned_agri_ha, burned_wild_ha, burned_urbn_ha),
    names_to = "burned_type",
    values_to = "hectares_value"
  ) %>%
  mutate(
    burned_type = gsub("_ha$", "", burned_type),
    percentile_type = factor(
      pct_map[as.character(percentile)],
      levels = c("Median (50th Percentile)", "90th Percentile",
                 "95th Percentile", "99th Percentile")
    )
  ) %>%
  select(-percentile)

summary_relative <- summary_data %>%
  group_by(burned_type, strategy, year, time, percentile_type) %>%
  mutate(area_relative = hectares_value / hectares_value[sal_int == 0] - 1) %>%
  ungroup()

saveRDS(summary_relative, file = "./data/rdat/summary_relative.rds")


# =========================================================
# 9. Build GAM dataset from perc_data_long
# =========================================================

# Helper functions
tidy_gam_parametric <- function(model) {
  sm <- summary(model)
  out <- as.data.frame(sm$p.table)
  out <- rownames_to_column(out, "term")
  names(out) <- c("term", "estimate", "std_error", "statistic", "p_value")
  out %>%
    mutate(
      conf_low  = estimate - 1.96 * std_error,
      conf_high = estimate + 1.96 * std_error
    )
}

tidy_gam_smooths <- function(model) {
  sm <- summary(model)
  out <- as.data.frame(sm$s.table)
  out <- rownames_to_column(out, "smooth")
  names(out) <- names(out) %>%
    str_replace("Ref.df", "ref_df") %>%
    str_replace("p-value", "p_value") %>%
    str_replace("^F$", "F_stat")
  out
}

label_factor_terms <- function(df) {
  df %>%
    mutate(
      group = case_when(
        term == "(Intercept)" ~ "Intercept",
        str_detect(term, "^strategy") ~ "Strategy",
        str_detect(term, "^burned_type") ~ "Avoid type",
        str_detect(term, "^year") ~ "Year",
        TRUE ~ "Other"
      ),
      level = case_when(
        str_detect(term, "^strategy") ~ str_remove(term, "^strategy"),
        str_detect(term, "^burned_type") ~ str_remove(term, "^burned_type"),
        str_detect(term, "^year") ~ str_remove(term, "^year"),
        TRUE ~ term
      )
    )
}

make_partial_curve <- function(model, data, focal_var,
                               n = 200,
                               strategy_ref = "Random",
                               burned_ref = NULL,
                               year_ref = NULL) {
  
  if (is.null(burned_ref)) burned_ref <- levels(data$burned_type)[1]
  if (is.null(year_ref)) year_ref <- names(sort(table(data$year), decreasing = TRUE))[1]
  
  nd <- tibble(
    time = median(data$time, na.rm = TRUE),
    percentile_num = median(data$percentile_num, na.rm = TRUE),
    sal_int = median(data$sal_int, na.rm = TRUE),
    strategy = factor(strategy_ref, levels = levels(data$strategy)),
    burned_type = factor(burned_ref, levels = levels(data$burned_type)),
    year = factor(year_ref, levels = levels(data$year))
  )
  
  vals <- seq(
    min(data[[focal_var]], na.rm = TRUE),
    max(data[[focal_var]], na.rm = TRUE),
    length.out = n
  )
  
  nd <- nd[rep(1, n), ]
  nd[[focal_var]] <- vals
  
  pr <- predict(model, newdata = nd, se.fit = TRUE, type = "link")
  
  nd %>%
    mutate(
      fit_link = pr$fit,
      se_link  = pr$se.fit,
      fit_resp = pmax(exp(fit_link) - 1, 0),
      lwr_resp = pmax(exp(fit_link - 1.96 * se_link) - 1, 0),
      upr_resp = pmax(exp(fit_link + 1.96 * se_link) - 1, 0)
    )
}

# --- Pivot avoided area columns to long ---
avoid_cols <- c("abs_avoid_tot_ha", "abs_avoid_urbn_ha")

hist(perc_data_long$abs_avoid_tot_ha)

gam_df <- perc_data_long %>%
  as_tibble() %>%
  select(year, strategy, sal_int, time, percentile,
         all_of(avoid_cols)) %>%
  pivot_longer(
    cols = all_of(avoid_cols),
    names_to = "burned_type",
    values_to = "avoided_ha"
  ) %>%
  mutate(
    burned_type = recode(
      burned_type,
      abs_avoid_tot_ha  = "Avoided Land (Total ha)",
      abs_avoid_urbn_ha = "Avoided Urban area (ha)"
    ),
    percentile_num = percentile * 100,
    time = as.numeric(time),
    sal_int = as.numeric(sal_int),
    year = as.character(year)
  ) %>%
  # Exclude BAU
  filter(sal_int > 0) %>%
  # Join intervention costs
  left_join(
    target_df_agg %>% mutate(year = as.character(year)),
    by = c("year" = "year",
           "strategy" = "strategy",
           "sal_int" = "threshold")
  ) %>%
  mutate(
    cost_eff     = sum_SAL_costs / avoided_ha,
    log_cost_eff = log1p(cost_eff)
  ) %>%
  # Clean / filter
  mutate(
    strategy    = fct_relevel(factor(strategy), "Random"),
    burned_type = factor(burned_type),
    year        = factor(year)
  ) %>%
  filter(
    is.finite(avoided_ha),
    is.finite(sum_SAL_costs),
    is.finite(cost_eff),
    is.finite(log_cost_eff),
    avoided_ha > 0,
    !is.na(percentile_num),
    !is.na(time),
    !is.na(sal_int)
  )

# Quick checks
print(gam_df %>% summarise(
  n = n(),
  min_time = min(time), max_time = max(time),
  min_sal  = min(sal_int), max_sal  = max(sal_int),
  min_pct  = min(percentile_num), max_pct  = max(percentile_num)
))


# =========================================================
# 10. Fit GAM + benchmark models
# =========================================================

nrow(gam_df)
head(gam_df$percentile_num)

m_gam <- gam(
  log_cost_eff ~
    s(time, k = min(8, n_distinct(gam_df$time) - 1), bs = "cr") +
    s(percentile_num, k = min(10, n_distinct(gam_df$percentile_num) - 1), bs = "cr") +
    s(sal_int, k = min(6, n_distinct(gam_df$sal_int) - 1), bs = "cr") +
    strategy +
    burned_type +
    year,
  data = gam_df |> sample_n(100000),
  method = "REML",
  select = TRUE
)

m_quad <- lm(
  log_cost_eff ~
    strategy +
    burned_type +
    year +
    poly(time, 2, raw = TRUE) +
    poly(percentile_num, 2, raw = TRUE) +
    poly(sal_int, 2, raw = TRUE),
  data = gam_df
)

model_comp <- tibble(
  model = c("GAM", "Quadratic OLS"),
  AIC = c(AIC(m_gam), AIC(m_quad)),
  BIC = c(BIC(m_gam), BIC(m_quad)),
  adj_r2 = c(summary(m_gam)$r.sq, summary(m_quad)$adj.r.squared),
  dev_explained = c(summary(m_gam)$dev.expl, NA_real_),
  n = c(nobs(m_gam), nobs(m_quad))
)

print(summary(m_gam))
print(summary(m_quad))
print(model_comp)


# =========================================================
# 11. Export supplementary tables
# =========================================================

gam_parametric <- tidy_gam_parametric(m_gam) %>%
  label_factor_terms()

gam_smooths <- tidy_gam_smooths(m_gam)

write.csv(gam_parametric, "./out/tab/gam_parametric_coefficients_runs.csv", row.names = FALSE)
write.csv(gam_smooths, "./out/tab/gam_smooth_terms_runs.csv", row.names = FALSE)
write.csv(model_comp, "./out/tab/model_comparison_runs.csv", row.names = FALSE)

gt(gam_parametric) %>%
  tab_header(title = "Supplementary Table S1. GAM parametric coefficients") %>%
  gtsave("./out/tab/gam_parametric_coefficients_runs.html")

gt(gam_smooths) %>%
  tab_header(title = "Supplementary Table S2. GAM smooth terms") %>%
  gtsave("./out/tab/gam_smooth_terms_runs.html")

gt(model_comp) %>%
  tab_header(title = "Supplementary Table S3. Model comparison") %>%
  gtsave("./out/tab/model_comparison_runs.html")


# =========================================================
# 12. Partial-effect plots (structural mechanism)
# =========================================================
strategy_ref <- if ("Random" %in% levels(gam_df$strategy)) "Random" else levels(gam_df$strategy)[1]
burned_ref   <- if ("Avoided Land (Total ha)" %in% levels(gam_df$burned_type)) {
  "Avoided Land (Total ha)"
} else {
  levels(gam_df$burned_type)[1]
}
year_ref <- names(sort(table(gam_df$year), decreasing = TRUE))[1]

curve_time <- make_partial_curve(
  model = m_gam, data = gam_df, focal_var = "time",
  strategy_ref = strategy_ref, burned_ref = burned_ref, year_ref = year_ref
)

curve_pct <- make_partial_curve(
  model = m_gam, data = gam_df, focal_var = "percentile_num",
  strategy_ref = strategy_ref, burned_ref = burned_ref, year_ref = year_ref
)

curve_sal <- make_partial_curve(
  model = m_gam, data = gam_df, focal_var = "sal_int",
  strategy_ref = strategy_ref, burned_ref = burned_ref, year_ref = year_ref
)


# =========================================================
# 12. Structural gradients (shared annotation, cleaner layout)
# =========================================================

# Optional: choose a more typical reference year for the structural panels
# rather than the outlier year, if desired. If you already created curve_* with
# a chosen year, this only affects the text shown in the subtitle.
year_ref_struct <- 2049
strategy_ref_struct <- "Random"
burned_ref_struct <- "Avoided Land (Total ha)"

struct_sub_shared <- paste0(
  "Partial GAM predictions for ", burned_ref_struct,
  "; holding strategy = ", strategy_ref_struct,
  " and year = ", year_ref_struct
)

common_curve_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(5, 8, 5, 5),
    panel.grid.minor = element_blank()
  )

p_time <- ggplot(curve_time, aes(x = time, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp), alpha = 0.15) +
  geom_line(linewidth = 1) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  common_curve_theme +
  labs(
    title = "A1) Duration",
    x = "Fire duration (time steps)",
    y = "Predicted cost per avoided ha"
  )

p_pct <- ggplot(curve_pct, aes(x = percentile_num, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp), alpha = 0.15) +
  geom_line(linewidth = 1) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "A2) Severity percentile",
    x = "Percentile of avoided burned area distribution"
  ) +
  common_curve_theme +
  theme(
    axis.title.y = element_blank()
  )

p_sal <- ggplot(curve_sal, aes(x = sal_int, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp), alpha = 0.15) +
  geom_line(linewidth = 1) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "A3) Intervention intensity",
    x = "Salvage logging intensity (%)"
  ) +
  common_curve_theme +
  theme(
    axis.title.y = element_blank()
  )

panel_A_curves <- p_time | p_pct | p_sal

panel_A <- panel_A_curves +
  plot_annotation(
    title = "A) Structural gradients from GAM",
    subtitle = struct_sub_shared,
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0)
    )
  )

# =========================================================
# 13. Adjusted year × strategy heatmap
# =========================================================
unique(gam_df$burned_type)
# choose fixed values for the continuous dimensions
time_ref <- 200
pct_ref  <- 95
sal_ref  <- 10

pred_grid <- expand.grid(
  year = levels(gam_df$year),
  strategy = levels(gam_df$strategy),
  burned_type = c("Avoided Land (Total ha)", "Avoided Urban area (ha)"),
  time = time_ref,
  percentile_num = pct_ref,
  sal_int = sal_ref,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    year = factor(year, levels = levels(gam_df$year)),
    strategy = factor(strategy, levels = levels(gam_df$strategy)),
    burned_type = factor(burned_type, levels = levels(gam_df$burned_type))
  )

pr <- predict(m_gam, newdata = pred_grid, se.fit = TRUE, type = "link")

pred_grid <- pred_grid %>%
  mutate(
    fit_link = pr$fit,
    se_link  = pr$se.fit,
    pred_cost_eff = pmax(exp(fit_link) - 1, 0),
    lwr_cost_eff  = pmax(exp(fit_link - 1.96 * se_link) - 1, 0),
    upr_cost_eff  = pmax(exp(fit_link + 1.96 * se_link) - 1, 0)
  )

p_heat_adj <- ggplot(pred_grid, aes(x = strategy, y = fct_rev(year), fill = pred_cost_eff)) +
  geom_tile() +
  facet_wrap(~ burned_type, ncol = 2) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    name = "Predicted cost\nper avoided ha"
  ) +
  labs(
    title = "B) Adjusted heterogeneity across years and strategies",
    subtitle = paste0(
      "Holding duration = ", round(time_ref, 0),
      ", percentile = ", round(pct_ref, 0),
      ", salvage logging = ", round(sal_ref, 0), "%"
    ),
    x = "Strategy",
    y = "Year"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.8, "cm"),
    plot.margin = margin(15, 10, 5, 10)
  )

# =========================================================
# 14. Combine GAM figure
# =========================================================

# Shared y-axis label for panel A using patchwork spacer trick
y_lab <- ggplot() +
  annotate(
    "text", x = 1, y = 1,
    label = "Predicted cost per avoided burned ha",
    angle = 90, fontface = "plain", size = 4
  ) +
  theme_void()

#panel_A <- p_time | p_pct | p_sal

p_main <- panel_A / p_heat_adj +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Cost-effectiveness: structural gradients versus factor heterogeneity",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    )
  )

print(p_main)

ggsave(
  filename = "./out/fig/gam_structure_heterogeneity_runs.png",
  plot = p_main,
  width = 12,
  height = 8,
  bg = "white"
)

# Diagnostics
png("./out/fig/gam_diagnostics_runs.png", width = 1400, height = 1200, res = 160)
par(mfrow = c(2, 2))
gam.check(m_gam)
dev.off()

conc <- concurvity(m_gam, full = TRUE)
write.csv(as.data.frame(conc$estimate), "./out/tables/gam_concurvity_runs.csv")

cat("\n============================\n")
cat("Main GAM summary\n")
cat("============================\n")
cat("Adjusted R^2:", round(summary(m_gam)$r.sq, 3), "\n")
cat("Deviance explained:", round(summary(m_gam)$dev.expl, 3), "\n")
cat("n:", nobs(m_gam), "\n")
cat("\nSmooth terms:\n")
print(gam_smooths)
cat("\nModel comparison:\n")
print(model_comp)




# =========================================================
# 15. Heatmap of cost per hectare (uses summary_relative)
# =========================================================

# join data
costeff_df <- summary_relative %>%
  filter(sal_int > 0) %>%
  mutate(year = as.character(year)) %>%
  left_join(target_df_agg,
            by = c("year" = "year",
                   "strategy" = "strategy",
                   "sal_int" = "threshold")) %>%
  mutate(
    cost_eff = sum_SAL_costs / hectares_value,
    log_cost_eff = log(cost_eff))

hist(costeff_df$log_cost_eff)

# --- Prepare Data for Heatmap ---
heatmap_data <- costeff_df %>%
  filter(burned_type %in% c("burned_land", "burned_urbn")) %>%
  mutate(sal_int_label = factor(sal_int, levels = sort(unique(sal_int)))) %>%
  mutate(burned_type = case_when(burned_type == "burned_land" ~ "Burned Land (Total ha)",
                                 burned_type == "burned_urbn" ~ "Burned Urban area (ha)",
                                 TRUE ~ burned_type),
         time = factor(time),
         percentile_type = factor(percentile_type))

heatmap_data <- heatmap_data %>%
  filter(percentile_type %in% c("99th Percentile", "Median (50th Percentile)"),
         time %in% c(30, 60, 120, 300, 900))

rect_data_final <- data.frame(
  xmin = 5.5, xmax = 15.5,
  ymin = 10.5, ymax = 33.5,
  percentile_type = "99th Percentile",
  burned_type = "Burned Land (Total ha)"
)

# ----------------------------------
# 15a. Overview (big picture)
# ----------------------------------
p.heatmap.ce <- ggplot(
  heatmap_data %>%
    filter(percentile_type %in% c("99th Percentile", "Median (50th Percentile)"),
           time %in% c(30, 60, 120, 300, 900)),
  aes(
    x = interaction(factor(strategy), factor(sal_int_label), sep = " - "),
    y = interaction(as.factor(year), as.factor(time), sep = " - "),
    fill = cost_eff
  )
) +
  geom_tile(color = NA) +
  geom_rect(
    data = rect_data_final,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, color = "black", linewidth = 1,
    inherit.aes = FALSE
  ) +
  facet_grid(burned_type ~ percentile_type, scales = "free_y") +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis_c(
    option = "plasma", direction = -1,
    name = "Cost-effectiveness",
    trans = "log10",
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = c(1e3, 1e5, 1e7, 1e9)
  ) +
  labs(
    title = "A) Cost-effectiveness (General tendency)",
    x = "Scenario - Salvage logging (%)",
    y = NULL
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# ----------------------------------
# 15b. Detailed (zoomed) plot
# ----------------------------------
p.heatmap.detail <- ggplot(
  heatmap_data %>%
    filter(percentile_type %in% c("99th Percentile"),
           burned_type == "Burned Land (Total ha)",
           time %in% c(120, 300),
           sal_int %in% c(20, 30)),
  aes(
    x = factor(strategy),
    y = as.factor(year),
    fill = cost_eff
  )
) +
  geom_tile(color = NA) +
  facet_grid(as.factor(time) ~ factor(sal_int_label), scales = "free_y") +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(
    option = "plasma", direction = -1,
    name = "Cost-effectiveness",
    trans = "log10",
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = c(1e3, 1e5, 1e7, 1e9)
  ) +
  labs(
    title = "B) Cost-effectiveness (Detail)",
    x = "Scenario",
    y = "Year"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0, "pt"),
    panel.spacing.y = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# ----------------------------------
# 15c. Combine heatmap
# ----------------------------------
p.heatmap.comb <- ggarrange(
  p.heatmap.ce + theme(plot.margin = margin(5, 10, 5, 5)),
  p.heatmap.detail + theme(plot.margin = margin(5, 5, 5, 10)),
  ncol = 2,
  widths = c(1, 1),
  common.legend = TRUE,
  legend = "bottom"
)

p.heatmap.comb
ggsave(p.heatmap.comb, file = "./out/fig/p.heatmap.comb.png",
       bg = "white", width = 12, height = 8)

# ----------------------------------
# 15d. Quick OLS on heatmap data
# ----------------------------------
heatmap_data$log_cost_eff <- log(1 + heatmap_data$cost_eff)
heatmap_data$time_step <- as.integer(heatmap_data$time) / 5
heatmap_data$time2 <- heatmap_data$time_step * as.integer(heatmap_data$time)
heatmap_data$sal_int2 <- heatmap_data$sal_int * heatmap_data$sal_int
heatmap_data$strategy <- relevel(factor(heatmap_data$strategy), ref = "Random")
summary(lm(log_cost_eff ~ strategy +
             factor(year) +
             factor(percentile_type) +
             factor(burned_type) +
             time_step +
             time2 +
             sal_int +
             sal_int2,
           data = heatmap_data %>%
             filter(!log_cost_eff == Inf)))

