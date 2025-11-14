# ==============================================
# 06_bcn_mc_sim.R
# Description  : Monte Carlo simulation of avoided costs with WUI structure losses.
# Inputs       : data/rdat/perc_data_long.rds, data/rdat/duration_params.rds
# Outputs      : data/rdat/mc_raw_dt.rds, data/rdat/mc_results_dt.rds, out/figures/net_benefit_distribution.png
# Author       : (auto-refactor)
# Created      : (auto)
# ==============================================

rm(list=ls())
setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

# ---------- Setup ----------
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(sf)
  library(terra)
  library(ggplot2)
  library(stringr)
  library(here)
  #library(decisionSupport)
  library(extRemes)
})

# --------------------------
# Helper: round duration to 5-min intervals
# --------------------------
round_to_5 <- function(x) 5 * round(x / 5)


set.seed(123)
options(stringsAsFactors = FALSE)
# Project root assumed to be Git repo top; use here::here()
data_raw_dir  <- here::here("data", "raw")
data_proc_dir <- here::here("data", "rdat")
out_fig_dir   <- here::here("out", "fig")
out_tab_dir   <- here::here("out", "tab")

# Source helper functions
funs_path <- here::here("src", "bcn_funs.R")
if (file.exists(funs_path)) source(funs_path)


# Example: read from CSV exported in Python
#dur_params <- read_csv("calibrated_fire_area_distribution.csv")   # columns: mu_dur, sigma_dur
#mu_dur    <- dur_params$value[dur_params$parameter == "mu_log_eff"]
#sigma_dur <- dur_params$value[dur_params$parameter == "sigma_log_eff"]
#mu_dur    <- dur_params$mu_log[1]
#sigma_dur <- dur_params$sigma_log[1]

# GEV parameters
library(evd)

xi     <- -0.205   # replace with your ξ from the plot
mu     <- 94.9     # replace with your μ
sigma  <- 46.5    # replace with your σ

n <- 10000  # number of duration samples

dur <- rgev(n, loc = mu, scale = sigma, shape = xi)

# Ensure strictly positive durations (Python also truncated)
dur <- round_to_5(dur[dur > 5])  



# TEST WUI
#perc_data_long <- readRDS(file = "./data/rdat/perc_data_long_WUI_TEST.rds")

perc_data_long <- readRDS(file = "./data/rdat/perc_data_long.rds")
perc_data_long$time <- round_to_5(perc_data_long$time)
perc_data_long$percentile <- round(perc_data_long$percentile, 2)

#summary_relative <- read_rds("./data/rdat/summary_relative.rds")
target_df_agg <- readRDS(file = "./data/rdat/target_df_agg.rds")


target_df_agg <- target_df_agg %>%
  mutate(strategy = case_when(
    strategy == "Highest artificial" ~ "High Built-up",
    strategy == "Highest biomass" ~ "High Biomass",
    strategy == "Highest connectivity" ~ "High Connectivity",
    strategy == "Low logging costs" ~ "Low SAL cost",
    TRUE ~ strategy
  ))


# function to return share of positive values
pos_share <- function(x) {
  sum(x > 0) / length(x)
}

suppress_console <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}




# ----------------------------------------
# Helper: linear interpolation for areas
# ----------------------------------------
# interp_area_dt <- function(perc_dt, sampled_duration) {
#   setorder(perc_dt, time)
#   lo_idx <- findInterval(sampled_duration, perc_dt$time)
#   lo_idx <- pmax(lo_idx, 1)
#   hi_idx <- pmin(lo_idx + 1, nrow(perc_dt))
#   t0 <- perc_dt$time[lo_idx]
#   t1 <- perc_dt$time[hi_idx]
#   w  <- (sampled_duration - t0) / (t1 - t0)
#   
#   out <- perc_dt[lo_idx, .(abs_avoid_agri_ha, abs_avoid_urbn_ha,
#                            abs_avoid_wild_ha, abs_avoid_tot_ha)] * (1 - w) +
#     perc_dt[hi_idx, .(abs_avoid_agri_ha, abs_avoid_urbn_ha,
#                       abs_avoid_wild_ha, abs_avoid_tot_ha)] * w
#   return(out)
# }

silent_random <- function(param, n) {
  suppressMessages(
    suppressWarnings(
      decisionSupport::random(param, n = n, verbosity = 0)
    )
  )
}

# ---------- Analysis ----------

# -----------------------------
# Parameters and inputs
# -----------------------------
num_simulations <- 5000
v_strategies <- unique(perc_data_long$strategy)
v_years <- unique(perc_data_long$year)
v_thresholds <- unique(perc_data_long$sal_int) %>% subset(.>0)

# Cost and value distributions
# rho_land_price_agr <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
#                            quantiles = c(4000, 8000))
# rho_houses_per_ha  <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
#                            quantiles = c(1, 20))
# rho_house_price    <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
#                            quantiles = c(10000, 100000))
# rho_fire_fight_cost_per_ha <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
#                                    quantiles = c(500, 5000))

setDT(target_df_agg)
setDT(perc_data_long)


# --------------------------
# Initialize result containers
# --------------------------
mc_raw_list_dt <- list()
mc_results_list_dt <- list()

# --------------------------
# Simulation loop
# --------------------------
for (s in seq_along(v_strategies)) {
  current_strategy <- v_strategies[s]
  cat("\nProcessing strategy:", current_strategy)
  
  target_strategy_dt <- target_df_agg[strategy == current_strategy]
  area_strategy_dt   <- perc_data_long[strategy == current_strategy]
  
  for (y in seq_along(v_years)) {
    current_year <- v_years[y]
    cat("\n - Year:", current_year)
    
    sal_costs_dt <- target_strategy_dt[year == current_year,
                                       .(sal_int = threshold, sum_SAL_costs)]
    area_year_dt <- area_strategy_dt[year == current_year]
    
    cat(" -- Threshold: ")
    
    for (t in seq_along(v_thresholds)) {
      current_threshold <- v_thresholds[t]
      cat(current_threshold, "..")
      
      # --- 1. Extract costs and areas
      current_sal_costs <- sal_costs_dt[sal_int == current_threshold, sum_SAL_costs]
      if (length(current_sal_costs) == 0) current_sal_costs <- 0
      
      # Subset area data for this threshold
      current_area_dt <- area_year_dt[sal_int == current_threshold,
                                      .(time, percentile,
                                        abs_avoid_tot_ha,
                                        abs_avoid_agri_ha,
                                        abs_avoid_urbn_ha,
                                        abs_avoid_wild_ha)]
      
      # --- 2. Sample number of fires per year for each simulation (1–10)
      number_fires_vec <- sample(1:10, num_simulations, replace = TRUE)
      total_fires <- sum(number_fires_vec)
      
      # --- 3. Mapping of fires to simulation IDs
      fire_map <- data.table(
        simulation_id = rep(seq_len(num_simulations), times = number_fires_vec)
      )
      
      # --- 4. Draw all random values
      sim_fire_dt <- data.table(
        simulation_id = fire_map$simulation_id,
        fire_id = sequence(number_fires_vec),
        #percentile = runif(total_fires, 0.01, 0.99),
        percentile = runif(total_fires, 0.95, 0.99),
        #duration_min = sample(dur, total_fires, replace = TRUE),
        duration_min = rep(200, total_fires), # 95 percentile of calibrated durations
        
        land_price_agr = rlnorm(total_fires, meanlog = log(6000), sdlog = 0.3),
        fire_fight_cost_per_ha = rlnorm(total_fires, meanlog = log(15000), sdlog = 0.5),
        eur_m2 = rlnorm(total_fires, meanlog = log(1600), sdlog = 0.25),
        damage_ratio = pmin(rbeta(total_fires, 2, 4), 1),
        
        dens_verylow  = runif(total_fires, 0.02, 0.06),
        dens_low      = runif(total_fires, 0.10, 0.50),
        dens_medium   = runif(total_fires, 0.50, 40.0),
        dens_intermix = runif(total_fires, 0.06, 7.50),
        dens_interface= runif(total_fires, 0.06, 40.0),
        
        fa_verylow  = rnorm(total_fires, 250, 20),
        fa_low      = rnorm(total_fires, 160, 20),
        fa_medium   = rnorm(total_fires, 100, 15),
        fa_intermix = rnorm(total_fires, 150, 15),
        fa_interface= rnorm(total_fires, 120, 15),
        
        vuln_verylow  = 0.6,
        vuln_low      = 0.7,
        vuln_medium   = 0.85,
        vuln_intermix = 1.0,
        vuln_interface= 1.0
      )
      
      # --- 5. Discretize to lookup resolution
      sim_fire_dt[, percentile := round(percentile, 2)]
      sim_fire_dt[, duration_min := round_to_5(duration_min)]
      
      
      # --- 6. Exact match join to attach avoided areas
      setkey(current_area_dt, percentile, time)
      setkey(sim_fire_dt, percentile, duration_min)
      
      sim_fire_dt <- current_area_dt[
        sim_fire_dt, on = .(percentile, time = duration_min), roll = "nearest"
      ]

      #check_join_dimensions(sim_fire_dt, current_area_dt, total_fires)
      
      # --- 7. Attach WUI avoided areas (static per threshold)
      this_wui <- area_year_dt[sal_int == current_threshold,
                               .(abs_avoid_wui_nohouse_ha,
                                 abs_avoid_wui_verylow_ha,
                                 abs_avoid_wui_low_ha,
                                 abs_avoid_wui_medium_ha,
                                 abs_avoid_wui_intermix_ha,
                                 abs_avoid_wui_interface_ha)][1]
      sim_fire_dt[, (names(this_wui)) := this_wui]
      
      # --- 8. Compute benefits
      sim_fire_dt[, benefit_agr := abs_avoid_agri_ha * land_price_agr]
      sim_fire_dt[, benefit_ff  := abs_avoid_tot_ha  * fire_fight_cost_per_ha]
      
      sim_fire_dt[, struct_loss :=
                    abs_avoid_wui_verylow_ha  * dens_verylow  * fa_verylow  * eur_m2 * damage_ratio * vuln_verylow  +
                    abs_avoid_wui_low_ha      * dens_low      * fa_low      * eur_m2 * damage_ratio * vuln_low      +
                    abs_avoid_wui_medium_ha   * dens_medium   * fa_medium   * eur_m2 * damage_ratio * vuln_medium   +
                    abs_avoid_wui_intermix_ha * dens_intermix * fa_intermix * eur_m2 * damage_ratio * vuln_intermix +
                    abs_avoid_wui_interface_ha* dens_interface* fa_interface* eur_m2 * damage_ratio * vuln_interface]
      
      sim_fire_dt[, benefit_per_event := benefit_agr + benefit_ff + struct_loss]
      
      # --- 9. Aggregate per year
      yearly_dt <- sim_fire_dt[, .(
        benefit_total = sum(benefit_per_event),
        number_fires_per_year = .N
      ), by = simulation_id]
      
      # --- 10. Compute net benefit
      yearly_dt[, net_benefit := benefit_total - current_sal_costs]
      yearly_dt[, `:=`(strategy = current_strategy,
                       year = current_year,
                       sal_int = current_threshold,
                       sal_cost = current_sal_costs)]
      
      # --- 11. Summarize probability of positive net benefit
      nbp <- mean(yearly_dt$net_benefit > 0)
      
      mc_results_list_dt[[length(mc_results_list_dt) + 1]] <- data.table(
        strategy = current_strategy,
        year = current_year,
        sal_int = current_threshold,
        net_benefit_positive = nbp,
        sal_cost = current_sal_costs
      )
      
      mc_raw_list_dt[[length(mc_raw_list_dt) + 1]] <- yearly_dt
    }
  }
}

# ----------------------------------------
# Combine and save
# ----------------------------------------
mc_results_dt <- rbindlist(mc_results_list_dt)
mc_raw_dt <- rbindlist(mc_raw_list_dt)


summary(mc_raw_dt$benefit_total)

saveRDS(mc_raw_dt, file = "mc_raw_dt.rds")




# ----------------------------------------
# Plot (same as before)
# ----------------------------------------
library(ggplot2)
p.heatmap.npv <- ggplot(mc_results_dt, aes(x = sal_int, y = factor(year), fill = net_benefit_positive)) +
  geom_tile() +
  facet_wrap(~ strategy) +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", midpoint = 0.5) +
  labs(
    x = "SAL Intensity", y = "Year", fill = "P(Net Benefit > 0)"
  ) +
  theme_minimal()
p.heatmap.npv


ggsave(plot = p.heatmap.npv, filename = "./out/fig/p_heatmap_npv_95tail.png", 
       width = 6, height = 3.5, bg = "white")




# make ridgeplot
library(ggridges)
library(scales)
library(viridisLite)
library(data.table)
library(ggplot2)
library(scales)
# Ensure mc_raw_dt is a data.table or data.frame with columns: year, net_benefit
mc_raw_df <- as.data.frame(mc_raw_dt)
mc_raw_df <- mc_raw_df[is.finite(mc_raw_df$net_benefit), ]

# --- helper ---
signed_log10 <- function(x) sign(x) * log10(1 + abs(x))

# --- prepare data ---
dt_plot <- rbindlist(list(
  data.table(type = "Benefit", value = mc_raw_df$benefit_total, sign = 1),
  data.table(type = "Cost", value = mc_raw_df$sal_cost, sign = -1)
))

dt_plot[, sign := as.numeric(sign)]

# --- pre-bin (so we can control Y signs and axis) ---
hdata <- dt_plot[, .(x = signed_log10(value)), by = .(type, sign)]
hdata <- hdata[, .N, by = .(type, sign, bin = cut(x, breaks = 100))]
# midpoints for plotting
get_mid <- function(b) {
  sapply(strsplit(gsub("\\(|\\]|\\[", "", b), ","), function(v) mean(as.numeric(v)))
}
hdata[, bin_mid := get_mid(bin)]
hdata[, count_pos := N] # always positive


# Load the scales package for custom transformations
library(scales) 

# --- Define break values and corresponding labels ---
# P: positive break values (in original currency scale)
P <- c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)

# Full set of values to label, including negative and zero
full_vals <- c(-rev(P), 0, P)

# Custom function is preserved from your snippet
signed_log10 <- function(x) sign(x) * log10(1 + abs(x))

library(scales)

signed_sqrt_trans <- function() {
  trans_new(
    name = "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse   = function(x) sign(x) * (x^2),
    breaks    = pretty_breaks()
  )
}


# --- plotting ---
p.mc.1 <- ggplot(hdata, aes(x = bin_mid, y = count_pos, fill = type)) +
  # Vertical line at X=0 (which is the median/zero point)
  geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.5) +
  geom_col(data = hdata[type == "Benefit"], alpha = 0.6, color = "grey40") +
  geom_col(data = hdata[type == "Cost"], aes(y = -count_pos), alpha = 0.6, color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Benefit" = "darkgreen", "Cost" = "firebrick")) +
  
  # --- Custom Axis Formatting ---
  scale_x_continuous(
    name = "Value (EUR)",
    # Apply signed_log10 to all break values (negative, zero, and positive)
    breaks = signed_log10(full_vals), 
    # Use scales::comma to format all labels
    labels = comma(full_vals) 
  ) +
  scale_y_continuous(
    # Keeping your original Y-axis definition, but note that 
    trans = "signed_sqrt", # would be the better choice if you need scaling.
    name = "Count",
    labels = label_number(accuracy = 1, big.mark = ",", style_negative = "parens"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(title = "Monte Carlo output: Benefits (up) vs Costs (down)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

p.mc.1
ggsave(plot = p.mc.1, "./out/fig/p_mc_cost_benefit_95tail.png", 
       width = 6, height = 4, bg = "white")


summary(lm(net_benefit~percentile + duration_min + number_fires_per_year + land_price_agr + houses_per_ha + house_price + strategy + sal_int, data = mc_raw_df))

mc_raw_df$net_benefit_log <- signed_log10(mc_raw_df$net_benefit)

# --- 2. ridge plot ---

# Assuming the original helper function is defined as:
# signed_log10 <- function(x) sign(x) * log10(1 + abs(x))

# ------------------------
# 1. Define the Inverse Transformation Function (Required for Labels)
# ------------------------

# This function maps the axis position (x) back to the original currency value.
# It is the inverse of the signed_log10 function.
signed_antilog10 <- function(x_transformed) {
  sign(x_transformed) * (10^(abs(x_transformed)) - 1)
}

# ------------------------
# 2. Ridge Plot with Corrected X-Axis Labels
# ------------------------
mc_raw_df$net_benefit_log <- signed_log10(mc_raw_df$net_benefit)

p.mc.2 <- ggplot(mc_raw_df %>%
         filter(!is.na(net_benefit_log)), # Simplified filter for clarity
       aes(x = net_benefit_log, y = factor(year), fill = ..density..)) +
  geom_density_ridges_gradient(
    scale = 2,
    rel_min_height = 0.01,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(name = "Density", option = "C") +
  
  # --- Corrected X-Axis Ticks ---
  scale_x_continuous(
    name = "Net Benefit (€)",
    # Define a clean set of break locations on the transformed axis
    breaks = seq(-8, 8, by = 2), # Adjusted breaks to cover common log10 range
    
    # Use the inverse function to calculate the true value for each break location (x)
    labels = function(x) {
      # 1. Back-transform x (the axis break location) to the original currency value (val)
      val <- signed_antilog10(x) 
      
      # 2. Format the true currency value (val) with suffixes
      sapply(val, function(v) {
        # Check for NA/Inf/non-finite values
        if (!is.finite(v) || abs(v) < 1) return("0") 
        
        # Use abs_v for formatting to handle negative signs correctly later
        abs_v <- abs(v)
        
        formatted_val <- if (abs_v < 1e3) {
          # < 1k
          sprintf("%.0f", v)
        } else if (abs_v < 1e6) {
          # 1k to 1M
          # Use scales::label_dollar to format with commas and 'k' suffix for cleanliness
          scales::label_number(scale = 1/1e3, suffix = "K", accuracy = 1)(abs_v)
        } else if (abs_v < 1e9) {
          # 1M to 1B
          scales::label_number(scale = 1/1e6, suffix = "M", accuracy = 1)(abs_v)
        } else {
          # 1B+
          scales::label_number(scale = 1/1e9, suffix = "B", accuracy = 1)(abs_v)
        }
        
        # Prepend the sign if the original value was negative, otherwise return as is
        if (v < 0) paste0("-", formatted_val) else formatted_val
      })
    }
  ) +
  
  labs(
    title = "Distribution of Simulated Net Benefits per Year",
    subtitle = "Monte Carlo results, symmetric log-transformed scale (back-transformed axis)",
    y = "Year"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(plot = p.mc.2, "./out/fig/p_mc_net_benefit_95tail.png", 
       height = 4, width = 6, bg = "white")



# --- 3. define facet variable (Year 2049 vs. other years) ---
mc_raw_df <- mc_raw_df %>%
  mutate(year_group = ifelse(year == 2049, "Year 2049", "Other years"))

# --- 4. Ridge plot by scenario (robust version) ---
ggplot(mc_raw_df,
       aes(x = net_benefit_log, y = strategy, fill = ..density..)) +
  geom_density_ridges_gradient(
    scale = 2,
    rel_min_height = 0.01,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(name = "Density", option = "C") +
  
  # --- Back-transform x-axis labels ---
  scale_x_continuous(
    name = "Net Benefit (€)",
    labels = function(x) {
      val <- inv_signed_log10(x)
      sapply(val, function(v) {
        if (is.na(v)) return(NA_character_)
        if (abs(v) < 1e3) {
          sprintf("%.0f", v)
        } else if (abs(v) < 1e6) {
          sprintf("%s%sk",
                  ifelse(v < 0, "-", ""),
                  formatC(abs(v) / 1e3, format = "f", digits = 0))
        } else {
          sprintf("%s%.1fM",
                  ifelse(v < 0, "-", ""),
                  abs(v) / 1e6)
        }
      })
    },
    breaks = seq(-10, 10, by = 2)
  ) +
  
  facet_wrap(~ year_group, nrow = 2, scales = "free_y") +
  labs(
    title = "Distribution of Simulated Net Benefits by Scenario",
    subtitle = "Comparison between Year 2049 (outlier) and all other years",
    y = "Management Scenario"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )






