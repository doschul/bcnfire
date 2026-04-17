# ==============================================
# 06_bcn_mc_sim.R
# Description  : Monte Carlo simulation of avoided costs with WUI structure losses.
# Inputs       : data/rdat/perc_data_long.rds
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
  library(evd)
  library(extRemes)
  library(ggridges)
  library(scales)
  library(viridisLite)
  library(data.table)
  library(ggplot2)
})



# --------------------------
# Helper functions
# --------------------------
round_to_5 <- function(x) 5 * round(x / 5)
# function to return share of positive values
pos_share <- function(x) {
  sum(x > 0) / length(x)
}

suppress_console <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}
silent_random <- function(param, n) {
  suppressMessages(
    suppressWarnings(
      decisionSupport::random(param, n = n, verbosity = 0)
    )
  )
}

signed_log10 <- function(x) sign(x) * log10(1 + abs(x))
inv_signed_log10 <- function(y) sign(y) * (10^abs(y) - 1)

# midpoints for plotting
get_mid <- function(b) {
  sapply(strsplit(gsub("\\(|\\]|\\[", "", b), ","), function(v) mean(as.numeric(v)))
}

signed_sqrt_trans <- function() {
  trans_new(
    name = "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse   = function(x) sign(x) * (x^2),
    breaks    = pretty_breaks()
  )
}

# ----------------------------------------
# Helper: linear interpolation for areas (not needed, I simply round instead)
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

set.seed(123)
options(stringsAsFactors = FALSE)
# Project root assumed to be Git repo top; use here::here()
data_raw_dir  <- here::here("data", "raw")
data_proc_dir <- here::here("data", "rdat")
out_fig_dir   <- here::here("out", "fig")
out_tab_dir   <- here::here("out", "tab")

# Source helper functions
funs_path <- here::here("src", "R", "bcn_funs.R")
if (file.exists(funs_path)) source(funs_path)


# set GEV parameters manually based on the plot (or read from a file if you have them stored)

xi     <- -0.205   # replace with your ξ from the plot
mu     <- 94.9     # replace with your μ
sigma  <- 46.5    # replace with your σ

# Example: read from CSV exported in Python
#dur_params <- read_csv("calibrated_fire_area_distribution.csv")   # columns: mu_dur, sigma_dur
#mu_dur    <- dur_params$value[dur_params$parameter == "mu_log_eff"]
#sigma_dur <- dur_params$value[dur_params$parameter == "sigma_log_eff"]
#mu_dur    <- dur_params$mu_log[1]
#sigma_dur <- dur_params$sigma_log[1]

n <- 10000  # number of duration samples

dur <- rgev(n, loc = mu, scale = sigma, shape = xi)

# Ensure strictly positive durations (Python also truncated)
dur <- round_to_5(dur[dur > 5])  


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


# ---------- Analysis ----------

# -----------------------------
# Parameters and inputs
# -----------------------------
num_simulations <- 10000
v_strategies <- unique(perc_data_long$strategy)
v_years <- unique(perc_data_long$year)
v_thresholds <- unique(perc_data_long$sal_int) %>% subset(.>0)


setDT(target_df_agg)
setDT(perc_data_long)


# --------------------------
# Initialize result containers
# --------------------------
mc_raw_list_dt <- list()
mc_results_list_dt <- list()

# Scenario settings: full distribution and single extreme event
scenario_params <- data.table(
  fire_domain = c("Full", "Tail"),
  percentile_start = c(0.01, 0.95),
  percentile_end = c(0.99, 0.99),
  freq_min = c(1L, 1L),
  freq_max = c(10L, 1L),
  fixed_duration_min = c(NA_real_, 200)
)

# --------------------------
# Simulation loop
# --------------------------
for (i in seq_len(nrow(scenario_params))) {
  scen <- scenario_params[i]
  cat("\n=== Fire domain:", scen$fire_domain, "===\n")

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
        
        # --- 2. Sample number of fires per year for each simulation
        number_fires_vec <- sample(scen$freq_min:scen$freq_max, num_simulations, replace = TRUE)
        total_fires <- sum(number_fires_vec)
        
        # --- 3. Mapping of fires to simulation IDs
        fire_map <- data.table(
          simulation_id = rep(seq_len(num_simulations), times = number_fires_vec)
        )
        
        # --- 4. Draw all random values
        sampled_duration <- if (is.na(scen$fixed_duration_min)) {
          sample(dur, total_fires, replace = TRUE)
        } else {
          rep(scen$fixed_duration_min, total_fires)
        }

        sim_fire_dt <- data.table(
          simulation_id = fire_map$simulation_id,
          fire_id = sequence(number_fires_vec),
          percentile = runif(total_fires, scen$percentile_start, scen$percentile_end),
          duration_min = sampled_duration,
          
          land_price_agr = rlnorm(total_fires, meanlog = log(3000), sdlog = 0.6),
          fire_fight_cost_per_ha = rlnorm(total_fires, meanlog = log(24000), sdlog = 0.15),
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
                         sal_cost = current_sal_costs,
                         fire_domain = scen$fire_domain)]
        
        # --- 11. Summarize probability of positive net benefit
        nbp <- mean(yearly_dt$net_benefit > 0)
        
        mc_results_list_dt[[length(mc_results_list_dt) + 1]] <- data.table(
          strategy = current_strategy,
          year = current_year,
          sal_int = current_threshold,
          net_benefit_positive = nbp,
          sal_cost = current_sal_costs,
          fire_domain = scen$fire_domain
        )
        
        mc_raw_list_dt[[length(mc_raw_list_dt) + 1]] <- yearly_dt
      }
    }
  }
}

# ----------------------------------------
# Combine and save
# ----------------------------------------
mc_results_dt <- rbindlist(mc_results_list_dt)
mc_raw_dt <- rbindlist(mc_raw_list_dt)


saveRDS(mc_raw_dt, file = "mc_raw_dt.rds")
saveRDS(mc_results_dt, file = "mc_results_dt.rds")
# Optional backwards-compatible exports
#saveRDS(mc_results_dt[fire_domain == "Full"], file = "mc_results_dt_full.rds")
#saveRDS(mc_results_dt[fire_domain == "Tail"], file = "mc_results_dt_tail.rds")

# reload
mc_raw_dt <- readRDS("mc_raw_dt.rds")
mc_results_dt <- readRDS("mc_results_dt.rds")


# ----------------------------------------
# Plot (same as before)
# ----------------------------------------
library(ggplot2)
library(patchwork)
library(scales)
library(data.table)
library(tidyverse)

mc_results_dt <- mc_results_dt %>%
  mutate(
    fire_domain = factor(
      fire_domain,
      levels = c("Full", "Tail"),
      labels = c(
        "Panel A: Full distribution",
        "Panel B: Single long-duration high-severity fire"
      )
    )
  )



# 1. Define the shared color scale as a reusable function or object
# This ensures the 0.2 'yellow' transition is identical in both
shared_fill_scale <- scale_fill_gradientn(
  colors = c("firebrick", "khaki2", "darkgreen"),
  values = rescale(c(0, 0.5, 1)), 
  limits = c(0, 1),
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  trans = "sqrt"
)

# 2. Create Plot for Panel A
p1 <- ggplot(mc_results_dt[fire_domain == "Panel A: Full distribution"], 
             aes(x = sal_int, y = factor(year), fill = net_benefit_positive)) +
  geom_tile() +
  facet_wrap(~ strategy, nrow = 1) +
  shared_fill_scale +
  labs(title = "Panel A: Full distribution", x = NULL, y = "Year", fill = "P(Net Benefit > 0)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) # Hide x-axis for the top plot

# 3. Create Plot for Panel B
p2 <- ggplot(mc_results_dt[fire_domain == "Panel B: Single long-duration high-severity fire"], 
             aes(x = sal_int, y = factor(year), fill = net_benefit_positive)) +
  geom_tile() +
  facet_wrap(~ strategy, nrow = 1) +
  shared_fill_scale +
  labs(title = "Panel B: Single long-duration high-severity fire", x = "SAL Intensity", y = "Year", fill = "P(Net Benefit > 0)") +
  theme_minimal() +
  theme(strip.text = element_blank()) # Hide strategy labels if they are redundant with p1

# 4. Combine with Patchwork
# '/' stacks them vertically, 'guides = "collect"' merges the legends
p.heatmap.npv <- (p1 / p2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(plot = p.heatmap.npv, filename = "./out/fig/p_heatmap_npv.png", 
       width = 6, height = 4.5, bg = "white")




# make plot

# --- prepare data ---
# Ensure mc_raw_dt is a data.table or data.frame with columns: year, net_benefit
mc_raw_df <- as.data.frame(mc_raw_dt)
mc_raw_df <- mc_raw_df[is.finite(mc_raw_df$net_benefit), ]

build_cost_benefit_hist <- function(df_in, panel_title, show_x_axis = TRUE) {
  dt_plot <- rbindlist(list(
    data.table(type = "Benefit", value = df_in$benefit_total, sign = 1),
    data.table(type = "Cost", value = df_in$sal_cost, sign = -1)
  ))
  dt_plot[, sign := as.numeric(sign)]

  hdata <- dt_plot[, .(x = signed_log10(value)), by = .(type, sign)]
  hdata <- hdata[, .N, by = .(type, sign, bin = cut(x, breaks = 100))]
  hdata[, bin_mid := get_mid(bin)]
  hdata[, count_pos := N]

  P <- c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
  full_vals <- c(-rev(P), 0, P)

  p <- ggplot(hdata, aes(x = bin_mid, y = count_pos, fill = type)) +
    geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.5) +
    geom_col(data = hdata[type == "Benefit"], alpha = 0.6, color = "grey40") +
    geom_col(data = hdata[type == "Cost"], aes(y = -count_pos), alpha = 0.6, color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Benefit" = "darkgreen", "Cost" = "firebrick")) +
    scale_x_continuous(
      name = "Value (EUR)",
      breaks = signed_log10(full_vals),
      labels = comma(full_vals)
    ) +
    scale_y_continuous(
      trans = "signed_sqrt",
      name = "Count",
      labels = label_number(accuracy = 1, big.mark = ",", style_negative = "parens"),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(title = panel_title, fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }

  p
}

p.mc.full <- build_cost_benefit_hist(
  df_in = mc_raw_df[mc_raw_df$fire_domain == "Full", ],
  panel_title = "Panel A: Full distribution",
  show_x_axis = FALSE
)

p.mc.tail <- build_cost_benefit_hist(
  df_in = mc_raw_df[mc_raw_df$fire_domain == "Tail", ],
  panel_title = "Panel B: Single long-duration high-severity fire",
  show_x_axis = TRUE
)

p.mc.combined <- (p.mc.full / p.mc.tail) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

p.mc.combined
ggsave(plot = p.mc.combined, "./out/fig/p_mc_cost_benefit_full_vs_extreme.png",
       width = 7, height = 8, bg = "white")



