# MC Sim with calibrated durations

rm(list = ls())

library(tidyverse)

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")


# Example: read from CSV exported in Python
dur_params <- read_csv("posterior_duration_params.csv")   # columns: mu_dur, sigma_dur
mu_dur    <- dur_params$mu_dur
sigma_dur <- dur_params$sigma_dur


# TEST WUI
perc_data_long <- readRDS(file = "./data/rdat/perc_data_long_WUI_TEST.rds")

perc_data_long <- readRDS(file = "./data/rdat/perc_data_long.rds")
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
interp_area_dt <- function(perc_dt, sampled_duration) {
  setorder(perc_dt, time)
  lo_idx <- findInterval(sampled_duration, perc_dt$time)
  lo_idx <- pmax(lo_idx, 1)
  hi_idx <- pmin(lo_idx + 1, nrow(perc_dt))
  t0 <- perc_dt$time[lo_idx]
  t1 <- perc_dt$time[hi_idx]
  w  <- (sampled_duration - t0) / (t1 - t0)
  
  out <- perc_dt[lo_idx, .(abs_avoid_agri_ha, abs_avoid_urbn_ha,
                           abs_avoid_wild_ha, abs_avoid_tot_ha)] * (1 - w) +
    perc_dt[hi_idx, .(abs_avoid_agri_ha, abs_avoid_urbn_ha,
                      abs_avoid_wild_ha, abs_avoid_tot_ha)] * w
  return(out)
}

silent_random <- function(param, n) {
  suppressMessages(
    suppressWarnings(
      decisionSupport::random(param, n = n, verbosity = 0)
    )
  )
}

# ----------------------------------------
# Setup
# ----------------------------------------
library(data.table)
library(decisionSupport)
library(extRemes)

# -----------------------------
# Parameters and inputs
# -----------------------------
num_simulations <- 2000
v_strategies <- unique(perc_data_long$strategy)
v_years <- unique(perc_data_long$year)
v_thresholds <- unique(perc_data_long$sal_int)

# Cost and value distributions
rho_land_price_agr <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
                           quantiles = c(4000, 8000))
rho_houses_per_ha  <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
                           quantiles = c(1, 20))
rho_house_price    <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
                           quantiles = c(10000, 100000))
rho_fire_fight_cost_per_ha <- list(distribution = "lnorm", probabilities = c(0.05, 0.95),
                                   quantiles = c(500, 5000))

setDT(target_df_agg)
setDT(perc_data_long)




# ----------------------------------------
# Initialize results
# ----------------------------------------
mc_raw_list_dt <- list()
mc_results_list_dt <- list()

# ----------------------------------------
# Simulation loop (vectorized version)
# ----------------------------------------

for (s in seq_along(v_strategies)) {
  current_strategy <- v_strategies[s]
  cat("\n", "Processing strategy:", current_strategy, "\n")
  
  target_strategy_dt <- target_df_agg[strategy == current_strategy]
  area_strategy_dt   <- perc_data_long[strategy == current_strategy]
  
  for (y in seq_along(v_years)) {
    current_year <- v_years[y]
    cat("- Year:", current_year, "\n")
    
    sal_costs_dt <- target_strategy_dt[year == current_year,
                                       .(sal_int = threshold, sum_SAL_costs)]
    area_year_dt <- area_strategy_dt[year == current_year]
    
    cat("-- Threshold: ")
    
    for (t in seq_along(v_thresholds)) {
      current_threshold <- v_thresholds[t]
      cat(current_threshold, "..")
      
      # --- 1. Extract costs and areas
      current_sal_costs <- sal_costs_dt[sal_int == current_threshold, sum_SAL_costs]
      if (length(current_sal_costs) == 0) current_sal_costs <- 0
      
      current_area_dt <- area_year_dt[sal_int == current_threshold,
                                      .(time, percentile, abs_avoid_tot_ha,
                                        abs_avoid_agri_ha, abs_avoid_urbn_ha,
                                        abs_avoid_wild_ha)]
      
      # --- 2. Sample number of fires per year for each simulation (1–10)
      number_fires_vec <- sample(1:10, num_simulations, replace = TRUE)
      total_fires <- sum(number_fires_vec)
      
      # --- 3. Mapping of fires to simulation IDs
      fire_map <- data.table(
        simulation_id = rep(seq_len(num_simulations), times = number_fires_vec)
      )
      
      # --- 4. Draw all random values at once (vectorized)
      sim_fire_dt <- data.table(
        simulation_id = fire_map$simulation_id,
        fire_id = sequence(number_fires_vec),
        percentile = runif(total_fires, 0.01, 0.99),
        duration_min = rlnorm(total_fires, meanlog = mu_dur, sdlog = sigma_dur),
        
        # farmland + firefighting have not changed
        land_price_agr = rlnorm(total_fires, meanlog = log(6000), sdlog = 0.3),
        fire_fight_cost_per_ha = rlnorm(total_fires, meanlog = log(1500), sdlog = 0.5),
        
        # structure cost parameters (€/m2 replacement cost)
        eur_m2 = rlnorm(total_fires, meanlog = log(1100), sdlog = 0.25),
        
        # damage ratio (fraction of structure value actually destroyed)
        damage_ratio = pmin(rbeta(total_fires, 2, 8), 0.7),
        
        # dwelling densities per ha (random by class)
        dens_verylow  = runif(total_fires, 0.03, 0.06),
        dens_low      = runif(total_fires, 0.10, 0.50),
        dens_medium   = runif(total_fires, 0.50, 0.90),
        dens_intermix = runif(total_fires, 0.60, 1.20),
        dens_interface= runif(total_fires, 0.40, 0.90),
        
        # dwelling floor areas (m2)
        fa_verylow  = rnorm(total_fires, 150, 20),
        fa_low      = rnorm(total_fires, 120, 20),
        fa_medium   = rnorm(total_fires, 100, 15),
        fa_intermix = rnorm(total_fires, 90, 15),
        fa_interface= rnorm(total_fires, 100, 20),
        
        # vulnerability multipliers
        vuln_verylow  = 0.5,
        vuln_low      = 0.7,
        vuln_medium   = 0.85,
        vuln_intermix = 1.0,
        vuln_interface= 0.9
      )
      
      # --- 5. Interpolate avoided TOTAL areas (as before)
      interp_areas <- interp_area_dt(current_area_dt, sim_fire_dt$duration_min)
      sim_fire_dt[, c("abs_avoid_agri_ha",
                      "abs_avoid_urbn_dummy",
                      "abs_avoid_wild_ha",
                      "abs_avoid_tot_ha") := interp_areas]
      
      # --- 5b. Attach WUI avoided areas for this strategy-year-threshold
      # (static slice repeated per fire)
      this_wui <- area_year_dt[sal_int == current_threshold,
                               .(abs_avoid_wui_nohouse_ha,
                                 abs_avoid_wui_verylow_ha,
                                 abs_avoid_wui_low_ha,
                                 abs_avoid_wui_medium_ha,
                                 abs_avoid_wui_intermix_ha,
                                 abs_avoid_wui_interface_ha)
      ][1]
      sim_fire_dt[, (names(this_wui)) := this_wui]
      
      # --- 6. Compute agricultural & firefighting benefits (unchanged)
      sim_fire_dt[, benefit_agr := abs_avoid_agri_ha * land_price_agr]
      sim_fire_dt[, benefit_ff  := abs_avoid_tot_ha  * fire_fight_cost_per_ha]
      
      # --- 6b. Compute structure loss avoided per WUI class
      sim_fire_dt[, struct_loss :=
                    # no-house area contributes zero
                    abs_avoid_wui_verylow_ha  * dens_verylow  * fa_verylow  * eur_m2 * damage_ratio * vuln_verylow  +
                    abs_avoid_wui_low_ha      * dens_low      * fa_low      * eur_m2 * damage_ratio * vuln_low      +
                    abs_avoid_wui_medium_ha   * dens_medium   * fa_medium   * eur_m2 * damage_ratio * vuln_medium   +
                    abs_avoid_wui_intermix_ha * dens_intermix * fa_intermix * eur_m2 * damage_ratio * vuln_intermix +
                    abs_avoid_wui_interface_ha* dens_interface* fa_interface* eur_m2 * damage_ratio * vuln_interface
      ]
      
      # --- 6c. Total benefit per fire
      sim_fire_dt[, benefit_per_event := benefit_agr + benefit_ff + struct_loss]
      
      # --- 7. Aggregate per year (2000 sims)
      yearly_dt <- sim_fire_dt[, .(
        benefit_total = sum(benefit_per_event),
        number_fires_per_year = .N
      ), by = simulation_id]
      
      
      # --- 8. Compute net benefits (subtract cost once per scenario)
      yearly_dt[, net_benefit := benefit_total - current_sal_costs]
      yearly_dt[, `:=`(strategy = current_strategy,
                       year = current_year,
                       sal_int = current_threshold)]
      
      # --- 9. Summarize probability of positive net benefit
      nbp <- mean(yearly_dt$net_benefit > 0)
      
      mc_results_list_dt[[length(mc_results_list_dt) + 1]] <- data.table(
        strategy = current_strategy,
        year = current_year,
        sal_int = current_threshold,
        net_benefit_positive = nbp
      )
      
      # store raw simulations
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
quantile(mc_raw_dt$benefit_total, c(0.5, 0.9, 0.99, 0.999))



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


ggsave(plot = p.heatmap.npv, filename = "./out/fig/p_heatmap_npv.png", 
       width = 6, height = 3.5, bg = "white")




# make ridgeplot
library(ggridges)
library(scales)
library(viridisLite)

# Ensure mc_raw_dt is a data.table or data.frame with columns: year, net_benefit
mc_raw_df <- as.data.frame(mc_raw_dt)
mc_raw_df <- mc_raw_df[is.finite(mc_raw_df$net_benefit), ]


# Symmetric pseudo-log transform (handles negative + positive values)
signed_log10 <- function(x) sign(x) * log10(1 + abs(x))
inv_signed_log10 <- function(y) sign(y) * (10^abs(y) - 1)
mc_raw_df$net_benefit_log <- signed_log10(mc_raw_df$net_benefit)

# Check benefit distribution before subtracting costs
hist(signed_log10(mc_raw_dt$benefit_total), breaks=100, main="Benefit_total before cost")


summary(lm(net_benefit~percentile + duration_min + number_fires_per_year + land_price_agr + houses_per_ha + house_price + strategy + sal_int, data = mc_raw_df))



# --- 2. ridge plot ---
ggplot(mc_raw_df %>%
         filter(sal_int  == 50), aes(x = net_benefit_log, y = factor(year), fill = ..density..)) +
  geom_density_ridges_gradient(
    scale = 2,
    rel_min_height = 0.01,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(name = "Density", option = "C") +
  facet_wrap(~ as.factor(strategy), scales = "free_y") +
  #--- 3. back-transform the x-axis ticks to show actual € values ---
  scale_x_continuous(
    name = "Net Benefit (€)",
    labels = function(x) {
      val <- inv_signed_log10(x)
      # format large values nicely, showing negatives correctly
      sapply(val, function(v) {
        if (abs(v) < 1e3) return(sprintf("%.0f", v))
        if (abs(v) < 1e6) return(sprintf("%s%sk",
                                         ifelse(v < 0, "-", ""), formatC(abs(v) / 1e3, format = "f", digits = 0)))
        return(sprintf("%s%.1fM",
                       ifelse(v < 0, "-", ""), abs(v) / 1e6))
      })
    },
    breaks = seq(-10, 10, by = 2)
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
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # <— rotate labels
  )


# --- 3. define facet variable (Year 2049 vs. other years) ---
mc_raw_df <- mc_raw_df %>%
  mutate(year_group = ifelse(year == 2049, "Year 2049", "Other years"))

# --- 4. ridge plot by scenario ---
ggplot(mc_raw_df, aes(x = net_benefit_log, y = strategy, fill = ..density..)) +
  geom_density_ridges_gradient(
    scale = 2,
    rel_min_height = 0.01,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(name = "Density", option = "C") +
  
  # back-transform x-axis labels
  scale_x_continuous(
    name = "Net Benefit (€)",
    labels = function(x) {
      val <- inv_signed_log10(x)
      sapply(val, function(v) {
        if (abs(v) < 1e3) return(sprintf("%.0f", v))
        if (abs(v) < 1e6) return(sprintf("%s%sk",
                                         ifelse(v < 0, "-", ""), formatC(abs(v) / 1e3, format = "f", digits = 0)))
        return(sprintf("%s%.1fM",
                       ifelse(v < 0, "-", ""), abs(v) / 1e6))
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


