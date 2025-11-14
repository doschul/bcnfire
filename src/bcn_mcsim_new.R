# make some test data

rm(list = ls())

library(tidyverse)
setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

perc_data_long <- readRDS(file = "./data/rdat/perc_data_long.rds")
summary_relative <- read_rds("./data/rdat/summary_relative.rds")
target_df_agg <- readRDS(file = "./data/rdat/target_df_agg.rds")


perc_data_long$percentile <- as.numeric(round(perc_data_long$percentile, 2))

target_df_agg <- target_df_agg %>%
  mutate(strategy = case_when(
    strategy == "Highest artificial" ~ "High Built-up",
    strategy == "Highest biomass" ~ "High Biomass",
    strategy == "Highest connectivity" ~ "High Connectivity",
    strategy == "Low logging costs" ~ "Low SAL cost",
    TRUE ~ strategy
  ))

# filter out BAU
summary_relative <- summary_relative %>%
  filter(sal_int > 0)


# iterate over different scenarios, different fixed costs

# function to return share of positive values
pos_share <- function(x) {
  sum(x > 0) / length(x)
}

suppress_console <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}

selected_time <- 180
num_simulations <- 2000
number_fires_per_year <- 5


v_strategies <- unique(summary_relative$strategy)
v_years <- unique(summary_relative$year)
v_thresholds <- unique(summary_relative$sal_int)


rho_land_price_agr <- list(distribution = "lnorm", 
                           probabilities = c(0.05, 0.95), 
                           quantiles = c(4000, 8000))

rho_houses_per_ha <- list(distribution = "lnorm", 
                          probabilities = c(0.05, 0.95), 
                          quantiles = c(1, 20))

rho_house_price <- list(distribution = "lnorm", 
                        probabilities = c(0.05, 0.95), 
                        quantiles = c(50000, 500000))

rho_fire_fight_cost_per_ha <- list(distribution = "lnorm", 
                                   probabilities = c(0.05, 0.95), 
                                   quantiles = c(500, 5000))




library(data.table)

setDT(target_df_agg)
setDT(perc_data_long)

# Create sim_dt directly as a data.table
sim_dt <- data.table(
  simulation_id = 1:num_simulations,
  percentile = as.numeric(round(sample(seq(0.01, 0.99, by = 0.01), num_simulations, replace = TRUE), 2)),
  land_price_agr = suppress_console(decisionSupport::random(rho_land_price_agr, n = num_simulations)),
  houses_per_ha = suppress_console(decisionSupport::random(rho_houses_per_ha, n = num_simulations)),
  house_price = suppress_console(decisionSupport::random(rho_house_price, n = num_simulations)),
  fire_fight_cost_per_ha = suppress_console(decisionSupport::random(rho_fire_fight_cost_per_ha, n = num_simulations))
)

# Initialize lists for results
mc_raw_list_dt <- list()
mc_results_list_dt <- list() # Store results for bind_rows at the end

# --- 2. Refactor the loop with data.table operations ---
# Instead of filtering inside the loop for each iteration, we can perform joins
# and calculations more efficiently.

# Pre-calculate `sal_costs` and merge `area_dat` outside the innermost loop
# for a given strategy and year, then iterate over thresholds.

# measure time of function call
start_time <- Sys.time()

for (s in 1:length(v_strategies)) {
  current_strategy <- v_strategies[s]
  cat("Processing strategy:", current_strategy, "\n")
  
  # Filter target_df_agg once per strategy to reduce repeated filtering
  target_strategy_dt <- target_df_agg[strategy == current_strategy]
  
  # Filter perc_data_long once per strategy for relevant 'time'
  area_strategy_dt <- perc_data_long[strategy == current_strategy & time == selected_time]
  
  for (y in 1:length(v_years)) {
    current_year <- v_years[y]
    cat("- Year:", current_year, "\n")
    
    # Filter target_df_agg for the current year
    sal_costs_dt <- target_strategy_dt[year == current_year, .(sal_int = threshold, sum_SAL_costs)]
    
    # Filter perc_data_long for the current year
    area_year_dt <- area_strategy_dt[year == current_year,
                                     .(percentile, abs_avoid_tot_ha, abs_avoid_agri_ha,
                                       abs_avoid_urbn_ha, abs_avoid_wild_ha, sal_int)]
    
    # Loop through thresholds. Here we can perform a non-equi join or a merge
    # that handles all thresholds at once if the structure allows, but given
    # that `sal_costs` and `area_dat` depend on `v_thresholds[t]`,
    # the inner loop for `t` might still be clearer for direct translation.
    
    for (t in 1:length(v_thresholds)) {
      current_threshold <- v_thresholds[t]
      cat("-- Threshold:", current_threshold, "\n")
      
      # Get sal_costs for the current threshold
      current_sal_costs <- sal_costs_dt[sal_int == current_threshold, sum_SAL_costs]
      if (length(current_sal_costs) == 0) {
        current_sal_costs <- 0 # Handle cases where no match is found, or decide how to proceed
        warning(paste("No sum_SAL_costs found for strategy", current_strategy, "year", current_year, "threshold", current_threshold))
      }
      
      # Get area_dat for the current threshold
      current_area_dt <- area_year_dt[sal_int == current_threshold,
                                      .(percentile, abs_avoid_tot_ha, abs_avoid_agri_ha,
                                        abs_avoid_urbn_ha, abs_avoid_wild_ha)]
      
      # Perform the join and calculations using data.table syntax
      # Join sim_dt with current_area_dt
      tmp_dt <- sim_dt[current_area_dt, on = "percentile", nomatch = 0] # nomatch=0 removes unmatched rows
      
      # Perform calculations directly on the data.table using := for in-place modification
      tmp_dt[, `:=`(
        benefit_agr = abs_avoid_agri_ha * land_price_agr,
        benefit_urb = abs_avoid_urbn_ha * houses_per_ha * house_price,
        benefit_ff = abs_avoid_tot_ha * fire_fight_cost_per_ha
      )]
      
      tmp_dt[, benefit_per_event := benefit_agr + benefit_urb + benefit_ff]
      tmp_dt[, benefit_total := benefit_per_event * number_fires_per_year] # Assuming number_fires_per_year is a single value or vector
      
      # Ensure sal_costs is correctly broadcasted if it's a single value
      tmp_dt[, net_benefit := benefit_total - current_sal_costs]
      
      # Store the raw results if needed
      mc_raw_list_dt <- append(mc_raw_list_dt, list(tmp_dt))
      
      # Calculate net_benefit_positive
      nbp <- pos_share(tmp_dt$net_benefit)
      
      # Store results
      mc_results_list_dt[[length(mc_results_list_dt) + 1]] <- data.table(
        strategy = current_strategy,
        year = current_year,
        sal_int = current_threshold,
        net_benefit_positive = nbp
      )
    }
  }
}

end_time <- Sys.time()
execution_time <- end_time - start_time
cat("Execution time:", execution_time, "seconds\n")

# --- 3. Combine results ---
# Use rbindlist for data.tables, which is faster than bind_rows for lists of data.tables
mc_results_dt <- rbindlist(mc_results_list_dt)
mc_raw_dt <- rbindlist(mc_raw_list_dt) # If you need the combined raw data



# same plot, color scale traffic light
ggplot(mc_results_dt, aes(x = sal_int, y = factor(year), fill = net_benefit_positive)) +
  geom_tile() +
  facet_wrap(~ strategy) +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", midpoint = 0.5) +
  labs(title = "Monte Carlo Simulation Results",
       x = "SAL Intensity",
       y = "Year",
       fill = "Chance of Benefits > Costs") +
  theme_minimal()






