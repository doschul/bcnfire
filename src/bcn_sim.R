# load raster data


library(sf)
library(tidyverse)

# get functions
source("bcn_funs.R")

# load grid
grd <- st_read("grd_bau.gpkg") %>%
  # transform cs to lon lat
  st_transform(., 4326)

# load neighbors
load("neighbor_idx.RData")


##### Monte Carlo Sim #####

# run 100 random ignition cells, store full result in list

rdm_ign_cells <- sample(grd$id, 200)
mc_time_step <- 10
mc_time_horizon <- 12


mc_res_bau <- lapply(rdm_ign_cells, 
                     get_burners, 
                     time_horizon = mc_time_horizon, 
                     time_step = mc_time_step,
                     full = TRUE, 
                     
                     scenario = "bau") %>%
  lapply(., function(x) apply(x, 2, sum)) %>%
  do.call(bind_rows, .)  %>%
  as.data.frame(.) %>%
  mutate(scenario = "BAU")

mc_res_sal <- lapply(rdm_ign_cells, 
                     get_burners, 
                     time_horizon = mc_time_horizon, 
                     time_step = mc_time_step,
                     full = TRUE, 
                     scenario = "sal") %>%
  lapply(., function(x) apply(x, 2, sum)) %>%
  do.call(rbind, .) %>%
  as.data.frame(.) %>%
  mutate(scenario = "SAL")

# merge and bring into long format
mc_res_long <- bind_rows(mc_res_bau, mc_res_sal) %>%
  pivot_longer(-scenario, names_to = "time", values_to = "burned") %>%
  mutate(time = as.numeric(gsub("V", "", time)) * mc_time_step)



# ggplot boxplot over time
p.mc.boxplot <- ggplot(mc_res_long, aes(x = factor(time), y = burned, fill = scenario)) +
  geom_boxplot() +
  labs(x = "Time (minutes)", y = "Burned cells") +
  theme_minimal()

# calculate difference in means between scenarios for each time
mc_res_long %>%
  group_by(time, scenario) %>%
  summarise(mean_burned = mean(burned)) %>%
  ungroup() %>%
  group_by(time) %>%
  summarise(burn_diff = diff(mean_burned))

# save plot
setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcn_fire")
ggsave(filename = "./out/fig/p.mc.boxplot.png", p.mc.boxplot, height = 5, width = 8, bg = "white")


