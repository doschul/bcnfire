# load raster data

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

rm(list = ls())

library(sf)
library(tidyverse)

# get functions
source("./src/bcn_funs.R")

# load grid
grd <- st_read("./data/rdat/grd.gpkg") %>%
  # transform cs to lon lat
  st_transform(., 4326)

# every pixel is 4 ha, use lc_wildland to get total wildland area
total_wildland <- sum(grd$lc_wildland[!grd$biomass_bau == grd$biomass_sal]) * 4 # in ha

# load neighbors
load("./data/rdat/neighbor_idx.RData")


##### Monte Carlo Sim #####

# run 100 random ignition cells, store full result in list

rdm_ign_cells <- sample(x = unique(grd$id), 
                        size = 100, 
                        replace = FALSE,
                        prob = grd$ignition_probability_surface)
mc_time_step <- 5
mc_time_horizon <- 36


# run lapply in parallel using future apply and show progress bar
library(future)
library(future.apply)
library(progressr)

handlers(handler_progress(format="[:bar] :percent :eta :message"))

plan(multisession, workers=7)

with_progress({
  p <- progressor(along=rdm_ign_cells)
  mc_res_bau <- future_lapply(rdm_ign_cells,  function(x) {
    p() # monitor progress
    get_burners(x, 
                time_horizon = mc_time_horizon, 
                time_step = mc_time_step,
                full = TRUE, 
                scenario = "bau",
                grd = grd, 
                neighbor_idx = neighbor_idx)})
}) 

agg_bau <- mc_res_bau  %>%
  lapply(., function(x) apply(x, 2, sum)) %>%
  do.call(bind_rows, .)  %>%
  as.data.frame(.) %>%
  mutate(scenario = "BAU")

with_progress({
  p <- progressor(along=rdm_ign_cells)
  mc_res_sal <- future_lapply(rdm_ign_cells,  function(x) {
    p() # monitor progress
    get_burners(x, 
                time_horizon = mc_time_horizon, 
                time_step = mc_time_step,
                full = TRUE, 
                scenario = "sal",
                grd = grd, 
                neighbor_idx = neighbor_idx)})
}) 

agg_sal <- mc_res_sal  %>%
  lapply(., function(x) apply(x, 2, sum)) %>%
  do.call(bind_rows, .)  %>%
  as.data.frame(.) %>%
  mutate(scenario = "SAL")



# merge and bring into long format
mc_res_long <- bind_rows(agg_bau, agg_sal) %>%
  mutate(id = NULL) %>%
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




### landcover multiplication test

# within each model run, multiply the landcover values with the burned cells, matching by cell id
# then sum the burned cells for each landcover class
# landcover shares per gridcell is stored in grd

lc_dat <- grd %>%
  st_set_geometry(NULL) %>%
  select(all_of(c("id", "lc_agriculture", "lc_wildland", "lc_rock",
                  "lc_artificial","lc_water")))



# multipy burned cells in mc_res list items with landcover values for that cell
# then sum the burned cells for each landcover class
# landcover shares per gridcell is stored in grd
mc_bau_lc <- lapply(mc_res_bau, function(x) {
  
  ignition_cell <- x$id[x$V1 == 1]
  
  x %>%
    # long format all columns starting with V to time, keep id
    # use column name pattern, e.g. V8 after "V" for time values 8
    pivot_longer(-id, names_to = "time", values_to = "burned", names_prefix = "V") %>%
    mutate(time = as.numeric(time)) %>%
    # join with landcover
    left_join(., lc_dat, by = "id") %>%
    group_by(time) %>%
    summarise(acres_burned = sum(burned),
              burned_agric = sum( burned * lc_agriculture),
              burned_wild = sum(burned * lc_wildland),
              burned_urban = sum(burned * lc_artificial)
    ) %>%
    ungroup() %>%
    mutate(ignition_cell = ignition_cell)
}) %>% do.call(bind_rows, .) %>%
  as.data.frame(.) %>%
  mutate(scenario = "BAU",
         id = 1:nrow(.))

mc_sal_lc <- lapply(mc_res_sal, function(x) {
  
  ignition_cell <- x$id[x$V1 == 1]
  
  x %>%
    # long format all columns starting with V to time, keep id
    # use column name pattern, e.g. V8 after "V" for time values 8
    pivot_longer(-id, names_to = "time", values_to = "burned", names_prefix = "V") %>%
    mutate(time = as.numeric(time)) %>%
    # join with landcover
    left_join(., lc_dat, by = "id") %>%
    group_by(time) %>%
    summarise(acres_burned = sum(burned),
              burned_agric = sum( burned * lc_agriculture),
              burned_wild = sum(burned * lc_wildland),
              burned_urban = sum(burned * lc_artificial)
    ) %>%
    ungroup() %>%
    mutate(ignition_cell = ignition_cell)
}) %>% do.call(bind_rows, .) %>%
  as.data.frame(.) %>%
  mutate(scenario = "SAL", 
         id = 1:nrow(.))

# boxplots for each landcover class
mc_bau_lc_long <- bind_rows(mc_bau_lc, mc_sal_lc) %>%
  pivot_longer(-c(scenario,time, id, ignition_cell), names_to = "landcover", values_to = "burned") %>%
  mutate(landcover = gsub("burned_", "", landcover))

# save result
save(mc_bau_lc_long, file = "./out/mc_bau_lc_long.RData")


# density plot of last timestep comparing scenarios
p.dens.mc.urb <- mc_bau_lc_long %>% 
  filter(landcover == "urban") %>%
  filter(time == max(time)) %>%
  ggplot(aes(x = burned, fill = scenario)) +
  geom_density(alpha = 0.5) +
  labs(x = "Burned urban surface", y = "Density") +
  theme_minimal()

p.dens.mc.agr <- mc_bau_lc_long %>% 
  filter(landcover == "agric") %>%
  filter(time == max(time)) %>%
  ggplot(aes(x = burned, fill = scenario)) +
  geom_density(alpha = 0.5) +
  labs(x = "Burned agricultural surface", y = "Density") +
  theme_minimal()

p.dens.mc.wld <- mc_bau_lc_long %>% 
  filter(landcover == "wild") %>%
  filter(time == max(time)) %>%
  ggplot(aes(x = burned, fill = scenario)) +
  geom_density(alpha = 0.5) +
  labs(x = "Burned shrub/forest surface", y = "Density") +
  theme_minimal()

# combine plots
library(ggpubr)

p.dens.mc <- ggarrange(p.dens.mc.urb, p.dens.mc.agr, p.dens.mc.wld, nrow = 3, 
          common.legend = TRUE, legend = "bottom", align = "hv")

# save
ggsave(filename = "./out/fig/p.dens.mc.png", p.dens.mc, height = 8, width = 5, bg = "white")


# calculate difference between scenarios for each landcover class and time step
mc_diff <- mc_bau_lc_long %>%
  group_by(time, landcover, id) %>%
  summarise(burn_diff = diff(burned))


# plot histogram for each landcover difference
mc_diff %>%
  filter(time == 25) %>%
  filter(landcover == "urban") %>%
  ggplot(aes(x = burn_diff)) +
  geom_histogram() +
  #facet_wrap(~landcover, scales = "free_x") +
  labs(x = "Difference in burned cells", y = "Density") +
  theme_minimal()


urb_diff25 <- mc_diff %>%
  filter(time == 25) %>%
  filter(landcover == "urban")

quantile(urb_diff25$burn_diff, probs = seq(0, 1, 0.1))

hist(mc_diff$burn_diff, breaks = 200)

# density plot for urban difference with time as color
mc_diff %>%
  filter(landcover == "urban", 
         time > 6) %>%
  ggplot(aes(x = burn_diff, fill = factor(time))) +
  # viridis color
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  # transform y axis to logscale
  #scale_y_sqrt() +
  labs(x = "Difference in burned urban cells", y = "Density") +
  theme_minimal() +
  # axis range -1 to 0
  coord_cartesian(xlim = c(-1, 0))
