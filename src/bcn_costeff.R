# Calculate cost-effectiveness of targeting strategies

##### Setup #####

rm(list=ls())

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

# libraries
library(tidyverse)
library(scales)


# load data
summary_relative <- read_rds("./data/rdat/summary_relative.rds")
target_df_agg <- readRDS(file = "./data/rdat/target_df_agg.rds")

# rename strtegies to match
target_df_agg <- target_df_agg %>%
  mutate(strategy = case_when(
    strategy == "Highest artificial" ~ "High Built-up",
    strategy == "Highest biomass" ~ "High Biomass",
    strategy == "Highest connectivity" ~ "High Connectivity",
    strategy == "Low logging costs" ~ "Low SAL cost",
    TRUE ~ strategy
  ))

table(unique(target_df_agg$strategy) %in% unique(summary_relative$strategy))

# join data
costeff_df <- summary_relative %>%
  # exclude BAU (no cost effectiveness)
  filter(sal_int > 0) %>%
  mutate(year = as.character(year)) %>%
  left_join(target_df_agg, 
            by = c("year" = "year", 
                   "strategy" = "strategy",
                   "sal_int" = "threshold"
            )) %>%
  mutate(
    cost_eff = sum_SAL_costs / hectares_value,
    log_cost_eff = log(cost_eff))

# log-costs normally distributed:
hist(costeff_df$log_cost_eff)


##### Plotting #####

# --- Prepare Data for Heatmap ---
# Filter summary_data to get only the median of burned_land_ha
heatmap_data <- costeff_df %>%
  filter(#percentile_type == "Median (50th Percentile)", # We want the median value
         burned_type %in% c("burned_land", "burned_urbn")) %>%
  mutate(sal_int_label = factor(sal_int, levels = sort(unique(sal_int)))) %>%
  mutate(burned_type = case_when(burned_type == "burned_land" ~ "Burned Land (Total ha)",
                                 burned_type == "burned_urbn" ~ "Burned Urban area (ha)", 
                                 TRUE ~ burned_type),
         time = factor(time),
         percentile_type = factor(percentile_type))

heatmap_data <-heatmap_data %>%
  filter(percentile_type %in%  c("99th Percentile", "Median (50th Percentile)"),
         time %in% c(30, 60, 120, 300, 900))

# Recalculate rect_data using these indices
rect_data_final <- data.frame(
  xmin = 5.5, xmax = 15.5,
  ymin = 10.5, ymax = 33.5,
  percentile_type = "99th Percentile",
  burned_type = "Burned Land (Total ha)"
)

# ----------------------------------
# 1. Overview (big picture)
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
    fill = NA, # No fill
    color = "black", # Highlight color
    linewidth = 1, # Thicker line
    inherit.aes = FALSE # Do not inherit main plot aesthetics
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
    panel.spacing = unit(0.5, "lines"),  # spacing only between facets
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# ----------------------------------
# 2. Detailed (zoomed) plot
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
    #panel.spacing = unit(0, "lines"),  # no space between facets
    panel.spacing.x = unit(0, "pt"),  # ensure zero spacing
    panel.spacing.y = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# ----------------------------------
# 3. Combine
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

options(scipen = 999)

heatmap_data$log_cost_eff <- log(1+heatmap_data$cost_eff)
heatmap_data$time_step = as.integer(heatmap_data$time)/5
heatmap_data$time2 = heatmap_data$time_step * as.integer(heatmap_data$time)
heatmap_data$sal_int2 = heatmap_data$sal_int * heatmap_data$sal_int
heatmap_data$strategy <- relevel(factor(heatmap_data$strategy), ref = "Random")
summary(lm(log_cost_eff ~ strategy+ 
             factor(year) + 
             factor(percentile_type) + 
             factor(burned_type) +
             time_step + 
             time2 + 
             sal_int + 
             sal_int2,  
           data = heatmap_data %>% 
             filter(!log_cost_eff == Inf) ))




