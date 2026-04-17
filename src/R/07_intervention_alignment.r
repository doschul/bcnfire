# Intervention alignment with AGP priority zones
# - For each scenario-year-intensity combination, compute:
#   - share of cells intervened (SAL share)
#   - mean AGP share among intervened cells (AGP share)

library(sf)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(here)

# ----------------------------
# 0) Read inputs
# ----------------------------
grd <- st_read("./data/rdat/grd_filt.gpkg", quiet = TRUE)
row.names(grd) <- paste0("cell_", seq_len(nrow(grd)))

target_df <- readRDS(file = "./data/rdat/target_df.rds")

agp_path <- here::here("data", "AGP_shp", "AGP.shp")
agp <- st_read(agp_path, quiet = TRUE)

# transform to same CRS as grd and clip to AOI
agp <- st_transform(agp, st_crs(grd))
agp_clipped <- st_intersection(agp, st_union(grd))

table(agp_clipped$TipusAGP)
#    AFG AFG_IUF     IUF     PEG PEG_IUF 
#    255      57     532     221      44 

# ----------------------------
# 1) Add cell_id and cell area
# ----------------------------
grd$cell_id <- row.names(grd)
grd$cell_area <- as.numeric(st_area(grd))

# optional sanity check: if all cells are equal-sized, weighting won't matter much
summary(grd$cell_area)

# ----------------------------
# 2) Intersect AGP with grid and compute area shares per cell x AGP class
# ----------------------------
# keep only needed columns
agp_clipped <- agp_clipped %>%
  select(TipusAGP)

# intersection grid x agp
agp_x_grd <- st_intersection(
  grd %>% select(cell_id, cell_area),
  agp_clipped
)

# area of overlap
agp_x_grd$int_area <- as.numeric(st_area(agp_x_grd))

# summarise to cell x AGP type
agp_share_long <- agp_x_grd %>%
  st_drop_geometry() %>%
  group_by(cell_id, TipusAGP, cell_area) %>%
  summarise(agp_share = sum(int_area) / first(cell_area), .groups = "drop")

# wide version: one column per AGP class
agp_share_wide <- agp_share_long %>%
  select(cell_id, TipusAGP, agp_share) %>%
  pivot_wider(
    names_from  = TipusAGP,
    values_from = agp_share,
    values_fill = 0,
    names_prefix = "agp_share_"
  )

# make sure all expected columns exist
expected_agp_cols <- c("agp_share_AFG", "agp_share_AFG_IUF", "agp_share_IUF",
                       "agp_share_PEG", "agp_share_PEG_IUF")
for (nm in expected_agp_cols) {
  if (!nm %in% names(agp_share_wide)) agp_share_wide[[nm]] <- 0
}

# total AGP share per cell
agp_share_wide <- agp_share_wide %>%
  mutate(
    agp_share_total = agp_share_AFG + agp_share_AFG_IUF + agp_share_IUF +
                      agp_share_PEG + agp_share_PEG_IUF
  )

# just in case of tiny topology/numeric issues
agp_share_wide <- agp_share_wide %>%
  mutate(
    agp_share_total = pmin(agp_share_total, 1)
  )

# merge back to grd if wanted
grd <- grd %>%
  left_join(agp_share_wide, by = "cell_id")

# ----------------------------
# 3) Prepare target_df
# ----------------------------
# enforce cell_id
if (!"cell_id" %in% names(target_df)) {
  target_df$cell_id <- row.names(target_df)
}

# keep only scenario mask columns
mask_cols <- grep("^new_mask_.*_t[0-9]+_y[0-9]+$", names(target_df), value = TRUE)

if (length(mask_cols) == 0) {
  stop("No scenario columns matching '^new_mask_.*_t[0-9]+_y[0-9]+$' found in target_df.")
}

# convert to data.table for speed
target_dt <- as.data.table(target_df[, c("cell_id", mask_cols), drop = FALSE])

# long format: one row per cell x scenario-year-intensity
target_long <- melt(
  target_dt,
  id.vars = "cell_id",
  measure.vars = mask_cols,
  variable.name = "scenario_name",
  value.name = "sal_bin"
)

# make sure binary is numeric 0/1
target_long[, sal_bin := as.numeric(sal_bin)]
target_long[is.na(sal_bin), sal_bin := 0]

# parse strategy / intensity / year
# example: new_mask_sal_costs_cell_t30_y2044
target_long[, c("strategy", "intensity", "year") := {
  m <- str_match(scenario_name, "^new_mask_(.+?)_t([0-9]+)_y([0-9]+)$")
  .(m[,2], as.integer(m[,3]), as.integer(m[,4]))
}]

# check parsing
if (any(is.na(target_long$strategy) | is.na(target_long$intensity) | is.na(target_long$year))) {
  bad_names <- unique(target_long[is.na(strategy) | is.na(intensity) | is.na(year), scenario_name])
  stop("Could not parse some scenario names:\n", paste(bad_names, collapse = "\n"))
}

# ----------------------------
# 4) Join AGP shares to each cell
# ----------------------------
agp_dt <- as.data.table(st_drop_geometry(grd)[, c(
  "cell_id",
  "agp_share_total",
  "agp_share_AFG",
  "agp_share_AFG_IUF",
  "agp_share_IUF",
  "agp_share_PEG",
  "agp_share_PEG_IUF"
)])

target_long <- merge(
  target_long,
  agp_dt,
  by = "cell_id",
  all.x = TRUE
)

# fill missing AGP shares with 0
agp_cols <- c("agp_share_total", "agp_share_AFG", "agp_share_AFG_IUF",
              "agp_share_IUF", "agp_share_PEG", "agp_share_PEG_IUF")

for (nm in agp_cols) {
  target_long[is.na(get(nm)), (nm) := 0]
}

# ----------------------------
# 5) Summarise to one point per strategy x intensity x year
# ----------------------------
# Interpretation:
# overlap_area = SAL ∩ AGP area
# x_agp_covered = share of total AGP area covered by SAL
# y_sal_in_agp  = share of SAL area falling inside AGP

# use explicit cell area if available
cell_dt <- as.data.table(st_drop_geometry(grd)[, c("cell_id", "cell_area")])

target_long <- merge(target_long, cell_dt, by = "cell_id", all.x = TRUE)

scenario_summary <- target_long[, {
  sal_ind <- as.numeric(sal_bin > 0)

  overlap_area <- sum(sal_ind * agp_share_total * cell_area, na.rm = TRUE)
  agp_area     <- sum(agp_share_total * cell_area, na.rm = TRUE)
  sal_area     <- sum(sal_ind * cell_area, na.rm = TRUE)

  # class-specific numerators and denominators
  num_AFG     <- sum(sal_ind * agp_share_AFG * cell_area, na.rm = TRUE)
  den_AFG     <- sum(agp_share_AFG * cell_area, na.rm = TRUE)

  num_AFG_IUF <- sum(sal_ind * agp_share_AFG_IUF * cell_area, na.rm = TRUE)
  den_AFG_IUF <- sum(agp_share_AFG_IUF * cell_area, na.rm = TRUE)

  num_IUF     <- sum(sal_ind * agp_share_IUF * cell_area, na.rm = TRUE)
  den_IUF     <- sum(agp_share_IUF * cell_area, na.rm = TRUE)

  num_PEG     <- sum(sal_ind * agp_share_PEG * cell_area, na.rm = TRUE)
  den_PEG     <- sum(agp_share_PEG * cell_area, na.rm = TRUE)

  num_PEG_IUF <- sum(sal_ind * agp_share_PEG_IUF * cell_area, na.rm = TRUE)
  den_PEG_IUF <- sum(agp_share_PEG_IUF * cell_area, na.rm = TRUE)

  list(
    n_cells = .N,
    n_sal_cells = sum(sal_ind, na.rm = TRUE),

    overlap_area = overlap_area,
    agp_area = agp_area,
    sal_area = sal_area,

    # global metrics
    x_agp_covered = if (agp_area > 0) overlap_area / agp_area else NA_real_,
    y_sal_in_agp  = if (sal_area > 0) overlap_area / sal_area else NA_real_,

    # class-specific X = share of AGP class covered by SAL
    x_AFG     = if (den_AFG > 0) num_AFG / den_AFG else NA_real_,
    x_AFG_IUF = if (den_AFG_IUF > 0) num_AFG_IUF / den_AFG_IUF else NA_real_,
    x_IUF     = if (den_IUF > 0) num_IUF / den_IUF else NA_real_,
    x_PEG     = if (den_PEG > 0) num_PEG / den_PEG else NA_real_,
    x_PEG_IUF = if (den_PEG_IUF > 0) num_PEG_IUF / den_PEG_IUF else NA_real_,

    # class-specific Y = share of SAL area falling in AGP class
    y_AFG     = if (sal_area > 0) num_AFG / sal_area else NA_real_,
    y_AFG_IUF = if (sal_area > 0) num_AFG_IUF / sal_area else NA_real_,
    y_IUF     = if (sal_area > 0) num_IUF / sal_area else NA_real_,
    y_PEG     = if (sal_area > 0) num_PEG / sal_area else NA_real_,
    y_PEG_IUF = if (sal_area > 0) num_PEG_IUF / sal_area else NA_real_
  )
}, by = .(strategy, intensity, year, scenario_name)]


# ----------------------------
# 6) Scale axes to [0,1]
# ----------------------------
# X scaling: min/max SAL share across all scenario-year-intensity combinations
#sal_min <- min(scenario_summary$sal_share_raw, na.rm = TRUE)
#sal_max <- max(scenario_summary$sal_share_raw, na.rm = TRUE)

# Y scaling: min/max AGP cell share across all cells (static baseline)
#agp_min <- min(agp_dt$agp_share_total, na.rm = TRUE)
#agp_max <- max(agp_dt$agp_share_total, na.rm = TRUE)

# robust min-max scaler
minmax_scale <- function(x, xmin, xmax) {
  if (isTRUE(all.equal(xmin, xmax))) {
    return(rep(0.5, length(x)))
  } else {
    return((x - xmin) / (xmax - xmin))
  }
}

#scenario_summary[, sal_scaled := minmax_scale(sal_share_raw, sal_min, sal_max)]
#scenario_summary[, agp_scaled := minmax_scale(agp_share_raw, agp_min, agp_max)]

# ----------------------------
# 7) Add labels / quadrant info
# ----------------------------
#scenario_summary[, quadrant := fifelse(
#  sal_scaled >= 0.5 & agp_scaled >= 0.5, "High SAL / High AGP",
#  fifelse(
#    sal_scaled < 0.5 & agp_scaled < 0.5, "Low SAL / Low AGP",
#    fifelse(
#      sal_scaled < 0.5 & agp_scaled >= 0.5, "High AGP / Low SAL",
#      "High SAL / Low AGP"
#    )
#  )
#)]

# prettier strategy names
scenario_summary[, strategy_pretty := fcase(
  grepl("bau_ROS_", strategy), "Rate of spread",
  grepl("FireConn", strategy), "Fire Connectivity",
  strategy == "sal_costs_cell", "Min. logging costs",
  strategy == "wui_risk_score", "WUI exposure",
  strategy == "random_target", "Random",
  default = strategy
)]

scenario_summary[, strategy_pretty := factor(
  strategy_pretty,
  levels = c("Random", "Min. logging costs", "WUI exposure", "Rate of spread", "Fire Connectivity")
)]

# alpha for year scaled within available years
yr_min <- min(scenario_summary$year, na.rm = TRUE)
yr_max <- max(scenario_summary$year, na.rm = TRUE)

scenario_summary[, year_alpha := minmax_scale(year, yr_min, yr_max)]
# avoid too faint points
scenario_summary[, year_alpha := 0.35 + 0.65 * year_alpha]

# intensity as ordered factor for size legend if preferred
scenario_summary[, intensity_f := factor(intensity, levels = sort(unique(intensity)))]

# save and reload
saveRDS(scenario_summary, "./data/rdat/scenario_agp_alignment_summary.rds")
scenario_summary <- readRDS("./data/rdat/scenario_agp_alignment_summary.rds")


# ----------------------------
# 8) Plot
# ----------------------------
p <- ggplot(scenario_summary, aes(x = x_agp_covered, y = y_sal_in_agp)) +
  #annotate("rect", xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1,
  #         alpha = 0.08, fill = "darkgreen") +
  #annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.5,
  #         alpha = 0.06, fill = "grey50") +
  #annotate("rect", xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1,
  #         alpha = 0.06, fill = "goldenrod") +
  #annotate("rect", xmin = 0.5, xmax = 1, ymin = 0, ymax = 0.5,
  #         alpha = 0.08, fill = "firebrick") +

  #geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.4, colour = "grey30") +
  #geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.4, colour = "grey30") +

  geom_point(
    aes(
      colour = strategy_pretty,
      size   = intensity,
      alpha  = year_alpha
    ),
    shape = 16
  ) +

  scale_x_continuous(
    #limits = c(0, 1),
    #breaks = seq(0, 1, 0.25),
    labels = number_format(accuracy = 0.01),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    #limits = c(0, 1),
    #breaks = seq(0, 1, 0.25),
    labels = number_format(accuracy = 0.01),
    expand = c(0.01, 0.01)
  ) +
  scale_alpha_identity(
    guide = "none"
  ) +
  scale_size_continuous(
    breaks = sort(unique(scenario_summary$intensity))
  ) +

  labs(
    x = "AGP coverage by SAL",
    y = "SAL share inside AGP",
    colour = "Targeting strategy",
    size = "SAL intensity",
    title = "Alignment of salvage logging with AGP priority zones",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p)

# ----------------------------
# 9) Optional: save outputs
# ----------------------------
#saveRDS(scenario_summary, "./data/rdat/scenario_agp_alignment_summary.rds")
ggsave(
  filename = "./out/fig/scenario_agp_alignment_pointplot.png",
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)



# exploratory plots
library(data.table)
library(ggplot2)
library(ggtern)

# --------------------------------------------------
# 1) Build ternary input from scenario_summary
# --------------------------------------------------
tern_dt <- as.data.table(scenario_summary)[, .(
  strategy_pretty,
  intensity,
  year,
  PEG = y_PEG + y_PEG_IUF,
  AFG = y_AFG + y_AFG_IUF,
  IUF = y_IUF
)]

tern_dt <- tern_dt[!is.na(PEG) & !is.na(AFG) & !is.na(IUF)]

# guard against tiny negative values from numeric noise before closing the simplex
tern_dt[, `:=`(
  PEG = pmax(PEG, 0),
  AFG = pmax(AFG, 0),
  IUF = pmax(IUF, 0)
)]

# 2049 vs all others
tern_dt[, year_group := ifelse(year == 2049, "2049", "Other years")]
tern_dt[, year_group := factor(year_group, levels = c("Other years", "2049"))]

# close composition to 1 within plotted components
tern_dt[, tern_sum := PEG + AFG + IUF]
tern_dt <- tern_dt[is.finite(tern_sum) & tern_sum > 0]

tern_dt[, `:=`(
  PEG_plot = PEG / tern_sum,
  AFG_plot = AFG / tern_sum,
  IUF_plot = IUF / tern_sum
)]

# force exact closure to 1 to avoid ggtern rounding warnings
tern_dt[, tern_rescale := PEG_plot + AFG_plot + IUF_plot]
tern_dt <- tern_dt[is.finite(tern_rescale) & tern_rescale > 0]
tern_dt[, `:=`(
  PEG_plot = PEG_plot / tern_rescale,
  AFG_plot = AFG_plot / tern_rescale,
  IUF_plot = IUF_plot / tern_rescale
)]
tern_dt[, c("tern_sum", "tern_rescale") := NULL]

if (nrow(tern_dt) == 0) {
  stop("No valid ternary plot rows after filtering and simplex normalisation.")
}

# --------------------------------------------------
# 2) Aggregate centroids
# --------------------------------------------------
centroids <- tern_dt[, .(
  PEG_plot = mean(PEG_plot, na.rm = TRUE),
  AFG_plot = mean(AFG_plot, na.rm = TRUE),
  IUF_plot = mean(IUF_plot, na.rm = TRUE)
), by = .(strategy_pretty, intensity, year_group)]

centroids <- centroids[
  is.finite(PEG_plot) & is.finite(AFG_plot) & is.finite(IUF_plot)
]

# --------------------------------------------------
# 3) Palette and symbols
# --------------------------------------------------
cols <- c(
  "Random"             = "#7F7F7F",
  "Min. logging costs" = "#E69F00",
  "WUI exposure"       = "#009E73",
  "Rate of spread"     = "#0072B2",
  "Fire Connectivity"  = "#CC79A7"
)

shape_vals <- c(
  "Other years" = 16,  # filled circle
  "2049"        = 17   # filled triangle
)

# --------------------------------------------------
# 4) Plot: put ternary mapping globally
# --------------------------------------------------
p_tern <- ggtern(
  data = tern_dt,
  aes(x = PEG_plot, y = AFG_plot, z = IUF_plot)
) +
  geom_point(
    aes(color = strategy_pretty, shape = year_group),
    size = 1.6,
    alpha = 0.18,
    stroke = 0.2,
    na.rm = TRUE
  ) +
  geom_point(
    data = centroids,
    mapping = aes(
      x = PEG_plot, y = AFG_plot, z = IUF_plot,
      color = strategy_pretty,
      shape = year_group,
      size  = intensity
    ),
    inherit.aes = FALSE,
    alpha = 0.75,
    stroke = 0.9,
    na.rm = TRUE
  ) +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_shape_manual(values = shape_vals, drop = FALSE) +
  scale_size_continuous(
    breaks = sort(unique(centroids$intensity)),
    range = c(1.5, 3)
  ) +
  labs(
    title = "Composition of SAL overlap across AGP functional classes",
    subtitle = paste(
      "Small points: scenario-year values.",
      "Large points: centroids by strategy x intensity.",
      "Triangles indicate 2049."
    ),
    T = "PEG (control)",
    L = "AFG (buffer)",
    R = "IUF (exposure)",
    color = "Targeting strategy",
    shape = "Year group",
    size  = "SAL intensity"
  ) +
  theme_bw(base_size = 12) +
  theme_showarrows() +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold")
  )

print(p_tern)

# --------------------------------------------------
# 7) Save
# --------------------------------------------------
ternary_outfile <- "./out/fig/scenario_agp_alignment_ternary_2049_vs_others.png"
dir.create(dirname(ternary_outfile), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = ternary_outfile,
  plot = p_tern,
  width = 10,
  height = 12,
  dpi = 300
)
