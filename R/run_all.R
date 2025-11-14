# ==============================================
# run_all.R
# Description  : One-click reproduction wrapper that runs the full pipeline end-to-end.
# Inputs       : N/A
# Outputs      : All pipeline outputs
# Author       : (auto-refactor)
# Created      : (auto)
# ==============================================

# ---------- Setup ----------
suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(ggplot2)
  library(stringr)
  library(here)
})

set.seed(123)
options(stringsAsFactors = FALSE)
# Project root assumed to be Git repo top; use here::here()
data_raw_dir  <- here::here("data", "raw")
data_proc_dir <- here::here("data", "rdat")
out_fig_dir   <- here::here("out", "figures")
out_tab_dir   <- here::here("out", "tables")

# Source helper functions
funs_path <- here::here("R", "bcn_funs.R")
if (file.exists(funs_path)) source(funs_path)

# ---------- Pipeline ----------
scripts <- c(
  "R/01_bcn_fire_history_and_validation.R",
  "R/02_data_prep.R",
  "R/03_bcn_target.R",
  "R/04_bcn_firesim.R",
  "R/05_retrieve_duration.R",
  "R/06_bcn_mc_sim.R"
)
for (s in scripts) { 
  message(">>> Running: ", s)
  source(here::here(s), local = new.env())
}
message("Pipeline complete.")
