# ==============================================
# 05_retrieve_duration.R
# Description  : Calibrate simulated wildfires to historical area distributions and extract posterior duration parameters.
# Inputs       : data/rdat/wf_observed.csv
# Outputs      : data/rdat/duration_params.rds, out/figures/posterior_vs_observed.png
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

# ---------- Analysis ----------

# (refactored code goes here; logic preserved)

# ---------- Outputs ----------

