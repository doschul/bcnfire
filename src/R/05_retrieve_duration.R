# ==============================================
# 05_retrieve_duration.R
# Description  : Calibrate simulated wildfires to historical area distributions and extract posterior duration parameters.
# Inputs       : data/rdat/wf_observed.csv
# Outputs      : data/rdat/duration_params.rds, out/figures/posterior_vs_observed.png
# ==============================================

# ---------- Analysis ----------

# run ./Pyt/GEV_calibration.py to get the posterior distribution of duration parameters

# ---------- Outputs ----------
xi     <- -0.205   # replace with your ξ from the plot
mu     <- 94.9     # replace with your μ
sigma  <- 46.5    # replace with your σ
