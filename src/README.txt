Overview of scripts (ordered, start by running the first)


0. bcn_funs: Contains set of helper functions (sourced at beginning of other scripts)

1. bcn_fire_history: Descriptive and predictive analysis of historical wildfire patterns.
	- Requires: Raw spatial wildfire observations, secondary data
	- Provides: wf_observed.csv (for later calibration), ignition_probability_layer
	- Output: Figure of historic fire distribution, Map of ignition probability

2. data_prep: Create regular grid, load and extract spatial data to grid, filter grid to AOI, create spatial neighbor matrix
	- Requires: Raw spatial datasets (connectivity, rate of spread, landcover, wildland urban interface, ignition_probability_layer)
	- Provides: grd_filt, neighbor_list
	- Output: Map of study region

3. bcn_target: Defines targeting variables and strategies, implements salvage logging values according to strategy.
	- Requires: grd_filt
	- Provides: target_df,summary_relative
	- Output: Figures of SAL costs, absolute and relative costs of SAL under targeting

4. bcn_firesim: implements fire spread simulation scenarios across ignition cells, years, strategies. Aggregates resulting affected lancover shares per time step.
	- Requires: target_df, neighbor_list
	- Provides: long_df_agg (also intermediate batches), perc_data_long (percentile outcome)
	- Output: Heatmap of SAL costs per avoided burned areas (cost-effectiveness per ha, independent of land cover/land use)

5. retrieve_duration: Calibrates simulated wildfires to match historical affected area patterns
	- Requires: wf_observed.csv
	- Provides: Calibrated duration distribution parameters
	- Output: Posterior vs. observed percentile plot

6. bcn_mc_sim: Monte Carlo Simulation of avoided costs.
	- Requires: perc_data_long, distribution parameters (duration, costs...)
	- Provides: 