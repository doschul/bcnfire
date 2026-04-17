# %%
import jax
import jax.numpy as jnp
from jax import random, vmap
import numpy as np
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
import jax.scipy as jsp
import pandas as pd

# %% Actual data

directory = "C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire"
sim_dat_path = directory + "/data/rdat/mc_res_agg_random_sal0.csv"
obs_dat_path = directory + "/data/wildfires/wf_observed.csv"

sim_dat = pd.read_csv(sim_dat_path)
obs_dat = pd.read_csv(obs_dat_path)

# select 100 random ignition cells
rdm_ign_cell_idx = np.random.choice(sim_dat["ignition_id"].unique(), size=100, replace=False)
sim_dat = sim_dat[sim_dat["ignition_id"].isin(rdm_ign_cell_idx)]

# subset simdat to year 2040
sim_dat = sim_dat[sim_dat["year"] == 2040]

# rename ignition_id to sim_id for clarity
sim_dat = sim_dat.rename(columns={"ignition_id": "sim_id"})

observed_areas = jnp.array(obs_dat["area_ha"].values)
times = jnp.array(sorted(sim_dat["time"].unique()))
sim_matrix = jnp.array(sim_dat.pivot(index="time", columns="sim_id", values="burned_land_ha").values)
print("sim_matrix shape:", sim_matrix.shape)  # (T, M)
print("times shape:", times.shape)  # (T,)
print("observed_areas shape:", observed_areas.shape)  # (N,)


# %%
# ------------------------
# Interpolation utilities
# ------------------------
# We'll build two vmapped helpers:
#  - interp_one_sim(t, sim_curve) -> area at t for that sim
#  - interp_all_sims(t, sim_matrix) -> vector length M of areas for all sims at t
# We then vmap interp_all_sims over N durations to get (N, M).

# interp for a single sim curve (1D fp) and scalar t -> scalar area
def interp_one_sim(t_scalar, sim_curve, xp=times):
    # jnp.interp supports scalar x and 1d fp
    return jnp.interp(t_scalar, xp, sim_curve)

# vectorized over sims: for a given scalar t, returns array (M,)
interp_all_sims_for_scalar_t = vmap(lambda fp: interp_one_sim(fp_t, fp), in_axes=(0), out_axes=0)  # placeholder

# We'll instead define a function using nested vmaps:
def sim_areas_at_t(t_scalar, sim_matrix_local, times_local):
    """Return simulated areas (M,) for all sims interpolated at scalar t_scalar."""
    # sim_matrix_local: shape (T, M)
    # We vmap interp_one_sim over sim curves (columns)
    sim_cols = jnp.transpose(sim_matrix_local)  # shape (M, T)
    # vmap over M: each fp is length T, we call jnp.interp(t_scalar, xp, fp)
    return vmap(lambda fp: jnp.interp(t_scalar, times_local, fp))(sim_cols)

# Vectorize over an array of durations t_vec -> returns (len(t_vec), M)
sim_areas_at_t_vmap = vmap(sim_areas_at_t, in_axes=(0, None, None), out_axes=0)


# ------------------------
# Log-mean-exp helper (stable marginalization over simulated ignitions)
# ------------------------
def logmeanexp(a, axis=-1):
    """Stable log(mean(exp(a), axis=axis))."""
    a_max = jnp.max(a, axis=axis, keepdims=True)
    return jnp.log(jnp.mean(jnp.exp(a - a_max), axis=axis)) + jnp.squeeze(a_max, axis=axis)

# ------------------------
# NumPyro model
# ------------------------
def continuous_duration_model(obs_areas, times, sim_matrix):
    """
    obs_areas: (N,)
    times: (T,)
    sim_matrix: (T, M)
    """
    N = obs_areas.shape[0]
    T, M = sim_matrix.shape

    # sigma: observation multiplicative noise (log-scale)
    sigma = numpyro.sample("sigma", dist.HalfNormal(0.5))

    # weakly-informative prior favoring shorter durations:
    # u ~ Beta(1, 2) (density ~ (1-u)), map to [times[0], times[-1]]
    with numpyro.plate("fires", N):
        u = numpyro.sample("u", dist.Beta(1.0, 2.0))  # favors lower u -> shorter times
        duration = times[0] + u * (times[-1] - times[0])  # (N,) continuous durations

        # For each duration, compute the vector of M simulated areas at that time:
        # sim_areas: shape (N, M)
        sim_areas = sim_areas_at_t_vmap(duration, sim_matrix, times)  # (N, M)

        # Avoid zeros / negative sims
        sim_areas = jnp.maximum(sim_areas, 1e-8)

        # compute log-likelihood for each (obs_i, sim_j): log p(A_obs_i | area_sim_ij(t), sigma)
        # dlnorm(A_obs | meanlog = log(area_sim), sdlog = sigma)
        # gives shape (N, M)
        loglike_matrix = dist.LogNormal(jnp.log(sim_areas), sigma).log_prob(obs_areas[:, None])

        # marginalize over ignitions (empirical mixture) -> log p(A_obs | t) = logmeanexp over M
        log_marginal_over_sims = logmeanexp(loglike_matrix, axis=1)  # shape (N,)

        # multiply by any prior on duration: prior implicit via u's Beta(1,2) already
        # register the contribution to the joint log-probability
        numpyro.factor("obs_loglike", jnp.sum(log_marginal_over_sims))




# %%
# ------------------------
# Run MCMC
# ------------------------
nuts_kernel = NUTS(continuous_duration_model)
mcmc = MCMC(nuts_kernel, num_warmup=300, num_samples=600, num_chains=2)
mcmc.run(random.PRNGKey(1), obs_areas=observed_areas, times=times, sim_matrix=sim_matrix)
posterior = mcmc.get_samples()

# %%
# posterior["u"] shape: (num_samples, N); transform to durations:
u_samples = posterior["u"]                     # (S, N)
duration_samples = times[0] + u_samples * (times[-1] - times[0])  # (S, N)
sigma_samples = posterior["sigma"]             # (S,)

# Posterior summaries
mean_duration_per_fire = jnp.mean(duration_samples, axis=0)
overall_mean_duration = jnp.mean(mean_duration_per_fire)
print("Overall posterior mean duration (minutes):", float(overall_mean_duration))
print("Posterior sigma (median):", float(jnp.median(sigma_samples)))

# %%
# ------------------------
# Posterior predictive sampling using simulations
# ------------------------
# For posterior predictive you want to sample simulated areas conditional on posterior durations.
# Procedure:
#  - For each posterior draw s and fire i:
#      - take duration d = duration_samples[s,i]
#      - interpolate the M sim areas at d -> vector sim_areas_j
#      - sample one ignition index j uniformly (or use the empirical mixture)
#      - sample an observation from LogNormal(log(sim_areas_j), sigma_samples[s])
#
# Below we produce K posterior predictive draws per fire (may be large):

S, N = duration_samples.shape  # S posterior samples, N observed fires
K = 1  # number of draws per posterior sample per fire (set >1 if you want more)
rng = np.random.default_rng(123)

# Flatten (s,i) to compute interpolation in batches
flat_durations = np.asarray(duration_samples).reshape(-1)  # length S*N

# Compute simulated areas for these durations
# sim_areas_at_t_vmap must return shape (len(flat_durations), M)
sim_areas_flat = sim_areas_at_t_vmap(jnp.array(flat_durations), sim_matrix, times)

# Pick one simulated ignition trajectory uniformly for each flattened entry
M = sim_matrix.shape[1]
rand_idx = rng.integers(0, M, size=sim_areas_flat.shape[0])
chosen_sim_areas = sim_areas_flat[jnp.arange(sim_areas_flat.shape[0]), rand_idx]  # (S*N,)

# --- Shape alignment fix ---
# Repeat sigma for each fire, tile simulated areas for each posterior sample
sigma_flat = np.repeat(np.asarray(sigma_samples), N)        # shape (S*N,)
sim_flat   = np.asarray(chosen_sim_areas).reshape(S * N)    # shape (S*N,)

# Sample observation noise: LogNormal(log(sim_area), sigma)
pp_draws = np.exp(np.random.normal(np.log(sim_flat), sigma_flat))  # shape (S*N,)
pp_draws = pp_draws.reshape(S, N)

# Posterior predictive summaries
pp_all = pp_draws.flatten()
calibrated_percentiles = np.percentile(pp_all, [1, 5, 10, 25, 50, 75, 90, 95, 99])
print("Calibrated percentiles of area (posterior predictive):")
print(calibrated_percentiles)

# %%
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(8,5))
sns.kdeplot(np.asarray(observed_areas), label="Observed", fill=True, alpha=0.4)
sns.kdeplot(pp_all, label="Posterior predictive", fill=True, alpha=0.4)
plt.xscale("log")
plt.xlabel("Burned area [ha] (log scale)")
plt.ylabel("Density")
plt.title("Posterior predictive fire size calibration")
plt.legend()
plt.tight_layout()
plt.show()

# %%

# Your existing objects:
# times, sim_matrix, observed_areas, duration_samples, sigma_samples

# Flatten posterior durations (S, N) -> (S*N,)
flat_durations = np.asarray(duration_samples).reshape(-1)

# Draw a random subset of posterior durations (to keep plotting light)
rng = np.random.default_rng(123)
subset_idx = rng.choice(len(flat_durations), size=2000, replace=False)
dur_subset = flat_durations[subset_idx]

# Draw the same number of random simulation columns
M = sim_matrix.shape[1]
sim_subset_idx = rng.integers(0, M, size=len(dur_subset))

# Interpolate simulated areas at these durations
def interp_one_curve(t, curve):
    return np.interp(t, np.asarray(times), np.asarray(curve))

# Vectorized interpolation for the sampled durations and simulations
sim_areas_subset = np.array([
    interp_one_curve(dur_subset[i], sim_matrix[:, sim_subset_idx[i]])
    for i in range(len(dur_subset))
])

# This is your *posterior-sampled simulated area* distribution.
# Convert to hectares if needed (already in ha here)
sim_calibrated_areas = sim_areas_subset

# Convert observed_areas to numpy
true_areas = np.asarray(observed_areas)

# Combine into a long dataframe for seaborn
df_plot = pd.DataFrame({
    "area_ha": np.concatenate([true_areas, sim_calibrated_areas]),
    "type": ["Observed"] * len(true_areas) + ["Simulated (Calibrated)"] * len(sim_calibrated_areas)
})

# Plot distributions (log-scale is helpful due to heavy tail)
plt.figure(figsize=(8, 5))
sns.kdeplot(data=df_plot, x="area_ha", hue="type", fill=True, common_norm=False, alpha=0.4)
plt.xscale("log")
plt.xlabel("Final burned area (ha, log scale)")
plt.ylabel("Density")
plt.title("Comparison of true vs calibrated simulated fire size distributions")
plt.legend(title="")
plt.tight_layout()
plt.show()

# %%
# histogram true area data
plt.figure(figsize=(8, 5))
plt.hist(true_areas, bins=30, alpha=0.5, label="Observed", density=True)
plt.hist(sim_calibrated_areas, bins=30, alpha=0.5, label="Simulated (Calibrated)", density=True)
#plt.xscale("log")
plt.xlabel("Final burned area (ha, log scale)")
plt.ylabel("Density")
plt.title("Histogram: true vs calibrated simulated fire size distributions")
plt.legend()
plt.tight_layout()
plt.show()
# %%
