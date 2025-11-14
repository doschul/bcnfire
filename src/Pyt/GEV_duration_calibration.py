# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import jax
import jax.numpy as jnp
from jax.scipy.stats import beta as jbeta
import os
from scipy.stats import gaussian_kde, ks_2samp, anderson_ksamp, wasserstein_distance, entropy
from scipy.stats import genextreme   # SciPy uses shape = -ξ convention!
from scipy.optimize import minimize
jax.config.update("jax_enable_x64", True)

# %%
# ------------------------
# Paths & data load
# ------------------------
directory = "C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire"
figure_dir = os.path.join(directory, "out/fig")
sim_dat_path = os.path.join(directory, "data/rdat/mc_res_agg_random_sal0.csv")
obs_dat_path = os.path.join(directory, "data/wildfires/wf_observed.csv")

sim_dat = pd.read_csv(sim_dat_path)
obs_dat = pd.read_csv(obs_dat_path)

# ------------------------
# Select ignition cells & subset year
# ------------------------
# Use up to 500 ignition cells (or fewer if fewer available)
max_M = 1000
all_ign = sim_dat["ignition_id"].unique()
if len(all_ign) > max_M:
    rdm_ign_cell_idx = np.random.default_rng(0).choice(all_ign, size=max_M, replace=False)
else:
    rdm_ign_cell_idx = all_ign

sim_dat = sim_dat[sim_dat["ignition_id"].isin(rdm_ign_cell_idx)]
sim_dat = sim_dat[sim_dat["year"] < 2042]  # subset years for speed

# new column sim_id from ignition_id and year
sim_dat["sim_id"] = sim_dat["ignition_id"].astype(str) + "_" + sim_dat["year"].astype(str)

# ------------------------
# Build tensors
# ------------------------
times = jnp.array(sorted(sim_dat["time"].unique()))                      # (T,)
sim_matrix = jnp.array(                                                   # (T, M)
    sim_dat.pivot(index="time", columns="sim_id", values="burned_land_ha").values
)

print("sim_matrix shape (T, M):", sim_matrix.shape)
print("times shape (T,):", times.shape)
print("observed_areas shape raw (N,):", obs_dat.shape)

# check min max ranges of obs and sim
print(f"Observed areas: min={obs_dat['area_ha'].min():.2f}, max={obs_dat['area_ha'].max():.2f}")
print(f"Simulated areas: min={sim_matrix.min():.2f}, max={sim_matrix.max():.2f}")


# share obs fire smaller 4
obs_smaller_4 = obs_dat[obs_dat["area_ha"] < 4]
print(f"Observed fires smaller than 4 ha: {obs_smaller_4.shape[0]}")
# share obs fire larger 1000
obs_larger_1000 = obs_dat[obs_dat["area_ha"] >= 1000]
print(f"Observed fires larger than 1000 ha: {obs_larger_1000.shape[0]}")

# remove obs fires smaller 4 ha for calibration
obs_dat = obs_dat[obs_dat["area_ha"] >= 4]
obs_dat = obs_dat[obs_dat["area_ha"] < 1000]  

observed_areas = jnp.array(obs_dat["area_ha"].values)                    # (N,)
print("observed_areas shape after filter (N,):", observed_areas.shape)

all_ids = sim_dat["sim_id"].unique()

# wide matrix: rows = time, columns = ignition cell
sim_matrix = sim_dat.pivot(index="time", columns="sim_id", values="burned_land_ha")
sim_matrix = sim_matrix.values  # T x M

T, M = sim_matrix.shape
print("Simulator grid:", sim_matrix.shape, "times:", len(times))


# %%
def interp_area_at_time(duration, times, sim_matrix, ign_idx):
    if duration <= times[0]:
        return sim_matrix[0, ign_idx]
    if duration >= times[-1]:
        return sim_matrix[-1, ign_idx]

    hi = np.searchsorted(times, duration, side='right')
    lo = hi - 1
    t0, t1 = times[lo], times[hi]
    w = (duration - t0) / (t1 - t0)

    return (1 - w) * sim_matrix[lo, ign_idx] + w * sim_matrix[hi, ign_idx]

def sample_gev_durations(xi, mu, sigma, n):
    """
    duration ~ GEV(xi, mu, sigma)
    SciPy uses c = -xi
    """
    c = -xi
    # raw samples
    raw = genextreme.rvs(c=c, loc=mu, scale=sigma, size=n)
    # durations must be > 0
    return raw[raw > 0]

def simulate_areas_gev(xi, mu, sigma, n_samples=30000):
    """
    Draw durations from GEV (xi, mu, sigma),
    resample until we get at least n_samples positive values.
    Then sample ignitions and compute areas.
    """
    # SciPy shape is c = -xi
    c = -xi
    
    # continuously sample until enough positive durations
    durs = []
    batch_size = n_samples  # large batch to reduce loops
    
    while len(durs) < n_samples:
        raw = genextreme.rvs(c=c, loc=mu, scale=sigma, size=batch_size)
        raw = raw[raw > 0]  # keep only valid durations
        durs.extend(raw)
    
    durs = np.array(durs[:n_samples])  # truncate to exactly n_samples
    
    # ignition sampling
    ign_idxs = np.random.randint(0, M, size=n_samples)

    areas = np.array([
        interp_area_at_time(durs[j], times, sim_matrix, ign_idxs[j])
        for j in range(n_samples)
    ])

    return areas[np.isfinite(areas) & (areas > 0)]


obs_areas = obs_dat["area_ha"].values
obs_areas = obs_areas[np.isfinite(obs_areas) & (obs_areas > 0)]

def calibration_loss(params):
    xi, mu, sigma = params
    if sigma <= 0:  # invalid
        return np.inf

    sim = simulate_areas_gev(xi, mu, sigma)

    # use log scale for stability
    return wasserstein_distance(np.log1p(obs_areas), np.log1p(sim))


# %%
xis    = np.linspace(-0.3, 0.3, 7)   # shape parameter
mus    = np.linspace(50, 400, 8)     # location parameter
sigmas = np.linspace(20, 200, 8)     # scale parameter

results = []
for xi in xis:
    for mu in mus:
        for sigma in sigmas:
            loss = calibration_loss((xi, mu, sigma))
            print(xi, mu, sigma, loss)
            results.append((loss, xi, mu, sigma))

best_loss, best_xi, best_mu, best_sigma = sorted(results)[0]
print("Grid best:", best_xi, best_mu, best_sigma, "loss:", best_loss)

res = minimize(
    calibration_loss,
    x0=[best_xi, best_mu, best_sigma],
    bounds=[(-0.5, 0.5),   # xi
            (10,  600),    # mu
            (5,   400)],   # sigma
    method='Nelder-Mead',
    options={'maxiter': 200}
)

xi_star, mu_star, sigma_star = res.x
t_max = mu_star - sigma_star / xi_star
print("Optimized GEV params:")
print("xi     =", xi_star)
print("mu     =", mu_star)
print("sigma  =", sigma_star)
print("Theoretical max duration from GEV:", t_max, "vs simulation max time:", times[-1])
print("Bounds:", res.bounds if hasattr(res, 'bounds') else 'see code')



# %% Simulate areas with calibrated GEV
sim_pred = simulate_areas_gev(xi_star, mu_star, sigma_star, n_samples=50000)

# --------------------------------------------------
# CLEAN INPUT VECTORS
# --------------------------------------------------
obs = obs_areas[np.isfinite(obs_areas) & (obs_areas > 0)]
sim = sim_pred[np.isfinite(sim_pred) & (sim_pred > 0)]

log_obs = np.log10(obs)
log_sim = np.log10(sim)

# KDEs
kde_obs = gaussian_kde(log_obs)
kde_sim = gaussian_kde(log_sim)

log_grid = np.linspace(min(log_obs.min(), log_sim.min()),
                       max(log_obs.max(), log_sim.max()), 400)
dens_obs = kde_obs(log_grid)
dens_sim = kde_sim(log_grid)
area_grid = 10**log_grid

# ECDF helper
def ecdf(x):
    xs = np.sort(x)
    ys = np.linspace(0, 1, len(xs))
    return xs, ys

obs_x, obs_y = ecdf(obs)
sim_x, sim_y = ecdf(sim)

# --------------------------------------------------
# STATISTICS
# --------------------------------------------------
wd = wasserstein_distance(log_obs, log_sim)
ks = ks_2samp(log_obs, log_sim)
ad = anderson_ksamp([log_obs, log_sim])

p = kde_obs(log_grid)
q = kde_sim(log_grid)
kl = entropy(p, q)

# --------------------------------------------------
# FIGURE
# --------------------------------------------------
fig, axs = plt.subplots(3, 2, figsize=(16, 10))
sns.set_style("whitegrid")

# --------------------------------------------------
# Title block with estimated GEV parameters on top
# --------------------------------------------------
fig.suptitle(
    f"GEV Calibration Diagnostics\n"
    f"Estimated parameters: ξ = {xi_star:.3f}, μ = {mu_star:.3f}, σ = {sigma_star:.3f}",
    fontsize=18,
    y=1.04
)

# --------------------------------------------------
# (1) Log-space KDE → real area axis
# --------------------------------------------------
axs[0,0].plot(area_grid, dens_obs, label="Observed (log-KDE)", linewidth=2)
axs[0,0].plot(area_grid, dens_sim, label="Simulated (log-KDE)", linewidth=2)
axs[0,0].set_xscale("log")
axs[0,0].set_title("1. Log-KDE (density in log-space)")
axs[0,0].set_xlabel("Burned area (ha)")
axs[0,0].legend()

# --------------------------------------------------
# (2) ECDF
# --------------------------------------------------
axs[1,0].step(obs_x, obs_y, where="post", label="Observed")
axs[1,0].step(sim_x, sim_y, where="post", label="Simulated")
axs[1,0].set_xscale("log")
axs[1,0].set_title("2. Empirical CDF")
axs[1,0].set_xlabel("Burned area (ha)")
axs[1,0].legend()

# --------------------------------------------------
# (3) Q–Q Plot
# --------------------------------------------------
q = np.linspace(0,1,len(log_obs))
qq_obs = np.quantile(log_obs, q)
qq_sim = np.quantile(log_sim, q)

axs[0,1].plot(qq_obs, qq_sim, lw=2)
axs[0,1].plot([qq_obs.min(), qq_obs.max()],
              [qq_obs.min(), qq_obs.max()], 'k--')
axs[0,1].set_xlabel("Observed log10(area)")
axs[0,1].set_ylabel("Simulated log10(area)")
axs[0,1].set_title("3. Q–Q plot (log-space)")
axs[0,1].grid(True)

# --------------------------------------------------
# (4) P–P plot (correct, common grid)
# --------------------------------------------------
pgrid = np.linspace(0.001, 0.999, 300)
q_obs = np.quantile(log_obs, pgrid)
cdf_sim_at_qobs = np.interp(q_obs, np.sort(log_sim), np.linspace(0,1,len(log_sim)))

axs[1,1].plot([0,1],[0,1],'k--', linewidth=1)
axs[1,1].scatter(pgrid, cdf_sim_at_qobs, s=12, alpha=0.5)
axs[1,1].set_title("4. P–P Plot (Common Grid)")
axs[1,1].set_xlabel("Observed CDF")
axs[1,1].set_ylabel("Simulated CDF")
axs[1,1].set_xlim(0,1)
axs[1,1].set_ylim(0,1)
axs[1,1].grid(True)

# --------------------------------------------------
# (5) Log-binned histogram
# --------------------------------------------------
bins = np.logspace(np.log10(min(obs.min(), sim.min())),
                   np.log10(max(obs.max(), sim.max())), 40)
axs[2,0].hist(obs, bins=bins, alpha=0.5, label="Observed", density=True)
axs[2,0].hist(sim, bins=bins, alpha=0.5, label="Simulated", density=True)
axs[2,0].set_xscale("log")
axs[2,0].set_title("5. Log-binned histogram")
axs[2,0].set_xlabel("Burned area (ha)")
axs[2,0].legend()

# --------------------------------------------------
# (6) Residuals
# --------------------------------------------------
res = qq_sim - qq_obs
axs[2,1].plot(qq_obs, res, lw=2)
axs[2,1].axhline(0, color='k', linestyle='--')
axs[2,1].set_title("6. Residuals (log-space)")
axs[2,1].set_xlabel("Observed quantiles (log10 area)")
axs[2,1].set_ylabel("Residual = Sim - Obs")

plt.tight_layout()

# --------------------------------------------------
# PRINT STATS BELOW PLOT (also render on figure)
# --------------------------------------------------
stats_txt = (
    f"GEV: ξ={xi_star:.3f}, μ={mu_star:.1f}, σ={sigma_star:.1f}    "
    f"Wasserstein (log): {wd:.3f}    "
    f"KS: {ks.statistic:.3f} (p={ks.pvalue:.3f})    "
    f"AD: {ad.statistic:.3f} (p={getattr(ad, 'significance_level', float('nan')):.3f})    "
    f"KL (log-KDE): {kl:.3f}"
)

# place centered under the suptitle (adjust y coordinate if needed)
fig.text(0.5, 0.98, stats_txt, ha='center', va='top', fontsize=11)

# also keep console output
print("\n===========================================")
print("         GOODNESS-OF-FIT STATISTICS        ")
print("===========================================")
print(stats_txt)
print("===========================================\n")

# save plot in figure directory
fig.savefig(os.path.join(figure_dir, "gev_calibration_diagnostics.png"), dpi=300, bbox_inches='tight')

plt.show()




# %% # Export calibrated GEV parameters
pd.DataFrame({
    "xi_star":    [xi_star],
    "mu_star":    [mu_star],
    "sigma_star": [sigma_star]
}).to_csv("calibrated_gev_duration_params.csv", index=False)

# %%
# print current working directory
import os
print("Current working directory:", os.getcwd())
# %%
