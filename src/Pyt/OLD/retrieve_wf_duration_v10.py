# %%
import numpy as np
import jax
import jax.numpy as jnp
from jax.scipy.stats import beta as jbeta
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
jax.config.update("jax_enable_x64", True)

# %%
# ------------------------
# Paths & data load
# ------------------------
directory = "C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire"
sim_dat_path = os.path.join(directory, "data/rdat/mc_res_agg_random_sal0.csv")
obs_dat_path = os.path.join(directory, "data/wildfires/wf_observed.csv")

sim_dat = pd.read_csv(sim_dat_path)
obs_dat = pd.read_csv(obs_dat_path)

# ------------------------
# Select ignition cells & subset year
# ------------------------
# Use up to 500 ignition cells (or fewer if fewer available)
max_M = 800
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



# %%

# ---------------------------------------------
import numpy as np, jax, jax.numpy as jnp
from jax.scipy.stats import beta as jbeta

# ---------- core utilities ----------
def interp_cols_at_t(durations, times, sim_rows):
    """Linear-interpolate (row-wise) sim_rows(T,M) at durations(N,) -> (N,M)."""
    idx_hi = jnp.searchsorted(times, durations, side='right')
    idx_lo = jnp.clip(idx_hi - 1, 0, times.shape[0] - 2)
    t0, t1 = times[idx_lo], times[idx_lo + 1]
    a = (durations - t0) / (t1 - t0)
    lo, hi = sim_rows[idx_lo, :], sim_rows[idx_lo + 1, :]
    return (1 - a)[:, None] * lo + a[:, None] * hi  # (N,M)

def logsumexp_unweighted(logp, axis=-1):
    m = jnp.max(logp, axis=axis, keepdims=True)
    return jnp.log(jnp.mean(jnp.exp(logp - m), axis=axis)) + jnp.squeeze(m, axis=axis)

def _trunc_lognorm_logpdf(logA_i, mu_tm, sigma_log, aL):
    """Left-truncated Normal on log-scale: logA_i ~ N(mu_tm, sigma), truncated at aL=log(4)."""
    z = (aL - mu_tm) / sigma_log
    # log(1 - Phi(z)) using erf for stability
    logZ = jnp.log(1 - 0.5*(1 + jax.lax.erf(z / jnp.sqrt(2))))
    logpdf = -0.5 * jnp.square((logA_i - mu_tm) / sigma_log) \
             - jnp.log(sigma_log) - 0.5*jnp.log(2*jnp.pi)
    return logpdf - logZ

# ---------- inner: compute per-fire log posterior over t for given (tau, phi) ----------
def _logpost_over_t(times, sim_matrix, logA, sigma_log, alpha_beta, beta_beta, truncation, tau, phi):
    eps = 1e-12
    T, M = sim_matrix.shape
    # warp: evaluate sims at t' = clip(phi * t)
    t_eff = jnp.clip(times * phi, times[0], times[-1])
    sim_warp_TM = interp_cols_at_t(t_eff, times, sim_matrix)  # (T,M)
    # scale tau
    log_sim_TM = jnp.log(jnp.maximum(sim_warp_TM, eps)) + jnp.log(tau + eps)

    aL = jnp.log(truncation)
    # prior over t (on original grid)
    u_grid = (times - times[0]) / (times[-1] - times[0])
    log_prior_t = jbeta.logpdf(u_grid, a=alpha_beta, b=beta_beta)  # (T,)

    def one_fire(logA_i):
        # truncated log-normal per ignition -> (T,M)
        logpdf_TM = _trunc_lognorm_logpdf(logA_i, log_sim_TM, sigma_log, aL)
        # mixture over ignitions (unweighted mean)
        loglik_T = logsumexp_unweighted(logpdf_TM, axis=1)  # (T,)
        return loglik_T + log_prior_t

    return jax.vmap(one_fire)(logA)  # (N,T)

# ---------- FAST calibration with (tau, phi) grid search ----------
def calibrate_durations_tau_phi(times, sim_matrix, observed_areas,
                                sigma_log=0.35, alpha_beta=1., beta_beta=2.,
                                truncation=4.0,
                                tau_grid=None, phi_grid=None):
    """
    Returns: post_mean_t, post_low, post_high, post, best_tau, best_phi, mask
    """
    eps = 1e-12
    T, M = sim_matrix.shape

    # keep only >= truncation
    mask = np.asarray(observed_areas) >= truncation
    obs = np.asarray(observed_areas)[mask]
    logA = jnp.log(jnp.maximum(jnp.asarray(obs), eps))  # (N,)

    if tau_grid is None:
        tau_grid = jnp.linspace(0.5, 2.4, 20)
    if phi_grid is None:
        phi_grid = jnp.linspace(0.05, 0.5, 10)  # <1 slows early growth

    # score each (tau,phi) by total log-evidence across fires:
    # sum_i logsumexp_t [ loglik_i(t; tau,phi) ]
    best = (-jnp.inf, None, None, None)  # (score, tau, phi, logpost_NT)

    for tau in tau_grid:
        for phi in phi_grid:
            logpost_NT = _logpost_over_t(times, sim_matrix, logA,
                                         sigma_log, alpha_beta, beta_beta,
                                         truncation, tau, phi)
            # log-evidence per fire: logsumexp over t
            def lse(v):
                m = jnp.max(v)
                return m + jnp.log(jnp.sum(jnp.exp(v - m)))
            score = jnp.sum(jax.vmap(lse)(logpost_NT))
            if score > best[0]:
                best = (score, float(tau), float(phi), logpost_NT)

    best_tau, best_phi = best[1], best[2]
    logpost_NT = best[3]

    # normalize to get posteriors over t
    def normalize_row(v):
        m = jnp.max(v)
        return jnp.exp(v - (m + jnp.log(jnp.sum(jnp.exp(v - m)))))
    post = jax.vmap(normalize_row)(logpost_NT)  # (N,T)

    # posterior summaries
    t_grid = times
    post_mean_t = jnp.sum(post * t_grid[None, :], axis=1)
    cdf = jnp.cumsum(post, axis=1)
    cdf = cdf / cdf[:, -1][:, None]

    def qfun(c, row):
        idx = jnp.searchsorted(row, c, side="left")
        idx = jnp.clip(idx, 0, T-1)
        return t_grid[idx]
    post_low  = jax.vmap(lambda r: qfun(0.025, r))(cdf)
    post_high = jax.vmap(lambda r: qfun(0.975, r))(cdf)

    return (np.asarray(post_mean_t), np.asarray(post_low), np.asarray(post_high),
            np.asarray(post), best_tau, best_phi, mask)

# ---------- posterior predictive updated to use tau & phi ----------
def posterior_predictive(times, sim_matrix, duration_samples, tau=1.0, phi=1.0):
    """
    duration_samples: (S,N) sampled *on the original t-grid*; we evaluate at t_eff=phi*t
    and multiply by tau.
    """
    eps = 1e-12
    T, M = sim_matrix.shape
    S, N = duration_samples.shape

    # warp durations
    t_eff = np.clip(phi * duration_samples, float(times[0]), float(times[-1]))
    sim_interp = interp_cols_at_t(jnp.array(t_eff).reshape(-1), times, sim_matrix)  # (S*N, M)

    rng = np.random.default_rng(0)
    idx = rng.integers(0, M, size=sim_interp.shape[0])  # uniform over ignitions (replace with weights if available)
    chosen = sim_interp[jnp.arange(sim_interp.shape[0]), idx]
    return np.asarray((tau * jnp.maximum(chosen, eps)).reshape(S, N))

# ---------- P–P plot ---------- 
def pp_plot(observed, simulated, ax=None): 
    observed = np.sort(observed) 
    sim_q = np.quantile(simulated, np.linspace(0,1,len(observed))) 
    obs_q = np.quantile(observed, np.linspace(0,1,len(observed))) 
    if ax is None: 
        fig, ax = plt.subplots(figsize=(5,5)) 
    ax.plot(obs_q, sim_q, lw=2) 
    lim = [min(obs_q.min(), sim_q.min()), max(obs_q.max(), sim_q.max())] 
    ax.plot(lim, lim, 'k--', lw=1) 
    ax.set_xscale('log'); ax.set_yscale('log') 
    ax.set_xlabel("Observed quantiles (ha)") 
    ax.set_ylabel("Posterior predictive quantiles (ha)") 
    ax.set_title("P–P plot: observed vs posterior predictive") 
    ax.grid(True, ls=":") 
    return ax


# %%
# Calibrate with τ and φ (still fast)
post_mean, post_low, post_high, post, best_tau, best_phi, mask = \
    calibrate_durations_tau_phi(times, sim_matrix, observed_areas,
                                sigma_log=0.35, alpha_beta=1., beta_beta=2.,
                                truncation=4.0)

print("Best tau:", best_tau, "Best phi:", best_phi)
print("Posterior mean duration (min):", post_mean.mean())

# %% ---------- posterior predictive checks ----------
# Posterior sampling of durations (S draws per fire)
rng = np.random.default_rng(123)
T = len(times); S = 100
post_np = (post.T / post.sum(axis=1)).T  # row-normalize (safety)
duration_samples = np.empty((S, post_np.shape[0]))
for i in range(post_np.shape[0]):
    idx = rng.choice(T, size=S, p=post_np[i])
    duration_samples[:, i] = times[idx]

# Posterior predictive with learned (tau, phi)
pp_draws = posterior_predictive(times, sim_matrix, duration_samples, tau=best_tau, phi=best_phi)
pp_all = pp_draws.ravel()
# ---------- plots ----------
obs = np.asarray(observed_areas)
obs_clean = obs[np.isfinite(obs) & (obs > 0)]
pp_clean = pp_all[np.isfinite(pp_all) & (pp_all > 0)]
plt.figure(figsize=(8,5))
sns.kdeplot(obs_clean, label="Observed", fill=True, alpha=0.4)
sns.kdeplot(pp_clean, label="Posterior predictive", fill=True, alpha=0.4)
plt.xscale("log")
plt.xlabel("Burned area [ha] (log scale)")
plt.ylabel("Density")
plt.title("Posterior predictive calibration with learned tau & phi")
plt.legend()
plt.tight_layout()
plt.show()

# ---------- P–P plot ----------
pp_plot(obs_clean, pp_clean)
plt.tight_layout()
plt.show()


# %% Export
import numpy as np, pandas as pd
from scipy.stats import genextreme

# pp_clean = vector of posterior predictive fire areas (ha), from your run
pp = np.array(pp_clean)
pp = pp[pp > 0]  # safety

mu = np.mean(np.log(pp))
sigma = np.std(np.log(pp))

# Fit GEV to fire areas (log scale recommended)
gev_params = genextreme.fit(np.log(pp))

shape, loc, scale = gev_params

params = pd.DataFrame({
    "distribution": ["lognormal", "gev_log"],
    "mu_log": [mu, np.nan],
    "sigma_log": [sigma, np.nan],
    "shape": [np.nan, shape],
    "loc": [np.nan, loc],
    "scale": [np.nan, scale]
})

params.to_csv("calibrated_fire_area_distribution.csv", index=False)

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import lognorm

# ----------------------------------------
# INPUT: posterior durations (physical t!)
# ----------------------------------------
dur = np.array(post_mean)   # or use posterior mean vector
dur = dur[np.isfinite(dur) & (dur > 0)]

# Fit LogNormal to durations
log_dur = np.log(dur)
mu_dur  = np.mean(log_dur)
sigma_dur = np.std(log_dur)

# SciPy parameterization: shape = sigma, scale = exp(mu)
shape = sigma_dur
scale = np.exp(mu_dur)

# PDF grid
x = np.linspace(dur.min(), dur.max(), 500)
pdf = lognorm.pdf(x, s=shape, scale=scale)

# Plot histogram + PDF
plt.figure(figsize=(8,5))
plt.hist(dur, bins=30, density=True, alpha=0.4, color='gray', label="Posterior durations")
plt.plot(x, pdf, 'r-', lw=2, label=f"Fitted LogNormal\nμ={mu_dur:.3f}, σ={sigma_dur:.3f}")
plt.xlabel("Duration (minutes)")
plt.ylabel("Density")
plt.title("Posterior Calibrated Durations + LogNormal Fit")
plt.legend()
plt.tight_layout()
plt.show()

# %%
phi = best_phi  # e.g. 0.2
dur_eff = phi * dur
dur_eff = dur_eff[dur_eff > 0]

# fit to warped durations
mu_eff  = np.mean(np.log(dur_eff))
sigma_eff = np.std(np.log(dur_eff))
shape_eff = sigma_eff
scale_eff = np.exp(mu_eff)

x2 = np.linspace(dur_eff.min(), dur_eff.max(), 500)
pdf2 = lognorm.pdf(x2, s=shape_eff, scale=scale_eff)

plt.figure(figsize=(8,5))
plt.hist(dur_eff, bins=30, density=True, alpha=0.4, color='steelblue', label="Effective durations φ·t")
plt.plot(x2, pdf2, 'k--', lw=2, label=f"Fitted LogNormal (φ·t)\nμ={mu_eff:.3f}, σ={sigma_eff:.3f}")
plt.xlabel("Effective duration (minutes)")
plt.ylabel("Density")
plt.title("Effective Durations After Time Warp φ")
plt.legend()
plt.tight_layout()
plt.show()

# %% Export effective duration parameters
import pandas as pd
params_dur = pd.DataFrame({
    "parameter": ["mu_log_eff", "sigma_log_eff", "shape_eff", "scale_eff"],
    "value": [mu_eff, sigma_eff, shape_eff, scale_eff]
})
params_dur.to_csv("calibrated_effective_duration_distribution.csv", index=False)

# %%
# print working directory
import os
print("Current working directory:", os.getcwd())
# %%
