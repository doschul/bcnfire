# %%
import os
import numpy as np
import jax
import jax.numpy as jnp
from jax.scipy.stats import t as jstudent_t
from jax.scipy.stats import beta as jbeta
from jax.scipy.special import logsumexp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import numpyro
import numpyro.distributions as dist
from numpyro.infer import NUTS, MCMC

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
max_M = 100
all_ign = sim_dat["ignition_id"].unique()
if len(all_ign) > max_M:
    rdm_ign_cell_idx = np.random.default_rng(0).choice(all_ign, size=max_M, replace=False)
else:
    rdm_ign_cell_idx = all_ign

sim_dat = sim_dat[sim_dat["ignition_id"].isin(rdm_ign_cell_idx)]
# Optional: subset year if desired
# sim_dat = sim_dat[sim_dat["year"] == 2040]

# new column sim_id from ignition_id and year
sim_dat["sim_id"] = sim_dat["ignition_id"].astype(str) + "_" + sim_dat["year"].astype(str)

# ------------------------
# Build tensors
# ------------------------
observed_areas = jnp.array(obs_dat["area_ha"].values)  # (N,)
times = jnp.array(sorted(sim_dat["time"].unique()))    # (T,)
sim_matrix = jnp.array(
    sim_dat.pivot(index="time", columns="sim_id", values="burned_land_ha").values
)  # (T, M)

print("sim_matrix shape (T, M):", sim_matrix.shape)
print("times shape (T,):", times.shape)
print("observed_areas shape (N,):", observed_areas.shape)

# %%
# ---------- core utilities ----------
def interp_cols_at_t(durations, times, sim_rows):
    """Linear-interpolate (row-wise) sim_rows(T,M) at durations(N,) -> (N,M)."""
    idx_hi = jnp.searchsorted(times, durations, side='right')
    idx_lo = jnp.clip(idx_hi - 1, 0, times.shape[0] - 2)
    t0, t1 = times[idx_lo], times[idx_lo + 1]
    a = (durations - t0) / (t1 - t0)
    lo, hi = sim_rows[idx_lo, :], sim_rows[idx_lo + 1, :]
    return (1 - a)[:, None] * lo + a[:, None] * hi  # (N,M)

def logmeanexp(x, axis=-1):
    """log( mean(exp(x), axis) )"""
    m = jnp.max(x, axis=axis, keepdims=True)
    return jnp.log(jnp.mean(jnp.exp(x - m), axis=axis)) + jnp.squeeze(m, axis=axis)

# %%
# ---------- MODEL (Student-t outcome; integrate over time grid & sim ensemble) ----------
def model_student_t(times, sim_matrix, observed_areas, alpha_beta=1., beta_beta=2.):
    """
    Observed log-areas are modeled as a mixture over:
      - unweighted sim ensemble at each time
      - discrete prior over times given by a Beta(a,b) on normalized grid

    For a given fire i:
      logA_i | t ~ mixture_m ( StudentT(df=nu_t, loc=log(sim[t,m]), scale=sigma_log) )

    We marginalize over t (and over ensemble m) analytically via log-mean-exp,
    so there are no discrete latents to sample.
    """
    eps = 1e-12
    T, M = sim_matrix.shape
    N = observed_areas.shape[0]

    # parameters
    sigma_log = numpyro.sample("sigma_log", dist.HalfNormal(1.0))  # wide but reasonable
    raw_nu = numpyro.sample("raw_nu", dist.Normal(0., 1.0))
    # enforce df > 2 (finite variance) and flexible tails
    nu_t = 2.0 + jax.nn.softplus(raw_nu)
    numpyro.deterministic("nu_t", nu_t)

    # precompute
    log_sim = jnp.log(jnp.maximum(sim_matrix, eps))  # (T,M)
    logA = jnp.log(jnp.maximum(observed_areas, eps)) # (N,)

    # prior over time grid (normalize across T to form a discrete prior)
    u_grid = (times - times[0]) / jnp.maximum(times[-1] - times[0], 1e-9)
    log_prior_t = jbeta.logpdf(u_grid, a=alpha_beta, b=beta_beta)  # (T,)
    log_prior_t = log_prior_t - logsumexp(log_prior_t)             # normalize across grid

    # per-fire marginal log-likelihood
    def fire_log_marginal(logA_i):
        # Student-t logpdf across T x M
        z = (logA_i - log_sim) / sigma_log          # (T,M)
        logpdf = jstudent_t.logpdf(z, df=nu_t) - jnp.log(sigma_log)  # (T,M)

        # mixture over M (unweighted)
        logpdf_mix_t = logmeanexp(logpdf, axis=1)   # (T,)

        # integrate over t with discrete prior
        return logsumexp(logpdf_mix_t + log_prior_t)  # scalar

    # vectorize over fires
    total_loglik = jnp.sum(jax.vmap(fire_log_marginal)(logA))
    numpyro.factor("loglik", total_loglik)

# %%
# ---------- run MCMC to infer (nu_t, sigma_log) ----------
rng_key = jax.random.PRNGKey(0)
kernel = NUTS(model_student_t, target_accept_prob=0.9)
mcmc = MCMC(kernel, num_warmup=100, num_samples=200, num_chains=1, progress_bar=True)
mcmc.run(rng_key, times=times, sim_matrix=sim_matrix, observed_areas=observed_areas)
mcmc.print_summary()

posterior = mcmc.get_samples()
nu_t_samps = np.asarray(posterior["nu_t"])
sigma_samps = np.asarray(posterior["sigma_log"])
nu_t_mean = float(nu_t_samps.mean())
sigma_mean = float(sigma_samps.mean())
print(f"Posterior means — nu_t: {nu_t_mean:.2f}, sigma_log: {sigma_mean:.3f}")

# %%
# ---------- Duration posterior (using posterior-mean parameters) ----------
def durations_posterior(times, sim_matrix, observed_areas,
                        sigma_log, nu_t, alpha_beta=1., beta_beta=2.):
    """
    Compute posterior p(t | data, params) on the discrete time grid for each fire.
    Returns: post_mean_t (N,), post_low (N,), post_high (N,), post (N,T)
    """
    eps = 1e-12
    T, M = sim_matrix.shape
    N = observed_areas.shape[0]
    log_sim = jnp.log(jnp.maximum(sim_matrix, eps))        # (T,M)
    logA = jnp.log(jnp.maximum(observed_areas, eps))       # (N,)

    # time prior (normalized discrete version of Beta(a,b) on grid)
    u_grid = (times - times[0]) / jnp.maximum(times[-1] - times[0], 1e-9)
    log_prior_t = jbeta.logpdf(u_grid, a=alpha_beta, b=beta_beta)  # (T,)
    log_prior_t = log_prior_t - logsumexp(log_prior_t)             # (T,)

    # Fire-wise loglik over t (mix over M)
    def loglik_over_t(logA_i):
        z = (logA_i - log_sim) / sigma_log                  # (T,M)
        logpdf = jstudent_t.logpdf(z, df=nu_t) - jnp.log(sigma_log)
        logpdf_mix_t = logmeanexp(logpdf, axis=1)           # (T,)
        return logpdf_mix_t                                  # (T,)

    loglik_all = jax.vmap(loglik_over_t)(logA)               # (N,T)
    log_post = loglik_all + log_prior_t[None, :]             # (N,T)

    # normalize rows
    def normalize_row(v):
        return v - logsumexp(v)
    log_post_norm = jax.vmap(normalize_row)(log_post)        # (N,T)
    post = jnp.exp(log_post_norm)                            # (N,T)

    # posterior summaries
    post_mean_t = jnp.sum(post * times[None, :], axis=1)     # (N,)

    cdf = jnp.cumsum(post, axis=1)
    cdf = cdf / cdf[:, -1][:, None]

    def qfun(c, row):
        idx = jnp.searchsorted(row, c, side="left")
        idx = jnp.clip(idx, 0, T-1)
        return times[idx]

    post_low  = jax.vmap(lambda r: qfun(0.025, r))(cdf)
    post_high = jax.vmap(lambda r: qfun(0.975, r))(cdf)

    return np.asarray(post_mean_t), np.asarray(post_low), np.asarray(post_high), np.asarray(post)

post_mean, post_low, post_high, post = durations_posterior(
    times, sim_matrix, observed_areas,
    sigma_log=sigma_mean, nu_t=nu_t_mean,
    alpha_beta=1., beta_beta=2.
)

print(f"Posterior mean duration (avg across fires): {post_mean.mean():.1f} min")

# %%
# ---------- posterior predictive ----------
def posterior_predictive(times, sim_matrix, duration_samples, rng_seed=0):
    """Draw simulated areas at sampled durations."""
    eps = 1e-12
    T, M = sim_matrix.shape
    S, N = duration_samples.shape
    flat_dur = jnp.array(duration_samples).reshape(-1)
    sim_interp = interp_cols_at_t(flat_dur, times, sim_matrix)  # (S*N,M)
    rng = np.random.default_rng(rng_seed)
    idx = rng.integers(0, M, size=sim_interp.shape[0])
    chosen = sim_interp[jnp.arange(sim_interp.shape[0]), idx]
    return np.asarray(chosen.reshape(S, N))

# Draw posterior durations (1 sample per fire for simplicity)
rng = np.random.default_rng(123)
Tlen = len(times)
# sample durations per-fire from its posterior over the grid
# here we sample 800 draws (S=800)
S = 800
# build categorical sampler per fire using its posterior probs over grid
post_probs = post / post.sum(axis=1, keepdims=True)  # (N,T)
cat_idx = np.vstack([
    rng.choice(Tlen, size=S, p=post_probs[i, :]) for i in range(post_probs.shape[0])
]).T  # (S,N) indices into time grid
duration_samples = np.asarray(times)[cat_idx]        # (S,N)

pp_draws = posterior_predictive(times, sim_matrix, duration_samples, rng_seed=123)
pp_all = pp_draws.flatten()

# %%
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
plt.title("Posterior predictive calibration (Student-t likelihood; nu inferred)")
plt.legend()
plt.tight_layout()
plt.show()

# %%
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

pp_plot(obs_clean, pp_clean)
plt.tight_layout()
plt.show()
