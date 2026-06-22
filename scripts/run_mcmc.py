"""
run_mcmc.py — observable-based MCMC for the windowed-IDE model (P1).

Fits (beta0, z_star, sigma) [or +z_c, dz with --ndim 5] to synthetic
H(z)/D_M(z)/fsigma8(z) data generated from a known truth, using emcee.
Reports convergence (autocorrelation time, Gelman-Rubin R-hat across walkers),
parameter recovery, and writes a corner plot + summary CSV.

Usage:
    python scripts/run_mcmc.py --fit beta0                 # clean pipeline check
    python scripts/run_mcmc.py --fit beta0,z_star,sigma    # full (shows degeneracy)
    python scripts/run_mcmc.py --precision current ...      # current-precision null
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import emcee  # noqa: E402
import corner  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from mtp_cosmology import Params  # noqa: E402
from mtp_cosmology.likelihood import (  # noqa: E402
    log_posterior, make_synthetic_data, theta_to_params, observables,
    PRIOR_BOUNDS,
)

# Full injected truth (all 5 free params); fits use a subset, others fixed here.
TRUTH = Params(beta0=0.10, z_star=0.40, sigma=0.30, z_c=0.50, dz=0.20)
ALL_NAMES = ["beta0", "z_star", "sigma", "z_c", "dz"]


def gelman_rubin(chain_2d_per_walker: np.ndarray) -> np.ndarray:
    """R-hat per parameter. chain shape: (nwalkers, nsteps, ndim)."""
    m, n, ndim = chain_2d_per_walker.shape
    rhats = np.zeros(ndim)
    for d in range(ndim):
        x = chain_2d_per_walker[:, :, d]               # (m walkers, n steps)
        chain_means = x.mean(axis=1)                   # (m,)
        grand_mean = chain_means.mean()
        B = n * np.var(chain_means, ddof=1)            # between-chain
        W = np.mean(np.var(x, axis=1, ddof=1))         # within-chain
        var_hat = (n - 1) / n * W + B / n
        rhats[d] = np.sqrt(var_hat / W) if W > 0 else np.nan
    return rhats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default="beta0,z_star,sigma",
                    help="comma-separated subset of beta0,z_star,sigma,z_c,dz to vary; "
                         "the rest are fixed at the injected truth")
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--walkers", type=int, default=32)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--precision", default="forecast", choices=["current", "forecast"],
                    help="synthetic-data error level (see likelihood.PRECISION_PRESETS)")
    ap.add_argument("--tag", default="", help="optional output filename suffix")
    ap.add_argument("--fig-dir", default=os.path.join(os.path.dirname(__file__), "..", "figures"))
    ap.add_argument("--res-dir", default=os.path.join(os.path.dirname(__file__), "..", "results"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    fit = [s.strip() for s in args.fit.split(",") if s.strip()]
    assert all(f in ALL_NAMES for f in fit), f"fit names must be in {ALL_NAMES}"
    ndim = len(fit)
    idx = [ALL_NAMES.index(f) for f in fit]
    truth_full = TRUTH.as_array()[:5]
    truth_vec = truth_full[idx]
    names = fit
    tag = args.tag or "-".join(fit)

    probes = make_synthetic_data(TRUTH, seed=args.seed, precision=args.precision)
    n_data = sum(len(pr.z) for pr in probes)
    print(f"Synthetic data: {n_data} points across H/D_M/fsigma8; "
          f"fit={fit}; precision={args.precision}")

    # Expand a sampled subset to a full 5-vector (fixed params from truth).
    def expand(theta_sub):
        full = truth_full.copy()
        full[idx] = theta_sub
        return full

    def logpost_sub(theta_sub, probes):
        return log_posterior(expand(theta_sub), probes)

    # Initialize walkers around a deliberately-offset guess (NOT truth).
    offset = {"beta0": 0.05, "z_star": 0.60, "sigma": 0.20, "z_c": 0.40, "dz": 0.25}
    guess = np.array([offset[f] for f in fit])
    p0 = guess + 1e-3 * rng.standard_normal((args.walkers, ndim))
    p0 = np.abs(p0)

    sampler = emcee.EnsembleSampler(args.walkers, ndim, logpost_sub, args=(probes,))
    print(f"Running emcee: {args.walkers} walkers x {args.steps} steps (ndim={ndim}) ...")
    sampler.run_mcmc(p0, args.steps, progress=False)

    # Autocorrelation -> burn-in / thinning
    try:
        tau = sampler.get_autocorr_time(tol=0)
        tau_max = float(np.nanmax(tau))
    except Exception:
        tau_max = args.steps / 50.0
    burn = int(min(max(2 * tau_max, args.steps * 0.2), args.steps * 0.5))
    thin = max(int(tau_max / 2), 1)

    rhat = gelman_rubin(sampler.get_chain(discard=burn))
    flat = sampler.get_chain(discard=burn, thin=thin, flat=True)
    accept = float(np.mean(sampler.acceptance_fraction))

    print(f"\nConvergence: tau_max≈{tau_max:.1f}, burn={burn}, thin={thin}, "
          f"accept={accept:.2f}, {flat.shape[0]} samples")
    print(f"R-hat: " + ", ".join(f"{n}={r:.4f}" for n, r in zip(names, rhat)))

    # Recovery table
    mean = flat.mean(axis=0)
    std = flat.std(axis=0)
    print(f"\n{'param':<8}{'truth':>9}{'mean':>11}{'std':>9}{'pull(sigma)':>13}")
    print("-" * 50)
    rows = []
    for i, nm in enumerate(names):
        pull = (mean[i] - truth_vec[i]) / std[i] if std[i] > 0 else np.nan
        print(f"{nm:<8}{truth_vec[i]:>9.4f}{mean[i]:>11.4f}{std[i]:>9.4f}{pull:>13.2f}")
        rows.append({"param": nm, "truth": truth_vec[i], "mean": mean[i],
                     "std": std[i], "pull_sigma": pull, "rhat": rhat[i]})
    within_1sigma = all(abs(r["pull_sigma"]) < 1.0 for r in rows)
    print(f"\nAll params recovered within 1 sigma: {within_1sigma}")

    # Derived: w_eff(0) and CPL (w0, wa) at the posterior mean
    from mtp_cosmology import Hz, w_eff, cpl_fit  # noqa: E402
    z = np.linspace(0.001, 2.5, 500)
    _, _, rde = Hz(z, theta_to_params(expand(mean)))
    we = w_eff(z, rde)
    w0, wa = cpl_fit(z, we)
    print(f"\nDerived @ posterior mean: w_eff(0)={we[0]:.4f}, "
          f"w_min={we.min():.4f}, CPL (w0,wa)=({w0:.4f},{wa:.4f})")

    # Outputs
    os.makedirs(args.fig_dir, exist_ok=True)
    os.makedirs(args.res_dir, exist_ok=True)
    fig_path = os.path.join(args.fig_dir, f"mcmc_corner_{tag}_{args.precision}.png")
    if ndim == 1:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.hist(flat[:, 0], bins=40, color="#455A64", alpha=0.85)
        ax.axvline(truth_vec[0], color="#2196F3", lw=2, label="truth")
        ax.axvline(mean[0], color="#FF5722", ls="--", lw=2, label="posterior mean")
        ax.set_xlabel(names[0]); ax.set_ylabel("samples"); ax.legend()
        ax.set_title(f"{names[0]} = {mean[0]:.4f} +/- {std[0]:.4f}")
    else:
        fig = corner.corner(flat, labels=names, truths=list(truth_vec),
                            show_titles=True, title_fmt=".4f")
    fig.savefig(fig_path, dpi=130, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(flat, columns=names).to_csv(
        os.path.join(args.res_dir, f"mcmc_chain_{tag}_{args.precision}.csv"), index=False)
    summary = pd.DataFrame(rows)
    summary["w_eff0"] = we[0]
    summary["w0"] = w0
    summary["wa"] = wa
    summary["accept_frac"] = accept
    summary.to_csv(os.path.join(args.res_dir, f"mcmc_summary_{tag}_{args.precision}.csv"), index=False)

    print(f"\nWrote: {os.path.relpath(fig_path)}, "
          f"results/mcmc_summary_{tag}_{args.precision}.csv, results/mcmc_chain_{tag}_{args.precision}.csv")


if __name__ == "__main__":
    main()
