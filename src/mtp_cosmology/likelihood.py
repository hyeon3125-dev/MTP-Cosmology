"""
likelihood.py — observable-based likelihood for the windowed-IDE model (P1).

Supersedes legacy/advanced_mcmc_analysis_toyw.py, which fit a standalone
w(z) = -1 - beta0*logistic_pdf(z) shape unrelated to the model output. Here the
likelihood is driven by the *actual* model observables:

    H(z)    — expansion rate           (cosmic chronometers / BAO)
    D_M(z)  — comoving transverse dist. (BAO)
    f*sigma8(z) — structure growth      (RSD)

chi^2 is summed over the three probes. Free parameters default to the 3-param
set (beta0, z_star, sigma); the 5-param set (+z_c, dz) is also supported.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .model import Params, Hz, D_M, growth_fs8

# Parameter order used by the sampler (3- or 5-D).
PARAM_NAMES_3 = ["beta0", "z_star", "sigma"]
PARAM_NAMES_5 = ["beta0", "z_star", "sigma", "z_c", "dz"]

# Flat-prior bounds. Physically motivated, deliberately wide.
PRIOR_BOUNDS = {
    "beta0": (0.0, 1.0),     # coupling amplitude > 0
    "z_star": (0.0, 2.0),    # window center within observable range
    "sigma": (0.05, 1.5),    # window width
    "z_c": (0.0, 2.0),       # sign-switch redshift
    "dz": (0.05, 1.0),       # switch smoothing
}


def theta_to_params(theta) -> Params:
    """Map a sampler vector (len 3 or 5) to a model Params."""
    theta = np.asarray(theta, dtype=float)
    if theta.size == 3:
        return Params(beta0=theta[0], z_star=theta[1], sigma=theta[2])
    if theta.size == 5:
        return Params(beta0=theta[0], z_star=theta[1], sigma=theta[2],
                      z_c=theta[3], dz=theta[4])
    raise ValueError(f"theta must have length 3 or 5, got {theta.size}")


def param_names(ndim: int):
    return PARAM_NAMES_3 if ndim == 3 else PARAM_NAMES_5


def observables(z_data, p: Params, n_grid: int = 400):
    """
    Model H(z), D_M(z), f*sigma8(z) at the requested data redshifts.

    Integrates the background on a dense grid from 0 to max(z_data) and
    interpolates, so D_M (a cumulative integral from 0) is accurate even when
    z_data is sparse/unordered.
    """
    z_data = np.atleast_1d(np.asarray(z_data, dtype=float))
    z_max = max(float(z_data.max()), 0.01)
    grid = np.linspace(1e-4, z_max * 1.001, n_grid)

    H_grid, rho_c, _ = Hz(grid, p)
    DM_grid = D_M(grid, H_grid)
    fs8_grid = growth_fs8(grid, H_grid, rho_c)

    return {
        "H": np.interp(z_data, grid, H_grid),
        "D_M": np.interp(z_data, grid, DM_grid),
        "fs8": np.interp(z_data, grid, fs8_grid),
    }


@dataclass
class ProbeData:
    """One observable's data vector: redshifts, measured values, 1-sigma errors."""
    z: np.ndarray
    value: np.ndarray
    sigma: np.ndarray
    key: str  # "H" | "D_M" | "fs8"


def log_prior(theta) -> float:
    names = param_names(np.asarray(theta).size)
    for name, val in zip(names, np.asarray(theta, dtype=float)):
        lo, hi = PRIOR_BOUNDS[name]
        if not (lo < val < hi):
            return -np.inf
    return 0.0


def log_likelihood(theta, probes) -> float:
    """Gaussian chi^2 summed over probes. Returns -inf on integration failure.

    Integrates the background ONCE on a grid covering all probe redshifts, then
    interpolates each probe — ~3x faster than per-probe integration.
    """
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    p = theta_to_params(theta)
    z_all = np.concatenate([pr.z for pr in probes])
    z_max = max(float(z_all.max()), 0.01)
    grid = np.linspace(1e-4, z_max * 1.001, 250)
    try:
        H_grid, rho_c, _ = Hz(grid, p)
        DM_grid = D_M(grid, H_grid)
        fs8_grid = growth_fs8(grid, H_grid, rho_c)
        model_grid = {"H": H_grid, "D_M": DM_grid, "fs8": fs8_grid}
        chi2 = 0.0
        for pr in probes:
            model = np.interp(pr.z, grid, model_grid[pr.key])
            chi2 += np.sum(((pr.value - model) / pr.sigma) ** 2)
    except (RuntimeError, FloatingPointError, ValueError):
        return -np.inf
    if not np.isfinite(chi2):
        return -np.inf
    return -0.5 * chi2


def log_posterior(theta, probes) -> float:
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, probes)


# ── Synthetic data at realistic redshift nodes / error scales ────────────────
# Nodes/errors loosely mimic current compilations (cosmic chronometers, DESI/
# BOSS BAO, RSD). They are illustrative, not a literal data replica.
#
# Two precision presets (see worklog / paper):
#   "current"  — present-day fractional errors. The windowed-IDE signal at
#                beta0~0.1 is a sub-percent effect (H 0.6%, D_M 0.25%, fsigma8
#                1.3% vs LCDM), BELOW these errors -> parameters are essentially
#                unconstrained. This is an honest null / detectability statement.
#   "forecast" — errors tightened to a futuristic level at which the signal
#                crosses the detection threshold, used to validate that the
#                inference pipeline recovers injected parameters.
PRECISION_PRESETS = {
    #            (frac H, frac D_M, abs fsigma8)
    "current":  (0.030, 0.020, 0.0250),
    "forecast": (0.0030, 0.0020, 0.0030),
}


def make_synthetic_data(p_true: Params, seed: int = 42, precision: str = "forecast"):
    if precision not in PRECISION_PRESETS:
        raise ValueError(f"precision must be one of {list(PRECISION_PRESETS)}")
    f_H, f_DM, a_fs8 = PRECISION_PRESETS[precision]
    rng = np.random.default_rng(seed)

    z_H = np.array([0.10, 0.20, 0.35, 0.50, 0.68, 0.88, 1.04, 1.30, 1.65, 1.96])
    z_DM = np.array([0.15, 0.38, 0.51, 0.61, 0.70, 0.85, 1.00, 1.48])
    z_fs8 = np.array([0.15, 0.25, 0.38, 0.51, 0.61, 0.78, 0.93, 1.40])

    truth = {
        "H": observables(z_H, p_true)["H"],
        "D_M": observables(z_DM, p_true)["D_M"],
        "fs8": observables(z_fs8, p_true)["fs8"],
    }

    sig_H = f_H * truth["H"]
    sig_DM = np.maximum(f_DM * truth["D_M"], 0.5)
    sig_fs8 = np.full_like(truth["fs8"], a_fs8)

    probes = [
        ProbeData(z_H, truth["H"] + rng.normal(0, sig_H), sig_H, "H"),
        ProbeData(z_DM, truth["D_M"] + rng.normal(0, sig_DM), sig_DM, "D_M"),
        ProbeData(z_fs8, truth["fs8"] + rng.normal(0, sig_fs8), sig_fs8, "fs8"),
    ]
    return probes
