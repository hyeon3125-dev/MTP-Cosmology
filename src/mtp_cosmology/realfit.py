"""
realfit.py — fit the windowed-IDE model to real DESI DR1 BAO data and compare
its Bayesian evidence against LambdaCDM (P2).

LambdaCDM here is the beta0 -> 0 limit of the model with the same fixed Planck
2018 background, i.e. it has ZERO free parameters. The windowed-IDE model adds
the coupling parameters. Bayesian evidence therefore applies an automatic Occam
penalty for the extra freedom, so a near-zero or negative Delta ln Z is the
expected honest outcome if sub-percent late-time coupling is not required by the
data.

Observables (fixed r_d = 147.09 Mpc, justified in data.py):
    D_M/r_d,  D_H/r_d = c/(H r_d),  D_V/r_d = [z D_M^2 D_H]^(1/3) / r_d
"""
from __future__ import annotations

import numpy as np

from .model import Params, Hz, Hz_LCDM, D_M, C_KM
from .data import desi_dr1, R_D
from .likelihood import PRIOR_BOUNDS

ALL_NAMES = ["beta0", "z_star", "sigma", "z_c", "dz"]
_TRUTH_FIXED = {"z_star": 0.40, "sigma": 0.30, "z_c": 0.50, "dz": 0.20}


def _grid_observables(p: Params, z_max: float, n_grid: int = 300):
    grid = np.linspace(1e-4, z_max * 1.001, n_grid)
    H, _, _ = Hz(grid, p)
    DM = D_M(grid, H)
    return grid, H, DM


def _grid_observables_lcdm(z_max: float, n_grid: int = 300):
    grid = np.linspace(1e-4, z_max * 1.001, n_grid)
    H = Hz_LCDM(grid)
    DM = D_M(grid, H)
    return grid, H, DM


def _bao_at(z, grid, H, DM):
    dm = np.interp(z, grid, DM)
    h = np.interp(z, grid, H)
    dh = C_KM / h
    dv = (z * dm ** 2 * dh) ** (1.0 / 3.0)
    return dm / R_D, dh / R_D, dv / R_D


def chi2_desi(grid, H, DM) -> float:
    """chi^2 of the DESI DR1 BAO data given model grids of H(z), D_M(z)."""
    aniso, iso = desi_dr1()
    chi2 = 0.0
    for b in aniso:
        dm_rd, dh_rd, _ = _bao_at(b.z, grid, H, DM)
        resid = np.array([dm_rd - b.DM_rd, dh_rd - b.DH_rd])
        chi2 += float(resid @ np.linalg.solve(b.cov, resid))
    for b in iso:
        _, _, dv_rd = _bao_at(b.z, grid, H, DM)
        chi2 += ((dv_rd - b.DV_rd) / b.sigma) ** 2
    return chi2


def make_params(theta, fit_names) -> Params:
    """Build a Params from a sampled subset; fixed params from _TRUTH_FIXED."""
    vals = dict(_TRUTH_FIXED)
    vals["beta0"] = 0.0
    for name, v in zip(fit_names, np.atleast_1d(theta)):
        vals[name] = float(v)
    return Params(beta0=vals["beta0"], z_star=vals["z_star"], sigma=vals["sigma"],
                  z_c=vals["z_c"], dz=vals["dz"])


def loglike_model(theta, fit_names, z_max: float = 2.4) -> float:
    p = make_params(theta, fit_names)
    try:
        grid, H, DM = _grid_observables(p, z_max)
        chi2 = chi2_desi(grid, H, DM)
    except (RuntimeError, FloatingPointError, ValueError):
        return -1e30
    if not np.isfinite(chi2):
        return -1e30
    return -0.5 * chi2


def loglike_lcdm(z_max: float = 2.4) -> float:
    grid, H, DM = _grid_observables_lcdm(z_max)
    return -0.5 * chi2_desi(grid, H, DM)


def prior_transform(u, fit_names):
    """Map unit-cube samples to flat priors over PRIOR_BOUNDS (dynesty)."""
    x = np.empty_like(u)
    for i, name in enumerate(fit_names):
        lo, hi = PRIOR_BOUNDS[name]
        x[i] = lo + (hi - lo) * u[i]
    return x
