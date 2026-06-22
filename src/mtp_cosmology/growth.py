"""
growth.py — linear growth and f*sigma8(z) for the comparison set (phase_3).

Solves the standard sub-horizon growth ODE in e-folds N = ln a,

    δ'' + (2 + d ln H/dN) δ' - (3/2) Ω_m(a) δ = 0,

using each model's actual H(a) and matter density Ω_m(a). For the IDE models
Ω_m(a) carries the coupling's effect on the background matter density (ρ_c is
modified by the energy exchange), so the growth differs from ΛCDM/CPL beyond the
pure expansion history.

APPROXIMATION (documented honestly): this is *background-modified* growth — it
omits the explicit dark-sector momentum-transfer perturbation and any scale
dependence (quasi-static, single-fluid δ_c). A full perturbed-IDE Boltzmann
treatment is phase_4. σ8(0) is fixed to the Planck value for all models so the
comparison is on shape, not normalization.
"""
from __future__ import annotations

import numpy as np

from .model import OM0, H0


def fsigma8(model, theta, z_eval, sigma8_0: float = 0.811, n_grid: int = 300):
    """Return f*sigma8 at the requested redshifts for `model` at `theta`."""
    z_eval = np.atleast_1d(np.asarray(z_eval, dtype=float))
    z_max = max(float(z_eval.max()) + 0.5, 3.0)

    # work on an a-grid from deep matter domination to today
    a = np.linspace(1.0 / (1.0 + z_max), 1.0, n_grid)
    z = 1.0 / a - 1.0
    E = model.H(z, theta) / H0
    rho_m = model.rho_m(z, theta)                 # normalized rho_c/rho_c0
    Om_a = OM0 * rho_m / E ** 2                    # fractional matter density

    N = np.log(a)
    dlnH_dN = np.gradient(np.log(E), N)

    # integrate δ(N), δ'(N) with matter-dominated init: δ ∝ a => δ=a, δ'=a
    delta = np.zeros(n_grid)
    dprime = np.zeros(n_grid)
    delta[0] = a[0]
    dprime[0] = a[0]
    for i in range(n_grid - 1):
        h = N[i + 1] - N[i]
        # RHS of the 2nd-order system at i
        d2 = -(2.0 + dlnH_dN[i]) * dprime[i] + 1.5 * Om_a[i] * delta[i]
        # simple RK2 (midpoint) for stability
        dmid = delta[i] + 0.5 * h * dprime[i]
        dpmid = dprime[i] + 0.5 * h * d2
        d2mid = -(2.0 + dlnH_dN[i]) * dpmid + 1.5 * Om_a[i] * dmid
        delta[i + 1] = delta[i] + h * dpmid
        dprime[i + 1] = dprime[i] + h * d2mid

    f = dprime / delta                              # dlnδ/dlna
    D = delta / delta[-1]                            # normalized to today
    fs8_grid = sigma8_0 * f * D                      # f * sigma8(z), sigma8(z)=sigma8_0 D

    # interpolate (grid is increasing in a => decreasing in z; sort by z)
    order = np.argsort(z)
    return np.interp(z_eval, z[order], fs8_grid[order])
