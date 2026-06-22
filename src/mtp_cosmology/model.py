"""
MTP Cosmology — Windowed Interacting Dark Energy (IDE) toy model.

Core, side-effect-free model module. Importable; no I/O, no plotting.

Background equations (see ARCHITECTURE.md §2 for derivation):

    Continuity (energy exchange Q between CDM and DE):
        rho_c'  * (1+z) = 3 rho_c          - Q/H
        rho_de' * (1+z) = 3(1+w) rho_de     + Q/H            (w = -1 here)

    Coupling kernel:
        Q/H = beta(z) * F_hier(z) * W_late(z) * rho_c
        F_hier(z) = [ log(l(z)/L_P) / log(l_0/L_P) ]^n      (hierarchy term)
        W_late(z) = exp[ -(z - z_star)^2 / (2 sigma^2) ]     (late-time window)
        beta(z)   = beta0 * tanh((z - z_c)/dz)               (sign/strength)

    Friedmann:
        H^2(z) = H0^2 [ Om0 rho_c/rho_c0 + Ode0 rho_de/rho_de0 + Or0 (1+z)^4 ]

Free parameters: beta0, z_star, sigma, z_c, dz  (n has negligible effect;
F_hier is ~constant over z<~2, which is the model's honest self-assessment).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

# ── Fixed cosmology (Planck 2018 TT,TE,EE+lowE+lensing base-LCDM) ────────────
H0 = 67.4          # km/s/Mpc
OM0 = 0.315        # matter (CDM here treated as the coupled component)
ODE0 = 0.685       # dark energy
OR0 = 9.2e-5       # radiation
C_KM = 299792.458  # km/s

# Planck length and a representative coarse-graining scale for the hierarchy term.
L_P = 1.616e-35    # m
L_0 = 1.0e26       # m  (~ Hubble-ish coarse scale)
LOG_RATIO = np.log(L_0 / L_P)


@dataclass(frozen=True)
class Params:
    """Free + nuisance parameters of the windowed-IDE kernel."""
    beta0: float
    z_star: float
    sigma: float
    z_c: float = 0.5
    dz: float = 0.2
    n: float = 2.0

    def as_array(self) -> np.ndarray:
        return np.array([self.beta0, self.z_star, self.sigma, self.z_c, self.dz, self.n])

    @staticmethod
    def from_iterable(it) -> "Params":
        vals = list(it)
        # Accept 3 (beta0,z_star,sigma), 5 (+z_c,dz) or 6 (+n) element forms.
        defaults = [None, None, None, 0.5, 0.2, 2.0]
        for i, v in enumerate(vals):
            defaults[i] = v
        return Params(*defaults)  # type: ignore[arg-type]


# ── Kernel components ────────────────────────────────────────────────────────
def F_hier(z, n: float = 2.0):
    """Logarithmic hierarchy term. l(z) = L_0/(1+z). Nearly constant for z<~2."""
    l_z = L_0 / (1.0 + z)
    return (np.log(l_z / L_P) / LOG_RATIO) ** n


def W_late(z, z_star: float, sigma: float):
    """Gaussian late-time activation window. Auto-suppresses CMB/BBN epochs."""
    return np.exp(-((z - z_star) ** 2) / (2.0 * sigma ** 2))


def beta_z(z, beta0: float, z_c: float, dz: float):
    """Sign-switching coupling amplitude."""
    return beta0 * np.tanh((z - z_c) / dz)


def Q_over_H_per_rhoc(z, p: Params):
    """(Q/H)/rho_c  — the dimensionless coupling kernel."""
    return beta_z(z, p.beta0, p.z_c, p.dz) * F_hier(z, p.n) * W_late(z, p.z_star, p.sigma)


# ── Background ODE integration ───────────────────────────────────────────────
def _ide_system(z, y, p: Params):
    rho_c, rho_de = y
    qh = Q_over_H_per_rhoc(z, p) * rho_c
    w = -1.0
    drho_c = (3.0 * rho_c - qh) / (1.0 + z)
    drho_de = (3.0 * (1.0 + w) * rho_de + qh) / (1.0 + z)
    return [drho_c, drho_de]


def solve_densities(z_arr, p: Params):
    """Integrate normalized densities rho_c(z), rho_de(z) with rho(0)=1."""
    sol = solve_ivp(
        _ide_system, (0.0, float(np.max(z_arr))), [1.0, 1.0],
        t_eval=z_arr, args=(p,), method="RK45", rtol=1e-9, atol=1e-11,
    )
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")
    return sol.y[0], sol.y[1]


def Hz(z_arr, p: Params):
    """Hubble rate H(z) [km/s/Mpc] plus the density solutions."""
    rho_c, rho_de = solve_densities(z_arr, p)
    e2 = OM0 * rho_c + ODE0 * rho_de + OR0 * (1.0 + z_arr) ** 4
    return H0 * np.sqrt(np.maximum(e2, 1e-12)), rho_c, rho_de


def Hz_LCDM(z):
    return H0 * np.sqrt(OM0 * (1.0 + z) ** 3 + ODE0 + OR0 * (1.0 + z) ** 4)


def w_eff(z_arr, rho_de):
    """Effective DE equation of state from d ln rho_de / d ln(1+z)."""
    ln_r = np.log(np.maximum(rho_de, 1e-12))
    return -1.0 + np.gradient(ln_r, np.log(1.0 + z_arr)) / 3.0


def D_M(z_arr, Hz_vals):
    """Comoving (transverse) distance via trapezoidal integration of c/H."""
    integ = C_KM / Hz_vals
    d = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(z_arr))])
    return d


def growth_fs8(z_arr, Hz_vals, rho_c, sigma8_0: float = 0.811):
    """Approximate f*sigma8(z) using the f ~ Omega_m(z)^0.55 growth-index ansatz."""
    a = 1.0 / (1.0 + z_arr)
    om_z = OM0 * rho_c / (Hz_vals / H0) ** 2
    f = om_z ** 0.55
    ln_d = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] / a[1:] + f[:-1] / a[:-1]) * np.diff(a))])
    d = np.exp(ln_d)
    d = d / (d[0] if d[0] != 0 else 1.0)
    sigma8_z = sigma8_0 * d / d[0]
    return f * sigma8_z


def cpl_fit(z_arr, w_arr, z_max_fit: float = 1.5):
    """Least-squares CPL fit  w(a) = w0 + wa (1-a)  over z < z_max_fit."""
    a = 1.0 / (1.0 + z_arr)
    design = np.vstack([np.ones_like(a), 1.0 - a]).T
    mask = z_arr < z_max_fit
    w0, wa = np.linalg.lstsq(design[mask], w_arr[mask], rcond=None)[0]
    return float(w0), float(wa)
