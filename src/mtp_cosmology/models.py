"""
models.py — unified model registry for fair comparison (see
docs/comparison_methodology.yaml).

Every model exposes the SAME interface so the likelihood, sampler, and metrics
are identical across models (fairness_rules):

    model.name
    model.param_names      list of free dark-sector parameter names
    model.bounds           {name: (lo, hi)} flat-prior bounds
    model.H(z, theta)      H(z) [km/s/Mpc] given the free params

Background-level comparison: the early-universe cosmology is fixed to Planck
2018 (Om0, Ode0, Or0, H0), so only dark-sector parameters are free. This keeps
the comparison tractable without a Boltzmann code; freeing H0/Om is a later
stage (phase_2+). LambdaCDM therefore has zero free parameters here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from .model import H0, OM0, ODE0, OR0, Hz_LCDM


# ── shared IDE background solver (arbitrary coupling kernel) ──────────────────
def solve_ide(z_arr, kernel: Callable[[float], float]):
    """
    Solve the IDE continuity system with DE w=-1 and coupling kernel
    k(z) = (Q/H)/rho_c.  Returns normalized rho_c(z), rho_de(z) (rho(0)=1).

    kernel = 0 reproduces LambdaCDM: rho_c=(1+z)^3, rho_de=1.
    """
    def sys(z, y):
        rho_c, rho_de = y
        qh = kernel(z) * rho_c
        drho_c = (3.0 * rho_c - qh) / (1.0 + z)
        drho_de = qh / (1.0 + z)          # 3(1+w)=0 for w=-1
        return [drho_c, drho_de]

    sol = solve_ivp(sys, (0.0, float(np.max(z_arr))), [1.0, 1.0],
                    t_eval=z_arr, method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0], sol.y[1]


def _H_from_densities(z, rho_c, rho_de):
    e2 = OM0 * rho_c + ODE0 * rho_de + OR0 * (1.0 + z) ** 4
    return H0 * np.sqrt(np.maximum(e2, 1e-12))


@dataclass
class Model:
    name: str
    param_names: list
    bounds: dict
    _H: Callable
    description: str = ""

    @property
    def k(self) -> int:
        return len(self.param_names)

    def H(self, z, theta):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        return self._H(z, np.atleast_1d(np.asarray(theta, dtype=float)))


# ── LambdaCDM (0 free params) ────────────────────────────────────────────────
def _H_lcdm(z, theta):
    return Hz_LCDM(z)


# ── CPL w0-wa ────────────────────────────────────────────────────────────────
def _H_cpl(z, theta):
    w0, wa = theta
    # rho_de(z)/rho_de0 = (1+z)^{3(1+w0+wa)} exp(-3 wa z/(1+z))
    f_de = (1.0 + z) ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * z / (1.0 + z))
    e2 = OM0 * (1.0 + z) ** 3 + ODE0 * f_de + OR0 * (1.0 + z) ** 4
    return H0 * np.sqrt(np.maximum(e2, 1e-12))


# ── constant-coupling IDE: Q = xi H rho_c ────────────────────────────────────
def _H_ide_const(z, theta):
    (xi,) = theta
    rho_c, rho_de = solve_ide(z, lambda zz: xi)
    return _H_from_densities(z, rho_c, rho_de)


# ── sign-switching IDE: xi(z) = xi0 tanh((z-zc)/dz) ──────────────────────────
def _H_ide_signswitch(z, theta):
    xi0, zc, dz = theta
    rho_c, rho_de = solve_ide(z, lambda zz: xi0 * np.tanh((zz - zc) / dz))
    return _H_from_densities(z, rho_c, rho_de)


# ── MTP windowed IDE, 3-param variant (zc=z*, dz=sigma fixed) ────────────────
def _mtp_kernel(beta0, z_star, sigma):
    from .model import F_hier
    def k(zz):
        beta = beta0 * np.tanh((zz - z_star) / sigma)     # zc=z*, dz=sigma
        w = np.exp(-((zz - z_star) ** 2) / (2.0 * sigma ** 2))
        return beta * F_hier(zz, 2.0) * w
    return k


def _H_mtp3p(z, theta):
    beta0, z_star, sigma = theta
    rho_c, rho_de = solve_ide(z, _mtp_kernel(beta0, z_star, sigma))
    return _H_from_densities(z, rho_c, rho_de)


# ── MTP windowed IDE, 4-param variant (zc free; dz=sigma) ────────────────────
def _H_mtp4p(z, theta):
    beta0, z_star, sigma, zc = theta
    from .model import F_hier
    def k(zz):
        beta = beta0 * np.tanh((zz - zc) / sigma)
        w = np.exp(-((zz - z_star) ** 2) / (2.0 * sigma ** 2))
        return beta * F_hier(zz, 2.0) * w
    rho_c, rho_de = solve_ide(z, k)
    return _H_from_densities(z, rho_c, rho_de)


# Flat-prior bounds follow docs/comparison_methodology.yaml.
REGISTRY = {
    "lcdm": Model("lcdm", [], {}, _H_lcdm, "base control (0 free DE params)"),
    "cpl_w0wa": Model("cpl_w0wa", ["w0", "wa"],
                      {"w0": (-2.0, -0.3), "wa": (-3.0, 3.0)},
                      _H_cpl, "CPL dynamical DE"),
    "standard_ide": Model("standard_ide", ["xi"],
                          {"xi": (-0.3, 0.3)},
                          _H_ide_const, "constant-coupling IDE"),
    "sign_switching_ide": Model("sign_switching_ide", ["xi0", "zc", "dz"],
                                {"xi0": (-0.3, 0.3), "zc": (0.0, 1.5), "dz": (0.02, 1.0)},
                                _H_ide_signswitch, "sign-switching IDE"),
    "mtp_3p": Model("mtp_3p", ["beta0", "z_star", "sigma"],
                    {"beta0": (0.0, 0.3), "z_star": (0.0, 1.5), "sigma": (0.05, 1.0)},
                    _H_mtp3p, "MTP windowed IDE (zc=z*, dz=sigma)"),
    "mtp_4p": Model("mtp_4p", ["beta0", "z_star", "sigma", "zc"],
                    {"beta0": (0.0, 0.3), "z_star": (0.0, 1.5),
                     "sigma": (0.05, 1.0), "zc": (0.0, 1.5)},
                    _H_mtp4p, "MTP windowed IDE (zc free)"),
}


def get(name: str) -> Model:
    if name not in REGISTRY:
        raise KeyError(f"unknown model '{name}'; have {list(REGISTRY)}")
    return REGISTRY[name]
