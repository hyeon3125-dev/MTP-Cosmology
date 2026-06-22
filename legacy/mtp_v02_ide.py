"""
mtp_v02_ide.py — v0.2 IDE transition (RECONSTRUCTED, kept for history).

Fixes v0.1 defects 1-3 but the effect is too small (1/10-1/30 of the DESI
requirement) because the log structure f(z) dies off at low z. v0.3 cures this
by factoring f(z) -> F_hier(z) * W_late(z) (late-time activation window).

Changes vs v0.1:
  - ODE -> IDE 4-vector coupling (satisfies nabla_mu T^{mu nu}=0).   [fix #1]
  - alpha(z) = alpha0 * tanh(z - z_c)  sign switching.               [enables w<-1]
  - alpha redefined as IDE coupling, not G_eff/G_N-1 (avoids Cassini). [fix #3]

Residual problem: amplitude shortfall at observable redshifts. See report §1.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

H0, OM0, ODE0, OR0 = 67.4, 0.315, 0.685, 9.2e-5
L_P, L_0 = 1.616e-35, 1.0e26
LOG_RATIO = np.log(L_0 / L_P)


def f_z(z, n=2.0):
    l_z = L_0 / (1.0 + z)
    return (np.log(l_z / L_P) / LOG_RATIO) ** n


def alpha_z(z, alpha0, z_c, dz):
    return alpha0 * np.tanh((z - z_c) / dz)


def _system(z, y, alpha0, z_c, dz, n):
    rho_c, rho_de = y
    qh = alpha_z(z, alpha0, z_c, dz) * f_z(z, n) * rho_c
    w = -1.0
    drho_c = (3.0 * rho_c - qh) / (1.0 + z)
    drho_de = (3.0 * (1.0 + w) * rho_de + qh) / (1.0 + z)
    return [drho_c, drho_de]


def w_eff_v02(z_arr, alpha0=0.1, z_c=0.5, dz=0.2, n=2.0):
    sol = solve_ivp(_system, (0, z_arr.max()), [1.0, 1.0], t_eval=z_arr,
                    args=(alpha0, z_c, dz, n), method="RK45", rtol=1e-9, atol=1e-11)
    rde = sol.y[1]
    ln_r = np.log(np.maximum(rde, 1e-12))
    return -1.0 + np.gradient(ln_r, np.log(1.0 + z_arr)) / 3.0


if __name__ == "__main__":
    z = np.linspace(0.001, 2.5, 500)
    w = w_eff_v02(z)
    print("v0.2 IDE — GR-consistent but amplitude too small")
    print(f"  min w_eff = {w.min():.5f}  (crosses but |w+1| << DESI requirement)")
    print("  -> v0.3 adds W_late(z) window to boost low-z amplitude.")
