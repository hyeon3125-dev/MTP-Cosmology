"""
mtp_v01_prototype.py — v0.1 prototype ansatz (RECONSTRUCTED, kept for history).

Documents the three structural defects that motivated the IDE rewrite. This
file is *meant* to fail the phantom-crossing test; do not use for results.

Original ansatz:
    Q(z) = alpha * H * rho_c * f(z),   f(z) = [log(l(z)/L_P)/log(l_0/L_P)]^n

Defects (see MTP_consolidated_report.md §1, ARCHITECTURE.md §1):
  1. Q passed directly into the DE ODE -> violates nabla_mu T^{mu nu}=0 (GR).
  2. w_eff = -1 + (positive term) -> phantom barrier, w<-1 impossible.
  3. alpha interpreted as G_eff/G_N - 1 -> exceeds Cassini bound ~100x.
  Plus: f(z->0)->0, so the transfer term vanishes today -> w pinned to -1.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

H0, OM0, ODE0 = 67.4, 0.315, 0.685
L_P, L_0 = 1.616e-35, 1.0e26
LOG_RATIO = np.log(L_0 / L_P)


def f_z(z, n=2.0):
    l_z = L_0 / (1.0 + z)
    return (np.log(l_z / L_P) / LOG_RATIO) ** n


def w_eff_v01(z, alpha=0.1, n=2.0):
    """DEFECT #2 demonstration: the correction is strictly >= 0 -> w >= -1."""
    return -1.0 + alpha * f_z(z, n)  # positive term => phantom barrier


if __name__ == "__main__":
    z = np.linspace(0.0, 2.5, 50)
    w = w_eff_v01(z)
    print("v0.1 prototype — phantom barrier demonstration")
    print(f"  min w_eff = {w.min():.4f}  (>= -1 always -> NO phantom crossing)")
    print("  -> superseded by v0.2 (IDE) and v0.3 (windowed IDE).")
