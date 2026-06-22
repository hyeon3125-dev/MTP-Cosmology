"""
screening.py — chameleon thin-shell screening estimate (P3).

Addresses report §4-2: the cosmological coupling beta(z) must be suppressed in
high-density environments (the Solar System) to satisfy the Cassini bound on
the PPN parameter, |gamma - 1| <= 2.3e-5.

We adopt a chameleon scalar field with a density-dependent effective mass. For a
compact body of radius R and Newtonian potential Phi_N, the field develops a
"thin shell": only a thin outer layer of thickness Delta R sources the external
fifth force. The thin-shell factor is (Khoury & Weltman 2004)

    Delta R / R  ≈  (phi_out - phi_in) / (6 beta_c M_Pl Phi_N)

and the effective external coupling is suppressed by this factor:

    beta_eff  =  beta_c * (Delta R / R)            (when Delta R/R << 1)

The measurable PPN deviation from a thin-shelled Sun is

    |gamma - 1|  ≈  2 beta_eff^2 / (1 + beta_eff^2)  ≈  2 beta_eff^2.

This module gives an order-of-magnitude screening estimate (NOT a full
scalar-tensor solve): enough to show the cosmological coupling is compatible
with Solar-System tests, while a complete micro-model remains future work.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

# Constants (SI)
G = 6.674e-11          # m^3 kg^-1 s^-2
C = 2.998e8            # m/s
M_SUN = 1.989e30       # kg
R_SUN = 6.957e8        # m

CASSINI_BOUND = 2.3e-5  # |gamma - 1| upper limit (Bertotti et al. 2003)


def newtonian_potential(mass: float = M_SUN, radius: float = R_SUN) -> float:
    """Dimensionless surface Newtonian potential Phi_N = G M / (R c^2)."""
    return G * mass / (radius * C ** 2)


@dataclass
class ThinShellResult:
    phi_N: float
    thin_shell_factor: float
    beta_eff: float
    gamma_minus_1: float
    passes_cassini: bool


def thin_shell(beta_c: float, delta_phi_over_Mpl: float,
               mass: float = M_SUN, radius: float = R_SUN) -> ThinShellResult:
    """
    Thin-shell screening for a body of given mass/radius.

    Parameters
    ----------
    beta_c : chameleon matter coupling (order unity in unscreened theories).
    delta_phi_over_Mpl : (phi_out - phi_in)/M_Pl, the field excursion between the
        cosmological and high-density minima, in Planck units. For a chameleon
        screened deep inside a dense body this is tiny.
    """
    phi_N = newtonian_potential(mass, radius)
    ts = delta_phi_over_Mpl / (6.0 * beta_c * phi_N)
    ts = min(ts, 1.0)  # unscreened bodies saturate at 1 (no thin shell)
    beta_eff = beta_c * ts
    gamma_m1 = 2.0 * beta_eff ** 2 / (1.0 + beta_eff ** 2)
    return ThinShellResult(
        phi_N=phi_N,
        thin_shell_factor=ts,
        beta_eff=beta_eff,
        gamma_minus_1=gamma_m1,
        passes_cassini=gamma_m1 <= CASSINI_BOUND,
    )


def max_delta_phi_for_cassini(beta_c: float,
                              mass: float = M_SUN, radius: float = R_SUN) -> float:
    """
    Largest field excursion (phi_out-phi_in)/M_Pl that still satisfies Cassini.

    Inverts |gamma-1| <= 2.3e-5 -> beta_eff <= sqrt(bound/2) -> thin-shell
    factor <= that / beta_c -> delta_phi <= 6 beta_c Phi_N * (factor).
    """
    phi_N = newtonian_potential(mass, radius)
    beta_eff_max = np.sqrt(CASSINI_BOUND / 2.0)
    ts_max = beta_eff_max / beta_c
    return 6.0 * beta_c * phi_N * ts_max


if __name__ == "__main__":
    print("Chameleon thin-shell screening — Solar-System (Cassini) check")
    print(f"  Phi_N(Sun) = {newtonian_potential():.3e}")
    print(f"  Cassini bound: |gamma-1| <= {CASSINI_BOUND:.1e}\n")

    # Example: order-unity coupling, small field excursion (screened chameleon).
    for dphi in (1e-6, 1e-7, 1e-8):
        r = thin_shell(beta_c=1.0, delta_phi_over_Mpl=dphi)
        flag = "PASS" if r.passes_cassini else "FAIL"
        print(f"  delta_phi/M_Pl={dphi:.0e}: thin-shell={r.thin_shell_factor:.3e}, "
              f"beta_eff={r.beta_eff:.3e}, |gamma-1|={r.gamma_minus_1:.3e}  [{flag}]")

    dmax = max_delta_phi_for_cassini(beta_c=1.0)
    print(f"\n  Max delta_phi/M_Pl satisfying Cassini (beta_c=1): {dmax:.3e}")
    print("  => any chameleon with field excursion below this is Solar-System safe.")
