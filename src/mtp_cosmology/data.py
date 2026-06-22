"""
data.py — real cosmological data tables (P2).

DESI DR1 BAO (DESI 2024 VI, arXiv:2404.03002, Table 1). Anisotropic tracers
report D_M/r_d and D_H/r_d with a correlation coefficient r; BGS and QSO report
the isotropic D_V/r_d only.

D_H(z) = c / H(z);  D_V(z) = [ z * D_M(z)^2 * D_H(z) ]^(1/3).

The windowed-IDE coupling is suppressed at high z by W_late, so the early-time
sound horizon is unchanged from LambdaCDM: we fix r_d = 147.09 Mpc (Planck 2018).
This is the standard assumption for late-time-only dark-energy fits.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

# Planck 2018 sound horizon at the drag epoch (Mpc). Fixed for late-time fits.
R_D = 147.09

# DESI DR1 BAO, Table 1 of arXiv:2404.03002.
# Anisotropic rows: (z_eff, DM/rd, sigma_DM, DH/rd, sigma_DH, corr_DM_DH)
DESI_ANISO = [
    ("LRG1",      0.510, 13.62, 0.25, 20.98, 0.61, -0.445),
    ("LRG2",      0.706, 16.85, 0.32, 20.08, 0.60, -0.420),
    ("LRG3+ELG1", 0.930, 21.71, 0.28, 17.88, 0.35, -0.389),
    ("ELG2",      1.317, 27.79, 0.69, 13.82, 0.42, -0.444),
    ("Lya",       2.330, 39.71, 0.94,  8.52, 0.17, -0.477),
]
# Isotropic rows: (z_eff, DV/rd, sigma_DV)
DESI_ISO = [
    ("BGS", 0.295,  7.93, 0.15),
    ("QSO", 1.491, 26.07, 0.67),
]


@dataclass
class BAOAniso:
    name: str
    z: float
    DM_rd: float
    DH_rd: float
    cov: np.ndarray  # 2x2 covariance of (DM/rd, DH/rd)


@dataclass
class BAOIso:
    name: str
    z: float
    DV_rd: float
    sigma: float


def desi_dr1():
    """Return (aniso_list, iso_list) of DESI DR1 BAO data objects."""
    aniso = []
    for name, z, dm, sdm, dh, sdh, r in DESI_ANISO:
        cov = np.array([[sdm ** 2, r * sdm * sdh],
                        [r * sdm * sdh, sdh ** 2]])
        aniso.append(BAOAniso(name, z, dm, dh, cov))
    iso = [BAOIso(name, z, dv, s) for name, z, dv, s in DESI_ISO]
    return aniso, iso
