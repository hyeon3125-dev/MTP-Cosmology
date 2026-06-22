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


# "Gold-2018" growth-rate compilation (Sagredo, Nesseris & Perivolaropoulos
# 2018, arXiv:1806.10822, Table III): (z, fsigma8, sigma). The per-point
# fiducial Omega_m (used for the Alcock-Paczynski correction) is omitted here;
# we fit fsigma8 directly (a standard first-order approximation).
RSD_GOLD2018 = [
    (0.02, 0.428, 0.0465), (0.02, 0.398, 0.065), (0.02, 0.314, 0.048),
    (0.10, 0.370, 0.130), (0.15, 0.490, 0.145), (0.17, 0.510, 0.060),
    (0.18, 0.360, 0.090), (0.38, 0.440, 0.060), (0.25, 0.3512, 0.0583),
    (0.37, 0.4602, 0.0378), (0.32, 0.384, 0.095), (0.59, 0.488, 0.060),
    (0.44, 0.413, 0.080), (0.60, 0.390, 0.063), (0.73, 0.437, 0.072),
    (0.60, 0.550, 0.120), (0.86, 0.400, 0.110), (1.40, 0.482, 0.116),
    (0.978, 0.379, 0.176), (1.23, 0.385, 0.099), (1.526, 0.342, 0.070),
    (1.944, 0.364, 0.106),
]


def rsd_gold2018():
    """Return (z, fsigma8, sigma) arrays for the Gold-2018 RSD compilation."""
    arr = np.array(RSD_GOLD2018)
    return arr[:, 0], arr[:, 1], arr[:, 2]


# Planck 2018 compressed-CMB distance priors (Chen, Huang & Wang 2019,
# arXiv:1808.05724, Table I; LambdaCDM, TT,TE,EE+lowE). For a late-time-only DE
# model these carry the full geometric CMB information (R and l_A depend on the
# expansion history out to decoupling). omega_b and n_s are early-time and held
# fixed, so we use the (R, l_A) 2-parameter prior our background predicts.
PLANCK18_R = 1.7502
PLANCK18_LA = 301.471
PLANCK18_SIG_R = 0.0046
PLANCK18_SIG_LA = 0.090
PLANCK18_CORR_R_LA = 0.46
PLANCK18_ZSTAR = 1089.92          # decoupling redshift (Planck 2018)
# Sound horizon at z_star [Mpc]. Planck 2018 gives ~144.43; we use 144.476,
# calibrated so the fixed fiducial Planck-LambdaCDM background reproduces the
# prior's central l_A exactly (absorbs the precise r_s/z* and the numerical
# D_M(z*) integration). Well within the ~0.3 Mpc r_s uncertainty; affects all
# models equally and so cancels in the differential comparison.
PLANCK18_RS_ZSTAR = 144.476


def planck18_distance_prior():
    """Return (data=[R, l_A], cov 2x2, z_star, r_s) for the Planck 2018 prior."""
    c = PLANCK18_CORR_R_LA * PLANCK18_SIG_R * PLANCK18_SIG_LA
    cov = np.array([[PLANCK18_SIG_R ** 2, c],
                    [c, PLANCK18_SIG_LA ** 2]])
    data = np.array([PLANCK18_R, PLANCK18_LA])
    return data, cov, PLANCK18_ZSTAR, PLANCK18_RS_ZSTAR


def desi_dr1():
    """Return (aniso_list, iso_list) of DESI DR1 BAO data objects."""
    aniso = []
    for name, z, dm, sdm, dh, sdh, r in DESI_ANISO:
        cov = np.array([[sdm ** 2, r * sdm * sdh],
                        [r * sdm * sdh, sdh ** 2]])
        aniso.append(BAOAniso(name, z, dm, dh, cov))
    iso = [BAOIso(name, z, dv, s) for name, z, dv, s in DESI_ISO]
    return aniso, iso
