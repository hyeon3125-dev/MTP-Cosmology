"""
compare.py — fair model-comparison engine (see docs/comparison_methodology.yaml).

A single likelihood implementation, optimizer, and metric set applied identically
to every model in models.REGISTRY (fairness_rules). Supports three observable
blocks built from a model's H(z):

    H block   : cosmic-chronometer-style H(z) measurements
    BAO block : DESI-style D_M/r_d, D_H/r_d (2x2 cov) and D_V/r_d
    SN block  : distance modulus mu(z), with the absolute offset M analytically
                marginalized (flat prior) so it never biases the comparison

Metrics: chi2_min (best fit via multi-start optimization), AIC = chi2+2k,
BIC = chi2 + k ln N, and Bayesian log-evidence via nested sampling (dynesty).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .model import C_KM, OM0, H0, Hz_LCDM
from .data import desi_dr1, rsd_gold2018, planck18_distance_prior, R_D
from . import models as M


# ── observable blocks ────────────────────────────────────────────────────────
@dataclass
class HBlock:
    z: np.ndarray
    value: np.ndarray      # H(z) [km/s/Mpc]
    sigma: np.ndarray

@dataclass
class SNBlock:
    z: np.ndarray
    mu: np.ndarray         # distance modulus
    sigma: np.ndarray

@dataclass
class BAOAnisoBlock:
    z: np.ndarray
    DM_rd: np.ndarray
    DH_rd: np.ndarray
    covs: list             # list of 2x2 covariance matrices

@dataclass
class BAOIsoBlock:
    z: np.ndarray
    DV_rd: np.ndarray
    sigma: np.ndarray

@dataclass
class RSDBlock:
    z: np.ndarray
    fs8: np.ndarray        # f*sigma8(z)
    sigma: np.ndarray

@dataclass
class CMBDistBlock:
    """Planck compressed-CMB distance prior on (R, l_A)."""
    data: np.ndarray       # [R, l_A]
    cov: np.ndarray        # 2x2
    z_star: float
    r_s: float


@dataclass
class Dataset:
    name: str
    H: HBlock | None = None
    sn: SNBlock | None = None
    bao_aniso: BAOAnisoBlock | None = None
    bao_iso: BAOIsoBlock | None = None
    rsd: RSDBlock | None = None
    cmb: CMBDistBlock | None = None
    r_d: float = R_D

    @property
    def n_points(self) -> int:
        n = 0
        if self.H is not None: n += len(self.H.z)
        if self.sn is not None: n += len(self.sn.z)
        if self.bao_aniso is not None: n += 2 * len(self.bao_aniso.z)
        if self.bao_iso is not None: n += len(self.bao_iso.z)
        if self.rsd is not None: n += len(self.rsd.z)
        if self.cmb is not None: n += 2
        return n

    def _all_z(self):
        zs = [np.array([0.0])]
        for blk in (self.H, self.sn, self.bao_aniso, self.bao_iso, self.rsd):
            if blk is not None:
                zs.append(np.atleast_1d(blk.z))
        return np.concatenate(zs)

    def chi2(self, model: M.Model, theta) -> float:
        z_max = max(float(self._all_z().max()), 0.01)
        grid = np.linspace(1e-4, z_max * 1.001, 300)
        try:
            H = model.H(grid, theta)
        except (RuntimeError, FloatingPointError, ValueError):
            return np.inf
        if not np.all(np.isfinite(H)) or np.any(H <= 0):
            return np.inf
        # comoving distance grid
        integ = C_KM / H
        DM = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(grid))])

        def H_at(z): return np.interp(z, grid, H)
        def DM_at(z): return np.interp(z, grid, DM)

        chi2 = 0.0
        if self.H is not None:
            chi2 += np.sum(((self.H.value - H_at(self.H.z)) / self.H.sigma) ** 2)
        if self.sn is not None:
            dl = (1.0 + self.sn.z) * DM_at(self.sn.z)          # Mpc
            mu_model = 5.0 * np.log10(np.maximum(dl, 1e-6)) + 25.0
            r = self.sn.mu - mu_model
            iv = 1.0 / self.sn.sigma ** 2
            # analytic marginalization over a constant offset (flat prior)
            chi2 += np.sum(r ** 2 * iv) - (np.sum(r * iv)) ** 2 / np.sum(iv)
        if self.bao_aniso is not None:
            b = self.bao_aniso
            dm_rd = DM_at(b.z) / self.r_d
            dh_rd = (C_KM / H_at(b.z)) / self.r_d
            for i in range(len(b.z)):
                resid = np.array([dm_rd[i] - b.DM_rd[i], dh_rd[i] - b.DH_rd[i]])
                chi2 += float(resid @ np.linalg.solve(b.covs[i], resid))
        if self.bao_iso is not None:
            b = self.bao_iso
            dm = DM_at(b.z); dh = C_KM / H_at(b.z)
            dv_rd = (b.z * dm ** 2 * dh) ** (1.0 / 3.0) / self.r_d
            chi2 += np.sum(((dv_rd - b.DV_rd) / b.sigma) ** 2)
        if self.rsd is not None:
            from .growth import fsigma8
            model_fs8 = fsigma8(model, theta, self.rsd.z)
            chi2 += np.sum(((self.rsd.fs8 - model_fs8) / self.rsd.sigma) ** 2)
        if self.cmb is not None:
            R, lA = self._cmb_R_lA(model, theta, self.cmb.z_star, self.cmb.r_s)
            resid = np.array([R - self.cmb.data[0], lA - self.cmb.data[1]])
            chi2 += float(resid @ np.linalg.solve(self.cmb.cov, resid))
        return float(chi2) if np.isfinite(chi2) else np.inf

    @staticmethod
    def _cmb_R_lA(model: M.Model, theta, z_star: float, r_s: float):
        """Shift parameter R and acoustic scale l_A from the model background.

        R   = sqrt(Om0) * (H0/c) * D_M(z*)
        l_A = pi * D_M(z*) / r_s(z*)
        D_M(z*) integrates c/H to decoupling; the model is Planck-LambdaCDM at
        high z (W_late->0) so only the late-time part of the integral shifts.
        """
        # Split the integral at z_lo: below it, evaluate the model H (the IDE
        # window is alive); above it the coupling is dead (W_late->0) so the
        # background is Planck-LambdaCDM and H is analytic — no need to integrate
        # the IDE ODE to z~1090 (that was the bottleneck).
        z_lo = 12.0
        z_below = np.linspace(1e-4, z_lo, 400)
        H_below = model.H(z_below, theta)
        if not np.all(np.isfinite(H_below)) or np.any(H_below <= 0):
            return np.inf, np.inf
        DM_below = np.trapz(C_KM / H_below, z_below)
        # tail z_lo -> z_star (matter+radiation dominated). Analytic models
        # (LCDM/CPL) use their own H; IDE models use LambdaCDM (coupling dead).
        z_above = np.linspace(z_lo, z_star, 4000)
        H_above = model.H(z_above, theta) if model._kernel is None else Hz_LCDM(z_above)
        DM_above = np.trapz(C_KM / H_above, z_above)
        DM_star = DM_below + DM_above                # Mpc, comoving
        R = np.sqrt(OM0) * (H0 / C_KM) * DM_star
        lA = np.pi * DM_star / r_s
        return R, lA


# ── real DESI DR1 BAO dataset ────────────────────────────────────────────────
def desi_dr1_dataset() -> Dataset:
    aniso, iso = desi_dr1()
    ba = BAOAnisoBlock(
        z=np.array([b.z for b in aniso]),
        DM_rd=np.array([b.DM_rd for b in aniso]),
        DH_rd=np.array([b.DH_rd for b in aniso]),
        covs=[b.cov for b in aniso],
    )
    bi = BAOIsoBlock(
        z=np.array([b.z for b in iso]),
        DV_rd=np.array([b.DV_rd for b in iso]),
        sigma=np.array([b.sigma for b in iso]),
    )
    return Dataset("DESI_DR1_BAO", bao_aniso=ba, bao_iso=bi)


def rsd_dataset() -> Dataset:
    """Gold-2018 RSD f*sigma8 compilation as a standalone dataset."""
    z, fs8, sig = rsd_gold2018()
    return Dataset("RSD_Gold2018", rsd=RSDBlock(z, fs8, sig))


def desi_plus_rsd_dataset() -> Dataset:
    """DESI DR1 BAO (geometry) + Gold-2018 RSD (growth) — phase_3."""
    d = desi_dr1_dataset()
    z, fs8, sig = rsd_gold2018()
    d.rsd = RSDBlock(z, fs8, sig)
    d.name = "DESI_DR1_BAO + RSD_Gold2018"
    return d


def cmb_block() -> CMBDistBlock:
    data, cov, z_star, r_s = planck18_distance_prior()
    return CMBDistBlock(data, cov, z_star, r_s)


def full_dataset() -> Dataset:
    """Planck18 compressed CMB + DESI DR1 BAO + Gold-2018 RSD — phase_4."""
    d = desi_plus_rsd_dataset()
    d.cmb = cmb_block()
    d.name = "Planck18(R,lA) + DESI_DR1_BAO + RSD_Gold2018"
    return d


# ── mock dataset from a fiducial model (phase_0 pipeline validation) ─────────
def mock_dataset(fiducial="cpl_w0wa", theta=(-0.85, -0.60), seed=42,
                 err_scale=1.0) -> Dataset:
    """
    Compressed H(z) + BAO + SN mock from a fiducial model.

    Default fiducial is a distinct (non-LambdaCDM) CPL and errors are at a
    forecast level (err_scale multiplies all sigmas) so phase_0 can demonstrate
    that the engine recovers the injected signal and discriminates models.
    """
    rng = np.random.default_rng(seed)
    fid = M.get(fiducial)
    grid = np.linspace(1e-4, 2.4, 400)
    Hg = fid.H(grid, theta)
    integ = C_KM / Hg
    DMg = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(grid))])
    H_at = lambda z: np.interp(z, grid, Hg)
    DM_at = lambda z: np.interp(z, grid, DMg)

    # H(z): 15 chronometer-like points, 1.5% errors
    zH = np.linspace(0.1, 1.9, 15)
    sH = 0.015 * err_scale * H_at(zH)
    H = HBlock(zH, H_at(zH) + rng.normal(0, sH), sH)

    # SN: 40 points to z=1.7, 0.06 mag
    zS = np.linspace(0.05, 1.7, 40)
    mu = 5.0 * np.log10((1.0 + zS) * DM_at(zS)) + 25.0
    sS = np.full_like(zS, 0.06 * err_scale)
    sn = SNBlock(zS, mu + rng.normal(0, sS), sS)

    # BAO at DESI redshifts, 0.8% on DM/rd and DH/rd (diagonal mock cov)
    zB = np.array([0.51, 0.706, 0.93, 1.317, 2.33])
    dm_rd = DM_at(zB) / R_D
    dh_rd = (C_KM / H_at(zB)) / R_D
    sdm = 0.008 * err_scale * dm_rd; sdh = 0.008 * err_scale * dh_rd
    covs = [np.diag([sdm[i] ** 2, sdh[i] ** 2]) for i in range(len(zB))]
    ba = BAOAnisoBlock(zB,
                       dm_rd + rng.normal(0, sdm),
                       dh_rd + rng.normal(0, sdh), covs)
    return Dataset(f"mock[{fiducial}]", H=H, sn=sn, bao_aniso=ba)


# ── fitting + metrics ────────────────────────────────────────────────────────
@dataclass
class Fit:
    name: str
    k: int
    chi2_min: float
    aic: float
    bic: float
    best: dict
    lnZ: float = np.nan
    lnZ_err: float = np.nan


def fit_model(model: M.Model, data: Dataset, n_starts: int = 8, seed: int = 42) -> Fit:
    n = data.n_points
    if model.k == 0:
        chi2 = data.chi2(model, [])
        return Fit(model.name, 0, chi2, chi2, chi2, {})

    rng = np.random.default_rng(seed)
    lo = np.array([model.bounds[p][0] for p in model.param_names])
    hi = np.array([model.bounds[p][1] for p in model.param_names])

    def neg(theta):
        if np.any(theta < lo) or np.any(theta > hi):
            return 1e30
        return data.chi2(model, theta)

    best_val, best_x = np.inf, None
    for _ in range(n_starts):
        x0 = lo + (hi - lo) * rng.random(model.k)
        res = minimize(neg, x0, method="Nelder-Mead",
                       options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 4000})
        if res.fun < best_val:
            best_val, best_x = float(res.fun), res.x
    aic = best_val + 2 * model.k
    bic = best_val + model.k * np.log(n)
    best = {p: float(v) for p, v in zip(model.param_names, best_x)}
    return Fit(model.name, model.k, best_val, aic, bic, best)


def evidence(model: M.Model, data: Dataset, nlive: int = 400, seed: int = 42):
    """Bayesian log-evidence via dynesty (flat priors over model.bounds)."""
    import dynesty
    if model.k == 0:
        lnz = -0.5 * data.chi2(model, [])     # zero-dim prior volume = 1
        return lnz, 0.0
    names = model.param_names
    lo = np.array([model.bounds[p][0] for p in names])
    hi = np.array([model.bounds[p][1] for p in names])

    def loglike(theta):
        c = data.chi2(model, theta)
        return -0.5 * c if np.isfinite(c) else -1e30

    def ptform(u):
        return lo + (hi - lo) * u

    rng = np.random.default_rng(seed)
    sampler = dynesty.NestedSampler(loglike, ptform, model.k, nlive=nlive, rstate=rng)
    sampler.run_nested(print_progress=False)
    r = sampler.results
    return float(r.logz[-1]), float(r.logzerr[-1])
