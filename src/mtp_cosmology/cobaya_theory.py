"""
cobaya_theory.py — Cobaya `Theory` wrapper for the windowed-IDE background.

This exposes the model's background expansion to Cobaya so it can be driven by
Cobaya's native likelihoods (DESI BAO, Pantheon+ SN, RSD) and combined with the
Planck information. It is a scaffold for the phase_4 full-stack run; `cobaya` is
imported lazily so this module imports fine without it installed.

Scope / honest note on "full Planck":
- This wrapper provides H(z) and comoving distances — enough for BAO, SN, RSD,
  and the **compressed CMB distance prior** (R, l_A), which for a *late-time-only*
  model like ours carries the full geometric CMB information. That route is
  implemented and run in `compare.py` (phase_4) and needs no Boltzmann code.
- The **full Cℓ-level Planck likelihood** additionally requires the CMB power
  spectra, i.e. solving the cosmological perturbations *with the IDE coupling*.
  Standard CAMB/CLASS do not include this custom coupling, so a literal Cℓ run
  for the windowed-IDE model needs the coupling added to a Boltzmann code's
  perturbation module. For ΛCDM and CPL (native to CAMB) a full Cℓ Planck run is
  possible directly in Cobaya; for late-time DE it returns the same dark-energy
  information as the compressed prior used here.

Example Cobaya block:

    theory:
      mtp_cosmology.cobaya_theory:MTPWindowedIDE: {}
    params:
      beta0:  {prior: {min: 0.0, max: 0.3}}
      z_star: {prior: {min: 0.0, max: 1.5}}
      sigma:  {prior: {min: 0.05, max: 1.0}}
    likelihood:
      bao.desi_dr1: ...          # consumes Hubble / comoving distances
"""
from __future__ import annotations

import numpy as np

from .model import H0, C_KM
from . import models as M

try:
    from cobaya.theory import Theory as _CobayaTheory
    _HAVE_COBAYA = True
except Exception:                       # cobaya not installed
    _CobayaTheory = object
    _HAVE_COBAYA = False


class MTPWindowedIDE(_CobayaTheory):
    """Background-level Cobaya theory for the MTP windowed-IDE model (mtp_3p)."""

    params = {"beta0": None, "z_star": None, "sigma": None}

    def initialize(self):
        self._model = M.get("mtp_3p")

    def get_requirements(self):
        return {}

    def get_can_provide(self):
        return ["Hubble", "angular_diameter_distance", "comoving_radial_distance"]

    def _theta(self, params_values):
        return [params_values["beta0"], params_values["z_star"], params_values["sigma"]]

    def calculate(self, state, want_derived=True, **params_values):
        theta = self._theta(params_values)
        # Cobaya passes the redshifts it needs via the requirement system; here we
        # store the model so the get_* methods can evaluate on demand.
        state["theta"] = theta

    # --- provider methods (evaluate the background on requested z) ---
    def get_Hubble(self, z, units="km/s/Mpc"):
        z = np.atleast_1d(z)
        theta = self.current_state["theta"]
        H = self._model.H(z, theta)
        return H if units == "km/s/Mpc" else H / C_KM   # 1/Mpc

    def comoving_radial_distance(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        grid = np.linspace(1e-4, float(z.max()) * 1.001, 400)
        H = self._model.H(grid, self.current_state["theta"])
        integ = C_KM / H
        DM = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(grid))])
        return np.interp(z, grid, DM)

    def get_angular_diameter_distance(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        return self.comoving_radial_distance(z) / (1.0 + z)


def cobaya_available() -> bool:
    return _HAVE_COBAYA
