"""MTP Cosmology — windowed interacting dark energy toy model."""
from .model import (
    Params,
    F_hier,
    W_late,
    beta_z,
    Q_over_H_per_rhoc,
    solve_densities,
    Hz,
    Hz_LCDM,
    w_eff,
    D_M,
    growth_fs8,
    cpl_fit,
    H0,
    OM0,
    ODE0,
    OR0,
    C_KM,
)

__version__ = "0.3.0"

__all__ = [
    "Params", "F_hier", "W_late", "beta_z", "Q_over_H_per_rhoc",
    "solve_densities", "Hz", "Hz_LCDM", "w_eff", "D_M", "growth_fs8",
    "cpl_fit", "H0", "OM0", "ODE0", "OR0", "C_KM", "__version__",
]
