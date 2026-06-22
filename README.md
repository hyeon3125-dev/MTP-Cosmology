# MTP Cosmology — Windowed Interacting Dark Energy Toy Model

A GR-consistent **windowed interacting dark energy (IDE)** model. A late-time
activation window (centered near the dark-matter/dark-energy coincidence epoch,
z* ≈ 0.3–0.5) drives an effective **phantom crossing** (w < −1) and produces a
w₀–wₐ time-variation signal consistent in direction with DESI DR1.

> **Honest scope.** The hierarchy term `F_hier(z)` is nearly constant over the
> observable range (z ≲ 2), so the model is, operationally, a windowed IDE
> model. The "MTP / hierarchy" naming is kept for historical/philosophical
> context. See [ARCHITECTURE.md](ARCHITECTURE.md) §1.

## Quick start

```bash
pip install -r requirements.txt
python scripts/run_v03.py          # reproduces figures/mtp_v03.png + the §3 table
```

The core model is an importable, side-effect-free package:

```python
import sys; sys.path.insert(0, "src")
from mtp_cosmology import Params, Hz, w_eff
import numpy as np

z = np.linspace(0.001, 2.5, 500)
p = Params(beta0=0.10, z_star=0.5, sigma=0.4)
H, rho_c, rho_de = Hz(z, p)
print("w_eff(0) =", w_eff(z, rho_de)[0])
```

## Repository layout

```
src/mtp_cosmology/      Importable core (model.py: kernel, ODE, observables)
scripts/                Drivers: run_v03.py (figures+table), run_mcmc.py (P1)
legacy/                 v0.1/v0.2 prototypes + flat v0.3 + toy-w MCMC (history)
figures/  results/      Generated artifacts
paper/                  Paper draft (P3)
docs/                   Notes
MTP_consolidated_report.md   Original consolidated dev report (v0.1→v0.3)
ARCHITECTURE.md         Physics derivation + code architecture + decisions
worklog.md              Dated work log
```

## Status & roadmap

| Phase | Scope | State |
|-------|-------|-------|
| **P0** | Reproducible landing: importable core, fixed paths, docs | ✅ done |
| **P1** | Observable-based likelihood (H/D_M/fσ8) + emcee MCMC | ✅ done |
| **P2** | Real DESI DR1 BAO fit + Bayesian evidence vs ΛCDM | ✅ done |
| **P3** | Scalar-tensor screening micro-model + paper draft | ✅ done |

**Headline results (all honest, reproducible):**
- The coupling at β₀≈0.1 is a **sub-percent** effect on H/D_M/fσ8 — below current
  precision. Pipeline recovers injected β₀ within 1σ at forecast precision;
  β₀–σ is degenerate from background+growth data (only β₀ amplitude is
  constrained).
- **DESI DR1 BAO**: β₀ < 0.27 (95%), Δln Z = −1.9 (1-param) / −1.2 (3-param) →
  weak preference for ΛCDM. The coupling is too weak to relieve the LRG1 D_H
  feature.
- **Screening**: a chameleon with Solar field excursion < 4.3×10⁻⁸ M_Pl satisfies
  the Cassini bound.

Run: `python scripts/run_v03.py`, `python scripts/run_mcmc.py --fit beta0`,
`python scripts/run_realfit.py --fit beta0`. See [paper/paper.md](paper/paper.md).

## Model comparison (primary objective)

The point of the project is **fair model comparison**: does windowed IDE explain
late-time phantom-like dynamics and growth with comparable or fewer effective
degrees of freedom than ΛCDM, CPL, and existing IDE variants? The methodology
(model set, datasets, priors, metrics, success/failure thresholds, run matrix)
is specified in [docs/comparison_methodology.yaml](docs/comparison_methodology.yaml).

```bash
python scripts/run_compare.py --stage phase_0   # mock pipeline validation
python scripts/run_compare.py --stage phase_1   # real DESI DR1 BAO (geometry)
python scripts/run_compare.py --stage phase_3   # + Gold-2018 RSD (growth)
python scripts/run_compare.py --stage phase_4   # + Planck18 compressed CMB
```

Outputs an AIC/BIC/Δln Z table over {ΛCDM, CPL, constant IDE, sign-switching IDE,
MTP-3p, MTP-4p}, all under identical priors/likelihood/sampler.

**Verdict across geometry → growth → CMB (real data):** CPL is the only model
that beats ΛCDM on both AIC and BIC, capturing the DESI evolving-DE signal
(w₀≈−0.7, wₐ≈−1). The windowed IDE never beats ΛCDM; adding Planck's (R, l_A)
**drives the IDE coupling to zero** (D_M(z\*) pins the late-time expansion), so
β₀ ends up consistent with 0. As a GR-perturbative late-time coupling the model
is *not economical* vs CPL — see [worklog.md](worklog.md), [paper/paper.md](paper/paper.md),
and [docs/phase4_cmb.md](docs/phase4_cmb.md).

See [worklog.md](worklog.md) for the running log and
[ARCHITECTURE.md](ARCHITECTURE.md) for the model.
