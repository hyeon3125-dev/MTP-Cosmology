# Worklog — MTP Cosmology

Dated, append-only log of development. Newest last within each phase.

---

## 2026-06-22 — P0: reproducible landing

**Goal.** Land the consolidated dev artifacts in a clean, reproducible repo so
the model can be run, reviewed, and built on (paper track).

**Done.**
- Restructured flat scripts into an importable package + drivers:
  - `src/mtp_cosmology/{__init__,model}.py` — pure, side-effect-free core
    (kernel, background ODE, observables). Keyed off a frozen `Params` dataclass.
  - `scripts/run_v03.py` — driver reproducing `figures/mtp_v03.png` and the §3
    numerical table. Headless backend; repo-relative output (no hardcoded
    `/mnt/user-data/...` paths).
- Reconstructed the v0.1/v0.2 prototypes referenced by the report but missing
  from disk, as self-demonstrating `legacy/` scripts:
  - `legacy/mtp_v01_prototype.py` — phantom-barrier defect (min w ≈ −0.90, no crossing).
  - `legacy/mtp_v02_ide.py` — IDE transition, GR-consistent.
  - originals preserved: `mtp_v03_flat.py`, `*_toyw.py`, `mtp_v03_original.png`.
- Authored `requirements.txt` (pinned lower bounds), `.gitignore`, `README.md`,
  `ARCHITECTURE.md` (physics + code + decisions + limitations).

**Verification.** `python scripts/run_v03.py` reproduces the report §3 numbers
exactly: w(0) = −1.008/−1.015/−1.027, CPL (w₀,wₐ) = (−1.030,+0.133)/(−1.060,+0.264),
fσ8(0.5) = 0.477/0.480, H(0)=67.43 — confirming the refactor preserved the physics.

**Key finding (carried into P1).** The report's MCMC "parameter recovery" was
produced by `legacy/advanced_mcmc_analysis_toyw.py`, whose likelihood fits a
standalone `w(z) = −1 − β₀·logistic_pdf(z)` shape — **not** the windowed-IDE
model output. So §3-4 validates a logistic curve, not the model. P1 must rebuild
the likelihood on actual observables (H/D_M/fσ8). Documented in ARCHITECTURE §4.

**Next.** P1 — observable-based likelihood + emcee, then re-derive §3 numbers
honestly.
