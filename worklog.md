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

---

## 2026-06-22 — P1: observable-based inference (fixes the core defect)

**Goal.** Replace the toy `w(z)=-1-beta0*logistic_pdf` likelihood with one driven
by the *actual* model observables H(z)/D_M(z)/fsigma8(z), then re-assess parameter
recovery honestly. (`src/mtp_cosmology/likelihood.py`, `scripts/run_mcmc.py`.)

**Done.**
- `likelihood.py`: `observables()` from the IDE background solve; Gaussian
  chi^2 over H/D_M/fsigma8 (single integration per call, ~4 ms); flat priors;
  emcee-ready `log_posterior`. Two synthetic-data precision presets (`current`,
  `forecast`).
- `run_mcmc.py`: emcee with autocorrelation-based burn-in/thinning, Gelman-Rubin
  R-hat across walkers, recovery table (pull in sigma), corner/posterior plot,
  derived w_eff(0)/CPL at the posterior mean. `--fit` selects any parameter
  subset (others fixed at truth).

**Key physical finding — the signal is sub-percent.** The coupling at beta0~0.1
moves observables by only H 0.6%, D_M 0.25%, fsigma8 1.3% vs LambdaCDM — BELOW
current measurement errors. So with current-precision data the parameters are
essentially unconstrained; this is an honest detectability statement, not a bug.

**Results (synthetic, seed 42).**
- **Pipeline validation** — fit beta0 alone (window fixed), forecast precision:
  beta0 = 0.0912 ± 0.0163 (truth 0.10, pull −0.54), R-hat 0.998 → recovered
  within 1 sigma. The likelihood + sampler are correct.
- **Current precision**, beta0 alone: beta0 = 0.132 ± 0.096 — error ≈ value,
  ~6x weaker than forecast; weakly constrained as expected.
- **Full 3-param** (beta0,z*,sigma), forecast: strong **beta0–sigma degeneracy**
  (only the *integrated* energy transfer ∝ amplitude×width is constrained); z*
  moderately recovered (~0.37 vs 0.40), sigma prior-dominated. This is physics,
  not a fitting failure — see `figures/mcmc_corner_3p_forecast.png`.

**Supersedes report §3-4.** The earlier "recovery within 1 sigma" validated a
logistic curve, not the model. The honest statement: the *amplitude* (beta0) is
recoverable; the window's location/width are degenerate from background+growth
data alone, motivating CMB / full-shape inputs (P2+).

Artifacts: `figures/mcmc_corner_{beta0_forecast,beta0_current,3p_forecast}.png`,
`results/mcmc_summary_*.csv` (raw chains gitignored, regenerable).

---

## 2026-06-22 — P2: real-data fit + Bayesian evidence (DESI DR1 BAO)

**Goal.** Confront the model with real data and compare its Bayesian evidence to
LambdaCDM. Full Planck+CMB requires a modified Boltzmann stack (out of scope in
this environment); the achievable, honest test is real **DESI DR1 BAO**
(`src/mtp_cosmology/{data,realfit}.py`, `scripts/run_realfit.py`).

**Setup.** DESI 2024 VI Table 1 (verified against the arXiv HTML source): 5
anisotropic tracers (D_M/r_d, D_H/r_d with correlation) + BGS/QSO (D_V/r_d) = 12
points. r_d = 147.09 Mpc fixed (justified: W_late suppresses Q at high z, so the
drag-epoch sound horizon is the LambdaCDM value). LambdaCDM = the beta0→0 limit
with the same fixed Planck background → **zero free parameters**, so its evidence
is the likelihood at that point. Model evidence via nested sampling (dynesty).

**Results.**
- LambdaCDM: chi^2 = 20.44 / 12 — dominated by the well-known LRG1 D_H/r_d point
  at z=0.51 (pred 22.73 vs obs 20.98), the same feature that drives DESI's
  evolving-DE preference.
- windowed-IDE, beta0 only: **beta0 < 0.27 (95%)**, Delta ln Z = **−1.87 ± 0.11**.
- windowed-IDE, (beta0,z*,sigma): Delta ln Z = **−1.18 ± 0.08**; posterior
  prior-dominated.
- Both: **weak preference for LambdaCDM** (Jeffreys). The coupling moves H by
  ~0.6% — far too little to relieve an ~8% / few-sigma LRG1 D_H feature — so the
  extra freedom is Occam-penalized.

**Honest conclusion.** DESI DR1 BAO neither requires nor excludes the windowed-IDE
coupling; it bounds beta0 from above and mildly disfavors the model on evidence.
A genuine test of the DESI w0–wa signal needs the full Planck+DESI+SN likelihood
with a modified Boltzmann code (real Step 2, external environment).

Artifacts: `results/realfit_desi_*.csv`.

---

## 2026-06-22 — P3: screening micro-model + paper draft

- `src/mtp_cosmology/screening.py`: chameleon thin-shell estimate. For an
  order-unity matter coupling, a Solar field excursion below ~4.3e-8 M_Pl yields
  |gamma−1| ≤ 2.3e-5 (Cassini). Demonstrates the cosmological coupling is
  compatible with Solar-System tests; a full scalar-tensor solve remains future
  work (honest caveat retained).
- `paper/paper.md`: draft assembling model → synthetic validation → DESI
  evidence → screening → honest limitations.

---

## 2026-06-22 — Reorientation: fair model-comparison framework

**Why.** The earlier P2 was too narrow — a single MTP-vs-ΛCDM evidence test on
DESI. The actual objective (per `docs/comparison_methodology.yaml`) is a *fair
comparison* of windowed IDE against ΛCDM, CPL, constant IDE and sign-switching
IDE under identical priors/likelihood/sampler/metrics. Also removed the internal
persona "collaboration roles" from README/paper/report (not appropriate in a
scientific repo).

**Built.**
- `src/mtp_cosmology/models.py` — unified model registry with one H(z, θ)
  interface: `lcdm` (0 DE params), `cpl_w0wa`, `standard_ide` (Q=ξHρ_c),
  `sign_switching_ide` (ξ tanh), `mtp_3p` (zc=z*, dz=σ), `mtp_4p`. Shared
  `solve_ide(kernel)`; cosmology fixed to Planck (background-level comparison).
- `src/mtp_cosmology/compare.py` — one likelihood for H(z)+BAO(D_M/D_H/D_V)+SN
  (SN offset analytically marginalized); χ²_min via multi-start Nelder-Mead;
  AIC/BIC; dynesty evidence. Real DESI DR1 dataset + mock generator.
- `scripts/run_compare.py` — runs a methodology stage, prints/saves an
  AIC/BIC/Δln Z table referenced to ΛCDM.
- `docs/comparison_methodology.yaml` — the spec (saved as requested).

**phase_0 (mock from CPL, pipeline validation).** Engine recovers the injected
fiducial: CPL best-fit (w0,wa)=(−0.836,−0.666) vs injected (−0.85,−0.60),
ΔAIC=−14.9, Δln Z=+2.6 (moderate→model). mtp_3p partially mimics it
(ΔAIC=−4.5). The comparison machinery works and discriminates.

**phase_1 (real DESI DR1 BAO).** Key result:

| model | k | χ² | ΔAIC | ΔBIC | Δln Z |
|-------|---|------|------|------|-------|
| ΛCDM | 0 | 20.44 | 0 | 0 | 0 |
| **CPL** | 2 | **11.77** | **−4.66** | **−3.69** | −0.45 |
| standard_ide | 1 | 20.23 | +1.80 | +2.28 | −1.21 |
| sign_switching_ide | 3 | 18.46 | +4.03 | +5.48 | −0.79 |
| mtp_3p | 3 | 20.03 | +5.60 | +7.05 | −0.22 |
| mtp_4p | 4 | 18.17 | +5.73 | +7.67 | −0.34 |

**Verdict (honest, answers the methodology's main_question).** On DESI DR1 BAO
geometry, **CPL is the only model that beats ΛCDM on both AIC and BIC** — it
captures the known evolving-DE signal (w0=−0.78, wa=−0.76) with 2 parameters.
The windowed IDE barely reduces χ² (20.44→20.03): a *perturbative* late-time
coupling (β0<0.3) moves H by ≲1%, far too little to reshape geometry as much as a
free w(a). It triggers the methodology's `reject_condition`/`failure_case`
(β0 effectively unconstrained, σ rails to its bound, no AIC/BIC gain).

**Implication / path forward.** Geometry alone does not favour MTP. Its only
plausible edge over CPL is the **growth sector** (fσ8/S8), where IDE imprints a
distinctive dark-sector energy-transfer signature that a pure w(a) model lacks —
i.e. `phase_3`. That, plus the full Planck+SN likelihood (`phase_4`), is where
the comparison must go for a decisive answer. Artifacts:
`results/model_comparison_{phase_0,phase_1}.csv`.

---

## 2026-06-22 — phase_3 (growth) and phase_4 (CMB): the decisive verdict

Extended the engine to the growth sector and folded in Planck.

**New machinery.**
- `growth.py` — linear growth ODE (δ'' + (2+dlnH/dN)δ' − 1.5 Ω_m(a)δ = 0) on each
  model's actual H(a) and Ω_m(a); IDE models carry the coupling's effect on the
  background matter density. fσ8(z) with σ8(0) fixed (Planck) for all models.
  Approximation (documented): background-modified growth, no explicit
  dark-sector perturbation; quasi-static. Full perturbed-IDE growth is a
  Boltzmann task.
- `data.py` — Gold-2018 RSD fσ8 compilation (22 pts, arXiv:1806.10822) and the
  Planck 2018 compressed-CMB distance prior (R, l_A; arXiv:1808.05724).
- `compare.py` — RSD block and a CMB (R, l_A) block. R, l_A use D_M(z\*) from the
  model background; the high-z tail (z>12, coupling dead) is analytic, so the CMB
  integral is cheap and exact for this late-time model class.
- `docs/phase4_cmb.md` explains why the compressed prior is the *correct* Planck
  tool here (not a shortcut), and what a literal Cℓ run would/wouldn't add.
  `cobaya_theory.py` + `configs/*.yaml` scaffold the full-stack route.

**phase_3 — DESI DR1 BAO + Gold-2018 RSD (34 pts):**

| model | k | χ² | ΔAIC | ΔBIC | Δln Z |
|-------|---|------|------|------|-------|
| ΛCDM | 0 | 39.27 | 0 | 0 | 0 |
| **CPL** | 2 | 28.56 | −6.71 | −3.65 | +0.56 |
| standard_ide | 1 | 36.48 | −0.79 | +0.74 | −0.16 |
| sign_switching_ide | 3 | 33.77 | +0.50 | +5.48 | −0.15 |
| mtp_3p | 3 | 36.59 | +3.33 | +7.90 | −0.12 |
| mtp_4p | 4 | 35.40 | +4.13 | +10.24 | −0.04 |

Adding growth does not rescue the IDE models: a perturbative coupling suppresses
fσ8 by only ~2–3% (≈ RSD errors). CPL still wins on AIC+BIC; MTP does not beat ΛCDM.

**phase_4 — Planck18 (R, l_A) + DESI DR1 BAO + Gold-2018 RSD (36 pts):**

| model | k | χ² | ΔAIC | ΔBIC | Δln Z |
|-------|---|------|------|------|-------|
| ΛCDM | 0 | 39.42 | 0 | 0 | 0 |
| **CPL** | 2 | 29.11 | −6.31 | −3.14 | −1.69 |
| standard_ide | 1 | 39.42 | +2.00 | +3.58 | −4.93 |
| sign_switching_ide | 3 | 38.83 | +5.41 | +10.16 | −3.09 |
| mtp_3p | 3 | 39.36 | +5.94 | +10.69 | −1.78 |
| mtp_4p | 4 | 37.93 | +6.50 | +12.84 | −2.60 |

(χ² shifts by ~0.1 vs the fit-only run because the growth grid was reduced
600→300 for speed; the best-fits and verdict are unchanged.)

**The CMB is decisive against the IDE coupling.** Best-fit couplings collapse:
standard_ide ξ → −0.0000, sign_switching ξ0 → −0.048, MTP pushed to
negligible-effect corners (σ rails to its lower bound 0.05). Reason: R and l_A
pin D_M(z\*) — the integrated late-time expansion — so any IDE shift of the
matter/DE balance is penalized unless the coupling vanishes.

Note the AIC/BIC vs evidence split: **CPL wins on AIC/BIC** (the information
criteria used in DESI-style analyses) by capturing the evolving-DE signal, but on
**Bayesian evidence every extended model — CPL included — is mildly disfavoured**
(Δln Z < 0), because the broad w0/wa and coupling priors incur an Occam penalty
the modest χ² gains don't overcome. The IDE models lose on *all* metrics.

**Verdict (answers the methodology's main_question).** Across geometry, growth,
and CMB, the windowed-IDE model is **not economical**: it never beats ΛCDM on
AIC/BIC, loses to CPL throughout, and with Planck folded in its coupling is driven
to zero — triggering the spec's `failure_case` ("β0 consistent with 0"). CPL is
the only model favoured over ΛCDM, capturing the known DESI evolving-DE signal
(w0≈−0.7, wa≈−1). The windowed-IDE framework, as a GR-perturbative late-time
coupling, is disfavoured by current data. A literal Cℓ-level Planck run for the
IDE model (modified Boltzmann code) is the only remaining refinement, and for a
late-time model would not change this conclusion (see `docs/phase4_cmb.md`).

Artifacts: `results/model_comparison_{phase_3,phase_4}.csv`.
