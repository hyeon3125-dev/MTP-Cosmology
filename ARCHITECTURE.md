# Architecture — MTP Cosmology

This document covers (1) the model's physics and its evolution, (2) the code
structure, and (3) the key design decisions and known limitations. It is the
canonical technical reference; `MTP_consolidated_report.md` is the original
narrative report it supersedes.

## 1. Model evolution (why the current form)

| Version | Form | Outcome |
|---------|------|---------|
| v0.1 | `Q = α·H·ρ_c·f(z)`, f = log-hierarchy ansatz | **Failed**: (1) ODE breaks ∇μTᵘᵛ=0; (2) w_eff = −1 + (positive) ⇒ phantom barrier; (3) α≡G_eff/G_N−1 violates Cassini ~100×; f(z→0)→0 pins w to −1 today. |
| v0.2 | IDE 4-vector coupling, `α(z)=α₀ tanh(z−z_c)` | **Partial**: GR-consistent, sign switching works, Cassini avoided by redefining α as a coupling. But log structure dies at low z ⇒ amplitude 1/10–1/30 of DESI requirement. |
| v0.3 | **Windowed IDE** (current) | **Works**: factor f → `F_hier·W_late`. The Gaussian window restores low-z amplitude ⇒ w<−1, w₀–wₐ signal, sanity checks pass. |

Runnable prototypes for v0.1/v0.2 live in `legacy/` and self-demonstrate their
failure modes (`python legacy/mtp_v01_prototype.py`, `…v02_ide.py`).

**Self-assessment.** `F_hier(z)` varies by <1% over z ≲ 2; setting it to 1
leaves results unchanged. The model is therefore operationally a *windowed IDE*
model. The hierarchy naming is retained for continuity, not because it carries
weight at observable redshifts. This is stated up front to avoid overclaiming.

## 2. Final equations (v0.3)

Energy exchange Q between CDM (ρ_c) and DE (ρ_de), with w = −1 for the DE fluid:

```
ρ̇_c  + 3Hρ_c        = +Q
ρ̇_de + 3H(1+w)ρ_de  = −Q
```

In redshift (using d/dt = −(1+z)H d/dz), with the dimensionless kernel
`(Q/H)/ρ_c`:

```
(1+z) dρ_c/dz  = 3ρ_c          − (Q/H)
(1+z) dρ_de/dz = 3(1+w)ρ_de    + (Q/H)
(Q/H)/ρ_c = β(z) · F_hier(z) · W_late(z)

F_hier(z) = [ log(ℓ(z)/L_P) / log(ℓ₀/L_P) ]ⁿ,   ℓ(z)=ℓ₀/(1+z)   (hierarchy)
W_late(z) = exp[ −(z−z*)² / (2σ²) ]                              (late window)
β(z)      = β₀ · tanh((z−z_c)/Δz)                                (sign/strength)
```

Friedmann (normalized densities, ρ(0)=1):

```
H²(z) = H₀² [ Ω_m ρ_c + Ω_de ρ_de + Ω_r (1+z)⁴ ]
```

Effective EoS: `w_eff(z) = −1 + (1/3) d ln ρ_de / d ln(1+z)`.

**Free parameters:** β₀, z*, σ, z_c, Δz. (n is a near-inert nuisance.)
**Fixed:** Planck 2018 base-ΛCDM — H₀=67.4, Ω_m=0.315, Ω_de=0.685, Ω_r=9.2e-5.

## 3. Code architecture

```
src/mtp_cosmology/model.py     Pure model. No I/O, no plotting, no globals beyond
                               fixed cosmology constants. Everything keyed off a
                               frozen `Params` dataclass.
  ├─ kernel:   F_hier, W_late, beta_z, Q_over_H_per_rhoc
  ├─ background: solve_densities (RK45), Hz, Hz_LCDM, w_eff
  └─ observables: D_M, growth_fs8, cpl_fit

scripts/run_v03.py             Driver: regenerates figures/ + the §3 table.
scripts/run_mcmc.py            (P1) emcee MCMC against observable likelihood.
```

**Design rules**
- The model module is import-safe and deterministic; all randomness (synthetic
  data, MCMC) lives in scripts and is seeded.
- No hardcoded absolute paths. Drivers default to repo-relative `figures/`,
  `results/`; overridable via CLI flags.
- Observables (H, D_M, fσ8) are the contract surface for likelihoods — P1/P2 fit
  these, **not** `w_eff` directly (w_eff is a derived diagnostic, not data).

## 4. Likelihood design (P1 — supersedes the toy version)

The original `legacy/advanced_mcmc_analysis_toyw.py` fit a standalone shape
`w(z) = −1 − β₀·logistic_pdf(z)`, which is **not** the v0.3 model output. Its
"parameter recovery" therefore validated a logistic curve, not the windowed-IDE
model. P1 replaces this:

```
data:        synthetic H(z), D_M(z), fσ8(z) from the model + observation-scale noise
likelihood:  −½ Σ_probe Σ_i [ (d_i − model_i(θ)) / σ_i ]²
sampler:     emcee (affine-invariant); convergence via R̂ / autocorrelation
params:      θ = (β₀, z*, σ[, z_c, Δz])
```

**P1 results (synthetic).** The β₀≈0.1 coupling is a sub-percent effect on
observables (H 0.6%, D_M 0.25%, fσ8 1.3% vs ΛCDM). Pipeline validated: fitting β₀
with the window fixed recovers it within 1σ at forecast precision
(0.091±0.016, R̂≈1.00). The full (β₀,z*,σ) fit shows a β₀–σ degeneracy — only the
integrated transfer (amplitude×width) is constrained.

## 4b. Real-data fit + evidence (P2)

`src/mtp_cosmology/{data,realfit}.py`, `scripts/run_realfit.py`. Fits real DESI
DR1 BAO (`data.py`, DESI 2024 VI Table 1) with fixed r_d=147.09 Mpc (W_late
leaves the drag-epoch sound horizon at its ΛCDM value). ΛCDM = the β₀→0 limit
(zero free params); model evidence via nested sampling (dynesty).

Result: β₀ < 0.27 (95%), Δln Z = −1.9 (β₀ only) / −1.2 (3-param) — weak
preference for ΛCDM. The coupling cannot relieve the LRG1 D_H/r_d feature at
z=0.51 (the driver of DESI's evolving-DE hint), so the extra freedom is
Occam-penalized. A real test of the w₀–wₐ signal needs the full Planck+DESI+SN
likelihood with a modified Boltzmann code (external environment).

## 4c. Fair model comparison (primary objective)

`docs/comparison_methodology.yaml` is the spec. `src/mtp_cosmology/models.py`
holds a unified registry (ΛCDM, CPL, constant IDE, sign-switching IDE, MTP-3p,
MTP-4p), each exposing one `H(z, θ)` interface; `compare.py` applies one
likelihood (H+BAO+SN; SN offset analytically marginalized), χ²_min via
multi-start optimization, AIC/BIC, and dynesty evidence identically to every
model; `scripts/run_compare.py` runs a stage and tabulates Δ vs ΛCDM.

Results (`results/model_comparison_*.csv`): phase_0 (mock) validates the engine
(CPL recovered, ΔAIC=−14.9); phase_1 (real DESI DR1 BAO) shows **CPL is the only
model beating ΛCDM on AIC+BIC**, while the windowed IDE does not improve on ΛCDM
on geometry — a perturbative coupling moves H by ≲1%, far less than a free w(a).
Cosmology is fixed to Planck (background-level); freeing H₀/Ω_m and adding
growth/CMB are later stages (phase_2–4).

## 5. Known limitations / open items

- **F_hier inert** at observable z (§1) — acknowledged, not hidden.
- **Screening is asserted, not derived** (report §4-2): a thin-shell factor
  S(ρ) is invoked to satisfy Cassini |γ−1|≲2.3×10⁻⁵, but a concrete
  scalar-tensor micro-model is P3 work.
- **fσ8 is approximate** — uses the f≈Ω_m(z)^0.55 growth-index ansatz, not a
  full perturbation solve. Adequate for sanity; revisit for real-data fits.
- **CMB/BBN** are assumed untouched because W_late suppresses Q at high z;
  this should be checked explicitly once Cobaya is wired up (P2).
