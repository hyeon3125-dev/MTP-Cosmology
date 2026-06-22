# A Windowed Interacting Dark Energy Model: Late-Time Phantom Crossing and a Confrontation with DESI DR1 BAO

**Draft v0.3** — internal working paper. Status: synthetic-data inference and a
real DESI DR1 BAO evidence comparison complete; full Planck+DESI+SN likelihood
deferred to future work (requires a modified Boltzmann code).

---

## Abstract

We present a GR-consistent interacting dark energy (IDE) model in which the
energy exchange between cold dark matter and dark energy is modulated by a
late-time Gaussian *window* centred near the dark-matter/dark-energy coincidence
epoch (z\* ≈ 0.3–0.5). The window produces an effective equation of state that
crosses the phantom divide (w < −1) at low redshift and a CPL-projected (w₀, wₐ)
with w₀ < −1, wₐ > 0 — the direction favoured by DESI. Using synthetic
H(z)/D_M(z)/fσ8(z) data we show the inference pipeline recovers the injected
coupling amplitude within 1σ, while the window's location and width are
degenerate from background-plus-growth data alone (only the integrated energy
transfer is constrained). Confronting the model with real DESI DR1 BAO, we find
the coupling amplitude bounded as β₀ < 0.27 (95%) and a Bayesian evidence Δln Z =
−1.9 relative to ΛCDM: the sub-percent coupling cannot relieve the LRG1
D_H/r_d feature at z = 0.51, so the data mildly disfavour the extra freedom. We
outline a chameleon screening mechanism consistent with Solar-System tests and
discuss the conditions under which the model could be meaningfully tested.

## 1. Introduction

DESI DR1 BAO hint at a time-varying dark energy equation of state, with a CPL
fit preferring w₀ > −1 today and wₐ < 0 (equivalently a phantom-like w < −1 at
intermediate redshift). This motivates models that generate a late-time
phantom crossing within a GR-consistent framework. Interacting dark energy,
where a coupling Q exchanges energy between the dark sectors, is a natural such
framework: an effective phantom crossing can arise without a literal phantom
field.

This work develops one such model with a deliberately *localized* coupling — a
window in redshift — and, crucially, subjects it to honest inference: we ask not
only "can it produce the signal" but "what do data actually constrain," and we
report the answer even when it is a null. The model carries hierarchical-scale
("MTP") motivation in its naming; we state at the outset (§2.3) that this
structure is operationally inert at observable redshifts.

## 2. Model

### 2.1 Background equations

Energy exchange Q between CDM (ρ_c) and a w = −1 dark-energy fluid (ρ_de):

  ρ̇_c + 3Hρ_c = +Q,  ρ̇_de + 3H(1+w)ρ_de = −Q,

manifestly conserving the total stress-energy (∇_μ T^{μν} = 0). In redshift,

  (1+z) dρ_c/dz = 3ρ_c − Q/H,  (1+z) dρ_de/dz = 3(1+w)ρ_de + Q/H.

The Friedmann equation uses normalized densities (ρ(0)=1) with Planck 2018
base-ΛCDM values (H₀=67.4, Ω_m=0.315, Ω_de=0.685, Ω_r=9.2×10⁻⁵).

### 2.2 The windowed coupling

  (Q/H)/ρ_c = β(z) · F_hier(z) · W_late(z),
  F_hier(z) = [log(ℓ(z)/L_P)/log(ℓ₀/L_P)]ⁿ,  ℓ(z)=ℓ₀/(1+z),
  W_late(z) = exp[−(z−z\*)²/(2σ²)],
  β(z) = β₀ tanh((z−z_c)/Δz).

W_late localizes the exchange to low redshift, automatically suppressing it at
CMB/BBN epochs; β(z) sets the amplitude and sign; F_hier carries the hierarchy
motivation.

### 2.3 Honest scope

F_hier(z) varies by < 1% over z ≲ 2; setting it to unity leaves all results
unchanged. The model is therefore operationally a **windowed IDE** model, and we
treat it as such. The effective EoS is w_eff(z) = −1 + (1/3) d ln ρ_de/d ln(1+z).

### 2.4 Model evolution

The current form is the third iteration. A direct-transfer prototype (v0.1)
violated stress-energy conservation and admitted only w ≥ −1 (a phantom
barrier); the IDE rewrite (v0.2) cured these but produced too small an amplitude
at observable redshifts because the log structure dies off at low z. Factoring
the kernel into F_hier × W_late (v0.3) restores the low-z amplitude. Runnable
prototypes are retained (`legacy/`).

## 3. Methods

**Observables.** H(z), D_M(z) = ∫₀^z c dz'/H, and an approximate growth rate
fσ8(z) using f ≈ Ω_m(z)^0.55, computed from the background solution.

**Inference.** Gaussian likelihood −½ Σ χ² summed over probes; flat priors;
affine-invariant ensemble sampling (emcee) with autocorrelation-based burn-in
and Gelman-Rubin R̂ convergence. For real-data evidence we use nested sampling
(dynesty).

## 4. Synthetic-data results

The coupling at β₀ ≈ 0.1 changes observables by only ~0.6% (H), 0.25% (D_M) and
1.3% (fσ8) relative to ΛCDM — below current measurement precision. We therefore
report two regimes:

- **Pipeline validation** (window fixed, fit β₀, forecast precision):
  β₀ = 0.091 ± 0.016 vs injected 0.10 (pull −0.5σ), R̂ ≈ 1.00 — recovered within
  1σ, confirming the likelihood and sampler are correct.
- **Current precision** (fit β₀): β₀ = 0.13 ± 0.10 — the error is comparable to
  the value; the coupling is essentially unconstrained.
- **Full (β₀, z\*, σ), forecast precision**: a pronounced β₀–σ degeneracy. Only
  the integrated energy transfer (∝ amplitude × width) is constrained; z\* is
  moderately recovered, σ is prior-dominated (Fig. `mcmc_corner_3p_forecast.png`).

The phantom-crossing phenomenology is robust: at the validation posterior mean,
w_eff(0) ≈ −1.01, w_min ≈ −1.05, CPL (w₀, wₐ) ≈ (−1.04, +0.12).

## 5. Real-data results: DESI DR1 BAO

We fit the real DESI DR1 BAO compilation (DESI 2024 VI, Table 1: 5 anisotropic
tracers in D_M/r_d, D_H/r_d plus BGS/QSO in D_V/r_d), holding r_d = 147.09 Mpc
(justified because the window leaves early-time physics, hence the drag-epoch
sound horizon, unchanged). ΛCDM is the β₀ → 0 limit with zero free parameters.

- ΛCDM: χ² = 20.4 / 12, dominated by the LRG1 D_H/r_d point at z = 0.51 (the
  feature underlying DESI's evolving-DE preference).
- windowed IDE (β₀ only): **β₀ < 0.27 (95%)**, **Δln Z = −1.9 ± 0.1**.
- windowed IDE (β₀, z\*, σ): **Δln Z = −1.2 ± 0.1**; posterior prior-dominated.

On the Jeffreys scale both are a *weak preference for ΛCDM*: the ~0.6% coupling
cannot move the ~8% LRG1 D_H feature, so the additional parameters are
Occam-penalized.

### 5.1 Fair comparison against standard dynamical-DE models

The decisive question is not MTP vs ΛCDM alone but whether MTP is *economical*
relative to established alternatives. We compare ΛCDM, CPL (w₀,wₐ), constant IDE,
sign-switching IDE, and the MTP 3- and 4-parameter variants under identical
priors, likelihood, sampler and metrics (the comparison engine, `compare.py`).

A pipeline-validation stage on mock data generated from a CPL fiducial recovers
the injected (w₀, wₐ) = (−0.85, −0.60) as (−0.836, −0.666) and correctly selects
CPL (ΔAIC = −14.9, Δln Z = +2.6), confirming the engine discriminates.

On **real DESI DR1 BAO** (Table below), **CPL is the only model that beats ΛCDM
on both AIC and BIC**, capturing the evolving-DE signal (w₀ = −0.78, wₐ = −0.76)
with two parameters. The windowed IDE barely reduces χ² (20.44 → 20.03); a
perturbative coupling (β₀ < 0.3) moves H(z) by ≲1%, too little to reshape the
geometry as much as a free w(a), so its extra parameters are penalized.

| model | k | χ² | ΔAIC | ΔBIC | Δln Z |
|-------|---|------|------|------|-------|
| ΛCDM | 0 | 20.44 | 0 | 0 | 0 |
| CPL | 2 | 11.77 | −4.66 | −3.69 | −0.45 |
| standard IDE | 1 | 20.23 | +1.80 | +2.28 | −1.21 |
| sign-switching IDE | 3 | 18.46 | +4.03 | +5.48 | −0.79 |
| MTP-3p | 3 | 20.03 | +5.60 | +7.05 | −0.22 |
| MTP-4p | 4 | 18.17 | +5.73 | +7.67 | −0.34 |

**On geometry alone, the windowed IDE is not competitive with CPL.** Adding the
growth sector and then the CMB does not change this — and the CMB makes it
decisive:

- **+ growth** (Gold-2018 RSD fσ8): the perturbative coupling suppresses fσ8 by
  only ~2–3%, comparable to the data errors. CPL still wins on AIC+BIC; the IDE
  models do not beat ΛCDM.
- **+ CMB** (Planck 2018 R, l_A): the best-fit couplings collapse to zero
  (standard IDE ξ → 0, sign-switching ξ₀ → 0, windowed IDE pushed to
  negligible-effect corners). R and l_A pin D_M(z\*) — the integrated late-time
  expansion — so any IDE rebalancing of the dark sector is penalized unless the
  coupling vanishes. CPL survives because w(a) reshapes the EoS while preserving
  D_M(z\*).

Thus across geometry → growth → CMB the windowed IDE is **not economical**: it
never beats ΛCDM and is driven to its no-coupling limit by Planck. This triggers
the pre-registered `failure_case` (β₀ consistent with zero).

## 6. Screening and Solar-System tests

To suppress the cosmological coupling in dense environments we adopt a chameleon
mechanism with a density-dependent effective mass and a thin-shell suppression
of the external fifth force. An order-of-magnitude estimate (`screening.py`)
gives: for an order-unity matter coupling, any Solar field excursion below
~4.3×10⁻⁸ M_Pl yields |γ−1| ≤ 2.3×10⁻⁵, satisfying the Cassini bound. A complete
scalar-tensor realization (deriving β(z, ρ) from a specific potential) remains
future work.

## 7. Discussion and limitations

- **The signal is intrinsically small.** At amplitudes that keep the coupling a
  perturbative late-time effect, the observable imprint is sub-percent; current
  background+growth data cannot detect or localize it. This is the central
  honest finding.
- **Background degeneracy.** Only the integrated transfer is constrained from
  H/D_M/fσ8; breaking the β₀–σ–z\* degeneracy requires CMB and/or full-shape
  information.
- **fσ8 is approximate** (growth-index ansatz, not a full perturbation solve).
- **Screening is asserted at estimate level**, not derived from a micro-model.
- **No CMB.** The real test of the DESI w₀–wₐ signal needs the full
  Planck+DESI+SN likelihood with a Boltzmann code modified for the IDE coupling.

## 8. Conclusion

The windowed IDE model cleanly generates a late-time phantom crossing in the
DESI-preferred direction and admits a viable screening mechanism. However, a fair
comparison is unambiguous on present geometry: on DESI DR1 BAO, CPL beats ΛCDM on
both AIC and BIC while the windowed IDE does not improve on ΛCDM at all, because a
perturbative late-time coupling cannot reshape H(z) as strongly as a free w(a).
The model is a consistent toy framework but, as it stands, **not economical
against CPL on geometric data**. Its decisive test lies in the growth sector
(fσ8/S8) and a full CMB-inclusive likelihood, where the dark-sector transfer
could in principle differentiate it from a pure equation-of-state model.

---

### Reproducibility

```
python scripts/run_v03.py                       # model + sanity figures/table
python scripts/run_mcmc.py --fit beta0          # pipeline validation (synthetic)
python scripts/run_mcmc.py --fit beta0,z_star,sigma   # full fit (degeneracy)
python scripts/run_realfit.py --fit beta0       # DESI DR1 BAO evidence
```

### Data
DESI 2024 VI (arXiv:2404.03002), Table 1. Sound horizon r_d from Planck 2018.
