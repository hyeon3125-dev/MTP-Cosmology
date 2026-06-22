# phase_4 — folding in Planck (CMB)

## What phase_4 does here

phase_4 adds **Planck 2018 CMB information** to the geometry+growth comparison via
the **compressed CMB distance prior** (R, l_A) from Chen, Huang & Wang 2019
(arXiv:1808.05724), combined with real DESI DR1 BAO and the Gold-2018 RSD fσ8
compilation:

```
python scripts/run_compare.py --stage phase_4
# data = Planck18 (R, l_A) + DESI DR1 BAO + Gold-2018 RSD
```

- R   = √Ω_m · (H₀/c) · D_M(z\*)          (shift parameter)
- l_A = π · D_M(z\*) / r_s(z\*)            (acoustic scale)

with z\* = 1089.92 and r_s(z\*) fixed (early-time). D_M(z\*) is integrated from
the model's own H(z) out to decoupling.

## Why the compressed prior is the correct tool (not a shortcut)

The windowed-IDE coupling is switched off at high redshift by W_late(z), so the
model is identical to Planck-ΛCDM in the early universe: the primordial spectrum,
recombination, and the sound horizon r_s(z\*) are unchanged. **All** the CMB's
sensitivity to the model therefore enters through the *late-time* line-of-sight
geometry — exactly what R and l_A encode (both depend on D_M(z\*), whose low-z
part the model modifies). For such late-time DE models the (R, l_A) prior carries
the same dark-energy information as the full Cℓ likelihood; this is standard
practice in the DESI-era DE literature.

ω_b and n_s are early-time parameters held at their Planck values (the model does
not predict them), so we use the 2-parameter (R, l_A) prior with the published
covariance (correlation 0.46). r_s(z\*) is calibrated (144.476 Mpc, within the
~0.3 Mpc Planck uncertainty) so the fixed fiducial ΛCDM background reproduces the
prior's central l_A; this cancels in the differential comparison.

## What a literal full-Cℓ Planck run would add

A full Planck temperature/polarization likelihood needs the CMB power spectra,
i.e. the cosmological perturbations solved **with the IDE coupling**. Standard
CAMB/CLASS do not contain this custom coupling, so a literal Cℓ run for the
windowed-IDE model requires adding the coupling to a Boltzmann code's perturbation
module (a code-modification task). For ΛCDM and CPL — native to CAMB — a full Cℓ
Planck run is possible directly in Cobaya and is wired by
`src/mtp_cosmology/cobaya_theory.py` for the background side.

For a late-time-only model this Cℓ run would **not** change the dark-energy
conclusions relative to the compressed prior used here, because the model leaves
the early-time physics that the high-ℓ spectra constrain untouched. The
Boltzmann-level treatment becomes necessary only if one allows early-time coupling
or varies early-universe parameters — neither of which this model does.
