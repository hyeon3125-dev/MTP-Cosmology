# Model equations (comparison set)

All models share the Planck 2018 background (H₀=67.4, Ω_m=0.315, Ω_de=0.685,
Ω_r=9.2e-5); only the dark-sector parameters listed are free (background-level
comparison, see `docs/comparison_methodology.yaml`). Implemented in
`src/mtp_cosmology/models.py`.

## ΛCDM (k=0)
H²(z) = H₀²[Ω_m(1+z)³ + Ω_de + Ω_r(1+z)⁴].

## CPL — w₀wₐ (k=2)
w(a) = w₀ + wₐ(1−a),  a = 1/(1+z).
ρ_de(z)/ρ_de0 = (1+z)^{3(1+w₀+wₐ)} · exp[−3 wₐ z/(1+z)].
H² = H₀²[Ω_m(1+z)³ + Ω_de · ρ_de/ρ_de0 + Ω_r(1+z)⁴].

## Interacting DE (continuity, w_DE = −1)
(1+z) dρ_c/dz = 3ρ_c − (Q/H)/ρ_c · ρ_c,
(1+z) dρ_de/dz = +(Q/H)/ρ_c · ρ_c,
H² = H₀²[Ω_m ρ_c + Ω_de ρ_de + Ω_r(1+z)⁴]  (normalized ρ(0)=1).

- **standard IDE** (k=1): (Q/H)/ρ_c = ξ.
- **sign-switching IDE** (k=3): (Q/H)/ρ_c = ξ₀ tanh((z−z_c)/Δz).
- **MTP windowed IDE**: (Q/H)/ρ_c = β(z)·F_hier(z)·W_late(z),
  β(z)=β₀ tanh((z−z_c)/Δz), W_late=exp[−(z−z\*)²/2σ²],
  F_hier=[log(ℓ(z)/L_P)/log(ℓ₀/L_P)]ⁿ (≈const for z≲2).
  - **MTP-3p** (k=3): z_c=z\*, Δz=σ → free {β₀, z\*, σ}.
  - **MTP-4p** (k=4): Δz=σ → free {β₀, z\*, σ, z_c}.

## Derived observables
D_M(z)=∫₀^z c dz'/H,  D_H=c/H,  D_V=[z D_M² D_H]^{1/3},  D_L=(1+z)D_M,
μ=5log₁₀(D_L/10pc),  fσ8 via f≈Ω_m(z)^0.55.  BAO ratios use r_d=147.09 Mpc.
