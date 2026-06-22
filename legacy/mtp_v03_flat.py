"""
MTP Cosmology v0.3 — Windowed IDE Toy Model
조비 권고 반영: f(z) 역할 분리 (F_hier × W_late), 후기우주 활성화 창 함수.
경로 C: DESI 직접 설명이 아닌 후기우주 toy model. w0-wa 시간변화 신호 타겟.

수식:
  Q = β(z)·H·ρ_c · F_hier(z)·W_late(z)
  F_hier(z) = [log(ℓ(z)/L_P)/log(ℓ0/L_P)]^n   (계층 철학 보존)
  W_late(z) = exp[-(z-z*)²/2σ²]                 (후기우주 활성화)
  β(z) = β0·tanh((z-z_c)/Δz)                     (방향 전환)

Sanity check: H(z), D_M(z), fσ8 (MCMC 전 단계)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad

# ── 고정 (Planck 2018) ───────────────────────────────────────────────────────
H0, Om0, Ode0, Or0 = 67.4, 0.315, 0.685, 9.2e-5
c_km = 299792.458  # km/s
L_P, l_0 = 1.616e-35, 1e26
LOG_RATIO = np.log(l_0 / L_P)

# ── 계층항 (로그 철학 보존) ──────────────────────────────────────────────────
def F_hier(z, n):
    """ℓ(z) = l_0/(1+z). 로그 계층 전달. 완만한 변화."""
    l_z = l_0 / (1 + z)
    return (np.log(l_z / L_P) / LOG_RATIO) ** n

# ── 후기우주 창 함수 ──────────────────────────────────────────────────────────
def W_late(z, z_star, sigma):
    """가우시안 창: z* 근방에서만 전달 활성화. CMB/BBN 자동 차단."""
    return np.exp(-((z - z_star)**2) / (2 * sigma**2))

# ── 방향 전환 ─────────────────────────────────────────────────────────────────
def beta_z(z, beta0, z_c, dz):
    return beta0 * np.tanh((z - z_c) / dz)

# ── 결합 Q/H ──────────────────────────────────────────────────────────────────
def Q_over_H(z, p):
    beta0, n, z_star, sigma, z_c, dz = p
    return beta_z(z, beta0, z_c, dz) * F_hier(z, n) * W_late(z, z_star, sigma)

# ── IDE 연립 ODE (rho_c 결합계수에 곱) ───────────────────────────────────────
def ide_system(z, y, p):
    rho_c, rho_de = y
    QH = Q_over_H(z, p) * rho_c
    w = -1.0
    drho_c  = (3*rho_c - QH) / (1 + z)
    drho_de = (3*(1+w)*rho_de + QH) / (1 + z)
    return [drho_c, drho_de]

def solve_densities(z_arr, p):
    sol = solve_ivp(ide_system, (0, z_arr.max()), [1.0, 1.0],
                    t_eval=z_arr, args=(p,), method='RK45', rtol=1e-9, atol=1e-11)
    return sol.y[0], sol.y[1]

def Hz(z_arr, p):
    rho_c, rho_de = solve_densities(z_arr, p)
    E2 = Om0*rho_c + Ode0*rho_de + Or0*(1+z_arr)**4
    return H0*np.sqrt(np.maximum(E2, 1e-12)), rho_c, rho_de

def Hz_LCDM(z):
    return H0*np.sqrt(Om0*(1+z)**3 + Ode0 + Or0*(1+z)**4)

def w_eff(z_arr, rho_de):
    ln_r = np.log(np.maximum(rho_de, 1e-12))
    return -1 + np.gradient(ln_r, np.log(1+z_arr))/3

# ── 공변거리 D_M(z) ───────────────────────────────────────────────────────────
def D_M(z_arr, Hz_vals):
    integ = c_km / Hz_vals
    D = np.concatenate([[0], np.cumsum(0.5*(integ[1:]+integ[:-1])*np.diff(z_arr))])
    return D

# ── 성장률 fσ8 (간이) ─────────────────────────────────────────────────────────
def growth_fs8(z_arr, Hz_vals, rho_c, sigma8_0=0.811):
    """δ'' + ... 간이 적분. Om(z) 기반 성장방정식."""
    a = 1/(1+z_arr)
    Om_z = Om0*rho_c / (Hz_vals/H0)**2
    # f ≈ Om(z)^0.55 (근사)
    f = Om_z**0.55
    # σ8(z) ∝ D(z), D normalized
    lnD = np.concatenate([[0], np.cumsum(0.5*(f[1:]/a[1:]+f[:-1]/a[:-1])*np.diff(a))])
    D = np.exp(lnD); D /= D[0] if D[0]!=0 else 1
    sigma8_z = sigma8_0 * D / D[0]
    return f * sigma8_z

# ── 시뮬레이션 ────────────────────────────────────────────────────────────────
z_arr = np.linspace(0.001, 2.5, 500)

# p = [beta0, n, z_star, sigma, z_c, dz]
configs = [
    ([0.05, 2.0, 0.5, 0.4, 0.5, 0.2], '#2196F3', 'β₀=0.05 z*=0.5'),
    ([0.10, 2.0, 0.5, 0.4, 0.5, 0.2], '#4CAF50', 'β₀=0.10 z*=0.5'),
    ([0.10, 2.0, 0.7, 0.5, 0.6, 0.3], '#FF9800', 'β₀=0.10 z*=0.7'),
    ([0.20, 2.0, 0.4, 0.3, 0.4, 0.2], '#9C27B0', 'β₀=0.20 z*=0.4'),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle('MTP Cosmology v0.3 — Windowed IDE Toy Model', fontsize=14, fontweight='bold')

# (1) w_eff(z)
ax = axes[0,0]
ax.axhline(-1, color='gray', ls='--', lw=1.5, label='Λ-CDM')
ax.axhspan(-1.1, -1.0, alpha=0.06, color='red')
for p, col, lbl in configs:
    _, _, rde = Hz(z_arr, p)
    ax.plot(z_arr, w_eff(z_arr, rde), color=col, lw=2, label=lbl)
ax.set_xlabel('z'); ax.set_ylabel('w_eff(z)'); ax.set_title('Equation of State')
ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_ylim(-1.1, -0.9)

# (2) H(z)/LCDM
ax = axes[0,1]
HL = Hz_LCDM(z_arr)
ax.axhline(1, color='gray', ls='--', lw=1.5, label='Λ-CDM')
for p, col, lbl in configs:
    h,_,_ = Hz(z_arr, p)
    ax.plot(z_arr, h/HL, color=col, lw=2, label=lbl)
ax.set_xlabel('z'); ax.set_ylabel('H(z)/H_ΛCDM'); ax.set_title('Hubble (sanity 1)')
ax.legend(fontsize=7); ax.grid(alpha=0.3)

# (3) D_M(z)/LCDM
ax = axes[1,0]
DL = D_M(z_arr, HL)
ax.axhline(1, color='gray', ls='--', lw=1.5, label='Λ-CDM')
for p, col, lbl in configs:
    h,_,_ = Hz(z_arr, p)
    dm = D_M(z_arr, h)
    ratio = np.ones_like(dm); ratio[1:] = dm[1:]/DL[1:]
    ax.plot(z_arr[1:], ratio[1:], color=col, lw=2, label=lbl)
ax.set_xlabel('z'); ax.set_ylabel('D_M(z)/D_M^ΛCDM'); ax.set_title('Comoving Distance (sanity 2)')
ax.legend(fontsize=7); ax.grid(alpha=0.3)

# (4) 창 함수 구조
ax = axes[1,1]
for p, col, lbl in configs[:2]:
    beta0,n,zs,sig,zc,dz = p
    ax.plot(z_arr, W_late(z_arr,zs,sig), color=col, lw=2, label=f'W_late {lbl}')
    ax.plot(z_arr, beta_z(z_arr,beta0,zc,dz), color=col, lw=1.2, ls='--', alpha=0.6)
ax.plot(z_arr, F_hier(z_arr,2.0), color='black', lw=1.5, ls=':', label='F_hier(n=2)')
ax.axhline(0, color='gray', lw=0.8)
ax.set_xlabel('z'); ax.set_ylabel('함수값'); ax.set_title('Components (solid=W, dash=β, dot=F)')
ax.legend(fontsize=7); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/mtp_v03.png', dpi=150, bbox_inches='tight')
plt.close()

# ── 수치 ──────────────────────────────────────────────────────────────────────
print("="*72)
print("MTP v0.3 — Windowed IDE Toy Model")
print("="*72)
print(f"{'config':<22}{'w(0)':>9}{'w_min':>9}{'z@min':>8}{'Phantom':>9}{'H(0)':>9}")
print("-"*72)
for p, col, lbl in configs:
    h, rc, rde = Hz(z_arr, p)
    we = w_eff(z_arr, rde)
    imin = np.argmin(we)
    phantom = "YES" if we.min() < -1.001 else "no"
    print(f"{lbl:<22}{we[0]:>9.4f}{we.min():>9.4f}{z_arr[imin]:>8.2f}{phantom:>9}{h[0]:>9.2f}")

print()
print("w0-wa 신호 (CPL 근사 피팅):")
for p, col, lbl in configs:
    _, _, rde = Hz(z_arr, p)
    we = w_eff(z_arr, rde)
    # w(a) = w0 + wa(1-a), a=1/(1+z)
    a = 1/(1+z_arr)
    A = np.vstack([np.ones_like(a), 1-a]).T
    mask = z_arr < 1.5
    w0, wa = np.linalg.lstsq(A[mask], we[mask], rcond=None)[0]
    print(f"  {lbl}: w0={w0:.4f}, wa={wa:.4f}")

print()
print("Sanity: fσ8(z=0.5) 비교")
for p, col, lbl in configs[:2]:
    h, rc, rde = Hz(z_arr, p)
    fs8 = growth_fs8(z_arr, h, rc)
    i = np.argmin(np.abs(z_arr-0.5))
    print(f"  {lbl}: fσ8(0.5)={fs8[i]:.4f}")
print("  (관측 fσ8(0.5)≈0.46 참고)")
