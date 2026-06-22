"""
run_v03.py — reproduce the v0.3 sanity figures and the §3 numerical table.

Usage:
    python scripts/run_v03.py [--out FIG_DIR]

Side-effect-free model lives in src/mtp_cosmology/model.py; this is the driver.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from mtp_cosmology import (  # noqa: E402
    Params, Hz, Hz_LCDM, w_eff, D_M, growth_fs8, cpl_fit,
    W_late, beta_z, F_hier,
)

# config = (Params, color, label)
CONFIGS = [
    (Params(0.05, 0.5, 0.4, 0.5, 0.2, 2.0), "#2196F3", "b0=0.05 z*=0.5"),
    (Params(0.10, 0.5, 0.4, 0.5, 0.2, 2.0), "#4CAF50", "b0=0.10 z*=0.5"),
    (Params(0.10, 0.7, 0.5, 0.6, 0.3, 2.0), "#FF9800", "b0=0.10 z*=0.7"),
    (Params(0.20, 0.4, 0.3, 0.4, 0.2, 2.0), "#9C27B0", "b0=0.20 z*=0.4"),
]

Z = np.linspace(0.001, 2.5, 500)


def make_figure(out_dir: str) -> str:
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle("MTP Cosmology v0.3 — Windowed IDE Toy Model",
                 fontsize=14, fontweight="bold")

    ax = axes[0, 0]
    ax.axhline(-1, color="gray", ls="--", lw=1.5, label="LCDM")
    ax.axhspan(-1.1, -1.0, alpha=0.06, color="red")
    for p, col, lbl in CONFIGS:
        _, _, rde = Hz(Z, p)
        ax.plot(Z, w_eff(Z, rde), color=col, lw=2, label=lbl)
    ax.set_xlabel("z"); ax.set_ylabel("w_eff(z)"); ax.set_title("Equation of State")
    ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_ylim(-1.1, -0.9)

    ax = axes[0, 1]
    hl = Hz_LCDM(Z)
    ax.axhline(1, color="gray", ls="--", lw=1.5, label="LCDM")
    for p, col, lbl in CONFIGS:
        h, _, _ = Hz(Z, p)
        ax.plot(Z, h / hl, color=col, lw=2, label=lbl)
    ax.set_xlabel("z"); ax.set_ylabel("H(z)/H_LCDM"); ax.set_title("Hubble (sanity 1)")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    ax = axes[1, 0]
    dl = D_M(Z, hl)
    ax.axhline(1, color="gray", ls="--", lw=1.5, label="LCDM")
    for p, col, lbl in CONFIGS:
        h, _, _ = Hz(Z, p)
        dm = D_M(Z, h)
        ratio = np.ones_like(dm); ratio[1:] = dm[1:] / dl[1:]
        ax.plot(Z[1:], ratio[1:], color=col, lw=2, label=lbl)
    ax.set_xlabel("z"); ax.set_ylabel("D_M(z)/D_M^LCDM"); ax.set_title("Comoving Distance (sanity 2)")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    for p, col, lbl in CONFIGS[:2]:
        ax.plot(Z, W_late(Z, p.z_star, p.sigma), color=col, lw=2, label=f"W_late {lbl}")
        ax.plot(Z, beta_z(Z, p.beta0, p.z_c, p.dz), color=col, lw=1.2, ls="--", alpha=0.6)
    ax.plot(Z, F_hier(Z, 2.0), color="black", lw=1.5, ls=":", label="F_hier(n=2)")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel("z"); ax.set_ylabel("value"); ax.set_title("Components (solid=W, dash=beta, dot=F)")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "mtp_v03.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    return path


def print_table() -> None:
    print("=" * 72)
    print("MTP v0.3 — Windowed IDE Toy Model")
    print("=" * 72)
    print(f"{'config':<18}{'w(0)':>9}{'w_min':>9}{'z@min':>8}{'Phantom':>9}{'H(0)':>9}")
    print("-" * 72)
    for p, _, lbl in CONFIGS:
        h, rc, rde = Hz(Z, p)
        we = w_eff(Z, rde)
        imin = int(np.argmin(we))
        phantom = "YES" if we.min() < -1.001 else "no"
        print(f"{lbl:<18}{we[0]:>9.4f}{we.min():>9.4f}{Z[imin]:>8.2f}{phantom:>9}{h[0]:>9.2f}")

    print("\nw0-wa signal (CPL fit, z<1.5):")
    for p, _, lbl in CONFIGS:
        _, _, rde = Hz(Z, p)
        w0, wa = cpl_fit(Z, w_eff(Z, rde))
        print(f"  {lbl}: w0={w0:.4f}, wa={wa:.4f}")

    print("\nSanity: f*sigma8(z=0.5)  (obs ~ 0.46)")
    for p, _, lbl in CONFIGS[:2]:
        h, rc, rde = Hz(Z, p)
        fs8 = growth_fs8(Z, h, rc)
        i = int(np.argmin(np.abs(Z - 0.5)))
        print(f"  {lbl}: f*sigma8(0.5)={fs8[i]:.4f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    default_out = os.path.join(os.path.dirname(__file__), "..", "figures")
    ap.add_argument("--out", default=default_out, help="figure output directory")
    args = ap.parse_args()
    path = make_figure(args.out)
    print_table()
    print(f"\nFigure written: {os.path.relpath(path)}")


if __name__ == "__main__":
    main()
