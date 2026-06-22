"""
run_compare.py — fair model comparison over a methodology stage.

Fits every model in the stage's model set to the same dataset with the same
likelihood/sampler/metrics, and prints a table of chi2, k, AIC, BIC and Bayesian
evidence, all referenced to LambdaCDM (Delta vs LCDM). See
docs/comparison_methodology.yaml.

Usage:
    python scripts/run_compare.py --stage phase_0   # mock pipeline validation
    python scripts/run_compare.py --stage phase_1   # real DESI DR1 BAO
    python scripts/run_compare.py --models lcdm,cpl_w0wa,mtp_3p --data desi
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from mtp_cosmology import models as M           # noqa: E402
from mtp_cosmology import compare as C           # noqa: E402

STAGES = {
    "phase_0": dict(
        models=["lcdm", "cpl_w0wa", "mtp_3p"],
        data="mock",
        goal="pipeline validation (mock from CPL fiducial)",
    ),
    "phase_1": dict(
        models=["lcdm", "cpl_w0wa", "standard_ide", "sign_switching_ide", "mtp_3p", "mtp_4p"],
        data="desi",
        goal="late-time geometry (real DESI DR1 BAO)",
    ),
    "phase_3": dict(
        models=["lcdm", "cpl_w0wa", "standard_ide", "sign_switching_ide", "mtp_3p", "mtp_4p"],
        data="desi+rsd",
        goal="geometry + growth (DESI DR1 BAO + Gold-2018 RSD fsigma8)",
    ),
    "phase_4": dict(
        models=["lcdm", "cpl_w0wa", "standard_ide", "sign_switching_ide", "mtp_3p", "mtp_4p"],
        data="full",
        goal="Planck18 compressed CMB + DESI DR1 BAO + Gold-2018 RSD",
    ),
}


def make_data(kind: str) -> C.Dataset:
    if kind == "mock":
        return C.mock_dataset(fiducial="cpl_w0wa", theta=(-0.85, -0.60), seed=42)
    if kind == "desi":
        return C.desi_dr1_dataset()
    if kind == "rsd":
        return C.rsd_dataset()
    if kind == "desi+rsd":
        return C.desi_plus_rsd_dataset()
    if kind == "full":
        return C.full_dataset()
    raise ValueError(kind)


def jeffreys(d):
    a = abs(d)
    if a < 1: return "~"
    if a < 2.5: return "weak"
    if a < 5: return "moderate"
    return "strong"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=list(STAGES))
    ap.add_argument("--models", help="comma list overriding the stage model set")
    ap.add_argument("--data", help="mock|desi, overriding the stage dataset")
    ap.add_argument("--no-evidence", action="store_true", help="skip nested sampling")
    ap.add_argument("--nlive", type=int, default=400)
    ap.add_argument("--res-dir", default=os.path.join(os.path.dirname(__file__), "..", "results"))
    args = ap.parse_args()

    if args.stage:
        spec = STAGES[args.stage]
        model_names = args.models.split(",") if args.models else spec["models"]
        data_kind = args.data or spec["data"]
        tag = args.stage
        goal = spec["goal"]
    else:
        assert args.models and args.data, "provide --stage or both --models and --data"
        model_names = args.models.split(",")
        data_kind = args.data
        tag = "custom"
        goal = "custom"

    data = make_data(data_kind)
    print(f"Stage: {tag} | {goal}")
    print(f"Data: {data.name} ({data.n_points} points) | models: {model_names}\n")

    fits = []
    for nm in model_names:
        mdl = M.get(nm)
        f = C.fit_model(mdl, data)
        if not args.no_evidence:
            f.lnZ, f.lnZ_err = C.evidence(mdl, data, nlive=args.nlive)
        fits.append(f)
        bf = " ".join(f"{k}={v:.3f}" for k, v in f.best.items()) or "(none)"
        print(f"  fit {nm:<20} chi2={f.chi2_min:7.2f} k={f.k}  lnZ={f.lnZ:8.3f}  best: {bf}")

    lcdm = next((f for f in fits if f.name == "lcdm"), None)
    aic0 = lcdm.aic if lcdm else 0.0
    bic0 = lcdm.bic if lcdm else 0.0
    lnz0 = lcdm.lnZ if lcdm else 0.0

    print(f"\n{'model':<20}{'k':>3}{'chi2':>9}{'AIC':>9}{'dAIC':>8}{'BIC':>9}{'dBIC':>8}"
          f"{'lnZ':>9}{'dlnZ':>8}  pref")
    print("-" * 96)
    rows = []
    for f in fits:
        dA = f.aic - aic0
        dB = f.bic - bic0
        dZ = f.lnZ - lnz0
        if f.name == "lcdm":
            pref = ""
        elif np.isfinite(dZ):
            pref = f"{jeffreys(dZ)}->{'model' if dZ > 0 else 'LCDM'}"
        else:                              # no evidence: read AIC instead
            pref = f"AIC->{'model' if dA < -2 else 'LCDM'}"
        print(f"{f.name:<20}{f.k:>3}{f.chi2_min:>9.2f}{f.aic:>9.2f}{dA:>8.2f}"
              f"{f.bic:>9.2f}{dB:>8.2f}{f.lnZ:>9.2f}{dZ:>8.2f}  {pref}")
        rows.append(dict(model=f.name, k=f.k, chi2_min=f.chi2_min, AIC=f.aic,
                         dAIC=dA, BIC=f.bic, dBIC=dB, lnZ=f.lnZ, lnZ_err=f.lnZ_err,
                         dlnZ=dZ, best=str(f.best)))

    print("\nReading: dAIC<-2 favors the model over LCDM; dBIC<=0 means it wins "
          "despite the parameter penalty; dlnZ>0 favors the model (Occam-weighted).")

    os.makedirs(args.res_dir, exist_ok=True)
    out = os.path.join(args.res_dir, f"model_comparison_{tag}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nWrote: results/model_comparison_{tag}.csv")


if __name__ == "__main__":
    main()
