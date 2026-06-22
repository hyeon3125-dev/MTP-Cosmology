"""
run_realfit.py — fit windowed-IDE to real DESI DR1 BAO + LambdaCDM evidence (P2).

Computes the Bayesian evidence of the windowed-IDE model (free coupling) and of
LambdaCDM (beta0=0, zero free parameters) on the real DESI DR1 BAO data, and
reports Delta ln Z (model - LambdaCDM) with a Jeffreys-scale reading.

Usage:
    python scripts/run_realfit.py [--fit beta0 | beta0,z_star,sigma] [--nlive N]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from mtp_cosmology.realfit import (  # noqa: E402
    loglike_model, loglike_lcdm, prior_transform, ALL_NAMES,
)
from mtp_cosmology.data import desi_dr1  # noqa: E402

import dynesty  # noqa: E402


def jeffreys(dlnz: float) -> str:
    a = abs(dlnz)
    fav = "windowed-IDE" if dlnz > 0 else "LambdaCDM"
    if a < 1.0:
        return "inconclusive (|Delta lnZ| < 1)"
    if a < 2.5:
        return f"weak preference for {fav}"
    if a < 5.0:
        return f"moderate preference for {fav}"
    return f"strong preference for {fav}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default="beta0",
                    help="subset of beta0,z_star,sigma,z_c,dz to vary")
    ap.add_argument("--nlive", type=int, default=400)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--res-dir", default=os.path.join(os.path.dirname(__file__), "..", "results"))
    args = ap.parse_args()

    fit = [s.strip() for s in args.fit.split(",") if s.strip()]
    assert all(f in ALL_NAMES for f in fit), f"fit names must be in {ALL_NAMES}"
    ndim = len(fit)

    aniso, iso = desi_dr1()
    ndata = 2 * len(aniso) + len(iso)
    print(f"DESI DR1 BAO: {len(aniso)} anisotropic (DM/rd,DH/rd) + {len(iso)} isotropic (DV/rd)"
          f" = {ndata} data points")
    print(f"Fitting: {fit}  (ndim={ndim})\n")

    # LambdaCDM: zero free parameters -> ln Z = log-likelihood at the fixed point.
    lnz_lcdm = loglike_lcdm()
    chi2_lcdm = -2.0 * lnz_lcdm
    print(f"LambdaCDM (fixed Planck): chi^2 = {chi2_lcdm:.2f} / {ndata} dof, "
          f"ln Z = {lnz_lcdm:.3f}")

    # Windowed-IDE: nested sampling for evidence + posterior.
    rng = np.random.default_rng(args.seed)
    sampler = dynesty.NestedSampler(
        loglike_model, prior_transform, ndim,
        logl_args=(fit,), ptform_args=(fit,),
        nlive=args.nlive, rstate=rng,
    )
    sampler.run_nested(print_progress=False)
    res = sampler.results
    lnz_model = float(res.logz[-1])
    lnz_model_err = float(res.logzerr[-1])

    # Posterior summary (importance-weighted).
    samples = res.samples
    weights = np.exp(res.logwt - res.logz[-1])
    mean = np.average(samples, axis=0, weights=weights)
    var = np.average((samples - mean) ** 2, axis=0, weights=weights)
    std = np.sqrt(var)
    # best-fit chi^2 at posterior mean
    chi2_model = -2.0 * loglike_model(mean, fit)

    print(f"windowed-IDE: ln Z = {lnz_model:.3f} +/- {lnz_model_err:.3f}, "
          f"chi^2(mean) = {chi2_model:.2f}\n")
    print(f"{'param':<8}{'mean':>11}{'std':>10}{'95% upper':>12}")
    print("-" * 42)
    rows = []
    for i, nm in enumerate(fit):
        # 95% one-sided upper limit from weighted samples
        order = np.argsort(samples[:, i])
        cw = np.cumsum(weights[order])
        cw /= cw[-1]
        upper95 = float(samples[order][np.searchsorted(cw, 0.95), i])
        print(f"{nm:<8}{mean[i]:>11.4f}{std[i]:>10.4f}{upper95:>12.4f}")
        rows.append({"param": nm, "mean": mean[i], "std": std[i], "upper95": upper95})

    dlnz = lnz_model - lnz_lcdm
    print(f"\nDelta ln Z (windowed-IDE - LambdaCDM) = {dlnz:+.3f} +/- {lnz_model_err:.3f}")
    print(f"  Jeffreys scale: {jeffreys(dlnz)}")
    print(f"  (positive favors the coupling; Occam-penalized for {ndim} extra param(s))")

    os.makedirs(args.res_dir, exist_ok=True)
    tag = "-".join(fit)
    summary = pd.DataFrame(rows)
    summary["lnZ_model"] = lnz_model
    summary["lnZ_model_err"] = lnz_model_err
    summary["lnZ_lcdm"] = lnz_lcdm
    summary["delta_lnZ"] = dlnz
    summary["chi2_model"] = chi2_model
    summary["chi2_lcdm"] = chi2_lcdm
    summary["ndata"] = ndata
    out = os.path.join(args.res_dir, f"realfit_desi_{tag}.csv")
    summary.to_csv(out, index=False)
    print(f"\nWrote: results/realfit_desi_{tag}.csv")


if __name__ == "__main__":
    main()
