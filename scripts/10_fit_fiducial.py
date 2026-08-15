#!/usr/bin/env python
"""Step 2: cross-validate and fit the fiducial relation; report intrinsic scatter.

    run.sh scripts/10_fit_fiducial.py --tag primary --nir-control

Writes results/fiducial_<tag>.json, results/cv_<tag>.csv and
data/cache/derived/<tag>_resid.parquet.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import fiducial as fid
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("fit")


def build_covs(d: pd.DataFrame, mh_degree: int, nir: bool):
    covs = [(d["mh_gspphot"].to_numpy(dtype=float), mh_degree)]
    if nir:
        covs.append((d["j_ks0"].to_numpy(dtype=float), 1))
    return covs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--knots", type=int, default=None)
    ap.add_argument("--mh-degree", type=int, default=None)
    ap.add_argument("--cv-max-n", type=int, default=400_000)
    ap.add_argument("--nir-control", action="store_true",
                    help="add the dereddened (J-Ks) colour as a control; it is "
                         "immune to optical harvesting (see fiducial.py)")
    ap.add_argument("--out-tag", default=None)
    args = ap.parse_args()
    out_tag = args.out_tag or args.tag

    df = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}.parquet")
    log.info("loaded %d stars", len(df))

    # Dereddened NIR colour. A_J and A_Ks come from the same law already used
    # for A_G, so the reddening correction is internally consistent.
    from pipeline import extinction as ext
    a0 = df["A_0"].to_numpy(dtype=float)
    bp_rp = df["bp_rp"].to_numpy(dtype=float)
    a_j = ext.deredden("J", a0, bp_rp, law="fitz19")
    a_ks = df["A_Ks"].to_numpy(dtype=float)
    df["j_ks0"] = (df["tmass_j_m"].to_numpy(dtype=float) - a_j) - \
                  (df["tmass_ks_m"].to_numpy(dtype=float) - a_ks)

    need = df["mh_gspphot"].notna() & np.isfinite(df["j_ks0"])
    log.info("GSP-Phot [M/H] present: %.1f%%; finite J-Ks: %.1f%%; both: %.1f%%",
             100 * df["mh_gspphot"].notna().mean(),
             100 * np.isfinite(df["j_ks0"]).mean(), 100 * need.mean())

    d = df[need].reset_index(drop=True)
    m_ks = d["M_Ks"].to_numpy(dtype=float)
    m_g = d["M_G"].to_numpy(dtype=float)
    mh = d["mh_gspphot"].to_numpy(dtype=float)
    jks = d["j_ks0"].to_numpy(dtype=float)

    if args.knots is not None and args.mh_degree is not None:
        best_k, best_p = args.knots, args.mh_degree
    else:
        cv = fid.cross_validate(m_ks, mh, m_g, max_fit_n=args.cv_max_n,
                                extra_covs=[jks] if args.nir_control else None)
        best_k, best_p = cv.best_knots, cv.best_mh_degree
        cv.to_frame().to_csv(cfg.RESULT_DIR / f"cv_{out_tag}.csv", index=False)

    log.info("fitting with knots=%d mh_degree=%d nir_control=%s",
             best_k, best_p, args.nir_control)
    covs = build_covs(d, best_p, args.nir_control)
    fit = fid.fit_fiducial(m_ks, covs, m_g, best_k)
    resid = fit.residuals(m_ks, covs, m_g)

    sig_meas = fid.residual_uncertainty(d, fit, covs)
    s_obs, s_int = fid.intrinsic_scatter(resid, sig_meas)
    slopes = fid.slope(fit, m_ks, covs)

    # Variant A: no metallicity, no NIR -- pure M_G(M_Ks).
    fit_bare = fid.fit_fiducial(m_ks, None, m_g, best_k)
    resid_bare = fit_bare.residuals(m_ks, None, m_g)
    s_obs_bare, _ = fid.intrinsic_scatter(
        resid_bare, fid.residual_uncertainty(d, fit_bare, None))

    shape = st.shape_stats(resid)

    out = {
        "tag": out_tag, "source_tag": args.tag,
        "nir_control": bool(args.nir_control),
        "n_stars_total": int(len(df)), "n_stars_fitted": int(len(d)),
        "n_interior_knots": int(best_k), "mh_degree": int(best_p),
        "n_params": int(fit.n_params), "converged": bool(fit.converged),
        "sigma_observed_mag": float(s_obs),
        "sigma_measurement_mag": float(np.sqrt(np.nanmean(sig_meas ** 2))),
        "sigma_intrinsic_mag": float(s_int),
        "sigma_observed_bare_mag": float(s_obs_bare),
        "slope_dMG_dMKs_median": float(np.nanmedian(slopes)),
        "slope_dMG_dMKs_p16": float(np.nanpercentile(slopes, 16)),
        "slope_dMG_dMKs_p84": float(np.nanpercentile(slopes, 84)),
        "mean_residual": float(shape.mean),
        "median_residual": float(shape.median),
        "bowley_skew": float(shape.skew),
        "naive_sigma_over_sqrtN": float(s_obs / np.sqrt(len(d))),
    }
    (cfg.RESULT_DIR / f"fiducial_{out_tag}.json").write_text(json.dumps(out, indent=2))

    d = d.assign(residual=resid, sigma_meas=sig_meas, slope_local=slopes,
                 residual_bare=resid_bare)
    d.to_parquet(cfg.DERIVED_DIR / f"{out_tag}_resid.parquet",
                 index=False, compression="zstd")

    print("\n=== fiducial fit ===")
    for k, v in out.items():
        print(f"  {k:30s} {v}")
    print("\n=== residual shape ===")
    print(shape.summary())
    print(shape.table.to_string(index=False))
    print("\nNOTE: the ROBUST centre of the residuals is ~0 by construction "
          "(free spline\n      intercept), so neither the mean residual nor a "
          "uniform harvesting\n      fraction is measurable here. "
          "sigma/sqrt(N) is NOT a sensitivity: see step 3.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
