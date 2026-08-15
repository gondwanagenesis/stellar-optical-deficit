#!/usr/bin/env python
"""The headline number this work can claim: a limit on the BEAMED class.

    run.sh scripts/30_beamed_class_limit.py --tag primary

Combines the two things this analysis does better than the single-star test:

  * clean wide pairs (theta > 10"), where the background is MEASURABLY
    symmetric and can therefore be subtracted (paper Sec 5.6);
  * a mid-infrared veto requiring both components to have a MEASURED bare
    photosphere in W1 and W2 (paper Sec 29 / this script's sibling).

The result is a limit on stars that intercept a large optical fraction AND
show no warm re-emission -- the beaming-consistent, cold, or
non-thermally-exporting class.  Suazo et al. (2022) model
L = (1-gamma)L_star + gamma*BB(T_DS) and so constrain the isotropic warm case
far better than we do; they are not designed to constrain this one.

This is the only number in the paper that is not superseded by prior work.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import extinction as ext
from pipeline import statistics as st

A_W1_OVER_AV, A_W2_OVER_AV = 0.039, 0.026   # Wang & Chen 2019, Table 3
EXCESS_NSIGMA = 3.0
CLEAN_SEP_ARCSEC = 10.0


def bare_photosphere_flag(d: pd.DataFrame) -> np.ndarray:
    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    ks0 = d["tmass_ks_m"].to_numpy(float) - ext.deredden("Ks", a0, bp_rp)
    j0 = d["tmass_j_m"].to_numpy(float) - ext.deredden("J", a0, bp_rp)
    jks0 = j0 - ks0
    w1 = d["wise_w1mpro"].to_numpy(float) - A_W1_OVER_AV * a0
    w2 = d["wise_w2mpro"].to_numpy(float) - A_W2_OVER_AV * a0
    e1 = d["wise_w1mpro_error"].to_numpy(float)
    e2 = d["wise_w2mpro_error"].to_numpy(float)
    have = (np.isfinite(w1) & np.isfinite(w2) & np.isfinite(e1) & np.isfinite(e2)
            & (e1 < 0.2) & (e2 < 0.2) & np.isfinite(jks0))

    def expected(x, y):
        edges = np.linspace(np.nanpercentile(x[have], 1),
                            np.nanpercentile(x[have], 99), 61)
        idx = np.clip(np.digitize(x, edges) - 1, 0, 59)
        med = np.full(60, np.nan)
        for b in range(60):
            m = (idx == b) & have & np.isfinite(y)
            if m.sum() > 100:
                med[b] = np.median(y[m])
        g = np.isfinite(med)
        c = 0.5 * (edges[:-1] + edges[1:])
        return np.interp(x, c[g], med[g])

    ex1 = (ks0 - w1) - expected(jks0, ks0 - w1)
    ex2 = (ks0 - w2) - expected(jks0, ks0 - w2)
    s1 = st.robust_sigma(ex1[have]); s2 = st.robust_sigma(ex2[have])
    return have & (np.abs(ex1) < EXCESS_NSIGMA * s1) & (np.abs(ex2) < EXCESS_NSIGMA * s2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    bare = bare_photosphere_flag(d)
    print(f"stars with MEASURED bare photosphere: {bare.sum():,} / {len(d):,} "
          f"({100*bare.mean():.1f}%)")

    p = pd.read_parquet(cfg.DERIVED_DIR / f"wide_binaries_{args.tag}.parquet")
    clean = p[p["theta_arcsec"] > CLEAN_SEP_ARCSEC].reset_index(drop=True)
    both = bare[clean["i"].to_numpy()] & bare[clean["j"].to_numpy()]
    sub = clean[both].reset_index(drop=True)
    print(f"clean pairs (theta>{CLEAN_SEP_ARCSEC:.0f}\")            : {len(clean):,}")
    print(f"  ...with BOTH components bare        : {len(sub):,} "
          f"({len(sub)*2:,} stars)")

    dr = sub["dr"].to_numpy(float)
    n_pairs, n_stars = len(sub), 2 * len(sub)
    sigma = st.robust_sigma(dr)
    med = float(np.median(dr))
    print(f"  sigma(dr)                           : {sigma:.5f} mag")

    rows = []
    for k in (3, 4, 5, 6):
        npos = int(np.count_nonzero(dr > med + k * sigma))
        nneg = int(np.count_nonzero(dr < med - k * sigma))
        asym, err = npos - nneg, np.sqrt(npos + nneg)
        ul = st.poisson_upper_limit(max(asym, 0)) + 1.645 * err
        delta = k * sigma
        f_det = float(st.fraction_from_delta(delta))
        rows.append({"k": k, "threshold_mag": delta, "f_detectable": f_det,
                     "n_pos": npos, "n_neg": nneg, "asymmetry": asym,
                     "asym_err": err, "excess_UL_95": ul,
                     "p_UL": ul / n_stars, "mean_f_UL": ul / n_stars * f_det})
    t = pd.DataFrame(rows)
    print("\n=== limit on the beaming-consistent class ===")
    print(t.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))

    best = t.loc[t["mean_f_UL"].idxmin()]
    print(f"\nBEST: p < {best['p_UL']:.3e} at f >= {best['f_detectable']:.3f}")
    print(f"      mean harvested fraction f_bar < {best['mean_f_UL']:.3e}")
    print("\n  This applies to stars with a large optical deficit AND a")
    print("  MEASURED absence of mid-infrared excess in W1 and W2 -- the class")
    print("  an IR-excess estimator is not designed to constrain.")
    print(f"\n  For scale: Suazo+22 give p < 1.9e-4 at gamma >= 0.5, but for")
    print("  spheres re-radiating isotropically at ~300 K.")

    out = {"tag": args.tag, "n_pairs": int(n_pairs), "n_stars": int(n_stars),
           "sigma_dr": float(sigma),
           "frac_bare": float(bare.mean()),
           "best_p_UL": float(best["p_UL"]),
           "best_f": float(best["f_detectable"]),
           "best_mean_f_UL": float(best["mean_f_UL"]),
           "table": t.to_dict(orient="records"),
           "claim": ("limit on optical deficit WITH measured absence of mid-IR "
                     "excess; complements rather than competes with IR-excess "
                     "searches")}
    (cfg.RESULT_DIR / f"beamed_class_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    t.to_csv(cfg.RESULT_DIR / f"beamed_class_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
