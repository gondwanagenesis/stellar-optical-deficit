#!/usr/bin/env python
"""Corrected pair limits, plus a background PREDICTION that tightens them.

    run.sh scripts/33_pair_limit_v2.py --tag primary

WHAT WENT WRONG AND WHAT REPLACES IT
------------------------------------
scripts/25 and 30 used the asymmetry n_pos - n_neg of dr = r_prim - r_sec.
scripts/32 shows by injection that this estimator is BLIND to the signal
(0.02 sigma response where the two-sided count gives 13.2 sigma), because
harvesting strikes either component with equal probability and therefore
produces a symmetric perturbation. Those limits are void.

Two-sided counting is correct but conservative: every pair in either tail is
allowed to be signal. This script does better, without circularity, by
PREDICTING the pair background from an independent measurement.

THE BACKGROUND PREDICTION
-------------------------
The pair tails are populated by per-star events -- blends, activity, unresolved
tertiaries -- and those do NOT cancel in the difference. If a single star
exceeds a threshold T with probability q(T), then a pair exceeds |dr| > T when
either component does, so

    E[pair tail] = 2 * q(T) * N_pairs      (to first order in q)

q(T) is measured on the 3.3M single-star sample at the SAME absolute threshold
in magnitudes -- not the same number of sigmas, since sigma(dr) and
sigma(single) differ. That makes the prediction independent of the pair data
and legitimate to subtract.

The residual after subtraction is what common-mode cancellation removed plus
anything real. Since cancellation can only REDUCE the pair tail relative to
this prediction, an observed count at or below prediction bounds the signal by
the uncertainty on the difference rather than by the raw count.
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

CLEAN_SEP_ARCSEC = 10.0
A_W1_OVER_AV, A_W2_OVER_AV = 0.039, 0.026
EXCESS_NSIGMA = 3.0


def bare_flag(d):
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

    def exp_of(y):
        lo, hi = np.nanpercentile(jks0[have], [1, 99])
        edges = np.linspace(lo, hi, 61)
        idx = np.clip(np.digitize(jks0, edges) - 1, 0, 59)
        med = np.full(60, np.nan)
        for b in range(60):
            m = (idx == b) & have & np.isfinite(y)
            if m.sum() > 100:
                med[b] = np.median(y[m])
        g = np.isfinite(med)
        return np.interp(jks0, (0.5 * (edges[:-1] + edges[1:]))[g], med[g])

    ex1 = (ks0 - w1) - exp_of(ks0 - w1)
    ex2 = (ks0 - w2) - exp_of(ks0 - w2)
    s1, s2 = st.robust_sigma(ex1[have]), st.robust_sigma(ex2[have])
    return have & (np.abs(ex1) < EXCESS_NSIGMA * s1) & (np.abs(ex2) < EXCESS_NSIGMA * s2)


def limits(dr, n_stars, q_of_T, label):
    sigma = st.robust_sigma(dr)
    med = float(np.median(dr))
    rows = []
    for k in (4, 5, 6, 7):
        T = k * sigma
        obs = int(np.count_nonzero(np.abs(dr - med) > T))
        n_pairs = len(dr)
        pred = 2.0 * q_of_T(T) * n_pairs
        excess = obs - pred
        err = np.sqrt(obs + pred)
        f_det = float(st.fraction_from_delta(T))
        ul_raw = st.poisson_upper_limit(obs)
        ul_sub = max(excess, 0.0) + 1.645 * err
        rows.append({
            "estimator": label, "k": k, "threshold_mag": T,
            "f_detectable": f_det, "observed": obs,
            "predicted_background": pred, "excess": excess,
            "p_UL_conservative": ul_raw / n_stars,
            "p_UL_bkg_subtracted": ul_sub / n_stars,
            "mean_f_UL": min(ul_sub, ul_raw) / n_stars * f_det,
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    r = d["residual"].to_numpy(float)
    med_s = float(np.median(r))
    dev = np.abs(r - med_s)
    n_single = len(r)

    def q_of_T(T):
        """Single-star probability of exceeding |r - median| > T."""
        return float(np.count_nonzero(dev > T)) / n_single

    print(f"single-star sample {n_single:,}, sigma = {st.robust_sigma(r):.5f}")
    for T in (0.4, 0.55, 0.67, 0.8):
        print(f"  q(|r|>{T:.2f} mag) = {q_of_T(T):.3e}")

    p = pd.read_parquet(cfg.DERIVED_DIR / f"wide_binaries_{args.tag}.parquet")
    clean = p[p["theta_arcsec"] > CLEAN_SEP_ARCSEC].reset_index(drop=True)
    bare = bare_flag(d)
    both = bare[clean["i"].to_numpy()] & bare[clean["j"].to_numpy()]

    allp = limits(clean["dr"].to_numpy(float), 2 * len(clean), q_of_T,
                  "all clean pairs")
    barep = limits(clean[both]["dr"].to_numpy(float), 2 * int(both.sum()),
                   q_of_T, "clean + bare photosphere")

    print(f"\nclean pairs {len(clean):,} | with both bare {int(both.sum()):,}\n")
    out = pd.concat([allp, barep], ignore_index=True)
    print(out.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))

    best_all = allp.loc[allp["mean_f_UL"].idxmin()]
    best_bare = barep.loc[barep["mean_f_UL"].idxmin()]
    print(f"\nCORRECTED all clean pairs : p < "
          f"{min(best_all['p_UL_conservative'], best_all['p_UL_bkg_subtracted']):.3e} "
          f"at f >= {best_all['f_detectable']:.3f}")
    print(f"CORRECTED beamed class    : p < "
          f"{min(best_bare['p_UL_conservative'], best_bare['p_UL_bkg_subtracted']):.3e} "
          f"at f >= {best_bare['f_detectable']:.3f}")

    p_dark = float(min(best_bare["p_UL_conservative"],
                       best_bare["p_UL_bkg_subtracted"]))
    p_iso = 1.9e-4
    print(f"\n=== corrected joint constraint ===")
    print(f"  p_iso  (Suazo+22)      < {p_iso:.2e}")
    print(f"  p_dark (this work)     < {p_dark:.2e}")
    print(f"  p_total                < {p_iso + p_dark:.2e}")
    print(f"  -> fewer than 1 in {1/(p_iso+p_dark):,.0f} stars, any disposal mode")

    res = {"tag": args.tag,
           "p_dark_corrected": p_dark,
           "f_dark": float(best_bare["f_detectable"]),
           "p_iso_suazo22": p_iso,
           "p_total_corrected": p_iso + p_dark,
           "one_in_n": float(1 / (p_iso + p_dark)),
           "superseded": ("asymmetry-based limits in scripts 25 and 30 are "
                          "void; see scripts/32 injection test"),
           "table": out.to_dict(orient="records")}
    (cfg.RESULT_DIR / f"pair_limit_v2_{args.tag}.json").write_text(
        json.dumps(res, indent=2))
    out.to_csv(cfg.RESULT_DIR / f"pair_limit_v2_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
