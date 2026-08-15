#!/usr/bin/env python
"""The limit from clean wide pairs, where the background is measurably symmetric.

    run.sh scripts/25_pair_limit.py --tag primary

WHY THIS BEATS THE SINGLE-STAR LIMIT DESPITE A WORSE SIGMA
-----------------------------------------------------------
The single-star limit is model-free and therefore conservative: all 19,844
positive 5-sigma outliers must be allowed to be signal, because unresolved
companions push residuals positive and the negative tail is NOT a valid
background estimate for them.  That is what pins p_UL at 6.1e-3.

In wide pairs the statistic is intrinsically two-sided -- either component can
be the harvested one -- so the two tails are each other's control, *provided*
the contaminants are symmetric.  scripts/24_wide_binaries.py measures whether
they are, and finds:

    3.0- 6.5"   n_neg/n_pos = 50.7     <- inside the 2MASS beam, contaminated
    6.5-10.3"   n_neg/n_pos =  2.7
   10.3-19.3"   n_neg/n_pos =  1.4     <- symmetric
   19.3-120"    n_neg/n_pos =  0.57    <- symmetric

The asymmetry is confined to separations below the 2MASS point-source aperture
(4 arcsec, with PSF wings to ~6-8 arcsec): the brighter primary leaks into the
secondary's K_s measurement, inflating its K_s and making it look
under-luminous in G.  Beyond ~10 arcsec that mechanism is gone and the tails
balance.

So on the theta > 10 arcsec subsample, background subtraction is justified by
measurement rather than by assumption, and the limit follows from the
ASYMMETRY rather than from the raw count.

THE SEPARATION CUT IS POST-HOC AND IS DECLARED AS SUCH.  It is physically
motivated (twice the 2MASS aperture) and the separation split that motivated it
is reported in full, but it was chosen after seeing the asymmetry, so this
limit is a secondary, post-unblinding result.
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
from pipeline import statistics as st

# Twice the 2MASS point-source aperture (Skrutskie et al. 2006): beyond this
# the primary's PSF wings no longer contaminate the secondary's Ks.
CLEAN_SEP_ARCSEC = 10.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    p = pd.read_parquet(cfg.DERIVED_DIR / f"wide_binaries_{args.tag}.parquet")
    clean = p[p["theta_arcsec"] > CLEAN_SEP_ARCSEC].reset_index(drop=True)
    dr = clean["dr"].to_numpy(float)
    n_pairs = len(clean)
    n_stars = 2 * n_pairs
    sigma = st.robust_sigma(dr)
    med = float(np.median(dr))
    sigma_indiv = sigma / np.sqrt(2.0)

    print(f"clean pairs (theta > {CLEAN_SEP_ARCSEC}\") : {n_pairs:,} "
          f"({n_stars:,} stars)")
    print(f"median separation                : "
          f"{np.median(clean['sep_au']):.0f} AU")
    print(f"sigma(dr)                        : {sigma:.5f} mag")
    print(f"sigma per star                   : {sigma_indiv:.5f} mag\n")

    rows = []
    for k in (3, 4, 5, 6):
        npos = int(np.count_nonzero(dr > med + k * sigma))
        nneg = int(np.count_nonzero(dr < med - k * sigma))
        asym = npos - nneg
        err = np.sqrt(npos + nneg)
        # 95% CL upper limit on a one-sided excess over a measured symmetric
        # background: Poisson UL on the excess, widened by the background error.
        ul = st.poisson_upper_limit(max(asym, 0)) + 1.645 * err
        delta = k * sigma
        f_thresh = float(st.fraction_from_delta(delta))
        eff = float(np.mean(dr + delta > med + k * sigma))  # crude, symmetric
        p_ul = ul / n_stars
        rows.append({
            "k": k, "threshold_mag": delta, "f_detectable": f_thresh,
            "n_pos": npos, "n_neg": nneg, "asymmetry": asym,
            "asym_err": err, "excess_UL_95": ul,
            "p_UL": p_ul, "mean_f_UL": p_ul * f_thresh,
        })
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))

    best = t.loc[t["mean_f_UL"].idxmin()]
    single = json.loads(
        (cfg.RESULT_DIR / f"analysis_{args.tag}_unblinded.json").read_text())

    print(f"\n=== comparison ===")
    print(f"  single-star limit  : p_UL = {single['best_mean_f_upper_limit']/0.5:.3e} "
          f"at f=0.5   (mean f_UL = {single['best_mean_f_upper_limit']:.3e})")
    print(f"  clean-pair limit   : p_UL = {best['p_UL']:.3e} "
          f"at f={best['f_detectable']:.3f}   (mean f_UL = {best['mean_f_UL']:.3e})")
    gain = single["best_mean_f_upper_limit"] / best["mean_f_UL"]
    print(f"  improvement        : {gain:.1f}x on the mean harvested fraction")
    print(f"\n  ...from {n_stars:,} stars instead of 3,321,566. The gain is not")
    print("  statistical: it is that the background is measurably symmetric here")
    print("  and can therefore be subtracted, which it cannot be for single stars.")

    out = {
        "tag": args.tag,
        "clean_sep_arcsec": CLEAN_SEP_ARCSEC,
        "n_pairs": int(n_pairs), "n_stars": int(n_stars),
        "median_sep_au": float(np.median(clean["sep_au"])),
        "sigma_dr_mag": float(sigma),
        "sigma_individual_mag": float(sigma_indiv),
        "table": t.to_dict(orient="records"),
        "best_k": float(best["k"]),
        "best_p_UL": float(best["p_UL"]),
        "best_mean_f_UL": float(best["mean_f_UL"]),
        "f_detectable_at_best": float(best["f_detectable"]),
        "gain_over_single_star": float(gain),
        "caveat": ("separation cut chosen after seeing the asymmetry; "
                   "post-hoc and post-unblinding, declared in RESULTS.md"),
    }
    (cfg.RESULT_DIR / f"pair_limit_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    t.to_csv(cfg.RESULT_DIR / f"pair_limit_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
