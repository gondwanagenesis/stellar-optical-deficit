#!/usr/bin/env python
"""What KIND of deficit is each outlier? Dust, blend, or neither?

    run.sh scripts/26_deficit_colour.py --tag primary

A positive residual (fainter in G at fixed M_Ks) has three possible causes, and
they make DIFFERENT predictions for the optical colour:

  (a) unmodelled dust, alpha ~ 2
        reddens the star.  Predicted excess follows the Fitz19 reddening
        vector, and it is large: dC ~ 0.65 * r, so r = 0.5 mag implies a
        0.32 mag colour excess.  Unmissable.

  (b) 2MASS blending (a neighbour inflating K_s)
        the OPTICAL photometry is untouched, so dC = 0 while r > 0.

  (c) an absorber that is neither
        dC somewhere else entirely, including negative.

So the ratio  dC_observed / dC_predicted-if-dust  separates them:
      ~1  -> dust
      ~0  -> blending (our known dominant contaminant)
    else  -> anomalous, and worth a second look

TEMPERATURE ANCHOR
------------------
The expected intrinsic colour is taken from (J-H)_0, a near-infrared colour.
This matters: (J-H) is essentially untouched by an optical absorber, so using
it to predict the expected optical colour does not launder the signal away.
Using M_Ks would be circular here because M_Ks is exactly what blending
corrupts.

WHY THIS IS THE RIGHT TEST FOR A HIDER
--------------------------------------
An engineered absorber optimised to evade the deficit test would sit at the
measured blind spot, alpha ~ 0.19.  Nothing natural sits there: interstellar
dust is alpha ~ 2, grey occultation is alpha = 0, and free-free/electron
scattering is alpha ~ 0.  A population piling up near alpha ~ 0.2 would be the
single most suspicious thing this dataset could contain, precisely because it
is the value that makes the star invisible to the primary statistic.
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

# Leverage of a dust-like absorber into the residual, from anchor.leverage(2, s)
DUST_LEVERAGE = 0.892


def running_median(x, y, xbins):
    idx = np.digitize(x, xbins)
    med = np.full(len(xbins) + 1, np.nan)
    for b in range(len(xbins) + 1):
        m = idx == b
        if m.sum() > 50:
            med[b] = np.median(y[m])
    return med[idx]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["source_id", "ra", "dec", "residual", "A_0",
                                 "bp_rp", "tmass_j_m", "tmass_h_m",
                                 "tmass_ks_m", "phot_g_mean_mag",
                                 "cstar_nsigma", "ruwe", "j_ks0"])
    r = d["residual"].to_numpy(float)
    sigma = st.robust_sigma(r)
    med_r = float(np.median(r))

    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_j = ext.deredden("J", a0, bp_rp)
    a_h = ext.deredden("H", a0, bp_rp)
    jh0 = ((d["tmass_j_m"].to_numpy(float) - a_j)
           - (d["tmass_h_m"].to_numpy(float) - a_h))
    bprp0 = ext.intrinsic_bp_rp(a0, bp_rp)

    ok = np.isfinite(jh0) & np.isfinite(bprp0) & (jh0 > 0.0) & (jh0 < 1.0)
    print(f"{ok.sum():,} of {len(d):,} stars have usable (J-H)_0")

    bins = np.linspace(0.0, 1.0, 60)
    expected = np.full(len(d), np.nan)
    expected[ok] = running_median(jh0[ok], bprp0[ok], bins)
    dC_obs = bprp0 - expected

    # Colour excess predicted if the whole residual were unmodelled dust
    kbp = ext.k_fitz19("BP", bprp0, a0)
    krp = ext.k_fitz19("RP", bprp0, a0)
    kg = ext.k_fitz19("G", bprp0, a0)
    dC_dust = (r - med_r) * (kbp - krp) / (kg * DUST_LEVERAGE)

    ratio = np.where(np.abs(dC_dust) > 1e-3, dC_obs / dC_dust, np.nan)

    pos = ok & (r > med_r + args.k * sigma)
    neg = ok & (r < med_r - args.k * sigma)
    bulk = ok & (np.abs(r - med_r) < sigma)
    print(f"\npositive {args.k}-sigma outliers with usable colours: {pos.sum():,}")
    print(f"negative {args.k}-sigma outliers with usable colours: {neg.sum():,}")

    def describe(mask, label):
        v = ratio[mask]
        v = v[np.isfinite(v)]
        if len(v) < 10:
            return None
        q = np.percentile(v, [10, 25, 50, 75, 90])
        frac_dust = float(np.mean((v > 0.6) & (v < 1.4)))
        frac_grey = float(np.mean(np.abs(v) < 0.25))
        frac_other = 1.0 - frac_dust - frac_grey
        print(f"\n  {label}  (n={len(v):,})")
        print(f"    dC_obs/dC_dust percentiles 10/25/50/75/90: "
              f"{q[0]:+.2f} {q[1]:+.2f} {q[2]:+.2f} {q[3]:+.2f} {q[4]:+.2f}")
        print(f"    consistent with DUST   (0.6-1.4) : {100*frac_dust:5.1f}%")
        print(f"    consistent with BLEND  (|x|<0.25): {100*frac_grey:5.1f}%")
        print(f"    NEITHER                          : {100*frac_other:5.1f}%")
        return {"label": label, "n": int(len(v)), "median": float(q[2]),
                "frac_dust": frac_dust, "frac_blend": frac_grey,
                "frac_other": float(frac_other)}

    print("\n=== what are the deficits made of? ===")
    res = [describe(pos, "positive outliers (deficit-like)"),
           describe(neg, "negative outliers (over-luminous control)")]

    # Split the positive outliers by blending proxy: if the 'blend' class is
    # real, it must concentrate at high C*.
    hi_c = pos & (d["cstar_nsigma"].to_numpy(float)
                  > np.nanquantile(d["cstar_nsigma"], 0.75))
    lo_c = pos & (d["cstar_nsigma"].to_numpy(float)
                  < np.nanquantile(d["cstar_nsigma"], 0.25))
    res.append(describe(hi_c, "positive outliers, HIGH C* (blend-prone)"))
    res.append(describe(lo_c, "positive outliers, LOW C* (clean photometry)"))

    out = {"tag": args.tag, "k_sigma": args.k, "sigma_mag": float(sigma),
           "dust_leverage": DUST_LEVERAGE,
           "classes": [x for x in res if x]}
    (cfg.RESULT_DIR / f"deficit_colour_{args.tag}.json").write_text(
        json.dumps(out, indent=2))

    # The anomalous set: strong deficit, clean photometry, colour excess that
    # is neither dust nor zero.
    anom = pos & np.isfinite(ratio) & (ratio > -0.5) & (ratio < 0.4) & \
        (d["cstar_nsigma"].to_numpy(float) < 0)
    print(f"\n=== anomalous set ===")
    print(f"  strong deficit + clean C* + non-dust, non-zero colour : "
          f"{int(anom.sum()):,}")
    if anom.sum():
        sub = d[anom].copy()
        sub["ratio"] = ratio[anom]
        sub["implied_f"] = st.fraction_from_delta(r[anom] - med_r)
        sub[["source_id", "ra", "dec", "residual", "ratio", "implied_f",
             "cstar_nsigma", "ruwe", "A_0"]].nlargest(15, "residual").to_csv(
            cfg.RESULT_DIR / f"anomalous_colour_{args.tag}.csv", index=False)
        print(f"  written to results/anomalous_colour_{args.tag}.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
