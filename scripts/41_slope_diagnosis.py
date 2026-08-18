#!/usr/bin/env python
"""Are the low-alpha absorbers anomalous, or is it grain growth?

    run.sh scripts/41_slope_diagnosis.py --tag primary

scripts/40 finds 1,994 stars whose absorber slope is significantly below the
diffuse-ISM value, with clean photometry and astrometry. Before treating any of
them as candidates, note where the strongest ones sit:

    RA 79-89, Dec -8..+1   Orion
    RA 244-249, Dec -19..-24  Ophiuchus
    RA 285-289, Dec -36..-37  Corona Australis

Those are star-forming regions, and there is a textbook explanation. Grain
growth in dense molecular cores FLATTENS the extinction law: R_V rises from
~3.1 in the diffuse ISM to 5-6 in dense cores, which is exactly a LOWER alpha.
So low alpha is the expected signature of dense-cloud dust, not of engineering.

The discriminating prediction: if it is grain growth, low alpha must track
extinction column and cloud membership. If it is not, low-alpha stars should be
distributed like the rest of the sample.
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

# Well-known nearby star-forming complexes, (l, b, radius_deg)
SFR = [("Orion", 209.0, -19.4, 12.0), ("Ophiuchus", 353.6, 16.9, 8.0),
       ("Corona Australis", 359.7, -17.8, 6.0), ("Taurus", 172.5, -15.5, 10.0),
       ("Chamaeleon", 297.2, -15.6, 8.0), ("Lupus", 339.0, 15.0, 8.0),
       ("Perseus", 159.4, -20.0, 8.0), ("Serpens", 31.5, 5.3, 6.0),
       ("Cepheus", 110.0, 15.0, 10.0), ("Vela", 265.0, 1.5, 10.0)]


def in_sfr(l, b):
    lab = np.array([""] * len(l), dtype=object)
    for name, l0, b0, rad in SFR:
        dl = np.abs((l - l0 + 180) % 360 - 180)
        d = np.hypot(dl * np.cos(np.radians(b)), b - b0)
        m = (d < rad) & (lab == "")
        lab[m] = name
    return lab


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    sl = pd.read_parquet(cfg.DERIVED_DIR / f"spectral_slope_{args.tag}.parquet")
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["source_id", "l", "b", "A_0"])
    m = sl.merge(d, on="source_id", how="left", suffixes=("", "_d"))

    strong = m[(m["deficit_G_mag"] > 0.10) & (m["deficit_G_mag"] < 3.0)
               & (m["dchi2"] > 25)].copy()
    strong["sfr"] = in_sfr(strong["l"].to_numpy(float),
                           strong["b"].to_numpy(float))
    low = strong["alpha_upper_2sig"] < 1.5
    dusty = strong["alpha"].between(1.8, 2.2)

    print(f"significant absorbers: {len(strong):,}")
    print(f"  alpha significantly < 1.5 : {int(low.sum()):,}")
    print(f"  alpha consistent with dust: {int(dusty.sum()):,}\n")

    print("=== TEST 1: does low alpha track extinction column? ===")
    for lab, sub in (("alpha < 1.5", strong[low]),
                     ("alpha ~ 2 (dust)", strong[dusty]),
                     ("all significant", strong)):
        a0 = sub["A_0"].to_numpy(float)
        print(f"  {lab:20s} n={len(sub):7,}  median A_0 = {np.nanmedian(a0):.4f}"
              f"   90th pct = {np.nanpercentile(a0, 90):.4f}"
              f"   frac A_0>0.3 = {np.mean(a0 > 0.3):.3f}")

    print("\n=== TEST 2: are they in known star-forming complexes? ===")
    for lab, sub in (("alpha < 1.5", strong[low]),
                     ("alpha ~ 2 (dust)", strong[dusty])):
        f_sfr = float((sub["sfr"] != "").mean())
        print(f"  {lab:20s} fraction in a known SFR = {f_sfr:.3f}")
    print("\n  breakdown for alpha < 1.5:")
    vc = strong[low]["sfr"].value_counts()
    for k, v in vc.items():
        nm = k if k else "(not in a listed SFR)"
        print(f"    {nm:24s} {v:6,}  ({100*v/max(int(low.sum()),1):5.1f}%)")

    print("\n=== TEST 3: Galactic latitude ===")
    for lab, sub in (("alpha < 1.5", strong[low]),
                     ("alpha ~ 2 (dust)", strong[dusty]),
                     ("whole sample", m)):
        b = np.abs(sub["b"].to_numpy(float))
        print(f"  {lab:20s} median |b| = {np.nanmedian(b):5.1f} deg   "
              f"frac |b|<15 = {np.mean(b < 15):.3f}")

    # The residual set that survives every mundane explanation. The fit-quality
    # and grid-edge cuts are essential and were missing on the first pass: a
    # star with chi2 = 300 on 5 degrees of freedom is not an absorber, it is bad
    # photometry, and alpha pinned at the grid edge (-1 or 6) means the fit
    # never converged inside the range at all.
    edge = (strong["alpha"] < -0.9) | (strong["alpha"] > 5.9)
    survivors = strong[low & (strong["sfr"] == "")
                       & (strong["A_0"] < 0.15)
                       & (np.abs(strong["b"]) > 20)
                       & (strong["chi2_fit"] < 11.07)   # p>0.05 for 5 dof
                       & ~edge]
    print(f"\n=== SURVIVORS: low alpha, NOT in an SFR, low extinction, "
          f"|b|>20 ===")
    print(f"  {len(survivors):,} stars")
    if len(survivors):
        cols = ["source_id", "ra", "dec", "l", "b", "alpha",
                "alpha_upper_2sig", "deficit_G_mag", "A_0", "chi2_fit"]
        print(survivors.nlargest(25, "deficit_G_mag")[cols].to_string(
            index=False, float_format=lambda v: f"{v:10.4g}"))
        survivors.to_csv(cfg.RESULT_DIR / f"slope_survivors_{args.tag}.csv",
                         index=False)

    frac_sfr_low = float((strong[low]["sfr"] != "").mean())
    frac_sfr_dust = float((strong[dusty]["sfr"] != "").mean())
    # Judge on Galactic latitude, not on the crude SFR circles: those cover only
    # a handful of named complexes and miss most of the dusty plane.
    b_low = float(np.nanmedian(np.abs(strong[low]["b"])))
    b_dust = float(np.nanmedian(np.abs(strong[dusty]["b"])))
    verdict = ("low alpha is DUST-ASSOCIATED -- it concentrates strongly toward "
               "the Galactic plane, consistent with grain growth flattening the "
               "extinction law in dense material"
               if b_low < 0.7 * b_dust else
               "low alpha does NOT concentrate toward the plane -- not "
               "obviously dust")
    print(f"\nVERDICT: {verdict}")

    out = {"n_significant": int(len(strong)), "n_low_alpha": int(low.sum()),
           "median_A0_low_alpha": float(np.nanmedian(strong[low]["A_0"])),
           "median_A0_dustlike": float(np.nanmedian(strong[dusty]["A_0"])),
           "frac_in_sfr_low_alpha": frac_sfr_low,
           "frac_in_sfr_dustlike": frac_sfr_dust,
           "n_survivors": int(len(survivors)), "verdict": verdict}
    (cfg.RESULT_DIR / f"slope_diagnosis_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
