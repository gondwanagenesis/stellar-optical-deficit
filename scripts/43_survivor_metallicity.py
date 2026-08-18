#!/usr/bin/env python
"""Are the Search-A survivors metal-poor halo subdwarfs?

    run.sh scripts/43_survivor_metallicity.py --tag primary

Search A (scripts/40) predicted each band from (M_Ks, (J-H)_0) and deliberately
used NO metallicity control, because the point was to measure a spectral slope
rather than to minimise scatter. That leaves an obvious loophole.

Metal-poor subdwarfs are underluminous in the optical at fixed M_Ks: lower
metal opacity shifts the whole SED, and the shift is smooth across wavelength,
so a power-law absorber will happily fit it and return a slope unlike dust.
The survivors sit at high Galactic latitude, which is where the metal-poor halo
lives. So the prediction is sharp: if this is the explanation, the survivors
must be systematically metal-poor.
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    s = pd.read_csv(cfg.RESULT_DIR / f"searchA_candidates_{args.tag}.csv")
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["source_id", "mh_gspphot", "b", "dist_pc",
                                 "M_Ks", "pmra", "pmdec", "parallax_corr"])
    m = s.merge(d, on="source_id", how="left", suffixes=("", "_d"))

    mh_s = m["mh_gspphot"].to_numpy(float)
    mh_all = d["mh_gspphot"].to_numpy(float)
    hi_b = np.abs(d["b"].to_numpy(float)) > 20

    print("=== metallicity ===")
    for lab, v in (("Search-A survivors", mh_s),
                   ("all stars, |b|>20", mh_all[hi_b]),
                   ("all stars", mh_all)):
        v = v[np.isfinite(v)]
        if not len(v):
            print(f"  {lab:24s} no GSP-Phot [M/H] available")
            continue
        print(f"  {lab:24s} n={len(v):8,}  median [M/H] = {np.median(v):+.3f}"
              f"   frac < -0.5 = {np.mean(v < -0.5):.3f}"
              f"   frac < -1.0 = {np.mean(v < -1.0):.3f}")

    # tangential speed: halo stars are fast
    k = 4.74047 / np.maximum(d["parallax_corr"].to_numpy(float), 1e-6)
    v_all = k * np.hypot(d["pmra"].to_numpy(float), d["pmdec"].to_numpy(float))
    k_s = 4.74047 / np.maximum(m["parallax_corr"].to_numpy(float), 1e-6)
    v_s = k_s * np.hypot(m["pmra"].to_numpy(float), m["pmdec"].to_numpy(float))
    print("\n=== kinematics (halo stars are FAST) ===")
    print(f"  survivors      median |v_tan| = {np.nanmedian(v_s):6.1f} km/s   "
          f"frac >100 km/s = {np.nanmean(v_s > 100):.3f}")
    print(f"  all |b|>20     median |v_tan| = {np.nanmedian(v_all[hi_b]):6.1f} km/s   "
          f"frac >100 km/s = {np.nanmean(v_all[hi_b] > 100):.3f}")

    frac_poor_s = float(np.nanmean(mh_s < -0.5))
    frac_poor_all = float(np.nanmean(mh_all[hi_b] < -0.5))
    fast_s = float(np.nanmean(v_s > 100))
    fast_all = float(np.nanmean(v_all[hi_b] > 100))

    verdict = ("SURVIVORS ARE METAL-POOR / HALO -- Search A's lack of a "
               "metallicity control explains them"
               if frac_poor_s > 1.5 * frac_poor_all or fast_s > 1.5 * fast_all
               else "survivors are NOT preferentially metal-poor or fast; "
                    "the subdwarf explanation does not account for them")
    print(f"\nVERDICT: {verdict}")

    out = {"n_survivors": int(len(m)),
           "median_mh_survivors": float(np.nanmedian(mh_s)),
           "median_mh_highb": float(np.nanmedian(mh_all[hi_b])),
           "frac_metalpoor_survivors": frac_poor_s,
           "frac_metalpoor_highb": frac_poor_all,
           "median_vtan_survivors": float(np.nanmedian(v_s)),
           "median_vtan_highb": float(np.nanmedian(v_all[hi_b])),
           "frac_fast_survivors": fast_s, "frac_fast_highb": fast_all,
           "verdict": verdict}
    (cfg.RESULT_DIR / f"survivor_metallicity_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
