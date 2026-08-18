#!/usr/bin/env python
"""Is Search A's 513-star survivor set real, or is it noise?

    run.sh scripts/44_searchA_null.py --tag primary

THE MIRROR CONTROL
------------------
An absorber can only DIM a star. So run the identical pipeline on the opposite
sign: require the fitted amplitude to be NEGATIVE, i.e. the star is anomalously
BRIGHT in G relative to K_s with a flat spectral slope. Nothing physical
produces that, so however many "survivors" the mirror yields is the pipeline's
own false-positive rate under identical cuts.

If the mirror returns a comparable number, the 513 are noise and Search A has
found nothing. If it returns far fewer, the positive set is at least asymmetric
and deserves follow-up.

This is the test that should have been written before the candidate list, not
after it.
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

SFR = [("Orion", 209.0, -19.4, 12.0), ("Ophiuchus", 353.6, 16.9, 8.0),
       ("Corona Australis", 359.7, -17.8, 6.0), ("Taurus", 172.5, -15.5, 10.0),
       ("Chamaeleon", 297.2, -15.6, 8.0), ("Lupus", 339.0, 15.0, 8.0),
       ("Perseus", 159.4, -20.0, 8.0), ("Serpens", 31.5, 5.3, 6.0),
       ("Cepheus", 110.0, 15.0, 10.0), ("Vela", 265.0, 1.5, 10.0)]


def in_sfr(l, b):
    lab = np.zeros(len(l), dtype=bool)
    for _, l0, b0, rad in SFR:
        dl = np.abs((l - l0 + 180) % 360 - 180)
        lab |= np.hypot(dl * np.cos(np.radians(b)), b - b0) < rad
    return lab


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    sl = pd.read_parquet(cfg.DERIVED_DIR / f"spectral_slope_{args.tag}.parquet")
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["source_id", "l", "b", "A_0"])
    m = sl.merge(d, on="source_id", how="left", suffixes=("", "_d"))
    m["insfr"] = in_sfr(m["l"].to_numpy(float), m["b"].to_numpy(float))
    edge = (m["alpha"] < -0.9) | (m["alpha"] > 5.9)

    common = (~m["insfr"] & (m["A_0"] < 0.15) & (np.abs(m["b"]) > 20)
              & (m["chi2_fit"] < 11.07) & ~edge
              & (m["cstar_nsigma"].abs() < 1.0) & (m["ruwe"] < 1.1)
              & (m["dchi2"] > 25) & (m["alpha_upper_2sig"] < 1.5))

    dim = common & (m["deficit_G_mag"] > 0.10) & (m["deficit_G_mag"] < 3.0)
    bright = common & (m["deficit_G_mag"] < -0.10) & (m["deficit_G_mag"] > -3.0)

    n_dim, n_bright = int(dim.sum()), int(bright.sum())
    print("identical cuts, both signs of the fitted amplitude:\n")
    print(f"  DIMMED  (physical, an absorber can do this) : {n_dim:6,}")
    print(f"  BRIGHTENED (mirror, nothing does this)      : {n_bright:6,}")
    ratio = n_dim / max(n_bright, 1)
    print(f"  ratio dim/bright                            : {ratio:6.2f}")

    a = m.loc[common, "deficit_G_mag"].to_numpy(float)
    a = a[np.isfinite(a)]
    print(f"\n  amplitude median {np.median(a):+.4f} mag, frac>0 = "
          f"{np.mean(a > 0):.3f}  (n={len(a):,})")

    if ratio < 1.5:
        verdict = ("NOISE. The mirror yields a comparable count, so the survivor "
                   "set is the pipeline's own false-positive rate and Search A "
                   "has found nothing.")
    elif ratio < 3:
        verdict = ("MARGINAL. Mild asymmetry only; the set is mostly false "
                   "positives and no individual object is compelling.")
    else:
        verdict = (f"ASYMMETRIC at {ratio:.1f}:1. The dimmed set exceeds the "
                   f"unphysical mirror, so it is not purely noise.")
    print(f"\nVERDICT: {verdict}")

    out = {"n_dimmed": n_dim, "n_brightened": n_bright, "ratio": float(ratio),
           "frac_amplitude_positive": float(np.mean(a > 0)), "verdict": verdict}
    (cfg.RESULT_DIR / f"searchA_null_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
