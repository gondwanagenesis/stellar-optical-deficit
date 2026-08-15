#!/usr/bin/env python
"""Itemise the systematic budget from the null-test outputs.

    run.sh scripts/22_systematic_budget.py --tag primary

Each line is a measured or propagated contribution to the mean-residual
uncertainty, in magnitudes.  The total is quoted two ways: in quadrature
(appropriate if the terms were independent, which several are not) and as the
largest single term (the honest floor, since correlated terms cannot be assumed
to cancel).
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

# Lindegren et al. 2021, A&A 649, A4 quote residual spatial/magnitude structure
# in the parallax zero point at the ~10 uas level after their correction.
ZP_RESIDUAL_UAS = 10.0


def jload(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()
    R = cfg.RESULT_DIR

    fid = jload(R / f"fiducial_{args.tag}.json") or {}
    floor = jload(R / f"floor_{args.tag}.json") or {}
    splits = pd.read_csv(R / f"nulls_{args.tag}_splits.csv")
    paired = pd.read_csv(R / f"nulls_{args.tag}_paired.csv")
    extslope = jload(R / f"nulls_{args.tag}_extslope.json") or {}
    scaling = pd.read_csv(R / f"nulls_{args.tag}_scaling.csv")

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["parallax", "A_G", "slope_local"])
    slope = float(fid.get("slope_dMG_dMKs_median", np.nan))
    n = int(floor.get("n_stars", len(d)))

    def split(name_contains: str) -> float:
        m = splits["test"].str.contains(name_contains, case=False, regex=False)
        return float(splits.loc[m, "difference"].abs().max()) if m.any() else np.nan

    rows = []

    rows.append({
        "term": "extinction: map vs Gaia per-star A_G",
        "value_mag": float(paired["mean_difference"].abs().max()),
        "how": "paired mean shift, same stars, two extinction treatments",
    })
    rows.append({
        "term": "extinction: band law (Fitz19 vs Wang&Chen19)",
        "value_mag": float(paired.loc[paired["test"].str.contains("wangchen"),
                                      "mean_difference"].abs().max())
        if paired["test"].str.contains("wangchen").any() else np.nan,
        "how": "paired mean shift between band laws",
    })
    rows.append({
        "term": "extinction: low vs high A_0 split",
        "value_mag": split("extinction quartile"),
        "how": "two-sided null split that must return zero",
    })
    if extslope:
        rows.append({
            "term": "extinction: implied fractional error x median A_G",
            "value_mag": abs(float(extslope["implied_fractional_extinction_error"])
                             * float(np.nanmedian(d["A_G"]))),
            "how": (f"residual-vs-A_0 slope implies the correction is off by "
                    f"{100*extslope['implied_fractional_extinction_error']:.1f}% of A_G"),
        })

    rows.append({
        "term": "astrophysical: main-sequence structure vs colour",
        "value_mag": split("colour split WITHIN"),
        "how": ("colour split inside the lowest-extinction quartile; survives "
                "where there is no dust, so it is not an extinction error"),
    })

    # Parallax zero point residual, propagated through (1 - s)
    plx = d["parallax"].to_numpy(float)
    dmu = 5.0 / np.log(10) * (ZP_RESIDUAL_UAS / 1000.0) / np.nanmedian(plx)
    rows.append({
        "term": "parallax zero-point residual",
        "value_mag": float(abs(dmu * (1.0 - slope))),
        "how": (f"{ZP_RESIDUAL_UAS:.0f} uas residual (Lindegren+2021) at median "
                f"parallax {np.nanmedian(plx):.2f} mas, x |1-s| = "
                f"{abs(1-slope):.2f}"),
    })

    rows.append({
        "term": "photometric: bright vs faint G split",
        "value_mag": split("apparent G"),
        "how": "two-sided null split",
    })
    rows.append({
        "term": "crowding: sparse vs crowded split",
        "value_mag": split("crowding"),
        "how": "two-sided null split",
    })
    rows.append({
        "term": "spatial coherence plateau",
        "value_mag": float(floor.get("floor_spatial_mag", np.nan)),
        "how": "plateau of sky-patch mean-residual scatter at large N",
    })
    rows.append({
        "term": "metallicity calibration",
        "value_mag": split("metallicity"),
        "how": "metal-poor vs metal-rich split; GSP-Phot [M/H] is not "
               "independent of the BP/RP photometry",
    })

    bud = pd.DataFrame(rows)
    bud["value_mag"] = bud["value_mag"].astype(float)
    quad = float(np.sqrt(np.nansum(bud["value_mag"] ** 2)))
    largest = float(np.nanmax(bud["value_mag"]))
    naive = float(floor.get("naive_sigma_over_sqrtN_mag", np.nan))

    bud.to_csv(R / f"systematic_budget_{args.tag}.csv", index=False)

    print(f"\n=== systematic budget [{args.tag}], N = {n:,} ===\n")
    print(bud[["term", "value_mag", "how"]].to_string(
        index=False, float_format=lambda v: f"{v:9.5f}"))
    print(f"\n  quadrature sum            {quad:.5f} mag")
    print(f"  largest single term       {largest:.5f} mag")
    print(f"  naive sigma/sqrt(N)       {naive:.6f} mag")
    print(f"  budget / naive            {largest/naive:.0f}x")
    print(f"\n  implied floor on a uniform harvested fraction: "
          f"f = {1 - 10 ** (-largest / 2.5):.4f}")
    print("  (uniform f is not measurable at all self-calibrated -- this is "
          "the\n   SCALE of the coherent systematics, not a sensitivity)")

    summary = {
        "tag": args.tag, "n_stars": n,
        "quadrature_sum_mag": quad,
        "largest_term_mag": largest,
        "largest_term": str(bud.loc[bud["value_mag"].idxmax(), "term"]),
        "naive_sigma_over_sqrtN_mag": naive,
        "ratio_largest_to_naive": largest / naive if naive else np.nan,
        "implied_uniform_f": float(1 - 10 ** (-largest / 2.5)),
    }
    (R / f"systematic_budget_{args.tag}.json").write_text(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
