#!/usr/bin/env python
"""Steps 5-6: blinded analysis, then unblinding.

    run.sh scripts/14_analyse.py --tag primary --blinded     # before freezing
    run.sh scripts/14_analyse.py --tag primary --unblind     # after committing

Under --blinded a secret constant is added to every residual, so the reported
mean and tail counts carry no information about the true answer until the
analysis is frozen and committed.
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

from pipeline import blind
from pipeline import config as cfg
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("analyse")

F_GRID = [1e-3, 3e-3, 1e-2, 3e-2, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--blinded", action="store_true")
    g.add_argument("--unblind", action="store_true")
    ap.add_argument("--k", type=float, default=cfg.FIT.outlier_nsigma)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    resid_true = d["residual"].to_numpy(dtype=float)

    if args.blinded:
        blind.create()
        s = blind.status()
        print(f"blind commitment : {s.commitment_sha256}")
        print(f"created          : {s.created_utc}")
        if s.unblinded:
            print("WARNING: this blind has already been opened. Any analysis "
                  "choice made from here is post-unblinding and must be "
                  "declared as such in RESULTS.md.")
        resid = blind.apply(resid_true)
        label = "BLINDED"
    else:
        offset = blind.unblind("I have frozen the analysis and committed it")
        print(f"unblinded. offset was {offset:+.6f} mag")
        resid = resid_true
        label = "UNBLINDED"

    shape = st.shape_stats(resid)
    print(f"\n=== {label} residual shape ===")
    print(shape.summary())
    print(shape.table.to_string(index=False))

    excl = st.exclusion_curve(resid, F_GRID, k=args.k)
    excl.to_csv(cfg.RESULT_DIR / f"exclusion_{args.tag}.csv", index=False)
    print(f"\n=== {label} 95% CL exclusion ===")
    print(excl[["f", "delta_mag", "delta_over_sigma", "efficiency",
                "n_pos_observed", "n_neg_observed", "p_upper_limit",
                "mean_f_upper_limit"]].to_string(
        index=False, float_format=lambda v: f"{v:12.6g}"))

    outl = st.find_outliers(d, resid, n_sigma=args.k)
    keep = ["source_id", "ra", "dec", "l", "b", "dist_pc", "phot_g_mean_mag",
            "M_G", "M_Ks", "bp_rp0", "A_0", "residual", "residual_sigma",
            "implied_f", "ruwe", "cstar_nsigma", "side",
            "wise_w1mpro", "wise_w2mpro", "wise_w3mpro", "wise_w4mpro",
            "allwise_designation", "tmass_designation"]
    keep = [c for c in keep if c in outl.columns]
    outl[keep].to_csv(cfg.RESULT_DIR / f"outliers_{args.tag}.csv", index=False)
    n_pos = int((outl["side"].str.startswith("positive")).sum())
    n_neg = len(outl) - n_pos
    print(f"\n=== {label} outliers at {args.k} sigma ===")
    print(f"  positive (deficit-like) : {n_pos}")
    print(f"  negative (control)      : {n_neg}")
    print(f"  written to results/outliers_{args.tag}.csv")

    # A row where p_upper_limit saturates at 1 is not a limit -- it only says
    # "every star could be harvested at this f and we would not see it".  Those
    # rows must not be quoted as a constraint on the mean fraction.
    real = excl[excl["p_upper_limit"] < 0.999]
    if len(real):
        best_mean_f = float(real["mean_f_upper_limit"].min())
        f_best = float(real.loc[real["mean_f_upper_limit"].idxmin(), "f"])
        f_min_constrained = float(real["f"].min())
    else:
        best_mean_f = f_best = f_min_constrained = float("nan")

    summary = {
        "tag": args.tag, "state": label,
        "n_stars": int(len(d)),
        "sigma_robust_mag": float(shape.sigma_robust),
        "mean_residual_mag": float(shape.mean),
        "median_residual_mag": float(shape.median),
        "bowley_skew": float(shape.skew),
        "n_pos_outliers": n_pos, "n_neg_outliers": n_neg,
        "k_sigma": float(args.k),
        "smallest_f_with_any_constraint": f_min_constrained,
        "best_mean_f_upper_limit": best_mean_f,
        "f_at_best_limit": f_best,
        "note": ("rows with p_upper_limit == 1 are unconstrained and excluded "
                 "from the best-limit search"),
    }
    (cfg.RESULT_DIR / f"analysis_{args.tag}_{label.lower()}.json").write_text(
        json.dumps(summary, indent=2))
    print(f"\nbest limit on mean harvested fraction: "
          f"{summary['best_mean_f_upper_limit']:.3g} "
          f"(at f = {summary['f_at_best_limit']:g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
