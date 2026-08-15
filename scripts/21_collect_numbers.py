#!/usr/bin/env python
"""Print every headline number in one place, straight from the result files.

    run.sh scripts/21_collect_numbers.py --tag primary

Exists so that RESULTS.md is transcribed from a single authoritative dump
rather than from a dozen scrollback fragments.
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

pd.set_option("display.width", 220)
pd.set_option("display.max_rows", 200)


def jload(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def csv(p: Path):
    return pd.read_csv(p) if p.exists() else None


def section(title: str):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--alt-tag", default="primary_optcol")
    args = ap.parse_args()
    R = cfg.RESULT_DIR

    section("1. SAMPLE")
    cf = csv(R / f"cutflow_{args.tag}.csv")
    if cf is not None:
        print(cf[["cut", "n_after", "frac_removed"]].to_string(
            index=False, float_format=lambda v: f"{v:8.5f}"))

    for tag in (args.tag, args.alt_tag):
        fid = jload(R / f"fiducial_{tag}.json")
        if not fid:
            continue
        section(f"2. FIDUCIAL FIT  [{tag}]")
        for k in ["n_stars_total", "n_stars_fitted", "nir_control",
                  "optical_colour_control", "sensitive_to",
                  "n_interior_knots", "mh_degree", "n_params", "converged",
                  "sigma_observed_mag", "sigma_measurement_mag",
                  "sigma_intrinsic_mag", "sigma_intrinsic_note",
                  "sigma_observed_bare_mag",
                  "slope_dMG_dMKs_median", "slope_dMG_dMKs_p16",
                  "slope_dMG_dMKs_p84", "mean_residual", "median_residual",
                  "bowley_skew", "naive_sigma_over_sqrtN"]:
            if k in fid:
                print(f"  {k:28s} {fid[k]}")

        section(f"3. SYSTEMATIC FLOOR  [{tag}]")
        fl = jload(R / f"floor_{tag}.json")
        if fl:
            for k, v in fl.items():
                print(f"  {k:32s} {v}")
            if fl.get("naive_sigma_over_sqrtN_mag"):
                print(f"  {'RATIO floor/naive':32s} "
                      f"{fl['floor_total_mag']/fl['naive_sigma_over_sqrtN_mag']:.1f}x")
        sp = csv(R / f"nulls_{tag}_splits.csv")
        if sp is not None:
            print("\n  -- split tests")
            print(sp[["test", "n_a", "n_b", "mean_a", "mean_b", "difference",
                      "difference_err", "n_sigma"]].to_string(
                index=False, float_format=lambda v: f"{v:10.5f}"))
        ex = jload(R / f"nulls_{tag}_extslope.json")
        if ex:
            print("\n  -- residual vs A_0")
            for k, v in ex.items():
                print(f"     {k:38s} {v}")
        pr = csv(R / f"nulls_{tag}_paired.csv")
        if pr is not None:
            print("\n  -- paired extinction treatments")
            print(pr.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))
        sc = csv(R / f"nulls_{tag}_scaling.csv")
        if sc is not None:
            print("\n  -- plateau per axis (groups with >=1% of the sample)")
            n = fl["n_stars"] if fl else 0
            big = sc[sc["mean_group_n"] > 0.01 * n]
            if len(big):
                print(big.groupby("axis")["excess"].median().to_string())
        rnd = csv(R / f"nulls_{tag}_random.csv")
        if rnd is not None:
            print("\n  -- random-subsample control")
            print(rnd.to_string(index=False, float_format=lambda v: f"{v:11.6f}"))

    section("4. INJECTION-RECOVERY")
    inj = csv(R / f"injection_{args.tag}.csv")
    if inj is not None:
        print("  -- uniform")
        u = inj[inj["mode"] == "uniform"]
        print(u[["f_injected", "delta_injected", "n_stars", "mean_residual",
                 "n_pos_5sig", "anchored_f"]].to_string(
            index=False, float_format=lambda v: f"{v:13.7f}"))
        print("\n  -- sparse")
        s = inj[inj["mode"] == "sparse"]
        print(s[["f_injected", "p_injected", "p_recovered", "p_recovered_std",
                 "p_recovered_over_injected", "n_pos_5sig"]].to_string(
            index=False, float_format=lambda v: f"{v:13.7f}"))
    thr = csv(R / f"injection_threshold_{args.tag}.csv")
    if thr is not None:
        print("\n  -- smallest recoverable p (>3 sigma)")
        print(thr.to_string(index=False))

    section("5. SPECTRAL LEVERAGE")
    lv = jload(R / "spectral_leverage.json")
    if lv:
        for k, v in lv.items():
            if k != "interpretation":
                print(f"  {k:28s} {v}")
    lt = csv(R / "spectral_leverage_table.csv")
    if lt is not None:
        print()
        print(lt.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))

    section("6. MASS ANCHOR")
    ma = jload(R / "mass_anchor.json")
    if ma:
        print(pd.DataFrame(ma["subsets"]).to_string(
            index=False, float_format=lambda v: f"{v:10.4f}"))
        print(f"\n  best: {ma['best_subset']}  "
              f"flat-absorber bound f < {ma['flat_absorber_bound_frac']:.4f}")

    section("7. BLINDED / UNBLINDED ANALYSIS")
    for state in ("blinded", "unblinded"):
        a = jload(R / f"analysis_{args.tag}_{state}.json")
        if a:
            print(f"  -- {state.upper()}")
            for k, v in a.items():
                print(f"     {k:34s} {v}")
    ec = csv(R / f"exclusion_{args.tag}.csv")
    if ec is not None:
        print("\n  -- exclusion curve")
        print(ec[["f", "delta_mag", "delta_over_sigma", "efficiency",
                  "n_pos_observed", "n_neg_observed", "p_upper_limit",
                  "mean_f_upper_limit"]].to_string(
            index=False, float_format=lambda v: f"{v:13.6g}"))

    section("8. OUTLIER FOLLOW-UP")
    fu = csv(R / f"outlier_followup_{args.tag}.csv")
    if fu is not None:
        for side, g in fu.groupby("side"):
            print(f"  {side}: n={len(g)}")
            for c in ["ir_measured", "any_ir_excess", "ir_no_excess_measured",
                      "ir_unknown", "known_contaminant", "clean_candidate"]:
                if c in g:
                    print(f"     {c:26s} {int(g[c].sum())}")
        if "simbad_otype" in fu:
            print("\n  SIMBAD types (positive side):")
            pos = fu[fu["side"].str.startswith("positive")]
            print(pos["simbad_otype"].value_counts(dropna=False).head(15).to_string())

    section("9. DISTANCE TRADE STUDY")
    dt = csv(R / "distance_trade.csv")
    if dt is not None:
        print(dt.to_string(index=False, float_format=lambda v: f"{v:12.5g}"))

    section("BLIND STATE")
    cm = jload(cfg.BLIND_DIR / "commitment.json")
    if cm:
        for k, v in cm.items():
            if k != "note":
                print(f"  {k:24s} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
