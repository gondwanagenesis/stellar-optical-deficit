#!/usr/bin/env python
"""Step 3: measure the systematic floor. MUST be run before looking at signal.

    run.sh scripts/11_null_tests.py --tag primary

Writes results/nulls_<tag>_*.csv and results/floor_<tag>.json
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

from pipeline import config as cfg
from pipeline import extinction as ext
from pipeline import fiducial as fid
from pipeline import nulls
from pipeline import sample as smp

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("nulls")


def _covs(d: pd.DataFrame, mh_degree: int, nir: bool, optical: bool = False):
    c = [(d["mh_gspphot"].to_numpy(float), mh_degree)]
    if nir and "j_ks0" in d:
        c.append((d["j_ks0"].to_numpy(float), 1))
    if optical and "bp_rp0" in d:
        c.append((d["bp_rp0"].to_numpy(float), 2))
    return c


def band_law_variant(df: pd.DataFrame, knots: int, mh_degree: int,
                     dust_map: str, band_law: str, nir: bool,
                     optical: bool = False):
    """Refit the whole chain under a different extinction treatment."""
    d = smp.add_absolute_magnitudes(df, dust_map, band_law)
    ok = (d["mh_gspphot"].notna() & np.isfinite(d["M_G"])
          & np.isfinite(d["M_Ks"])).to_numpy()
    d = d[ok].reset_index(drop=True)
    covs = _covs(d, mh_degree, nir, optical)
    fit = fid.fit_fiducial(d["M_Ks"].to_numpy(float), covs,
                           d["M_G"].to_numpy(float), knots)
    return d, fit.residuals(d["M_Ks"].to_numpy(float), covs,
                            d["M_G"].to_numpy(float))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--dust-map", default="edenhofer23")
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    meta = json.loads((cfg.RESULT_DIR / f"fiducial_{args.tag}.json").read_text())
    knots, mh_degree = meta["n_interior_knots"], meta["mh_degree"]
    nir = bool(meta.get("nir_control", False))
    optical = bool(meta.get("optical_colour_control", False))
    resid = d["residual"].to_numpy(dtype=float)
    sigma = float(meta["sigma_observed_mag"])
    log.info("%d stars, sigma=%.4f mag", len(d), sigma)

    # ---- A. two-sided splits -------------------------------------------
    splits = nulls.standard_splits(d, resid)
    splits.to_csv(cfg.RESULT_DIR / f"nulls_{args.tag}_splits.csv", index=False)
    print("\n=== A. two-sided split tests (all must be consistent with zero) ===")
    print(splits[["test", "n_a", "n_b", "mean_a", "mean_b",
                  "difference", "difference_err", "n_sigma"]].to_string(
        index=False, float_format=lambda v: f"{v:9.5f}"))

    # ---- A2. extinction slope, expressed as a physical error ------------
    exs = nulls.extinction_residual_slope(d, resid)
    if exs:
        print("\n=== A2. residual vs A_0 ===")
        for k, v in exs.items():
            print(f"  {k:38s} {v}")
        (cfg.RESULT_DIR / f"nulls_{args.tag}_extslope.json").write_text(
            json.dumps(exs, indent=2))

    # ---- B. extinction-treatment paired tests --------------------------
    print("\n=== B. paired extinction-treatment differences ===")
    paired = []
    base_map, base_law = args.dust_map, "fitz19"
    d_base, r_base = band_law_variant(d, knots, mh_degree, base_map, base_law,
                                      nir, optical)
    key = d_base["source_id"].to_numpy()
    for law in ["wangchen19"]:
        d_alt, r_alt = band_law_variant(d, knots, mh_degree, base_map, law,
                                        nir, optical)
        common, ia, ib = np.intersect1d(key, d_alt["source_id"].to_numpy(),
                                        return_indices=True)
        paired.append(nulls.dust_map_paired_test(
            d_base.iloc[ia], r_base[ia], r_alt[ib],
            f"{base_map}/{base_law}", f"{base_map}/{law}"))
    # Gaia's own per-star extinction, a genuinely different systematic
    if "ag_gspphot" in d.columns and d["ag_gspphot"].notna().any():
        d2 = d.copy()
        good = d2["ag_gspphot"].notna() & d2["ebpminrp_gspphot"].notna()
        d2 = d2[good].reset_index(drop=True)
        mu = d2["dist_mod"].to_numpy(float)
        # A_Ks from GSP-Phot A_0 with the Fitz19 law, for consistency
        a0 = d2["azero_gspphot"].to_numpy(float)
        a_ks = ext.deredden("Ks", a0, d2["bp_rp"].to_numpy(float), law="fitz19")
        d2["M_G"] = d2["phot_g_mean_mag"].to_numpy(float) - mu - d2["ag_gspphot"].to_numpy(float)
        d2["M_Ks"] = d2["tmass_ks_m"].to_numpy(float) - mu - a_ks
        covs2 = _covs(d2, mh_degree, nir, optical)
        f2 = fid.fit_fiducial(d2["M_Ks"].to_numpy(float), covs2,
                              d2["M_G"].to_numpy(float), knots)
        r2 = f2.residuals(d2["M_Ks"].to_numpy(float), covs2,
                          d2["M_G"].to_numpy(float))
        common, ia, ib = np.intersect1d(key, d2["source_id"].to_numpy(),
                                        return_indices=True)
        paired.append(nulls.dust_map_paired_test(
            d_base.iloc[ia], r_base[ia], r2[ib],
            f"{base_map}/fitz19", "GSP-Phot per-star A_G"))
    pdf = pd.DataFrame(paired)
    pdf.to_csv(cfg.RESULT_DIR / f"nulls_{args.tag}_paired.csv", index=False)
    print(pdf.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))

    # ---- C. structured group-mean scaling ------------------------------
    print("\n=== C. group-mean scatter vs group size (the floor measurement) ===")
    frames = [nulls.spatial_scaling(d, resid)]
    for col in ["A_0", "phot_g_mean_mag", "dist_pc", "sky_density"]:
        frames.append(nulls.binned_scaling(d, resid, col))
    scaling = pd.concat(frames, ignore_index=True)
    scaling.to_csv(cfg.RESULT_DIR / f"nulls_{args.tag}_scaling.csv", index=False)
    for axis, sub in scaling.groupby("axis"):
        print(f"\n  -- {axis}")
        print(sub[["n_groups", "mean_group_n", "rms_group_means",
                   "expected_if_noise", "excess"]].to_string(
            index=False, float_format=lambda v: f"{v:10.5f}"))

    # ---- D. random-subsample control -----------------------------------
    rnd = nulls.random_subsample_scaling(resid)
    rnd.to_csv(cfg.RESULT_DIR / f"nulls_{args.tag}_random.csv", index=False)
    print("\n=== D. random-subsample CONTROL (expected to track sigma/sqrt(N)) ===")
    print(rnd.to_string(index=False, float_format=lambda v: f"{v:11.6f}"))

    # ---- floor ----------------------------------------------------------
    # The floor is the plateau of `excess` at large group size, taken as the
    # median excess over groups containing at least 1% of the sample.
    big = scaling[scaling["mean_group_n"] > 0.01 * len(d)]
    floor_spatial = float(np.nanmedian(
        scaling[(scaling["axis"] == "sky (HEALPix)")
                & (scaling["mean_group_n"] > 0.01 * len(d))]["excess"]))
    floor_all = float(np.nanmedian(big["excess"])) if len(big) else np.nan
    worst_split = splits.loc[splits["n_sigma"].abs().idxmax()]

    floor = {
        "tag": args.tag,
        "n_stars": int(len(d)),
        "sigma_mag": sigma,
        "naive_sigma_over_sqrtN_mag": sigma / np.sqrt(len(d)),
        "floor_spatial_mag": floor_spatial,
        "floor_all_axes_mag": floor_all,
        "worst_split_test": str(worst_split["test"]),
        "worst_split_difference_mag": float(worst_split["difference"]),
        "worst_split_n_sigma": float(worst_split["n_sigma"]),
        "max_paired_rms_mag": float(pdf["rms"].max()) if len(pdf) else np.nan,
        "max_paired_mean_shift_mag": float(pdf["mean_difference"].abs().max())
        if len(pdf) else np.nan,
    }
    floor["floor_total_mag"] = float(np.nanmax([
        floor.get("floor_all_axes_mag") or 0.0,
        abs(floor["worst_split_difference_mag"]),
        floor.get("max_paired_mean_shift_mag") or 0.0,
    ]))
    floor["floor_implied_f"] = float(1 - 10 ** (-floor["floor_total_mag"] / 2.5))
    (cfg.RESULT_DIR / f"floor_{args.tag}.json").write_text(json.dumps(floor, indent=2))

    print("\n=== SYSTEMATIC FLOOR ===")
    for k, v in floor.items():
        print(f"  {k:32s} {v}")
    print(f"\n  naive sigma/sqrt(N) = {floor['naive_sigma_over_sqrtN_mag']:.3e} mag")
    print(f"  measured floor      = {floor['floor_total_mag']:.3e} mag")
    print(f"  ratio               = "
          f"{floor['floor_total_mag']/floor['naive_sigma_over_sqrtN_mag']:.1f}x worse")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
