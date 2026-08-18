#!/usr/bin/env python
"""Search C: radio counterparts via CDS XMatch (replaces TAP cone queries).

    run.sh scripts/47b_searchC_xmatch.py --tag primary
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from astropy.table import Table
import astropy.units as u

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchC_radio")

MATCH_RADIUS_ARCSEC = 5.0
MAX_CANDIDATES = 5000
CONTROL_N = 5000
RNG_SEED = 47_2026

RADIO_CATALOGS = {
    "NVSS": {
        "vizier_id": "vizier:VIII/65/nvss",
        "flux_col": "S1.4",
        "freq_ghz": 1.4,
        "cite": "Condon+1998",
    },
    "FIRST": {
        "vizier_id": "vizier:VIII/92/first14",
        "flux_col": "Fpeak",
        "freq_ghz": 1.4,
        "cite": "Helfand+2015",
    },
}


def xmatch_catalog(stars_df, cat_name, cat_info):
    """Cross-match stars against a VizieR radio catalog via CDS XMatch."""
    from astroquery.xmatch import XMatch

    t = Table.from_pandas(stars_df[["source_id", "ra", "dec"]])
    log.info("  querying %s (%d stars)...", cat_name, len(t))
    t0 = time.time()

    try:
        result = XMatch.query(
            cat1=t,
            cat2=cat_info["vizier_id"],
            max_distance=MATCH_RADIUS_ARCSEC * u.arcsec,
            colRA1="ra",
            colDec1="dec",
        )
    except Exception as exc:
        log.error("  %s xmatch failed: %s", cat_name, exc)
        return None

    dt = time.time() - t0
    log.info("  %s: %d matches in %.1fs", cat_name, len(result), dt)

    if len(result) == 0:
        return pd.DataFrame(columns=["source_id", "angDist", "flux_mJy"])

    df = result.to_pandas()
    flux_col = cat_info["flux_col"]
    flux_candidates = [flux_col]
    if "." in flux_col:
        flux_candidates.append(flux_col.replace(".", "_"))
    found_flux = None
    for fc in flux_candidates:
        if fc in df.columns:
            found_flux = fc
            break
    if found_flux is None:
        for c in df.columns:
            if "flux" in c.lower() or c.startswith("S") or c == "Fpeak":
                found_flux = c
                break

    out = pd.DataFrame({
        "source_id": df["source_id"],
        "angDist": df["angDist"],
        "flux_mJy": df[found_flux].astype(float) if found_flux else np.nan,
    })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()
    tag = args.tag

    log.info("loading %s_resid.parquet", tag)
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{tag}_resid.parquet")
    n_all = len(d)
    log.info("loaded %d stars", n_all)

    significance = np.where(d.sigma_meas > 0, d.residual / d.sigma_meas, 0.0)
    clean = (d.ruwe < 1.4) & (d.cstar_nsigma.abs() < 3.0)
    dimmed = (d.residual > 0) & (significance > 3.0) & clean

    n_dimmed = int(dimmed.sum())
    log.info("optically dimmed (>3sigma, clean): %d (%.3f%%)",
             n_dimmed, 100 * n_dimmed / n_all)

    d_dimmed = d[dimmed].copy()
    d_dimmed["significance"] = significance[dimmed]
    d_dimmed = d_dimmed.nlargest(MAX_CANDIDATES, "significance").reset_index(drop=True)
    n_cand = len(d_dimmed)
    log.info("top %d candidates (sig range %.1f -- %.1f)",
             n_cand, d_dimmed.significance.min(), d_dimmed.significance.max())

    control_mask = (np.abs(significance) < 1.0) & clean
    d_ctrl_pool = d[control_mask].copy()
    g_min = d_dimmed.phot_g_mean_mag.quantile(0.01)
    g_max = d_dimmed.phot_g_mean_mag.quantile(0.99)
    d_ctrl_pool = d_ctrl_pool[
        (d_ctrl_pool.phot_g_mean_mag >= g_min) &
        (d_ctrl_pool.phot_g_mean_mag <= g_max)]
    rng = np.random.default_rng(RNG_SEED)
    n_ctrl = min(CONTROL_N, len(d_ctrl_pool))
    d_control = d_ctrl_pool.iloc[
        rng.choice(len(d_ctrl_pool), size=n_ctrl, replace=False)
    ].reset_index(drop=True)
    log.info("control sample: %d stars (G in [%.1f, %.1f])", n_ctrl, g_min, g_max)

    results = {}
    all_cand_matches = []

    for cat_name, cat_info in RADIO_CATALOGS.items():
        log.info("=== %s (%s, %.1f GHz) ===",
                 cat_name, cat_info["cite"], cat_info["freq_ghz"])

        cand_m = xmatch_catalog(d_dimmed, cat_name, cat_info)
        ctrl_m = xmatch_catalog(d_control, cat_name, cat_info)

        if cand_m is None and ctrl_m is None:
            results[cat_name] = {"error": "CDS XMatch unreachable"}
            continue

        n_cand_m = len(cand_m) if cand_m is not None else 0
        n_ctrl_m = len(ctrl_m) if ctrl_m is not None else 0

        cand_rate = n_cand_m / max(n_cand, 1)
        ctrl_rate = n_ctrl_m / max(n_ctrl, 1)
        excess = (cand_rate / ctrl_rate) if ctrl_rate > 0 else (
            float("inf") if n_cand_m > 0 else 0.0)

        if cand_m is not None and len(cand_m) > 0:
            for _, row in cand_m.iterrows():
                all_cand_matches.append({
                    "source_id": int(row.source_id),
                    "angDist": float(row.angDist),
                    "flux_mJy": float(row.flux_mJy),
                    "catalog": cat_name,
                    "freq_ghz": cat_info["freq_ghz"],
                })

        results[cat_name] = {
            "freq_ghz": cat_info["freq_ghz"],
            "cite": cat_info["cite"],
            "n_candidates": n_cand,
            "n_candidate_matches": n_cand_m,
            "candidate_match_rate": float(cand_rate),
            "n_control": n_ctrl,
            "n_control_matches": n_ctrl_m,
            "control_match_rate": float(ctrl_rate),
            "excess_factor": float(excess) if np.isfinite(excess) else None,
        }

        log.info("  candidates: %d/%d = %.4f", n_cand_m, n_cand, cand_rate)
        log.info("  control   : %d/%d = %.4f", n_ctrl_m, n_ctrl, ctrl_rate)
        log.info("  excess: %s",
                 f"{excess:.2f}" if np.isfinite(excess) else "inf/undef")

    unique_ids = set(m["source_id"] for m in all_cand_matches)
    n_unique = len(unique_ids)
    n_multi = sum(
        1 for sid in unique_ids
        if sum(1 for m in all_cand_matches if m["source_id"] == sid) > 1)

    total_cm = sum(r.get("n_candidate_matches", 0) for r in results.values()
                   if isinstance(r, dict) and "n_candidate_matches" in r)
    total_ctm = sum(r.get("n_control_matches", 0) for r in results.values()
                    if isinstance(r, dict) and "n_control_matches" in r)
    if total_ctm > 0:
        overall_excess = (total_cm / n_cand) / (total_ctm / n_ctrl)
    else:
        overall_excess = float("inf") if total_cm > 0 else 0.0

    if total_cm == 0 and total_ctm == 0:
        interpretation = (
            "No radio counterparts found for either dimmed candidates or "
            "controls. Main-sequence stars below the NVSS/FIRST detection "
            "threshold (~2.5 mJy) as expected. No anomalous GHz emission "
            "correlated with the optical deficit. This is the expected "
            "null for a photometric absorber that does not radiate at "
            "radio frequencies.")
    elif total_cm == 0:
        interpretation = (
            "No radio counterparts for dimmed candidates but some for controls. "
            "Controls may include chance alignments with background AGN.")
    elif overall_excess > 3.0:
        interpretation = (
            f"Radio detection rate {overall_excess:.1f}x above control. "
            f"Check individual matches against SIMBAD for AGN/active stars.")
    elif overall_excess > 1.5:
        interpretation = (
            f"Mild excess ({overall_excess:.1f}x). Likely residual contamination.")
    else:
        interpretation = (
            f"Radio detection rate ({overall_excess:.2f}x control) consistent "
            f"with chance alignments. No anomalous radio emission.")

    print(f"\n{'='*60}")
    print("SEARCH C: RADIO COUNTERPARTS (CDS XMatch)")
    print(f"{'='*60}")
    for cat_name, r in results.items():
        if "error" in r:
            print(f"  {cat_name:8s}: {r['error']}")
        else:
            print(f"  {cat_name:8s}: {r['n_candidate_matches']:3d}/{n_cand} dimmed   "
                  f"{r['n_control_matches']:3d}/{n_ctrl} control   "
                  f"excess={r.get('excess_factor', 'N/A')}")
    print(f"\n  unique dimmed with radio: {n_unique}")
    print(f"  detected in >1 survey:   {n_multi}")
    print(f"\nINTERPRETATION: {interpretation}")

    summary = {
        "tag": tag,
        "match_radius_arcsec": MATCH_RADIUS_ARCSEC,
        "n_all_stars": n_all,
        "n_optically_dimmed": n_dimmed,
        "n_candidates_queried": n_cand,
        "n_control": n_ctrl,
        "catalogs": results,
        "n_unique_radio_matches": n_unique,
        "n_multi_survey_detections": n_multi,
        "overall_excess_factor": (
            float(overall_excess) if np.isfinite(overall_excess) else None),
        "interpretation": interpretation,
    }
    json_path = cfg.RESULT_DIR / f"searchC_radio_{tag}.json"
    json_path.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", json_path)

    if all_cand_matches:
        match_df = pd.DataFrame(all_cand_matches)
        merge_cols = [c for c in ["source_id", "ra", "dec", "l", "b",
                                   "phot_g_mean_mag", "M_Ks", "residual",
                                   "sigma_meas"] if c in d_dimmed.columns]
        match_df = match_df.merge(
            d_dimmed[merge_cols].drop_duplicates(subset="source_id"),
            on="source_id", how="left", suffixes=("", "_star"))
        csv_path = cfg.RESULT_DIR / f"searchC_radio_candidates_{tag}.csv"
        match_df.to_csv(csv_path, index=False)
        log.info("wrote %d radio-matched candidates to %s",
                 len(match_df), csv_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
