#!/usr/bin/env python
"""Diagnose the Search-A 28.6:1 asymmetry: what population produces it?

    run.sh scripts/49_searchA_diagnosis.py --tag primary

THE QUESTION
------------
Search A (scripts 40-44) found 372 dimmed vs 13 brightened stars under identical
cuts, a 28.6:1 asymmetry.  The survivors are not metal-poor subdwarfs (refuted
by script 43) and are not obviously grain growth (high latitude, low extinction).

Three hypotheses remain:
  1. Unresolved NIR companions inflating K_s (aperture-mismatch family)
  2. Debris discs adding NIR flux (less common than circumstellar material)
  3. Something actually absorbing optical light with a non-dust spectrum

This script tests each hypothesis against the data we already have, without
needing spectroscopy.

TESTS
-----
  A. Cross-match quality: if K_s is inflated by a companion, the 2MASS/WISE
     positional match should be systematically offset.
  B. W1-W2 colour: a cool companion contributes in W2 > W1, making W1-W2 red.
  C. BP-RP excess factor: optical blending inflates this quantity.
  D. Teff distribution: M dwarfs are more susceptible to unresolved binarity.
  E. Correlation structure: if the mechanism is blend-related, the deficit should
     correlate with cross-match distance but NOT with IR colour.
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
    clean = s[~s["known_contaminant"].astype(bool)][
        ["source_id", "alpha", "deficit_G_mag"]
    ].copy()

    cols = [
        "source_id", "l", "b", "dist_pc", "M_G", "M_Ks", "residual",
        "A_0", "ruwe", "sigma_meas", "teff_gspphot", "mh_gspphot",
        "tmass_xm_nmates", "tmass_xm_nnb", "tmass_xm_dist",
        "wise_xm_nmates", "wise_xm_nnb", "wise_xm_dist",
        "non_single_star", "astrometric_excess_noise_sig",
        "ipd_frac_multi_peak", "bp_rp0",
        "wise_w1mpro", "wise_w2mpro", "wise_w3mpro",
        "tmass_j_m", "tmass_h_m", "tmass_ks_m",
        "phot_bp_rp_excess_factor",
    ]
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet", columns=cols)
    m = clean.merge(d, on="source_id", how="left")

    np.random.seed(42)
    ctrl = d.sample(n=min(len(m) * 10, len(d)), random_state=42)

    results = {}

    # ---- A. Cross-match quality ----
    print("=== A. Cross-match quality ===")
    for name, col in [("tmass_xm_dist", "tmass_xm_dist"),
                       ("wise_xm_dist", "wise_xm_dist")]:
        ms = float(m[col].median())
        mc = float(ctrl[col].median())
        ratio = ms / max(mc, 1e-9)
        print(f"  {name:18s} survivors={ms:.3f}\"  control={mc:.3f}\"  ratio={ratio:.2f}")
        results[f"{name}_surv"] = ms
        results[f"{name}_ctrl"] = mc
        results[f"{name}_ratio"] = ratio

    # ---- B. W1-W2 colour ----
    print("\n=== B. W1-W2 colour (companion diagnostic) ===")
    w12_s = (m.wise_w1mpro - m.wise_w2mpro).dropna()
    w12_c = (ctrl.wise_w1mpro - ctrl.wise_w2mpro).dropna()
    print(f"  W1-W2 survivors median: {w12_s.median():.3f} mag")
    print(f"  W1-W2 control   median: {w12_c.median():.3f} mag")
    print(f"  frac W1-W2 > 0.1: surv={float((w12_s > 0.1).mean()):.3f}  "
          f"ctrl={float((w12_c > 0.1).mean()):.3f}")
    results["w12_surv"] = float(w12_s.median())
    results["w12_ctrl"] = float(w12_c.median())

    # ---- C. BP-RP excess factor ----
    print("\n=== C. BP-RP excess factor (optical blending) ===")
    bpe_s = float(m.phot_bp_rp_excess_factor.median())
    bpe_c = float(ctrl.phot_bp_rp_excess_factor.median())
    print(f"  survivors median: {bpe_s:.4f}")
    print(f"  control   median: {bpe_c:.4f}")
    results["bprp_excess_surv"] = bpe_s
    results["bprp_excess_ctrl"] = bpe_c

    # ---- D. Teff distribution ----
    print("\n=== D. Temperature distribution ===")
    ts = float(m.teff_gspphot.median())
    tc = float(ctrl.teff_gspphot.median())
    print(f"  Teff survivors median: {ts:.0f} K")
    print(f"  Teff control   median: {tc:.0f} K")
    frac_m_s = float(m.teff_gspphot.between(3000, 4000).mean())
    frac_m_c = float(ctrl.teff_gspphot.between(3000, 4000).mean())
    print(f"  frac 3000-4000K: surv={frac_m_s:.3f}  ctrl={frac_m_c:.3f}")
    results["teff_surv"] = ts
    results["teff_ctrl"] = tc
    results["frac_mdwarf_surv"] = frac_m_s
    results["frac_mdwarf_ctrl"] = frac_m_c

    # ---- E. Correlation structure ----
    print("\n=== E. Correlation structure ===")
    deficit = m.deficit_G_mag.to_numpy(float)

    for name, col in [("wise_xm_dist", "wise_xm_dist"),
                       ("W1-W2", None),
                       ("tmass_xm_dist", "tmass_xm_dist")]:
        if col is not None:
            x = m[col].to_numpy(float)
        else:
            x = (m.wise_w1mpro - m.wise_w2mpro).to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(deficit)
        if mask.sum() > 10:
            r = np.corrcoef(deficit[mask], x[mask])[0, 1]
            print(f"  deficit vs {name:18s}: r = {r:+.3f}  (n={mask.sum()})")
            results[f"corr_deficit_{name}"] = float(r)

    # ---- Verdict ----
    print("\n" + "=" * 60)
    wise_ratio = results.get("wise_xm_dist_ratio", 1.0)
    w12_excess = results["w12_surv"] - results["w12_ctrl"]

    evidence = []
    if wise_ratio > 2.0:
        evidence.append(f"WISE positional offset {wise_ratio:.1f}x larger")
    if w12_excess > 0.03:
        evidence.append(f"W1-W2 {w12_excess:.3f} mag redder")
    if results["frac_mdwarf_surv"] > results["frac_mdwarf_ctrl"] * 1.3:
        evidence.append("concentrated in M dwarfs")
    if results.get("bprp_excess_surv", 0) > results.get("bprp_excess_ctrl", 0) * 1.02:
        evidence.append("slightly elevated BP-RP excess factor")

    if len(evidence) >= 3:
        verdict = ("CONSISTENT WITH UNRESOLVED NIR COMPANIONS. "
                   + "; ".join(evidence)
                   + ". The asymmetry is real (something physical dims these stars "
                   "with a flat slope), but the multi-indicator pattern matches "
                   "the aperture-mismatch/companion family rather than engineering. "
                   "Spectroscopy remains the decisive test.")
    elif len(evidence) >= 1:
        verdict = ("PARTIALLY CONSISTENT WITH NIR COMPANIONS. "
                   + "; ".join(evidence)
                   + ". Some blend indicators are elevated but the pattern is "
                   "not conclusive. Spectroscopic follow-up needed.")
    else:
        verdict = ("NO CLEAR BLEND SIGNATURE. The survivors do not show elevated "
                   "cross-match offsets, IR colour anomalies, or blending indicators. "
                   "The mechanism remains unexplained.")

    print(f"VERDICT: {verdict}")
    results["evidence"] = evidence
    results["verdict"] = verdict
    results["n_survivors"] = len(m)

    (cfg.RESULT_DIR / f"searchA_diagnosis_{args.tag}.json").write_text(
        json.dumps(results, indent=2))
    print(f"\nwrote searchA_diagnosis_{args.tag}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
