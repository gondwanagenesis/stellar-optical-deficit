#!/usr/bin/env python
"""SEARCH B: optical deficit + cold IR excess -- the intercept-and-re-emit channel.

    run.sh scripts/46_searchB_cold_excess.py --tag primary

THE PHYSICS
-----------
A Dyson-type structure intercepts optical starlight and must do something with
the energy.  Search A (script 40) looks for the interception; it is agnostic
about what happens to the waste.  This script looks for the complementary
signature: the waste heat.

WISE W3 (12 um) and W4 (22 um) probe thermal emission from ~100-300 K material.
Wien's law gives T_peak = 2898/lambda_um:

    W3 (12 um) -> 242 K      W4 (22 um) -> 132 K

A partial sphere at 200-300 K (the equilibrium temperature at ~1 AU) radiates
most strongly in W3; at 100-150 K (farther out, or a cooler structure) the peak
shifts into W4.  A star that is BOTH optically dimmed (positive residual in our
M_G vs M_Ks fiducial) AND infrared-bright (negative residual in M_W3 or M_W4
relative to its photosphere) carries the full intercept-and-re-emit signature.

This is the channel Suazo et al. (2022) exploit with L = (1-gamma)L_star +
gamma*BB(T_DS), and it is where the strongest existing limits come from.  Our
contribution is (a) applying it to the same sample that Search A used, (b)
requiring the optical deficit to be independently significant, not just fitting a
single SED model, and (c) looking at W3-W4 colour anomalies for temperature
information.

WHAT THIS IS NOT
----------------
This is not a limit on the beaming-consistent class (that is script 29).  Any
candidate that survives here IS radiating thermally, so the beaming hypothesis
does not apply to them.  The two searches are complementary, not competing.

QUALITY NOTE ON W3/W4
---------------------
WISE W3 and W4 are far noisier than W1/W2.  Typical W3 SNR ~ 5-20 for our
G < 19 sample; W4 often has SNR < 5 or no detection at all.  We require
sigma < 0.3 mag (SNR > 3.6) to include a measurement, and treat each band
independently since many stars have W3 but not W4.
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

# WISE effective wavelengths (micron) from Wright et al. (2010), AJ 140, 1868.
LAM_W3 = 11.56
LAM_W4 = 22.09

# A_band / A_V for WISE W3 and W4 (Wang & Chen 2019, ApJ 877, 116, Table 3).
# Both are tiny: at A_V = 0.5 the correction is < 0.01 mag, but included for
# internal consistency with the optical/NIR dereddening.
A_W3_OVER_AV = 0.015
A_W4_OVER_AV = 0.013

# Quality and selection thresholds.
W_ERR_MAX = 0.3          # mag; W3/W4 photometric error ceiling
EXCESS_NSIGMA = 3.0      # significance threshold for IR excess (in scatter units)
DEFICIT_NSIGMA = 3.0     # significance threshold for optical deficit (in sigma)
DEFICIT_MAG_MIN = 0.05   # minimum optical residual to count as a deficit (mag)
DEFICIT_MAG_MAX = 3.0    # sanity cap on optical residual (mag)

# Star-forming regions to mask (same list as 44_searchA_null.py).
SFR = [("Orion", 209.0, -19.4, 12.0), ("Ophiuchus", 353.6, 16.9, 8.0),
       ("Corona Australis", 359.7, -17.8, 6.0), ("Taurus", 172.5, -15.5, 10.0),
       ("Chamaeleon", 297.2, -15.6, 8.0), ("Lupus", 339.0, 15.0, 8.0),
       ("Perseus", 159.4, -20.0, 8.0), ("Serpens", 31.5, 5.3, 6.0),
       ("Cepheus", 110.0, 15.0, 10.0), ("Vela", 265.0, 1.5, 10.0)]


def in_sfr(l: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Galactic-coordinate mask for known star-forming regions."""
    lab = np.zeros(len(l), dtype=bool)
    for _, l0, b0, rad in SFR:
        dl = np.abs((l - l0 + 180) % 360 - 180)
        lab |= np.hypot(dl * np.cos(np.radians(b)), b - b0) < rad
    return lab


def running_fiducial(x: np.ndarray, y: np.ndarray, nbins: int = 60,
                     lo: float | None = None, hi: float | None = None
                     ) -> np.ndarray:
    """Running median of y as a function of x -- the empirical photospheric locus.

    Returns a y-prediction for every input x via linear interpolation between
    bin medians.  Bins with < 50 stars are dropped (W3/W4 coverage is sparse
    at the faint end).
    """
    lo = np.nanpercentile(x, 1) if lo is None else lo
    hi = np.nanpercentile(x, 99) if hi is None else hi
    edges = np.linspace(lo, hi, nbins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, nbins - 1)
    med = np.full(nbins, np.nan)
    for b in range(nbins):
        m = (idx == b) & np.isfinite(y)
        if m.sum() > 50:
            med[b] = np.median(y[m])
    good = np.isfinite(med)
    if good.sum() < 5:
        return np.full(len(x), np.nan)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return np.interp(x, centres[good], med[good])


def polyfit_fiducial(x: np.ndarray, y: np.ndarray, valid: np.ndarray,
                     degree: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Robust polynomial fit of y vs x for valid points.

    Returns (predicted_y_for_all, coefficients).  Uses iterative sigma-clipping
    to suppress outliers (debris discs, YSOs) that would bias the photospheric
    relation.
    """
    xv, yv = x[valid], y[valid]
    mask = np.ones(len(xv), dtype=bool)
    coeffs = None
    for _iteration in range(5):
        xm, ym = xv[mask], yv[mask]
        if len(xm) < degree + 10:
            break
        coeffs = np.polyfit(xm, ym, degree)
        pred = np.polyval(coeffs, xv)
        resid = yv - pred
        sigma = st.robust_sigma(resid)
        mask = np.abs(resid) < 3.0 * sigma
    if coeffs is None:
        return np.full(len(x), np.nan), np.array([])
    return np.polyval(coeffs, x), coeffs


def main() -> int:
    ap = argparse.ArgumentParser(description="Search B: optical deficit + cold IR excess")
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--deficit-nsigma", type=float, default=DEFICIT_NSIGMA)
    ap.add_argument("--excess-nsigma", type=float, default=EXCESS_NSIGMA)
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Load sample
    # ------------------------------------------------------------------
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    n_all = len(d)
    print(f"loaded {n_all:,} stars from {args.tag}_resid.parquet")

    # Unpack the columns we need.
    r = d["residual"].to_numpy(float)
    sigma_r = st.robust_sigma(r)
    med_r = float(np.median(r))
    m_ks = d["M_Ks"].to_numpy(float)
    dist_pc = d["dist_pc"].to_numpy(float)
    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float) if "bp_rp" in d.columns else np.full(n_all, np.nan)

    # Absolute WISE magnitudes.  Distance modulus uses the same dist_pc that
    # produced M_G and M_Ks, so the calibration is internally consistent.
    dm = 5.0 * np.log10(np.clip(dist_pc, 1.0, None)) - 5.0

    w3 = d["wise_w3mpro"].to_numpy(float)
    w4 = d["wise_w4mpro"].to_numpy(float)
    e_w3 = d["wise_w3mpro_error"].to_numpy(float) if "wise_w3mpro_error" in d.columns else np.full(n_all, np.nan)
    e_w4 = d["wise_w4mpro_error"].to_numpy(float) if "wise_w4mpro_error" in d.columns else np.full(n_all, np.nan)

    # Deredden W3 and W4.
    w3_dered = w3 - A_W3_OVER_AV * a0
    w4_dered = w4 - A_W4_OVER_AV * a0
    M_W3 = w3_dered - dm
    M_W4 = w4_dered - dm

    # ------------------------------------------------------------------
    # Quality masks
    # ------------------------------------------------------------------
    have_w3 = np.isfinite(w3) & np.isfinite(e_w3) & (e_w3 < W_ERR_MAX) & (e_w3 > 0)
    have_w4 = np.isfinite(w4) & np.isfinite(e_w4) & (e_w4 < W_ERR_MAX) & (e_w4 > 0)
    have_either = have_w3 | have_w4

    print(f"\nWISE coverage:")
    print(f"  W3 usable (finite, sigma < {W_ERR_MAX}) : {have_w3.sum():,} ({100*have_w3.mean():.1f}%)")
    print(f"  W4 usable (finite, sigma < {W_ERR_MAX}) : {have_w4.sum():,} ({100*have_w4.mean():.1f}%)")
    print(f"  either W3 or W4                    : {have_either.sum():,}")
    print(f"  both W3 and W4                     : {(have_w3 & have_w4).sum():,}")

    # Star-forming region mask.
    l_arr = d["l"].to_numpy(float)
    b_arr = d["b"].to_numpy(float)
    sfr_mask = in_sfr(l_arr, b_arr)
    print(f"  in star-forming regions            : {sfr_mask.sum():,}")

    # ------------------------------------------------------------------
    # Fit M_W3 and M_W4 vs M_Ks: the photospheric fiducial
    # ------------------------------------------------------------------
    # Clean training sample: not in SFR, low extinction, off the plane,
    # good RUWE, and not an optical outlier (so the fiducial is trained on
    # normal photospheres, not dimmed stars).
    clean = (~sfr_mask & (a0 < 0.15) & (np.abs(b_arr) > 15)
             & (np.abs(r - med_r) < 3.0 * sigma_r))
    if "ruwe" in d.columns:
        clean &= d["ruwe"].to_numpy(float) < 1.4

    print(f"\n--- W3 fiducial fit (M_W3 vs M_Ks) ---")
    train_w3 = clean & have_w3 & np.isfinite(m_ks)
    print(f"  training sample: {train_w3.sum():,} stars")
    # Use polynomial fit (degree 4): M_W3 vs M_Ks is smooth and nearly linear
    # for main-sequence stars, with mild curvature at the ends.
    pred_w3, coeff_w3 = polyfit_fiducial(m_ks, M_W3, train_w3, degree=4)
    resid_w3 = M_W3 - pred_w3   # negative = brighter than expected = IR excess
    sigma_w3 = st.robust_sigma(resid_w3[have_w3 & np.isfinite(pred_w3)])
    print(f"  scatter of W3 residual: {sigma_w3:.4f} mag")
    print(f"  poly coefficients: {coeff_w3}")

    print(f"\n--- W4 fiducial fit (M_W4 vs M_Ks) ---")
    train_w4 = clean & have_w4 & np.isfinite(m_ks)
    print(f"  training sample: {train_w4.sum():,} stars")
    pred_w4, coeff_w4 = polyfit_fiducial(m_ks, M_W4, train_w4, degree=4)
    resid_w4 = M_W4 - pred_w4
    sigma_w4 = st.robust_sigma(resid_w4[have_w4 & np.isfinite(pred_w4)])
    print(f"  scatter of W4 residual: {sigma_w4:.4f} mag")
    print(f"  poly coefficients: {coeff_w4}")

    # Also compute a running-median fiducial as a cross-check.
    run_w3 = running_fiducial(m_ks, np.where(train_w3, M_W3, np.nan))
    run_w4 = running_fiducial(m_ks, np.where(train_w4, M_W4, np.nan))

    # ------------------------------------------------------------------
    # Flag IR excess and optical deficit
    # ------------------------------------------------------------------
    # IR EXCESS: star is brighter than expected in W3 or W4.
    # In magnitudes, brighter = smaller number = negative residual.
    ir_excess_w3 = have_w3 & (resid_w3 < -args.excess_nsigma * sigma_w3)
    ir_excess_w4 = have_w4 & (resid_w4 < -args.excess_nsigma * sigma_w4)
    ir_excess_any = ir_excess_w3 | ir_excess_w4

    # OPTICAL DEFICIT: star is dimmer than expected in G (positive residual).
    opt_deficit = (r > med_r + args.deficit_nsigma * sigma_r) & \
                  (r > DEFICIT_MAG_MIN) & (r < DEFICIT_MAG_MAX)

    print(f"\n--- selection ---")
    print(f"  optical deficit (>{args.deficit_nsigma:.0f} sigma, >{DEFICIT_MAG_MIN} mag) : {opt_deficit.sum():,}")
    print(f"  W3 excess (< -{args.excess_nsigma:.0f} sigma)                   : {ir_excess_w3.sum():,}")
    print(f"  W4 excess (< -{args.excess_nsigma:.0f} sigma)                   : {ir_excess_w4.sum():,}")
    print(f"  any IR excess (W3 or W4)                      : {ir_excess_any.sum():,}")

    # ------------------------------------------------------------------
    # THE INTERCEPT-AND-RE-EMIT CANDIDATES
    # ------------------------------------------------------------------
    # Stars that are BOTH optically dimmed AND infrared-bright.
    intercept_w3 = opt_deficit & ir_excess_w3 & ~sfr_mask
    intercept_w4 = opt_deficit & ir_excess_w4 & ~sfr_mask
    intercept_any = opt_deficit & ir_excess_any & ~sfr_mask
    # Without SFR cut (for comparison).
    intercept_any_raw = opt_deficit & ir_excess_any

    n_int_w3 = int(intercept_w3.sum())
    n_int_w4 = int(intercept_w4.sum())
    n_int_any = int(intercept_any.sum())
    n_int_raw = int(intercept_any_raw.sum())
    print(f"\n=== INTERCEPT-AND-RE-EMIT CANDIDATES ===")
    print(f"  optical deficit + W3 excess (no SFR) : {n_int_w3:,}")
    print(f"  optical deficit + W4 excess (no SFR) : {n_int_w4:,}")
    print(f"  optical deficit + any IR excess (no SFR) : {n_int_any:,}")
    print(f"  (before SFR cut: {n_int_raw:,})")

    # ------------------------------------------------------------------
    # W3-W4 colour anomaly: warm dust/shell temperature diagnostic
    # ------------------------------------------------------------------
    have_both = have_w3 & have_w4
    w3w4 = w3_dered - w4_dered   # colour; positive = W4 brighter = cooler excess
    w3w4_exp = running_fiducial(m_ks, np.where(clean & have_both, w3w4, np.nan))
    w3w4_excess = w3w4 - w3w4_exp   # positive = redder than photosphere
    sigma_w3w4 = st.robust_sigma(w3w4_excess[have_both & np.isfinite(w3w4_exp)])

    colour_red = have_both & (w3w4_excess > args.excess_nsigma * sigma_w3w4)
    colour_red_deficit = colour_red & opt_deficit & ~sfr_mask
    print(f"\n--- W3-W4 colour anomaly ---")
    print(f"  W3-W4 scatter: {sigma_w3w4:.4f} mag")
    print(f"  W3-W4 red (> {args.excess_nsigma:.0f} sigma): {colour_red.sum():,}")
    print(f"  W3-W4 red + optical deficit (no SFR): {int(colour_red_deficit.sum()):,}")

    # Estimate implied dust temperature from W3-W4 colour for candidates.
    if n_int_any > 0:
        cands = d.loc[intercept_any].copy()
        cands["resid_optical"] = r[intercept_any]
        cands["resid_optical_nsigma"] = (r[intercept_any] - med_r) / sigma_r
        cands["resid_w3"] = resid_w3[intercept_any]
        cands["resid_w4"] = resid_w4[intercept_any]
        cands["resid_w3_nsigma"] = resid_w3[intercept_any] / sigma_w3
        cands["resid_w4_nsigma"] = resid_w4[intercept_any] / sigma_w4
        cands["w3w4_colour"] = w3w4[intercept_any]
        cands["w3w4_excess"] = w3w4_excess[intercept_any]
        cands["implied_f_optical"] = st.fraction_from_delta(
            np.clip(r[intercept_any] - med_r, 0, None))
    else:
        cands = pd.DataFrame()

    # ------------------------------------------------------------------
    # MIRROR CONTROL: bright-in-G + IR-deficit
    # ------------------------------------------------------------------
    # Nothing physical makes a star anomalously bright in G AND faint in the
    # mid-IR.  This mirror count is the false-positive rate under identical cuts.
    opt_bright = (r < med_r - args.deficit_nsigma * sigma_r) & \
                 (r < -DEFICIT_MAG_MIN) & (r > -DEFICIT_MAG_MAX)
    ir_deficit_w3 = have_w3 & (resid_w3 > args.excess_nsigma * sigma_w3)
    ir_deficit_w4 = have_w4 & (resid_w4 > args.excess_nsigma * sigma_w4)
    mirror = opt_bright & (ir_deficit_w3 | ir_deficit_w4) & ~sfr_mask
    n_mirror = int(mirror.sum())
    ratio = n_int_any / max(n_mirror, 1)
    print(f"\n=== MIRROR CONTROL (bright-G + IR-faint) ===")
    print(f"  mirror count : {n_mirror:,}")
    print(f"  signal/mirror ratio : {ratio:.2f}")

    if ratio < 1.5:
        mirror_verdict = ("NOISE. The mirror yields a comparable count; the "
                          "intercept-and-re-emit candidates are consistent with "
                          "the pipeline's false-positive rate.")
    elif ratio < 3:
        mirror_verdict = ("MARGINAL. Mild asymmetry only; most candidates are "
                          "likely false positives.")
    else:
        mirror_verdict = (f"ASYMMETRIC at {ratio:.1f}:1. The physical-sign set "
                          f"exceeds the unphysical mirror, so the signal is at "
                          f"least not purely noise. Astrophysical contaminants "
                          f"(YSOs, debris discs) are the likely source.")
    print(f"  VERDICT: {mirror_verdict}")

    # ------------------------------------------------------------------
    # Prevalence limit: what fraction of stars could host cold re-emitters?
    # ------------------------------------------------------------------
    # Among stars with usable W3 or W4, the fraction showing BOTH the optical
    # deficit and the IR excess is the upper bound on the intercept-and-re-emit
    # population.  Apply the Poisson upper limit for a proper statistical bound.
    n_wise_clean = int((have_either & ~sfr_mask).sum())
    p_raw = n_int_any / max(n_wise_clean, 1)
    p_ul = st.poisson_upper_limit(n_int_any) / max(n_wise_clean, 1)
    print(f"\n=== PREVALENCE ===")
    print(f"  stars with usable W3/W4 (no SFR) : {n_wise_clean:,}")
    print(f"  intercept-and-re-emit candidates  : {n_int_any:,}")
    print(f"  raw fraction                      : {p_raw:.3e}")
    print(f"  95% CL Poisson upper limit        : {min(p_ul, 1.0):.3e}")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    # Candidate table.
    if len(cands) > 0:
        keep = [c for c in [
            "source_id", "ra", "dec", "l", "b", "dist_pc",
            "phot_g_mean_mag", "M_G", "M_Ks", "A_0",
            "wise_w3mpro", "wise_w4mpro", "wise_w3mpro_error", "wise_w4mpro_error",
            "resid_optical", "resid_optical_nsigma",
            "resid_w3", "resid_w4", "resid_w3_nsigma", "resid_w4_nsigma",
            "w3w4_colour", "w3w4_excess", "implied_f_optical",
            "allwise_designation", "tmass_designation",
        ] if c in cands.columns]
        out_csv = cfg.RESULT_DIR / f"searchB_cold_candidates_{args.tag}.csv"
        cands[keep].sort_values("resid_optical_nsigma", ascending=False).to_csv(
            out_csv, index=False)
        print(f"\n  candidates written to {out_csv.name}")

        # Show the top candidates.
        show = [c for c in [
            "source_id", "dist_pc", "M_Ks",
            "resid_optical_nsigma", "resid_w3_nsigma", "resid_w4_nsigma",
            "w3w4_colour", "implied_f_optical",
        ] if c in cands.columns]
        top = cands[keep].nlargest(20, "resid_optical_nsigma")
        print(f"\n  top 20 candidates:")
        print(top[show].to_string(index=False))

    # Summary JSON.
    summary = {
        "tag": args.tag,
        "n_all": n_all,
        "n_with_w3": int(have_w3.sum()),
        "n_with_w4": int(have_w4.sum()),
        "n_with_either": int(have_either.sum()),
        "n_with_both": int((have_w3 & have_w4).sum()),
        "sigma_optical": float(sigma_r),
        "sigma_w3": float(sigma_w3),
        "sigma_w4": float(sigma_w4),
        "sigma_w3w4_colour": float(sigma_w3w4),
        "w3_fiducial_poly_degree": 4,
        "w3_fiducial_coefficients": coeff_w3.tolist() if len(coeff_w3) else [],
        "w4_fiducial_poly_degree": 4,
        "w4_fiducial_coefficients": coeff_w4.tolist() if len(coeff_w4) else [],
        "deficit_nsigma": args.deficit_nsigma,
        "excess_nsigma": args.excess_nsigma,
        "n_optical_deficit": int(opt_deficit.sum()),
        "n_ir_excess_w3": int(ir_excess_w3.sum()),
        "n_ir_excess_w4": int(ir_excess_w4.sum()),
        "n_ir_excess_any": int(ir_excess_any.sum()),
        "n_intercept_reemit_w3": n_int_w3,
        "n_intercept_reemit_w4": n_int_w4,
        "n_intercept_reemit_any": n_int_any,
        "n_intercept_reemit_before_sfr_cut": n_int_raw,
        "n_mirror_control": n_mirror,
        "signal_mirror_ratio": float(ratio),
        "mirror_verdict": mirror_verdict,
        "n_w3w4_colour_red_deficit": int(colour_red_deficit.sum()),
        "n_wise_clean_sample": n_wise_clean,
        "p_raw": float(p_raw),
        "p_UL_95": float(min(p_ul, 1.0)),
        "interpretation": (
            "Search B looks for the intercept-and-re-emit signature: optical "
            "dimming (positive residual in M_G vs M_Ks) paired with mid-IR "
            "excess (negative residual in M_W3 or M_W4 vs M_Ks). This is the "
            "channel that constrains structures radiating waste heat "
            "isotropically at 100-300 K. p_UL_95 is the 95% CL Poisson upper "
            "limit on the fraction of stars hosting such structures, among "
            "those with usable WISE W3/W4 photometry outside star-forming "
            "regions. Astrophysical contaminants (debris discs, YSOs) dominate "
            "any real signal at this sensitivity."),
    }
    out_json = cfg.RESULT_DIR / f"searchB_cold_{args.tag}.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"\n  summary written to {out_json.name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
