#!/usr/bin/env python
"""AUDIT 3: is channel 5's spectral slope actually measured, or invented?

    run.sh scripts/91_audit_spectral_slope.py --tag primary

WHY THIS CHANNEL
----------------
Channel 5 (scripts/40) is the largest-N channel in the project -- 2,991,398
stars, a rate upper limit of 1.0e-6 -- and the only one that ever produced a
survivor list (513 stars, 372 after the mirror). Every downstream claim rests
on one quantity, alpha, the power-law index of the fitted absorber:

    "median alpha = 1.95, squarely on the dust value"
    "5,239 stars have alpha significantly below 1.5"
    "372 asymmetric survivors"

None of those mean anything unless the fitter can RECOVER alpha. That has never
been tested. Two bugs of exactly this shape have already been caught in this
project: a hand-rolled C* formula that rejected 91% of the sample, and a
cross-match that ignored proper motion. Both were reimplementations of
something whose behaviour nobody had measured.

WHAT THIS AUDIT DOES
--------------------
T0  REPRODUCTION. Refit the real data with the fitter block transcribed from
    scripts/40 and compare against the stored spectral_slope parquet. If this
    does not match, the audit is testing something other than the pipeline and
    every later number here is void. This is the guard against auditing a
    reimplementation instead of the code.

T1  INJECTION-RECOVERY. Add a synthetic absorber of known (alpha, amplitude)
    to real residuals and refit. Measures three things nobody has:
      - bias:   does recovered alpha equal injected alpha?
      - power:  what fraction of a genuinely flat (alpha=0.5) absorber is
                flagged by the "alpha_upper_2sig < 1.5" criterion?
      - TYPE-S: what fraction of a pure DUST absorber (alpha=2.0) is
                MISCLASSIFIED as significantly flat? If that rate times 79,426
                is of order 5,239, the flat population is dust plus noise and
                the grain-growth diagnosis in scripts/41 was diagnosing the
                fitter, not the Galaxy.

T2  THE MISSING SLOPE MIRROR. This project's own standard says an unphysical
    opposite-sign tail must be measured -- Search T carries beta > 3 as the
    mirror for beta ~ 0. Channel 5 never did the analogous thing. Grain physics
    cannot make alpha much steeper than ~2 any more than it can make it much
    flatter than ~1.5, so "alpha significantly ABOVE 2.6" is an equally
    unphysical population under identical cuts, and its count is the slope
    axis's own false-positive rate. scripts/44's mirror flips the AMPLITUDE
    sign, which tests whether something dims; it does not test the slope, and
    the slope is the entire specificity claim.

T3  MARK PERMUTATION, two kinds.
      A. Whole residual 7-vectors shuffled over stars. The fit is unchanged by
         construction, so this isolates the sky-and-quality cuts: does the
         survivor set's association with A_0, |b| and SFR membership survive?
      B. Each band shuffled INDEPENDENTLY over stars. This preserves every
         band's marginal distribution and the survey geometry but destroys the
         cross-band coherence that a real SED shape has. Refit. Whatever
         survives is what incoherent noise alone can manufacture.
    Neither generates synthetic positions, per the standing rule.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import extinction as ext
from pipeline import statistics as st

# Import the CHANNEL'S OWN module rather than restating its constants. The
# numeric filename prefix is not a legal identifier, so go through importlib.
_spec = importlib.util.spec_from_file_location(
    "ch5", Path(__file__).resolve().parent / "40_spectral_slope_search.py")
ch5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ch5)

LAM, FITBANDS, AV_RATIO = ch5.LAM, ch5.FITBANDS, ch5.AV_RATIO
ALPHA_GRID = ch5.ALPHA_GRID
running_pred = ch5.running_pred

# Star-forming complexes, verbatim from scripts/41 and scripts/44.
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


def build_basis():
    """The normalised absorber basis, as scripts/40 constructs it."""
    lam = np.array([LAM[b] for b in FITBANDS])
    iG = FITBANDS.index("G")
    basis = np.array([(lam / LAM["Ks"]) ** (-a) - 1.0 for a in ALPHA_GRID])
    return basis / basis[:, iG][:, None]


def fit_grid(R, sel, sig, basis, want_lo=False, chunk=200_000):
    """The fitter block of scripts/40, transcribed verbatim, plus the LOWER
    edge of the profile interval (which scripts/40 never computes -- that is
    T2's whole point).  T0 checks this transcription against the stored
    parquet before any conclusion is drawn from it."""
    w = 1.0 / sig ** 2
    den = (basis ** 2 * w).sum(axis=1)
    bw = (basis * w).T.astype(np.float32)
    n = len(sel)
    alpha = np.empty(n, dtype=np.float32)
    Abest = np.empty(n, dtype=np.float32)
    chi2best = np.empty(n, dtype=np.float32)
    dchi2 = np.empty(n, dtype=np.float32)
    a_hi = np.empty(n, dtype=np.float32)
    a_lo = np.empty(n, dtype=np.float32) if want_lo else None
    idx = np.arange(len(ALPHA_GRID))
    for s0 in range(0, n, chunk):
        s1 = min(s0 + chunk, n)
        Rs = R[sel[s0:s1]].astype(np.float32)
        num = Rs @ bw
        c0 = (Rs ** 2 * w.astype(np.float32)).sum(axis=1)
        chi2 = c0[:, None] - num ** 2 / den.astype(np.float32)
        ib = np.argmin(chi2, axis=1)
        r_ = np.arange(s1 - s0)
        alpha[s0:s1] = ALPHA_GRID[ib]
        Abest[s0:s1] = num[r_, ib] / den[ib]
        chi2best[s0:s1] = chi2[r_, ib]
        dchi2[s0:s1] = c0 - chi2[r_, ib]
        ok_a = chi2 <= (chi2[r_, ib] + 4.0)[:, None]
        a_hi[s0:s1] = np.where(ok_a.any(axis=1),
                               ALPHA_GRID[(ok_a * idx).argmax(axis=1)], np.nan)
        if want_lo:
            a_lo[s0:s1] = np.where(ok_a.any(axis=1),
                                   ALPHA_GRID[ok_a.argmax(axis=1)], np.nan)
        del Rs, num, chi2, ok_a
    out = {"alpha": alpha, "deficit_G_mag": Abest, "chi2_fit": chi2best,
           "dchi2": dchi2, "alpha_upper_2sig": a_hi}
    if want_lo:
        out["alpha_lower_2sig"] = a_lo
    return out


def build_residual_matrix(d):
    """Residual construction of scripts/40, lines 100-137, transcribed."""
    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    mu = d["dist_mod"].to_numpy(float)
    absmag = {}
    for b in ("BP", "G", "RP"):
        col = {"BP": "phot_bp_mean_mag", "G": "phot_g_mean_mag",
               "RP": "phot_rp_mean_mag"}[b]
        absmag[b] = d[col].to_numpy(float) - mu - ext.deredden(b, a0, bp_rp)
    for b in ("J", "H", "Ks"):
        col = {"J": "tmass_j_m", "H": "tmass_h_m", "Ks": "tmass_ks_m"}[b]
        absmag[b] = d[col].to_numpy(float) - mu - ext.deredden(b, a0, bp_rp)
    for b in ("W1", "W2"):
        col = {"W1": "wise_w1mpro", "W2": "wise_w2mpro"}[b]
        absmag[b] = d[col].to_numpy(float) - mu - AV_RATIO[b] * a0

    jh = absmag["J"] - absmag["H"]
    mks = absmag["Ks"]
    ok = np.isfinite(jh) & np.isfinite(mks)
    for b in FITBANDS:
        ok &= np.isfinite(absmag[b])
    ok &= ((d["wise_w1mpro_error"].to_numpy(float) < 0.15)
           & (d["wise_w2mpro_error"].to_numpy(float) < 0.15))

    resid = {b: absmag[b] - running_pred(mks, jh, absmag[b], ok)
             for b in FITBANDS}
    R = np.column_stack([resid[b] for b in FITBANDS])
    sig = np.array([st.robust_sigma(resid[b][ok]) for b in FITBANDS])
    R = R - np.nanmedian(R[ok], axis=0)
    return R, sig, ok


# ---- the downstream cut stack, exactly as scripts/41 and scripts/44 apply it
def survivor_cuts(f, l, b, a0, low=True):
    edge = (f["alpha"] < -0.9) | (f["alpha"] > 5.9)
    common = (~in_sfr(l, b) & (a0 < 0.15) & (np.abs(b) > 20)
              & (f["chi2_fit"] < 11.07) & ~edge
              & (np.abs(f["cstar_nsigma"]) < 1.0) & (f["ruwe"] < 1.1)
              & (f["dchi2"] > 25))
    if low:
        common = common & (f["alpha_upper_2sig"] < 1.5)
    else:
        common = common & (f["alpha_lower_2sig"] > 2.6)
    dim = common & (f["deficit_G_mag"] > 0.10) & (f["deficit_G_mag"] < 3.0)
    bright = common & (f["deficit_G_mag"] < -0.10) & (f["deficit_G_mag"] > -3.0)
    return int(np.nansum(dim)), int(np.nansum(bright))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--ninject", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=20260821)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    findings, out = [], {"tag": args.tag, "seed": args.seed}

    cols = ["source_id", "ra", "dec", "l", "b", "A_0", "bp_rp", "dist_mod",
            "phot_bp_mean_mag", "phot_g_mean_mag", "phot_rp_mean_mag",
            "tmass_j_m", "tmass_h_m", "tmass_ks_m",
            "wise_w1mpro", "wise_w1mpro_error", "wise_w2mpro",
            "wise_w2mpro_error", "cstar_nsigma", "ruwe"]
    print("loading residual sample ...", flush=True)
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=cols)
    print(f"  {len(d):,} rows", flush=True)

    R, sig, ok = build_residual_matrix(d)
    sel = np.flatnonzero(ok)
    print(f"  {len(sel):,} fitted; band sigma "
          + " ".join(f"{b}={s:.4f}" for b, s in zip(FITBANDS, sig)), flush=True)
    basis = build_basis()

    # ------------------------------------------------------------------ T0
    print("\n=== T0  reproduction against the stored channel-5 output ===",
          flush=True)
    fit = fit_grid(R, sel, sig, basis, want_lo=True)
    fit["cstar_nsigma"] = d["cstar_nsigma"].to_numpy(float)[sel]
    fit["ruwe"] = d["ruwe"].to_numpy(float)[sel]
    sid = d["source_id"].to_numpy()[sel]
    L = d["l"].to_numpy(float)[sel]
    B = d["b"].to_numpy(float)[sel]
    A0 = d["A_0"].to_numpy(float)[sel]

    stored_path = cfg.DERIVED_DIR / f"spectral_slope_{args.tag}.parquet"
    if not stored_path.exists():
        # A smoke run on the testrun tag has no stored channel-5 output. Say so
        # rather than silently skipping: an audit whose reproduction check did
        # not run has not established that it is auditing the pipeline.
        print(f"  {stored_path.name} absent -- reproduction check NOT RUN")
        out["T0_reproduction"] = {"ran": False, "why": "no stored output for "
                                  "this tag; T1-T3 are unvalidated here"}
        findings.append("T0 DID NOT RUN (no stored channel-5 parquet for this "
                        "tag). Nothing below is validated against the pipeline.")
        stored = None
    else:
        stored = pd.read_parquet(
            stored_path,
            columns=["source_id", "alpha", "alpha_upper_2sig", "deficit_G_mag",
                     "dchi2", "chi2_fit"])
    mine = pd.DataFrame({"source_id": sid, "alpha_m": fit["alpha"],
                         "ahi_m": fit["alpha_upper_2sig"],
                         "A_m": fit["deficit_G_mag"],
                         "dchi2_m": fit["dchi2"]})
    cmp = stored.merge(mine, on="source_id", how="inner") if stored is not None \
        else mine.iloc[:0].assign(alpha=[], deficit_G_mag=[], dchi2=[],
                                  alpha_upper_2sig=[])
    t0 = {"ran": stored is not None,
          "n_stored": int(len(stored)) if stored is not None else 0,
          "n_refit": int(len(sel)), "n_matched": int(len(cmp))}
    for a, b_, nm in (("alpha", "alpha_m", "alpha"),
                      ("deficit_G_mag", "A_m", "deficit_G_mag"),
                      ("dchi2", "dchi2_m", "dchi2"),
                      ("alpha_upper_2sig", "ahi_m", "alpha_upper_2sig")):
        dv = cmp[b_].to_numpy(float) - cmp[a].to_numpy(float)
        dv = dv[np.isfinite(dv)]
        t0[f"max_abs_diff_{nm}"] = float(np.max(np.abs(dv))) if len(dv) else None
        t0[f"frac_exact_{nm}"] = float(np.mean(np.abs(dv) < 1e-4)) if len(dv) else None
        if len(dv):
            print(f"  {nm:18s} max|diff| = {t0[f'max_abs_diff_{nm}']:.3e}   "
                  f"frac identical = {t0[f'frac_exact_{nm}']:.5f}")
    if t0["ran"]:
        t0["reproduced"] = bool(t0["frac_exact_alpha"] is not None
                                and t0["frac_exact_alpha"] > 0.999
                                and t0["max_abs_diff_deficit_G_mag"] < 1e-3)
        print(f"  REPRODUCED: {t0['reproduced']}")
        if not t0["reproduced"]:
            findings.append(
                "T0 FAILED: the transcribed fitter does not reproduce the "
                "stored channel-5 output, so T1-T3 below describe a different "
                "estimator and must not be cited.")
        out["T0_reproduction"] = t0

    # ------------------------------------------------------------------ T1
    print("\n=== T1  injection-recovery of alpha ===", flush=True)
    pool = rng.choice(len(sel), size=min(args.ninject, len(sel)), replace=False)
    sub = sel[pool]
    Rsub = R[sub].copy()
    cn = d["cstar_nsigma"].to_numpy(float)[sub]
    rw = d["ruwe"].to_numpy(float)[sub]
    quiet = fit["dchi2"][pool] < 9.0   # no significant absorber before injection
    print(f"  canvas {len(sub):,} random fitted stars "
          f"({quiet.mean()*100:.1f}% carry no significant absorber already)")

    grid_alpha = [-0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    grid_amp = [0.05, 0.10, 0.20, 0.50, 1.00]
    rows = []
    ident = np.arange(len(sub))
    for at in grid_alpha:
        ia = int(np.argmin(np.abs(ALPHA_GRID - at)))
        vec = basis[ia]
        for A in grid_amp:
            Rinj = Rsub + A * vec[None, :]
            f = fit_grid(Rinj, ident, sig, basis, want_lo=True)
            sigmask = ((f["dchi2"] > 25) & (f["deficit_G_mag"] > 0.10)
                       & (f["deficit_G_mag"] < 3.0))
            clean = (sigmask & (np.abs(cn) < 1.0) & (rw < 1.1)
                     & (f["chi2_fit"] < 14.0))
            flagged_low = clean & (f["alpha_upper_2sig"] < 1.5)
            flagged_hi = clean & (f["alpha_lower_2sig"] > 2.6)
            q = quiet
            rows.append({
                "alpha_true": float(ALPHA_GRID[ia]), "amp_true": A,
                "median_alpha_rec": float(np.median(f["alpha"][q])),
                "median_amp_rec": float(np.median(f["deficit_G_mag"][q])),
                "frac_significant": float(sigmask[q].mean()),
                "frac_flagged_flat": float(flagged_low[q].mean()),
                "frac_flagged_steep": float(flagged_hi[q].mean()),
                "n_quiet": int(q.sum()),
                "n_flagged_flat": int(flagged_low[q].sum()),
                "n_flagged_steep": int(flagged_hi[q].sum()),
            })
            r = rows[-1]
            print(f"  a_true={r['alpha_true']:+5.2f} A={A:4.2f}  ->  "
                  f"a_rec(med)={r['median_alpha_rec']:+6.3f}  "
                  f"A_rec(med)={r['median_amp_rec']:+6.3f}  "
                  f"sig={r['frac_significant']:.3f}  "
                  f"flat-flag={r['frac_flagged_flat']:.4f}  "
                  f"steep-flag={r['frac_flagged_steep']:.4f}", flush=True)
            del Rinj, f
    out["T1_injection"] = rows

    dust = sorted([r for r in rows if abs(r["alpha_true"] - 2.0) < 0.06],
                  key=lambda r: r["amp_true"])
    flat = sorted([r for r in rows if abs(r["alpha_true"] - 0.5) < 0.06],
                  key=lambda r: r["amp_true"])
    # Below what amplitude does the fitter stop measuring alpha at all? Read it
    # off the injection table rather than asserting it.
    floor = {}
    for at in grid_alpha:
        rs = sorted([r for r in rows if abs(r["alpha_true"] - at) < 0.06],
                    key=lambda r: r["amp_true"])
        good = [r["amp_true"] for r in rs
                if abs(r["median_alpha_rec"] - r["alpha_true"]) < 0.25]
        floor[f'{rs[0]["alpha_true"]:.2f}'] = min(good) if good else None
    print("\n  smallest injected amplitude at which alpha is recovered "
          "to +/-0.25 mag:")
    for k, v in floor.items():
        print(f"    alpha_true={k:>6s}  ->  "
              f"{f'{v:.2f} mag' if v else 'never, within the injected grid'}")
    out["T1_alpha_recovery_floor_mag"] = floor

    out["T1_summary"] = {
        "type_s_rate_dust_called_flat": {str(r["amp_true"]): r["frac_flagged_flat"]
                                         for r in dust},
        "power_flat_called_flat": {str(r["amp_true"]): r["frac_flagged_flat"]
                                   for r in flat},
        "alpha_bias_at_dust": {str(r["amp_true"]):
                               r["median_alpha_rec"] - r["alpha_true"]
                               for r in dust},
    }

    # --- T1b. Fold the measured type-S rate through the REAL amplitude
    # distribution. "5,239 stars have alpha significantly below 1.5" only means
    # something if it exceeds what pure dust plus this fitter's noise would
    # produce at the amplitudes those stars actually have. The comparison has
    # never been made because the type-S rate was never measured.
    print("\n=== T1b  what pure dust would produce at the REAL amplitudes ===",
          flush=True)
    real_sig = ((fit["dchi2"] > 25) & (fit["deficit_G_mag"] > 0.10)
                & (fit["deficit_G_mag"] < 3.0))
    amp_real = fit["deficit_G_mag"][real_sig].astype(float)
    xa = np.log10([r["amp_true"] for r in dust])
    ya = np.array([r["frac_flagged_flat"] for r in dust])
    # clamp outside the injected range rather than extrapolate a rate
    rate = np.interp(np.log10(np.clip(amp_real, 10 ** xa[0], 10 ** xa[-1])),
                     xa, ya)
    exp_flat = float(rate.sum())
    # the same fold for the ANOMALOUS cut stack, which is what 5,239 counts
    real_low = real_sig & (fit["alpha_upper_2sig"] < 1.5) \
        & (np.abs(fit["cstar_nsigma"]) < 1.0) & (fit["ruwe"] < 1.1)
    n_low = int(np.nansum(real_low))
    print(f"  real significant absorbers            : {int(real_sig.sum()):,}")
    print(f"  median amplitude                      : "
          f"{np.median(amp_real):.3f} mag")
    print(f"  observed flagged flat (clean phot/ast): {n_low:,}")
    print(f"  expected if ALL of them were pure dust: {exp_flat:.0f}")
    excess = n_low / max(exp_flat, 1e-9)
    print(f"  observed / dust-expectation           : {excess:.2f}")
    out["T1b_dust_expectation"] = {
        "n_significant_real": int(real_sig.sum()),
        "median_amplitude_real": float(np.median(amp_real)),
        "n_flagged_flat_observed": n_low,
        "n_flagged_flat_expected_if_pure_dust": exp_flat,
        "observed_over_dust_expectation": float(excess),
    }
    if excess < 2.0:
        findings.append(
            f"T1b: the flat-slope population ({n_low:,}) is within a factor "
            f"{excess:.2f} of what an ALL-DUST sample would produce through "
            f"this fitter ({exp_flat:.0f}). The 'flat absorber' population is "
            f"consistent with misclassified dust and carries no independent "
            f"information.")

    # --- T1c. Is the steep tail a FAIR mirror for the flat tail? Compare the
    # two misclassification rates on the SAME injected dust population. If the
    # estimator scatters alpha downward more readily than upward, a small steep
    # count does not certify a large flat count -- the same one-signed-
    # contaminant logic already applied to the 2MASS aperture mismatch.
    r05 = next(r for r in dust if abs(r["amp_true"] - 0.5) < 1e-9)
    skew = r05["frac_flagged_flat"] / max(r05["frac_flagged_steep"], 1e-9)
    print(f"\n  mirror fairness at alpha_true=2, A=0.5 mag: "
          f"flat-flag {r05['frac_flagged_flat']:.5f} vs steep-flag "
          f"{r05['frac_flagged_steep']:.5f}  (ratio {skew:.1f})")
    out["T1c_mirror_fairness"] = {
        "flat_flag_rate_at_dust": r05["frac_flagged_flat"],
        "steep_flag_rate_at_dust": r05["frac_flagged_steep"],
        "flat_over_steep_misclassification": float(skew),
    }
    if skew > 3.0:
        findings.append(
            f"T1c: on injected PURE DUST the fitter calls a star significantly "
            f"flat at least {skew:.0f}x more often than significantly steep "
            f"({r05['n_flagged_flat']} vs {r05['n_flagged_steep']} of "
            f"{r05['n_quiet']} injected). The alpha axis is skewed, so "
            f"the steep tail is NOT a valid mirror for the flat tail and a "
            f"small steep count certifies nothing. Use the injection-measured "
            f"type-S rate (T1b) instead.")

    # ------------------------------------------------------------------ T2
    print("\n=== T2  the slope mirror channel 5 never ran "
          "(alpha significantly > 2.6) ===", flush=True)
    n_flat_dim, n_flat_bright = survivor_cuts(fit, L, B, A0, low=True)
    n_steep_dim, n_steep_bright = survivor_cuts(fit, L, B, A0, low=False)
    print(f"  FLAT  (alpha_2sig_upper < 1.5) survivors : "
          f"{n_flat_dim:6,} dimmed   {n_flat_bright:5,} brightened")
    print(f"  STEEP (alpha_2sig_lower > 2.6) survivors : "
          f"{n_steep_dim:6,} dimmed   {n_steep_bright:5,} brightened")
    ratio = n_flat_dim / max(n_steep_dim, 1)
    print(f"  flat/steep = {ratio:.2f}")
    out["T2_slope_mirror"] = {
        "n_flat_dimmed": n_flat_dim, "n_flat_brightened": n_flat_bright,
        "n_steep_dimmed": n_steep_dim, "n_steep_brightened": n_steep_bright,
        "flat_over_steep": float(ratio),
        "published_flat_dimmed": 372,
    }
    if n_steep_dim >= 0.5 * n_flat_dim:
        findings.append(
            f"T2: the unphysical STEEP tail returns {n_steep_dim:,} survivors "
            f"against {n_flat_dim:,} flat ones under identical cuts. Grain "
            f"physics produces neither, so the flat set is dominated by the "
            f"fitter's own slope-axis noise and 'alpha away from dust' carries "
            f"little specificity.")

    # ------------------------------------------------------------------ T3
    print("\n=== T3A  mark permutation of whole residual vectors "
          "(fit fixed, cuts move) ===", flush=True)
    # 372 has never been compared to a null. The amplitude mirror in scripts/44
    # (372 vs 13) only asks whether something dims; it does not ask whether the
    # survivors sit on clean sky more often than chance. Shuffle the fitted SED
    # vectors -- with the photometric and astrometric quality marks that travel
    # with the star -- over sky positions and extinctions. Values over
    # positions, never synthetic positions.
    n_perm = 200
    counts = np.empty(n_perm, dtype=int)
    for i in range(n_perm):
        perm = rng.permutation(len(sel))
        fitA = {k: v[perm] for k, v in fit.items()}
        counts[i], _ = survivor_cuts(fitA, L, B, A0, low=True)
    mu_p, sd_p = float(counts.mean()), float(counts.std(ddof=1))
    z = (n_flat_dim - mu_p) / max(sd_p, 1e-9)
    # analytic cross-check: the permutation is estimating
    # N(passes fit-side cuts) x fraction(passes sky-side cuts)
    edge_a = (fit["alpha"] < -0.9) | (fit["alpha"] > 5.9)
    fitside = ((fit["chi2_fit"] < 11.07) & ~edge_a
               & (np.abs(fit["cstar_nsigma"]) < 1.0) & (fit["ruwe"] < 1.1)
               & (fit["dchi2"] > 25) & (fit["alpha_upper_2sig"] < 1.5)
               & (fit["deficit_G_mag"] > 0.10) & (fit["deficit_G_mag"] < 3.0))
    skyside = ~in_sfr(L, B) & (A0 < 0.15) & (np.abs(B) > 20)
    analytic = float(np.nansum(fitside)) * float(np.mean(skyside))
    print(f"  real survivors                     : {n_flat_dim:,}")
    print(f"  mark-permutation null ({n_perm} draws)  : "
          f"{mu_p:.1f} +/- {sd_p:.1f}")
    print(f"  analytic factorisation cross-check : {analytic:.1f}")
    print(f"  z = {z:+.1f}   real/null = {n_flat_dim / max(mu_p, 1e-9):.2f}")
    out["T3A_vector_permutation"] = {
        "n_real": n_flat_dim, "n_permutations": n_perm,
        "null_mean": mu_p, "null_sd": sd_p,
        "analytic_factorisation": analytic,
        "z": float(z), "real_over_null": float(n_flat_dim / max(mu_p, 1e-9)),
    }
    if z < -3:
        findings.append(
            f"T3A: against a mark-permutation null the survivor set is a "
            f"DEFICIT, not an excess: {n_flat_dim:,} observed against "
            f"{mu_p:.0f} +/- {sd_p:.0f} expected (z = {z:+.1f}). The flat-slope "
            f"population avoids clean, high-latitude, low-extinction sky, "
            f"which is what dust does. The published framing -- 372 survivors, "
            f"asymmetric at 28.6:1 over the amplitude mirror -- was never "
            f"compared to this null, and the amplitude mirror does not test "
            f"sky association at all.")

    print("\n=== T3B  band-wise permutation (SED coherence destroyed) ===",
          flush=True)
    Rperm = R.copy()
    for j in range(Rperm.shape[1]):
        Rperm[sel, j] = Rperm[sel[rng.permutation(len(sel))], j]
    fitB = fit_grid(Rperm, sel, sig, basis, want_lo=True)
    fitB["cstar_nsigma"] = fit["cstar_nsigma"]
    fitB["ruwe"] = fit["ruwe"]
    sigB = ((fitB["dchi2"] > 25) & (fitB["deficit_G_mag"] > 0.10)
            & (fitB["deficit_G_mag"] < 3.0))
    sigR = ((fit["dchi2"] > 25) & (fit["deficit_G_mag"] > 0.10)
            & (fit["deficit_G_mag"] < 3.0))
    nB_dim, nB_bright = survivor_cuts(fitB, L, B, A0, low=True)
    nB_steep, _ = survivor_cuts(fitB, L, B, A0, low=False)
    print(f"  'significant absorbers' real {int(sigR.sum()):,}  "
          f"vs band-permuted {int(sigB.sum()):,}")
    print(f"  median alpha real {np.median(fit['alpha'][sigR]):.3f}  "
          f"vs band-permuted {np.median(fitB['alpha'][sigB]):.3f}")
    print(f"  flat survivors  real {n_flat_dim:,}  vs band-permuted {nB_dim:,}")
    print(f"  steep survivors real {n_steep_dim:,}  vs band-permuted {nB_steep:,}")
    out["T3B_bandwise_permutation"] = {
        "n_significant_real": int(sigR.sum()),
        "n_significant_permuted": int(sigB.sum()),
        "median_alpha_real": float(np.median(fit["alpha"][sigR])),
        "median_alpha_permuted": float(np.median(fitB["alpha"][sigB])),
        "n_flat_survivors_real": n_flat_dim,
        "n_flat_survivors_permuted": nB_dim,
        "n_steep_survivors_permuted": nB_steep,
    }
    if sigB.sum() > 0.5 * sigR.sum():
        findings.append(
            f"T3B: incoherent band-wise noise alone manufactures "
            f"{int(sigB.sum()):,} 'significant absorbers' against "
            f"{int(sigR.sum()):,} real ones. The significance threshold "
            f"dchi2>25 is not selecting coherent SED shape.")

    out["findings"] = findings
    (cfg.RESULT_DIR / f"audit_spectral_slope_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    print("\n=== FINDINGS ===")
    if findings:
        for f_ in findings:
            print("  * " + f_)
    else:
        print("  none of the pre-registered failure conditions fired")
    print(f"\nwrote results/audit_spectral_slope_{args.tag}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
