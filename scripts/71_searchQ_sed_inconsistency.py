#!/usr/bin/env python
"""Search Q: does the SED disagree with the spectrum?

    run.sh scripts/71_searchQ_sed_inconsistency.py --tag primary

WHY THIS IS INDEPENDENT OF EVERYTHING ELSE HERE
-----------------------------------------------
Every optical channel in this project measures absolute G against absolute Ks.
That makes them all share two weaknesses: they are anchored on 2MASS, whose
4 arcsec beam is the dominant contaminant, and they need a distance, so the
parallax and the extinction map enter every number.

This channel needs none of that. Gaia derives stellar parameters twice, by
methods that fail differently:

    GSP-Phot   fits the low-resolution BP/RP SED, i.e. the SHAPE OF THE LIGHT.
    GSP-Spec   fits the RVS spectrum, i.e. the ATOMIC LINES.

An absorber sitting in front of a star changes the first and not the second.
Intervening material reddens and dims the energy distribution, so a
photometric fit reads a cooler star; the photosphere's line ratios are
untouched, so the spectroscopic fit reads the true temperature. The signature
is therefore

    Teff(spec) > Teff(phot)

with no distance, no 2MASS, and no dust map involved anywhere.

THE CONTAMINANT THAT DOES THE SAME THING
----------------------------------------
Interstellar dust. It reddens the SED exactly as an absorber would, which is
precisely why GSP-Phot fits extinction as a free parameter -- but it fits it
imperfectly, and any residual error lands in Teff. So the test is repeated at
low extinction, where there is almost nothing for the fit to get wrong.

THE MIRROR CONTROL
------------------
Absorbing material can only make the photometric temperature too COOL. A star
whose photometric temperature is too HOT has no absorber explanation and can
only be fitting error, so the negative tail measures this estimator's own
false-positive rate directly.

THE CROSS-CHECK THAT WOULD MATTER
---------------------------------
If a star is genuinely veiled, the two independent estimators must agree about
it: it should show BOTH a spectroscopic-minus-photometric temperature excess
AND an optical deficit in M_G at fixed M_Ks. Neither estimator shares a
failure mode with the other, so their overlap is far more informative than
either tail alone.
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
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchQ")

COLS = ["source_id", "l", "b", "A_0", "dist_pc", "residual", "bp_rp",
        "phot_g_mean_mag", "teff_gspphot", "logg_gspphot", "mh_gspphot",
        "teff_gspspec", "logg_gspspec", "mh_gspspec", "flags_gspspec",
        "ruwe", "cstar_nsigma", "M_Ks"]

K_SIGMA = 3.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=K_SIGMA)
    args = ap.parse_args()

    import pyarrow.parquet as pq
    path = cfg.DERIVED_DIR / f"{args.tag}_resid.parquet"
    have = set(pq.read_schema(path).names)
    use = [c for c in COLS if c in have]
    missing = [c for c in COLS if c not in have]
    if missing:
        log.warning("absent from the derived table: %s", ", ".join(missing))
    d = pd.read_parquet(path, columns=use)
    log.info("loaded %d stars", len(d))

    if "teff_gspspec" not in d.columns:
        log.error("teff_gspspec is not in the derived table; "
                  "this channel needs the GSP-Spec columns")
        return 1

    tp = d["teff_gspphot"].to_numpy(float)
    ts = d["teff_gspspec"].to_numpy(float)
    both = np.isfinite(tp) & np.isfinite(ts) & (tp > 2500) & (tp < 8000) \
        & (ts > 2500) & (ts < 8000)
    log.info("with BOTH photometric and spectroscopic Teff: %d",
             int(both.sum()))
    if both.sum() < 1000:
        log.error("too few stars have both temperatures")
        return 1

    # ---- GSP-Spec quality ------------------------------------------------
    # flags_gspspec is a 41-character quality string; the leading characters
    # encode the reliability of the parameter fit, where '0' is best.
    if "flags_gspspec" in d.columns:
        fl = d["flags_gspspec"].astype(str)
        good_spec = fl.str[:3].isin(["000", "001", "010", "100",
                                     "011", "101", "110", "111"]).to_numpy().copy()
        good_spec &= (fl.str.len() > 3).to_numpy()
        log.info("with usable GSP-Spec flags: %d", int((both & good_spec).sum()))
    else:
        good_spec = np.ones(len(d), dtype=bool)

    clean = (both & good_spec
             & (d["ruwe"].to_numpy(float) < 1.4)
             & (np.abs(d["cstar_nsigma"].to_numpy(float)) < 3.0))
    log.info("clean sample for this channel: %d", int(clean.sum()))

    dteff = ts - tp
    # Work in a relative measure so the scatter is comparable across the
    # temperature range rather than dominated by the hot end.
    rel = dteff / np.maximum(tp, 1.0)

    r_clean = rel[clean]
    med = float(np.median(r_clean))
    sig = float(st.robust_sigma(r_clean))
    log.info("Teff(spec)/Teff(phot) - 1 : median %+.4f, robust sigma %.4f",
             med, sig)

    thr = args.k * sig
    veiled = clean & (rel > med + thr)      # photometric Teff too COOL
    inverse = clean & (rel < med - thr)     # too HOT: no absorber explanation
    n_v, n_i = int(veiled.sum()), int(inverse.sum())
    ratio = n_v / n_i if n_i else np.inf
    log.info("")
    log.info("SED cooler than the spectrum (absorber-like) : %d", n_v)
    log.info("SED hotter than the spectrum (mirror)        : %d", n_i)
    log.info("ratio                                        : %s",
             f"{ratio:.2f}" if np.isfinite(ratio) else "inf")

    # ---- extinction control ----------------------------------------------
    a0 = d["A_0"].to_numpy(float)
    log.info("")
    log.info("=== extinction control ===")
    ext_rows = []
    for lo, hi in [(0.0, 0.02), (0.02, 0.05), (0.05, 0.10),
                   (0.10, 0.20), (0.20, 1.0)]:
        m = clean & (a0 >= lo) & (a0 < hi)
        if m.sum() < 500:
            continue
        rr = rel[m]
        mm = float(np.median(rr))
        nv = int(np.count_nonzero(rr > med + thr))
        ni = int(np.count_nonzero(rr < med - thr))
        ext_rows.append({"a0_lo": lo, "a0_hi": hi, "n": int(m.sum()),
                         "median_rel": mm, "veiled": nv, "inverse": ni,
                         "ratio": (nv / ni) if ni else None})
        log.info("  A_0 %.2f-%-5.2f n=%7d  median=%+.4f  veiled=%5d "
                 "inverse=%5d  ratio=%s", lo, hi, int(m.sum()), mm, nv, ni,
                 f"{nv/ni:6.2f}" if ni else "   inf")

    # ---- the cross-check with the optical deficit ------------------------
    log.info("")
    log.info("=== agreement with the independent M_G deficit ===")
    resid = d["residual"].to_numpy(float)
    rs = st.robust_sigma(resid[np.isfinite(resid)])
    dim = np.isfinite(resid) & (resid > args.k * rs)
    log.info("optically dim (> %.0f sigma of %.4f mag): %d",
             args.k, rs, int((clean & dim).sum()))

    n_both = int((veiled & dim).sum())
    p_v = veiled[clean].mean()
    p_d = dim[clean].mean()
    expected = p_v * p_d * clean.sum()
    log.info("veiled AND optically dim: %d observed, %.1f expected if "
             "independent  (enrichment %.2fx)",
             n_both, expected, n_both / max(expected, 1e-9))

    # mirror of the joint test: the two estimators agreeing in the
    # unphysical direction is the joint false-positive rate
    n_both_mirror = int((inverse & np.isfinite(resid)
                         & (resid < -args.k * rs)).sum())
    log.info("inverse AND optically bright (joint mirror): %d", n_both_mirror)

    # ---- verdict ----------------------------------------------------------
    low_ext = [r for r in ext_rows if r["a0_hi"] <= 0.05 and r["ratio"]]
    low_ratio = np.mean([r["ratio"] for r in low_ext]) if low_ext else np.nan

    if n_v <= n_i:
        verdict = (
            f"NULL. The photometric temperature is no more often too cool than "
            f"too hot ({n_v} vs {n_i}). The BP/RP energy distribution agrees "
            f"with the RVS line spectrum, so no population carries an SED "
            f"distorted in the way intervening material would distort it. This "
            f"is the one optical result here that uses neither 2MASS, nor a "
            f"distance, nor the dust map.")
    elif np.isfinite(low_ratio) and low_ratio < 1.5:
        verdict = (
            f"EXTINCTION. The absorber-like tail exceeds its mirror overall "
            f"({ratio:.2f}) but the excess disappears at low extinction "
            f"(ratio {low_ratio:.2f} for A_0 < 0.05). The asymmetry is residual "
            f"error in the extinction that GSP-Phot fits alongside "
            f"temperature, not intervening material.")
    elif n_both > 3 * expected:
        verdict = (
            f"TWO INDEPENDENT ESTIMATORS AGREE. {n_v} stars have a "
            f"photometric temperature too cool for their spectrum against "
            f"{n_i} in the unphysical mirror, the excess survives low "
            f"extinction (ratio {low_ratio:.2f}), and {n_both} of them are "
            f"ALSO optically dim in M_G at fixed M_Ks -- {n_both/max(expected,1e-9):.1f}x "
            f"the chance overlap, against {n_both_mirror} in the joint mirror. "
            f"These two estimators share no failure mode: one uses BP/RP shape "
            f"versus RVS lines, the other uses G versus 2MASS Ks with a "
            f"distance and a dust map. Follow-up is warranted.")
    else:
        verdict = (
            f"ASYMMETRIC but unconfirmed. The absorber-like tail exceeds its "
            f"mirror ({n_v} vs {n_i}, ratio {ratio:.2f}) and survives low "
            f"extinction, but the overlap with the independent optical deficit "
            f"is {n_both} against {expected:.1f} expected by chance, i.e. no "
            f"real enrichment. Two estimators that should both respond to a "
            f"veiled star do not pick out the same stars, which argues the "
            f"asymmetry is a property of the GSP-Phot fit rather than of the "
            f"stars.")

    print(f"\n{'='*74}")
    print("SEARCH Q: DOES THE ENERGY DISTRIBUTION DISAGREE WITH THE SPECTRUM?")
    print(f"{'='*74}")
    print(f"  stars with both Teff estimates      : {int(both.sum()):,}")
    print(f"  clean sample                        : {int(clean.sum()):,}")
    print(f"  relative Teff scatter               : {sig:.4f}")
    print()
    print(f"  SED too COOL (absorber-like)        : {n_v}")
    print(f"  SED too HOT  (mirror, unphysical)   : {n_i}")
    print(f"  ratio                               : "
          f"{ratio:.2f}" if np.isfinite(ratio) else "  ratio: inf")
    print(f"  ratio at A_0 < 0.05                 : "
          f"{low_ratio:.2f}" if np.isfinite(low_ratio) else "  n/a")
    print()
    print(f"  veiled AND optically dim            : {n_both} "
          f"(expected {expected:.1f})")
    print(f"  joint mirror                        : {n_both_mirror}")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "n_both_teff": int(both.sum()),
        "n_clean": int(clean.sum()),
        "median_rel_dteff": med,
        "robust_sigma_rel": sig,
        "k_sigma": args.k,
        "n_veiled": n_v,
        "n_inverse_mirror": n_i,
        "ratio": float(ratio) if np.isfinite(ratio) else None,
        "extinction_bins": ext_rows,
        "low_extinction_ratio": (
            float(low_ratio) if np.isfinite(low_ratio) else None),
        "n_veiled_and_dim": n_both,
        "n_expected_by_chance": float(expected),
        "n_joint_mirror": n_both_mirror,
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchQ_sed_inconsistency_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)

    if n_both:
        cand = d[veiled & dim].copy()
        cand["rel_dteff"] = rel[veiled & dim]
        cand.sort_values("rel_dteff", ascending=False).to_csv(
            cfg.RESULT_DIR / f"searchQ_candidates_{args.tag}.csv", index=False)
        log.info("wrote %d joint candidates", n_both)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
