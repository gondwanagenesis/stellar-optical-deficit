#!/usr/bin/env python
"""Is the discard-pile deficit aperture mismatch, or something else?

    run.sh scripts/67_searchM_blend_test.py --tag primary

THE OPEN QUESTION FROM SEARCH M
-------------------------------
Search M found every auditable discard class to be one-sidedly dim far beyond
the retained sample: variables 180:1, C* failures 98:1, extragalactic 60:1,
poor-Ks 70:1, against a reference of 5:1. All survived a low-extinction cut.

That looked like a signal and it cannot be read as one, because the test's own
premise was wrong. Search M assumed contamination scatters a residual both
ways. This project's dominant contaminant does not: 2MASS resolves a 4 arcsec
beam where Gaia resolves sub-arcsecond, so a neighbour inside the beam adds
flux to Ks and not to G. That makes the star look intrinsically bright in Ks,
which the estimator reads as a G-band deficit. The bias is ONE-SIGNED by
construction, and three of the four flagged classes select for exactly the
conditions that produce it.

So one-sidedness cannot discriminate. This can.

THE TEST
--------
2MASS publishes, per cross-match, how many Gaia sources fell in the beam
(xm_nnb), how many other 2MASS entries competed for the match (xm_nmates), and
how far the match sat from the Gaia position (xm_dist). Aperture mismatch is
a direct, monotonic function of those three: it is largest when the beam is
crowded and the match is displaced, and it vanishes when a single 2MASS source
sits precisely on a single Gaia source.

If the dim excess is aperture mismatch, it must weaken sharply in the
blend-clean subsample and strengthen in the blend-suspect one. If it is
present at the same strength in stars with a single, well-centred, uncontested
2MASS match, then the beam is not the explanation and the class needs another.
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
from pipeline import sample as smp
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("blendtest")

COLS = [
    "source_id", "l", "b", "parallax", "parallax_error",
    "phot_g_mean_mag", "bp_rp", "phot_bp_rp_excess_factor",
    "phot_variable_flag", "classprob_dsc_combmod_star",
    "in_qso_candidates", "in_galaxy_candidates",
    "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat",
    "astrometric_params_solved",
    "tmass_ks_m", "tmass_ph_qual",
    "tmass_xm_dist", "tmass_xm_nnb", "tmass_xm_nmates",
]

K_SIGMA = 3.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=K_SIGMA)
    args = ap.parse_args()

    import pyarrow.parquet as pq
    files = sorted(cfg.RAW_DIR.glob("sample_d500_p*.parquet"))
    frames = []
    for f in files:
        have = set(pq.read_schema(f).names)
        frames.append(pd.read_parquet(f, columns=[c for c in COLS if c in have]))
    d = pd.concat(frames, ignore_index=True)
    log.info("raw rows: %d", len(d))

    # ---- geometry, extinction, absolute magnitudes -----------------------
    d = smp.add_astrometry(d)
    a0 = ext.query_a0("edenhofer23", d["l"].to_numpy(float),
                      d["b"].to_numpy(float), d["dist_pc"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = d["dist_mod"].to_numpy(float)
    M_G = d["phot_g_mean_mag"].to_numpy(float) - mu - a_g
    M_Ks = d["tmass_ks_m"].to_numpy(float) - mu - a_ks

    cstar = smp.corrected_excess_factor(
        bp_rp, d["phot_bp_rp_excess_factor"].to_numpy(float))
    cstar_n = cstar / np.maximum(
        smp.excess_factor_sigma(d["phot_g_mean_mag"].to_numpy(float)), 1e-9)

    is_var = (d["phot_variable_flag"].astype(str).str.upper() == "VARIABLE").to_numpy()
    dsc = d["classprob_dsc_combmod_star"].fillna(1.0).to_numpy()
    qso = d["in_qso_candidates"].fillna(False).to_numpy(bool)
    gal = d["in_galaxy_candidates"].fillna(False).to_numpy(bool)
    ks_q = d["tmass_ph_qual"].astype(str)

    base = (np.isfinite(M_G) & np.isfinite(M_Ks) & np.isfinite(a0)
            & (M_Ks > 3.0) & (M_Ks < 8.0)
            & (bp_rp > 0.7) & (bp_rp < 3.6)
            & (d["dist_pc"].to_numpy() > 10) & (d["dist_pc"].to_numpy() < 500)
            & (a_g < 0.5))

    clean = (base & ~is_var & ~qso & ~gal & (dsc > 0.5)
             & (np.abs(cstar_n) < 3.0)
             & (ks_q.str[2] == "A").to_numpy())

    x, y = M_Ks[clean], M_G[clean]
    best = None
    for deg in (3, 4, 5):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, sigma, deg = best
    resid = M_G - np.polyval(coef, M_Ks)
    thr = args.k * sigma
    log.info("fiducial degree %d sigma %.4f; threshold %.3f mag", deg, sigma, thr)

    # ---- the blend axis ---------------------------------------------------
    nnb = d["tmass_xm_nnb"].fillna(1).to_numpy()
    nmates = d["tmass_xm_nmates"].fillna(1).to_numpy()
    xmd = d["tmass_xm_dist"].fillna(0.0).to_numpy()

    blend_clean = (nnb <= 1) & (nmates <= 1) & (xmd < 0.3)
    blend_bad = (nnb > 1) | (nmates > 1) | (xmd > 1.0)
    log.info("blend-clean: %d   blend-suspect: %d",
             int(blend_clean.sum()), int(blend_bad.sum()))

    classes = {
        "RETAINED (reference)": clean,
        "variable": base & is_var,
        "C* beyond 3 sigma": base & (np.abs(cstar_n) >= 3.0),
        "classified extragalactic": base & ((dsc <= 0.5) | qso | gal),
        "2MASS Ks not 'A'": base & (ks_q.str[2] != "A").to_numpy(),
    }

    def ratio_of(mask):
        r = resid[mask & np.isfinite(resid)]
        if len(r) < 100:
            return None
        med = float(np.median(r))
        nd = int(np.count_nonzero(r > med + thr))
        nb = int(np.count_nonzero(r < med - thr))
        return {"n": len(r), "median": med, "dim": nd, "bright": nb,
                "ratio": (nd / nb) if nb else None}

    log.info("")
    log.info("%-28s %-22s %-22s", "class", "blend-CLEAN", "blend-SUSPECT")
    rows = []
    for name, m in classes.items():
        # For the C* class the C* indicator is the selector, so only the
        # 2MASS beam indicators are used to define cleanliness there.
        c_ok = ratio_of(m & blend_clean)
        c_bad = ratio_of(m & blend_bad)
        f_clean = c_ok["ratio"] if c_ok and c_ok["ratio"] else np.nan
        f_bad = c_bad["ratio"] if c_bad and c_bad["ratio"] else np.nan
        rows.append({"class": name, "clean": c_ok, "suspect": c_bad})
        log.info("%-28s n=%-7s r=%-8s  n=%-7s r=%-8s",
                 name,
                 f"{c_ok['n']}" if c_ok else "-",
                 f"{f_clean:.2f}" if np.isfinite(f_clean) else "-",
                 f"{c_bad['n']}" if c_bad else "-",
                 f"{f_bad:.2f}" if np.isfinite(f_bad) else "-")

    # ---- does the deficit scale with beam crowding? ----------------------
    log.info("")
    log.info("=== deficit vs 2MASS cross-match distance (all sources) ===")
    dist_rows = []
    for lo, hi in [(0.0, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0),
                   (1.0, 2.0), (2.0, 100.0)]:
        m = base & (xmd >= lo) & (xmd < hi) & np.isfinite(resid)
        if m.sum() < 500:
            continue
        r = resid[m]
        dist_rows.append({"lo": lo, "hi": hi, "n": int(m.sum()),
                          "median_resid": float(np.median(r))})
        log.info("  xm_dist %4.1f-%-5.1f  n=%8d  median resid = %+.4f",
                 lo, hi, int(m.sum()), float(np.median(r)))

    # ---- verdict ----------------------------------------------------------
    ref = next(r for r in rows if r["class"].startswith("RETAINED"))
    ref_clean = ref["clean"]["ratio"] if ref["clean"] and ref["clean"]["ratio"] else 1.0

    survivors = []
    for r in rows:
        if r["class"].startswith("RETAINED"):
            continue
        rc = r["clean"]["ratio"] if r["clean"] and r["clean"]["ratio"] else None
        rb = r["suspect"]["ratio"] if r["suspect"] and r["suspect"]["ratio"] else None
        if rc is None:
            continue
        # survives if it is still far above the reference AFTER cleaning,
        # and is not simply much worse when blended
        if rc > 3.0 * ref_clean:
            survivors.append({"class": r["class"], "clean_ratio": rc,
                              "suspect_ratio": rb})

    if not survivors:
        verdict = (
            "APERTURE MISMATCH. Every flagged discard class collapses toward "
            "the reference once the 2MASS cross-match is required to be "
            "single, uncontested and well-centred. The one-sided dim excess "
            "found in Search M is the 4 arcsec beam picking up flux that Gaia "
            "resolves out, which inflates Ks and is read as a G deficit. "
            "Search M's asymmetry is explained and carries no signal.")
    else:
        names = ", ".join(s["class"] for s in survivors)
        verdict = (
            f"SURVIVES the blend test in: {names}. These remain far more "
            f"one-sidedly dim than the reference even when the 2MASS match is "
            f"single, uncontested and centred to better than 0.3 arcsec, so "
            f"the beam does not explain them. Each needs its own astrophysical "
            f"account before anything else is claimed -- for variables the "
            f"first candidate is starspot coverage, which suppresses optical "
            f"flux more than infrared and is genuinely one-signed.")

    print(f"\n{'='*76}")
    print("SEARCH M FOLLOW-UP: IS THE DISCARD DEFICIT APERTURE MISMATCH?")
    print(f"{'='*76}")
    print(f"{'class':<28} {'clean n':>9} {'clean r':>9} "
          f"{'susp n':>9} {'susp r':>9}")
    for r in rows:
        c, b = r["clean"], r["suspect"]
        print(f"{r['class']:<28} "
              f"{(c['n'] if c else 0):>9} "
              f"{(f'{c['ratio']:.2f}' if c and c['ratio'] else '-'):>9} "
              f"{(b['n'] if b else 0):>9} "
              f"{(f'{b['ratio']:.2f}' if b and b['ratio'] else '-'):>9}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchM_blendtest_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag, "threshold_mag": float(thr),
        "fiducial_sigma": float(sigma),
        "reference_clean_ratio": float(ref_clean),
        "classes": rows,
        "deficit_vs_xm_dist": dist_rows,
        "survivors": survivors,
        "verdict": verdict,
    }, indent=2, default=str))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
