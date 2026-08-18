#!/usr/bin/env python
"""Diagnose the Search E anomaly box before believing any of it.

    run.sh scripts/62_searchE_diagnosis.py --tag primary

WHY THIS EXISTS
---------------
Search E found 486 systems with a 0.20-0.45 Msun dark companion at P > 1500 d
and e > 0.15, against ZERO in the circular control at the same mass and
period. A 486:0 asymmetry looks decisive.

The zero is what makes it suspicious. If fitted eccentricity carried real
scatter, some long-period systems would land at low e by chance. Finding none
says that beyond the baseline the fit essentially never returns a small
eccentricity -- which would be a property of the fitting, not of the
population, and would manufacture the entire signal.

Three tests, each of which can kill it:

1. THE BASELINE TEST. DR3's astrometric baseline is ~1000 d. Every one of the
   486 has P > 1500 d, so every orbit is extrapolated. If the fraction of
   eccentric systems jumps across the baseline, the anomaly is an artefact of
   fitting periods the data cannot constrain.

2. THE MASS-INDEPENDENCE TEST. This is the decisive one. The physical claim is
   specific to LOW-mass companions: a helium white dwarf must have been
   stripped, and stripping welds mass to period. Companions above ~0.5 Msun
   are CO white dwarfs, which have their own wide-orbit channel and are NOT
   anomalous there. So if wide eccentric systems are equally common at high
   companion mass, the selection is not testing the helium-white-dwarf
   relation at all -- it is just selecting long-period Gaia solutions.

3. THE SOLUTION-QUALITY TEST. If the 486 have systematically worse fits than
   the parent sample, they are the tail of the error distribution.
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchE_diag")

DR3_BASELINE_D = 1000.0
ECC_CUT = 0.15


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    d = pd.read_csv(cfg.RESULT_DIR / f"searchE_all_amrf_{args.tag}.csv")
    log.info("clean AMRF sample: %d", len(d))

    P = d["period"].to_numpy(float)
    e = d["eccentricity"].to_numpy(float)
    m2 = d["m2_min"].to_numpy(float)

    # ---- 1. the baseline test --------------------------------------------
    log.info("")
    log.info("=== 1. ECCENTRICITY ACROSS THE DR3 BASELINE ===")
    bins = [0, 200, 400, 700, 1000, 1500, 2500, 1e9]
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (P >= lo) & (P < hi)
        if m.sum() < 20:
            continue
        frac = float((e[m] > ECC_CUT).mean())
        rows.append({"P_lo": lo, "P_hi": hi, "n": int(m.sum()),
                     "median_e": float(np.median(e[m])),
                     "frac_ecc": frac})
        log.info("  P %5.0f-%-7.0f n=%6d  median e=%.3f  frac(e>%.2f)=%.3f",
                 lo, hi, int(m.sum()), float(np.median(e[m])), ECC_CUT, frac)

    within = (P < DR3_BASELINE_D)
    beyond = (P > 1500.0)
    frac_within = float((e[within] > ECC_CUT).mean())
    frac_beyond = float((e[beyond] > ECC_CUT).mean())
    log.info("  frac(e>%.2f) within baseline : %.3f (n=%d)",
             ECC_CUT, frac_within, int(within.sum()))
    log.info("  frac(e>%.2f) beyond 1500 d   : %.3f (n=%d)",
             ECC_CUT, frac_beyond, int(beyond.sum()))
    log.info("  ratio                        : %.2fx",
             frac_beyond / max(frac_within, 1e-9))

    # ---- 2. the mass-independence test -----------------------------------
    log.info("")
    log.info("=== 2. IS THE EFFECT SPECIFIC TO LOW-MASS COMPANIONS? ===")
    mass_bins = [(0.20, 0.45), (0.45, 0.60), (0.60, 0.90),
                 (0.90, 1.40), (1.40, 10.0)]
    mrows = []
    for lo, hi in mass_bins:
        m = beyond & (m2 >= lo) & (m2 < hi)
        if m.sum() < 5:
            log.info("  M2 %.2f-%.2f : n=%d (too few)", lo, hi, int(m.sum()))
            continue
        frac = float((e[m] > ECC_CUT).mean())
        mrows.append({"m2_lo": lo, "m2_hi": hi, "n": int(m.sum()),
                      "frac_ecc": frac,
                      "median_e": float(np.median(e[m]))})
        log.info("  M2 %.2f-%.2f Msun, P>1500d: n=%5d  frac(e>%.2f)=%.3f  "
                 "median e=%.3f", lo, hi, int(m.sum()), ECC_CUT, frac,
                 float(np.median(e[m])))

    if len(mrows) >= 2:
        fr = np.array([r["frac_ecc"] for r in mrows])
        spread = float(fr.max() - fr.min())
        log.info("  spread in frac(e>%.2f) across mass bins: %.3f",
                 ECC_CUT, spread)
        mass_independent = spread < 0.15
    else:
        spread = float("nan")
        mass_independent = False

    # ---- 3. solution quality ---------------------------------------------
    log.info("")
    log.info("=== 3. SOLUTION QUALITY OF THE ANOMALY BOX ===")
    box = beyond & (m2 >= 0.20) & (m2 <= 0.45) & (e > ECC_CUT)
    rest = ~box
    for col in ("goodness_of_fit", "significance", "eccentricity_error",
                "period_error", "ruwe"):
        if col not in d.columns:
            continue
        a = d.loc[box, col].median()
        b = d.loc[rest, col].median()
        log.info("  %-20s box=%9.3f   rest=%9.3f   ratio=%.2f",
                 col, a, b, a / b if b else np.nan)

    # relative period error is the sharpest quality indicator here
    if "period_error" in d.columns:
        rel_p = d["period_error"].to_numpy(float) / np.maximum(P, 1e-9)
        log.info("  relative period error: box=%.3f  rest=%.3f",
                 float(np.median(rel_p[box])), float(np.median(rel_p[rest])))

    # ---- verdict ----------------------------------------------------------
    log.info("")
    all_beyond_baseline = bool((P[box] > DR3_BASELINE_D).all())

    if mass_independent:
        verdict = (
            f"ARTEFACT. The fraction of eccentric wide systems varies by only "
            f"{spread:.3f} across companion mass from 0.2 to >1.4 Msun. The "
            f"selection is therefore not testing the helium-white-dwarf "
            f"mass-period relation -- it is selecting long-period Gaia "
            f"solutions, which are eccentric regardless of what orbits them. "
            f"Combined with frac(e>{ECC_CUT}) rising {frac_beyond/max(frac_within,1e-9):.1f}x "
            f"across the {DR3_BASELINE_D:.0f} d baseline and every box member "
            f"having an extrapolated period, the 486:0 asymmetry is a property "
            f"of fitting orbits longer than the data, not of the companions.")
    else:
        verdict = (
            f"SURVIVES the mass-independence test: the eccentric fraction "
            f"varies by {spread:.3f} across companion mass, so the effect is "
            f"not purely a long-period fitting artefact. It still requires "
            f"radial-velocity confirmation, since all {int(box.sum())} box "
            f"members have periods beyond the DR3 baseline.")

    print(f"\n{'='*70}")
    print("SEARCH E DIAGNOSIS: IS THE 486:0 ASYMMETRY REAL?")
    print(f"{'='*70}")
    print(f"  frac(e>{ECC_CUT}) within the {DR3_BASELINE_D:.0f} d baseline : "
          f"{frac_within:.3f}")
    print(f"  frac(e>{ECC_CUT}) beyond 1500 d                : {frac_beyond:.3f}")
    print(f"  ratio across the baseline                  : "
          f"{frac_beyond/max(frac_within,1e-9):.2f}x")
    print(f"  every box member extrapolated              : {all_beyond_baseline}")
    print()
    for r in mrows:
        print(f"  M2 {r['m2_lo']:.2f}-{r['m2_hi']:.2f} Msun, P>1500d : "
              f"n={r['n']:5d}  frac(e>{ECC_CUT})={r['frac_ecc']:.3f}")
    print(f"  spread across mass bins                    : {spread:.3f}")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "n_sample": int(len(d)),
        "ecc_cut": ECC_CUT,
        "dr3_baseline_d": DR3_BASELINE_D,
        "period_bins": rows,
        "frac_ecc_within_baseline": frac_within,
        "frac_ecc_beyond_1500d": frac_beyond,
        "baseline_ratio": float(frac_beyond / max(frac_within, 1e-9)),
        "mass_bins_beyond_1500d": mrows,
        "frac_ecc_spread_across_mass": spread,
        "mass_independent": bool(mass_independent),
        "all_box_beyond_baseline": all_beyond_baseline,
        "n_box": int(box.sum()),
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchE_diagnosis_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
