#!/usr/bin/env python
"""Joint constraint across both waste-heat disposal modes.

    run.sh scripts/31_joint_constraint.py

THE STATEMENT NOBODY HAS MADE
-----------------------------
Every published Dyson-sphere limit constrains ONE disposal mode. Infrared
searches (Carrigan 2009; Wright et al. 2014; Suazo et al. 2022, 2024) constrain
spheres that re-radiate isotropically at a warm temperature. Optical-deficit
searches (Annis 1999; Zackrisson et al. 2015, 2018; this work) constrain
interception irrespective of disposal, but weakly.

Partition the population of stars intercepting a fraction >= f:

    p_total = p_iso + p_dark

  p_iso  : re-radiates isotropically in the WISE bands  -> IR searches see it
  p_dark : beams, radiates too cold for WISE, or exports non-thermally
           -> IR searches do not see it; only a deficit search does

Section 5.7 of the paper measures p_dark DIRECTLY, by requiring a measured bare
photosphere in W1 and W2 alongside the optical deficit. Suazo et al. (2022)
measure p_iso. Together they bound p_total for the first time without assuming
a disposal mode, and they bound the beamed fraction

    beta = p_dark / p_total

which is the quantity the entire "infrared searches are evadable" argument
turns on, and which has never been measured.

WHAT THIS CANNOT DO
-------------------
Both inputs are upper limits with zero detections, so beta is bounded only
trivially from above (beta <= 1) and not at all from below: with no detections
in either channel we cannot say what fraction of a population that may not
exist chooses to beam. What IS new is (a) the two-sided bound on p_total, and
(b) an explicit statement of how much better the deficit channel must get
before beta becomes measurable. That number is the useful output.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg

# Suazo et al. 2022, MNRAS 512, 2988, Table: upper limits on the fraction of
# stars with a partial Dyson sphere, within 100 pc, T_DS = 300 K.
SUAZO22 = {0.1: 6.6e-3, 0.5: 1.9e-4, 0.9: 1.8e-5}


def main() -> int:
    R = cfg.RESULT_DIR
    beamed = json.loads((R / "beamed_class_primary.json").read_text())
    pair = json.loads((R / "pair_limit_primary.json").read_text())

    p_dark = float(beamed["best_p_UL"])
    f_dark = float(beamed["best_f"])
    n_dark = int(beamed["n_stars"])

    # Match Suazo's covering fraction to ours as closely as their grid allows.
    f_match = min(SUAZO22, key=lambda g: abs(g - f_dark))
    p_iso = SUAZO22[f_match]

    p_total = p_iso + p_dark

    print("=" * 72)
    print("JOINT CONSTRAINT ON STELLAR ENERGY HARVESTING, BOTH DISPOSAL MODES")
    print("=" * 72)
    print(f"\n  covering fraction        f >= {f_dark:.2f} (this work), "
          f"gamma >= {f_match:.1f} (Suazo+22)\n")
    print(f"  p_iso   (isotropic warm re-emission, Suazo+22, 100 pc)")
    print(f"          < {p_iso:.2e}")
    print(f"  p_dark  (beamed / cold / non-thermal, THIS WORK, 500 pc)")
    print(f"          < {p_dark:.2e}      [{n_dark:,} stars, measured bare "
          f"photosphere]")
    print(f"\n  p_total = p_iso + p_dark")
    print(f"          < {p_total:.2e}")
    print(f"\n  -> Fewer than 1 in {1/p_total:,.0f} nearby lower-main-sequence")
    print(f"     stars intercepts >={100*f_dark:.0f}% of its optical output by ANY")
    print(f"     means, warm or dark. This is the first bound that does not")
    print(f"     assume how the waste heat is disposed of.")

    # How much better must the deficit channel get for beta to be measurable?
    print("\n" + "-" * 72)
    print("WHEN DOES THE BEAMED FRACTION BECOME MEASURABLE?")
    print("-" * 72)
    print("\n  beta = p_dark / p_total is bounded only trivially today (both")
    print("  channels are non-detections). It becomes informative once the")
    print("  deficit limit drops below the IR limit, because then a detection")
    print("  in one and not the other is decisive.\n")
    need = p_iso
    factor = p_dark / need
    print(f"  deficit limit now          : {p_dark:.2e}")
    print(f"  must reach                 : {need:.2e}  (the IR limit)")
    print(f"  improvement required       : {factor:.1f}x")

    # The pair estimator is background-subtracted, so its limit scales ~1/N
    # rather than 1/sqrt(N) while the asymmetry stays consistent with zero.
    rows = []
    for mult in (1, 2, 5, 10, 20, 50):
        n = n_dark * mult
        # counts scale with N; the 95% UL on a zero asymmetry grows as sqrt(N)
        scaled_ul = float(beamed["table"][-1]["excess_UL_95"]) * np.sqrt(mult)
        p = scaled_ul / n
        rows.append({"pair_sample_multiple": mult, "stars": n,
                     "projected_p_UL": p,
                     "beats_IR_limit": p < need})
    t = pd.DataFrame(rows)
    print()
    print(t.to_string(index=False, float_format=lambda v: f"{v:12.4g}"))

    first = t[t["beats_IR_limit"]]
    if len(first):
        m = int(first.iloc[0]["pair_sample_multiple"])
        print(f"\n  => {m}x the current clean-pair sample "
              f"({int(first.iloc[0]['stars']):,} stars) puts the deficit")
        print(f"     channel below the isotropic limit, at which point beta")
        print(f"     becomes a measurable quantity rather than a rhetorical one.")
        print(f"     El-Badry+21 hold ~1.3e6 pairs; we currently use 6,431.")
    else:
        print("\n  => not reachable by scaling pairs alone within the grid tested.")

    out = {
        "f_this_work": f_dark, "gamma_suazo": f_match,
        "p_iso_upper": p_iso, "p_dark_upper": p_dark,
        "p_total_upper": p_total,
        "one_in_n_stars": float(1 / p_total),
        "improvement_needed_for_beta": float(factor),
        "projection": t.to_dict(orient="records"),
        "caveat": ("both inputs are non-detections, so beta is bounded only "
                   "trivially from above; the new content is the "
                   "disposal-agnostic bound on p_total and the explicit "
                   "threshold at which beta becomes measurable"),
    }
    (R / "joint_constraint.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {R / 'joint_constraint.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
