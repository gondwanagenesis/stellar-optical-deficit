#!/usr/bin/env python
"""Search U: does the energy actually balance?

    run.sh scripts/75_searchU_energy_balance.py --tag primary

WHAT CHANNEL 9 DID NOT ASK
--------------------------
Channel 9 selected stars that are BOTH optically dimmed and infrared-bright,
found 2,632 of them, and concluded debris discs. That selection is a
coincidence requirement: each quantity had to be individually significant.

It never checked whether the two amounts MATCH.

Energy conservation is the whole physics of the intercept-and-re-emit
hypothesis. A structure that absorbs a fraction of a star's optical output has
to re-radiate that much and no more. So the discovery statistic is not the
deficit, and not the excess, but their RATIO -- and a real absorber sits on a
locus, not merely in a quadrant.

WHY THAT IS FAR MORE SELECTIVE
------------------------------
The two dominant natural populations fail the balance test in opposite,
diagnostic ways.

  debris discs        optically thin and rarely edge-on, so they are heated by
                      the star and re-radiate WITHOUT dimming it. Infrared
                      gain with little or no optical loss: the ratio runs away
                      to infinity.

  extinction, blends  dim the optical through material that is cold, distant
                      or simply elsewhere, contributing nothing at 12 microns.
                      Optical loss with no infrared gain: the ratio goes to
                      zero.

  a real absorber     gain equals loss. The ratio is of order unity.

So the natural populations occupy the two ends and the signal occupies the
middle. Searching either quantity alone throws away half the information and
all of the specificity.

THE HONEST CAVEAT ON THE BOLOMETRIC CORRECTION
----------------------------------------------
W3 is a narrow band, not a bolometer. Material at 100-300 K peaks between 10
and 29 microns, so W3 at 12 microns and W4 at 22 microns capture a large but
temperature-dependent share of the re-radiated energy. The ratio computed here
is therefore a band-limited proxy, and its absolute normalisation depends on
the assumed dust temperature. That is why the test is run as a DISTRIBUTION
comparison against controls rather than as an absolute threshold at unity: the
question asked is whether any population sits at a systematically different
balance point from the debris discs, not whether the ratio equals 1.000.
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
log = logging.getLogger("searchU")

# Effective wavelengths (micron) and Vega zero points (Jy)
LAM = {"G": 0.6230, "W3": 11.56, "W4": 22.09}
ZP = {"G": 3229.0, "W3": 31.674, "W4": 8.363}
C_UM_HZ = 2.998e14


def nu_fnu(mag, band):
    """nu*F_nu in arbitrary consistent units, from a Vega magnitude."""
    f_jy = ZP[band] * 10.0 ** (-0.4 * np.asarray(mag, float))
    return (C_UM_HZ / LAM[band]) * f_jy


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    src = cfg.RESULT_DIR / f"searchB_cold_candidates_{args.tag}.csv"
    if not src.exists():
        log.error("missing %s -- run scripts/46_searchB_cold_excess.py first",
                  src)
        return 1

    d = pd.read_csv(src)
    log.info("channel 9 candidates: %d", len(d))
    log.info("columns: %s", ", ".join(d.columns[:24]))

    need = ["resid_optical", "resid_w3", "M_G", "wise_w3mpro"]
    missing = [c for c in need if c not in d.columns]
    if missing:
        log.error("missing columns: %s", ", ".join(missing))
        return 1

    r_opt = d["resid_optical"].to_numpy(float)     # >0 means dim
    r_w3 = d["resid_w3"].to_numpy(float)           # <0 means IR-bright
    r_w4 = d["resid_w4"].to_numpy(float) if "resid_w4" in d else np.full(len(d), np.nan)

    # ---- convert residuals to fractional flux changes --------------------
    # observed / expected = 10^(-0.4 * residual)
    frac_lost_opt = 1.0 - 10.0 ** (-0.4 * r_opt)        # >0 when dimmed
    frac_gain_w3 = 10.0 ** (-0.4 * r_w3) - 1.0          # >0 when brightened
    frac_gain_w4 = 10.0 ** (-0.4 * r_w4) - 1.0

    # ---- convert to energy, using the star's own measured bands ----------
    # The expected (unabsorbed) optical energy is the observed one divided by
    # what survived; the expected photospheric infrared is the observed one
    # divided by what was added.
    nfn_g_obs = nu_fnu(d["M_G"].to_numpy(float), "G")
    nfn_w3_obs = nu_fnu(d["wise_w3mpro"].to_numpy(float), "W3")

    with np.errstate(divide="ignore", invalid="ignore"):
        e_opt_expected = nfn_g_obs / np.maximum(10.0 ** (-0.4 * r_opt), 1e-9)
        e_lost = e_opt_expected * frac_lost_opt
        w3_photosphere = nfn_w3_obs / np.maximum(1.0 + frac_gain_w3, 1e-9)
        e_gained = w3_photosphere * frac_gain_w3
        balance = e_gained / np.maximum(e_lost, 1e-30)

    ok = (np.isfinite(balance) & (e_lost > 0) & np.isfinite(e_gained)
          & (frac_lost_opt > 0))
    log.info("with a computable energy balance: %d", int(ok.sum()))
    if ok.sum() < 50:
        log.error("too few usable objects")
        return 1

    b = balance[ok]
    logb = np.log10(np.maximum(b, 1e-12))
    log.info("")
    log.info("energy balance  (IR gained / optical lost), log10:")
    for q in (1, 5, 16, 50, 84, 95, 99):
        log.info("   %2d th pct : %+.2f  (ratio %.3g)",
                 q, float(np.percentile(logb, q)),
                 10 ** float(np.percentile(logb, q)))

    frac_gg1 = float((b > 1).mean())
    log.info("fraction with IR gain EXCEEDING optical loss: %.3f", frac_gg1)

    # ---- where does a real absorber sit, and is anything there? ----------
    # Balanced means the re-radiated band-limited energy is within a factor
    # of a few of the intercepted optical energy.
    balanced = ok.copy()
    balanced[ok] = (b > 0.2) & (b < 5.0)
    runaway = ok.copy()
    runaway[ok] = b >= 5.0            # disc-like: gain without loss
    starved = ok.copy()
    starved[ok] = b <= 0.2            # extinction-like: loss without gain

    log.info("")
    log.info("balanced   (0.2 < ratio < 5)  : %d", int(balanced.sum()))
    log.info("runaway    (ratio >= 5)       : %d   disc-like", int(runaway.sum()))
    log.info("starved    (ratio <= 0.2)     : %d   extinction-like",
             int(starved.sum()))

    # ---- what distinguishes the balanced group? --------------------------
    log.info("")
    log.info("=== properties by balance class ===")
    rows = []
    for name, m in (("balanced", balanced), ("runaway", runaway),
                    ("starved", starved)):
        if m.sum() < 20:
            continue
        rec = {"class": name, "n": int(m.sum())}
        for col, lab in (("w3w4_colour", "W3-W4"), ("A_0", "A_0"),
                         ("b", "|b|"), ("resid_optical", "deficit"),
                         ("implied_f_optical", "implied f")):
            if col in d.columns:
                v = d.loc[m, col]
                v = v.abs() if col == "b" else v
                rec[lab] = float(v.median())
        rows.append(rec)
        log.info("  %-9s n=%5d  %s", name, int(m.sum()),
                 "  ".join(f"{k}={v:+.3f}" for k, v in rec.items()
                           if k not in ("class", "n")))

    # ---- the control that makes it a measurement -------------------------
    # Reverse the roles: stars that are optically BRIGHT and IR-FAINT cannot
    # host an absorber at all, so running the identical arithmetic on them
    # gives the rate at which this estimator manufactures a balanced ratio.
    r_opt_m = -r_opt
    r_w3_m = -r_w3
    fl_m = 1.0 - 10.0 ** (-0.4 * r_opt_m)
    fg_m = 10.0 ** (-0.4 * r_w3_m) - 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        e_lost_m = (nfn_g_obs / np.maximum(10.0 ** (-0.4 * r_opt_m), 1e-9)) * fl_m
        w3p_m = nfn_w3_obs / np.maximum(1.0 + fg_m, 1e-9)
        bal_m = (w3p_m * fg_m) / np.maximum(e_lost_m, 1e-30)
    ok_m = np.isfinite(bal_m) & (e_lost_m > 0) & (fl_m > 0)
    n_bal_m = int(((bal_m > 0.2) & (bal_m < 5.0) & ok_m).sum())
    log.info("")
    log.info("mirror control (signs reversed): %d usable, %d 'balanced'",
             int(ok_m.sum()), n_bal_m)

    ratio_vs_mirror = (int(balanced.sum()) / n_bal_m) if n_bal_m else np.inf
    log.info("balanced : mirror-balanced = %s",
             f"{ratio_vs_mirror:.2f}" if np.isfinite(ratio_vs_mirror) else "inf")

    # ---- verdict ----------------------------------------------------------
    n_bal = int(balanced.sum())
    if n_bal == 0:
        verdict = (
            "NULL. No channel-9 candidate re-radiates in the mid-infrared "
            "anything close to what it appears to lose in the optical. The "
            "population splits into disc-like objects that gain without "
            "losing and extinction-like objects that lose without gaining, "
            "which is what the two natural explanations predict and what an "
            "intercept-and-re-emit structure cannot look like.")
    elif np.isfinite(ratio_vs_mirror) and ratio_vs_mirror < 2.0:
        verdict = (
            f"NULL. {n_bal} candidates fall in the energy-balanced band, but "
            f"the sign-reversed mirror produces {n_bal_m} by the same "
            f"arithmetic ({ratio_vs_mirror:.2f}:1). The balanced group is "
            f"what this estimator generates from photometric scatter, not a "
            f"population that conserves energy.")
    else:
        verdict = (
            f"{n_bal} candidates sit in the energy-balanced band against "
            f"{n_bal_m} in the sign-reversed mirror ({ratio_vs_mirror:.1f}:1). "
            f"These lose optical energy and regain a comparable amount at 12 "
            f"microns, which neither a debris disc (gain without loss) nor "
            f"extinction (loss without gain) does. Before this means "
            f"anything: check the W3-W4 colour for a disc temperature, "
            f"confirm the optical deficit is not aperture mismatch by "
            f"requiring a single centred 2MASS match, and note that the "
            f"absolute normalisation depends on the assumed re-emission "
            f"temperature because W3 is a band, not a bolometer.")

    print(f"\n{'='*72}")
    print("SEARCH U: DOES THE ENERGY BALANCE?")
    print(f"{'='*72}")
    print(f"  channel-9 candidates                : {len(d):,}")
    print(f"  with a computable balance           : {int(ok.sum()):,}")
    print(f"  median log10(IR gain / optical loss): {float(np.median(logb)):+.2f}")
    print()
    print(f"  balanced  (0.2 < ratio < 5)         : {n_bal}")
    print(f"  runaway   (>= 5, disc-like)         : {int(runaway.sum())}")
    print(f"  starved   (<= 0.2, extinction-like) : {int(starved.sum())}")
    print(f"  mirror-control balanced             : {n_bal_m}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchU_energy_balance_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "n_candidates": int(len(d)),
        "n_usable": int(ok.sum()),
        "log10_balance_percentiles": {
            str(q): float(np.percentile(logb, q))
            for q in (1, 5, 16, 50, 84, 95, 99)},
        "frac_gain_exceeds_loss": frac_gg1,
        "n_balanced": n_bal,
        "n_runaway_disclike": int(runaway.sum()),
        "n_starved_extinctionlike": int(starved.sum()),
        "n_mirror_balanced": n_bal_m,
        "balanced_over_mirror": (
            float(ratio_vs_mirror) if np.isfinite(ratio_vs_mirror) else None),
        "class_properties": rows,
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)

    if n_bal:
        cand = d[balanced].copy()
        cand["energy_balance"] = balance[balanced]
        cand.sort_values("energy_balance").to_csv(
            cfg.RESULT_DIR / f"searchU_balanced_{args.tag}.csv", index=False)
        log.info("wrote %d balanced candidates", n_bal)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
