#!/usr/bin/env python
"""Quantify how spectrally selective an absorber must be for this test to see it.

    run.sh scripts/12_spectral_leverage.py --tag testrun_nir

This answers the brief's question -- "for a greybody absorber with
wavelength-dependent optical depth, what spectral slope is required before the
test has any leverage?" -- and, in doing so, corrects the premise behind it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import anchor
from pipeline import config as cfg

pd.set_option("display.width", 200)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="testrun_nir")
    ap.add_argument("--slope", type=float, default=None,
                    help="override dM_G/dM_Ks; default reads the fitted value")
    args = ap.parse_args()

    if args.slope is not None:
        slope = args.slope
        slope_src = "command line"
        s16 = s84 = np.nan
    else:
        meta = json.loads((cfg.RESULT_DIR / f"fiducial_{args.tag}.json").read_text())
        slope = meta["slope_dMG_dMKs_median"]
        s16, s84 = meta["slope_dMG_dMKs_p16"], meta["slope_dMG_dMKs_p84"]
        slope_src = f"fitted, {args.tag}"

    print(f"fiducial slope s = dM_G/dM_Ks = {slope:.3f}  ({slope_src})")
    if np.isfinite(s16):
        print(f"  16-84 percentile across the sample: {s16:.3f} .. {s84:.3f}")

    print("\nThe measured residual for an attenuated star is")
    print("    r = dm_G - s * dm_Ks")
    print("so the test is blind when dm_G / dm_Ks == s, NOT when the absorber "
          "is flat.\n")

    tab = anchor.leverage_table(slope)
    tab_disp = tab.rename(columns={
        "alpha": "alpha", "dm_G": "dm_G(mag)", "dm_Ks": "dm_Ks(mag)",
        "ratio_dmG_dmKs": "dm_G/dm_Ks", "residual_per_dmG": "leverage",
        "residual_mag": "residual(mag)"})
    print(tab_disp.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))

    a_blind_an = anchor.blind_slope_analytic(slope)          # SED-weighted
    a_blind_vega = anchor.blind_slope_analytic(slope, teff=None)
    a_blind_num = anchor.blind_slope_numeric(slope)
    grey_lev = anchor.leverage(0.0, slope)
    dust_lev = anchor.leverage(2.0, slope)

    lam_g_sed = anchor.effective_wavelength("G")
    lam_ks_sed = anchor.effective_wavelength("Ks")
    print(f"\neffective wavelengths for a 4500 K photosphere:")
    print(f"  G  : {lam_g_sed:.4f} um   (Vega-referenced catalogue value "
          f"{anchor.LAM_G_VEGA:.4f} um)")
    print(f"  Ks : {lam_ks_sed:.4f} um   (catalogue {anchor.LAM_KS:.4f} um)")
    print("  G is 400-950 nm wide, so a cool star's flux weighting pushes its "
          "effective\n  wavelength far to the red of the Vega value. Using the "
          "catalogue number here\n  would be wrong by ~10% in alpha_blind.")

    print(f"\nblind spot (analytic, SED-weighted) : alpha = {a_blind_an:.3f}")
    print(f"blind spot (analytic, Vega lambdas) : alpha = {a_blind_vega:.3f}  "
          f"<- not the right choice, shown for contrast")
    print(f"blind spot (full numeric)           : alpha = {a_blind_num:.3f}")
    print(f"leverage for a GREY absorber (a=0)  : {grey_lev:+.3f}")
    print(f"leverage for dust-like (a=2)        : {dust_lev:+.3f}")

    # Sensitivity of the blind spot to the slope uncertainty
    if np.isfinite(s16):
        print(f"\nblind alpha across the sample slope range: "
              f"{anchor.blind_slope_analytic(s16):.3f} .. "
              f"{anchor.blind_slope_analytic(s84):.3f}")

    out = {
        "slope_dMG_dMKs": float(slope),
        "alpha_blind_analytic_sed": float(a_blind_an),
        "alpha_blind_analytic_vega": float(a_blind_vega),
        "alpha_blind_numeric": float(a_blind_num),
        "leverage_grey_alpha0": float(grey_lev),
        "leverage_dustlike_alpha2": float(dust_lev),
        "lambda_eff_G_um_sed_4500K": float(lam_g_sed),
        "lambda_eff_Ks_um_sed_4500K": float(lam_ks_sed),
        "lambda_eff_G_um_vega": anchor.LAM_G_VEGA,
        "lambda_eff_Ks_um_vega": anchor.LAM_KS,
        "interpretation": {
            "flat_absorber": ("NOT invisible: leverage is negative, so a grey "
                              "absorber makes a star look OVER-luminous in G "
                              "at fixed M_Ks"),
            "true_blind_spot": ("alpha ~ ln(s)/ln(lambda_Ks/lambda_G); an "
                                "absorber already moderately selective"),
            "dust": ("interstellar dust has alpha ~ 2, far on the positive "
                     "side, so under-corrected reddening mimics the signal "
                     "directly"),
        },
    }
    (cfg.RESULT_DIR / "spectral_leverage.json").write_text(json.dumps(out, indent=2))
    tab.to_csv(cfg.RESULT_DIR / "spectral_leverage_table.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
