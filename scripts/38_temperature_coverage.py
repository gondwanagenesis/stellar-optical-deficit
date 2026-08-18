#!/usr/bin/env python
"""What re-radiation temperatures does this work actually cover?

    run.sh scripts/38_temperature_coverage.py

Our mid-IR veto (paper Sec 5.7) uses WISE W1 and W2 at 3.4 and 4.6 micron.
Suazo et al. use W1-W4, out to 22 micron. Both are blind below some
temperature, and "runs cold" is not an exotic design choice -- it is what a
computation-limited civilisation is driven to, because the Landauer cost of a
bit erasure is kT ln2 and the value of energy is set by the SINK temperature.

This computes, for a star enclosed at covering fraction f and re-radiating at
temperature T, the flux density in each available far-infrared and submillimetre
band, and compares it to that survey's detection limit. The output is the
temperature interval that is genuinely constrained versus the interval where
nobody has looked.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg

H, C, KB = 6.62607015e-34, 2.99792458e8, 1.380649e-23
LSUN = 3.828e26          # W
PC = 3.0856775814913673e16

# band centre (micron), 5-sigma point-source sensitivity (Jy), survey
BANDS = [
    (3.4,   5e-5,   "WISE W1"),
    (4.6,   1e-4,   "WISE W2"),
    (12.0,  1e-3,   "WISE W3"),
    (22.0,  6e-3,   "WISE W4"),
    (25.0,  0.5,    "IRAS 25um"),
    (60.0,  0.6,    "IRAS 60um"),
    (65.0,  2.4,    "AKARI FIS N60"),
    (90.0,  0.55,   "AKARI FIS WIDE-S"),
    (100.0, 1.0,    "IRAS 100um"),
    (140.0, 1.4,    "AKARI FIS WIDE-L"),
    (160.0, 6.3,    "AKARI FIS N160"),
    (250.0, 0.02,   "Herschel SPIRE 250 (pointed)"),
    (350.0, 0.02,   "Herschel SPIRE 350 (pointed)"),
    (850.0, 0.05,   "SCUBA-2 850 (pointed)"),
]


def flux_density_jy(T, f, L_star_lsun, d_pc, lam_um):
    """Fnu in Jy for a blackbody shell re-radiating f * L_star at temperature T."""
    L = f * L_star_lsun * LSUN
    lam = lam_um * 1e-6
    nu = C / lam
    x = np.clip(H * nu / (KB * T), 1e-8, 700.0)
    # fraction of a blackbody's power per unit frequency, normalised
    bnu = (2 * H * nu ** 3 / C ** 2) / np.expm1(x)
    sigma_t4 = 5.670374419e-8 * T ** 4
    # total emitted power = L, distributed as B_nu / (sigma T^4 / pi)
    fnu = L * (np.pi * bnu / sigma_t4) / (4 * np.pi * (d_pc * PC) ** 2)
    return fnu * 1e26     # W/m^2/Hz -> Jy


def main() -> int:
    rows = []
    for T in (20, 30, 50, 75, 100, 150, 200, 300, 600):
        for d_pc, f, L in ((100, 1.0, 1.0), (500, 1.0, 1.0), (500, 0.5, 1.0)):
            best, best_snr = None, 0.0
            for lam, lim, name in BANDS:
                fnu = flux_density_jy(T, f, L, d_pc, lam)
                snr = fnu / lim
                if snr > best_snr:
                    best, best_snr = name, snr
            rows.append({"T_K": T, "d_pc": d_pc, "f": f,
                         "best_band": best, "best_SNR": best_snr,
                         "detectable": best_snr > 5})
    t = pd.DataFrame(rows)
    print("=== detectability of a re-radiating shell, 1 Lsun star ===\n")
    for (d_pc, f), g in t.groupby(["d_pc", "f"]):
        print(f"  d = {d_pc} pc, covering fraction f = {f}")
        print(g[["T_K", "best_band", "best_SNR", "detectable"]].to_string(
            index=False, float_format=lambda v: f"{v:12.4g}"))
        print()

    print("=== what THIS PAPER's veto actually covers ===")
    print("  Sec 5.7 uses W1 (3.4um) and W2 (4.6um) only.")
    for T in (20, 30, 50, 100, 200, 400, 800):
        f34 = flux_density_jy(T, 1.0, 1.0, 500, 3.4)
        f46 = flux_density_jy(T, 1.0, 1.0, 500, 4.6)
        det = max(f34 / 5e-5, f46 / 1e-4)
        print(f"    T={T:4d} K  best W1/W2 SNR = {det:10.3g}  "
              f"{'covered' if det > 5 else 'INVISIBLE to our veto'}")

    # The gap: coldest temperature any all-sky survey reaches at 500 pc
    print("\n=== the coverage gap ===")
    allsky = [(lam, lim, n) for lam, lim, n in BANDS if "pointed" not in n]
    for T in (10, 15, 20, 25, 30, 40, 60):
        best = max(flux_density_jy(T, 1.0, 1.0, 500, lam) / lim
                   for lam, lim, n in allsky)
        print(f"    T={T:3d} K, f=1, 500 pc : best all-sky SNR = {best:9.3g}  "
              f"{'detectable' if best > 5 else 'NO ALL-SKY SURVEY REACHES THIS'}")

    t.to_csv(cfg.RESULT_DIR / "temperature_coverage.csv", index=False)
    out = {"bands": [{"lam_um": l, "limit_jy": q, "survey": n}
                     for l, q, n in BANDS],
           "note": ("W1/W2 veto in Sec 5.7 covers only hot shells; the far-IR "
                    "all-sky surveys (IRAS, AKARI FIS) extend coverage down to "
                    "~20-30 K at 500 pc; below that no all-sky survey reaches "
                    "and the CMB floor at 2.7 K is unreachable in principle")}
    (cfg.RESULT_DIR / "temperature_coverage.json").write_text(
        json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
