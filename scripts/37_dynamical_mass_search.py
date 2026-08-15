#!/usr/bin/env python
"""I2: optical deficit against DYNAMICAL mass, not against a photometric proxy.

    run.sh scripts/37_dynamical_mass_search.py

WHY THIS ESCAPES THE FRAME
--------------------------
Every limit in this paper so far measures M_G against M_Ks. That inherits two
structural blind spots:

  * a GREY absorber dims G and Ks together and produces a residual of the
    WRONG SIGN (paper Sec 3.3), so grey harvesting is invisible or anti-visible;
  * the anchor is 2MASS, whose 4-arcsec beam is the dominant contaminant
    (Sec 5.5).

A dynamical mass has neither problem. Gravity does not care what an absorber
does to photons, and the mass comes from Gaia astrometry rather than from a
low-resolution infrared survey. Regressing M_G directly on log(M_dyn) therefore
probes absorbers of ANY spectral slope, including exactly the alpha ~ 0.19
blind spot and the grey case.

We already pulled 99,174 systems with dynamical masses (gaiadr3.binary_masses)
and used them only as a flat-absorber control. This runs them as a search.

THE COST, STATED UP FRONT
-------------------------
Dynamical masses come from binaries, and in an unresolved binary the observed G
contains BOTH components while m1 is the primary alone. That inflates the
luminosity and is the dominant systematic. It biases toward OVER-luminosity,
i.e. against the signal, so the search is conservative -- but it is noisy, and
that noise sets the sensitivity.
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
from pipeline import sample as smp
from pipeline import statistics as st


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    df = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    print(f"systems with dynamical masses: {len(df):,}")

    df = smp.add_astrometry(df)
    a0 = ext.query_a0("edenhofer23", df["l"].to_numpy(float),
                      df["b"].to_numpy(float), df["dist_pc"].to_numpy(float))
    bp_rp = df["bp_rp"].to_numpy(float)
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = df["dist_mod"].to_numpy(float)
    df["A_0"] = a0
    df["M_G"] = df["phot_g_mean_mag"].to_numpy(float) - mu - a_g
    df["M_Ks"] = df["tmass_ks_m"].to_numpy(float) - mu - a_ks

    good = (np.isfinite(df["M_G"]) & np.isfinite(a0) & df["m1"].notna()
            & (df["m1"] > 0.2) & (df["m1"] < 1.6)
            & (df["tmass_ph_qual"].str[2] == "A")
            & (a0 < 0.6))
    d = df[good].reset_index(drop=True)
    print(f"usable: {len(d):,}  (0.2-1.6 Msun, A_0<0.6, clean Ks)")

    x = np.log10(d["m1"].to_numpy(float))
    for band, col in (("M_G", "M_G"), ("M_Ks", "M_Ks")):
        y = d[col].to_numpy(float)
        best = None
        for deg in (2, 3, 4):
            c = np.polyfit(x, y, deg)
            s = st.robust_sigma(y - np.polyval(c, x))
            if best is None or s < best[1]:
                best = (c, s, deg)
        c, s, deg = best
        d[f"res_{band}"] = y - np.polyval(c, x)
        print(f"  {band} vs log M_dyn : deg {deg}, robust scatter {s:.4f} mag")

    # The flux-ratio subset removes most of the secondary-light systematic.
    if "fluxratio" in d and d["fluxratio"].notna().any():
        faint = d["fluxratio"].notna() & (d["fluxratio"] < 0.1)
        print(f"  systems with fluxratio<0.1 (faint secondary): "
              f"{int(faint.sum()):,}")
    else:
        faint = pd.Series(False, index=d.index)

    rows = []
    for label, m in (("all systems", pd.Series(True, index=d.index)),
                     ("faint secondary only", faint)):
        if m.sum() < 500:
            continue
        for band in ("M_G", "M_Ks"):
            r = d.loc[m, f"res_{band}"].to_numpy(float)
            s = st.robust_sigma(r); med = float(np.median(r))
            n = len(r)
            npos = int(np.count_nonzero(r > med + args.k * s))
            nneg = int(np.count_nonzero(r < med - args.k * s))
            thr = args.k * s
            f_det = float(st.fraction_from_delta(thr))
            ul = st.poisson_upper_limit(npos)
            rows.append({"subset": label, "band": band, "N": n,
                         "sigma": s, "threshold_mag": thr,
                         "f_detectable": f_det, "n_pos": npos, "n_neg": nneg,
                         "p_UL": ul / n, "mean_f_UL": ul / n * f_det})
    t = pd.DataFrame(rows)
    print("\n=== limits from dynamical mass ===")
    print(t.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))

    print("\n  M_G is the band that matters: it probes absorbers of ANY")
    print("  spectral slope, including the alpha~0.19 blind spot and the grey")
    print("  case that the M_G-vs-M_Ks estimator cannot see at all.")

    best_g = t[(t["band"] == "M_G")].loc[lambda z: z["mean_f_UL"].idxmin()]
    print(f"\n  best M_G limit: p < {float(best_g['p_UL']):.3e} at "
          f"f >= {float(best_g['f_detectable']):.3f}  "
          f"({best_g['subset']}, N={int(best_g['N']):,})")
    print(f"  photometric estimator for comparison: p < 4.3e-4 at f >= 0.507")

    out = {"n_systems": int(len(d)), "k": args.k,
           "table": t.to_dict(orient="records"),
           "best_MG_p_UL": float(best_g["p_UL"]),
           "best_MG_f": float(best_g["f_detectable"]),
           "note": ("M_G vs dynamical mass is immune to the grey-absorber and "
                    "alpha~0.19 blind spots and does not use 2MASS as the "
                    "anchor; unresolved secondary light is the dominant "
                    "systematic and biases toward over-luminosity, i.e. "
                    "against the signal")}
    (cfg.RESULT_DIR / "dynamical_mass_search.json").write_text(
        json.dumps(out, indent=2, default=float))
    t.to_csv(cfg.RESULT_DIR / "dynamical_mass_search.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
