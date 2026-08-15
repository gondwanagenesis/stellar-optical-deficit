#!/usr/bin/env python
"""Optical deficit WITHOUT infrared excess: a limit on the beamed class.

    run.sh scripts/29_ir_veto.py --tag primary

THE GAP THIS FILLS
------------------
Suazo et al. (2022) obtain p < 1.9e-4 at gamma >= 0.5 by modelling
L = (1-gamma)L_star + gamma*BB(T_DS): obscuration AND isotropic re-emission at
an assumed temperature.  That is a strong limit on spheres that radiate their
waste heat into 4*pi at ~300 K.

It is not a limit on a sphere that beams its waste heat, radiates far outside
the WISE bands, or exports the energy non-thermally.  That class is the entire
motivation for the optical-deficit channel (Section 1 of the paper), and nobody
has put a number on it -- including us, so far.

This script does.  It asks: among stars whose mid-infrared is MEASURED and
consistent with a bare photosphere, how many show an optical deficit?

The veto cuts both ways, which is the point:
  * it removes YSOs, debris discs and dusty contaminants, which are a large
    part of our 19,844-strong positive tail -- so the limit should IMPROVE;
  * and what survives is, by construction, the beaming-consistent population.

REQUIRING A MEASUREMENT, NOT AN ABSENCE
---------------------------------------
Section 5.4 of the paper found 278 of the 300 strongest outliers had no AllWISE
photometry at all.  Counting those as "no excess" would be the single easiest
way to fabricate a candidate list.  Every star here must have a real W1 and W2
detection with finite uncertainties; absence of data is never treated as
absence of excess.
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

# Wang & Chen 2019, ApJ 877, 116, Table 3: A_band / A_V for WISE.
# Both are <0.04 A_V, i.e. <0.02 mag over our A_V < 0.6 sample, but included
# so the dereddening is internally consistent with the optical bands.
A_W1_OVER_AV = 0.039
A_W2_OVER_AV = 0.026

EXCESS_NSIGMA = 3.0     # |excess| below this counts as "bare photosphere"


def running_expected(x, y, nbins=60, lo=None, hi=None):
    lo = np.nanpercentile(x, 1) if lo is None else lo
    hi = np.nanpercentile(x, 99) if hi is None else hi
    edges = np.linspace(lo, hi, nbins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, nbins - 1)
    med = np.full(nbins, np.nan)
    for b in range(nbins):
        m = (idx == b) & np.isfinite(y)
        if m.sum() > 100:
            med[b] = np.median(y[m])
    good = np.isfinite(med)
    if good.sum() < 5:
        return np.full(len(x), np.nan)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return np.interp(x, centres[good], med[good])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    r = d["residual"].to_numpy(float)
    sigma = st.robust_sigma(r)
    med_r = float(np.median(r))
    n_all = len(d)
    print(f"sample {n_all:,}   sigma = {sigma:.5f} mag")

    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_ks = ext.deredden("Ks", a0, bp_rp)
    a_j = ext.deredden("J", a0, bp_rp)
    ks0 = d["tmass_ks_m"].to_numpy(float) - a_ks
    j0 = d["tmass_j_m"].to_numpy(float) - a_j
    jks0 = j0 - ks0

    w1 = d["wise_w1mpro"].to_numpy(float) - A_W1_OVER_AV * a0
    w2 = d["wise_w2mpro"].to_numpy(float) - A_W2_OVER_AV * a0
    e_w1 = d["wise_w1mpro_error"].to_numpy(float)
    e_w2 = d["wise_w2mpro_error"].to_numpy(float)

    have = (np.isfinite(w1) & np.isfinite(w2) & np.isfinite(e_w1)
            & np.isfinite(e_w2) & np.isfinite(jks0)
            & (e_w1 < 0.2) & (e_w2 < 0.2))
    print(f"with usable W1 AND W2 photometry : {have.sum():,} "
          f"({100*have.mean():.1f}%)")

    # Photospheric colour, calibrated on the sample itself as a function of a
    # NIR temperature proxy. Main-sequence stars are Rayleigh-Jeans here so
    # (Ks-W1) is close to zero and only weakly colour dependent.
    c1 = ks0 - w1
    c2 = ks0 - w2
    exp1 = running_expected(jks0, np.where(have, c1, np.nan))
    exp2 = running_expected(jks0, np.where(have, c2, np.nan))
    ex1 = c1 - exp1          # positive = IR excess (W1 brighter than photosphere)
    ex2 = c2 - exp2

    s1 = st.robust_sigma(ex1[have])
    s2 = st.robust_sigma(ex2[have])
    print(f"scatter of (Ks-W1) excess : {s1:.4f} mag")
    print(f"scatter of (Ks-W2) excess : {s2:.4f} mag")

    bare = have & (np.abs(ex1) < EXCESS_NSIGMA * s1) & \
        (np.abs(ex2) < EXCESS_NSIGMA * s2)
    excess = have & ((ex1 > EXCESS_NSIGMA * s1) | (ex2 > EXCESS_NSIGMA * s2))
    print(f"\nbare photosphere (|excess| < {EXCESS_NSIGMA:.0f} sigma in BOTH) : "
          f"{bare.sum():,}")
    print(f"mid-IR EXCESS detected                              : "
          f"{excess.sum():,}")
    print(f"no usable WISE (never counted as 'no excess')        : "
          f"{(~have).sum():,}")

    # --- does the IR veto clean the optical tail? ------------------------
    pos = r > med_r + args.k * sigma
    print(f"\n=== does the mid-IR veto clean the optical deficit tail? ===")
    for lab, m in [("all stars", np.ones(n_all, bool)),
                   ("with WISE", have),
                   ("bare photosphere", bare),
                   ("IR excess present", excess)]:
        n = int(m.sum())
        npos = int((pos & m).sum())
        nneg = int(((r < med_r - args.k * sigma) & m).sum())
        rate = npos / n if n else np.nan
        print(f"  {lab:22s} N={n:9,}  pos={npos:7,}  neg={nneg:6,}  "
              f"rate={rate:.5f}  pos/neg={npos/max(nneg,1):6.1f}")

    # --- the limit on the beaming-consistent class -----------------------
    rows = []
    for lab, m in [("all stars (paper Sec 5.3)", np.ones(n_all, bool)),
                   ("bare photosphere only", bare)]:
        rr = r[m]
        n = len(rr)
        s = st.robust_sigma(rr)
        mm = float(np.median(rr))
        npos = int(np.count_nonzero(rr > mm + args.k * s))
        for f in (0.3, 0.5, 0.9):
            eff = st.detection_efficiency(rr, f, args.k, sigma=s, median=mm)
            if eff <= 0:
                continue
            p_ul = st.poisson_upper_limit(npos) / (n * eff)
            rows.append({"subsample": lab, "N": n, "f": f, "sigma": s,
                         "n_pos": npos, "efficiency": eff,
                         "p_UL": min(p_ul, 1.0),
                         "mean_f_UL": min(p_ul, 1.0) * f})
    t = pd.DataFrame(rows)
    print(f"\n=== limits ===")
    print(t.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))

    base = t[(t["subsample"].str.startswith("all")) & (t["f"] == 0.5)]
    veto = t[(t["subsample"].str.startswith("bare")) & (t["f"] == 0.5)]
    if len(base) and len(veto):
        gain = float(base["p_UL"].iloc[0] / veto["p_UL"].iloc[0])
        print(f"\n  IR veto changes p_UL at f=0.5 by {gain:.2f}x "
              f"({float(base['p_UL'].iloc[0]):.3e} -> "
              f"{float(veto['p_UL'].iloc[0]):.3e})")

    t.to_csv(cfg.RESULT_DIR / f"ir_veto_{args.tag}.csv", index=False)
    out = {
        "tag": args.tag,
        "n_all": int(n_all), "n_with_wise": int(have.sum()),
        "n_bare_photosphere": int(bare.sum()),
        "n_ir_excess": int(excess.sum()),
        "n_no_wise": int((~have).sum()),
        "excess_scatter_w1": float(s1), "excess_scatter_w2": float(s2),
        "table": t.to_dict(orient="records"),
        "interpretation": (
            "p_UL in the 'bare photosphere' rows is a limit on stars with an "
            "optical deficit AND a measured absence of mid-IR excess, i.e. the "
            "beaming-consistent / non-thermally-exporting class that "
            "IR-excess searches such as Suazo+22 are not designed to constrain."),
    }
    (cfg.RESULT_DIR / f"ir_veto_{args.tag}.json").write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
