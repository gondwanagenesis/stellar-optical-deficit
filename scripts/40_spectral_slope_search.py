#!/usr/bin/env python
"""SEARCH A: measure the spectral slope of every anomaly across 7 bands.

    run.sh scripts/40_spectral_slope_search.py --tag primary

SPECIFICITY, NOT FLUX
---------------------
Broadband dimming has an appalling natural background: dust, blends, spots and
starspots all do it, which is why 50 years of deficit and excess searches keep
converging on "reddened dusty object". But dust has a SHAPE. Interstellar
extinction follows alpha ~ 2 over the optical-NIR, set by grain physics. An
absorber with alpha far from 2 -- and especially a needle-sharp or inverted
slope -- has no natural producer.

So instead of asking "is this star faint?" we ask "faint in a pattern nothing
natural makes?".

MODEL
-----
Anchor on M_Ks and a NIR temperature proxy (J-H)_0, predict absolute magnitude
in each of BP, G, RP, J, H, W1, W2 from the bulk of the sample, then fit the
7 residuals with a two-parameter power-law absorber:

    dm(lambda) = A * [ (lambda/lambda_Ks)^(-alpha) - 1 ]

The -1 is because K_s is the anchor, so everything is measured differentially
against it. Consequence, stated plainly: a GREY absorber (alpha = 0) gives
dm = 0 in every band and is invisible here by construction. Grey is covered by
the dynamical-mass channel (scripts/37) instead, and the two must be read
together.

What this channel CAN do is measure alpha for every non-grey anomaly and ask
whether any of them sit away from the dust value.
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

# effective wavelengths, micron
LAM = {"BP": 0.5036, "G": 0.5822, "RP": 0.7620, "J": 1.235, "H": 1.662,
       "Ks": 2.159, "W1": 3.353, "W2": 4.603}
FITBANDS = ["BP", "G", "RP", "J", "H", "W1", "W2"]
# A_band/A_V for dereddening the WISE bands (Wang & Chen 2019, Table 3)
AV_RATIO = {"W1": 0.039, "W2": 0.026}
# alpha = 0 must be EXCLUDED: the basis vector (lam/lam_Ks)^0 - 1 is identically
# zero, so the normal equation divides by zero and NaN propagates through every
# star's fit. (First run returned "0 significant absorbers" for exactly this
# reason, not because there were none.) Grey is unreachable in this
# differential parameterisation anyway -- see the module docstring.
ALPHA_GRID = np.concatenate([np.linspace(-1.0, -0.05, 20),
                             np.linspace(0.05, 6.0, 120)])


def running_pred(x1, x2, y, ok, n1=40, n2=12):
    """Median y on a 2D grid of (x1, x2), interpolated back to every star."""
    e1 = np.nanpercentile(x1[ok], np.linspace(1, 99, n1 + 1))
    e2 = np.nanpercentile(x2[ok], np.linspace(1, 99, n2 + 1))
    i1 = np.clip(np.digitize(x1, e1) - 1, 0, n1 - 1)
    i2 = np.clip(np.digitize(x2, e2) - 1, 0, n2 - 1)
    grid = np.full((n1, n2), np.nan)
    key = i1 * n2 + i2
    okk = ok & np.isfinite(y)
    df = pd.DataFrame({"k": key[okk], "y": y[okk]})
    med = df.groupby("k")["y"].median()
    cnt = df.groupby("k")["y"].size()
    good = med[cnt >= 25]
    flat = np.full(n1 * n2, np.nan)
    flat[good.index.to_numpy()] = good.to_numpy()
    grid = flat.reshape(n1, n2)
    # fill gaps along the temperature axis
    for a in range(n1):
        col = grid[a]
        m = np.isfinite(col)
        if m.sum() >= 2:
            grid[a] = np.interp(np.arange(n2), np.flatnonzero(m), col[m])
    out = grid.reshape(-1)[key]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    n = len(d)
    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    mu = d["dist_mod"].to_numpy(float)

    # dereddened absolute magnitudes in every band
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
    e1 = d["wise_w1mpro_error"].to_numpy(float)
    e2 = d["wise_w2mpro_error"].to_numpy(float)
    ok &= (e1 < 0.15) & (e2 < 0.15)
    print(f"{n:,} stars, {ok.sum():,} with all 7 bands + Ks anchor and clean WISE")

    resid = {}
    for b in FITBANDS:
        pred = running_pred(mks, jh, absmag[b], ok)
        resid[b] = absmag[b] - pred
        s = st.robust_sigma(resid[b][ok])
        print(f"  {b:3s} residual scatter {s:.4f} mag")

    R = np.column_stack([resid[b] for b in FITBANDS])
    sig = np.array([st.robust_sigma(resid[b][ok]) for b in FITBANDS])
    # centre each band on its own median so the fit sees only excursions
    R = R - np.nanmedian(R[ok], axis=0)

    lam = np.array([LAM[b] for b in FITBANDS])
    lam_ks = LAM["Ks"]
    iG = FITBANDS.index("G")

    # REPARAMETERISED. The raw basis A*[(lam/lamKs)^-alpha - 1] is degenerate:
    # as alpha -> 0 the basis vanishes, so A diverges to compensate and the fit
    # runs to the grid corner with |A| ~ 9 mag, which is physical nonsense. The
    # first run did exactly that. Normalising each basis vector by its G-band
    # entry makes the free amplitude the G-band differential deficit itself --
    # a bounded, physical quantity -- and the ratio stays finite as alpha -> 0
    # because both numerator and denominator go as alpha * ln(lam/lamKs).
    basis = np.array([(lam / lam_ks) ** (-a) - 1.0 for a in ALPHA_GRID])
    basis = basis / basis[:, iG][:, None]
    w = 1.0 / sig ** 2

    sel = np.flatnonzero(ok)
    # A full (nstar, nalpha) array is 3e6 x 140 x 8 bytes = 3.3 GB per matrix and
    # we would need three of them. Chunk over stars and keep only the argmin.
    den = (basis ** 2 * w).sum(axis=1)            # (na,)
    bw = (basis * w).T.astype(np.float32)         # (nb, na)
    alpha = np.empty(len(sel), dtype=np.float32)
    Abest = np.empty(len(sel), dtype=np.float32)
    chi2best = np.empty(len(sel), dtype=np.float32)
    dchi2 = np.empty(len(sel), dtype=np.float32)
    alpha_hi = np.empty(len(sel), dtype=np.float32)
    CH = 200_000
    for s0 in range(0, len(sel), CH):
        s1 = min(s0 + CH, len(sel))
        Rs = R[sel[s0:s1]].astype(np.float32)
        num = Rs @ bw                              # (chunk, na)
        c0 = (Rs ** 2 * w.astype(np.float32)).sum(axis=1)
        chi2 = c0[:, None] - num ** 2 / den.astype(np.float32)
        ib = np.argmin(chi2, axis=1)
        r_ = np.arange(s1 - s0)
        alpha[s0:s1] = ALPHA_GRID[ib]
        Abest[s0:s1] = num[r_, ib] / den[ib]
        chi2best[s0:s1] = chi2[r_, ib]
        dchi2[s0:s1] = c0 - chi2[r_, ib]
        # Is alpha significantly BELOW the dust value? Profile the grid: the
        # upper edge of the delta-chi2 < 4 interval (2 sigma, one parameter).
        ok_a = chi2 <= (chi2[r_, ib] + 4.0)[:, None]
        hi = np.where(ok_a.any(axis=1),
                      ALPHA_GRID[(ok_a * np.arange(len(ALPHA_GRID))).argmax(axis=1)],
                      np.nan)
        alpha_hi[s0:s1] = hi
        del Rs, num, chi2, ok_a

    out = pd.DataFrame({
        "source_id": d["source_id"].to_numpy()[sel],
        "ra": d["ra"].to_numpy()[sel], "dec": d["dec"].to_numpy()[sel],
        "residual": d["residual"].to_numpy()[sel],
        "alpha": alpha, "alpha_upper_2sig": alpha_hi,
        "deficit_G_mag": Abest, "dchi2": dchi2,
        "chi2_fit": chi2best,
        "cstar_nsigma": d["cstar_nsigma"].to_numpy()[sel],
        "ruwe": d["ruwe"].to_numpy()[sel],
        "A_0": d["A_0"].to_numpy()[sel],
    })
    # A real absorber must DIM the star in G relative to Ks, be significant,
    # and stay physically bounded. The amplitude is now the G-band deficit, so
    # requiring it below 3 mag rejects the runaway fits directly.
    strong = out[(out["deficit_G_mag"] > 0.10) & (out["deficit_G_mag"] < 3.0)
                 & (out["dchi2"] > 25)]
    print(f"\nstars with a significant, physically bounded dimming absorber "
          f"(0.1 < dm_G < 3 mag, dchi2>25): {len(strong):,}")

    print("\n=== alpha distribution of significant absorbers ===")
    bins = [-1, 0.5, 1.0, 1.5, 1.8, 2.2, 2.6, 3.5, 6.01]
    h = pd.cut(strong["alpha"], bins).value_counts().sort_index()
    tot = max(len(strong), 1)
    for iv, c in h.items():
        flag = ""
        if iv.right <= 1.5:
            flag = "  <-- NOT dust-like"
        if iv.left >= 2.6:
            flag = "  <-- steeper than dust"
        print(f"  alpha in {str(iv):16s} {c:8,}  ({100*c/tot:5.2f}%){flag}")

    print(f"\n  median alpha : {strong['alpha'].median():.3f}   "
          f"(interstellar dust ~ 2.0)")

    # The anomalous set. alpha must be SIGNIFICANTLY below the dust value --
    # the 2-sigma upper edge of its confidence interval below 1.5 -- not merely
    # best-fit low, which any noisy star can manage.
    anom = strong[(strong["alpha_upper_2sig"] < 1.5) & (strong["chi2_fit"] < 14.0)
                  & (strong["cstar_nsigma"].abs() < 1.0)
                  & (strong["ruwe"] < 1.1)]
    print(f"\n=== ANOMALOUS: alpha significantly < 1.5 (2sig), good fit, "
          f"clean photometry + astrometry ===")
    print(f"  {len(anom):,} stars")
    if len(anom):
        show = anom.nlargest(20, "deficit_G_mag")[
            ["source_id", "ra", "dec", "alpha", "alpha_upper_2sig",
             "deficit_G_mag", "chi2_fit", "dchi2", "cstar_nsigma", "ruwe"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:11.4g}"))

    out.to_parquet(cfg.DERIVED_DIR / f"spectral_slope_{args.tag}.parquet",
                   index=False)
    anom.to_csv(cfg.RESULT_DIR / f"slope_anomalies_{args.tag}.csv", index=False)
    res = {"n_fitted": int(ok.sum()), "n_significant": int(len(strong)),
           "median_alpha": float(strong["alpha"].median()) if len(strong) else None,
           "n_anomalous_lowalpha": int(len(anom)),
           "frac_alpha_below_1p5": float((strong["alpha"] < 1.5).mean())
           if len(strong) else None,
           "note": "grey (alpha=0) is invisible by construction; see scripts/37"}
    (cfg.RESULT_DIR / f"spectral_slope_{args.tag}.json").write_text(
        json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
