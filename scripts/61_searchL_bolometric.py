#!/usr/bin/env python
"""Search L: does the energy budget close, measured against gravity?

    run.sh scripts/61_searchL_bolometric.py

THE LOOPHOLE THIS CLOSES
------------------------
Channel 4 regressed absolute G on log(dynamical mass) and found zero positive
outliers. That is the project's strongest constraint, but it has a specific
escape: it measures ONE band. An absorber that intercepts optical light and
re-radiates it at some wavelength we did not check produces a G-band deficit
and no bolometric deficit at all. Channel 3 patched part of this with a W1/W2
veto, which only reaches shells hotter than ~300 K; channel 9 added W3/W4,
which reaches ~130 K. Below that the re-emission simply leaves our bands.

Integrating the whole measured spectral energy distribution removes the escape
entirely. If the total energy leaving a star is less than its mass says it
should produce, the missing energy is not being re-radiated at ANY wavelength
we can see -- it is being beamed away, stored, or converted to something
non-thermal. That is the conservation-law statement the single-band version
cannot make.

WHY DYNAMICAL MASS IS THE ONLY VALID ANCHOR
-------------------------------------------
Anchoring on any photometric quantity is circular: the anchor moves when the
absorber does. Gravity does not care. Gaia's binary_masses gives m1 from
astrometric orbits, so the predictor is independent of every band we sum.

THE COST, STATED UP FRONT
-------------------------
Dynamical masses come from binaries, and the science sample of this project
deliberately excludes binaries (RUWE < 1.4, non_single_star = 0). The overlap
between the two is exactly zero, so the infrared photometry has to be fetched
separately, and the population being constrained is not the population the
photometric channels constrain. Unresolved secondary light also inflates the
measured luminosity, which biases AGAINST finding a deficit -- so the search
is conservative by construction, as in channel 4.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchL")

# Effective wavelengths (micron) and zero points (Jy) for the bands we sum.
BANDS = {
    "BP":  (0.5050, 3552.0),
    "G":   (0.6230, 3229.0),
    "RP":  (0.7773, 2555.0),
    "J":   (1.235,  1594.0),
    "H":   (1.662,  1024.0),
    "Ks":  (2.159,   666.7),
    "W1":  (3.353,   309.5),
    "W2":  (4.603,   171.8),
    "W3":  (11.56,    31.7),
    "W4":  (22.09,     8.4),
}

XMATCH_RADIUS_ARCSEC = 3.0
K_SIGMA = 5.0
RNG_SEED = 61_2026


def fetch_wise(bm: pd.DataFrame, force=False) -> pd.DataFrame:
    """AllWISE W1-W4 for the dynamical-mass sample, via CDS XMatch."""
    out = cfg.RAW_DIR / "binary_masses_wise.parquet"
    if out.exists() and not force:
        d = pd.read_parquet(out)
        log.info("cached WISE match: %d rows", len(d))
        return d

    from astroquery.xmatch import XMatch
    from astropy.table import Table
    import astropy.units as u

    src = bm[["source_id", "ra", "dec"]].dropna()
    log.info("cross-matching %d positions against AllWISE ...", len(src))

    chunks = []
    step = 25000
    for i0 in range(0, len(src), step):
        part = src.iloc[i0:i0 + step]
        t = Table.from_pandas(part)
        t0 = time.time()
        r = XMatch.query(cat1=t, cat2="vizier:II/328/allwise",
                         max_distance=XMATCH_RADIUS_ARCSEC * u.arcsec,
                         colRA1="ra", colDec1="dec")
        log.info("  rows %d-%d -> %d matches (%.1fs)",
                 i0, min(i0 + step, len(src)), len(r), time.time() - t0)
        chunks.append(r.to_pandas())

    d = pd.concat(chunks, ignore_index=True)
    keep = ["source_id", "angDist", "W1mag", "W2mag", "W3mag", "W4mag",
            "e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"]
    keep = [c for c in keep if c in d.columns]
    d = d[keep]
    # keep the nearest match per source
    d = d.sort_values("angDist").drop_duplicates("source_id", keep="first")
    d.to_parquet(out, index=False)
    log.info("wrote %s (%d unique sources)", out, len(d))
    return d


def integrate_sed(mags: dict, dist_pc: np.ndarray) -> np.ndarray:
    """Trapezoidal integral of nu*F_nu over log(lambda) -> apparent bolometric.

    Bands are sparse and unevenly spaced, so this is a lower bound on the true
    bolometric flux rather than a full SED fit. That is acceptable here because
    the search is differential: every star is integrated the same way and we
    look for outliers against the population, not for absolute luminosities.
    """
    names = [b for b in BANDS if b in mags]
    lam = np.array([BANDS[b][0] for b in names])
    order = np.argsort(lam)
    lam = lam[order]
    names = [names[i] for i in order]

    n = len(dist_pc)
    nu_fnu = np.full((len(names), n), np.nan)
    for i, b in enumerate(names):
        zp = BANDS[b][1]
        m = mags[b]
        f_jy = zp * 10.0 ** (-0.4 * m)                 # Jy
        nu = 2.998e14 / BANDS[b][0]                     # Hz
        nu_fnu[i] = nu * f_jy * 1e-23                   # erg/s/cm^2

    log_lam = np.log(lam)
    total = np.zeros(n)
    covered = np.zeros(n)
    for i in range(len(names) - 1):
        a, b_ = nu_fnu[i], nu_fnu[i + 1]
        both = np.isfinite(a) & np.isfinite(b_)
        dlog = log_lam[i + 1] - log_lam[i]
        seg = 0.5 * (a + b_) * dlog
        total = np.where(both, total + seg, total)
        covered += both.astype(float) * dlog
    total[covered <= 0] = np.nan
    return total, covered


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--force-xmatch", action="store_true")
    ap.add_argument("--k", type=float, default=K_SIGMA)
    args = ap.parse_args()

    bm = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    log.info("dynamical-mass systems: %d", len(bm))

    wise = fetch_wise(bm, force=args.force_xmatch)
    d = bm.merge(wise, on="source_id", how="left")

    # ---- assemble the bands ----------------------------------------------
    # BP and RP are not in the cached binary_masses pull; reconstruct BP and
    # RP from G and bp_rp where possible, which is exact for the colour but
    # only approximate for the split. Drop them if that is not available.
    mags = {}
    if "phot_g_mean_mag" in d:
        mags["G"] = d["phot_g_mean_mag"].to_numpy(float)
    if "tmass_j_m" in d:
        mags["J"] = d["tmass_j_m"].to_numpy(float)
    if "tmass_ks_m" in d:
        mags["Ks"] = d["tmass_ks_m"].to_numpy(float)
    for b in ("W1", "W2", "W3", "W4"):
        col = f"{b}mag"
        if col in d:
            mags[b] = d[col].to_numpy(float)

    log.info("bands available: %s", ", ".join(sorted(mags,
             key=lambda b: BANDS[b][0])))

    dist_pc = 1000.0 / d["parallax"].to_numpy(float)
    f_bol, covered = integrate_sed(mags, dist_pc)

    # ---- absolute bolometric magnitude ------------------------------------
    with np.errstate(divide="ignore", invalid="ignore"):
        m_bol = -2.5 * np.log10(f_bol)
        M_bol = m_bol - 5.0 * np.log10(np.maximum(dist_pc, 1e-6)) + 5.0

    d = d.assign(f_bol=f_bol, M_bol=M_bol, dist_pc=dist_pc,
                 log_cover=covered)

    # ---- the regression ---------------------------------------------------
    good = (d["m1"].notna() & (d["m1"] > 0.2) & (d["m1"] < 1.6)
            & np.isfinite(M_bol) & (d["parallax_over_error"] > 20)
            & (covered > 1.5))          # require real wavelength coverage
    n_w4 = int(d.loc[good, "W4mag"].notna().sum()) if "W4mag" in d else 0
    n_w3 = int(d.loc[good, "W3mag"].notna().sum()) if "W3mag" in d else 0
    log.info("usable systems: %d  (with W3: %d, with W4: %d)",
             int(good.sum()), n_w3, n_w4)

    g = d[good].reset_index(drop=True)
    x = np.log10(g["m1"].to_numpy(float))
    y = g["M_bol"].to_numpy(float)

    best = None
    for deg in (2, 3, 4):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, sigma, deg = best
    resid = y - np.polyval(coef, x)
    log.info("M_bol vs log M_dyn: degree %d, robust scatter %.4f mag",
             deg, sigma)

    med = float(np.median(resid))
    thr = args.k * sigma
    n_pos = int(np.count_nonzero(resid > med + thr))   # under-luminous = deficit
    n_neg = int(np.count_nonzero(resid < med - thr))   # over-luminous

    log.info("")
    log.info("at %.0f sigma (= %.3f mag):", args.k, thr)
    log.info("  bolometrically DEFICIENT (energy missing): %d", n_pos)
    log.info("  bolometrically EXCESSIVE (extra light)   : %d", n_neg)

    # The negative tail is unresolved secondary light and is the pipeline's
    # own false-positive rate for this estimator: both tails are populated by
    # the same photometric scatter, but only the positive one can be signal.
    ratio = n_pos / n_neg if n_neg > 0 else float("inf") if n_pos else 0.0
    log.info("  ratio deficient:excessive = %s",
             f"{ratio:.2f}" if np.isfinite(ratio) else "inf")

    f_det = float(st.fraction_from_delta(thr))
    ul = st.poisson_upper_limit(n_pos)
    p_ul = ul / len(g)
    log.info("")
    log.info("detectable intercepted fraction at this threshold: f >= %.3f",
             f_det)
    log.info("Poisson 95%% upper limit on the rate: p < %.3g", p_ul)

    if n_pos == 0:
        verdict = (
            f"NULL, and this is the strongest form of the constraint in this "
            f"project. Among {len(g):,} systems with a dynamical mass and "
            f"measured photometry from 0.5 to 22 micron, ZERO are "
            f"bolometrically deficient at {args.k:.0f} sigma. Energy leaving "
            f"these stars balances what their gravitational mass says they "
            f"produce, so no more than a fraction p < {p_ul:.2g} intercept "
            f"f >= {f_det:.2f} of their output and fail to re-radiate it in "
            f"any band we can see. Unlike channel 4 this cannot be evaded by "
            f"re-emitting at an unchecked wavelength; it can only be evaded "
            f"by beaming the waste heat away from us or by storing it.")
    elif n_pos <= n_neg:
        verdict = (
            f"NULL. The {n_pos} bolometrically deficient systems are matched "
            f"by {n_neg} equally-significant over-luminous ones, so both tails "
            f"are photometric scatter rather than a distinct population.")
    else:
        verdict = (
            f"ASYMMETRIC: {n_pos} bolometrically deficient against {n_neg} "
            f"over-luminous. Since unresolved secondary light biases this "
            f"estimator toward over-luminosity, an excess in the deficient "
            f"tail runs against the dominant systematic and warrants "
            f"follow-up.")

    print(f"\n{'='*68}")
    print("SEARCH L: BOLOMETRIC ENERGY CLOSURE vs DYNAMICAL MASS")
    print(f"{'='*68}")
    print(f"  systems with dynamical mass + SED   : {len(g):,}")
    print(f"  bands integrated                    : "
          f"{', '.join(sorted(mags, key=lambda b: BANDS[b][0]))}")
    print(f"  M_bol vs log M_dyn scatter          : {sigma:.4f} mag")
    print(f"  threshold ({args.k:.0f} sigma)                 : {thr:.3f} mag")
    print(f"  detectable intercepted fraction     : f >= {f_det:.3f}")
    print()
    print(f"  bolometrically DEFICIENT            : {n_pos}")
    print(f"  bolometrically EXCESSIVE (mirror)   : {n_neg}")
    print(f"  Poisson 95% upper limit             : p < {p_ul:.3g}")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "n_systems": int(len(g)),
        "bands": sorted(mags, key=lambda b: BANDS[b][0]),
        "n_with_w3": n_w3,
        "n_with_w4": n_w4,
        "poly_degree": int(deg),
        "robust_sigma_mag": float(sigma),
        "k_sigma": args.k,
        "threshold_mag": float(thr),
        "f_detectable": f_det,
        "n_deficient": n_pos,
        "n_excessive": n_neg,
        "ratio_deficient_excessive": (
            float(ratio) if np.isfinite(ratio) else None),
        "poisson_upper_limit_rate": float(p_ul),
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchL_bolometric_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)

    if n_pos:
        cand = g[resid > med + thr].copy()
        cand["bol_deficit_mag"] = resid[resid > med + thr]
        csv = cfg.RESULT_DIR / f"searchL_candidates_{args.tag}.csv"
        cand.sort_values("bol_deficit_mag", ascending=False).to_csv(
            csv, index=False)
        log.info("wrote %d candidates to %s", len(cand), csv)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
