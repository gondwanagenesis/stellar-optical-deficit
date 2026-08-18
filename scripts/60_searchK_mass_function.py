#!/usr/bin/env python
"""Search K: holes in the stellar population, not dim stars.

    run.sh scripts/60_searchK_mass_function.py --tag primary

THE BLIND SPOT THIS TARGETS
---------------------------
Every channel in this project so far examines stars that are IN the catalogue
and asks whether their light is anomalous. That is structurally blind to the
limiting case it was built to find. A star intercepted at f = 0.5 is 0.75 mag
faint and sits in our residual tail; a star intercepted at f -> 1 is not faint,
it is ABSENT. It falls below the magnitude limit, fails the quality cuts, and
leaves the sample entirely.

So the signature of complete harvesting is not a dim star. It is a hole.

WHY THE MASS FUNCTION RATHER THAN THE NUMBER DENSITY
----------------------------------------------------
Counting stars per unit volume is the obvious test and the wrong one: stellar
density has enormous natural gradients from Galactic structure, and it is
degraded by extinction, crowding and Gaia's scanning-law-dependent
completeness. A deficit means nothing against that background.

The SHAPE of the mass function is a far better observable, because it is close
to universal. Whatever sets the initial mass function does not vary much from
place to place in the thin disc, so a spatial cell whose mass-function slope
differs from the global one is anomalous in a way a low-density cell is not.

And the shape carries the physics: harvesting is not mass-blind. A civilisation
choosing targets would choose by stellar type -- long-lived K and M dwarfs for
persistence, or bright F and G stars for power -- and either preference removes
stars from one end of the mass function and tilts the slope.

THE NULL
--------
Masses are permuted over positions. This preserves the spatial selection
function, the density field, the extinction pattern and the global mass
distribution exactly, and destroys only the association between where a star
is and how massive it is. Drawing synthetic positions instead would attribute
every survey-geometry artefact to the signal.

THE DISCRIMINANT AGAINST DUST
-----------------------------
Extinction also removes stars, preferentially the faint red ones, which tilts
the slope in a cell. But extinction is measured per star here, so the test is
repeated in low-extinction cells only: a dust-driven tilt weakens, an
engineering-driven one does not.
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
log = logging.getLogger("searchK")

MAX_DIST_PC = 300.0        # completeness is defensible here
MASS_LO, MASS_HI = 0.30, 0.80    # the range our sample actually populates
MIN_PER_CELL = 120
N_SHUFFLE = 300
RNG_SEED = 60_2026


def fit_mass_relation(bm: pd.DataFrame):
    """log10(mass) as a cubic in M_Ks, from Gaia's own dynamical masses."""
    d = bm[bm["m1"].notna() & bm["tmass_ks_m"].notna()
           & (bm["parallax"] > 0)]
    dm = 5.0 * np.log10(1000.0 / d["parallax"].to_numpy(float)) - 5.0
    m_ks = d["tmass_ks_m"].to_numpy(float) - dm
    mass = d["m1"].to_numpy(float)
    ok = (np.isfinite(m_ks) & np.isfinite(mass) & (m_ks > 1.0) & (m_ks < 9.5)
          & (mass > 0.1) & (mass < 2.0))
    coef = np.polyfit(m_ks[ok], np.log10(mass[ok]), 3)
    return coef, int(ok.sum())


def slope_mle(masses, lo, hi):
    """Maximum-likelihood power-law index for dN/dM ~ M^-alpha on [lo, hi].

    The truncated-Pareto MLE has no closed form, so solve the score equation
    by bisection. Using the MLE rather than a binned fit matters here: cells
    hold a few hundred stars and binning throws away most of the information.
    """
    m = masses[(masses >= lo) & (masses <= hi)]
    n = len(m)
    if n < 10:
        return np.nan, n
    s = np.log(m).sum()
    ln_lo, ln_hi = np.log(lo), np.log(hi)

    # Solve by matching the observed mean of log(m) to its expectation under
    # the truncated Pareto, which is monotonic in alpha and so bisects safely.
    target = s / n

    def mean_logm(a):
        if abs(a - 1.0) < 1e-6:
            a = 1.0 + 1e-6
        p = 1.0 - a
        num = (hi ** p * ln_hi - lo ** p * ln_lo) / p - \
              (hi ** p - lo ** p) / p ** 2
        den = (hi ** p - lo ** p) / p
        return num / den

    lo_a, hi_a = -6.0, 12.0
    if (mean_logm(lo_a) - target) * (mean_logm(hi_a) - target) > 0:
        return np.nan, n
    for _ in range(80):
        mid = 0.5 * (lo_a + hi_a)
        if (mean_logm(lo_a) - target) * (mean_logm(mid) - target) <= 0:
            hi_a = mid
        else:
            lo_a = mid
    return 0.5 * (lo_a + hi_a), n


def cell_slopes(cell_id, masses, n_cells):
    """Slope per cell; NaN where the cell is too sparse."""
    out = np.full(n_cells, np.nan)
    counts = np.zeros(n_cells, dtype=int)
    order = np.argsort(cell_id)
    cid_s = cell_id[order]
    m_s = masses[order]
    edges = np.searchsorted(cid_s, np.arange(n_cells + 1))
    for c in range(n_cells):
        i0, i1 = edges[c], edges[c + 1]
        if i1 - i0 < MIN_PER_CELL:
            counts[c] = i1 - i0
            continue
        a, n = slope_mle(m_s[i0:i1], MASS_LO, MASS_HI)
        out[c] = a
        counts[c] = n
    return out, counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--nside-cells", type=int, default=12,
                    help="cells per axis on the sky grid (l x b)")
    ap.add_argument("--max-a0", type=float, default=None,
                    help="keep only stars below this extinction, and rebuild "
                         "the null on that subsample. Comparing a "
                         "low-extinction maximum against a null computed on "
                         "the full grid is not a fair test: the expected "
                         "maximum grows with the number of cells.")
    ap.add_argument("--min-absb", type=float, default=None,
                    help="keep only stars above this |b|")
    ap.add_argument("--max-dist", type=float, default=MAX_DIST_PC,
                    help="distance limit. The dominant systematic is "
                         "magnitude-limited completeness: at fixed distance a "
                         "0.3 Msun star is far fainter than a 0.8 Msun one, so "
                         "any cell-to-cell variation in effective survey depth "
                         "tilts the fitted slope. Shrinking this until the "
                         "sample is volume-complete for the whole mass window "
                         "is the test.")
    args = ap.parse_args()
    max_dist = args.max_dist

    rng = np.random.default_rng(RNG_SEED)

    bm = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    coef, n_cal = fit_mass_relation(bm)
    log.info("mass scale from %d Gaia dynamical masses", n_cal)

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["source_id", "l", "b", "dist_pc",
                                 "M_Ks", "A_0", "phot_g_mean_mag"])
    log.info("sample: %d", len(d))

    ok = (d["dist_pc"] > 0) & (d["dist_pc"] < max_dist) & d["M_Ks"].notna()
    if args.max_a0 is not None:
        ok &= d["A_0"] < args.max_a0
    if args.min_absb is not None:
        ok &= d["b"].abs() > args.min_absb
    d = d[ok].reset_index(drop=True)
    log.info("within %.0f pc with good M_Ks%s%s: %d", max_dist,
             f", A_0 < {args.max_a0}" if args.max_a0 is not None else "",
             f", |b| > {args.min_absb}" if args.min_absb is not None else "",
             len(d))

    mass = 10.0 ** np.polyval(coef, np.clip(d["M_Ks"].to_numpy(float), 1.0, 9.5))
    d = d.assign(mass=mass)

    in_range = (mass >= MASS_LO) & (mass <= MASS_HI)
    d = d[in_range].reset_index(drop=True)
    mass = d["mass"].to_numpy(float)
    log.info("in the mass window %.2f-%.2f Msun: %d", MASS_LO, MASS_HI, len(d))

    global_alpha, _ = slope_mle(mass, MASS_LO, MASS_HI)
    log.info("global mass-function index alpha = %.3f", global_alpha)

    # ---- sky cells (equal-area in l, sin b) ------------------------------
    nl = args.nside_cells * 2
    nb = args.nside_cells
    li = np.clip((d["l"].to_numpy(float) / 360.0 * nl).astype(int), 0, nl - 1)
    sb = np.sin(np.radians(d["b"].to_numpy(float)))
    bi = np.clip(((sb + 1.0) / 2.0 * nb).astype(int), 0, nb - 1)
    cell = li * nb + bi
    n_cells = nl * nb
    log.info("sky grid: %d x %d = %d cells", nl, nb, n_cells)

    alphas, counts = cell_slopes(cell, mass, n_cells)
    good = np.isfinite(alphas)
    log.info("cells with >= %d stars: %d of %d",
             MIN_PER_CELL, int(good.sum()), n_cells)

    if good.sum() < 10:
        log.error("too few populated cells")
        return 1

    dev = alphas[good] - global_alpha
    log.info("cell-to-cell scatter in alpha: %.4f (rms)", float(np.std(dev)))

    # ---- mark-permutation null -------------------------------------------
    max_abs_obs = float(np.nanmax(np.abs(dev)))
    null_max = np.empty(N_SHUFFLE)
    for i in range(N_SHUFFLE):
        shuffled = mass[rng.permutation(len(mass))]
        a_n, _ = cell_slopes(cell, shuffled, n_cells)
        g = np.isfinite(a_n)
        null_max[i] = float(np.nanmax(np.abs(a_n[g] - global_alpha)))
        if (i + 1) % 100 == 0:
            log.info("  null %d/%d", i + 1, N_SHUFFLE)

    mu, sd = float(null_max.mean()), float(null_max.std())
    z = (max_abs_obs - mu) / sd if sd > 0 else 0.0
    p_emp = float((null_max >= max_abs_obs).mean())

    log.info("")
    log.info("max |alpha - alpha_global| observed : %.4f", max_abs_obs)
    log.info("shuffled null                       : %.4f +/- %.4f", mu, sd)
    log.info("excess                              : %+.2f sigma", z)
    log.info("empirical p                         : %.4f", p_emp)

    # ---- the extinction discriminant -------------------------------------
    # A dust-driven tilt should weaken when only low-extinction cells are kept.
    a0_cell = np.full(n_cells, np.nan)
    for c in np.unique(cell):
        a0_cell[c] = float(np.median(d["A_0"].to_numpy(float)[cell == c]))
    low_ext = good & np.isfinite(a0_cell) & (a0_cell < 0.05)
    if low_ext.sum() >= 5:
        dev_low = alphas[low_ext] - global_alpha
        max_low = float(np.nanmax(np.abs(dev_low)))
        log.info("low-extinction cells (A_0 < 0.05): %d, max |dev| = %.4f",
                 int(low_ext.sum()), max_low)
    else:
        max_low = float("nan")
        log.info("too few low-extinction cells for the dust discriminant")

    corr_alpha_a0 = float(pd.Series(alphas[good]).corr(
        pd.Series(a0_cell[good])))
    log.info("corr(alpha, cell extinction) = %+.3f", corr_alpha_a0)

    # ---- worst cell ------------------------------------------------------
    idx = np.where(good)[0]
    worst = idx[int(np.argmax(np.abs(dev)))]
    wl = (worst // nb + 0.5) / nl * 360.0
    wsb = ((worst % nb) + 0.5) / nb * 2.0 - 1.0
    wb = np.degrees(np.arcsin(np.clip(wsb, -1, 1)))
    log.info("most deviant cell: l=%.0f b=%.0f, alpha=%.3f, N=%d, A_0=%.3f",
             wl, wb, alphas[worst], counts[worst], a0_cell[worst])

    # ---- verdict ----------------------------------------------------------
    if z > 5.0 and p_emp < 0.01:
        verdict = (
            f"The mass-function slope varies across the sky by more than "
            f"permutation allows ({z:+.1f} sigma, p = {p_emp:.3g}). Before "
            f"this means anything, the extinction correlation is "
            f"{corr_alpha_a0:+.3f} and the low-extinction maximum is "
            f"{max_low:.4f} against {max_abs_obs:.4f} overall -- if the signal "
            f"survives there it is not dust.")
    else:
        verdict = (
            f"NULL. The stellar mass function has the same shape everywhere "
            f"within {max_dist:.0f} pc: the largest cell-to-cell deviation "
            f"is {max_abs_obs:.4f} against a permutation null of "
            f"{mu:.4f} +/- {sd:.4f} ({z:+.2f} sigma, p = {p_emp:.3f}). No "
            f"region is missing stars of any particular mass. This constrains "
            f"the f -> 1 case that every photometric channel is blind to, "
            f"because a completely enshrouded star leaves the catalogue and "
            f"would deplete its cell's mass function.")

    print(f"\n{'='*68}")
    print("SEARCH K: SPATIAL ANOMALIES IN THE STELLAR MASS FUNCTION")
    print(f"{'='*68}")
    print(f"  stars within {max_dist:.0f} pc, {MASS_LO}-{MASS_HI} Msun : {len(d):,}")
    print(f"  populated sky cells (>= {MIN_PER_CELL})        : {int(good.sum())} of {n_cells}")
    print(f"  global mass-function index alpha        : {global_alpha:.3f}")
    print(f"  cell-to-cell rms in alpha               : {float(np.std(dev)):.4f}")
    print(f"  max |deviation|                         : {max_abs_obs:.4f}")
    print(f"  permutation null                        : {mu:.4f} +/- {sd:.4f}")
    print(f"  excess                                  : {z:+.2f} sigma")
    print(f"  empirical p                             : {p_emp:.4f}")
    print(f"  corr(alpha, extinction)                 : {corr_alpha_a0:+.3f}")
    print(f"  max |dev| in low-extinction cells       : {max_low:.4f}")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "max_dist_pc": max_dist,
        "mass_window": [MASS_LO, MASS_HI],
        "n_stars": int(len(d)),
        "n_cells": int(n_cells),
        "n_populated_cells": int(good.sum()),
        "min_per_cell": MIN_PER_CELL,
        "global_alpha": float(global_alpha),
        "cell_alpha_rms": float(np.std(dev)),
        "max_abs_deviation": max_abs_obs,
        "null_max_mean": mu,
        "null_max_std": sd,
        "z": float(z),
        "p_empirical": p_emp,
        "corr_alpha_extinction": corr_alpha_a0,
        "max_dev_low_extinction": max_low,
        "worst_cell": {"l_deg": float(wl), "b_deg": float(wb),
                       "alpha": float(alphas[worst]),
                       "n": int(counts[worst]),
                       "A_0": float(a0_cell[worst])},
        "n_shuffle": N_SHUFFLE,
        "verdict": verdict,
    }
    summary["cuts"] = {"max_a0": args.max_a0, "min_absb": args.min_absb}
    suffix = ""
    if args.max_a0 is not None:
        suffix += f"_a0lt{args.max_a0:g}"
    if args.min_absb is not None:
        suffix += f"_blt{args.min_absb:g}"
    out = cfg.RESULT_DIR / f"searchK_mass_function_{args.tag}{suffix}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
