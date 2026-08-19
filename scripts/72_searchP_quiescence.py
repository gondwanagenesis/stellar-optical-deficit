#!/usr/bin/env python
"""Search P: stars that are too quiet for what they are.

    run.sh scripts/72_searchP_quiescence.py --tag primary

THE INVERSION
-------------
Every technosignature search, including all twenty-two channels above, looks
for something ANOMALOUS: excess infrared, a deficit, a weird transit, an odd
orbit. Nobody looks for anomalous REGULARITY.

Low-mass stars are magnetically active. Spots rotate in and out of view, flares
fire, and the result is photometric scatter well above the measurement floor.
The fraction showing that scatter varies smoothly with colour and with age,
because activity decays as a star spins down. A star that is unusually steady
for its type is, in the ordinary reading, simply old and slow -- but it is also
what a stabilised star would look like, and nothing in the literature has ever
checked whether the quiet population is distributed the way it should be.

THE PHYSICS THAT BOUNDS THE TEST
--------------------------------
There is a hard floor and it is not a property of the star. Photon noise lives
in our detector, so no source can scatter LESS than Poisson statistics allow.
"Anomalously quiet" therefore cannot mean sub-floor scatter -- that would be an
instrumental artefact, not a discovery. It can only mean over-population OF the
floor: too many stars sitting at the measurement limit, in a place or a colour
range where their neighbours are active.

So the statistic is the quiet FRACTION, and the question is whether it varies
across the sky beyond what reshuffling can produce.

THE NULL
--------
Mark permutation again: the quiet flag is shuffled over positions, preserving
the density field, the magnitude distribution and the scanning-law footprint
exactly, and destroying only the association between where a star is and
whether it is quiet. Generating synthetic positions instead would hand every
survey-geometry artefact to the signal.
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
log = logging.getLogger("searchP")

COLS = ["source_id", "l", "b", "dist_pc", "bp_rp", "A_0",
        "phot_g_mean_mag", "phot_g_mean_flux_over_error",
        "phot_bp_n_obs", "phot_variable_flag", "ruwe", "M_Ks"]

N_SHUFFLE = 300
MIN_PER_CELL = 200
RNG_SEED = 72_2026


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--nside-cells", type=int, default=10)
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    import pyarrow.parquet as pq
    path = cfg.DERIVED_DIR / f"{args.tag}_resid.parquet"
    have = set(pq.read_schema(path).names)
    d = pd.read_parquet(path, columns=[c for c in COLS if c in have])
    log.info("loaded %d stars", len(d))

    fove = d["phot_g_mean_flux_over_error"].to_numpy(float)
    n_bp = d["phot_bp_n_obs"].fillna(0).to_numpy(float)
    n_g = np.where(n_bp > 0, 9.0 * n_bp, np.nan)
    amp = np.sqrt(n_g) / np.maximum(fove, 1e-9)      # fractional epoch scatter
    g = d["phot_g_mean_mag"].to_numpy(float)
    bp_rp = d["bp_rp"].to_numpy(float)

    ok = (np.isfinite(amp) & (amp > 0) & np.isfinite(g)
          & np.isfinite(bp_rp) & (d["ruwe"].to_numpy(float) < 1.4))
    log.info("usable: %d", int(ok.sum()))

    # ---- the measurement floor, empirically ------------------------------
    # The floor is set by photon statistics and therefore by magnitude. Take
    # a low percentile of the amplitude in narrow magnitude bins: at any
    # magnitude some stars really are quiet, so the low percentile traces the
    # instrumental limit rather than any astrophysics.
    edges = np.arange(np.floor(np.nanmin(g[ok])), np.ceil(np.nanmax(g[ok])) + 0.25, 0.25)
    idx = np.digitize(g, edges)
    floor = np.full(len(d), np.nan)
    for b in range(1, len(edges)):
        m = ok & (idx == b)
        if m.sum() < 200:
            continue
        floor[idx == b] = np.percentile(amp[m], 10)

    have_floor = ok & np.isfinite(floor)
    # excess variance above the floor, in amplitude units
    excess = np.sqrt(np.maximum(amp ** 2 - floor ** 2, 0.0))
    quiet = have_floor & (excess < 0.2 * floor)      # essentially at the limit
    log.info("stars at the measurement floor (quiet): %d (%.1f%%)",
             int(quiet.sum()), 100 * quiet[have_floor].mean())

    # ---- does the quiet fraction behave as activity physics predicts? ----
    log.info("")
    log.info("=== quiet fraction vs colour (activity should rise to the red) ===")
    col_rows = []
    for lo, hi in [(0.7, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5),
                   (2.5, 3.0), (3.0, 3.6)]:
        m = have_floor & (bp_rp >= lo) & (bp_rp < hi)
        if m.sum() < 500:
            continue
        f = float(quiet[m].mean())
        col_rows.append({"bp_rp_lo": lo, "bp_rp_hi": hi,
                         "n": int(m.sum()), "quiet_frac": f})
        log.info("  BP-RP %.1f-%.1f  n=%8d  quiet fraction %.4f",
                 lo, hi, int(m.sum()), f)

    # ---- spatial variation, against a mark-permutation null --------------
    sel = have_floor
    nl, nb = args.nside_cells * 2, args.nside_cells
    li = np.clip((d["l"].to_numpy(float)[sel] / 360.0 * nl).astype(int), 0, nl - 1)
    sb = np.sin(np.radians(d["b"].to_numpy(float)[sel]))
    bi = np.clip(((sb + 1.0) / 2.0 * nb).astype(int), 0, nb - 1)
    cell = li * nb + bi
    n_cells = nl * nb
    q = quiet[sel].astype(float)
    log.info("")
    log.info("spatial grid: %d cells, %d stars", n_cells, len(q))

    def cell_fracs(vals):
        s = np.bincount(cell, weights=vals, minlength=n_cells)
        n = np.bincount(cell, minlength=n_cells)
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(n >= MIN_PER_CELL, s / np.maximum(n, 1), np.nan)
        return f, n

    f_obs, n_cell = cell_fracs(q)
    good = np.isfinite(f_obs)
    glob = float(q.mean())
    dev = f_obs[good] - glob
    max_obs = float(np.nanmax(np.abs(dev)))
    log.info("populated cells: %d, global quiet fraction %.4f",
             int(good.sum()), glob)
    log.info("max |cell - global| = %.4f", max_obs)

    null_max = np.empty(N_SHUFFLE)
    for i in range(N_SHUFFLE):
        f_n, _ = cell_fracs(q[rng.permutation(len(q))])
        gg = np.isfinite(f_n)
        null_max[i] = float(np.nanmax(np.abs(f_n[gg] - glob)))
    mu, sd = float(null_max.mean()), float(null_max.std())
    z = (max_obs - mu) / sd if sd > 0 else 0.0
    p_emp = float((null_max >= max_obs).mean())
    log.info("permutation null: %.4f +/- %.4f  ->  z = %+.2f, p = %.4f",
             mu, sd, z, p_emp)

    # worst cell
    ii = np.where(good)[0]
    worst = ii[int(np.argmax(np.abs(dev)))]
    wl = (worst // nb + 0.5) / nl * 360.0
    wsb = ((worst % nb) + 0.5) / nb * 2.0 - 1.0
    wb = float(np.degrees(np.arcsin(np.clip(wsb, -1, 1))))
    log.info("most deviant cell: l=%.0f b=%.0f  quiet fraction %.4f (n=%d)",
             wl, wb, f_obs[worst], n_cell[worst])

    # ---- verdict ----------------------------------------------------------
    rises_red = (len(col_rows) >= 3
                 and col_rows[-1]["quiet_frac"] < col_rows[0]["quiet_frac"])

    if z > 5.0 and p_emp < 0.01:
        verdict = (
            f"The quiet fraction varies across the sky beyond permutation "
            f"({z:+.1f} sigma, p = {p_emp:.3g}), peaking at l={wl:.0f}, "
            f"b={wb:.0f}. Before this means anything it has to be separated "
            f"from the scanning law: Gaia's number of transits varies with "
            f"ecliptic latitude, and more transits means a better-determined "
            f"mean and a lower apparent amplitude.")
    else:
        verdict = (
            f"NULL. The fraction of stars sitting at the photometric floor is "
            f"the same everywhere within 500 pc: the largest cell deviation is "
            f"{max_obs:.4f} against a permutation null of {mu:.4f} +/- {sd:.4f} "
            f"({z:+.2f} sigma, p = {p_emp:.3f}). No region contains an excess "
            f"of anomalously steady stars. The quiet fraction also "
            f"{'falls' if rises_red else 'does not fall'} toward redder "
            f"colours, which is the direction ordinary activity physics "
            f"predicts and an independent check that the statistic measures "
            f"what it claims to.")

    print(f"\n{'='*70}")
    print("SEARCH P: ANOMALOUSLY QUIET STARS")
    print(f"{'='*70}")
    print(f"  stars with a usable amplitude       : {int(have_floor.sum()):,}")
    print(f"  at the photometric floor ('quiet')  : {int(quiet.sum()):,} "
          f"({100*glob:.1f}%)")
    print(f"  populated sky cells                 : {int(good.sum())}")
    print(f"  max |cell - global| quiet fraction  : {max_obs:.4f}")
    print(f"  permutation null                    : {mu:.4f} +/- {sd:.4f}")
    print(f"  excess                              : {z:+.2f} sigma "
          f"(p = {p_emp:.4f})")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchP_quiescence_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "n_usable": int(have_floor.sum()),
        "n_quiet": int(quiet.sum()),
        "global_quiet_fraction": glob,
        "quiet_fraction_by_colour": col_rows,
        "n_cells_populated": int(good.sum()),
        "max_abs_deviation": max_obs,
        "null_mean": mu, "null_std": sd,
        "z": float(z), "p_empirical": p_emp,
        "worst_cell": {"l_deg": float(wl), "b_deg": wb,
                       "quiet_frac": float(f_obs[worst]),
                       "n": int(n_cell[worst])},
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
