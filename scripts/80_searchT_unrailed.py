#!/usr/bin/env python
"""Search T, third pass: unrail the fit at BOTH edges, then ask what is left.

    run.sh scripts/80_searchT_unrailed.py --tag primary

WHAT WENT WRONG TWICE
---------------------
Script 74 capped beta at 2.99999 and then defined its unphysical mirror as
beta > 3.0, so the mirror was identically zero for any input. Script 78 fixed
that by extending the grid to beta = 5 and making the mirror symmetric about
dust -- and got 1,023 cold candidates against 233 mirror, a 4.4x asymmetry,
which script 79 diagnosed as CMB confusion on the strength of the fitted
temperatures clustering at 2.75 K and a 56x plane-to-pole rate gradient.

That diagnosis was built on unconstrained fits. 88.6% of those 1,023
candidates sit at beta = -0.50 EXACTLY -- the lower grid edge, which script 78
inherited from script 74 and never questioned because the bug it was hunting
was at the other end. A fit pinned to a boundary is not a measurement of
anything; it is the model announcing it cannot describe the data. The fitted
temperature clustering that the CMB argument rested on is then a nuisance
parameter absorbing a slope the grid would not let beta reach.

The general lesson, which is the reusable part: extending a grid at one end
does not validate the other end. Check for boundary occupancy at EVERY edge of
every fit, as a matter of routine, and treat a railed fit as missing data.

WHAT THE LOW EDGE WAS HIDING
----------------------------
beta below -0.5 is not exotic. In this parameterisation the flux density goes
as nu^(2+beta) in the Rayleigh-Jeans limit, so:

    beta = +1.6   dust, S ~ nu^3.6
    beta =  0.0   cold blackbody, S ~ nu^2
    beta = -2.0   FLAT spectrum, S ~ nu^0
    beta = -2.7   S ~ nu^-0.7, the canonical optically-thin synchrotron slope

The high-latitude millimetre sky is dominated by blazars and flat-spectrum
radio quasars, which live in exactly that range and which a grid floored at
-0.5 cannot represent. They are extragalactic and therefore isotropic, so they
concentrate at high |b| purely because Galactic confusion suppresses detection
in the plane -- which reproduces the same 56x gradient that was read as CMB.

So this pass runs beta over [-4.0, +7.2], symmetric about dust at 1.6 so that
both tails of the mirror have identical grid room, requires the best fit to lie
strictly INSIDE the grid on both axes, and then re-asks every question script
79 asked.

The two hypotheses now make opposite predictions, which is the point:

    CMB confusion     interior fits pile at beta ~ 0, T ~ 2.7 K
    radio sources     interior fits pile at beta ~ -2, T unconstrained/hot
"""
from __future__ import annotations

import argparse
import importlib.util
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
log = logging.getLogger("searchT3")

T_CMB = 2.72548
BETA_DUST, BETA_OFFSET, T_MAX, DCHI2 = 1.6, 1.1, 10.0, 4.0
BETA_LO, BETA_HI = -4.0, 7.2          # symmetric about BETA_DUST: 1.6 -/+ 5.6
T_GRID_LO, T_GRID_HI = 2.0, 39.5


def _load(name, fname):
    p = Path(__file__).resolve().parent / fname
    sp = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    s78 = _load("s78", "78_searchT_pccs2e_fixed.py")
    s78.BETAS = np.arange(BETA_LO, BETA_HI + 1e-9, 0.05)
    log.info("beta grid: %.2f to %.2f, %d points (symmetric about dust %.1f)",
             s78.BETAS.min(), s78.BETAS.max(), len(s78.BETAS), BETA_DUST)

    d = s78.assemble_seds(s78.load_pccs2e())
    groups = s78.group_by_position(d)
    s = s78.run_fits(d, groups)
    log.info("sources with >= 3 bands: %d", len(s))

    beta, T, ab = (s["beta"].to_numpy(float), s["T"].to_numpy(float),
                   np.abs(s["glat"].to_numpy(float)))

    # ---- boundary occupancy at EVERY edge --------------------------------
    edges = {
        "beta_floor": float(np.mean(beta <= BETA_LO + 1e-9)),
        "beta_ceiling": float(np.mean(beta >= BETA_HI - 1e-9)),
        "T_floor": float(np.mean(T <= T_GRID_LO + 1e-9)),
        "T_ceiling": float(np.mean(T >= T_GRID_HI - 1e-9)),
    }
    log.info("")
    log.info("boundary occupancy on the widened grid:")
    for k, v in edges.items():
        log.info("  %-14s %.2f%%", k, 100 * v)

    interior = ((beta > BETA_LO + 1e-9) & (beta < BETA_HI - 1e-9)
                & (T > T_GRID_LO + 1e-9) & (T < T_GRID_HI - 1e-9))
    log.info("interior (usable) fits: %d / %d = %.1f%%",
             int(interior.sum()), len(s), 100 * interior.mean())

    # ---- what the previously-railed population actually is ---------------
    prev_railed = beta < -0.5 - 1e-9      # unreachable on script 78's grid
    log.info("")
    log.info("sources now fitting below script 78's floor (beta < -0.5): %d",
             int(prev_railed.sum()))
    if prev_railed.any():
        bb = beta[prev_railed & interior]
        log.info("  their interior beta: median %.2f, 16-84%% [%.2f, %.2f]",
                 float(np.median(bb)), float(np.quantile(bb, 0.16)),
                 float(np.quantile(bb, 0.84)))
        log.info("  fraction within 0.3 of beta = -2 (flat spectrum): %.2f",
                 float(np.mean(np.abs(bb + 2.0) < 0.3)))
        log.info("  their |b|<10 fraction: %.3f (all sources %.3f)",
                 float(np.mean(ab[prev_railed] < 10)), float(np.mean(ab < 10)))

    # ---- selection, interior only ----------------------------------------
    prefers_bb = (s["chi2_blackbody"] < s["chi2_dust"] - DCHI2).to_numpy()
    prefers_best = (s["chi2"] < s["chi2_dust"] - DCHI2).to_numpy()
    cold = interior & (np.abs(beta - 0.0) < BETA_OFFSET) & (T < T_MAX) & prefers_bb
    mirror = (interior & (np.abs(beta - 2 * BETA_DUST) < BETA_OFFSET)
              & (T < T_MAX) & prefers_best)
    log.info("")
    log.info("INTERIOR-ONLY, beta within %.1f of 0.0 (cold) or %.1f (mirror), "
             "T < %.0f K, dchi2 > %.0f", BETA_OFFSET, 2 * BETA_DUST, T_MAX, DCHI2)
    log.info("  cold   : %d", int(cold.sum()))
    log.info("  mirror : %d", int(mirror.sum()))

    ratio = (float(cold.sum()) / float(mirror.sum())
             if mirror.sum() else float("inf"))

    # ---- re-ask script 79's CMB questions on interior fits ---------------
    if cold.sum():
        near = float(np.mean(np.abs(T[cold] - T_CMB) <= 0.5))
        medT = float(np.median(T[cold]))
    else:
        near, medT = float("nan"), float("nan")
    log.info("")
    log.info("cold tail (interior): median T %.2f K, %.1f%% within 0.5 K of T_CMB",
             medT, 100 * near)

    prof = []
    for lo, hi in [(0, 5), (5, 10), (10, 20), (20, 40), (40, 90)]:
        m = (ab >= lo) & (ab < hi)
        n = int(m.sum())
        prof.append({"b_lo": lo, "b_hi": hi, "n": n,
                     "cold_rate": float(cold[m].mean()) if n else None,
                     "mirror_rate": float(mirror[m].mean()) if n else None,
                     "sub_minus_half_rate": float(prev_railed[m].mean()) if n else None})
        log.info("  |b| %2d-%2d n=%6d  cold %.4f  mirror %.4f  beta<-0.5 %.4f",
                 lo, hi, n, prof[-1]["cold_rate"] or 0,
                 prof[-1]["mirror_rate"] or 0, prof[-1]["sub_minus_half_rate"] or 0)

    n_c, n_m = int(cold.sum()), int(mirror.sum())
    if n_c <= n_m:
        verdict = (
            f"NULL, and the earlier CMB diagnosis is withdrawn. On a grid "
            f"widened at BOTH edges (beta from {BETA_LO} to {BETA_HI}, "
            f"symmetric about dust) and restricted to fits lying strictly "
            f"inside it, the cold-blackbody tail holds {n_c} sources against "
            f"{n_m} in the symmetric unphysical mirror, a ratio of "
            f"{ratio:.2f}. The 4.4x asymmetry reported by script 78 was 88.6% "
            f"fits pinned at beta = -0.50, its grid floor: the model could not "
            f"reach the slope the data wanted, and the temperature clustering "
            f"near 2.7 K that the CMB argument rested on was a nuisance "
            f"parameter absorbing that. Widening the floor lets those sources "
            f"find their real slope, and {int(prev_railed.sum())} of them land "
            f"below -0.5 where flat-spectrum radio sources live "
            f"(beta = -2 is S ~ nu^0). Both the cold tail and its mirror are "
            f"now consistent, so this channel is a null.")
    else:
        verdict = (
            f"{n_c} interior cold-blackbody fits against {n_m} in the "
            f"symmetric mirror, ratio {ratio:.2f}. This survives the "
            f"boundary-occupancy cut that killed the script 78 result and "
            f"needs per-source validation against HI4PI, the Solar System "
            f"object registry and the Planck artefact flags before it means "
            f"anything.")

    print(f"\n{'=' * 74}")
    print("SEARCH T, THIRD PASS: UNRAILED AT BOTH EDGES")
    print(f"{'=' * 74}")
    print(f"  beta grid                           : "
          f"{BETA_LO} to {BETA_HI} ({len(s78.BETAS)} points)")
    print(f"  sources with >= 3 bands             : {len(s):,}")
    print(f"  interior (non-railed) fits          : "
          f"{int(interior.sum()):,} ({100*interior.mean():.1f}%)")
    print(f"  now fitting below the old floor     : {int(prev_railed.sum()):,}")
    print(f"  cold (interior)                     : {n_c}")
    print(f"  symmetric mirror (interior)         : {n_m}")
    print(f"  ratio                               : {ratio:.2f}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchT_unrailed_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "beta_grid": [BETA_LO, BETA_HI, len(s78.BETAS)],
        "n_sources": int(len(s)),
        "boundary_occupancy": edges,
        "n_interior": int(interior.sum()),
        "frac_interior": float(interior.mean()),
        "n_below_old_floor": int(prev_railed.sum()),
        "median_beta_below_old_floor": (
            float(np.median(beta[prev_railed & interior]))
            if (prev_railed & interior).any() else None),
        "n_cold_interior": n_c,
        "n_mirror_interior": n_m,
        "cold_over_mirror": ratio,
        "cold_median_T": medT,
        "cold_frac_within_0p5K_of_TCMB": near,
        "latitude_profile": prof,
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    s.to_csv(cfg.RESULT_DIR / f"searchT_sed_fits_unrailed_{args.tag}.csv",
             index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
