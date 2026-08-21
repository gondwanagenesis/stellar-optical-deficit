#!/usr/bin/env python
"""Search V diagnosis: is the 2.7 sigma excess sub-cell survey geometry?

    run.sh scripts/87_searchV_resolution_scan.py --tag primary

WHAT SEARCH V FOUND, AND WHY IT IS NOT YET A RESULT
----------------------------------------------------
Channel 22 measured the alignment of `ipd_gof_harmonic_phase` with the
proper-motion position angle over 224,081 high-RUWE sources. Observed
R = 0.04091. A global mark permutation puts that at the p floor, but so does
the local one: within nside=8 HEALPix cells the null mean is 0.03747 +- 0.00126,
so 92% of the apparent alignment is survey geometry and what remains is an
excess of +0.00344, or 2.7 sigma, at a preferred dphi of 3.9 deg.

Two features of that result argue against taking it at face value.

1. THE EXCESS IS FLAT. Every subset the channel cut -- four quartiles of
   harmonic amplitude, four of proper motion, four of multi-peak fraction, the
   dim tail, the bright tail -- returns an excess between +0.0005 and +0.0086,
   with no subset near zero. A physical per-star alignment should concentrate
   where the mechanism says it should: in high-amplitude, high-proper-motion,
   multi-peak sources, and vanish in the clean ones. A flat pedestal across
   every cut is the signature of something the estimator does to all of them
   equally.

2. TWO OF THE THREE POSITIVE CONTROLS FAIL. The proper-motion gradient, which
   the stationary-blend mechanism requires most directly -- the drift term is
   what swings the separation vector onto the PM axis -- comes out at 1.6
   sigma. The multi-peak gradient, on sources where Gaia itself saw a second
   peak, comes out at 0.55 sigma.

THE HYPOTHESIS THIS SCRIPT TESTS
--------------------------------
An nside=8 cell is 7.3 deg across and PA_PM varies by 49.5 deg inside one.
Both PA_PM and the scanning-law phase field are smooth on the sky, so they are
still correlated WITHIN a cell. A within-cell shuffle destroys that residual
sub-cell correlation while the observed R keeps it, and the difference appears
as a positive excess in every subset alike. If that is what the 2.7 sigma is,
the excess must fall toward zero as the cells shrink and the surviving
correlation is squeezed out.

WHY A FALLING EXCESS IS NOT AUTOMATICALLY THE ANSWER
-----------------------------------------------------
It falls for a trivial reason too. As nside rises the cells empty out, groups
of one cannot be shuffled, and a null that shuffles nothing reproduces the
observed R exactly -- driving the excess to zero with no geometry involved at
all. So the scan is only readable next to a measurement of the null's
remaining POWER, and this script measures it three ways at every nside:

* the fraction of sources sitting in a cell with at least one companion, which
  is the fraction the shuffle can touch at all;
* the within-cell circular dispersion of PA_PM, which is the variation the
  shuffle has to work with;
* and directly, by injecting an aligned population into the real phases and
  asking whether that nside's own null still recovers it.

An excess that falls while injected signal is still recovered is geometry.
An excess that falls only once the injection stops being recovered is a
statement about the test, not about the sky.

CORRECTNESS OF THE FAST SHUFFLE
-------------------------------
The scan needs ~10^4 within-cell permutations, which the group-loop version in
script 85 is too slow for at fine nside, so the shuffle here is vectorised with
a lexsort. That is a reimplementation of the estimator the verdict rests on,
which is exactly the kind of thing that has produced two silent errors in this
project already, so it is checked against the original at nside=8 before any
of the scan is believed.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
_v = __import__("85_searchV_harmonic_phase")

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchV_scan")

NSIDES = (2, 4, 8, 16, 32, 64, 128)
N_PERM = 500
RNG_SEED = 20260820


def resultant(phase, pa):
    th = np.deg2rad(2.0 * ((phase - pa) % 180.0))
    return float(np.hypot(np.mean(np.cos(th)), np.mean(np.sin(th))))


def permute_local_fast(phase, pa, cell, rng, n_perm):
    """Within-cell mark permutation, vectorised.

    `np.lexsort((key, cell))` orders by cell and randomly inside each cell, so
    writing those values into the cell-ordered slots is a uniform independent
    permutation within every group at once. Cells of one map to themselves,
    which is what the group-loop version does by skipping them.
    """
    dest = np.argsort(cell, kind="stable")
    out = np.empty(n_perm)
    pa_s = np.empty_like(pa)
    for i in range(n_perm):
        src = np.lexsort((rng.random(pa.size), cell))
        pa_s[dest] = pa[src]
        out[i] = resultant(phase, pa_s)
    return out


def shufflable_fraction(cell) -> float:
    _, n = np.unique(cell, return_counts=True)
    return float(n[n > 1].sum() / n.sum())


def within_cell_dispersion(pa_deg, cell, min_n=20):
    """Mean circular sd of PA_PM inside a cell, degrees. Power of the shuffle."""
    th = np.deg2rad(2.0 * np.asarray(pa_deg, float))
    df = pd.DataFrame({"c": cell, "x": np.cos(th), "y": np.sin(th)})
    g = df.groupby("c")[["x", "y"]].agg(["mean", "count"])
    r = np.hypot(g[("x", "mean")], g[("y", "mean")]).to_numpy()
    n = g[("x", "count")].to_numpy()
    r = r[n >= min_n]
    if r.size < 5:
        return None
    r = np.clip(r, 1e-6, 1 - 1e-9)
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(r))).mean() / 2.0)


def _trend(rows):
    if len(rows) < 2:
        return {"status": "too few bins"}
    lo, hi = rows[0], rows[-1]
    diff = hi["excess_R"] - lo["excess_R"]
    sd = float(np.hypot(hi["null_sd"], lo["null_sd"]))
    e = [r["excess_R"] for r in rows]
    rho = float(np.corrcoef(np.argsort(np.argsort(e)),
                            np.arange(len(e)))[0, 1]) if len(e) > 2 else None
    return {"top_minus_bottom": float(diff), "sigma": sd,
            "nsigma": float(diff / sd) if sd > 0 else None,
            "spearman_rho_over_bins": rho,
            "passes": bool(sd > 0 and diff / sd > 2.0)}


def _verdict(rows, ruwe_trend, r_obs):
    """Read the scan only where the null still demonstrably has power."""
    live = [r for r in rows
            if any(p["recovered"] for p in r["injection_power"])]
    if not live:
        return ("INCONCLUSIVE: the local null recovered no injected population "
                "at any nside, so the scan measures the test's power and not "
                "the sky.")
    base = next((r for r in rows if r["nside"] == 8), rows[0])
    finest = live[-1]
    drop = base["excess_R"] - finest["excess_R"]
    frac = drop / base["excess_R"] if base["excess_R"] else float("nan")
    tail = ""
    if ruwe_trend.get("passes"):
        tail = (" The excess also rises with RUWE at "
                f"{ruwe_trend['nsigma']:.1f} sigma, which points at Gaia's own "
                "window placement -- a mis-centred source drifts off centre "
                "ALONG its proper motion and degrades the PSF fit "
                "anisotropically -- rather than at anything on the sky.")
    if finest["nsigma"] is not None and finest["nsigma"] < 2.0:
        return (f"SUB-CELL SURVEY GEOMETRY. The nside=8 excess of "
                f"{base['excess_R']:+.5f} ({base['nsigma']:.1f} sigma) falls to "
                f"{finest['excess_R']:+.5f} ({finest['nsigma']:.1f} sigma) at "
                f"nside={finest['nside']} ({finest['cell_deg']:.2f} deg cells), "
                f"a loss of {100 * frac:.0f}%, while that same null still "
                f"recovers an injected aligned population. The null retains "
                f"power and the signal does not survive, so the excess "
                f"reported by channel 22 is residual correlation between two "
                f"smooth sky fields inside a 7.3 deg cell. Channel 22 is NULL "
                f"and its verdict string must be corrected.{tail}")
    return (f"THE EXCESS SURVIVES. At nside={finest['nside']} "
            f"({finest['cell_deg']:.2f} deg cells) the excess is still "
            f"{finest['excess_R']:+.5f} ({finest['nsigma']:.1f} sigma, "
            f"p={finest['p_local']:.3g}) with injected populations still "
            f"recovered, so it is not a smooth-field correlation on any scale "
            f"the scan reaches. It is a per-star pairing of sub-PSF "
            f"orientation with the proper-motion axis.{tail}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--inject", type=float, nargs="*", default=[0.02, 0.005])
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    import healpy as hp

    d = pd.read_parquet(cfg.RAW_DIR / "high_ruwe_500pc.parquet")
    ok = (d["ipd_gof_harmonic_phase"].notna()
          & d["ipd_gof_harmonic_amplitude"].notna()
          & d["pmra"].notna() & d["pmdec"].notna())
    d = d[ok].reset_index(drop=True)
    pmra = d["pmra"].to_numpy(float)
    pmdec = d["pmdec"].to_numpy(float)
    pa = np.rad2deg(np.arctan2(pmra, pmdec)) % 180.0
    phase = d["ipd_gof_harmonic_phase"].to_numpy(float) % 180.0
    theta = np.radians(90.0 - d["dec"].to_numpy(float))
    phi = np.radians(d["ra"].to_numpy(float))

    r_obs = resultant(phase, pa)
    out = {"tag": args.tag, "n": int(len(d)), "R_observed": r_obs,
           "n_perm": int(args.n_perm), "seed": RNG_SEED}
    log.info("n=%d  observed R=%.5f", len(d), r_obs)

    # ---- the sample is 200 source_id stripes, and source_id is sky ---------
    # DR3 source_id embeds the level-12 HEALPix index in its high bits, so the
    # every-5th-chunk grid scripts 80/81/84 share is a set of 200 contiguous
    # SKY stripes, not a random 20% of the population. It is why only 356 of
    # 768 nside=8 cells hold any source at all. This does not bias a
    # within-cell shuffle -- a partly covered cell is a smaller region, over
    # which the smooth fields vary less, which makes the null more
    # conservative rather than less -- but it is a property of the sample that
    # anyone reading a spatial statistic off it has to know.
    c8 = hp.ang2pix(8, theta, phi)
    out["sample_geometry"] = {
        "nside8_cells_occupied": int(np.unique(c8).size),
        "nside8_cells_total": int(hp.nside2npix(8)),
        "note": ("source_id chunks are sky stripes: DR3 source_id carries the "
                 "level-12 HEALPix index in its high bits, so taking every "
                 "5th of 1000 id chunks samples 200 contiguous sky stripes. "
                 "Partial cell coverage makes the within-cell null more "
                 "conservative, not less.")}
    log.info("sample occupies %d/%d nside=8 cells (source_id chunks are sky "
             "stripes)", np.unique(c8).size, hp.nside2npix(8))

    # ---- the fast shuffle must reproduce the estimator the verdict used ----
    log.info("cross-checking the vectorised shuffle against script 85's "
             "group loop at nside=8 ...")
    rng_a = np.random.default_rng(1234)
    rng_b = np.random.default_rng(1234)
    n_chk = 60
    a = _v.permute_local(phase, pa, c8, rng_a, n_chk)
    b = permute_local_fast(phase, pa, c8, rng_b, n_chk)
    # Different RNG consumption patterns, so the draws differ; the null
    # DISTRIBUTION is what has to agree.
    dm = abs(float(a.mean() - b.mean()))
    tol = 3.0 * float(np.hypot(a.std(), b.std())) / np.sqrt(n_chk)
    out["estimator_crosscheck"] = {
        "n_draws": n_chk,
        "group_loop_mean_R": float(a.mean()), "group_loop_sd": float(a.std()),
        "vectorised_mean_R": float(b.mean()), "vectorised_sd": float(b.std()),
        "mean_difference": dm, "tolerance_3se": tol,
        "agrees": bool(dm < tol)}
    log.info("  group loop %.5f+-%.5f vs vectorised %.5f+-%.5f  (diff %.5f, "
             "3se tol %.5f) -> %s", a.mean(), a.std(), b.mean(), b.std(),
             dm, tol, "OK" if dm < tol else "MISMATCH")
    if dm >= tol:
        log.error("the vectorised shuffle does not reproduce the original "
                  "null; refusing to run a scan on it")
        out["verdict"] = ("ABORTED: the fast within-cell shuffle does not "
                          "reproduce script 85's null distribution.")
        (cfg.RESULT_DIR / f"searchV_resolution_scan_{args.tag}.json"
         ).write_text(json.dumps(out, indent=2))
        return 1

    # ---- the scan ---------------------------------------------------------
    rows = []
    for ns in NSIDES:
        cell = hp.ang2pix(ns, theta, phi)
        nb = permute_local_fast(phase, pa, cell, rng, args.n_perm)
        exc = r_obs - nb.mean()
        row = {"nside": int(ns),
               "cell_deg": float(np.degrees(np.sqrt(hp.nside2pixarea(ns)))),
               "cells_occupied": int(np.unique(cell).size),
               "shufflable_fraction": shufflable_fraction(cell),
               "within_cell_pa_dispersion_deg": within_cell_dispersion(pa, cell),
               "null_R": float(nb.mean()), "null_sd": float(nb.std()),
               "excess_R": float(exc),
               "nsigma": float(exc / nb.std()) if nb.std() > 0 else None,
               "p_local": _v.perm_p(r_obs, nb)}

        # power at this nside, measured on the real sky
        power = []
        for f in args.inject:
            ph = phase.copy()
            m = rng.random(ph.size) < f
            ph[m] = (pa[m] + rng.normal(0, 15, int(m.sum()))) % 180.0
            ri = resultant(ph, pa)
            nbi = permute_local_fast(ph, pa, cell, rng, args.n_perm)
            p = _v.perm_p(ri, nbi)
            power.append({"injected_fraction": float(f), "R": ri,
                          "null_R": float(nbi.mean()),
                          "excess_R": float(ri - nbi.mean()),
                          "p_local": p, "recovered": bool(p < 0.01)})
        row["injection_power"] = power
        rows.append(row)
        log.info("nside %3d (%.2f deg): shufflable %.3f  disp %s  null "
                 "%.5f+-%.5f  excess %+.5f (%.2f sig, p=%.4g)  injected "
                 "recovery %s", ns, row["cell_deg"], row["shufflable_fraction"],
                 ("%.1f" % row["within_cell_pa_dispersion_deg"]
                  if row["within_cell_pa_dispersion_deg"] else "n/a"),
                 nb.mean(), nb.std(), exc, row["nsigma"] or float("nan"),
                 row["p_local"],
                 {p["injected_fraction"]: p["recovered"] for p in power})

    out["scan"] = rows

    # ---- RUWE gradient: the astrometric-mismodelling alternative -----------
    # There is a per-star artefact that is neither a blend nor a
    # technosignature. Gaia predicts where each source will be and reads out a
    # small window there; a source whose astrometric model is wrong sits off
    # the window centre, and the offset accumulates ALONG the proper-motion
    # direction. A mis-centred source degrades the PSF fit anisotropically, so
    # its harmonic phase would prefer the PM axis for reasons internal to the
    # instrument. That mechanism scales with how badly the model fits, which
    # is what RUWE measures, so an excess rising steeply with RUWE points here
    # rather than at anything on the sky. The whole sample is RUWE >= 1.4 by
    # construction, so this is a gradient within the disturbed population.
    ruwe = d["ruwe"].to_numpy(float)
    q = np.nanpercentile(ruwe, [25, 50, 75])
    bins, rrows = np.digitize(ruwe, q), []
    for b in range(4):
        m = bins == b
        if m.sum() < 200:
            continue
        ri = resultant(phase[m], pa[m])
        nbi = permute_local_fast(phase[m], pa[m], c8[m], rng,
                                 max(200, args.n_perm // 2))
        rrows.append({"bin": int(b), "n": int(m.sum()),
                      "lo": float(ruwe[m].min()), "hi": float(ruwe[m].max()),
                      "R": ri, "null_R": float(nbi.mean()),
                      "null_sd": float(nbi.std()),
                      "excess_R": float(ri - nbi.mean()),
                      "p_local": _v.perm_p(ri, nbi)})
    out["control_ruwe"] = {"bins": rrows, "trend": _trend(rrows)}
    log.info("RUWE gradient: excess by bin = %s -> %s",
             ", ".join(f"{r['excess_R']:+.5f}" for r in rrows),
             json.dumps(out["control_ruwe"]["trend"]))

    out["verdict"] = _verdict(rows, out["control_ruwe"]["trend"], r_obs)

    print("\n" + "=" * 78)
    print("SEARCH V DIAGNOSIS: does the excess survive shrinking the cell?")
    print("=" * 78)
    print(f"  observed R  {r_obs:.5f}   (fixed; only the null changes below)\n")
    print(f"  {'nside':>6} {'deg':>6} {'shuf':>6} {'disp':>6} {'null R':>9} "
          f"{'excess':>9} {'sig':>6} {'p':>8}  injected recovery")
    for r in rows:
        rec = " ".join(f"{p['injected_fraction']:g}:"
                       f"{'yes' if p['recovered'] else 'NO'}"
                       for p in r["injection_power"])
        disp = (f"{r['within_cell_pa_dispersion_deg']:.1f}"
                if r["within_cell_pa_dispersion_deg"] else "  n/a")
        print(f"  {r['nside']:>6} {r['cell_deg']:>6.2f} "
              f"{r['shufflable_fraction']:>6.3f} {disp:>6} {r['null_R']:>9.5f} "
              f"{r['excess_R']:>+9.5f} {r['nsigma']:>6.2f} {r['p_local']:>8.4g}"
              f"  {rec}")
    print(f"\nVERDICT: {out['verdict']}\n")

    p = cfg.RESULT_DIR / f"searchV_resolution_scan_{args.tag}.json"
    p.write_text(json.dumps(out, indent=2))
    log.info("wrote %s", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
