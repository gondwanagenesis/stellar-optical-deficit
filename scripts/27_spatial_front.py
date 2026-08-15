#!/usr/bin/env python
"""Is there a SPREADING FRONT? Spatial structure in the residual field.

    run.sh scripts/27_spatial_front.py --tag primary

Every other test in this repository treats stars one at a time, and therefore
pays the full single-star background for every candidate.  But the hypothesis
is *spreading*, and spreading is a spatial process: a front is a connected
region, not a scatter of unrelated points.  Correlation statistics detect
clustered weak signals far below the per-object threshold, and unlike the tail
statistic their sensitivity keeps improving with N, because the astrophysical
background (blends, activity, binaries) is spatially RANDOM while a front is
not.

It also makes a prediction nothing else here makes: an EDGE.  A discontinuity
in mean residual across a surface in 3D, with ordinary stars on the far side.

METHOD, AND THE ONE THING THAT WOULD RUIN IT
--------------------------------------------
Interstellar dust is spatially correlated and is our dominant systematic, so a
raw correlation function of residuals measures dust, not aliens.  The residuals
are therefore first cleaned by regressing out every systematic we can measure
-- extinction, apparent magnitude, crowding, colour, metallicity, blending and
astrometric-quality proxies -- using a model that is NEVER shown (l, b,
distance).  Because the model cannot see position, it cannot absorb a spatial
front; it can only absorb structure that is explained by the systematics.

A_0 is itself positional, so removing its effect does remove genuinely
spatial structure -- but only the part that follows the dust, which a front has
no reason to do.

The null is measured, not assumed: residuals are shuffled among stars, which
destroys any position association while preserving the residual distribution
and the wildly uneven star counts per cell.  Everything is quoted against that
shuffled null.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

from pipeline import config as cfg
from pipeline import statistics as st

CELL_PC = 25.0            # 3D cell size; ~10-100 stars per cell at 500 pc
MIN_STARS_PER_CELL = 20
N_SHUFFLE = 20
# Pair counting grows as r^3: at 90 pc with 400k stars it is 4.6e8 pairs and
# the process is killed. Small scales are done by direct pair counting on a
# subsample; large scales come from the cell means instead, which is O(N).
CORR_SUBSAMPLE = 200_000
CORR_BINS_PC = np.array([1, 2, 4, 8, 15, 25])

SYSTEMATIC_FEATURES = ["A_0", "phot_g_mean_mag", "sky_density", "M_Ks",
                       "bp_rp0", "mh_gspphot", "cstar_nsigma", "ruwe",
                       "astrometric_excess_noise", "parallax_over_error"]


def clean_residuals(d: pd.DataFrame, r: np.ndarray) -> np.ndarray:
    """Remove everything explainable by NON-POSITIONAL systematics."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    X = d[SYSTEMATIC_FEATURES].to_numpy(dtype=np.float32)
    m = HistGradientBoostingRegressor(
        max_iter=200, learning_rate=0.1, max_depth=6,
        early_stopping=True, validation_fraction=0.1, random_state=0)
    m.fit(X, r)
    pred = m.predict(X)
    print(f"  systematics model explains "
          f"{100*(1 - np.var(r - pred)/np.var(r)):.2f}% of the residual variance")
    return r - pred


def cell_scan(xyz: np.ndarray, rc: np.ndarray, rng: np.random.Generator,
              n_shuffle: int = N_SHUFFLE) -> dict:
    """Mean cleaned residual per 3D cell, against a shuffled null."""
    ijk = np.floor((xyz + 1300.0) / CELL_PC).astype(np.int32)
    key = (ijk[:, 0].astype(np.int64) * 4_000_000
           + ijk[:, 1].astype(np.int64) * 2000 + ijk[:, 2])
    order = np.argsort(key, kind="stable")
    k_sorted, r_sorted = key[order], rc[order]
    edges = np.flatnonzero(np.diff(k_sorted)) + 1
    starts = np.concatenate([[0], edges])
    ends = np.concatenate([edges, [len(k_sorted)]])
    counts = ends - starts
    keep = counts >= MIN_STARS_PER_CELL
    sums = np.add.reduceat(r_sorted, starts)[keep]
    n = counts[keep]
    means = sums / n
    sigma = float(np.std(rc, ddof=1))
    z = means / (sigma / np.sqrt(n))
    cell_keys = k_sorted[starts][keep]

    # shuffled null: same cells, same counts, residuals reassigned at random
    z_null = []
    for _ in range(n_shuffle):
        rs = rng.permutation(rc)[order]
        s = np.add.reduceat(rs, starts)[keep]
        z_null.append(s / n / (sigma / np.sqrt(n)))
    z_null = np.concatenate(z_null)

    out = {"n_cells": int(keep.sum()), "median_stars_per_cell": int(np.median(n)),
           "z_std_real": float(np.std(z)), "z_std_null": float(np.std(z_null))}
    for thr in (3.0, 4.0, 5.0):
        obs = int(np.count_nonzero(z > thr))
        exp = float(np.count_nonzero(z_null > thr)) / n_shuffle
        out[f"cells_above_{thr:.0f}sigma"] = obs
        out[f"expected_above_{thr:.0f}sigma"] = exp
        out[f"excess_above_{thr:.0f}sigma"] = obs - exp
    return out, z, cell_keys, ijk, keep, starts, k_sorted


def cell_aux(xyz: np.ndarray, d: pd.DataFrame, cell_keys: np.ndarray):
    """Per-cell median |b| and A_0, aligned with the cell_keys ordering."""
    ijk = np.floor((xyz + 1300.0) / CELL_PC).astype(np.int32)
    key = (ijk[:, 0].astype(np.int64) * 4_000_000
           + ijk[:, 1].astype(np.int64) * 2000 + ijk[:, 2])
    df = pd.DataFrame({"key": key, "b": d["b"].to_numpy(float),
                       "a0": d["A_0"].to_numpy(float)})
    g = df.groupby("key").agg(b=("b", "median"), a0=("a0", "median"),
                              n=("b", "size"))
    g = g.reindex(cell_keys)
    return (g["b"].to_numpy(float), g["a0"].to_numpy(float),
            g["n"].to_numpy(float))


def connected_regions(cell_keys: np.ndarray, z: np.ndarray,
                      thr: float = 3.0) -> dict:
    """Largest CONNECTED clump of anomalous cells -- the front signature."""
    sel = z > thr
    if sel.sum() < 2:
        return {"n_anomalous_cells": int(sel.sum()), "largest_clump": int(sel.sum())}
    k = cell_keys[sel]
    i = (k // 4_000_000).astype(np.int64)
    j = ((k % 4_000_000) // 2000).astype(np.int64)
    l = (k % 2000).astype(np.int64)
    pts = np.column_stack([i, j, l]).astype(float)
    tree = cKDTree(pts)
    # 26-connectivity: neighbours within sqrt(3) in cell units
    pairs = tree.query_pairs(1.75, output_type="ndarray")
    parent = np.arange(len(pts))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in pairs:
        ra_, rb = find(a), find(b)
        if ra_ != rb:
            parent[rb] = ra_
    roots = np.array([find(a) for a in range(len(pts))])
    sizes = np.bincount(roots)
    return {"n_anomalous_cells": int(sel.sum()),
            "largest_clump": int(sizes.max()),
            "n_clumps": int((sizes > 0).sum())}


def correlation_function(xyz: np.ndarray, rc: np.ndarray,
                         rng: np.random.Generator) -> pd.DataFrame:
    """xi(s) = <dr_i dr_j> in 3D separation bins, against a shuffled null."""
    n = len(rc)
    if n > CORR_SUBSAMPLE:
        idx = rng.choice(n, CORR_SUBSAMPLE, replace=False)
        xyz, rc = xyz[idx], rc[idx]
    tree = cKDTree(xyz)
    rows = []
    prev_sum = prev_cnt = 0.0
    shuf = rng.permutation(rc)
    prev_sum_n = prev_cnt_n = 0.0
    for s in CORR_BINS_PC:
        pairs = tree.query_pairs(float(s), output_type="ndarray")
        if len(pairs) == 0:
            continue
        prod = rc[pairs[:, 0]] * rc[pairs[:, 1]]
        tot_sum, tot_cnt = prod.sum(), float(len(prod))
        prod_n = shuf[pairs[:, 0]] * shuf[pairs[:, 1]]
        tot_sum_n, tot_cnt_n = prod_n.sum(), float(len(prod_n))

        d_sum, d_cnt = tot_sum - prev_sum, tot_cnt - prev_cnt
        d_sum_n, d_cnt_n = tot_sum_n - prev_sum_n, tot_cnt_n - prev_cnt_n
        if d_cnt > 100:
            xi = d_sum / d_cnt
            xi_n = d_sum_n / d_cnt_n
            err = np.std(prod) / np.sqrt(d_cnt)
            rows.append({"s_max_pc": float(s), "n_pairs": int(d_cnt),
                         "xi": xi, "xi_shuffled": xi_n, "err": err,
                         "n_sigma": (xi - xi_n) / err if err > 0 else np.nan})
        prev_sum, prev_cnt = tot_sum, tot_cnt
        prev_sum_n, prev_cnt_n = tot_sum_n, tot_cnt_n
        if len(pairs) > 3e8:
            break
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--min-abs-b", type=float, default=0.0,
                    help="restrict to |b| above this. The decisive test: the "
                         "full-sky structure is concentrated in the plane and "
                         "at high A_0, i.e. it is dust. At high |b| there is "
                         "almost no dust, so any structure surviving there is "
                         "NOT explicable that way.")
    ap.add_argument("--max-a0", type=float, default=None)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()
    rng = np.random.default_rng(20260815)
    label = args.label or args.tag

    cols = ["l", "b", "dist_pc", "residual", "pmra", "pmdec"] + SYSTEMATIC_FEATURES
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=sorted(set(cols)))
    d = d.dropna(subset=["l", "b", "dist_pc", "residual"]).reset_index(drop=True)
    if args.min_abs_b > 0:
        d = d[d["b"].abs() >= args.min_abs_b].reset_index(drop=True)
        print(f"restricted to |b| >= {args.min_abs_b:.0f} deg")
    if args.max_a0 is not None:
        d = d[d["A_0"] <= args.max_a0].reset_index(drop=True)
        print(f"restricted to A_0 <= {args.max_a0}")
    for c in SYSTEMATIC_FEATURES:
        d[c] = d[c].fillna(d[c].median())
    r = d["residual"].to_numpy(float)
    print(f"{len(d):,} stars, raw sigma = {st.robust_sigma(r):.5f} mag")

    print("\n=== removing non-positional systematics ===")
    rc = clean_residuals(d, r)
    print(f"  cleaned sigma = {st.robust_sigma(rc):.5f} mag")

    lam, bet = np.radians(d["l"].to_numpy(float)), np.radians(d["b"].to_numpy(float))
    dist = d["dist_pc"].to_numpy(float)
    xyz = np.column_stack([dist * np.cos(bet) * np.cos(lam),
                           dist * np.cos(bet) * np.sin(lam),
                           dist * np.sin(bet)])

    print(f"\n=== 3D cell scan ({CELL_PC:.0f} pc cells) ===")
    stats, z, cell_keys, *_ = cell_scan(xyz, rc, rng)
    for k, v in stats.items():
        print(f"  {k:32s} {v}")

    # --- IS THE STRUCTURE JUST DUST? -------------------------------------
    # Decisive discriminant: interstellar dust is concentrated toward the
    # galactic plane and toward high A_0.  A spreading front has no reason to
    # be.  If the anomalous cells track |b| and A_0, this is unmodelled
    # extinction and not a technosignature.
    print(f"\n=== is the structure dust, or not? ===")
    cell_b, cell_a0, cell_n = cell_aux(xyz, d, cell_keys)
    hi = z > 3.0
    print(f"  anomalous cells (z>3)  : median |b| = "
          f"{np.median(np.abs(cell_b[hi])):5.1f} deg,  median A_0 = "
          f"{np.median(cell_a0[hi]):.4f}")
    print(f"  all cells              : median |b| = "
          f"{np.median(np.abs(cell_b)):5.1f} deg,  median A_0 = "
          f"{np.median(cell_a0):.4f}")
    rho_b = float(np.corrcoef(z, np.abs(cell_b))[0, 1])
    rho_a = float(np.corrcoef(z, cell_a0)[0, 1])
    print(f"  corr(z, |b|)           : {rho_b:+.3f}")
    print(f"  corr(z, A_0)           : {rho_a:+.3f}")
    frac_lowb = float(np.mean(np.abs(cell_b[hi]) < 20))
    frac_lowb_all = float(np.mean(np.abs(cell_b) < 20))
    print(f"  fraction at |b|<20 deg : {100*frac_lowb:.1f}% of anomalous vs "
          f"{100*frac_lowb_all:.1f}% of all cells")
    stats.update({"corr_z_absb": rho_b, "corr_z_A0": rho_a,
                  "frac_anomalous_lowb": frac_lowb,
                  "frac_all_lowb": frac_lowb_all,
                  "median_absb_anomalous": float(np.median(np.abs(cell_b[hi]))),
                  "median_absb_all": float(np.median(np.abs(cell_b)))})

    print(f"\n=== connected anomalous regions (the front signature) ===")
    for thr in (3.0, 4.0):
        cc = connected_regions(cell_keys, z, thr)
        # Null: the same NUMBER of anomalous cells, placed at random among the
        # occupied cells. Without this, a large clump means nothing -- cells
        # are not uniformly distributed, so clumping happens for free.
        null_sizes = []
        n_anom = cc["n_anomalous_cells"]
        for _ in range(20):
            fake_z = np.full(len(z), -9.0)
            fake_z[rng.choice(len(z), n_anom, replace=False)] = 9.0
            null_sizes.append(connected_regions(cell_keys, fake_z, 0.0)["largest_clump"])
        cc["largest_clump_null_mean"] = float(np.mean(null_sizes))
        cc["largest_clump_null_max"] = int(np.max(null_sizes))
        print(f"  at z > {thr:.0f}: {n_anom} cells in {cc['n_clumps']} clumps; "
              f"largest = {cc['largest_clump']} cells vs null "
              f"{cc['largest_clump_null_mean']:.1f} (max {cc['largest_clump_null_max']})")
        stats[f"connected_z{thr:.0f}"] = cc

    print(f"\n=== 3D two-point correlation of cleaned residuals ===")
    corr = correlation_function(xyz, rc, rng)
    print(corr.to_string(index=False, float_format=lambda v: f"{v:12.6g}"))
    corr.to_csv(cfg.RESULT_DIR / f"spatial_corr_{label}.csv", index=False)

    worst = corr.loc[corr["n_sigma"].abs().idxmax()] if len(corr) else None
    verdict = "NO spatial structure beyond the shuffled null"
    if worst is not None and abs(worst["n_sigma"]) > 5:
        # Detecting structure is NOT detecting a front. Dust is spatially
        # correlated too, and the extinction correction is known to be ~19%
        # short. Only call it interesting if the anomalous cells are NOT
        # preferentially dusty and NOT preferentially in the plane.
        dusty = (np.median(cell_a0[hi]) > 1.10 * np.median(cell_a0)
                 or np.median(np.abs(cell_b[hi])) < 0.90 * np.median(np.abs(cell_b)))
        if dusty:
            verdict = (f"structure at {worst['n_sigma']:+.1f} sigma, but the "
                       f"anomalous cells are preferentially dusty and/or "
                       f"low-latitude -- consistent with UNMODELLED EXTINCTION, "
                       f"not a front")
        else:
            verdict = (f"structure at s < {worst['s_max_pc']:.0f} pc "
                       f"({worst['n_sigma']:+.1f} sigma) and NOT dust-aligned "
                       f"-- INVESTIGATE")
    print(f"\nVERDICT: {verdict}")

    stats["verdict"] = verdict
    stats["cell_pc"] = CELL_PC
    stats["n_stars"] = int(len(d))
    stats["sigma_clean"] = float(st.robust_sigma(rc))
    (cfg.RESULT_DIR / f"spatial_front_{label}.json").write_text(
        json.dumps(stats, indent=2, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
