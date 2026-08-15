#!/usr/bin/env python
"""I1 v2: kinematic clustering with the controls the v1 result demanded.

    run.sh scripts/36_velocity_clustering_v2.py --tag primary

WHAT v1 FOUND AND WHY IT WAS NOT A DETECTION
---------------------------------------------
v1 reported a 4301-sigma excess of low-dv anomaly pairs. Three controls were
missing, and each removed a large part of it:

  1. bound systems -- dv = 0 by construction. Removed (tiny effect, 66 pairs).
  2. |v| matching  -- anomalies are KINEMATICALLY COLD (median 22.8 vs 35.1
     km/s; 43% below 20 km/s against 23% for the field). Low |v| makes low
     pairwise dv automatic. Matching |v| cut the excess from 160x to 43x.
  3. SEPARATION matching -- this script. v is tangential only, so two stars
     far apart on the sky have their velocity vectors projected differently
     and get a spuriously large dv. Anomalies are more spatially concentrated
     than any control set (3.9M pairs vs 2.0M), so their pairs sit at smaller
     separations and therefore at systematically smaller dv. Comparing whole
     samples confounds this with the signal.

The fix is to compare dv distributions ONLY WITHIN MATCHED SEPARATION BINS.

Plus a youth veto: kinematically cold stars are young, young stars sit in
associations, and associations cluster in velocity for entirely mundane
reasons. Restricting to |v| > 40 km/s selects the dynamically heated thin/thick
disc, where coherent young associations are absent. If the excess survives
there, youth is not the explanation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from pipeline import config as cfg
from pipeline import statistics as st

K_SIGMA = 5.0
R_PC = 100.0
MIN_SEP_PC = 1.0
SEP_EDGES = np.array([1, 5, 10, 20, 35, 55, 75, 100])
DV_THRESH = 5.0          # km/s


def velocity_vectors(d):
    ra = np.radians(d["ra"].to_numpy(float))
    dec = np.radians(d["dec"].to_numpy(float))
    k = 4.74047 / np.maximum(d["parallax_corr"].to_numpy(float), 1e-6)
    va = k * d["pmra"].to_numpy(float)
    vd = k * d["pmdec"].to_numpy(float)
    e_a = np.column_stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)])
    e_d = np.column_stack([-np.sin(dec) * np.cos(ra),
                           -np.sin(dec) * np.sin(ra), np.cos(dec)])
    return va[:, None] * e_a + vd[:, None] * e_d


def xyz_of(d):
    l = np.radians(d["l"].to_numpy(float)); b = np.radians(d["b"].to_numpy(float))
    r = d["dist_pc"].to_numpy(float)
    return np.column_stack([r * np.cos(b) * np.cos(l),
                            r * np.cos(b) * np.sin(l), r * np.sin(b)])


def frac_by_sep(xyz, v, rng=None, max_pairs=4_000_000):
    """Fraction of pairs with dv < DV_THRESH, in bins of spatial separation."""
    pairs = cKDTree(xyz).query_pairs(R_PC, output_type="ndarray")
    if len(pairs) == 0:
        return np.full(len(SEP_EDGES) - 1, np.nan), np.zeros(len(SEP_EDGES) - 1)
    if len(pairs) > max_pairs and rng is not None:
        pairs = pairs[rng.choice(len(pairs), max_pairs, replace=False)]
    sep = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
    keep = sep > MIN_SEP_PC
    pairs, sep = pairs[keep], sep[keep]
    dv = np.linalg.norm(v[pairs[:, 0]] - v[pairs[:, 1]], axis=1)
    idx = np.digitize(sep, SEP_EDGES) - 1
    fr = np.full(len(SEP_EDGES) - 1, np.nan)
    n = np.zeros(len(SEP_EDGES) - 1)
    for b in range(len(SEP_EDGES) - 1):
        m = idx == b
        n[b] = m.sum()
        if n[b] > 50:
            fr[b] = np.mean(dv[m] < DV_THRESH)
    return fr, n


def run(d, r, label, args, rng):
    sig = st.robust_sigma(r); med = float(np.median(r))
    anom = r > med + K_SIGMA * sig
    if anom.sum() < 200:
        print(f"  {label}: only {anom.sum()} anomalies, skipping")
        return None
    xyz, v = xyz_of(d), velocity_vectors(d)
    speed = np.linalg.norm(v, axis=1)
    dist = d["dist_pc"].to_numpy(float)

    import healpy as hp
    pix = hp.ang2pix(32, np.radians(90 - d["b"].to_numpy(float)),
                     np.radians(d["l"].to_numpy(float)), nest=True)
    dens = np.bincount(pix, minlength=hp.nside2npix(32))[pix].astype(float)

    f_obs, n_obs = frac_by_sep(xyz[anom], v[anom], rng)

    d_e = np.percentile(dist[anom], np.linspace(0, 100, 7))
    r_e = np.percentile(dens[anom], np.linspace(0, 100, 4))
    v_e = np.percentile(speed[anom], np.linspace(0, 100, 7))
    idx_all = np.arange(len(d))
    ctrl = []
    for _ in range(args.n_control):
        pick = []
        for a in range(len(d_e) - 1):
            for c_ in range(len(r_e) - 1):
                for w in range(len(v_e) - 1):
                    want = int(((dist[anom] >= d_e[a]) & (dist[anom] < d_e[a+1])
                                & (dens[anom] >= r_e[c_]) & (dens[anom] < r_e[c_+1])
                                & (speed[anom] >= v_e[w])
                                & (speed[anom] < v_e[w+1])).sum())
                    if not want:
                        continue
                    pool = idx_all[(dist >= d_e[a]) & (dist < d_e[a+1])
                                   & (dens >= r_e[c_]) & (dens < r_e[c_+1])
                                   & (speed >= v_e[w]) & (speed < v_e[w+1]) & ~anom]
                    if len(pool):
                        pick.append(rng.choice(pool, min(want, len(pool)), False))
        sel = np.concatenate(pick)
        fr, _ = frac_by_sep(xyz[sel], v[sel], rng)
        ctrl.append(fr)
    ctrl = np.array(ctrl)

    print(f"\n  {label}: {int(anom.sum()):,} anomalies")
    rows = []
    for b in range(len(SEP_EDGES) - 1):
        mu, sd = np.nanmean(ctrl[:, b]), np.nanstd(ctrl[:, b], ddof=1)
        z = (f_obs[b] - mu) / sd if sd > 0 and np.isfinite(f_obs[b]) else np.nan
        rows.append({"sep_pc": f"{SEP_EDGES[b]}-{SEP_EDGES[b+1]}",
                     "n_pairs": int(n_obs[b]), "frac_anom": f_obs[b],
                     "frac_ctrl": mu, "ratio": f_obs[b]/mu if mu > 0 else np.nan,
                     "n_sigma": z})
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))
    return t


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-control", type=int, default=12)
    args = ap.parse_args()
    rng = np.random.default_rng(20260816)

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["ra", "dec", "l", "b", "dist_pc", "residual",
                                 "parallax_corr", "pmra", "pmdec",
                                 "wise_w1mpro", "wise_w2mpro",
                                 "wise_w1mpro_error", "wise_w2mpro_error"])
    ok = (np.isfinite(d["pmra"]) & np.isfinite(d["pmdec"])
          & np.isfinite(d["dist_pc"])).to_numpy()
    d = d[ok].reset_index(drop=True)
    r = d["residual"].to_numpy(float)
    print(f"{len(d):,} stars")

    out = {}
    t_all = run(d, r, "ALL (separation-matched)", args, rng)
    if t_all is not None:
        out["all"] = t_all.to_dict(orient="records")

    v = velocity_vectors(d); speed = np.linalg.norm(v, axis=1)
    fast = speed > 40.0
    print(f"\n  youth veto: |v| > 40 km/s keeps {fast.sum():,} stars "
          f"({100*fast.mean():.0f}%) -- dynamically heated disc, "
          f"coherent young associations absent")
    t_fast = run(d[fast].reset_index(drop=True), r[fast],
                 "KINEMATICALLY OLD only", args, rng)
    if t_fast is not None:
        out["kinematically_old"] = t_fast.to_dict(orient="records")

    zz = []
    for k, v_ in out.items():
        zz += [row["n_sigma"] for row in v_ if np.isfinite(row.get("n_sigma", np.nan))]
    verdict = ("no kinematic clustering survives the controls"
               if not zz or max(np.abs(zz)) < 5 else
               f"residual excess up to {max(np.abs(zz)):.1f} sigma")
    print(f"\nVERDICT: {verdict}")
    out["verdict"] = verdict
    (cfg.RESULT_DIR / f"velocity_clustering_v2_{args.tag}.json").write_text(
        json.dumps(out, indent=2, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
