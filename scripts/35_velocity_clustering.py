#!/usr/bin/env python
"""I1: do anomalies cluster in VELOCITY space rather than position space?

    run.sh scripts/35_velocity_clustering.py --tag primary

THE IDEA
--------
scripts/27 searched for a spreading front in 3D position and found only dust.
But interstellar transfer cost is set by DELTA-V, not by distance. Two stars
50 pc apart co-moving at 2 km/s are adjacent in transfer cost; two stars 5 pc
apart differing by 60 km/s are not. A propagating process spreads along the
kinematic network, so position-space clustering is the wrong metric and its
null result is not informative about the hypothesis.

STATISTIC
---------
For stars within R pc of each other, compare the distribution of |dv| for
anomaly-anomaly pairs against a control set matched in distance and sky
position. An excess of low-|dv| anomaly pairs is the signature.

VELOCITY IS TANGENTIAL ONLY
---------------------------
Radial velocities were not pulled, so v is built from proper motion alone:
    v = 4.74047 / parallax * (mu_alpha* e_alpha + mu_delta e_delta)
as a 3D Cartesian vector with zero radial component. This UNDERESTIMATES true
|dv| and adds scatter, but it is computed identically for anomalies and
controls, so the comparison is fair. It can only dilute a real signal, never
manufacture one.
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
# Bound systems have dv ~ 0 BY CONSTRUCTION, and anomalies are enriched in
# multiplicity (paper Sec 5.5), so including them manufactures an enormous
# fake excess at low dv -- 4301 sigma on the first run. The widest bound wide
# binaries reach ~1 pc, so pairs closer than this are excluded outright.
MIN_SEP_PC = 1.0
DV_BINS = np.array([2, 5, 10, 20, 40, 80])


def velocity_vectors(d: pd.DataFrame) -> np.ndarray:
    ra = np.radians(d["ra"].to_numpy(float))
    dec = np.radians(d["dec"].to_numpy(float))
    plx = d["parallax_corr"].to_numpy(float)
    k = 4.74047 / np.maximum(plx, 1e-6)
    va = k * d["pmra"].to_numpy(float)
    vd = k * d["pmdec"].to_numpy(float)
    e_a = np.column_stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)])
    e_d = np.column_stack([-np.sin(dec) * np.cos(ra),
                           -np.sin(dec) * np.sin(ra), np.cos(dec)])
    return va[:, None] * e_a + vd[:, None] * e_d


def xyz_of(d: pd.DataFrame) -> np.ndarray:
    l = np.radians(d["l"].to_numpy(float))
    b = np.radians(d["b"].to_numpy(float))
    r = d["dist_pc"].to_numpy(float)
    return np.column_stack([r * np.cos(b) * np.cos(l),
                            r * np.cos(b) * np.sin(l), r * np.sin(b)])


def dv_profile(xyz: np.ndarray, v: np.ndarray, r_pc: float) -> tuple:
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(r_pc, output_type="ndarray")
    if len(pairs) == 0:
        return np.zeros(len(DV_BINS)), 0
    sep = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
    keep = sep > MIN_SEP_PC          # drop gravitationally bound systems
    pairs = pairs[keep]
    if len(pairs) == 0:
        return np.zeros(len(DV_BINS)), 0
    dv = np.linalg.norm(v[pairs[:, 0]] - v[pairs[:, 1]], axis=1)
    counts = np.array([np.count_nonzero(dv < b) for b in DV_BINS], dtype=float)
    return counts, len(pairs)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-control", type=int, default=40)
    args = ap.parse_args()
    rng = np.random.default_rng(20260816)

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["ra", "dec", "l", "b", "dist_pc", "residual",
                                 "parallax_corr", "pmra", "pmdec", "A_0",
                                 "M_Ks", "wise_w1mpro", "wise_w2mpro",
                                 "wise_w1mpro_error", "wise_w2mpro_error"])
    r = d["residual"].to_numpy(float)
    sig = st.robust_sigma(r)
    med = float(np.median(r))
    ok = np.isfinite(d["pmra"]) & np.isfinite(d["pmdec"]) & np.isfinite(d["dist_pc"])
    d = d[ok].reset_index(drop=True)
    r = r[ok.to_numpy()]

    anom = r > med + K_SIGMA * sig
    n_anom = int(anom.sum())
    print(f"{len(d):,} stars, {n_anom:,} anomalies at +{K_SIGMA}sigma")

    xyz = xyz_of(d)
    v = velocity_vectors(d)
    speed = np.linalg.norm(v, axis=1)
    print(f"tangential speed: median {np.median(speed):.1f} km/s, "
          f"16-84% {np.percentile(speed,16):.1f}-{np.percentile(speed,84):.1f}")

    obs_counts, n_obs_pairs = dv_profile(xyz[anom], v[anom], R_PC)
    print(f"\nanomaly-anomaly pairs within {R_PC:.0f} pc: {n_obs_pairs:,}")

    # Controls matched in distance AND local sky density. Distance alone is
    # not enough: anomalies sit preferentially in crowded fields (that is what
    # produces them), so a distance-only control has a different clustering
    # amplitude and the comparison is not like for like.
    dist = d["dist_pc"].to_numpy(float)
    import healpy as hp
    pix = hp.ang2pix(32, np.radians(90 - d["b"].to_numpy(float)),
                     np.radians(d["l"].to_numpy(float)), nest=True)
    dens = np.bincount(pix, minlength=hp.nside2npix(32))[pix].astype(float)

    # THE DIAGNOSTIC THAT MATTERS. If anomalies are kinematically COLD -- low
    # |v| -- then low pairwise dv follows trivially and means nothing. Young
    # stars are kinematically cold (they have not been disc-heated) and our
    # positive tail is enriched in young objects (YSOs, pre-MS, active stars,
    # paper Sec 5.4). So |v| must be matched, not just distance and density.
    print(f"\n  |v| anomalies : median {np.median(speed[anom]):6.1f} km/s")
    print(f"  |v| field     : median {np.median(speed[~anom]):6.1f} km/s")
    print(f"  fraction |v| < 20 km/s: anomalies "
          f"{np.mean(speed[anom] < 20):.3f}, field {np.mean(speed[~anom] < 20):.3f}")

    d_edges = np.percentile(dist[anom], np.linspace(0, 100, 9))
    rho_edges = np.percentile(dens[anom], np.linspace(0, 100, 5))
    v_edges = np.percentile(speed[anom], np.linspace(0, 100, 9))
    idx_all = np.arange(len(d))
    ctrl_counts, ctrl_pairs = [], []
    for _ in range(args.n_control):
        pick = []
        for a in range(len(d_edges) - 1):
            for c_ in range(len(rho_edges) - 1):
                for w in range(len(v_edges) - 1):
                    sel_a = ((dist[anom] >= d_edges[a]) & (dist[anom] < d_edges[a + 1])
                             & (dens[anom] >= rho_edges[c_])
                             & (dens[anom] < rho_edges[c_ + 1])
                             & (speed[anom] >= v_edges[w])
                             & (speed[anom] < v_edges[w + 1]))
                    want = int(sel_a.sum())
                    if not want:
                        continue
                    pool = idx_all[(dist >= d_edges[a]) & (dist < d_edges[a + 1])
                                   & (dens >= rho_edges[c_])
                                   & (dens < rho_edges[c_ + 1])
                                   & (speed >= v_edges[w]) & (speed < v_edges[w + 1])
                                   & ~anom]
                    if len(pool):
                        pick.append(rng.choice(pool, min(want, len(pool)),
                                               replace=False))
        sel = np.concatenate(pick)
        c, npr = dv_profile(xyz[sel], v[sel], R_PC)
        ctrl_counts.append(c)
        ctrl_pairs.append(npr)
    ctrl_counts = np.array(ctrl_counts)
    ctrl_pairs = np.array(ctrl_pairs, dtype=float)

    print(f"control sets: {args.n_control}, mean pairs "
          f"{ctrl_pairs.mean():,.0f}\n")

    rows = []
    for i, b in enumerate(DV_BINS):
        f_obs = obs_counts[i] / max(n_obs_pairs, 1)
        f_ctl = ctrl_counts[:, i] / np.maximum(ctrl_pairs, 1)
        mu, sd = f_ctl.mean(), f_ctl.std(ddof=1)
        z = (f_obs - mu) / sd if sd > 0 else np.nan
        rows.append({"dv_max_kms": b, "n_anom_pairs_below": int(obs_counts[i]),
                     "frac_anomaly": f_obs, "frac_control": mu,
                     "control_sd": sd, "ratio": f_obs / mu if mu > 0 else np.nan,
                     "n_sigma": z})
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda v: f"{v:12.6g}"))

    worst = t.loc[t["n_sigma"].abs().idxmax()]
    verdict = ("NO kinematic clustering of anomalies beyond the matched control"
               if abs(worst["n_sigma"]) < 4 else
               f"EXCESS at dv < {worst['dv_max_kms']:.0f} km/s "
               f"({worst['n_sigma']:+.1f} sigma) -- INVESTIGATE")
    print(f"\nVERDICT: {verdict}")

    out = {"tag": args.tag, "k_sigma": K_SIGMA, "r_pc": R_PC,
           "n_stars": int(len(d)), "n_anomalies": n_anom,
           "n_anomaly_pairs": int(n_obs_pairs),
           "table": t.to_dict(orient="records"), "verdict": verdict}
    (cfg.RESULT_DIR / f"velocity_clustering_{args.tag}.json").write_text(
        json.dumps(out, indent=2, default=float))
    t.to_csv(cfg.RESULT_DIR / f"velocity_clustering_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
