#!/usr/bin/env python
"""Why was Search F's coherence BELOW its own null?

    run.sh scripts/76_searchF_deficit_diagnosis.py

THE LOOSE END
-------------
Search F measured the mean local resultant of proper-motion-anomaly directions
at 0.14074 against a mark-permutation null of 0.14573 +/- 0.00128, i.e.
-3.89 sigma. It was reported as a null and it is not one.

Under the null hypothesis the observed statistic and the permuted statistic
have the same expectation. A significant DEFICIT of alignment is exactly as
anomalous as an excess, and an unexplained 3.9-sigma deficit in the headline
statistic of a channel is not something to leave sitting there.

THE HYPOTHESIS
--------------
Resolved wide binaries. If both components of a pair are in the Hipparcos-Gaia
catalogue, each one's proper-motion anomaly is dominated by the pull of the
other, so the two vectors point TOWARD each other -- anti-aligned. A
neighbourhood containing such a pair therefore has LOWER resultant length than
random, because the pair partly cancels.

Mark permutation destroys that pairing: it reassigns anomaly vectors to
positions at random, so the anti-alignment disappears and the permuted
statistic sits higher. The observed value being below the null is then not a
deficit of signal but the presence of real, mundane, anti-correlated physics
that the null deliberately erases.

THE TEST
--------
Directly measure the alignment of the CLOSEST pairs as a function of
separation. If the hypothesis holds:

  * pairs separated by less than a few thousand AU should be systematically
    anti-aligned, with mean cos(angle) < 0;
  * the anti-alignment must weaken with separation and vanish for unbound
    neighbours;
  * removing one component of each close pair should raise the global
    resultant back toward the null.

If instead the closest pairs show no anti-alignment, the hypothesis is wrong
and the deficit is a bug in the statistic, which would put Search F's headline
number in question.
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
log = logging.getLogger("searchF_diag")

KMS_PER_MASYR_KPC = 4.740470446
SNR_MIN = 3.0
PLX_MIN = 2.0
K_NEIGHBOURS = 40
N_SHUFFLE = 200
RNG_SEED = 76_2026


def radec_basis(ra_deg, dec_deg):
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    ca, sa, cd, sd = np.cos(a), np.sin(a), np.cos(d), np.sin(d)
    r = np.stack([cd * ca, cd * sa, sd], axis=1)
    e_a = np.stack([-sa, ca, np.zeros_like(sa)], axis=1)
    e_d = np.stack([-sd * ca, -sd * sa, cd], axis=1)
    return r, e_a, e_d


def icrs_to_gal():
    ra_ngp, dec_ngp, l_ncp = (np.radians(192.85948), np.radians(27.12825),
                              np.radians(122.93192))
    sd, cd = np.sin(dec_ngp), np.cos(dec_ngp)
    sa, ca = np.sin(ra_ngp), np.cos(ra_ngp)
    sl, cl = np.sin(l_ncp), np.cos(l_ncp)
    m1 = np.array([[ca, sa, 0], [-sa, ca, 0], [0, 0, 1]])
    m2 = np.array([[sd, 0, -cd], [0, 1, 0], [cd, 0, sd]])
    m3 = np.array([[cl, sl, 0], [-sl, cl, 0], [0, 0, 1]])
    return m3 @ m2 @ m1


def fit_rotation_glide(r_hat, e_a, e_d, da, dd, sa_, sd_):
    n = len(r_hat)
    design = np.zeros((2 * n, 6))
    for k in range(3):
        ek = np.zeros(3); ek[k] = 1.0
        rot = np.cross(np.broadcast_to(ek, r_hat.shape), r_hat)
        design[0::2, k] = np.einsum("ij,ij->i", rot, e_a)
        design[1::2, k] = np.einsum("ij,ij->i", rot, e_d)
        design[0::2, 3 + k] = e_a[:, k]
        design[1::2, 3 + k] = e_d[:, k]
    y = np.empty(2 * n); y[0::2] = da; y[1::2] = dd
    w = np.empty(2 * n)
    w[0::2] = 1.0 / np.maximum(sa_, 1e-6)
    w[1::2] = 1.0 / np.maximum(sd_, 1e-6)
    p, *_ = np.linalg.lstsq(design * w[:, None], y * w, rcond=None)
    m = design @ p
    return p, m[0::2], m[1::2]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=K_NEIGHBOURS)
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    d = pd.read_parquet(cfg.RAW_DIR / "hgca_brandt2021.parquet")
    need = ["plx", "pmRA", "pmDE", "pmRAhg", "pmDEhg",
            "e_pmRA", "e_pmDE", "e_pmRAhg", "e_pmDEhg", "RA_ICRS", "DE_ICRS"]
    ok = d[need].notna().all(axis=1) & (d["plx"] > PLX_MIN)
    d = d[ok].reset_index(drop=True)
    log.info("HGCA within 500 pc: %d", len(d))

    dpa = d["dpmRA"].fillna(0).to_numpy(float)
    dpd = d["dpmDE"].fillna(0).to_numpy(float)
    da = d["pmRA"].to_numpy(float) - (d["pmRAhg"].to_numpy(float) + dpa)
    dd = d["pmDE"].to_numpy(float) - (d["pmDEhg"].to_numpy(float) + dpd)
    sa_ = np.hypot(d["e_pmRA"].to_numpy(float), d["e_pmRAhg"].to_numpy(float))
    sd_ = np.hypot(d["e_pmDE"].to_numpy(float), d["e_pmDEhg"].to_numpy(float))

    r_hat, e_a, e_d = radec_basis(d["RA_ICRS"].to_numpy(float),
                                  d["DE_ICRS"].to_numpy(float))
    _, ma, md = fit_rotation_glide(r_hat, e_a, e_d, da, dd, sa_, sd_)
    ra_, rd_ = da - ma, dd - md
    snr = np.hypot(ra_ / sa_, rd_ / sd_)
    sel = snr > SNR_MIN
    log.info("significant PMa: %d", int(sel.sum()))

    dist_kpc = 1.0 / d["plx"].to_numpy(float)[sel]
    scale = KMS_PER_MASYR_KPC * dist_kpc
    dv = (ra_[sel, None] * e_a[sel] + rd_[sel, None] * e_d[sel]) * scale[:, None]
    R = icrs_to_gal()
    dv_gal = dv @ R.T
    xyz = (r_hat[sel] * (dist_kpc * 1000.0)[:, None]) @ R.T
    nrm = np.linalg.norm(dv_gal, axis=1)
    good = nrm > 0
    unit = dv_gal[good] / nrm[good, None]
    xyz = xyz[good]
    log.info("directional sample: %d", len(unit))

    # ---- the decisive measurement: alignment vs pair separation ----------
    from scipy.spatial import cKDTree
    tree = cKDTree(xyz)
    dist_nn, idx_nn = tree.query(xyz, k=2)
    sep_pc = dist_nn[:, 1]
    partner = idx_nn[:, 1]
    cosang = np.einsum("ij,ij->i", unit, unit[partner])

    log.info("")
    log.info("=== alignment with the NEAREST neighbour, vs separation ===")
    log.info("(anti-alignment means cos < 0; random means cos ~ 0)")
    rows = []
    for lo, hi in [(0, 0.005), (0.005, 0.02), (0.02, 0.05), (0.05, 0.1),
                   (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 1e9)]:
        m = (sep_pc >= lo) & (sep_pc < hi)
        if m.sum() < 30:
            continue
        mc = float(np.mean(cosang[m]))
        se = float(np.std(cosang[m]) / np.sqrt(m.sum()))
        rows.append({"sep_lo_pc": lo, "sep_hi_pc": hi, "n": int(m.sum()),
                     "mean_cos": mc, "sem": se, "z": mc / se if se else 0.0})
        log.info("  sep %7.3f-%-8.3f pc  n=%6d  mean cos = %+.4f +/- %.4f "
                 "(%+.1f sigma)", lo, hi, int(m.sum()), mc, se,
                 mc / se if se else 0.0)

    # AU scale for the closest bin, which is where binaries live
    close = sep_pc < 0.02
    log.info("")
    log.info("pairs closer than 0.02 pc (~4100 AU): %d, mean cos = %+.4f",
             int(close.sum()), float(np.mean(cosang[close]))
             if close.sum() else np.nan)

    # ---- does removing one of each close pair restore the null? ----------
    def coherence(vecs, positions, k):
        t = cKDTree(positions)
        _, ii = t.query(positions, k=k)
        return float((np.linalg.norm(vecs[ii].sum(axis=1), axis=1) / k).mean())

    k = min(args.k, len(unit) - 1)
    obs_all = coherence(unit, xyz, k)
    null_all = []
    for _ in range(N_SHUFFLE):
        null_all.append(coherence(unit[rng.permutation(len(unit))], xyz, k))
    mu, sd = float(np.mean(null_all)), float(np.std(null_all))
    z_all = (obs_all - mu) / sd if sd else 0.0
    log.info("")
    log.info("ALL stars      : observed %.5f  null %.5f +/- %.5f  z=%+.2f",
             obs_all, mu, sd, z_all)

    # drop the higher-index member of every close pair
    drop = np.zeros(len(unit), dtype=bool)
    for i in np.where(close)[0]:
        j = partner[i]
        if not drop[i] and not drop[j]:
            drop[max(i, j)] = True
    keep = ~drop
    log.info("removing one member of each close pair: dropping %d",
             int(drop.sum()))

    u2, x2 = unit[keep], xyz[keep]
    k2 = min(args.k, len(u2) - 1)
    obs2 = coherence(u2, x2, k2)
    null2 = []
    for _ in range(N_SHUFFLE):
        null2.append(coherence(u2[rng.permutation(len(u2))], x2, k2))
    mu2, sd2 = float(np.mean(null2)), float(np.std(null2))
    z2 = (obs2 - mu2) / sd2 if sd2 else 0.0
    log.info("close pairs cut: observed %.5f  null %.5f +/- %.5f  z=%+.2f",
             obs2, mu2, sd2, z2)

    # ---- verdict ----------------------------------------------------------
    closest = rows[0] if rows else None
    anti = closest and closest["mean_cos"] < 0 and closest["z"] < -2

    if anti and abs(z2) < abs(z_all) * 0.6:
        verdict = (
            f"EXPLAINED, and it is physics rather than a bug. The nearest-"
            f"neighbour pairs are significantly ANTI-aligned at the closest "
            f"separations (mean cos = {closest['mean_cos']:+.4f}, "
            f"{closest['z']:+.1f} sigma), which is what resolved wide binaries "
            f"must do: each component's proper-motion anomaly is dominated by "
            f"the pull of the other, so the two vectors point toward each "
            f"other and partly cancel inside a neighbourhood. Mark permutation "
            f"destroys the pairing, so the permuted statistic sits higher and "
            f"the observed value falls below its own null. Removing one member "
            f"of each close pair moves the global statistic from "
            f"{z_all:+.2f} to {z2:+.2f} sigma. Search F's headline null "
            f"stands; the deficit was real, mundane, anti-correlated binary "
            f"physics that the null erases by construction.")
    elif anti:
        verdict = (
            f"PARTLY EXPLAINED. Close pairs are anti-aligned "
            f"({closest['mean_cos']:+.4f}, {closest['z']:+.1f} sigma) as wide "
            f"binaries require, but removing them only moves the global "
            f"statistic from {z_all:+.2f} to {z2:+.2f} sigma, so binary "
            f"cancellation is not the whole effect and something else also "
            f"suppresses local alignment.")
    else:
        verdict = (
            f"NOT EXPLAINED, and this is now a concern. The closest pairs show "
            f"no significant anti-alignment (mean cos = "
            f"{closest['mean_cos'] if closest else float('nan'):+.4f}), so "
            f"resolved binaries do not account for the deficit. An observed "
            f"statistic {z_all:+.2f} sigma BELOW its own permutation null "
            f"without a physical cause points at the statistic itself. Search "
            f"F's coherence estimator should be re-derived before its limit is "
            f"quoted.")

    print(f"\n{'='*72}")
    print("SEARCH F DIAGNOSIS: WHY WAS COHERENCE BELOW THE NULL?")
    print(f"{'='*72}")
    for r in rows:
        print(f"  sep {r['sep_lo_pc']:7.3f}-{r['sep_hi_pc']:<8.3f} pc  "
              f"n={r['n']:6d}  mean cos = {r['mean_cos']:+.4f} "
              f"({r['z']:+.1f} sigma)")
    print()
    print(f"  all stars       : z = {z_all:+.2f}")
    print(f"  close pairs cut : z = {z2:+.2f}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / "searchF_deficit_diagnosis.json"
    out.write_text(json.dumps({
        "n_directional": int(len(unit)),
        "separation_bins": rows,
        "z_all_stars": float(z_all),
        "z_close_pairs_removed": float(z2),
        "n_dropped": int(drop.sum()),
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
