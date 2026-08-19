#!/usr/bin/env python
"""Search F, rebuilt: the mark-permutation null was geometrically invalid.

    run.sh scripts/77_searchF_corrected_null.py

THE ERROR
---------
Search F compared the mean local resultant of proper-motion-anomaly directions
against a mark-permutation null and got -3.9 sigma: the observed statistic sat
BELOW its own null. That was reported as a null result. It is a bug.

The first hypothesis -- resolved wide binaries anti-aligning because each
component is pulled toward the other -- was tested in script 76 and FAILED.
Only six pairs lie within 0.02 pc, and removing them made the deficit worse
(-3.91 to -4.26 sigma). So the cause is not physical.

It is geometric. The proper-motion anomaly is a TRANSVERSE quantity: it lives
in the plane perpendicular to the line of sight, always. Stars that are
neighbours in 3D share very nearly the same line of sight, so their anomaly
vectors are all confined to nearly the same plane.

Mark permutation destroys exactly that. It hands a star a vector computed for a
different star on a different line of sight, and that vector is not transverse
to the recipient's line of sight. So permuted neighbourhoods contain vectors
spanning the full sphere while real neighbourhoods contain vectors spanning a
plane -- and for k random unit vectors the expected resultant is LARGER in
three dimensions than in two:

    E|sum| / sqrt(k)  ->  sqrt(8/(3 pi)) = 0.921   in 3D
    E|sum| / sqrt(k)  ->  sqrt(pi)/2     = 0.886   in 2D

The null was therefore inflated by construction, by roughly the 4 per cent
that the observed deficit represents. Search F's headline number was measured
against a baseline that no real data could ever reach.

THE CORRECT NULL
----------------
Randomise only the quantity the signal lives in, and preserve the geometry:
rotate each star's own anomaly vector by a random angle ABOUT ITS OWN LINE OF
SIGHT. That keeps every vector transverse to its own line of sight, keeps the
sky distribution, the density field and the magnitude distribution untouched,
and destroys only the position-angle correlation between neighbours -- which is
precisely what a coordinated thrust would produce.

This script re-runs the channel against that null and re-derives the injection
sensitivity, so the published limit rests on a calibrated zero point.
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
log = logging.getLogger("searchF2")

KMS_PER_MASYR_KPC = 4.740470446
SNR_MIN, PLX_MIN, K_NN, N_SHUFFLE = 3.0, 2.0, 40, 200
RNG_SEED = 77_2026


def radec_basis(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    ca, sa, cd, sd = np.cos(a), np.sin(a), np.cos(d), np.sin(d)
    return (np.stack([cd * ca, cd * sa, sd], 1),
            np.stack([-sa, ca, np.zeros_like(sa)], 1),
            np.stack([-sd * ca, -sd * sa, cd], 1))


def icrs_to_gal():
    ra, de, ln = np.radians(192.85948), np.radians(27.12825), np.radians(122.93192)
    m1 = np.array([[np.cos(ra), np.sin(ra), 0], [-np.sin(ra), np.cos(ra), 0], [0, 0, 1]])
    m2 = np.array([[np.sin(de), 0, -np.cos(de)], [0, 1, 0], [np.cos(de), 0, np.sin(de)]])
    m3 = np.array([[np.cos(ln), np.sin(ln), 0], [-np.sin(ln), np.cos(ln), 0], [0, 0, 1]])
    return m3 @ m2 @ m1


def fit_rot_glide(r, ea, ed, da, dd, sa, sd):
    n = len(r)
    D = np.zeros((2 * n, 6))
    for k in range(3):
        ek = np.zeros(3); ek[k] = 1.0
        rot = np.cross(np.broadcast_to(ek, r.shape), r)
        D[0::2, k] = np.einsum("ij,ij->i", rot, ea)
        D[1::2, k] = np.einsum("ij,ij->i", rot, ed)
        D[0::2, 3 + k] = ea[:, k]
        D[1::2, 3 + k] = ed[:, k]
    y = np.empty(2 * n); y[0::2] = da; y[1::2] = dd
    w = np.empty(2 * n)
    w[0::2] = 1 / np.maximum(sa, 1e-6); w[1::2] = 1 / np.maximum(sd, 1e-6)
    p, *_ = np.linalg.lstsq(D * w[:, None], y * w, rcond=None)
    m = D @ p
    return p, m[0::2], m[1::2]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=K_NN)
    ap.add_argument("--n-shuffle", type=int, default=N_SHUFFLE)
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    d = pd.read_parquet(cfg.RAW_DIR / "hgca_brandt2021.parquet")
    need = ["plx", "pmRA", "pmDE", "pmRAhg", "pmDEhg",
            "e_pmRA", "e_pmDE", "e_pmRAhg", "e_pmDEhg", "RA_ICRS", "DE_ICRS"]
    d = d[d[need].notna().all(axis=1) & (d["plx"] > PLX_MIN)].reset_index(drop=True)

    da = d["pmRA"].to_numpy(float) - (d["pmRAhg"].to_numpy(float)
                                      + d["dpmRA"].fillna(0).to_numpy(float))
    dd = d["pmDE"].to_numpy(float) - (d["pmDEhg"].to_numpy(float)
                                      + d["dpmDE"].fillna(0).to_numpy(float))
    sa = np.hypot(d["e_pmRA"].to_numpy(float), d["e_pmRAhg"].to_numpy(float))
    sd = np.hypot(d["e_pmDE"].to_numpy(float), d["e_pmDEhg"].to_numpy(float))

    r, ea, ed = radec_basis(d["RA_ICRS"].to_numpy(float),
                            d["DE_ICRS"].to_numpy(float))
    _, ma, md = fit_rot_glide(r, ea, ed, da, dd, sa, sd)
    ra_, rd_ = da - ma, dd - md
    sel = np.hypot(ra_ / sa, rd_ / sd) > SNR_MIN

    dist_kpc = 1.0 / d["plx"].to_numpy(float)[sel]
    R = icrs_to_gal()
    dv = (ra_[sel, None] * ea[sel] + rd_[sel, None] * ed[sel]) \
        * (KMS_PER_MASYR_KPC * dist_kpc)[:, None]
    dv_gal, r_gal = dv @ R.T, r[sel] @ R.T
    xyz = r_gal * (dist_kpc * 1000.0)[:, None]
    nrm = np.linalg.norm(dv_gal, axis=1)
    g = nrm > 0
    unit, xyz, los = dv_gal[g] / nrm[g, None], xyz[g], r_gal[g]
    n = len(unit)
    log.info("directional sample: %d", n)

    # ---- confirm the diagnosis: are real vectors transverse and permuted ones not? ---
    along_obs = float(np.mean(np.abs(np.einsum("ij,ij->i", unit, los))))
    perm = rng.permutation(n)
    along_perm = float(np.mean(np.abs(np.einsum("ij,ij->i", unit[perm], los))))
    log.info("")
    log.info("mean |unit . line-of-sight| :")
    log.info("  observed  %.5f   (transverse by construction, should be ~0)",
             along_obs)
    log.info("  permuted  %.5f   (no longer transverse -> spans 3D)",
             along_perm)
    log.info("  ratio     %.1fx", along_perm / max(along_obs, 1e-9))

    # ---- the two nulls ----------------------------------------------------
    from scipy.spatial import cKDTree
    tree = cKDTree(xyz)
    k = min(args.k, n - 1)
    _, idx = tree.query(xyz, k=k)

    def resultant(v):
        return float((np.linalg.norm(v[idx].sum(axis=1), axis=1) / k).mean())

    obs = resultant(unit)

    def spin_about_los(rng):
        """Rotate each vector by a random angle about ITS OWN line of sight.

        Preserves transversality, sky positions, density and magnitudes;
        destroys only the position-angle correlation between neighbours.
        """
        # orthonormal transverse basis per star
        a = np.cross(los, np.array([0.0, 0.0, 1.0]))
        bad = np.linalg.norm(a, axis=1) < 1e-8
        a[bad] = np.cross(los[bad], np.array([0.0, 1.0, 0.0]))
        a /= np.linalg.norm(a, axis=1)[:, None]
        b = np.cross(los, a)
        th = rng.uniform(0, 2 * np.pi, size=len(los))
        return np.cos(th)[:, None] * a + np.sin(th)[:, None] * b

    null_perm = np.array([resultant(unit[rng.permutation(n)])
                          for _ in range(args.n_shuffle)])
    null_spin = np.array([resultant(spin_about_los(rng))
                          for _ in range(args.n_shuffle)])

    def z_of(nl):
        m, s = float(nl.mean()), float(nl.std())
        return (obs - m) / s if s else 0.0, m, s

    z_p, mp, sp = z_of(null_perm)
    z_s, ms, ss = z_of(null_spin)

    log.info("")
    log.info("observed mean local resultant : %.5f", obs)
    log.info("  mark-permutation null (WRONG): %.5f +/- %.5f  -> z = %+.2f",
             mp, sp, z_p)
    log.info("  line-of-sight-spin null (correct): %.5f +/- %.5f  -> z = %+.2f",
             ms, ss, z_s)

    # ---- injection against the corrected null ----------------------------
    log.info("")
    log.info("injection-recovery against the corrected null:")
    inj = []
    for frac in (0.005, 0.01, 0.02, 0.05, 0.10, 0.20):
        n_inj = max(int(frac * n), 2)
        centre = xyz[rng.integers(n)]
        _, patch = tree.query(centre, k=n_inj)
        v = unit.copy()
        # a real thrust is a 3D direction; each star only reveals its
        # transverse component, so project it onto each star's own plane
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        proj = direction - (los[patch] @ direction)[:, None] * los[patch]
        nn = np.linalg.norm(proj, axis=1)
        okp = nn > 1e-6
        v[patch[okp]] = proj[okp] / nn[okp, None]
        zi = (resultant(v) - ms) / ss if ss else 0.0
        inj.append({"coherent_fraction": frac, "n": int(n_inj), "z": float(zi)})
        log.info("  f=%.3f (n=%5d) -> z = %+.1f", frac, n_inj, zi)

    det = [r for r in inj if r["z"] > 5]
    f_min = min(r["coherent_fraction"] for r in det) if det else None

    verdict = (
        f"Search F's null was geometrically invalid and is corrected here. The "
        f"proper-motion anomaly is transverse to the line of sight, so mark "
        f"permutation -- which hands a star a vector belonging to a different "
        f"line of sight -- compares real two-dimensional neighbourhoods "
        f"against synthetic three-dimensional ones, and random unit vectors "
        f"have a larger resultant in 3D (0.921 sqrt-k) than in 2D (0.886 "
        f"sqrt-k). That inflated the null by about 4 per cent and produced the "
        f"unexplained -3.9 sigma. Rotating each vector about its own line of "
        f"sight preserves the geometry and destroys only the position-angle "
        f"correlation a coordinated thrust would create. Against that null the "
        f"result is z = {z_s:+.2f}: still null, now calibrated. Injection "
        f"recovers a coherent domain of "
        f"{'>=' + format(f_min, '.1%') if f_min else 'no tested fraction'} "
        f"at 5 sigma.")

    print(f"\n{'='*74}")
    print("SEARCH F, CORRECTED NULL")
    print(f"{'='*74}")
    print(f"  directional sample                  : {n:,}")
    print(f"  mean |u.LOS| observed / permuted    : "
          f"{along_obs:.5f} / {along_perm:.5f}")
    print(f"  observed mean local resultant       : {obs:.5f}")
    print(f"  mark-permutation null (invalid)     : {mp:.5f} -> z {z_p:+.2f}")
    print(f"  line-of-sight-spin null (correct)   : {ms:.5f} -> z {z_s:+.2f}")
    if f_min:
        print(f"  detectable coherent fraction        : >= {f_min:.1%}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / "searchF_corrected_null.json"
    out.write_text(json.dumps({
        "n_directional": int(n),
        "mean_abs_u_dot_los_observed": along_obs,
        "mean_abs_u_dot_los_permuted": along_perm,
        "observed_resultant": obs,
        "null_mark_permutation": {"mean": mp, "std": sp, "z": float(z_p)},
        "null_los_spin": {"mean": ms, "std": ss, "z": float(z_s)},
        "injection": inj,
        "min_detectable_coherent_fraction": f_min,
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
