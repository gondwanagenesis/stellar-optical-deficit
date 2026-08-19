#!/usr/bin/env python
"""Search R: normal in every observable, impossible in combination.

    run.sh scripts/73_searchR_joint_anomaly.py --tag primary

THE ADVERSARIAL ARGUMENT
------------------------
Every channel in this project -- all twenty-five of them -- thresholds on one
or two observables. Optical deficit. Infrared excess. Spectral slope. RUWE.
Emissivity index. Each asks whether a single number is extreme.

Anything that tuned each number individually to sit inside the normal range
would pass all of them, and would have to be found some other way.

Making a MARGINAL distribution look ordinary is a tractable problem: it is one
number at a time, and the target is public. Making the JOINT distribution look
ordinary is a different problem. Stellar observables are correlated by physics
nobody designed -- temperature sets colour sets luminosity sets radius sets
surface gravity, and the infrared colours follow from the same atmosphere that
produced the optical ones. To pass a joint test you must reproduce every one of
those correlations simultaneously, including ones the observer has not thought
to check.

So this search deliberately DISCARDS every star that is already a marginal
outlier -- those are what the other channels are for -- and asks which of the
remaining, individually unremarkable stars sit where the correlation structure
says nothing should be.

METHOD
------
Robust standardisation, then the Mahalanobis distance from the bulk covariance,
which is the natural metric for "far from the manifold" when the manifold is
locally linear. Two guards make the result meaningful:

  the marginal veto   any star beyond 3 sigma on ANY single feature is removed
                      BEFORE scoring, so a high score cannot be inherited from
                      one extreme axis. What survives is anomalous only in
                      combination.

  the held-out null   the covariance is fitted on one random half and the score
                      computed on the other, so the score distribution of
                      ordinary stars is measured rather than assumed. Stellar
                      data is not Gaussian and a chi-squared expectation would
                      over-trigger badly.

WHAT A DETECTION WOULD LOOK LIKE
--------------------------------
Not a long tail -- that is what non-Gaussianity produces. A genuine population
would be a CLUMP: objects at similar joint-anomaly scores, sharing which
features drive the score, and ideally clustered in space. A smooth tail with
scattered feature attributions is the manifold being curved.
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
log = logging.getLogger("searchR")

# Features spanning photometry, astrometry, atmosphere and data quality.
# Deliberately broad: the point is to include correlations nobody checks.
BASE_COLS = [
    "source_id", "l", "b", "dist_pc", "A_0", "M_Ks", "residual",
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "phot_bp_rp_excess_factor", "cstar_nsigma",
    "ruwe", "astrometric_excess_noise", "parallax_over_error",
    "ipd_frac_multi_peak", "ipd_frac_odd_win", "visibility_periods_used",
    "phot_g_mean_flux_over_error", "phot_bp_n_obs", "phot_rp_n_obs",
    "teff_gspphot", "logg_gspphot", "mh_gspphot", "ebpminrp_gspphot",
    "tmass_j_m", "tmass_h_m", "tmass_ks_m", "tmass_xm_dist",
    "wise_w1mpro", "wise_w2mpro", "wise_w3mpro", "wise_w4mpro",
    "wise_xm_dist",
]

MARGINAL_VETO_SIGMA = 3.0
RNG_SEED = 73_2026


def build_features(d: pd.DataFrame):
    """Physically meaningful combinations, not raw magnitudes.

    Raw apparent magnitudes encode distance and would make the search find
    nearby stars. Colours and absolute quantities are what carry the physics.
    """
    f = {}
    g = d["phot_g_mean_mag"].to_numpy(float)
    bp = d["phot_bp_mean_mag"].to_numpy(float)
    rp = d["phot_rp_mean_mag"].to_numpy(float)
    j = d["tmass_j_m"].to_numpy(float)
    h = d["tmass_h_m"].to_numpy(float)
    ks = d["tmass_ks_m"].to_numpy(float)
    w1 = d["wise_w1mpro"].to_numpy(float)
    w2 = d["wise_w2mpro"].to_numpy(float)
    w3 = d["wise_w3mpro"].to_numpy(float)

    f["BP_RP"] = bp - rp
    f["BP_G"] = bp - g
    f["G_RP"] = g - rp
    f["G_J"] = g - j
    f["J_H"] = j - h
    f["H_Ks"] = h - ks
    f["Ks_W1"] = ks - w1
    f["W1_W2"] = w1 - w2
    f["W2_W3"] = w2 - w3
    f["M_Ks"] = d["M_Ks"].to_numpy(float)
    f["deficit"] = d["residual"].to_numpy(float)
    f["cstar"] = d["cstar_nsigma"].to_numpy(float)
    f["ruwe"] = d["ruwe"].to_numpy(float)
    f["ast_noise"] = np.log10(
        np.maximum(d["astrometric_excess_noise"].to_numpy(float), 1e-3))
    f["teff"] = d["teff_gspphot"].to_numpy(float)
    f["logg"] = d["logg_gspphot"].to_numpy(float)
    f["mh"] = d["mh_gspphot"].to_numpy(float)

    nbp = d["phot_bp_n_obs"].fillna(0).to_numpy(float)
    fove = d["phot_g_mean_flux_over_error"].to_numpy(float)
    f["amp"] = np.log10(np.maximum(
        np.sqrt(np.where(nbp > 0, 9 * nbp, np.nan)) / np.maximum(fove, 1e-9),
        1e-6))
    return pd.DataFrame(f)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--top", type=int, default=300)
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    import pyarrow.parquet as pq
    path = cfg.DERIVED_DIR / f"{args.tag}_resid.parquet"
    have = set(pq.read_schema(path).names)
    use = [c for c in BASE_COLS if c in have]
    log.info("using %d of %d requested columns", len(use), len(BASE_COLS))
    d = pd.read_parquet(path, columns=use)
    log.info("loaded %d stars", len(d))

    X = build_features(d)
    feat = list(X.columns)
    log.info("feature space: %d dimensions -> %s", len(feat), ", ".join(feat))

    complete = X.notna().all(axis=1).to_numpy()
    log.info("with every feature measured: %d (%.1f%%)",
             int(complete.sum()), 100 * complete.mean())
    if complete.sum() < 5000:
        log.error("too few complete rows")
        return 1

    Xc = X[complete].to_numpy(float)
    dc = d[complete].reset_index(drop=True)

    # ---- robust standardisation ------------------------------------------
    med = np.median(Xc, axis=0)
    mad = np.median(np.abs(Xc - med), axis=0) * 1.4826
    mad[mad <= 0] = 1.0
    Z = (Xc - med) / mad

    # ---- the marginal veto -----------------------------------------------
    # Anything extreme on a single axis is already the business of another
    # channel. Removing it guarantees the survivors are anomalous ONLY in
    # combination.
    max_abs = np.max(np.abs(Z), axis=1)
    inlier = max_abs < MARGINAL_VETO_SIGMA
    log.info("marginally unremarkable on ALL %d axes (|z| < %.0f): %d (%.1f%%)",
             len(feat), MARGINAL_VETO_SIGMA, int(inlier.sum()),
             100 * inlier.mean())

    Zi = Z[inlier]
    di = dc[inlier].reset_index(drop=True)

    # ---- train / test split so the null is measured, not assumed ---------
    n = len(Zi)
    perm = rng.permutation(n)
    half = n // 2
    tr, te = perm[:half], perm[half:]

    cov = np.cov(Zi[tr], rowvar=False)
    cov += np.eye(cov.shape[0]) * 1e-6
    inv = np.linalg.pinv(cov)
    mu = np.median(Zi[tr], axis=0)

    def mahal(A):
        D = A - mu
        return np.einsum("ij,jk,ik->i", D, inv, D)

    d2_tr = mahal(Zi[tr])
    d2_te = mahal(Zi[te])
    log.info("Mahalanobis^2 on held-out data: median %.1f, 99.9th pct %.1f, "
             "max %.1f", float(np.median(d2_te)),
             float(np.percentile(d2_te, 99.9)), float(d2_te.max()))

    # The training half defines the null; anything in the test half beyond
    # the training maximum is unprecedented rather than merely rare.
    thr_999 = float(np.percentile(d2_tr, 99.9))
    thr_max = float(d2_tr.max())
    n_beyond_999 = int((d2_te > thr_999).sum())
    n_beyond_max = int((d2_te > thr_max).sum())
    expected_999 = 0.001 * len(te)
    log.info("held-out beyond the training 99.9th pct: %d (expected %.0f)",
             n_beyond_999, expected_999)
    log.info("held-out beyond the training MAXIMUM   : %d", n_beyond_max)

    # ---- score everything and inspect the extreme tail -------------------
    d2_all = mahal(Zi)
    di = di.assign(joint_d2=d2_all)
    order = np.argsort(-d2_all)
    top = di.iloc[order[:args.top]].copy()

    # which features drive each extreme score?
    Dtop = Zi[order[:args.top]] - mu
    contrib = (Dtop @ inv) * Dtop          # per-feature contribution to d^2
    driver = np.array(feat)[np.argmax(contrib, axis=1)]
    top["top_driver"] = driver
    drv = pd.Series(driver).value_counts()
    log.info("")
    log.info("features driving the %d most jointly-anomalous stars:", args.top)
    for k, v in drv.head(8).items():
        log.info("  %-12s %4d (%.0f%%)", k, v, 100 * v / len(top))

    # ---- clump or tail? ---------------------------------------------------
    # A real population is a clump at similar scores with a shared driver; a
    # curved manifold gives a smooth tail with scattered drivers.
    top_frac = float(drv.iloc[0] / len(top))
    ratio_1_10 = float(d2_all[order[0]] / d2_all[order[9]]) if len(order) > 10 else 1.0

    # spatial clustering of the extreme tail, permutation-tested
    from scipy.spatial import cKDTree
    lr = np.radians(di["l"].to_numpy(float))
    br = np.radians(di["b"].to_numpy(float))
    xyz = np.column_stack([np.cos(br) * np.cos(lr),
                           np.cos(br) * np.sin(lr), np.sin(br)])
    sel_top = order[:args.top]
    tree = cKDTree(xyz)
    r = 2 * np.sin(np.radians(2.0) / 2)          # 2 degree radius
    obs_pairs = int(sum(len(p) - 1 for p in
                        tree.query_ball_point(xyz[sel_top], r=r)) / 2)
    null_pairs = []
    for _ in range(200):
        s = rng.choice(len(di), size=args.top, replace=False)
        null_pairs.append(int(sum(len(p) - 1 for p in
                                  tree.query_ball_point(xyz[s], r=r)) / 2))
    npm, nps = float(np.mean(null_pairs)), float(np.std(null_pairs))
    z_clump = (obs_pairs - npm) / nps if nps > 0 else 0.0
    log.info("")
    log.info("spatial pairs within 2 deg among the top %d: %d "
             "(null %.1f +/- %.1f -> %+.1f sigma)",
             args.top, obs_pairs, npm, nps, z_clump)

    # ---- verdict ----------------------------------------------------------
    if n_beyond_max > 0 and z_clump > 5:
        verdict = (
            f"{n_beyond_max} held-out stars exceed anything in the training "
            f"half, and the extreme tail is spatially clustered at "
            f"{z_clump:+.1f} sigma. That combination -- unprecedented joint "
            f"scores AND spatial structure -- is what a real population looks "
            f"like rather than a curved manifold. Inspect the candidate list.")
    elif z_clump > 5:
        verdict = (
            f"The most jointly-anomalous stars are spatially clustered "
            f"({z_clump:+.1f} sigma) but none exceeds the training-half "
            f"maximum, so their scores are rare rather than unprecedented. "
            f"Clustering with ordinary scores usually means a real but known "
            f"population -- a cluster, an association, or a region of "
            f"systematically different data quality. Check what is at those "
            f"coordinates before anything else.")
    else:
        verdict = (
            f"NULL. Among {int(inlier.sum()):,} stars that are unremarkable on "
            f"every one of {len(feat)} observables individually, the joint "
            f"distribution contains no population sitting where the "
            f"correlation structure forbids: {n_beyond_999} exceed the "
            f"training 99.9th percentile against {expected_999:.0f} expected, "
            f"the extreme tail is not spatially clustered ({z_clump:+.1f} "
            f"sigma), and its scores are driven by {top_frac:.0%} a single "
            f"feature, which is the signature of a curved manifold rather "
            f"than a distinct class. Nothing here is hiding by being average "
            f"everywhere at once.")

    print(f"\n{'='*74}")
    print("SEARCH R: NORMAL IN EVERY AXIS, ANOMALOUS IN COMBINATION")
    print(f"{'='*74}")
    print(f"  feature dimensions                  : {len(feat)}")
    print(f"  stars with all features measured    : {int(complete.sum()):,}")
    print(f"  unremarkable on every axis          : {int(inlier.sum()):,}")
    print(f"  held-out beyond training 99.9 pct   : {n_beyond_999} "
          f"(expected {expected_999:.0f})")
    print(f"  held-out beyond training maximum    : {n_beyond_max}")
    print(f"  spatial clustering of extreme tail  : {z_clump:+.1f} sigma")
    print(f"  dominant score driver               : {drv.index[0]} "
          f"({top_frac:.0%})")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchR_joint_anomaly_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "features": feat,
        "n_complete": int(complete.sum()),
        "n_marginal_inliers": int(inlier.sum()),
        "marginal_veto_sigma": MARGINAL_VETO_SIGMA,
        "d2_median_heldout": float(np.median(d2_te)),
        "d2_train_p999": thr_999,
        "d2_train_max": thr_max,
        "n_heldout_beyond_p999": n_beyond_999,
        "n_expected_beyond_p999": float(expected_999),
        "n_heldout_beyond_train_max": n_beyond_max,
        "spatial_clump_z": float(z_clump),
        "top_drivers": {str(k): int(v) for k, v in drv.head(10).items()},
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)

    cols = [c for c in ["source_id", "l", "b", "dist_pc", "joint_d2",
                        "top_driver", "M_Ks", "residual", "bp_rp", "ruwe",
                        "teff_gspphot", "phot_g_mean_mag"]
            if c in top.columns]
    top[cols].to_csv(cfg.RESULT_DIR / f"searchR_candidates_{args.tag}.csv",
                     index=False)
    log.info("wrote the top %d joint anomalies", len(top))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
