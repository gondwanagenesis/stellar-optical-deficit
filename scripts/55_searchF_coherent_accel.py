#!/usr/bin/env python
"""Search F: coherent acceleration domains from the Hipparcos-Gaia PMa.

    run.sh scripts/55_searchF_coherent_accel.py

THE HYPOTHESIS
--------------
Light can be beamed, cooled, or spectrally shaped to hide. Momentum cannot.
A star being deliberately moved carries an acceleration that no photometric
concealment strategy removes.

The catch is that the proper-motion anomaly is already the standard way to
find unseen companions, so a large anomaly on one star means "binary", not
"engine". Magnitude is not the discriminant.

DIRECTION IS. An unseen companion pulls its primary along whatever orbit it
happens to occupy, so across a population the anomaly vectors are isotropic.
A coordinated programme of stellar migration is, by construction, not
isotropic: it would leave a spatially localised patch of ALIGNED acceleration
vectors. So the search statistic is local directional coherence, and the whole
game is removing the other things that produce coherence.

THE THREE CONFOUNDS, AND WHY ONLY TWO SURVIVE
---------------------------------------------
1. Hipparcos-to-Gaia frame rotation. The two catalogues' reference frames have
   a residual relative spin, which imprints a rigid-rotation pattern on every
   PMa in the sky. This is GLOBAL and coherent -- exactly our signal shape, at
   the largest scale. Removed here by fitting and subtracting the six-parameter
   rotation + glide (Mignard & Klioner 2012) before any coherence is measured.

2. Perspective acceleration. A star with radial velocity changes its proper
   motion secularly even with no companion. Brandt supplies this per star as
   (dpmRA, dpmDE); we subtract it.

3. The Galactic potential itself. Coherent, but at ~2e-10 m/s^2 it moves a
   star by 0.16 m/s over the 25 yr baseline, which at 100 pc is 3e-4 mas/yr --
   three orders of magnitude below the per-star noise. Not a confound at this
   precision, which is worth stating rather than assuming.

WHY THIS BEATS THE VELOCITY-SPACE SEARCH (channel 7)
----------------------------------------------------
Channel 7 found a 4301-sigma velocity-clustering excess that dissolved into
young comoving associations. Associations share a VELOCITY because they formed
together. They do not share an ACCELERATION: the cluster's internal potential
produces nothing measurable on a 25 yr baseline. Moving to the acceleration
domain removes the contaminant that defeated channel 7 outright.
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
log = logging.getLogger("searchF")

K_NEIGHBOURS = 40          # local neighbourhood for the coherence statistic
SNR_MIN = 3.0              # PMa significance to enter the directional sample
N_SHUFFLE = 200            # null realisations
PLX_MIN = 2.0              # 500 pc
RNG_SEED = 55_2026

# 1 mas/yr at 1 kpc == 4.7405 km/s
KMS_PER_MASYR_KPC = 4.740470446


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------

def radec_to_cartesian(ra_deg, dec_deg):
    """Unit vector and the local (east, north) tangent basis, all ICRS."""
    a = np.radians(ra_deg)
    d = np.radians(dec_deg)
    ca, sa, cd, sd = np.cos(a), np.sin(a), np.cos(d), np.sin(d)

    r_hat = np.stack([cd * ca, cd * sa, sd], axis=1)
    e_alpha = np.stack([-sa, ca, np.zeros_like(sa)], axis=1)
    e_delta = np.stack([-sd * ca, -sd * sa, cd], axis=1)
    return r_hat, e_alpha, e_delta


def icrs_to_galactic_matrix():
    """ICRS -> Galactic rotation matrix (standard, Hipparcos convention)."""
    ra_ngp = np.radians(192.85948)
    dec_ngp = np.radians(27.12825)
    l_ncp = np.radians(122.93192)

    sd, cd = np.sin(dec_ngp), np.cos(dec_ngp)
    sa, ca = np.sin(ra_ngp), np.cos(ra_ngp)
    sl, cl = np.sin(l_ncp), np.cos(l_ncp)

    # R = R_z(-l_ncp) R_y(dec_ngp - 90) R_z(ra_ngp)
    m1 = np.array([[ca, sa, 0], [-sa, ca, 0], [0, 0, 1]])
    m2 = np.array([[sd, 0, -cd], [0, 1, 0], [cd, 0, sd]])
    m3 = np.array([[cl, sl, 0], [-sl, cl, 0], [0, 0, 1]])
    return m3 @ m2 @ m1


# ---------------------------------------------------------------------------
# global rotation + glide removal
# ---------------------------------------------------------------------------

def fit_rotation_glide(r_hat, e_alpha, e_delta, dmu_a, dmu_d, sig_a, sig_d):
    """Fit the six-parameter rigid rotation + glide to the PMa field.

    A frame spin omega produces apparent proper motion omega x r_hat; a glide
    g produces the tangential projection of g. Both are built numerically from
    the three Cartesian basis directions rather than written out in spherical
    trigonometry, which is where sign errors live.

    Returns (params, model_a, model_d).
    """
    n = len(r_hat)
    design = np.zeros((2 * n, 6))

    for k in range(3):
        e_k = np.zeros(3)
        e_k[k] = 1.0
        # rotation: omega_k x r_hat
        rot = np.cross(np.broadcast_to(e_k, r_hat.shape), r_hat)
        design[0::2, k] = np.einsum("ij,ij->i", rot, e_alpha)
        design[1::2, k] = np.einsum("ij,ij->i", rot, e_delta)
        # glide: tangential projection of g_k
        design[0::2, 3 + k] = e_alpha[:, k]
        design[1::2, 3 + k] = e_delta[:, k]

    y = np.empty(2 * n)
    y[0::2] = dmu_a
    y[1::2] = dmu_d

    w = np.empty(2 * n)
    w[0::2] = 1.0 / np.maximum(sig_a, 1e-6)
    w[1::2] = 1.0 / np.maximum(sig_d, 1e-6)

    dw = design * w[:, None]
    yw = y * w
    params, *_ = np.linalg.lstsq(dw, yw, rcond=None)

    model = design @ params
    return params, model[0::2], model[1::2]


# ---------------------------------------------------------------------------
# coherence statistic
# ---------------------------------------------------------------------------

def coherence_field(xyz, unit_vecs, k, rng=None, shuffle=False):
    """Mean local resultant length of the unit vectors over k-NN neighbourhoods.

    Returns (per-star resultant lengths, global mean).
    """
    from scipy.spatial import cKDTree

    v = unit_vecs
    if shuffle:
        v = v[rng.permutation(len(v))]

    tree = cKDTree(xyz)
    _, idx = tree.query(xyz, k=k)

    # idx: (n, k). Resultant of each neighbourhood's unit vectors.
    summed = v[idx].sum(axis=1)
    r = np.linalg.norm(summed, axis=1) / k
    return r, float(r.mean())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=int, default=K_NEIGHBOURS)
    ap.add_argument("--snr-min", type=float, default=SNR_MIN)
    ap.add_argument("--n-shuffle", type=int, default=N_SHUFFLE)
    args = ap.parse_args()

    rng = np.random.default_rng(RNG_SEED)

    df = pd.read_parquet(cfg.RAW_DIR / "hgca_brandt2021.parquet")
    log.info("HGCA: %d stars", len(df))

    # ---- quality and volume ------------------------------------------------
    need = ["plx", "pmRA", "pmDE", "pmRAhg", "pmDEhg",
            "e_pmRA", "e_pmDE", "e_pmRAhg", "e_pmDEhg",
            "RA_ICRS", "DE_ICRS"]
    ok = df[need].notna().all(axis=1) & (df["plx"] > PLX_MIN)
    d = df[ok].reset_index(drop=True)
    log.info("within 500 pc with complete astrometry: %d", len(d))

    # ---- the proper-motion anomaly ----------------------------------------
    # Perspective acceleration is supplied by Brandt as a correction to the
    # long-baseline proper motion; apply it before differencing.
    dpm_a = d["dpmRA"].fillna(0.0).to_numpy(float)
    dpm_d = d["dpmDE"].fillna(0.0).to_numpy(float)

    dmu_a = d["pmRA"].to_numpy(float) - (d["pmRAhg"].to_numpy(float) + dpm_a)
    dmu_d = d["pmDE"].to_numpy(float) - (d["pmDEhg"].to_numpy(float) + dpm_d)

    sig_a = np.hypot(d["e_pmRA"].to_numpy(float), d["e_pmRAhg"].to_numpy(float))
    sig_d = np.hypot(d["e_pmDE"].to_numpy(float), d["e_pmDEhg"].to_numpy(float))

    log.info("raw PMa: median |dmu| = %.4f mas/yr, median sigma = %.4f",
             float(np.median(np.hypot(dmu_a, dmu_d))),
             float(np.median(np.hypot(sig_a, sig_d))))

    # ---- remove the frame rotation + glide --------------------------------
    r_hat, e_alpha, e_delta = radec_to_cartesian(
        d["RA_ICRS"].to_numpy(float), d["DE_ICRS"].to_numpy(float))

    params, mod_a, mod_d = fit_rotation_glide(
        r_hat, e_alpha, e_delta, dmu_a, dmu_d, sig_a, sig_d)

    log.info("global frame terms removed (mas/yr):")
    log.info("  rotation omega = [%+.4f %+.4f %+.4f]", *params[:3])
    log.info("  glide    g     = [%+.4f %+.4f %+.4f]", *params[3:])

    res_a = dmu_a - mod_a
    res_d = dmu_d - mod_d

    # ---- significance ------------------------------------------------------
    snr = np.hypot(res_a / sig_a, res_d / sig_d)
    sel = snr > args.snr_min
    log.info("PMa significant at >%.1f sigma: %d (%.1f%%)",
             args.snr_min, int(sel.sum()), 100 * sel.mean())

    dd = d[sel].reset_index(drop=True)
    ra_s = res_a[sel]
    rd_s = res_d[sel]
    rhat_s = r_hat[sel]
    ea_s = e_alpha[sel]
    ed_s = e_delta[sel]

    # ---- 3D acceleration direction ----------------------------------------
    # The PMa is a transverse velocity change; promote it to a full 3D vector
    # so that stars on different lines of sight can be compared honestly.
    dist_kpc = 1.0 / dd["plx"].to_numpy(float)          # plx in mas -> kpc
    scale = KMS_PER_MASYR_KPC * dist_kpc                 # mas/yr -> km/s
    dv_icrs = (ra_s[:, None] * ea_s + rd_s[:, None] * ed_s) * scale[:, None]

    R_gal = icrs_to_galactic_matrix()
    dv_gal = dv_icrs @ R_gal.T
    pos_icrs = rhat_s * (dist_kpc * 1000.0)[:, None]     # pc
    xyz = pos_icrs @ R_gal.T

    norm = np.linalg.norm(dv_gal, axis=1)
    good = norm > 0
    unit = dv_gal[good] / norm[good, None]
    xyz = xyz[good]
    dd = dd[good].reset_index(drop=True)
    log.info("directional sample: %d stars", len(unit))

    dv_kms = norm[good]
    log.info("velocity-anomaly magnitude: median %.3f km/s, 90th pct %.3f",
             float(np.median(dv_kms)), float(np.percentile(dv_kms, 90)))

    # ---- coherence vs shuffled null ---------------------------------------
    k = min(args.k, len(unit) - 1)
    r_obs, mean_obs = coherence_field(xyz, unit, k)
    log.info("observed mean local resultant (k=%d): %.5f", k, mean_obs)

    null_means = np.empty(args.n_shuffle)
    null_maxes = np.empty(args.n_shuffle)
    for i in range(args.n_shuffle):
        r_n, m_n = coherence_field(xyz, unit, k, rng=rng, shuffle=True)
        null_means[i] = m_n
        null_maxes[i] = r_n.max()
        if (i + 1) % 50 == 0:
            log.info("  null %d/%d", i + 1, args.n_shuffle)

    mu_n, sd_n = float(null_means.mean()), float(null_means.std())
    z_mean = (mean_obs - mu_n) / sd_n if sd_n > 0 else 0.0

    max_obs = float(r_obs.max())
    mu_mx, sd_mx = float(null_maxes.mean()), float(null_maxes.std())
    z_max = (max_obs - mu_mx) / sd_mx if sd_mx > 0 else 0.0

    log.info("null mean resultant: %.5f +/- %.5f  ->  z = %+.2f",
             mu_n, sd_n, z_mean)
    log.info("observed max local resultant: %.4f", max_obs)
    log.info("null max: %.4f +/- %.4f  ->  z = %+.2f", mu_mx, sd_mx, z_max)

    # ---- where is the most coherent patch, and what is in it? -------------
    j = int(np.argmax(r_obs))
    patch = {
        "resultant": max_obs,
        "l_deg": float(np.degrees(np.arctan2(xyz[j, 1], xyz[j, 0])) % 360),
        "b_deg": float(np.degrees(np.arcsin(
            xyz[j, 2] / max(np.linalg.norm(xyz[j]), 1e-9)))),
        "dist_pc": float(np.linalg.norm(xyz[j])),
    }
    log.info("most coherent neighbourhood: l=%.1f b=%.1f d=%.0f pc, R=%.3f",
             patch["l_deg"], patch["b_deg"], patch["dist_pc"], max_obs)

    # ---- injection-recovery: what coherent fraction would we have seen? ---
    log.info("injection-recovery ...")
    inj_rows = []
    for frac in (0.02, 0.05, 0.10, 0.20, 0.40):
        # Coherently align a random localised patch of this fraction.
        n_inj = max(int(frac * len(unit)), 2)
        from scipy.spatial import cKDTree
        tree = cKDTree(xyz)
        centre = xyz[rng.integers(len(xyz))]
        _, patch_idx = tree.query(centre, k=n_inj)

        v_inj = unit.copy()
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        v_inj[patch_idx] = direction

        _, m_inj = coherence_field(xyz, v_inj, k)
        z_inj = (m_inj - mu_n) / sd_n if sd_n > 0 else 0.0
        inj_rows.append({"coherent_fraction": frac, "n_injected": n_inj,
                         "mean_resultant": m_inj, "z": z_inj})
        log.info("  f=%.2f (n=%d): mean R = %.5f, z = %+.1f",
                 frac, n_inj, m_inj, z_inj)

    inj = pd.DataFrame(inj_rows)
    detectable = inj[inj["z"] > 5.0]
    f_min = float(detectable["coherent_fraction"].min()) if len(detectable) else None

    # ---- verdict -----------------------------------------------------------
    if z_mean > 5.0:
        verdict = (f"COHERENT ACCELERATION FIELD detected at {z_mean:.1f} sigma "
                   f"above the shuffled null. Investigate: check for residual "
                   f"frame systematics and for a single dominant cluster.")
    elif z_max > 5.0:
        verdict = (f"No global coherence ({z_mean:+.1f} sigma) but a localised "
                   f"patch reaches {z_max:.1f} sigma above the null maximum. "
                   f"Most likely a residual cluster or a noise excursion; "
                   f"inspect the patch membership.")
    else:
        verdict = (f"NULL. Acceleration directions are isotropic: mean local "
                   f"resultant {z_mean:+.2f} sigma, max patch {z_max:+.2f} "
                   f"sigma from the shuffled null. Consistent with unseen "
                   f"binaries pointing at random, which is the expected "
                   f"background. No coordinated stellar migration in this "
                   f"volume at the injected sensitivity.")

    print(f"\n{'='*64}")
    print("SEARCH F: COHERENT ACCELERATION DOMAINS")
    print(f"{'='*64}")
    print(f"  stars within 500 pc          : {len(d):,}")
    print(f"  PMa significant (>{args.snr_min:.0f} sigma)   : {len(unit):,}")
    print(f"  frame rotation removed       : "
          f"|omega| = {np.linalg.norm(params[:3]):.4f} mas/yr")
    print(f"  mean local resultant         : {mean_obs:.5f}")
    print(f"  shuffled null                : {mu_n:.5f} +/- {sd_n:.5f}")
    print(f"  excess                       : {z_mean:+.2f} sigma")
    print(f"  max local patch              : {max_obs:.4f}  ({z_max:+.2f} sigma)")
    if f_min is not None:
        print(f"  detectable coherent fraction : >= {f_min:.0%}")
    else:
        print(f"  detectable coherent fraction : none below 40%")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "n_hgca": int(len(df)),
        "n_within_500pc": int(len(d)),
        "n_significant_pma": int(len(unit)),
        "snr_min": args.snr_min,
        "k_neighbours": k,
        "frame_rotation_masyr": [float(v) for v in params[:3]],
        "frame_glide_masyr": [float(v) for v in params[3:]],
        "median_dv_kms": float(np.median(dv_kms)),
        "mean_resultant_observed": mean_obs,
        "mean_resultant_null": mu_n,
        "mean_resultant_null_std": sd_n,
        "z_mean": float(z_mean),
        "max_resultant_observed": max_obs,
        "max_resultant_null": mu_mx,
        "z_max": float(z_max),
        "most_coherent_patch": patch,
        "injection": inj_rows,
        "min_detectable_coherent_fraction": f_min,
        "n_shuffle": args.n_shuffle,
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchF_coherent_accel_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)

    dd = dd.assign(local_resultant=r_obs[good] if len(r_obs) == len(good) else r_obs,
                   dv_kms=dv_kms)
    top = dd.nlargest(200, "local_resultant")
    csv = cfg.RESULT_DIR / f"searchF_coherent_patches_{args.tag}.csv"
    top.to_csv(csv, index=False)
    log.info("wrote %s", csv)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
