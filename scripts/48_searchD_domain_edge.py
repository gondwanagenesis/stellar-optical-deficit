#!/usr/bin/env python
"""Search D: domain-edge discontinuity in the 3D residual field.

    run.sh scripts/48_searchD_domain_edge.py --tag primary

THE HYPOTHESIS
--------------
The "grabby aliens" model (Hanson, Martin & McCarter 2021) predicts expanding
domains -- civilisations converting stars at near-c expansion speed.  A domain
boundary is a SURFACE in 3D space where the stellar population transitions from
normal to harvested.  Stars on the harvested side are systematically dimmer in
the optical (positive residual); stars on the normal side are not.

This script searches for that boundary by looking for DISCONTINUITIES in the
3D distribution of optical residuals, using three complementary methods:

  1. PLANAR SCAN -- fit oriented half-spaces at many orientations and offsets.
     For each plane, compare the mean residual on the two sides.  A domain edge
     intersecting our volume produces a plane orientation where the two sides
     differ by much more than the shuffled null.

  2. RADIAL SCAN -- test whether the mean residual changes discontinuously at
     some distance from the Sun.  An expanding domain centred elsewhere
     produces a transition radius if the boundary crosses our volume.

  3. LOCAL GRADIENT -- for each star, compare its local neighbourhood residual
     to a surrounding shell.  A domain boundary produces a coherent gradient
     vector field; random noise does not.  The ALIGNMENT of the gradient vectors
     is the discriminant.

THE CAVEAT
----------
Our 500 pc volume is tiny on cosmological scales.  A domain expanding at
near-c for 1 Gyr is ~300 Mpc across, so if the boundary passed through our
neighbourhood it would appear FLAT to us.  The planar scan is therefore the
geometrically correct search for a cosmological-speed domain.  The radial scan
catches a domain centred on a nearby star cluster that started expanding more
recently (within ~1000 yr for a 500 pc boundary).

This is the weakest prior of all our searches, but it completes the channel
space: if the boundary is here, one of these three scans will find it.

CONTROLS
--------
- Shuffle residuals to measure the null distribution of every statistic.
- Check whether any detected edge correlates with Galactic latitude or
  extinction, which would identify it as dust rather than engineering.
- Restrict to clean, low-extinction stars to avoid dust-driven edges.
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

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

# Planar scan: number of plane orientations (HEALPix-like icosphere sampling)
N_DIRECTIONS = 192            # 192 directions = HEALPix nside=4
N_OFFSETS_PER_DIR = 20        # plane offsets along each normal
MIN_STARS_PER_SIDE = 500      # require enough stars on each side of the plane

# Radial scan
RADIAL_BINS = 25              # number of radial shells
RADIAL_ANNULUS_N = 5          # sliding window half-width in bins

# Local gradient
K_INNER = 30                  # neighbours for the inner sphere
SHELL_FACTOR = 2.5            # shell radius = SHELL_FACTOR * inner radius
K_SHELL = 60                  # neighbours to sample in the shell
GRADIENT_SUBSAMPLE = 50_000   # subsample for gradient computation (memory)

# Null calibration
N_SHUFFLE = 50

# Quality cuts for this search
A0_MAX = 0.3
CSTAR_NSIGMA_MAX = 3.0
RUWE_MAX = 1.4


def galactic_to_xyz(l_deg: np.ndarray, b_deg: np.ndarray,
                    dist_pc: np.ndarray) -> np.ndarray:
    """Galactic (l, b, dist) to Cartesian (X, Y, Z) with Sun at origin."""
    l_rad = np.radians(l_deg)
    b_rad = np.radians(b_deg)
    cb = np.cos(b_rad)
    return np.column_stack([dist_pc * cb * np.cos(l_rad),
                            dist_pc * cb * np.sin(l_rad),
                            dist_pc * np.sin(b_rad)])


def fibonacci_sphere(n: int) -> np.ndarray:
    """Generate n approximately uniformly spaced unit vectors on the sphere.

    Uses the Fibonacci-lattice method (Gonzalez 2010).  This is faster and
    simpler than building a HEALPix grid and gives adequate angular uniformity
    for the plane-scan.
    """
    golden = (1 + np.sqrt(5)) / 2
    i = np.arange(n)
    theta = 2 * np.pi * i / golden
    phi = np.arccos(1 - 2 * (i + 0.5) / n)
    return np.column_stack([np.sin(phi) * np.cos(theta),
                            np.sin(phi) * np.sin(theta),
                            np.cos(phi)])


# --------------------------------------------------------------------------
# Method 1: planar scan
# --------------------------------------------------------------------------

def _scan_directions(projections: np.ndarray, resid: np.ndarray,
                     sigma: float, dir_indices: np.ndarray) -> float:
    """Vectorised scan over directions: sort once, cumsum to sweep offsets.

    Returns the maximum |z| found across all directions and offsets.
    """
    n_stars = len(resid)
    total_sum = resid.sum()
    best_abs_z = 0.0

    for d_idx in dir_indices:
        proj = projections[:, d_idx]
        order = np.argsort(proj)
        r_sorted = resid[order]
        p_sorted = proj[order]
        cumsum = np.cumsum(r_sorted)

        lo_idx = int(n_stars * 0.05)
        hi_idx = int(n_stars * 0.95)
        split_points = np.linspace(lo_idx, hi_idx, N_OFFSETS_PER_DIR,
                                   dtype=int)
        split_points = np.clip(split_points, MIN_STARS_PER_SIDE,
                               n_stars - MIN_STARS_PER_SIDE)
        split_points = np.unique(split_points)

        left_sum = cumsum[split_points - 1]
        left_n = split_points.astype(float)
        right_n = float(n_stars) - left_n
        right_sum = total_sum - left_sum

        left_mean = left_sum / left_n
        right_mean = right_sum / right_n
        delta = right_mean - left_mean
        se = sigma * np.sqrt(1.0 / left_n + 1.0 / right_n)
        z = np.abs(delta / se)
        z_max = float(z.max())
        if z_max > best_abs_z:
            best_abs_z = z_max

    return best_abs_z


def planar_scan(xyz: np.ndarray, resid: np.ndarray,
                rng: np.random.Generator) -> dict:
    """Test every plane orientation for a mean-residual discontinuity.

    For each normal direction n, project all stars onto n and slide a dividing
    plane along the projected axis.  At each offset, compute the difference in
    mean residual between the two half-spaces.  The maximum difference over all
    orientations and offsets is the test statistic.

    Vectorised: sort once per direction, cumulative-sum sweep over offsets.
    The null is calibrated by shuffling the residuals.
    """
    normals = fibonacci_sphere(N_DIRECTIONS)
    sigma = float(np.std(resid, ddof=1))
    n_stars = len(resid)
    total_sum = float(resid.sum())

    projections = xyz @ normals.T  # (n_stars, N_DIRECTIONS)

    best_delta = 0.0
    best_z = 0.0
    best_dir_idx = -1
    best_offset_idx = 0

    all_dirs = np.arange(N_DIRECTIONS)
    for d_idx in all_dirs:
        proj = projections[:, d_idx]
        order = np.argsort(proj)
        r_sorted = resid[order]
        cumsum = np.cumsum(r_sorted)

        lo_idx = int(n_stars * 0.05)
        hi_idx = int(n_stars * 0.95)
        split_points = np.linspace(lo_idx, hi_idx, N_OFFSETS_PER_DIR,
                                   dtype=int)
        split_points = np.clip(split_points, MIN_STARS_PER_SIDE,
                               n_stars - MIN_STARS_PER_SIDE)
        split_points = np.unique(split_points)

        left_sum = cumsum[split_points - 1]
        left_n = split_points.astype(float)
        right_n = float(n_stars) - left_n
        right_sum = total_sum - left_sum

        left_mean = left_sum / left_n
        right_mean = right_sum / right_n
        delta = right_mean - left_mean
        se = sigma * np.sqrt(1.0 / left_n + 1.0 / right_n)
        z = delta / se
        k = int(np.argmax(np.abs(z)))
        if abs(z[k]) > abs(best_z):
            best_z = float(z[k])
            best_delta = float(delta[k])
            best_dir_idx = d_idx
            best_offset_idx = int(split_points[k])

    best_offset_frac = best_offset_idx / max(n_stars, 1)

    # --- null calibration: shuffle residuals N_SHUFFLE times ----------------
    null_dirs = np.arange(0, N_DIRECTIONS, 4)  # every 4th for speed
    null_max_z = np.empty(N_SHUFFLE)
    for s in range(N_SHUFFLE):
        resid_s = rng.permutation(resid)
        null_max_z[s] = _scan_directions(projections, resid_s, sigma,
                                         null_dirs)
        if (s + 1) % 10 == 0:
            print(f"    null shuffle {s + 1}/{N_SHUFFLE}", flush=True)

    trials_correction = np.sqrt(
        np.log(N_DIRECTIONS * N_OFFSETS_PER_DIR)
        / max(np.log(len(null_dirs) * N_OFFSETS_PER_DIR), 1.0))
    null_max_z_corrected = null_max_z * trials_correction

    excess_sigma = ((abs(best_z) - np.mean(null_max_z_corrected))
                    / max(np.std(null_max_z_corrected, ddof=1), 1e-10))

    best_normal = normals[best_dir_idx] if best_dir_idx >= 0 else [0, 0, 0]
    # Convert normal to Galactic l, b
    x, y, z_n = best_normal
    best_l = float(np.degrees(np.arctan2(y, x)) % 360)
    best_b = float(np.degrees(np.arcsin(np.clip(z_n, -1, 1))))

    return {
        "method": "planar_scan",
        "n_directions": N_DIRECTIONS,
        "n_offsets_per_dir": N_OFFSETS_PER_DIR,
        "best_delta_mag": float(best_delta),
        "best_z": float(best_z),
        "best_normal_lbdeg": [best_l, best_b],
        "best_offset_frac": float(best_offset_frac),
        "null_max_z_mean": float(np.mean(null_max_z_corrected)),
        "null_max_z_std": float(np.std(null_max_z_corrected, ddof=1)),
        "excess_over_null_sigma": float(excess_sigma),
        "n_shuffle": N_SHUFFLE,
    }


# --------------------------------------------------------------------------
# Method 2: radial scan
# --------------------------------------------------------------------------

def radial_scan(dist_pc: np.ndarray, resid: np.ndarray,
                rng: np.random.Generator) -> dict:
    """Test for a discontinuity in mean residual vs heliocentric distance.

    Bin stars radially, compute the running mean, and look for a step -- a
    sudden change in mean residual at some distance.  The maximum absolute
    difference between adjacent running means is the test statistic.
    """
    edges = np.linspace(dist_pc.min(), dist_pc.max(), RADIAL_BINS + 1)
    bin_idx = np.clip(np.digitize(dist_pc, edges) - 1, 0, RADIAL_BINS - 1)
    bin_means = np.array([np.mean(resid[bin_idx == i])
                          if np.sum(bin_idx == i) > 10 else np.nan
                          for i in range(RADIAL_BINS)])
    bin_counts = np.array([np.sum(bin_idx == i) for i in range(RADIAL_BINS)])
    bin_centres = 0.5 * (edges[:-1] + edges[1:])

    # Running mean with window = 2*RADIAL_ANNULUS_N + 1
    valid = np.isfinite(bin_means) & (bin_counts > 30)
    if valid.sum() < 5:
        return {"method": "radial_scan", "verdict": "insufficient data"}

    # Step statistic: max |difference| between the mean of bins [0..k] and
    # the mean of bins [k+1..end], weighted by count
    best_step = 0.0
    best_step_z = 0.0
    best_k = -1
    sigma = float(np.std(resid, ddof=1))
    for k in range(2, RADIAL_BINS - 3):
        left = (bin_idx <= k)
        right = (bin_idx > k)
        n_l, n_r = left.sum(), right.sum()
        if n_l < 100 or n_r < 100:
            continue
        delta = np.mean(resid[right]) - np.mean(resid[left])
        se = sigma * np.sqrt(1.0 / n_l + 1.0 / n_r)
        z = delta / se if se > 0 else 0.0
        if abs(z) > abs(best_step_z):
            best_step_z = z
            best_step = delta
            best_k = k

    # Null calibration
    null_max_z = []
    for _ in range(N_SHUFFLE):
        rs = rng.permutation(resid)
        null_best = 0.0
        for k in range(2, RADIAL_BINS - 3):
            left = (bin_idx <= k)
            right = (bin_idx > k)
            n_l, n_r = left.sum(), right.sum()
            if n_l < 100 or n_r < 100:
                continue
            delta = np.mean(rs[right]) - np.mean(rs[left])
            se = sigma * np.sqrt(1.0 / n_l + 1.0 / n_r)
            z_null = abs(delta / se) if se > 0 else 0.0
            null_best = max(null_best, z_null)
        null_max_z.append(null_best)
    null_max_z = np.array(null_max_z)

    excess_sigma = ((abs(best_step_z) - np.mean(null_max_z))
                    / max(np.std(null_max_z, ddof=1), 1e-10))

    transition_pc = float(edges[best_k + 1]) if best_k >= 0 else np.nan

    return {
        "method": "radial_scan",
        "n_bins": RADIAL_BINS,
        "best_step_mag": float(best_step),
        "best_step_z": float(best_step_z),
        "transition_dist_pc": transition_pc,
        "null_max_z_mean": float(np.mean(null_max_z)),
        "null_max_z_std": float(np.std(null_max_z, ddof=1)),
        "excess_over_null_sigma": float(excess_sigma),
        "n_shuffle": N_SHUFFLE,
        "bin_centres_pc": [float(x) for x in bin_centres[valid]],
        "bin_means_mag": [float(x) for x in bin_means[valid]],
        "bin_counts": [int(x) for x in bin_counts[valid]],
    }


# --------------------------------------------------------------------------
# Method 3: local gradient
# --------------------------------------------------------------------------

def local_gradient(xyz: np.ndarray, resid: np.ndarray,
                   rng: np.random.Generator) -> dict:
    """Detect coherent gradient in the residual field.

    For a subsample of stars, compute the vector pointing from "inner
    neighbourhood mean residual" to "outer shell mean residual".  If a domain
    boundary exists, these gradient vectors are ALIGNED -- they all point
    from the harvested side to the normal side.  Random noise produces
    randomly oriented gradients whose mean direction is zero.

    The test statistic is the RESULTANT LENGTH of the unit gradient vectors,
    compared against a shuffled null.
    """
    n = len(resid)
    if n > GRADIENT_SUBSAMPLE:
        idx = rng.choice(n, GRADIENT_SUBSAMPLE, replace=False)
    else:
        idx = np.arange(n)

    xyz_sub = xyz[idx]
    # Build tree on ALL stars for querying
    tree = cKDTree(xyz)

    # For each subsample star, query K_INNER nearest neighbours to define
    # the inner sphere, then query K_INNER + K_SHELL neighbours to get the
    # shell stars
    k_total = K_INNER + K_SHELL
    dists, indices = tree.query(xyz_sub, k=min(k_total + 1, n))
    # indices[:,0] is the point itself (distance 0); skip it
    inner_idx = indices[:, 1:K_INNER + 1]
    outer_idx = indices[:, K_INNER + 1:k_total + 1]

    # Inner mean residual
    inner_mean = np.mean(resid[inner_idx], axis=1)

    # Outer mean residual
    outer_mean = np.mean(resid[outer_idx], axis=1)

    # Gradient direction: vector from this star toward the mean position of
    # the outer shell, weighted by (outer_mean - inner_mean)
    # Actually: the gradient vector is the mean displacement of outer stars
    # minus the mean displacement of inner stars, weighted by the residual
    # difference.  Simpler: use the centroid of the outer shell relative to
    # the star as the gradient direction.
    outer_centroids = np.mean(xyz[outer_idx], axis=1) - xyz_sub
    norms = np.linalg.norm(outer_centroids, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    directions = outer_centroids / norms

    # Sign: positive gradient magnitude = outer is dimmer than inner
    gradient_mag = outer_mean - inner_mean

    # Weighted unit vectors: sign(gradient) * direction
    signed_dirs = directions * np.sign(gradient_mag)[:, None]

    # Resultant length of unit vectors -- tests for alignment
    resultant = np.linalg.norm(np.mean(signed_dirs, axis=0))
    # Under random orientations, E[R] = 1/sqrt(N), so the z-score is:
    n_sub = len(idx)
    expected_R = 1.0 / np.sqrt(n_sub)
    rayleigh_z = (resultant - expected_R) / expected_R  # rough z-score

    # Also compute: mean absolute gradient magnitude and the fraction of
    # stars where outer > inner (would be >50% if a boundary makes the far
    # side dimmer)
    frac_outer_dimmer = float(np.mean(gradient_mag > 0))

    # Mean gradient direction in Galactic coordinates
    mean_dir = np.mean(signed_dirs, axis=0)
    mean_dir_norm = mean_dir / max(np.linalg.norm(mean_dir), 1e-10)
    grad_l = float(np.degrees(np.arctan2(mean_dir_norm[1], mean_dir_norm[0])) % 360)
    grad_b = float(np.degrees(np.arcsin(np.clip(mean_dir_norm[2], -1, 1))))

    # --- null calibration: shuffle residuals --------------------------------
    null_resultants = []
    for _ in range(N_SHUFFLE):
        rs = rng.permutation(resid)
        inner_m = np.mean(rs[inner_idx], axis=1)
        outer_m = np.mean(rs[outer_idx], axis=1)
        gm = outer_m - inner_m
        sd = directions * np.sign(gm)[:, None]
        null_resultants.append(np.linalg.norm(np.mean(sd, axis=0)))
    null_resultants = np.array(null_resultants)

    excess_sigma = ((resultant - np.mean(null_resultants))
                    / max(np.std(null_resultants, ddof=1), 1e-10))

    return {
        "method": "local_gradient",
        "n_subsample": n_sub,
        "K_inner": K_INNER,
        "K_shell": K_SHELL,
        "resultant_length": float(resultant),
        "null_resultant_mean": float(np.mean(null_resultants)),
        "null_resultant_std": float(np.std(null_resultants, ddof=1)),
        "excess_over_null_sigma": float(excess_sigma),
        "frac_outer_dimmer": frac_outer_dimmer,
        "mean_gradient_direction_lbdeg": [grad_l, grad_b],
        "mean_gradient_magnitude_mag": float(np.mean(np.abs(gradient_mag))),
        "n_shuffle": N_SHUFFLE,
    }


# --------------------------------------------------------------------------
# Dust / structure cross-checks
# --------------------------------------------------------------------------

def dust_correlation_check(xyz: np.ndarray, resid: np.ndarray,
                           b: np.ndarray, A0: np.ndarray,
                           best_normal: list, best_offset_frac: float) -> dict:
    """Check whether the best planar edge aligns with the dust plane."""
    normal = np.array(best_normal)
    norm = np.linalg.norm(normal)
    if norm < 1e-10:
        return {"dust_check": "no edge found"}
    normal = normal / norm

    # Angle between the best-fit plane normal and the Galactic Z axis
    z_axis = np.array([0.0, 0.0, 1.0])
    cos_angle = abs(float(np.dot(normal, z_axis)))
    angle_from_plane_deg = float(np.degrees(np.arcsin(np.clip(cos_angle, 0, 1))))

    # Residual correlation with |b| and A_0
    proj = xyz @ normal
    median_proj = np.median(proj)
    side_A = proj < median_proj
    side_B = ~side_A

    mean_absb_A = float(np.mean(np.abs(b[side_A])))
    mean_absb_B = float(np.mean(np.abs(b[side_B])))
    mean_A0_A = float(np.mean(A0[side_A]))
    mean_A0_B = float(np.mean(A0[side_B]))

    is_dust = (angle_from_plane_deg < 30.0
               and abs(mean_A0_A - mean_A0_B) > 0.01)

    return {
        "plane_normal_angle_from_gal_pole_deg": float(90.0 - angle_from_plane_deg),
        "plane_normal_angle_from_gal_plane_deg": angle_from_plane_deg,
        "mean_absb_side_A_deg": mean_absb_A,
        "mean_absb_side_B_deg": mean_absb_B,
        "mean_A0_side_A": mean_A0_A,
        "mean_A0_side_B": mean_A0_B,
        "consistent_with_dust_plane": is_dust,
    }


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Search D: domain-edge discontinuity in the 3D residual field")
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()
    rng = np.random.default_rng(20260818)

    # ------------------------------------------------------------------
    # Load and clean
    # ------------------------------------------------------------------
    cols = ["source_id", "l", "b", "dist_pc", "residual", "A_0",
            "cstar_nsigma", "ruwe"]
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=cols)
    n_raw = len(d)
    d = d.dropna(subset=["l", "b", "dist_pc", "residual"]).reset_index(drop=True)

    # Quality cuts: low extinction, clean photometry, good astrometry
    mask = ((d["A_0"].fillna(9.0) < A0_MAX)
            & (d["cstar_nsigma"].fillna(9.0).abs() < CSTAR_NSIGMA_MAX)
            & (d["ruwe"].fillna(9.0) < RUWE_MAX))
    d = d[mask].reset_index(drop=True)
    print(f"loaded {n_raw:,} stars -> {len(d):,} after quality cuts "
          f"(A_0<{A0_MAX}, |c*|<{CSTAR_NSIGMA_MAX}, ruwe<{RUWE_MAX})")

    resid = d["residual"].to_numpy(float)
    l_deg = d["l"].to_numpy(float)
    b_deg = d["b"].to_numpy(float)
    dist_pc = d["dist_pc"].to_numpy(float)
    A0 = d["A_0"].to_numpy(float)
    sigma = st.robust_sigma(resid)
    print(f"residual robust sigma = {sigma:.5f} mag, mean = {np.mean(resid):.6f}")

    xyz = galactic_to_xyz(l_deg, b_deg, dist_pc)

    # ------------------------------------------------------------------
    # Method 1: planar scan
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"METHOD 1: PLANAR SCAN ({N_DIRECTIONS} orientations x {N_OFFSETS_PER_DIR} offsets)")
    print(f"{'='*60}")
    planar = planar_scan(xyz, resid, rng)
    print(f"  best delta          = {planar['best_delta_mag']:+.6f} mag")
    print(f"  best z              = {planar['best_z']:+.2f}")
    print(f"  best plane normal   = l={planar['best_normal_lbdeg'][0]:.1f}, "
          f"b={planar['best_normal_lbdeg'][1]:.1f}")
    print(f"  null max z          = {planar['null_max_z_mean']:.2f} "
          f"+/- {planar['null_max_z_std']:.2f}")
    print(f"  excess over null    = {planar['excess_over_null_sigma']:+.2f} sigma")

    # ------------------------------------------------------------------
    # Method 2: radial scan
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"METHOD 2: RADIAL SCAN ({RADIAL_BINS} distance bins)")
    print(f"{'='*60}")
    radial = radial_scan(dist_pc, resid, rng)
    if "best_step_z" in radial:
        print(f"  best step           = {radial['best_step_mag']:+.6f} mag")
        print(f"  best step z         = {radial['best_step_z']:+.2f}")
        print(f"  transition distance = {radial['transition_dist_pc']:.0f} pc")
        print(f"  null max z          = {radial['null_max_z_mean']:.2f} "
              f"+/- {radial['null_max_z_std']:.2f}")
        print(f"  excess over null    = {radial['excess_over_null_sigma']:+.2f} sigma")
    else:
        print(f"  {radial.get('verdict', 'no result')}")

    # ------------------------------------------------------------------
    # Method 3: local gradient
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"METHOD 3: LOCAL GRADIENT (K_inner={K_INNER}, K_shell={K_SHELL}, "
          f"n_sub={GRADIENT_SUBSAMPLE})")
    print(f"{'='*60}")
    gradient = local_gradient(xyz, resid, rng)
    print(f"  resultant length    = {gradient['resultant_length']:.6f}")
    print(f"  null resultant      = {gradient['null_resultant_mean']:.6f} "
          f"+/- {gradient['null_resultant_std']:.6f}")
    print(f"  excess over null    = {gradient['excess_over_null_sigma']:+.2f} sigma")
    print(f"  frac outer dimmer   = {gradient['frac_outer_dimmer']:.3f} (0.5 = no edge)")
    print(f"  mean gradient dir   = l={gradient['mean_gradient_direction_lbdeg'][0]:.1f}, "
          f"b={gradient['mean_gradient_direction_lbdeg'][1]:.1f}")

    # ------------------------------------------------------------------
    # Dust cross-check on the planar result
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("DUST CROSS-CHECK")
    print(f"{'='*60}")
    # Convert the best normal from (l,b) back to Cartesian for the check
    best_l_rad = np.radians(planar["best_normal_lbdeg"][0])
    best_b_rad = np.radians(planar["best_normal_lbdeg"][1])
    best_normal_xyz = [float(np.cos(best_b_rad) * np.cos(best_l_rad)),
                       float(np.cos(best_b_rad) * np.sin(best_l_rad)),
                       float(np.sin(best_b_rad))]
    dust = dust_correlation_check(xyz, resid, b_deg, A0, best_normal_xyz,
                                  planar["best_offset_frac"])
    for k, v in dust.items():
        print(f"  {k:42s} {v}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("VERDICT")
    print(f"{'='*60}")
    methods = [planar, radial, gradient]
    max_excess = max(m.get("excess_over_null_sigma", 0.0) for m in methods)
    any_significant = max_excess > 5.0

    if not any_significant:
        verdict = (f"NO DOMAIN EDGE DETECTED. All three methods are consistent "
                   f"with the shuffled null (max excess = {max_excess:+.1f} sigma). "
                   f"Within 500 pc the stellar residual field has no "
                   f"discontinuity surface above the noise floor.")
    elif dust.get("consistent_with_dust_plane", False):
        who = max(methods, key=lambda m: m.get("excess_over_null_sigma", 0.0))
        verdict = (f"Structure detected by {who['method']} at "
                   f"{who.get('excess_over_null_sigma', 0):.1f} sigma, but the "
                   f"plane normal aligns with the Galactic plane and correlates "
                   f"with extinction -- consistent with UNMODELLED DUST, not a "
                   f"domain boundary.")
    else:
        who = max(methods, key=lambda m: m.get("excess_over_null_sigma", 0.0))
        verdict = (f"ANOMALOUS EDGE detected by {who['method']} at "
                   f"{who.get('excess_over_null_sigma', 0):.1f} sigma, NOT "
                   f"aligned with the dust plane -- INVESTIGATE.")
    print(f"\n  {verdict}")

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    result = {
        "tag": args.tag,
        "n_stars_raw": n_raw,
        "n_stars_clean": len(d),
        "quality_cuts": {"A0_max": A0_MAX,
                         "cstar_nsigma_max": CSTAR_NSIGMA_MAX,
                         "ruwe_max": RUWE_MAX},
        "residual_robust_sigma": float(sigma),
        "planar_scan": planar,
        "radial_scan": radial,
        "local_gradient": gradient,
        "dust_check": dust,
        "verdict": verdict,
    }

    out_path = cfg.RESULT_DIR / f"searchD_domain_edge_{args.tag}.json"
    out_path.write_text(json.dumps(result, indent=2, default=float))
    print(f"\nwrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
