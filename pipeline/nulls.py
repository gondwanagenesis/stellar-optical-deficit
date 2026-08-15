"""Step 3: measure the systematic floor.

The point of this module
------------------------
The naive sensitivity sigma/sqrt(N) assumes the residuals are independent and
unbiased.  They are not: extinction errors, photometric calibration, and the
selection function all vary coherently across the sky and across magnitude, so
averaging more stars stops helping once the coherent component dominates.

A subtlety that is easy to get wrong: if you draw *random* subsamples from a
single catalogue and look at the scatter of their means, you will recover
sigma/sqrt(N) essentially perfectly all the way to N = N_total, because random
subsampling cannot manufacture a systematic.  That test is a control, not a
measurement of the floor.

The floor is measured with **structured** subsamples -- groups of stars that
share a systematic.  Here that means contiguous sky patches (HEALPix pixels at
a range of nside), and bins in extinction, apparent magnitude, distance and
crowding.  The scatter of group means as a function of group size N follows
sigma/sqrt(N) while statistics dominate and flattens at the coherent-systematic
level.  That plateau is the reportable sensitivity.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import config as cfg

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Two-sided split tests
# --------------------------------------------------------------------------

def _mean_and_err(r: np.ndarray) -> tuple[float, float, int]:
    r = r[np.isfinite(r)]
    n = len(r)
    if n == 0:
        return np.nan, np.nan, 0
    return float(np.mean(r)), float(np.std(r, ddof=1) / np.sqrt(n)), n


def split_test(resid: np.ndarray, mask_a: np.ndarray, mask_b: np.ndarray,
               label: str) -> dict:
    """Mean residual in two disjoint halves that must agree if unbiased."""
    ma, ea, na = _mean_and_err(resid[mask_a])
    mb, eb, nb = _mean_and_err(resid[mask_b])
    diff = ma - mb
    err = np.hypot(ea, eb)
    return {
        "test": label,
        "mean_a": ma, "err_a": ea, "n_a": na,
        "mean_b": mb, "err_b": eb, "n_b": nb,
        "difference": diff,
        "difference_err": err,
        "n_sigma": diff / err if err > 0 else np.nan,
    }


def standard_splits(df: pd.DataFrame, resid: np.ndarray) -> pd.DataFrame:
    """The battery of splits that must all return zero if the pipeline is clean."""
    rows = []

    b = df["b"].to_numpy(dtype=float)
    rows.append(split_test(resid, b > 0, b < 0, "galactic hemisphere (N vs S)"))

    a0 = df["A_0"].to_numpy(dtype=float)
    q1, q3 = np.nanquantile(a0, [0.25, 0.75])
    rows.append(split_test(resid, a0 <= q1, a0 >= q3,
                           "extinction quartile (low vs high A_0)"))

    g = df["phot_g_mean_mag"].to_numpy(dtype=float)
    gmed = np.nanmedian(g)
    rows.append(split_test(resid, g <= gmed, g > gmed, "apparent G (bright vs faint)"))

    d = df["dist_pc"].to_numpy(dtype=float)
    dmed = np.nanmedian(d)
    rows.append(split_test(resid, d <= dmed, d > dmed, "distance (near vs far)"))

    rho = df["sky_density"].to_numpy(dtype=float)
    rmed = np.nanmedian(rho)
    rows.append(split_test(resid, rho <= rmed, rho > rmed,
                           "crowding (sparse vs crowded)"))

    if "mh_gspphot" in df:
        z = df["mh_gspphot"].to_numpy(dtype=float)
        zmed = np.nanmedian(z)
        rows.append(split_test(resid, z <= zmed, z > zmed,
                               "metallicity (metal-poor vs metal-rich)"))

    absb = np.abs(b)
    rows.append(split_test(resid, absb < 20, absb >= 20,
                           "galactic latitude (|b|<20 vs |b|>=20)"))

    colour = df["bp_rp0"].to_numpy(dtype=float)
    cmed = np.nanmedian(colour)
    rows.append(split_test(resid, colour <= cmed, colour > cmed,
                           "colour (blue vs red half)"))

    # --- conditional splits ------------------------------------------------
    # The colour split is usually the worst offender, and the obvious suspect
    # is the colour-dependent extinction coefficient.  Repeating it inside the
    # lowest-extinction quartile separates the two hypotheses: if the bias
    # survives where there is almost no dust, it is not an extinction-law
    # error but a genuine inadequacy of the M_Ks + [M/H] parametrisation.
    lowdust = a0 <= q1
    if lowdust.sum() > 1000:
        c_lo = np.nanmedian(colour[lowdust])
        rows.append(split_test(
            resid, lowdust & (colour <= c_lo), lowdust & (colour > c_lo),
            "colour split WITHIN lowest-extinction quartile"))

    # Same idea for distance: a distance-dependent bias at low extinction
    # points at the parallax zero point or the selection function rather than
    # at the dust map.
    if lowdust.sum() > 1000:
        d_lo = np.nanmedian(d[lowdust])
        rows.append(split_test(
            resid, lowdust & (d <= d_lo), lowdust & (d > d_lo),
            "distance split WITHIN lowest-extinction quartile"))

    return pd.DataFrame(rows)


def dust_map_paired_test(df: pd.DataFrame, resid_a: np.ndarray,
                         resid_b: np.ndarray, label_a: str,
                         label_b: str) -> dict:
    """Paired difference between two extinction treatments on the same stars.

    Paired, so the intrinsic scatter cancels and what is left is purely the
    systematic difference between the two treatments.
    """
    d = resid_a - resid_b
    d = d[np.isfinite(d)]
    n = len(d)
    return {
        "test": f"paired: {label_a} - {label_b}",
        "mean_difference": float(np.mean(d)),
        "std_difference": float(np.std(d, ddof=1)),
        "err_of_mean": float(np.std(d, ddof=1) / np.sqrt(n)),
        "n": n,
        "rms": float(np.sqrt(np.mean(d ** 2))),
    }


# --------------------------------------------------------------------------
# Mean-residual scatter versus group size
# --------------------------------------------------------------------------

# An RMS estimated from a handful of groups is itself so noisy that it can
# fake a plateau ending.  Below this many surviving groups the point is dropped
# rather than plotted, because the alternative is a downturn at large N that
# looks like the systematic going away.
MIN_GROUPS_FOR_RMS = 8


def _group_mean_scatter(resid: np.ndarray, group_id: np.ndarray,
                        min_group: int = 20) -> tuple[float, float, int, int]:
    """Scatter of per-group mean residuals.

    Returns (rms_of_group_means, expected_rms_if_pure_noise, n_groups, mean_n).
    """
    ok = np.isfinite(resid)
    resid = resid[ok]
    group_id = group_id[ok]

    order = np.argsort(group_id, kind="stable")
    r = resid[order]
    g = group_id[order]
    edges = np.flatnonzero(np.diff(g)) + 1
    starts = np.concatenate([[0], edges])
    ends = np.concatenate([edges, [len(g)]])
    counts = ends - starts
    keep = counts >= min_group
    if keep.sum() < MIN_GROUPS_FOR_RMS:
        return np.nan, np.nan, int(keep.sum()), 0

    sums = np.add.reduceat(r, starts)[keep]
    counts = counts[keep]
    means = sums / counts

    # Expected scatter of group means under pure independent noise.
    var_all = float(np.var(r, ddof=1))
    expected = float(np.sqrt(np.mean(var_all / counts)))
    observed = float(np.sqrt(np.mean((means - np.mean(means)) ** 2)))
    return observed, expected, int(len(means)), int(np.mean(counts))


def spatial_scaling(df: pd.DataFrame, resid: np.ndarray,
                    nsides=(64, 32, 16, 8, 4, 2, 1)) -> pd.DataFrame:
    """Scatter of sky-patch mean residuals vs patch size.

    Statistical noise falls as 1/sqrt(N).  A coherent spatial systematic does
    not, so the observed curve peels away from the expected curve and flattens.
    The plateau is the spatial systematic floor.
    """
    import healpy as hp
    theta = np.radians(90.0 - df["b"].to_numpy(dtype=float))
    phi = np.radians(df["l"].to_numpy(dtype=float))
    rows = []
    for nside in nsides:
        pix = hp.ang2pix(nside, theta, phi, nest=True)
        obs, exp, ngrp, meann = _group_mean_scatter(resid, pix)
        rows.append({"axis": "sky (HEALPix)", "nside": nside,
                     "n_groups": ngrp, "mean_group_n": meann,
                     "rms_group_means": obs, "expected_if_noise": exp,
                     "excess": np.sqrt(max(obs ** 2 - exp ** 2, 0.0))
                     if np.isfinite(obs) and np.isfinite(exp) else np.nan})
        log.info("nside=%-4d n_groups=%5d  <n>=%7d  rms=%.5f  expected=%.5f",
                 nside, ngrp, meann, obs, exp)
    return pd.DataFrame(rows)


def binned_scaling(df: pd.DataFrame, resid: np.ndarray, column: str,
                   n_bins_grid=(2000, 1000, 400, 200, 100, 40, 20, 10, 5)) -> pd.DataFrame:
    """Same diagnostic along a non-spatial axis (extinction, magnitude, ...)."""
    x = df[column].to_numpy(dtype=float)
    ok = np.isfinite(x)
    rows = []
    for nb in n_bins_grid:
        q = np.linspace(0, 1, nb + 1)
        edges = np.nanquantile(x[ok], q)
        edges = np.unique(edges)
        if len(edges) < 3:
            continue
        gid = np.digitize(x, edges[1:-1])
        obs, exp, ngrp, meann = _group_mean_scatter(resid, gid)
        rows.append({"axis": column, "n_bins": nb, "n_groups": ngrp,
                     "mean_group_n": meann, "rms_group_means": obs,
                     "expected_if_noise": exp,
                     "excess": np.sqrt(max(obs ** 2 - exp ** 2, 0.0))
                     if np.isfinite(obs) and np.isfinite(exp) else np.nan})
    return pd.DataFrame(rows)


def random_subsample_scaling(resid: np.ndarray, n_grid=None,
                             n_realisations: int | None = None,
                             seed: int | None = None) -> pd.DataFrame:
    """Control: random subsamples, which cannot manufacture a systematic.

    This is expected to track sigma/sqrt(N) exactly.  It is reported so that
    the contrast with the structured tests is explicit, and so that anyone
    quoting sigma/sqrt(N) can see why that number is not a sensitivity.
    """
    n_grid = n_grid or cfg.NULLS.n_grid
    n_realisations = n_realisations or cfg.NULLS.n_realisations
    seed = seed if seed is not None else cfg.NULLS.rng_seed

    rng = np.random.default_rng(seed)
    r = resid[np.isfinite(resid)]
    sigma = float(np.std(r, ddof=1))
    n_tot = len(r)
    rows = []
    for n in n_grid:
        if n > n_tot:
            continue
        means = np.empty(n_realisations)
        for i in range(n_realisations):
            idx = rng.integers(0, n_tot, size=n)
            means[i] = r[idx].mean()
        rows.append({"axis": "random subsample", "n": n,
                     "rms_of_means": float(np.std(means, ddof=1)),
                     "sigma_over_sqrt_n": sigma / np.sqrt(n)})
    return pd.DataFrame(rows)
