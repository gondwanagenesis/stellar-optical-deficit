"""Co-natal wide-pair identification by self-join of the analysis sample.

Extracted from scripts/24 so that several analyses can share one definition.
Criteria follow El-Badry, Rix & Heintz (2021, MNRAS 506, 2269): small projected
separation, consistent parallax, and a proper-motion difference small enough to
be a bound orbit. Chance alignment is measured by a scramble test, never
assumed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

MAX_ANGSEP_ARCSEC = 120.0
MAX_PHYS_SEP_AU = 20_000.0
MIN_PHYS_SEP_AU = 100.0
PARALLAX_CONSISTENCY_NSIGMA = 3.0
MU_ERR_PER_PARALLAX_ERR = 1.2      # crude DR3 scaling; PM errors were not pulled
MU_MARGIN_NSIGMA = 2.0
M_TOT_MSUN = 2.0                   # generous, so the criterion errs toward keeping
CLEAN_SEP_ARCSEC = 10.0            # beyond 2x the 2MASS aperture; see paper 5.5


def orbital_mu_max(sep_au: np.ndarray, dist_pc: np.ndarray) -> np.ndarray:
    """Maximum proper-motion difference for a bound orbit, mas/yr."""
    v_esc = 29.78 * np.sqrt(2.0 * M_TOT_MSUN / np.maximum(sep_au, 1.0))
    return v_esc / (4.74 * np.maximum(dist_pc, 1.0) / 1000.0)


def find_pairs(d: pd.DataFrame, scramble: bool = False,
               seed: int = 7) -> pd.DataFrame:
    """Self-join for co-moving, co-distant, close pairs.

    ``scramble=True`` permutes the astrometry while keeping positions, which
    destroys real pairs and leaves the chance-alignment rate.
    """
    ra = np.radians(d["ra"].to_numpy(float))
    dec = np.radians(d["dec"].to_numpy(float))
    xyz = np.column_stack([np.cos(dec) * np.cos(ra),
                           np.cos(dec) * np.sin(ra),
                           np.sin(dec)])

    plx = d["parallax_corr"].to_numpy(float)
    plx_err = d["parallax_error"].to_numpy(float)
    pmra = d["pmra"].to_numpy(float)
    pmdec = d["pmdec"].to_numpy(float)
    dist = d["dist_pc"].to_numpy(float)

    if scramble:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(d))
        plx, plx_err = plx[perm], plx_err[perm]
        pmra, pmdec, dist = pmra[perm], pmdec[perm], dist[perm]

    r_chord = 2.0 * np.sin(np.radians(MAX_ANGSEP_ARCSEC / 3600.0) / 2.0)
    pairs = cKDTree(xyz).query_pairs(r_chord, output_type="ndarray")
    if len(pairs) == 0:
        return pd.DataFrame()

    i, j = pairs[:, 0], pairs[:, 1]
    dot = np.clip(np.einsum("ij,ij->i", xyz[i], xyz[j]), -1.0, 1.0)
    theta = np.degrees(np.arccos(dot)) * 3600.0
    dmean = 0.5 * (dist[i] + dist[j])
    sep_au = theta * dmean

    dplx = np.abs(plx[i] - plx[j])
    sig_dplx = np.hypot(plx_err[i], plx_err[j])
    dmu = np.hypot(pmra[i] - pmra[j], pmdec[i] - pmdec[j])
    sig_dmu = MU_ERR_PER_PARALLAX_ERR * sig_dplx
    mu_max = orbital_mu_max(sep_au, dmean)

    keep = ((sep_au > MIN_PHYS_SEP_AU) & (sep_au < MAX_PHYS_SEP_AU)
            & (dplx < PARALLAX_CONSISTENCY_NSIGMA * sig_dplx)
            & (dmu < mu_max + MU_MARGIN_NSIGMA * sig_dmu))

    return pd.DataFrame({"i": i[keep], "j": j[keep],
                         "theta_arcsec": theta[keep], "sep_au": sep_au[keep],
                         "dist_pc": dmean[keep]})


def drop_shared_components(pairs: pd.DataFrame) -> pd.DataFrame:
    """Keep only pairs whose members appear in no other pair, so pairs are
    statistically independent (removes triples and chance overlaps)."""
    i, j = pairs["i"].to_numpy(), pairs["j"].to_numpy()
    counts = pd.Series(np.concatenate([i, j])).value_counts()
    multi = set(counts[counts > 1].index)
    solo = ~(pd.Series(i).isin(multi).to_numpy()
             | pd.Series(j).isin(multi).to_numpy())
    return pairs[solo].reset_index(drop=True)
