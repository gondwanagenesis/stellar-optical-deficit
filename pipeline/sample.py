"""Turn raw TAP chunks into the analysis sample.

Order of operations matters and is fixed here:
  1. parallax zero-point correction   (Lindegren et al. 2021)
  2. distance and distance modulus
  3. A_0 from each 3D dust map
  4. A_band from each band law, iterating for the intrinsic colour
  5. absolute magnitudes
  6. photometric-quality cuts that need corrected quantities
  7. main-sequence box

Every stage records how many stars it removed, and the cut-flow table is a
reported result, not a debug print.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import config as cfg
from . import extinction as ext

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Riello et al. 2021, A&A 649, A3 -- BP/RP flux excess factor
# --------------------------------------------------------------------------

def corrected_excess_factor(bp_rp: np.ndarray, excess: np.ndarray) -> np.ndarray:
    """C* -- the colour-corrected BP/RP flux excess factor, Riello+2021 eq. 6.

    A well-behaved single star has C* ~ 0 regardless of colour.  Positive C*
    means BP+RP flux exceeds what the G aperture sees, i.e. a blend, a nearby
    neighbour, or an extended source.  Because such contamination makes a star
    look *brighter*, not fainter, it works against the deficit signal -- but it
    also inflates the scatter, which is what actually limits us.
    """
    x = np.asarray(bp_rp, dtype=float)
    c = np.asarray(excess, dtype=float)
    out = np.empty_like(x)
    lo = x < 0.5
    mid = (x >= 0.5) & (x < 4.0)
    hi = x >= 4.0
    out[lo] = c[lo] - (1.154360 + 0.033772 * x[lo] + 0.032277 * x[lo] ** 2)
    out[mid] = c[mid] - (1.162004 + 0.011464 * x[mid] + 0.049255 * x[mid] ** 2
                         - 0.005879 * x[mid] ** 3)
    out[hi] = c[hi] - (1.057572 + 0.140537 * x[hi])
    return out


def excess_factor_sigma(g_mag: np.ndarray) -> np.ndarray:
    """1-sigma width of the C* distribution vs G, Riello+2021 eq. 18."""
    g = np.asarray(g_mag, dtype=float)
    return 0.0059898 + 8.817481e-12 * g ** 7.618399


# --------------------------------------------------------------------------
# Parallax zero point
# --------------------------------------------------------------------------

def parallax_zero_point(df: pd.DataFrame) -> np.ndarray:
    """Lindegren et al. 2021, A&A 649, A4 parallax bias, in mas.

    The published correction is calibrated for 6 < G < 21.  Our sample floor is
    G = 4, so for G < 6 we hold the correction at its G = 6 value rather than
    extrapolating a fitted spline outside its support.  Fewer than 0.1% of the
    sample is affected (checked in the cut-flow report).

    Why this matters at our precision: the zero point is -20 to -40 uas, which
    at parallax 2 mas is a 1-2% distance error, i.e. 0.02-0.04 mag of distance
    modulus.  That is an order of magnitude above the sensitivity we are trying
    to reach, so it cannot be neglected.  Note also that a distance-modulus
    error moves M_G and M_Ks *together*, so it enters the residual only through
    (1 - dM_G/dM_Ks) -- a partial but not complete cancellation.
    """
    from zero_point import zpt
    zpt.load_tables()

    g = df["phot_g_mean_mag"].to_numpy(dtype=float)
    g_clipped = np.clip(g, 6.0, 21.0)
    nu_eff = df["nu_eff_used_in_astrometry"].to_numpy(dtype=float)
    pseudo = df["pseudocolour"].to_numpy(dtype=float)
    ecl_lat = df["ecl_lat"].to_numpy(dtype=float)
    solved = df["astrometric_params_solved"].to_numpy(dtype=int)

    # zpt needs a finite colour proxy in the branch it will actually use.
    nu_eff = np.where(np.isfinite(nu_eff), nu_eff, 1.45)
    pseudo = np.where(np.isfinite(pseudo), pseudo, 1.45)

    zp = zpt.get_zpt(g_clipped, nu_eff, pseudo, ecl_lat, solved, _warnings=False)
    zp = np.asarray(zp, dtype=float)
    # Sources with astrometric_params_solved not in {31, 95} get no correction.
    zp = np.where(np.isfinite(zp), zp, 0.0)

    # gaiadr3-zeropoint returns the offset in MILLIarcseconds, the same units
    # as `parallax` -- not microarcseconds.  Dividing by 1000 here silently
    # reduced the correction by three orders of magnitude, which looked like
    # "the zero point is negligible" instead of like a bug.  Typical values are
    # -0.02 to -0.05 mas, i.e. -20 to -50 uas, matching Lindegren et al. 2021.
    med_uas = float(np.median(zp)) * 1000.0
    if not (-150.0 < med_uas < -1.0):
        raise ValueError(
            f"parallax zero point median {med_uas:.2f} uas is outside the "
            f"plausible range (-150, -1); the units returned by "
            f"gaiadr3-zeropoint have probably changed. Refusing to continue "
            f"with a silently wrong distance scale.")
    return zp


# --------------------------------------------------------------------------
# Cut-flow bookkeeping
# --------------------------------------------------------------------------

@dataclass
class CutFlow:
    rows: list = None

    def __post_init__(self):
        if self.rows is None:
            self.rows = []

    def record(self, label: str, mask: np.ndarray, note: str = "") -> None:
        n_before = len(mask)
        n_after = int(np.count_nonzero(mask))
        self.rows.append({
            "cut": label,
            "n_before": n_before,
            "n_after": n_after,
            "n_removed": n_before - n_after,
            "frac_removed": (n_before - n_after) / max(n_before, 1),
            "note": note,
        })
        log.info("%-42s %9d -> %9d  (-%6.3f%%)  %s",
                 label, n_before, n_after, 100 * (n_before - n_after) / max(n_before, 1),
                 note)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


# --------------------------------------------------------------------------
# Derived quantities
# --------------------------------------------------------------------------

def add_astrometry(df: pd.DataFrame) -> pd.DataFrame:
    zp = parallax_zero_point(df)
    df = df.copy()
    df["parallax_zp"] = zp
    df["parallax_corr"] = df["parallax"].to_numpy(dtype=float) - zp
    df["dist_pc"] = 1000.0 / df["parallax_corr"]
    df["dist_mod"] = 5.0 * np.log10(df["dist_pc"]) - 5.0
    return df


def add_extinction(df: pd.DataFrame, maps=("edenhofer23",),
                   laws=("fitz19", "wangchen19")) -> pd.DataFrame:
    """Attach A_0 per map and A_band per (map, law) combination."""
    df = df.copy()
    l = df["l"].to_numpy(dtype=float)
    b = df["b"].to_numpy(dtype=float)
    d = df["dist_pc"].to_numpy(dtype=float)
    bp_rp = df["bp_rp"].to_numpy(dtype=float)

    for m in maps:
        log.info("querying dust map %s for %d stars ...", m, len(df))
        a0 = ext.query_a0(m, l, b, d)
        df[f"a0_{m}"] = a0
        for law in laws:
            good = np.isfinite(a0)
            a0f = np.where(good, a0, 0.0)
            for band in ("G", "Ks", "BP", "RP"):
                a = ext.deredden(band, a0f, bp_rp, law=law)
                df[f"a_{band}_{m}_{law}"] = np.where(good, a, np.nan)
            df[f"bp_rp0_{m}_{law}"] = np.where(
                good, ext.intrinsic_bp_rp(a0f, bp_rp, law=law), np.nan)
    return df


def add_absolute_magnitudes(df: pd.DataFrame, dust_map: str, band_law: str) -> pd.DataFrame:
    """M_G and M_Ks for one (map, law) choice, written to fixed column names."""
    df = df.copy()
    mu = df["dist_mod"].to_numpy(dtype=float)
    df["M_G"] = (df["phot_g_mean_mag"].to_numpy(dtype=float) - mu
                 - df[f"a_G_{dust_map}_{band_law}"].to_numpy(dtype=float))
    df["M_Ks"] = (df["tmass_ks_m"].to_numpy(dtype=float) - mu
                  - df[f"a_Ks_{dust_map}_{band_law}"].to_numpy(dtype=float))
    df["bp_rp0"] = df[f"bp_rp0_{dust_map}_{band_law}"]
    df["A_0"] = df[f"a0_{dust_map}"]
    df["A_G"] = df[f"a_G_{dust_map}_{band_law}"]
    df["A_Ks"] = df[f"a_Ks_{dust_map}_{band_law}"]
    return df


def add_sky_density(df: pd.DataFrame, nside: int = 64) -> pd.DataFrame:
    """Crowding proxy: sample stars per HEALPix pixel (nside 64 ~ 0.84 deg^2).

    This is the density of the *selected* sample, not the true stellar density.
    It is used only to split the sample into crowded and sparse halves for the
    null test, where a monotone proxy is all that is required.
    """
    import healpy as hp
    df = df.copy()
    theta = np.radians(90.0 - df["b"].to_numpy(dtype=float))
    phi = np.radians(df["l"].to_numpy(dtype=float))
    pix = hp.ang2pix(nside, theta, phi, nest=True)
    counts = np.bincount(pix, minlength=hp.nside2npix(nside))
    df["hpx64"] = pix
    df["sky_density"] = counts[pix]
    return df


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------

def apply_quality_cuts(df: pd.DataFrame, flow: CutFlow) -> pd.DataFrame:
    c = cfg.CUTS
    n0 = len(df)
    flow.record("0. raw joined rows", np.ones(n0, dtype=bool), "after server-side ADQL cuts")

    m = np.ones(n0, dtype=bool)

    def step(cond, label, note=""):
        nonlocal m
        m = m & np.asarray(cond, dtype=bool)
        flow.record(label, m, note)

    step(df["tmass_ph_qual"].str[2] == c.tmass_ph_qual_ks,
         "1. 2MASS Ks ph_qual == 'A'", "SNR>10 and good profile fit")

    if c.reject_variable:
        step(df["phot_variable_flag"].fillna("") != "VARIABLE",
             "2. not photometrically variable", "Gaia phot_variable_flag")

    if c.reject_non_single_star:
        step(df["non_single_star"].fillna(0) == 0,
             "3. non_single_star == 0", "astrometric/spectro/eclipsing binary flags")

    step(~df["in_qso_candidates"].fillna(False) & ~df["in_galaxy_candidates"].fillna(False),
         "4. not QSO/galaxy candidate")

    step(df["classprob_dsc_combmod_star"].fillna(1.0) > 0.5,
         "5. DSC star probability > 0.5", "extragalactic/contaminant rejection")

    cstar = corrected_excess_factor(df["bp_rp"].to_numpy(dtype=float),
                                   df["phot_bp_rp_excess_factor"].to_numpy(dtype=float))
    sig = excess_factor_sigma(df["phot_g_mean_mag"].to_numpy(dtype=float))
    df = df.copy()
    df["cstar"] = cstar
    df["cstar_nsigma"] = cstar / sig
    step(np.abs(cstar) < c.bp_rp_excess_nsigma_max * sig,
         "6. |C*| < 3 sigma_C*(G)", "Riello+2021 BP/RP excess locus")

    return df[m].reset_index(drop=True)


def apply_extinction_and_ms_cuts(df: pd.DataFrame, flow: CutFlow,
                                 distance_max_pc: float | None = None,
                                 a_g_max: float | None = None) -> pd.DataFrame:
    c = cfg.CUTS
    dmax = distance_max_pc if distance_max_pc is not None else c.distance_max_pc
    agmax = a_g_max if a_g_max is not None else c.a_g_max_primary

    m = np.ones(len(df), dtype=bool)

    def step(cond, label, note=""):
        nonlocal m
        m = m & np.asarray(cond, dtype=bool)
        flow.record(label, m, note)

    step(np.isfinite(df["A_0"]), "7. dust map coverage", "finite A_0 at this l,b,d")
    step((df["dist_pc"] > c.distance_min_pc) & (df["dist_pc"] < dmax),
         f"8. {c.distance_min_pc:.0f} < d < {dmax:.0f} pc", "zero-point-corrected distance")
    step(df["A_G"] < agmax, f"9. A_G < {agmax}", "extinction-law error stays sub-threshold")
    step(np.isfinite(df["M_G"]) & np.isfinite(df["M_Ks"]) & np.isfinite(df["bp_rp0"]),
         "10. finite M_G, M_Ks, (BP-RP)_0")
    step((df["M_Ks"] > c.m_ks_min) & (df["M_Ks"] < c.m_ks_max),
         f"11. {c.m_ks_min} < M_Ks < {c.m_ks_max}", "lower main sequence box")
    step((df["bp_rp0"] > c.bp_rp0_min) & (df["bp_rp0"] < c.bp_rp0_max),
         f"12. {c.bp_rp0_min} < (BP-RP)_0 < {c.bp_rp0_max}")

    return df[m].reset_index(drop=True)
