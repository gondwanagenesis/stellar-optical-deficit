"""Extinction: 3D dust maps and band extinction coefficients.

Two independent axes of choice are exposed, and the pipeline runs all four
combinations so that the spread can be reported as a systematic rather than
hidden inside a single "best" answer:

  dust map   : bayestar19 | edenhofer23
  band law   : fitz19     | wangchen19

Why the band law matters so much
--------------------------------
The Gaia/Fitzpatrick-2019 coefficient k_G = A_G/A_0 runs from 1.00 at
(BP-RP)_0 = 0 to 0.65 at (BP-RP)_0 = 3.  A constant A_G/A_0 would therefore
imprint a strong colour-dependent residual on the main sequence, and a
colour-dependent residual is exactly what a mass-dependent harvesting signal
would look like.  Using the colour-dependent law is not optional.

The two laws also disagree substantially in the near-infrared:
Fitzpatrick et al. (2019) gives A_Ks/A_0 = 0.194 while Wang & Chen (2019)
gives A_Ks/A_V = 0.078.  That factor-of-2.5 disagreement is a real, published
disagreement about the NIR slope of the extinction curve, and at A_0 = 0.3 it
is a 0.035 mag shift in M_Ks -- far above the target sensitivity.  It is
propagated, not averaged away.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as cfg

log = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Band-extinction laws
# --------------------------------------------------------------------------

# Column order of the ESA Fitz19_EDR3_*.csv coefficient tables.
_FITZ_TERMS = ["Intercept", "X", "X2", "X3", "A", "A2", "A3", "XA", "AX2", "XA2"]

# Map our band names onto the table's Kname values.
_FITZ_BAND = {"G": "kG", "BP": "kBP", "RP": "kRP",
              "J": "kJ", "H": "kH", "Ks": "kK"}


@lru_cache(maxsize=1)
def _fitz_table() -> pd.DataFrame:
    path = cfg.DATA_DIR / "extinction_coeffs" / "Fitz19_EDR3_MainSequence.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing -- run scripts/04_fetch_extinction_coeffs.py")
    df = pd.read_csv(path)
    return df[df["Xname"] == "BPRP"].set_index("Kname")


def k_fitz19(band: str, bp_rp0: np.ndarray, a0: np.ndarray) -> np.ndarray:
    """A_band / A_0 from the official Gaia EDR3 (Fitzpatrick+2019) polynomial.

    k = a1 + a2 X + a3 X^2 + a4 X^3 + a5 A0 + a6 A0^2 + a7 A0^3
          + a8 X A0 + a9 A0 X^2 + a10 X A0^2
    with X = (BP-RP)_0 (the *intrinsic*, de-reddened colour).

    Source: ESA Gaia EDR3 extinction-law release,
    Fitz19_EDR3_MainSequence.csv, log g = 4.5 main-sequence variant.
    Valid 3500-10000 K and 0.01 < A0 < 20 mag.
    """
    row = _fitz_table().loc[_FITZ_BAND[band]]
    x = np.asarray(bp_rp0, dtype=float)
    a = np.asarray(a0, dtype=float)
    return (row["Intercept"]
            + row["X"] * x + row["X2"] * x**2 + row["X3"] * x**3
            + row["A"] * a + row["A2"] * a**2 + row["A3"] * a**3
            + row["XA"] * x * a + row["AX2"] * a * x**2 + row["XA2"] * x * a**2)


def k_wangchen19(band: str, bp_rp0: np.ndarray, a0: np.ndarray) -> np.ndarray:
    """A_band / A_V, constant, from Wang & Chen 2019, ApJ 877, 116, Table 3.

    Colour-independent by construction.  Used as the systematic alternative;
    A_V is identified with A_0 here, which is accurate to a few per cent for
    R_V = 3.1 and is itself part of the systematic being probed.
    """
    r = cfg.EXTINCTION_RATIOS_WANG_CHEN_2019[band]
    return np.full(np.shape(np.asarray(bp_rp0, dtype=float)), r, dtype=float)


BAND_LAWS = {"fitz19": k_fitz19, "wangchen19": k_wangchen19}


def deredden(band: str, a0: np.ndarray, bp_rp_obs: np.ndarray,
             law: str = "fitz19", n_iter: int = 4) -> np.ndarray:
    """Return A_band given A_0 and the *observed* BP-RP.

    The Fitz19 law wants the intrinsic colour, which we do not know until we
    have dereddened, so iterate.  Four iterations converges to < 1e-4 mag over
    the whole colour and A_0 range used here (checked in tests).
    """
    fn = BAND_LAWS[law]
    a0 = np.asarray(a0, dtype=float)
    bp_rp = np.asarray(bp_rp_obs, dtype=float)
    bp_rp0 = bp_rp.copy()
    for _ in range(n_iter):
        e_bp_rp = (fn("BP", bp_rp0, a0) - fn("RP", bp_rp0, a0)) * a0
        bp_rp0 = bp_rp - e_bp_rp
    return fn(band, bp_rp0, a0) * a0


def intrinsic_bp_rp(a0: np.ndarray, bp_rp_obs: np.ndarray,
                    law: str = "fitz19", n_iter: int = 4) -> np.ndarray:
    fn = BAND_LAWS[law]
    a0 = np.asarray(a0, dtype=float)
    bp_rp = np.asarray(bp_rp_obs, dtype=float)
    bp_rp0 = bp_rp.copy()
    for _ in range(n_iter):
        e_bp_rp = (fn("BP", bp_rp0, a0) - fn("RP", bp_rp0, a0)) * a0
        bp_rp0 = bp_rp - e_bp_rp
    return bp_rp0


# --------------------------------------------------------------------------
# 3D dust maps
# --------------------------------------------------------------------------

def _configure_dustmaps() -> None:
    from dustmaps.config import config as dm_config
    dm_config["data_dir"] = str(cfg.DUSTMAPS_DATA_DIR)


@lru_cache(maxsize=4)
def _query_object(which: str):
    _configure_dustmaps()
    if which == "bayestar19":
        from dustmaps.bayestar import BayestarQuery
        return BayestarQuery(version="bayestar2019", max_samples=1)
    if which == "edenhofer23":
        from dustmaps.edenhofer2023 import Edenhofer2023Query
        # integrated=True is REQUIRED.  The default returns extinction *density*
        # per parsec, which is ~1e-4 and silently produces an almost-zero A_0;
        # that bug is invisible in a summary table because the numbers merely
        # look small rather than wrong.  Building the integrated map takes a
        # couple of minutes on first construction.
        return Edenhofer2023Query(integrated=True)
    if which == "sfd":
        from dustmaps.sfd import SFDQuery
        return SFDQuery()
    raise ValueError(which)


# Native-unit -> A_0 conversion for each map.
MAP_TO_A0 = {
    "bayestar19": cfg.BAYESTAR19_TO_AV,     # Green et al. 2019, ApJ 887, 93, Tab. 1
    "edenhofer23": cfg.EDENHOFER23_TO_AV,   # Edenhofer et al. 2024, A&A 685, A82
    "sfd": 3.1,                             # nominal R_V; 2D, sanity check only
}

# Approximate validity limits, used to set a coverage mask rather than to
# silently extrapolate.
MAP_DISTANCE_MAX_PC = {"bayestar19": 63_000.0, "edenhofer23": 1_250.0, "sfd": np.inf}
MAP_DISTANCE_MIN_PC = {"bayestar19": 63.0, "edenhofer23": 69.0, "sfd": 0.0}


def query_a0(which: str, l_deg: np.ndarray, b_deg: np.ndarray,
             dist_pc: np.ndarray, chunk: int = 200_000) -> np.ndarray:
    """A_0 from a 3D map at (l, b, distance).  NaN where out of coverage."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    q = _query_object(which)
    scale = MAP_TO_A0[which]
    n = len(l_deg)
    out = np.full(n, np.nan)

    for i0 in range(0, n, chunk):
        i1 = min(i0 + chunk, n)
        sl = slice(i0, i1)
        if which == "sfd":
            coords = SkyCoord(l_deg[sl] * u.deg, b_deg[sl] * u.deg, frame="galactic")
        else:
            coords = SkyCoord(l_deg[sl] * u.deg, b_deg[sl] * u.deg,
                              distance=dist_pc[sl] * u.pc, frame="galactic")
        with np.errstate(invalid="ignore"):
            vals = q(coords)
        out[sl] = np.asarray(vals, dtype=float) * scale
        log.debug("%s: %d/%d", which, i1, n)

    lo = MAP_DISTANCE_MIN_PC[which]
    hi = MAP_DISTANCE_MAX_PC[which]
    bad = (dist_pc < lo) | (dist_pc > hi)
    out[bad] = np.nan
    # Bayestar returns NaN outside its declination footprint already.
    out[out < 0] = 0.0
    return out
