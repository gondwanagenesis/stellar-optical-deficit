"""Spectral leverage, the flat-absorber case, and the model-anchored analysis.

THE FLAT-ABSORBER CAVEAT IS NOT WHAT THE BRIEF SAYS IT IS
---------------------------------------------------------
The brief states that "a spectrally FLAT absorber suppresses numerator and
denominator together and is invisible to this test".  Working it through, that
is not right, and the true blind spot sits somewhere else entirely.

The measured residual for a star whose light is attenuated by dm_G and dm_Ks is

    r = dm_G - s * dm_Ks           where s = dM_G/dM_Ks along the fiducial

and s is NOT 1.  On the lower main sequence s ~ 1.5-2.0, because G is a bluer
and far more temperature-sensitive band than Ks.  Therefore:

  * flat absorber (dm_G = dm_Ks = d):   r = d (1 - s) < 0.
    A grey absorber makes a star look *over*-luminous in G at fixed M_Ks.  It
    is not invisible; it is anti-correlated with the signal, and roughly
    |1 - s| ~ 0.5-1.0 times as large as the naive expectation.

  * the true blind spot is  dm_G / dm_Ks = s,  i.e. an absorber that is already
    moderately selective.  For a power-law optical depth
    tau(lambda) ~ lambda^-alpha this is

        (lambda_G / lambda_Ks)^-alpha = s   =>   alpha_blind = ln(s) / ln(lambda_Ks/lambda_G)

    which for s = 1.8 and the Gaia G / 2MASS Ks effective wavelengths gives
    alpha_blind ~ 0.5.

  * alpha > alpha_blind gives a positive (deficit-like) residual;
    alpha < alpha_blind gives a negative one.

A consequence worth stating loudly: ordinary interstellar dust has
alpha ~ 2 (A_G/A_Ks ~ 10), far on the positive side.  Under-corrected reddening
therefore mimics the harvesting signal directly, which is why the extinction
treatment carries so much of the systematic budget.

None of this removes the need for an independent mass anchor -- it changes what
the anchor is for.  It is no longer "the flat case is invisible, bound it"; it
is "the alpha ~ 0.5 case is invisible, and every other alpha has a sign and an
amplitude the anchor can calibrate".
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg

# Vega-referenced effective wavelengths, micron (config.LAMBDA_EFF_UM).
# These are NOT the right numbers for this sample -- see effective_wavelength()
# below and the caveat in config.py -- and are kept only so the Vega-based
# analytic estimate can be reported alongside the SED-weighted one.
LAM_G_VEGA = cfg.LAMBDA_EFF_UM["G"]
LAM_KS = cfg.LAMBDA_EFF_UM["Ks"]

# Approximate rectangular band limits, micron, used for the SED-weighted
# integration cross-check.
# Gaia EDR3 passbands: Riello et al. 2021, A&A 649, A3, Fig. 2 (G spans the
# full 330-1050 nm range; the 10% points are ~400-950 nm).
# 2MASS: Cohen, Wheaton & Megeath 2003, AJ 126, 1090, Table 2.
BAND_LIMITS_UM = {
    "G": (0.400, 0.950),
    "Ks": (2.028, 2.290),
}


# --------------------------------------------------------------------------
# Absorber model
# --------------------------------------------------------------------------

def tau_of_lambda(lam_um, alpha: float, tau_ref: float,
                  lam_ref_um: float = LAM_G_VEGA) -> np.ndarray:
    """Power-law optical depth, normalised to tau_ref at lam_ref.

    The reference wavelength only sets the normalisation of tau; every
    quantity reported here is a *ratio* between bands, so the choice cancels.
    """
    return tau_ref * (np.asarray(lam_um, dtype=float) / lam_ref_um) ** (-alpha)


def _planck_lambda(lam_um, teff):
    """Planck B_lambda up to a constant; only the SED *shape* matters here."""
    lam_m = np.asarray(lam_um, dtype=float) * 1e-6
    h, c, k_b = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    x = h * c / (lam_m * k_b * teff)
    return lam_m ** -5 / np.expm1(np.clip(x, 1e-8, 700.0))


def band_deficit(band: str, alpha: float, tau_ref: float, teff: float = 4500.0,
                 n_grid: int = 400, sed_weighted: bool = True) -> float:
    """Magnitude deficit in ``band`` for a power-law absorber.

    ``sed_weighted=True`` integrates a blackbody through a rectangular band;
    ``False`` uses the effective wavelength only.  The two agree to better than
    0.01 mag for Ks at any tau, and for G they diverge once tau is large
    because G is very broad -- which is itself worth knowing.
    """
    if not sed_weighted:
        lam = cfg.LAMBDA_EFF_UM[band]
        return float(1.0857362 * tau_of_lambda(lam, alpha, tau_ref))
    lo, hi = BAND_LIMITS_UM[band]
    lam = np.linspace(lo, hi, n_grid)
    w = _planck_lambda(lam, teff)
    t = np.exp(-tau_of_lambda(lam, alpha, tau_ref))
    return float(-2.5 * np.log10(np.trapezoid(w * t, lam) / np.trapezoid(w, lam)))


def deficit_ratio(alpha: float, tau_ref: float = 0.05, teff: float = 4500.0,
                  sed_weighted: bool = True) -> float:
    """dm_G / dm_Ks for a given spectral slope."""
    dg = band_deficit("G", alpha, tau_ref, teff, sed_weighted=sed_weighted)
    dk = band_deficit("Ks", alpha, tau_ref, teff, sed_weighted=sed_weighted)
    return dg / dk if dk > 0 else np.inf


def leverage(alpha: float, slope: float, tau_ref: float = 0.05,
             teff: float = 4500.0, sed_weighted: bool = True) -> float:
    """Fraction of the raw G-band deficit that survives into the residual.

        leverage = (dm_G - slope * dm_Ks) / dm_G

    leverage = 1 would mean Ks is a perfect, unaffected mass anchor.
    leverage = 0 is the blind spot.  leverage < 0 means the residual has the
    opposite sign to the deficit.
    """
    dg = band_deficit("G", alpha, tau_ref, teff, sed_weighted=sed_weighted)
    dk = band_deficit("Ks", alpha, tau_ref, teff, sed_weighted=sed_weighted)
    return (dg - slope * dk) / dg if dg != 0 else np.nan


def effective_wavelength(band: str, teff: float = 4500.0,
                         n_grid: int = 2000) -> float:
    """Flux-weighted effective wavelength of a band for a blackbody of ``teff``.

    lambda_eff = int(lambda * B_lambda * R) / int(B_lambda * R)

    This is the wavelength at which the band actually samples the absorber's
    optical depth *for these stars*, and it differs substantially from the
    Vega-referenced catalogue value because Gaia G spans 400-950 nm and the
    sample is dominated by K and M dwarfs.  At 4500 K the G band's effective
    wavelength lands near 0.67 um rather than the Vega value of 0.582 um.
    """
    lo, hi = BAND_LIMITS_UM[band]
    lam = np.linspace(lo, hi, n_grid)
    w = _planck_lambda(lam, teff)
    return float(np.trapezoid(lam * w, lam) / np.trapezoid(w, lam))


def blind_slope_analytic(slope: float, teff: float | None = 4500.0) -> float:
    """alpha at which the test has exactly zero leverage (thin-absorber limit).

    With ``teff`` given, the SED-weighted effective wavelengths are used, which
    is the physically correct choice.  Pass ``teff=None`` to use the
    Vega-referenced catalogue values instead; the two differ by ~10% in the
    resulting alpha, and the difference is entirely the broadness of G.
    """
    if teff is None:
        lam_g, lam_ks = LAM_G_VEGA, LAM_KS
    else:
        lam_g = effective_wavelength("G", teff)
        lam_ks = effective_wavelength("Ks", teff)
    return float(np.log(slope) / np.log(lam_ks / lam_g))


def blind_slope_numeric(slope: float, tau_ref: float = 0.05,
                        teff: float = 4500.0,
                        bracket=(-2.0, 6.0)) -> float:
    from scipy.optimize import brentq
    fn = lambda a: leverage(a, slope, tau_ref, teff)      # noqa: E731
    lo, hi = bracket
    try:
        return float(brentq(fn, lo, hi, xtol=1e-4))
    except ValueError:
        return np.nan


def leverage_table(slope: float, alphas=None, tau_ref: float = 0.05,
                   teff: float = 4500.0) -> pd.DataFrame:
    alphas = alphas if alphas is not None else np.array(
        [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0])
    rows = []
    for a in alphas:
        dg = band_deficit("G", a, tau_ref, teff)
        dk = band_deficit("Ks", a, tau_ref, teff)
        rows.append({
            "alpha": a,
            "dm_G": dg, "dm_Ks": dk,
            "ratio_dmG_dmKs": dg / dk if dk else np.inf,
            "residual_per_dmG": (dg - slope * dk) / dg if dg else np.nan,
            "residual_mag": dg - slope * dk,
        })
    df = pd.DataFrame(rows)
    df.attrs["slope"] = slope
    df.attrs["tau_ref"] = tau_ref
    return df


# --------------------------------------------------------------------------
# Model-anchored offset
# --------------------------------------------------------------------------

def anchored_offset(m_g: np.ndarray, m_ks: np.ndarray, mh: np.ndarray,
                    reference) -> dict:
    """Mean offset of the sample from a FIXED external reference relation.

    Unlike the self-calibrated residual, this one retains sensitivity to a
    uniform offset -- at the price of inheriting the reference's systematic
    error, which is the dominant term and must be quoted with the result.
    """
    pred = reference(m_ks, mh)
    d = np.asarray(m_g, dtype=float) - pred
    d = d[np.isfinite(d)]
    n = len(d)
    mean = float(np.mean(d))
    return {
        "n": n,
        "mean_offset_mag": mean,
        "stat_err_mag": float(np.std(d, ddof=1) / np.sqrt(n)),
        "median_offset_mag": float(np.median(d)),
        "implied_f": float(1.0 - 10.0 ** (-max(mean, 0.0) / 2.5)),
    }
