"""Test statistics for the optical-deficit search.

What is actually measurable
---------------------------
The fiducial spline has a free intercept in every knot interval, so the
weighted mean residual is zero *by construction* at every M_Ks.  This is not a
detail: it means

    the mean residual is not a measurable quantity in the self-calibrated
    analysis, and neither is a uniform harvesting fraction.

Even a sparse injection -- a fraction p of stars each dimmed by Delta -- shifts
the local mean by p*Delta, and the spline absorbs that too.  What the spline
*cannot* absorb is the **shape** of the residual distribution at fixed M_Ks.
Harvesting adds a one-sided positive (fainter) tail and nothing else.

So the self-calibrated analysis is a distribution-shape test:

  * tail asymmetry  A(k) = [N(r > +k s) - N(r < -k s)] / N
  * the count of individual >5-sigma positive outliers
  * an exclusion contour in the (f, p) plane, where f is the harvested
    fraction per star and p is the fraction of stars harvested

The null for A(k) is *not* zero.  Unresolved binaries and blends make stars
over-luminous, so the intrinsic residual distribution has a genuine negative
tail.  Any limit that estimates the positive-tail background by reflecting the
negative tail therefore over-subtracts and is anti-conservative.  The primary
limit here uses the positive tail alone -- every positive outlier is treated as
a potential signal -- and is conservative by construction.  The reflected
version is reported alongside, labelled as such.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


def delta_mag(f: float | np.ndarray) -> np.ndarray:
    """Magnitude deficit for an intercepted optical fraction f (positive = fainter)."""
    f = np.asarray(f, dtype=float)
    return -2.5 * np.log10(1.0 - f)


def fraction_from_delta(delta: float | np.ndarray) -> np.ndarray:
    """Inverse of delta_mag."""
    return 1.0 - 10.0 ** (-np.asarray(delta, dtype=float) / 2.5)


def robust_sigma(r: np.ndarray) -> float:
    r = r[np.isfinite(r)]
    return 1.4826 * float(np.median(np.abs(r - np.median(r))))


# --------------------------------------------------------------------------
# Shape statistics
# --------------------------------------------------------------------------

@dataclass
class ShapeStats:
    n: int
    sigma_robust: float
    mean: float
    median: float
    skew: float
    table: pd.DataFrame        # per-threshold tail counts

    def summary(self) -> str:
        return (f"N={self.n}  sigma_robust={self.sigma_robust:.4f} mag  "
                f"mean={self.mean:+.5f}  median={self.median:+.5f}  "
                f"skew={self.skew:+.4f}")


def shape_stats(resid: np.ndarray, k_values=(2, 3, 4, 5, 6)) -> ShapeStats:
    r = np.asarray(resid, dtype=float)
    r = r[np.isfinite(r)]
    n = len(r)
    s = robust_sigma(r)
    med = float(np.median(r))
    rows = []
    for k in k_values:
        n_pos = int(np.count_nonzero(r > med + k * s))
        n_neg = int(np.count_nonzero(r < med - k * s))
        asym = (n_pos - n_neg) / n
        # Poisson error on the difference of two counts
        err = np.sqrt(n_pos + n_neg) / n
        rows.append({
            "k": k, "n_pos": n_pos, "n_neg": n_neg,
            "frac_pos": n_pos / n, "frac_neg": n_neg / n,
            "asymmetry": asym, "asymmetry_err": err,
            "n_sigma": asym / err if err > 0 else np.nan,
        })
    # Robust skew (Bowley/quartile skewness): insensitive to the extreme tails
    q1, q2, q3 = np.percentile(r, [25, 50, 75])
    bowley = (q3 + q1 - 2 * q2) / (q3 - q1) if q3 > q1 else np.nan
    return ShapeStats(n=n, sigma_robust=s, mean=float(np.mean(r)),
                      median=med, skew=float(bowley), table=pd.DataFrame(rows))


# --------------------------------------------------------------------------
# Detection efficiency and exclusion in the (f, p) plane
# --------------------------------------------------------------------------

def detection_efficiency(resid: np.ndarray, f: float, k: float,
                         sigma: float | None = None,
                         median: float | None = None) -> float:
    """Fraction of harvested stars that would land above the +k sigma threshold.

    Estimated by shifting the *observed* residual distribution by Delta(f) and
    counting, so the true non-Gaussian shape of the intrinsic scatter is used
    rather than a Gaussian assumption.
    """
    r = np.asarray(resid, dtype=float)
    r = r[np.isfinite(r)]
    s = sigma if sigma is not None else robust_sigma(r)
    med = median if median is not None else float(np.median(r))
    thresh = med + k * s
    return float(np.mean(r + delta_mag(f) > thresh))


def poisson_upper_limit(n_obs: int, cl: float = 0.95) -> float:
    """One-sided upper limit on a Poisson mean given n_obs, no background.

    Uses the exact relation UL = 0.5 * chi2.ppf(cl, 2*(n+1)).
    """
    from scipy.stats import chi2
    return float(0.5 * chi2.ppf(cl, 2 * (n_obs + 1)))


def exclusion_curve(resid: np.ndarray, f_grid, k: float = 5.0,
                    cl: float = 0.95, subtract_reflected: bool = False) -> pd.DataFrame:
    """95% CL upper limit on the harvested fraction p, for each per-star f.

    p_UL(f) = N_UL / (N_total * efficiency(f, k))

    With ``subtract_reflected=False`` (the default and the reported primary),
    every star above the threshold is counted as a candidate signal: the limit
    is conservative and model-free.  With ``subtract_reflected=True`` the
    negative tail is used as a background estimate, which is *anti*-conservative
    because unresolved binaries populate the negative tail only.
    """
    r = np.asarray(resid, dtype=float)
    r = r[np.isfinite(r)]
    n = len(r)
    s = robust_sigma(r)
    med = float(np.median(r))
    n_pos = int(np.count_nonzero(r > med + k * s))
    n_neg = int(np.count_nonzero(r < med - k * s))

    if subtract_reflected:
        excess = max(n_pos - n_neg, 0)
        n_ul = poisson_upper_limit(excess, cl) + np.sqrt(n_neg)
    else:
        n_ul = poisson_upper_limit(n_pos, cl)

    rows = []
    for f in f_grid:
        eff = detection_efficiency(r, f, k, sigma=s, median=med)
        p_ul = n_ul / (n * eff) if eff > 0 else np.inf
        rows.append({
            "f": f,
            "delta_mag": float(delta_mag(f)),
            "delta_over_sigma": float(delta_mag(f) / s),
            "efficiency": eff,
            "n_pos_observed": n_pos,
            "n_neg_observed": n_neg,
            "n_upper_limit": n_ul,
            "p_upper_limit": min(p_ul, 1.0),
            "mean_f_upper_limit": min(p_ul, 1.0) * f,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Individual outliers
# --------------------------------------------------------------------------

def find_outliers(df: pd.DataFrame, resid: np.ndarray,
                  n_sigma: float = 5.0) -> pd.DataFrame:
    """Individual positive (fainter-than-fiducial) outliers.

    Also returns the matching negative-side list, which is the control: any
    contaminant that produces a symmetric tail should appear in both.
    """
    r = np.asarray(resid, dtype=float)
    s = robust_sigma(r)
    med = float(np.median(r))
    out = df.copy()
    out["residual"] = r
    out["residual_sigma"] = (r - med) / s
    out["implied_f"] = fraction_from_delta(r - med)
    pos = out[out["residual_sigma"] > n_sigma].copy()
    neg = out[out["residual_sigma"] < -n_sigma].copy()
    pos["side"] = "positive (deficit-like)"
    neg["side"] = "negative (over-luminous control)"
    return pd.concat([pos, neg], ignore_index=True).sort_values(
        "residual_sigma", ascending=False).reset_index(drop=True)
