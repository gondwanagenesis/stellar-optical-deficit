"""The fiducial main-sequence relation M_G = f(M_Ks, [M/H]) and its residuals.

Model
-----
    M_G = sum_k sum_j  c_{jk} B_k(M_Ks) * mh^j        j = 0 .. p

B_k is a cubic B-spline basis with knots at quantiles of M_Ks (so knot density
follows sample density and no knot interval is starved), and p is the
metallicity polynomial degree.  Both K and p are chosen by k-fold
cross-validation.

A DEGENERACY YOU MUST UNDERSTAND BEFORE READING ANY RESULT
----------------------------------------------------------
This fiducial is fit *to the data itself*.  A harvesting fraction f applied
uniformly to every star in the sample shifts M_G by a constant
-2.5 log10(1-f) and leaves M_Ks untouched; the fit simply absorbs that constant
into its intercept and the residuals return to zero.  The self-calibrated test
is therefore **identically blind to a uniform population offset**, and its
sensitivity to one is exactly zero, not sigma/sqrt(N).

What it *can* see is differential harvesting: a subpopulation offset from the
rest, or a deficit that varies with M_Ks or metallicity in a way the basis
cannot absorb.  A separate model-anchored analysis (pipeline/anchor.py) is
required to bound the uniform case, and that one is limited by stellar-model
systematics rather than by counting statistics.

Robustness
----------
The main sequence has a genuine equal-mass-binary sequence 0.75 mag above it.
Least squares would drag the fiducial upward toward it, so the fit uses Huber
IRLS.  Since binaries are over-luminous, any residual pull is toward *negative*
(brighter) residuals, which suppresses a deficit signal rather than creating
one.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.interpolate import BSpline

from . import config as cfg

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Design matrix
# --------------------------------------------------------------------------

def quantile_knots(x: np.ndarray, n_interior: int, degree: int = 3) -> np.ndarray:
    """Full knot vector with ``n_interior`` interior knots at data quantiles."""
    x = np.asarray(x, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    if n_interior > 0:
        qs = np.linspace(0, 1, n_interior + 2)[1:-1]
        interior = np.quantile(x, qs)
        # Collapse duplicate knots that arise if the sample piles up.
        interior = np.unique(interior)
    else:
        interior = np.array([])
    return np.concatenate([np.repeat(lo, degree + 1), interior,
                           np.repeat(hi, degree + 1)])


def _normalise_covs(covs) -> list[tuple[np.ndarray, int]]:
    """Accept ``None``, ``(values, degree)``, or a list of such pairs."""
    if covs is None:
        return []
    if isinstance(covs, tuple) and len(covs) == 2 and np.ndim(covs[1]) == 0:
        covs = [covs]
    out = []
    for values, deg in covs:
        if values is None or deg == 0:
            continue
        out.append((np.asarray(values, dtype=float), int(deg)))
    return out


def design_matrix(m_ks: np.ndarray, covs, knots: np.ndarray,
                  degree: int = 3) -> sp.csr_matrix:
    """Sparse design matrix: spline basis in M_Ks, tensored with covariate powers.

    ``covs`` is a list of ``(values, polynomial_degree)`` pairs -- metallicity,
    and optionally a near-infrared colour (see NIR CONTROL below).

    Sparse is not an optimisation detail here.  A cubic B-spline has only 4
    non-zero basis functions per point, so the matrix has 4*(1 + sum of degrees)
    non-zeros per row instead of (K+4)*(1 + sum of degrees).  At N = 5e6 and
    K = 20 that is the difference between fitting on the whole sample and not.

    NIR CONTROL
    -----------
    Adding the *optical* colour (BP-RP) as a control would be self-defeating:
    optical harvesting changes BP-RP, so the control would absorb the signal.
    A near-infrared colour such as (J-Ks) is safe for exactly the reason this
    whole method works -- the test is only sensitive to absorbers that are
    spectrally selective in the optical, and such an absorber leaves J-Ks
    alone.  (J-Ks) therefore controls for temperature, age and abundance
    structure at fixed M_Ks without touching the signal.  The cost is the
    2MASS photometric error it injects, which is why it is optional and its
    effect on the scatter is reported.
    """
    x = np.clip(np.asarray(m_ks, dtype=float), knots[0], knots[-1])
    B = BSpline.design_matrix(x, knots, degree, extrapolate=False).tocsr()
    pairs = _normalise_covs(covs)
    if not pairs:
        return B
    blocks = [B]
    for values, deg in pairs:
        for j in range(1, deg + 1):
            blocks.append(sp.diags(values ** j) @ B)
    return sp.hstack(blocks, format="csr")


# --------------------------------------------------------------------------
# Robust fit
# --------------------------------------------------------------------------

@dataclass
class FiducialFit:
    coef: np.ndarray
    knots: np.ndarray
    cov_degrees: tuple
    degree: int
    sigma_robust: float
    n_used: int
    n_params: int
    converged: bool
    ridge: float = 0.0

    @property
    def mh_degree(self) -> int:
        return self.cov_degrees[0] if self.cov_degrees else 0

    def predict(self, m_ks: np.ndarray, covs) -> np.ndarray:
        A = design_matrix(m_ks, covs, self.knots, self.degree)
        return np.asarray(A @ self.coef).ravel()

    def residuals(self, m_ks: np.ndarray, covs,
                  m_g: np.ndarray) -> np.ndarray:
        return np.asarray(m_g, dtype=float) - self.predict(m_ks, covs)


def _robust_sigma(r: np.ndarray) -> float:
    """Normalised MAD -- insensitive to the binary sequence and to outliers."""
    return 1.4826 * float(np.nanmedian(np.abs(r - np.nanmedian(r))))


def fit_fiducial(m_ks: np.ndarray, covs, m_g: np.ndarray,
                 n_interior: int,
                 degree: int | None = None,
                 huber_delta: float | None = None,
                 max_iter: int | None = None,
                 tol: float | None = None,
                 ridge: float = 1e-8) -> FiducialFit:
    degree = degree if degree is not None else cfg.FIT.spline_degree
    huber_delta = huber_delta if huber_delta is not None else cfg.FIT.huber_delta
    max_iter = max_iter if max_iter is not None else cfg.FIT.max_irls_iter
    tol = tol if tol is not None else cfg.FIT.irls_tol

    m_ks = np.asarray(m_ks, dtype=float)
    m_g = np.asarray(m_g, dtype=float)
    pairs = _normalise_covs(covs)

    knots = quantile_knots(m_ks, n_interior, degree)
    A = design_matrix(m_ks, covs, knots, degree)
    n, p = A.shape

    w = np.ones(n)
    coef = np.zeros(p)
    prev = None
    converged = False

    for _ in range(max_iter):
        Aw = sp.diags(w) @ A
        AtA = np.asarray((Aw.T @ A).todense())
        AtA.flat[:: p + 1] += ridge * n          # tiny ridge for conditioning
        Atb = np.asarray(Aw.T @ m_g).ravel()
        coef = np.linalg.solve(AtA, Atb)
        r = m_g - np.asarray(A @ coef).ravel()
        s = _robust_sigma(r)
        if s <= 0 or not np.isfinite(s):
            break
        u = np.abs(r) / s
        w = np.where(u <= huber_delta, 1.0, huber_delta / np.maximum(u, 1e-12))
        if prev is not None and np.max(np.abs(coef - prev)) < tol:
            converged = True
            break
        prev = coef.copy()

    r = m_g - np.asarray(A @ coef).ravel()
    return FiducialFit(coef=coef, knots=knots,
                       cov_degrees=tuple(d for _, d in pairs), degree=degree,
                       sigma_robust=_robust_sigma(r), n_used=n, n_params=p,
                       converged=converged, ridge=ridge)


# --------------------------------------------------------------------------
# Cross-validation over knot count and metallicity degree
# --------------------------------------------------------------------------

@dataclass
class CVResult:
    table: list = field(default_factory=list)
    best_knots: int = 0
    best_mh_degree: int = 0

    def to_frame(self):
        import pandas as pd
        return pd.DataFrame(self.table)


def cross_validate(m_ks: np.ndarray, mh: np.ndarray, m_g: np.ndarray,
                   knot_grid=None, mh_degree_grid=None, folds: int | None = None,
                   seed: int = 12345, max_fit_n: int = 400_000,
                   extra_covs: list | None = None) -> CVResult:
    """k-fold CV on robust (Huber) out-of-fold loss.

    The loss is the mean Huber loss, not the mean square, because the binary
    sequence is a real feature of the data and a squared loss would select the
    knot count that best chases it.

    For speed the CV runs on a random subsample of at most ``max_fit_n`` stars;
    the selected complexity is then refit on the full sample.  Knot count is
    only weakly N-dependent here because knots are placed at quantiles.
    """
    knot_grid = knot_grid or cfg.FIT.knot_grid
    mh_degree_grid = mh_degree_grid or cfg.FIT.mh_degree_grid
    folds = folds or cfg.FIT.cv_folds

    rng = np.random.default_rng(seed)
    n = len(m_ks)
    extra_covs = list(extra_covs or [])
    if n > max_fit_n:
        idx = rng.choice(n, max_fit_n, replace=False)
        m_ks, mh, m_g = m_ks[idx], mh[idx], m_g[idx]
        extra_covs = [c[idx] for c in extra_covs]
        n = max_fit_n

    fold_id = rng.integers(0, folds, size=n)
    res = CVResult()
    best = (np.inf, None, None)

    def build(sel, p):
        return [(mh[sel], p)] + [(c[sel], 1) for c in extra_covs]

    for k in knot_grid:
        for p in mh_degree_grid:
            losses = []
            for f in range(folds):
                tr = fold_id != f
                te = ~tr
                try:
                    fit = fit_fiducial(m_ks[tr], build(tr, p), m_g[tr], k)
                    r = fit.residuals(m_ks[te], build(te, p), m_g[te])
                except np.linalg.LinAlgError:
                    losses.append(np.inf)
                    continue
                s = fit.sigma_robust
                u = np.abs(r) / s
                d = cfg.FIT.huber_delta
                loss = np.where(u <= d, 0.5 * u ** 2, d * (u - 0.5 * d))
                losses.append(float(np.nanmean(loss)))
            mloss = float(np.mean(losses))
            res.table.append({"n_interior_knots": k, "mh_degree": p,
                              "cv_huber_loss": mloss,
                              "n_params": (k + cfg.FIT.spline_degree + 1) * (p + 1)})
            log.info("CV knots=%3d mh_deg=%d  loss=%.6f", k, p, mloss)
            if mloss < best[0]:
                best = (mloss, k, p)

    res.best_knots, res.best_mh_degree = best[1], best[2]
    log.info("CV best: knots=%d mh_degree=%d", res.best_knots, res.best_mh_degree)
    return res


# --------------------------------------------------------------------------
# Error budget
# --------------------------------------------------------------------------

def slope(fit: FiducialFit, m_ks: np.ndarray, covs,
          h: float = 0.01) -> np.ndarray:
    """dM_G/dM_Ks along the fitted relation, by central difference.

    Note this holds the covariates fixed, which is the relevant derivative:
    it is the slope a star moves along when its M_Ks changes at fixed
    metallicity and NIR colour.
    """
    return (fit.predict(m_ks + h, covs) - fit.predict(m_ks - h, covs)) / (2 * h)


def residual_uncertainty(df, fit: FiducialFit, covs=None) -> np.ndarray:
    """Per-star measurement contribution to the residual, in magnitudes.

        r = M_G - f(M_Ks)
        dr = dG - f' dKs - (1 - f') dmu - dA_G + f' dA_Ks

    Note the (1 - f') factor on the distance-modulus error: because M_G and
    M_Ks share the same distance, a parallax error moves a star nearly along
    the relation and is *partially* cancelled.  On the lower main sequence
    f' ~ 1.1-1.4, so the cancellation is strong but incomplete and the sign of
    the residual leakage flips across the colour range.
    """
    m_ks = df["M_Ks"].to_numpy(dtype=float)
    if covs is None:
        covs = [(df["mh_gspphot"].to_numpy(dtype=float), fit.mh_degree)]
    fprime = slope(fit, m_ks, covs)

    # Gaia G: sigma_mag = 1.0857 / (flux/flux_error), plus the calibration floor.
    # Calibration floor 2.0 mmag from Riello et al. 2021 Sect. 8.1.
    sig_g = np.sqrt((1.0857 / df["phot_g_mean_flux_over_error"].to_numpy(dtype=float)) ** 2
                    + 0.0020 ** 2)
    sig_ks = df["tmass_ks_msigcom"].to_numpy(dtype=float)
    # sigma_mu = (5/ln10) * sigma_plx/plx
    sig_mu = 2.1715 / df["parallax_over_error"].to_numpy(dtype=float)

    var = (sig_g ** 2
           + (fprime * sig_ks) ** 2
           + ((1.0 - fprime) * sig_mu) ** 2)
    return np.sqrt(var)


def intrinsic_scatter(resid: np.ndarray, sig_meas: np.ndarray) -> tuple[float, float]:
    """Deconvolve the measurement contribution from the observed scatter.

    Returns (sigma_observed, sigma_intrinsic), both robust (MAD-based).
    """
    s_obs = _robust_sigma(resid)
    s_meas = float(np.sqrt(np.nanmean(sig_meas ** 2)))
    s_int = float(np.sqrt(max(s_obs ** 2 - s_meas ** 2, 0.0)))
    return s_obs, s_int
