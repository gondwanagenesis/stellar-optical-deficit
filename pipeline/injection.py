"""Step 4: inject synthetic deficits and try to recover them with the identical
pipeline.

Two injection modes, and they behave completely differently:

  uniform : every star dimmed by Delta(f).  The fiducial is refit afterwards,
            and because a constant offset lies exactly in the span of the
            spline basis, the fit absorbs it and the recovered signal is zero
            for every f.  This is not a bug in the injection -- it is a
            demonstration of the degeneracy described in DECISIONS.md D11, and
            it is the reason the brief's "nonzero MEAN Delta_M_G" statistic is
            not measurable self-calibrated.

  sparse  : a fraction p of stars dimmed by Delta(f).  Here the mean shift
            p*Delta is NOT absorbed -- the robust (Huber) fit down-weights the
            injected outliers instead of following them, so ~77% of p*Delta
            survives into the residual mean.  Recovery is nevertheless measured
            in the tail, because the mean is also where the coherent
            systematics live and the tail limit is ~17x tighter.

Rule enforced here: the spline hyperparameters (knot count, metallicity degree)
are chosen once on the real data and then held FIXED for every injection.
Re-running cross-validation per injection would let the model adapt to the
signal and would inflate recovery.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from . import config as cfg
from . import fiducial as fid
from . import statistics as st

log = logging.getLogger(__name__)


@dataclass
class Recovery:
    mode: str
    f_injected: float
    p_injected: float
    delta_injected: float
    n_stars: int
    # self-calibrated observables
    mean_residual: float
    sigma_robust: float
    n_pos_5sig: int
    n_neg_5sig: int
    n_pos_excess_over_baseline: float
    p_recovered: float
    p_recovered_err: float
    # model-anchored observable (mean offset against a fixed reference curve)
    anchored_mean_offset: float
    anchored_f_recovered: float


def inject(m_g: np.ndarray, f: float, p: float,
           rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Return (injected M_G, boolean mask of which stars were harvested)."""
    delta = float(st.delta_mag(f))
    if p >= 1.0:
        mask = np.ones(len(m_g), dtype=bool)
    else:
        mask = rng.random(len(m_g)) < p
    out = np.asarray(m_g, dtype=float).copy()
    out[mask] += delta
    return out, mask


def _p_from_tail(n_pos: int, n_pos_baseline: float, n_tot: int,
                 efficiency: float) -> tuple[float, float]:
    if efficiency <= 0:
        return np.nan, np.nan
    excess = n_pos - n_pos_baseline
    err = np.sqrt(max(n_pos, 0) + max(n_pos_baseline, 0))
    return excess / (n_tot * efficiency), err / (n_tot * efficiency)


def run_one(m_ks: np.ndarray, covs, m_g: np.ndarray,
            *, mode: str, f: float, p: float,
            n_interior: int,
            baseline_resid: np.ndarray,
            baseline_n_pos: float,
            anchor_curve=None,
            k: float = 5.0,
            rng: np.random.Generator | None = None) -> Recovery:
    rng = rng or np.random.default_rng(0)
    delta = float(st.delta_mag(f))

    m_g_inj, mask = inject(m_g, f, p, rng)

    # Refit with FIXED hyperparameters -- see module docstring.
    fit = fid.fit_fiducial(m_ks, covs, m_g_inj, n_interior)
    resid = fit.residuals(m_ks, covs, m_g_inj)

    s = st.robust_sigma(resid)
    med = float(np.median(resid))
    n_pos = int(np.count_nonzero(resid > med + k * s))
    n_neg = int(np.count_nonzero(resid < med - k * s))

    eff = st.detection_efficiency(baseline_resid, f, k)
    p_rec, p_err = _p_from_tail(n_pos, baseline_n_pos, len(m_ks), eff)

    # Model-anchored: compare to a FIXED reference curve that does not adapt.
    if anchor_curve is not None:
        offset = float(np.mean(m_g_inj - anchor_curve))
        f_anch = float(st.fraction_from_delta(max(offset, 0.0)))
    else:
        offset, f_anch = np.nan, np.nan

    return Recovery(
        mode=mode, f_injected=f, p_injected=p, delta_injected=delta,
        n_stars=len(m_ks),
        mean_residual=float(np.mean(resid)), sigma_robust=s,
        n_pos_5sig=n_pos, n_neg_5sig=n_neg,
        n_pos_excess_over_baseline=n_pos - baseline_n_pos,
        p_recovered=p_rec, p_recovered_err=p_err,
        anchored_mean_offset=offset, anchored_f_recovered=f_anch,
    )


def injection_campaign(df: pd.DataFrame, fit0: fid.FiducialFit,
                       resid0: np.ndarray,
                       n_interior: int, mh_degree: int,
                       cov_columns=("mh_gspphot",),
                       f_grid=None, sparse_p_grid=None, sparse_f_grid=None,
                       n_realisations: int | None = None,
                       k: float = 5.0,
                       seed: int | None = None,
                       max_n: int | None = None) -> pd.DataFrame:
    """Full injection-recovery campaign.

    ``max_n`` subsamples for speed; when set, it is recorded in the output so
    that no recovery curve can be quoted at an N it was not measured at.
    """
    f_grid = f_grid or cfg.INJECT.f_grid
    sparse_p_grid = sparse_p_grid or cfg.INJECT.sparse_p_grid
    sparse_f_grid = sparse_f_grid or cfg.INJECT.sparse_f_grid
    n_realisations = n_realisations or cfg.INJECT.n_realisations
    seed = seed if seed is not None else cfg.INJECT.rng_seed

    rng = np.random.default_rng(seed)

    m_ks = df["M_Ks"].to_numpy(dtype=float)
    m_g = df["M_G"].to_numpy(dtype=float)
    cov_vals = [df[c].to_numpy(dtype=float) for c in cov_columns]
    degrees = [mh_degree] + [1] * (len(cov_vals) - 1)

    if max_n is not None and len(m_ks) > max_n:
        idx = rng.choice(len(m_ks), max_n, replace=False)
        m_ks, m_g = m_ks[idx], m_g[idx]
        cov_vals = [c[idx] for c in cov_vals]
        resid0 = resid0[idx]
        log.info("injection campaign on a %d-star subsample", max_n)

    covs = list(zip(cov_vals, degrees))
    s0 = st.robust_sigma(resid0)
    med0 = float(np.median(resid0))
    baseline_n_pos = float(np.count_nonzero(resid0 > med0 + k * s0))
    anchor = fit0.predict(m_ks, covs)        # frozen reference curve

    rows = []

    # --- uniform ---------------------------------------------------------
    for f in f_grid:
        for i in range(max(1, n_realisations // 10)):
            r = run_one(m_ks, covs, m_g, mode="uniform", f=f, p=1.0,
                        n_interior=n_interior,
                        baseline_resid=resid0, baseline_n_pos=baseline_n_pos,
                        anchor_curve=anchor, k=k, rng=rng)
            rows.append(asdict(r) | {"realisation": i})
        log.info("uniform  f=%.0e done", f)

    # --- sparse ----------------------------------------------------------
    for f in sparse_f_grid:
        for p in sparse_p_grid:
            for i in range(n_realisations):
                r = run_one(m_ks, covs, m_g, mode="sparse", f=f, p=p,
                            n_interior=n_interior,
                            baseline_resid=resid0, baseline_n_pos=baseline_n_pos,
                            anchor_curve=anchor, k=k, rng=rng)
                rows.append(asdict(r) | {"realisation": i})
            log.info("sparse   f=%.2f p=%.0e done", f, p)

    out = pd.DataFrame(rows)
    out["subsampled_to"] = max_n if max_n else len(df)
    return out


def summarise(campaign: pd.DataFrame) -> pd.DataFrame:
    """Collapse realisations into mean +/- std per injection setting."""
    g = campaign.groupby(["mode", "f_injected", "p_injected"], as_index=False)
    agg = g.agg(
        n_realisations=("p_recovered", "size"),
        delta_injected=("delta_injected", "first"),
        n_stars=("n_stars", "first"),
        mean_residual=("mean_residual", "mean"),
        mean_residual_std=("mean_residual", "std"),
        p_recovered=("p_recovered", "mean"),
        p_recovered_std=("p_recovered", "std"),
        n_pos_5sig=("n_pos_5sig", "mean"),
        anchored_f=("anchored_f_recovered", "mean"),
        anchored_f_std=("anchored_f_recovered", "std"),
    )
    agg["p_recovered_over_injected"] = np.where(
        agg["p_injected"] > 0, agg["p_recovered"] / agg["p_injected"], np.nan)
    return agg
