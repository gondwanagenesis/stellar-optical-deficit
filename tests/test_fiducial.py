"""Fiducial-fit behaviour, including the degeneracy that defines the method."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import fiducial as fid          # noqa: E402
from pipeline import statistics as st         # noqa: E402


def synth(n=20_000, sigma=0.09, seed=1):
    """A toy main sequence: M_G rises faster than M_Ks, plus a metallicity term."""
    rng = np.random.default_rng(seed)
    m_ks = rng.uniform(3.0, 8.0, n)
    mh = rng.normal(-0.1, 0.3, n)
    truth = 1.1 + 1.35 * m_ks + 0.02 * (m_ks - 5.5) ** 2 - 0.45 * mh
    m_g = truth + rng.normal(0, sigma, n)
    return m_ks, mh, m_g, truth


def test_recovers_a_known_relation():
    m_ks, mh, m_g, truth = synth()
    fit = fid.fit_fiducial(m_ks, [(mh, 1)], m_g, 10)
    pred = fit.predict(m_ks, [(mh, 1)])
    assert np.abs(pred - truth).max() < 0.05
    assert fit.sigma_robust == pytest.approx(0.09, abs=0.01)


def test_slope_is_measured_correctly():
    m_ks, mh, m_g, _ = synth()
    fit = fid.fit_fiducial(m_ks, [(mh, 1)], m_g, 10)
    s = fid.slope(fit, m_ks, [(mh, 1)])
    # true dM_G/dM_Ks = 1.35 + 0.04*(m_ks - 5.5), median at m_ks=5.5 -> 1.35
    assert np.median(s) == pytest.approx(1.35, abs=0.03)


def test_uniform_offset_is_absorbed_essentially_exactly():
    """THE core degeneracy. A uniform deficit is invisible, by construction.

    The absorption is not bit-exact: the tiny ridge term used for conditioning
    slightly penalises the intercept, leaving a leak of ~1.5e-7 of the injected
    offset. That is seven orders of magnitude below the systematic floor, so
    the correct assertion is on the leak *fraction*, not on an absolute
    tolerance that would scale with the injected value.

    If this test ever fails, either the basis lost its constant term or someone
    fixed the intercept -- and in either case every 'uniform f' number in
    RESULTS.md has to be recomputed.
    """
    m_ks, mh, m_g, _ = synth()
    covs = [(mh, 1)]
    base = fid.fit_fiducial(m_ks, covs, m_g, 10)
    r_base = base.residuals(m_ks, covs, m_g)

    for f in [1e-5, 1e-3, 1e-2, 0.1, 0.5]:
        delta = st.delta_mag(f)
        shifted = m_g + delta
        fit = fid.fit_fiducial(m_ks, covs, shifted, 10)
        r = fit.residuals(m_ks, covs, shifted)
        leak = abs(np.mean(r) - np.mean(r_base)) / delta
        assert leak < 1e-5, f"f={f} leaked {leak:.2e} of its offset"
        assert np.abs(r - r_base).max() / delta < 1e-5


def test_sparse_injection_moves_BOTH_the_tail_and_the_mean():
    """A sparse deficit is NOT absorbed the way a uniform one is.

    This corrects an assumption that is easy to make and wrong: because the fit
    is robust (Huber), injected outliers sitting several sigma out get
    down-weighted rather than followed, so the fiducial barely moves and most
    of the p*delta mean shift SURVIVES into the residuals. Measured here at
    ~77% of p*delta.

    The tail is still the preferred statistic -- not because the mean is
    insensitive, but because the mean is where the coherent extinction and
    photometric systematics also live (see step 3).
    """
    m_ks, mh, m_g, _ = synth(n=60_000, seed=7)
    covs = [(mh, 1)]
    base = fid.fit_fiducial(m_ks, covs, m_g, 10)
    r_base = base.residuals(m_ks, covs, m_g)
    s = st.robust_sigma(r_base)
    n_pos_base = int((r_base > np.median(r_base) + 5 * s).sum())

    rng = np.random.default_rng(3)
    p, f = 0.01, 0.5
    delta = st.delta_mag(f)
    mask = rng.random(len(m_g)) < p
    inj = m_g.copy()
    inj[mask] += delta

    fit = fid.fit_fiducial(m_ks, covs, inj, 10)
    r = fit.residuals(m_ks, covs, inj)
    s2 = st.robust_sigma(r)
    n_pos = int((r > np.median(r) + 5 * s2).sum())

    assert n_pos > n_pos_base + 100, (
        f"tail did not respond: {n_pos_base} -> {n_pos}")

    retained = (np.mean(r) - np.mean(r_base)) / (p * delta)
    assert 0.5 < retained < 1.2, (
        f"sparse mean shift retained {retained:.2f} of p*delta; the robust fit "
        f"should preserve most of it")


def test_huber_resists_an_overluminous_sequence():
    """A binary sequence 0.75 mag brighter must not drag the fiducial."""
    m_ks, mh, m_g, truth = synth(n=40_000, seed=11)
    rng = np.random.default_rng(5)
    is_bin = rng.random(len(m_g)) < 0.15
    contaminated = m_g.copy()
    contaminated[is_bin] -= 0.75

    covs = [(mh, 1)]
    fit = fid.fit_fiducial(m_ks, covs, contaminated, 10)
    pred = fit.predict(m_ks, covs)
    bias = float(np.median(pred - truth))
    # least squares would shift by about -0.15*0.75 = -0.11 mag
    assert abs(bias) < 0.05, f"fiducial dragged by {bias:+.3f} mag"


def test_design_matrix_is_sparse_and_right_shape():
    import scipy.sparse as sp
    m_ks, mh, m_g, _ = synth(n=1000)
    knots = fid.quantile_knots(m_ks, 8)
    A = fid.design_matrix(m_ks, [(mh, 2)], knots)
    assert sp.issparse(A)
    n_basis = len(knots) - 3 - 1
    assert A.shape == (1000, n_basis * 3)
    # cubic B-spline: at most 4 non-zeros per basis block per row
    assert A.nnz <= 1000 * 4 * 3


def test_extra_covariate_reduces_scatter_when_it_matters():
    """A covariate that genuinely drives M_G should tighten the fit."""
    rng = np.random.default_rng(2)
    n = 30_000
    m_ks = rng.uniform(3, 8, n)
    mh = rng.normal(0, 0.3, n)
    jks = rng.normal(0.6, 0.15, n)
    m_g = 1.0 + 1.3 * m_ks - 0.4 * mh + 1.8 * jks + rng.normal(0, 0.05, n)

    without = fid.fit_fiducial(m_ks, [(mh, 1)], m_g, 10).sigma_robust
    with_ = fid.fit_fiducial(m_ks, [(mh, 1), (jks, 1)], m_g, 10).sigma_robust
    assert with_ < 0.5 * without


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
