"""Null-test machinery: the floor diagnostics must behave on known inputs."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import nulls                    # noqa: E402


def toy(n=50_000, seed=0, ext_error=0.0):
    rng = np.random.default_rng(seed)
    a0 = rng.exponential(0.2, n)
    df = pd.DataFrame({
        "A_0": a0,
        "A_G": 0.78 * a0,
        "b": rng.uniform(-80, 80, n),
        "l": rng.uniform(0, 360, n),
        "dist_pc": rng.uniform(50, 500, n),
        "phot_g_mean_mag": rng.uniform(10, 19, n),
        "sky_density": rng.uniform(1, 100, n),
        "bp_rp0": rng.uniform(0.7, 3.5, n),
        "mh_gspphot": rng.normal(0, 0.3, n),
    })
    # residual = pure noise + a fractional extinction error leaking through
    resid = rng.normal(0, 0.1, n) + 0.89 * ext_error * df["A_G"].to_numpy()
    return df, resid


def test_extinction_slope_recovers_a_known_error():
    for eps in [0.05, 0.20, -0.10]:
        df, r = toy(ext_error=eps, seed=abs(int(eps * 100)))
        out = nulls.extinction_residual_slope(df, r)
        assert out["implied_fractional_extinction_error"] == pytest.approx(
            eps, abs=0.02), f"eps={eps} -> {out['implied_fractional_extinction_error']}"


def test_extinction_slope_is_zero_for_clean_residuals():
    df, r = toy(ext_error=0.0)
    out = nulls.extinction_residual_slope(df, r)
    assert abs(out["n_sigma"]) < 4


def test_random_subsamples_track_sigma_over_sqrt_n():
    """The control must reproduce the naive scaling; if it does not, the
    floor measurement is meaningless."""
    rng = np.random.default_rng(1)
    r = rng.normal(0, 0.1, 300_000)
    tab = nulls.random_subsample_scaling(r, n_grid=(100, 1000, 10_000),
                                         n_realisations=200)
    ratio = tab["rms_of_means"] / tab["sigma_over_sqrt_n"]
    assert np.all(np.abs(ratio - 1) < 0.15), ratio.tolist()


def test_group_scatter_detects_a_coherent_offset():
    """Structured groups must show excess when a per-group offset exists."""
    rng = np.random.default_rng(2)
    n, ngrp = 60_000, 60
    gid = rng.integers(0, ngrp, n)
    offsets = rng.normal(0, 0.02, ngrp)      # coherent 0.02 mag per group
    r = rng.normal(0, 0.1, n) + offsets[gid]
    obs, exp, ngroups, meann = nulls._group_mean_scatter(r, gid)
    assert obs > exp
    excess = np.sqrt(obs ** 2 - exp ** 2)
    assert excess == pytest.approx(0.02, abs=0.006)


def test_group_scatter_refuses_too_few_groups():
    rng = np.random.default_rng(3)
    r = rng.normal(0, 0.1, 10_000)
    gid = rng.integers(0, 3, 10_000)         # below MIN_GROUPS_FOR_RMS
    obs, exp, ngroups, meann = nulls._group_mean_scatter(r, gid)
    assert np.isnan(obs), "an RMS from 3 groups must not be reported"


def test_split_test_signs():
    r = np.concatenate([np.full(1000, 0.1), np.full(1000, -0.1)])
    a = np.zeros(2000, dtype=bool); a[:1000] = True
    out = nulls.split_test(r, a, ~a, "toy")
    assert out["difference"] == pytest.approx(0.2)
    assert out["mean_a"] == pytest.approx(0.1)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
