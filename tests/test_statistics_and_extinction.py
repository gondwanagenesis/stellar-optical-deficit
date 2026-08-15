"""Statistics, extinction law, blinding."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import anchor                    # noqa: E402
from pipeline import extinction as ext         # noqa: E402
from pipeline import statistics as st          # noqa: E402


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------

def test_delta_mag_roundtrip():
    for f in [1e-6, 1e-3, 0.1, 0.5, 0.9]:
        assert st.fraction_from_delta(st.delta_mag(f)) == pytest.approx(f)


def test_delta_mag_signs_and_limits():
    assert st.delta_mag(0.0) == 0.0
    assert st.delta_mag(0.5) == pytest.approx(0.7525, abs=1e-3)
    assert st.delta_mag(0.9) == pytest.approx(2.5, abs=1e-3)
    # harvesting makes stars FAINTER, i.e. positive magnitude change
    assert st.delta_mag(0.1) > 0


def test_detection_efficiency_is_monotone_in_f():
    rng = np.random.default_rng(0)
    r = rng.normal(0, 0.1, 200_000)
    effs = [st.detection_efficiency(r, f, k=5.0) for f in
            [0.01, 0.05, 0.1, 0.3, 0.5, 0.9]]
    assert all(b >= a for a, b in zip(effs, effs[1:]))
    assert effs[0] < 1e-3          # a 1% deficit never reaches 5 sigma
    assert effs[-1] > 0.99         # a 90% deficit always does


def test_poisson_upper_limit_known_values():
    # Standard one-sided 95% Poisson upper limits
    assert st.poisson_upper_limit(0) == pytest.approx(2.996, abs=0.01)
    assert st.poisson_upper_limit(1) == pytest.approx(4.744, abs=0.01)
    assert st.poisson_upper_limit(10) == pytest.approx(17.0, abs=0.1)


def test_exclusion_curve_saturates_and_is_conservative():
    rng = np.random.default_rng(1)
    r = rng.normal(0, 0.1, 100_000)
    ex = st.exclusion_curve(r, [0.01, 0.1, 0.5, 0.9], k=5.0)
    # tiny f -> no constraint
    assert ex.iloc[0]["p_upper_limit"] >= 0.99
    # large f -> a real constraint
    assert ex.iloc[-1]["p_upper_limit"] < 1e-3
    # limits must tighten as f grows
    assert list(ex["p_upper_limit"]) == sorted(ex["p_upper_limit"], reverse=True)


def test_shape_stats_detects_an_injected_positive_tail():
    rng = np.random.default_rng(4)
    r = rng.normal(0, 0.1, 200_000)
    clean = st.shape_stats(r)
    r2 = r.copy()
    r2[:400] += 0.8
    dirty = st.shape_stats(r2)
    a_clean = clean.table.set_index("k").loc[5, "asymmetry"]
    a_dirty = dirty.table.set_index("k").loc[5, "asymmetry"]
    assert a_dirty > a_clean
    assert dirty.table.set_index("k").loc[5, "n_pos"] > 300


# --------------------------------------------------------------------------
# extinction
# --------------------------------------------------------------------------

def test_fitz19_kG_is_strongly_colour_dependent():
    """If this collapses to a constant, the whole colour systematic changes."""
    a0 = np.array([0.3, 0.3, 0.3])
    x = np.array([0.0, 1.5, 3.0])
    k = ext.k_fitz19("G", x, a0)
    assert k[0] > k[-1]
    assert k[0] / k[-1] > 1.3, f"k_G barely varies: {k}"


def test_fitz19_kKs_is_much_smaller_than_kG():
    a0 = np.full(3, 0.3)
    x = np.array([0.7, 2.0, 3.5])
    kg = ext.k_fitz19("G", x, a0)
    kks = ext.k_fitz19("Ks", x, a0)
    assert np.all(kks < 0.35)
    assert np.all(kg / kks > 2.0)


def test_four_iterations_is_not_enough():
    """Documents why the default was raised from 4 to 12."""
    a0 = np.linspace(0.01, 1.5, 40)
    bp_rp = np.linspace(0.6, 3.8, 40)
    a4 = ext.deredden("G", a0, bp_rp, law="fitz19", n_iter=4)
    a12 = ext.deredden("G", a0, bp_rp, law="fitz19", n_iter=12)
    assert np.abs(a4 - a12).max() > 1e-3, (
        "if this now converges by iteration 4 the default can be lowered again")


def test_deredden_iteration_converges_at_the_default():
    a0 = np.linspace(0.01, 1.5, 40)
    bp_rp = np.linspace(0.6, 3.8, 40)
    a12 = ext.deredden("G", a0, bp_rp, law="fitz19", n_iter=12)
    a30 = ext.deredden("G", a0, bp_rp, law="fitz19", n_iter=30)
    assert np.abs(a12 - a30).max() < 1e-4


def test_wangchen_law_is_colour_independent():
    a0 = np.full(5, 0.4)
    x = np.linspace(0.5, 3.5, 5)
    k = ext.k_wangchen19("G", x, a0)
    assert np.allclose(k, k[0])


def test_the_two_laws_disagree_in_the_nir():
    """The factor ~2.5 disagreement in A_Ks is a real published one."""
    a0 = np.full(3, 0.3)
    x = np.array([1.0, 2.0, 3.0])
    kf = ext.k_fitz19("Ks", x, a0).mean()
    kw = ext.k_wangchen19("Ks", x, a0).mean()
    assert kf / kw > 2.0


# --------------------------------------------------------------------------
# spectral leverage
# --------------------------------------------------------------------------

def test_grey_absorber_has_negative_leverage():
    """The brief's 'flat absorber is invisible' claim, tested."""
    for s in [1.2, 1.5, 1.9]:
        lev = anchor.leverage(0.0, s)
        assert lev < 0, f"slope {s} gave leverage {lev:+.3f}"
        assert abs(lev) > 0.1, "grey absorber is not even close to invisible"


def test_blind_spot_is_where_the_formula_says():
    for s in [1.2, 1.5, 1.9]:
        a_an = anchor.blind_slope_analytic(s)
        a_num = anchor.blind_slope_numeric(s)
        assert a_an == pytest.approx(a_num, abs=0.03)
        assert anchor.leverage(a_num, s) == pytest.approx(0.0, abs=0.02)


def test_dustlike_absorber_has_high_leverage():
    lev = anchor.leverage(2.0, 1.25)
    assert lev > 0.8, "under-corrected reddening should mimic the signal"


def test_leverage_is_monotone_in_alpha():
    s = 1.25
    a = np.linspace(0.0, 4.0, 25)
    lev = [anchor.leverage(x, s) for x in a]
    assert all(b >= a_ - 1e-9 for a_, b in zip(lev, lev[1:]))


def test_grey_absorber_dims_both_bands_equally():
    dg = anchor.band_deficit("G", 0.0, 0.05)
    dk = anchor.band_deficit("Ks", 0.0, 0.05)
    assert dg == pytest.approx(dk, rel=1e-6)


# --------------------------------------------------------------------------
# blinding
# --------------------------------------------------------------------------

def test_blind_commitment_verifies(tmp_path, monkeypatch):
    from pipeline import blind
    from pipeline import config as cfg
    monkeypatch.setattr(cfg, "BLIND_DIR", tmp_path)
    monkeypatch.setattr(blind, "SECRET_PATH", tmp_path / "secret_offset.json")
    monkeypatch.setattr(blind, "COMMIT_PATH", tmp_path / "commitment.json")

    blind.create(force=True)
    assert blind.verify()
    r = np.zeros(100)
    blinded = blind.apply(r)
    assert not np.allclose(blinded, r)          # it actually shifted something
    assert np.allclose(blinded - blinded[0], 0)  # by a constant

    with pytest.raises(ValueError):
        blind.unblind("nope")
    off = blind.unblind("I have frozen the analysis and committed it")
    assert blinded[0] == pytest.approx(off)
    assert blind.status().unblinded


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
