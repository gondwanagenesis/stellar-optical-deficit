"""Does Search V's alignment machinery actually work, and does its null hold?

Three things have to be true before a null from channel 22 means anything, and
none of them is obvious enough to assume.

1. The estimator has to RECOVER an injected alignment. A test that cannot see
   a signal it was handed reports a null for the wrong reason. This project has
   already published one null that a positive control later showed was broken.

2. The local mark permutation has to REJECT a pure survey-geometry confound. The
   whole design rests on the claim that shuffling within sky cells preserves the
   smooth positional fields in both angles. If it does not, the channel's honest
   null is not honest.

3. The global mark permutation has to FAIL on that same confound. If both nulls
   behaved identically the local one would be pointless machinery, and the
   argument for preferring it would be decoration rather than method.

The confound is simulated exactly as the physics describes it: both the IPD
harmonic phase and the proper-motion position angle are given a smooth
dependence on sky position -- the scanning law for one, solar reflex motion for
the other -- with no per-star relationship whatsoever. Any test that calls that
a detection would have manufactured a discovery out of survey geometry.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "searchV", ROOT / "scripts" / "85_searchV_harmonic_phase.py")
sv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sv)

N = 20_000
NCELL = 400
NPERM = 200


def _sky(rng):
    """Cell id per star, plus a smooth angle field indexed by cell."""
    cell = rng.integers(0, NCELL, N)
    field = rng.uniform(0, 180, NCELL)
    return cell, field


def test_recovers_injected_alignment():
    """An injected concentration at a known dphi must come out significant."""
    rng = np.random.default_rng(7)
    cell, _ = _sky(rng)
    pa = rng.uniform(0, 180, N)
    # 20% of sources have their phase locked to the PM axis within +-15 deg;
    # the rest are uniform. This is a deliberately weak, realistic injection.
    aligned = rng.random(N) < 0.20
    phase = np.where(aligned,
                     (pa + rng.normal(0, 15, N)) % 180.0,
                     rng.uniform(0, 180, N))

    r_obs = sv.resultant(phase, pa)
    null = sv.permute_local(phase, pa, cell, rng, NPERM)
    assert r_obs > null.mean() + 5 * null.std(), (
        f"injected 20% alignment not recovered: R={r_obs:.4f} vs "
        f"null {null.mean():.4f}+-{null.std():.4f}")
    assert sv.perm_p(r_obs, null) < 0.01

    ray = sv.rayleigh((phase - pa) % 180.0)
    assert ray["p"] < 1e-6
    # The recovered axis must be the injected one, near dphi = 0 (mod 180).
    pref = ray["preferred_dphi_deg"]
    assert min(pref, 180.0 - pref) < 10.0


def test_local_null_rejects_a_pure_geometry_confound():
    """Two smooth sky fields, no per-star link. The local null must not care."""
    rng = np.random.default_rng(11)
    cell, field = _sky(rng)
    # PA_PM: smooth field (solar reflex) plus per-star scatter.
    pa = (field[cell] + rng.normal(0, 20, N)) % 180.0
    # IPD phase: a DIFFERENT smooth field that is correlated with the first
    # only through sky position -- here, the same field rotated and rescaled,
    # which is the worst realistic case.
    phase = (0.8 * field[cell] + 30.0 + rng.normal(0, 25, N)) % 180.0

    r_obs = sv.resultant(phase, pa)
    local = sv.permute_local(phase, pa, cell, rng, NPERM)
    p_local = sv.perm_p(r_obs, local)
    assert p_local > 0.01, (
        f"local null flagged pure geometry as signal: p={p_local:.4g}, "
        f"R={r_obs:.4f} vs null {local.mean():.4f}+-{local.std():.4f}")


def test_global_null_would_have_manufactured_that_discovery():
    """The same confound, judged by the global shuffle, looks like a detection.

    This is the failure the channel is designed around, so it is asserted
    rather than described.
    """
    rng = np.random.default_rng(11)
    cell, field = _sky(rng)
    pa = (field[cell] + rng.normal(0, 20, N)) % 180.0
    phase = (0.8 * field[cell] + 30.0 + rng.normal(0, 25, N)) % 180.0

    r_obs = sv.resultant(phase, pa)
    glob = sv.permute_global(phase, pa, rng, NPERM)
    assert sv.perm_p(r_obs, glob) < 0.01, (
        "the global shuffle failed to be fooled, which would mean the "
        "geometry injected here is too weak to exercise the distinction")
    assert glob.mean() < sv.permute_local(phase, pa, cell, rng, NPERM).mean()


def test_perpendicular_mirror_is_degenerate():
    """Documented in the script and asserted here: R is rotation invariant.

    The perpendicular tail cannot serve as a mirror control for this statistic,
    so the channel must not be allowed to quote one by accident.
    """
    rng = np.random.default_rng(3)
    pa = rng.uniform(0, 180, 5000)
    phase = (pa + rng.normal(0, 20, 5000)) % 180.0
    assert sv.resultant(phase, pa) == pytest.approx(
        sv.resultant(phase, (pa + 90.0) % 180.0), rel=1e-12)


def test_permutation_absorbs_nonuniform_marginals_and_rayleigh_does_not():
    """The reason the analytic Rayleigh p is not this channel's test.

    In the real data neither marginal is uniform: the harmonic phase piles up
    on the 0/180 axis because the scanning law does not sample scan angles
    evenly, and PA_PM is shaped by solar reflex motion. Two INDEPENDENT axial
    variables with such marginals have an expected resultant near the product
    of their marginal resultants rather than zero. A test that assumes
    uniformity reads that product as significance.

    Asserted here on independent data with deliberately non-uniform marginals:
    the permutation null tracks the product, so the permutation p is
    unbothered, while the Rayleigh p claims a detection that does not exist.
    """
    rng = np.random.default_rng(41)
    n = 20_000
    cell = rng.integers(0, NCELL, n)
    # Both concentrated toward 0/180, and drawn completely independently.
    phase = (rng.normal(0, 42, n)) % 180.0
    pa = (rng.normal(0, 55, n)) % 180.0

    r_phase = sv.rayleigh(phase)["R"]
    r_pa = sv.rayleigh(pa)["R"]
    r_obs = sv.resultant(phase, pa)
    null = sv.permute_local(phase, pa, cell, rng, NPERM)

    # The null lands on the product of the marginals, which is what it is for.
    assert null.mean() == pytest.approx(r_phase * r_pa, abs=0.02)
    # So the honest test sees nothing, because there is nothing.
    assert sv.perm_p(r_obs, null) > 0.01
    # And the analytic p, which assumes uniformity, calls it a detection.
    assert sv.rayleigh((phase - pa) % 180.0)["p"] < 1e-6


def test_uniform_data_gives_a_calibrated_p_value():
    """No signal, no geometry: the local-null p must be roughly uniform."""
    rng = np.random.default_rng(23)
    ps = []
    for _ in range(20):
        cell = rng.integers(0, NCELL, 5000)
        pa = rng.uniform(0, 180, 5000)
        phase = rng.uniform(0, 180, 5000)
        ps.append(sv.perm_p(sv.resultant(phase, pa),
                            sv.permute_local(phase, pa, cell, rng, 100)))
    ps = np.array(ps)
    assert 0.2 < ps.mean() < 0.8, f"p-values not calibrated: mean {ps.mean():.3f}"
    assert np.count_nonzero(ps < 0.05) <= 3
