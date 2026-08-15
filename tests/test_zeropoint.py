"""The parallax zero point must be applied in the right units.

This exists because it was wrong: `gaiadr3-zeropoint` returns milliarcseconds
(the same units as `parallax`), and an extra /1000 made the correction 1000x
too small. The failure mode is nasty because nothing crashes -- the pipeline
just quietly runs with a ~1.4% distance-scale error, which is 0.03 mag of
distance modulus and comparable to the whole systematic floor.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config as cfg          # noqa: E402
from pipeline import sample as smp          # noqa: E402


def _one_chunk() -> pd.DataFrame | None:
    files = sorted(cfg.RAW_DIR.glob("sample_d*_p*.parquet"))
    return pd.read_parquet(files[0]) if files else None


needs_data = pytest.mark.skipif(_one_chunk() is None,
                                reason="no cached Gaia chunk available")


@needs_data
def test_zero_point_magnitude_is_literature_scale():
    df = _one_chunk()
    zp = smp.parallax_zero_point(df)          # mas
    med_uas = float(np.median(zp)) * 1000.0
    # Lindegren et al. 2021, A&A 649, A4: the DR3 parallax bias is of order
    # -17 uas globally and reaches -50 uas for faint red sources.
    assert -80.0 < med_uas < -5.0, f"median zero point {med_uas:.2f} uas"


@needs_data
def test_zero_point_units_are_mas_not_uas():
    """The units bug directly: zp must be comparable to parallax/100, not /1e5."""
    df = _one_chunk()
    zp = smp.parallax_zero_point(df)
    plx = df["parallax"].to_numpy(float)
    frac = abs(np.median(zp) / np.median(plx))
    # ~1% for a 3 mas parallax and a 40 uas offset; the bug gave ~1e-5
    assert 1e-3 < frac < 1e-1, f"|zp/parallax| = {frac:.2e}"


@needs_data
def test_correction_moves_the_distance_modulus_measurably():
    df = _one_chunk()
    out = smp.add_astrometry(df)
    naive_mu = 5 * np.log10(1000.0 / df["parallax"].to_numpy(float)) - 5
    shift = float(np.median(out["dist_mod"].to_numpy(float) - naive_mu))
    # a negative zero point means true parallaxes are LARGER, so distances and
    # distance moduli are SMALLER
    assert shift < 0, f"zero point moved the distance modulus by {shift:+.4f} mag"
    assert 0.005 < abs(shift) < 0.15, f"shift {shift:+.4f} mag is implausible"


@needs_data
def test_bad_units_are_rejected_loudly(monkeypatch):
    """If the package ever switches to uas, we must fail rather than proceed."""
    df = _one_chunk()
    import pipeline.sample as s
    from zero_point import zpt

    real = zpt.get_zpt

    def scaled(*a, **k):
        return np.asarray(real(*a, **k), dtype=float) / 1000.0

    monkeypatch.setattr(zpt, "get_zpt", scaled)
    with pytest.raises(ValueError, match="plausible range"):
        s.parallax_zero_point(df)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
