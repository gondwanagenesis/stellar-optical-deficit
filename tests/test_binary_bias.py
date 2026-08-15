"""Does unresolved binarity bias the search conservative, or does it MIMIC the signal?

The brief asserted that "an unresolved companion makes a star OVERluminous,
which is the opposite of the signal, so contamination biases you conservative",
and asked for that to be verified rather than assumed.

It is false for this particular observable, and the sign is the dangerous one.

The residual is  r = M_G - f(M_Ks),  and  f' = dM_G/dM_Ks > 1 on the lower main
sequence (G is a bluer, more temperature-sensitive band than Ks, so M_G runs
faster).  Adding an unresolved companion brightens both magnitudes:

    r = dM_G - f' * dM_Ks

A *cool* companion contributes far more Ks flux than G flux, so
|dM_Ks| > |dM_G|, and with f' > 1 both terms push r positive.  Positive means
fainter-than-fiducial, which is exactly the harvesting signature.

Even the equal-mass case is adverse: both magnitudes drop by 0.753, giving
r = 0.753 * (f' - 1) > 0.

So unresolved binaries are a *primary* contaminant that manufactures a false
deficit, and the RUWE cut is load-bearing in the opposite sense to the one
stated.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def flux(mag):
    return 10.0 ** (-0.4 * np.asarray(mag, dtype=float))


def combine(mag_a, mag_b):
    return -2.5 * np.log10(flux(mag_a) + flux(mag_b))


# An illustrative monotone lower-main-sequence relation, shaped like the
# Pecaut & Mamajek (2013, ApJS 208, 9) dwarf sequence.  These are NOT precision
# values and nothing quantitative in the analysis uses them; the test depends
# only on the *shape* -- M_G rising faster than M_Ks -- which is a robust
# feature of any main sequence because G is bluer and more
# temperature-sensitive than Ks.  The measured slope from our own fitted
# fiducial is reported in RESULTS.md and is what the science uses.
#
# The table is deliberately extended well below the analysis box (down to
# M_Ks = 1.6) so that combining two stars never pushes the summed magnitude off
# the end of the table, where np.interp would silently clamp and fabricate a
# zero slope.
MS_M_KS = np.array([1.60, 2.10, 2.60, 3.00, 3.60, 4.20, 4.80,
                    5.40, 6.00, 6.60, 7.20, 7.80, 8.40])
MS_M_G = np.array([1.95, 2.62, 3.42, 4.10, 5.20, 6.10, 7.05,
                   8.05, 9.10, 10.20, 11.35, 12.55, 13.80])

MS_MIN, MS_MAX = MS_M_KS[0], MS_M_KS[-1]


def ms_m_g(m_ks):
    m_ks = np.asarray(m_ks, dtype=float)
    if np.any(m_ks < MS_MIN - 1e-9) or np.any(m_ks > MS_MAX + 1e-9):
        raise ValueError(
            f"M_Ks outside tabulated range [{MS_MIN}, {MS_MAX}]; np.interp "
            f"would clamp and fabricate a zero slope")
    return np.interp(m_ks, MS_M_KS, MS_M_G)


def ms_slope(m_ks, h=0.05):
    return (ms_m_g(m_ks + h) - ms_m_g(m_ks - h)) / (2 * h)


def test_main_sequence_slope_exceeds_unity():
    """The whole sign argument rests on dM_G/dM_Ks > 1."""
    m_ks = np.linspace(3.2, 7.6, 40)
    s = ms_slope(m_ks)
    assert np.all(s > 1.0), f"slope dips to {s.min():.3f}"


def binary_residual(m_ks_primary, m_ks_secondary):
    """Residual of an unresolved pair against the single-star relation."""
    g1, g2 = ms_m_g(m_ks_primary), ms_m_g(m_ks_secondary)
    g_tot = combine(g1, g2)
    ks_tot = combine(m_ks_primary, m_ks_secondary)
    return g_tot - ms_m_g(ks_tot)


def test_equal_mass_binary_looks_fainter_not_brighter():
    """Equal-mass binary: r = 0.753 (f' - 1) > 0, i.e. deficit-like."""
    for m_ks in [3.5, 4.5, 5.5, 6.5, 7.5]:
        r = binary_residual(m_ks, m_ks)
        assert r > 0.05, f"equal-mass binary at M_Ks={m_ks} gave r={r:+.4f}"


def test_cool_companion_looks_fainter():
    """Unequal pairs with a cooler secondary are worse, not better."""
    worst = 0.0
    for m_ks in [3.5, 4.0, 4.5, 5.0]:
        for dm in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
            sec = min(m_ks + dm, MS_M_KS[-1])
            r = binary_residual(m_ks, sec)
            assert r > -1e-6, (f"primary M_Ks={m_ks}, secondary M_Ks={sec} "
                               f"gave r={r:+.4f}")
            worst = max(worst, r)
    # The contaminant is not marginal: it reaches several tenths of a magnitude.
    assert worst > 0.15, f"peak binary residual only {worst:.3f} mag"


def test_binary_contamination_mimics_a_large_harvested_fraction():
    """Translate the worst-case binary residual into an apparent f."""
    from pipeline.statistics import fraction_from_delta
    rs = [binary_residual(4.0, 4.0 + dm) for dm in np.arange(0.25, 3.25, 0.25)]
    f_apparent = fraction_from_delta(max(rs))
    # If this is small the contaminant would be harmless; it is not.
    assert f_apparent > 0.10, f"apparent f only {f_apparent:.3f}"


def test_the_brief_claim_is_false():
    """Explicitly record that the assumed sign is wrong, so a regression
    that silently 'fixes' the sign is caught."""
    residuals = [binary_residual(4.5, 4.5 + dm) for dm in [0.0, 1.0, 2.0, 3.0]]
    assert all(r > 0 for r in residuals), (
        "if this ever passes with negative residuals, the main-sequence shape "
        "assumption has changed and the contamination argument must be redone")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
