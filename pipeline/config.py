"""Central configuration for the stellar optical-deficit search.

Every threshold used anywhere in the pipeline lives here as a named constant with
a justification and, where it is an empirical quantity, a literature citation.
Nothing downstream may hard-code a number.

The measurement
---------------
For main-sequence stars we compare the absolute G magnitude to the absolute Ks
magnitude.  If a fraction ``f`` of the star's *optical* output is intercepted and
never leaves the system, then

    Delta_M_G = -2.5 * log10(1 - f)          [positive = fainter]

M_Ks acts as the "how big is this star really" anchor.  This is the load-bearing
assumption and it is also the principal weakness of the method: a spectrally
*flat* absorber dims G and Ks together and is completely invisible here.  See
LIMITATIONS.md and ``pipeline/flat_absorber.py``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
CACHE_DIR = DATA_DIR / "cache"
RAW_DIR = CACHE_DIR / "raw_gaia"
DERIVED_DIR = CACHE_DIR / "derived"
FIG_DIR = PROJECT_ROOT / "figures"
RESULT_DIR = PROJECT_ROOT / "results"
BLIND_DIR = PROJECT_ROOT / "blind"
MANIFEST_PATH = CACHE_DIR / "manifest.jsonl"

# dustmaps data is multi-GB and randomly accessed; keep it on native ext4 rather
# than the /mnt/c DrvFs mount, which is roughly an order of magnitude slower for
# the HDF5 random reads Bayestar does.
DUSTMAPS_DATA_DIR = Path(os.environ.get("DUSTMAPS_DATA_DIR", Path.home() / "dustmaps_data"))

for _d in (DATA_DIR, CACHE_DIR, RAW_DIR, DERIVED_DIR, FIG_DIR, RESULT_DIR, BLIND_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# TAP service
# --------------------------------------------------------------------------

TAP_URL = "https://gea.esac.esa.int/tap-server/tap"

# The ESA Gaia archive caps anonymous asynchronous jobs at 3e6 rows.  We chunk
# the sky far below that so a single chunk never truncates silently.
TAP_ROW_CAP = 3_000_000

# Sky partition.  Gaia source_id encodes a level-12 HEALPix (nested, nside=4096)
# index in its high bits:  source_id = hpx_level12 * 2**35 + ...
# so  hpx_level_N = source_id / 2**(59 - 2N).   Partitioning on source_id ranges
# lets the archive use its primary-key index, and gives us contiguous sky
# patches for free -- which the hemisphere/crowding null tests in step 3 need.
HEALPIX_PARTITION_LEVEL = 2          # nside=4 -> 192 chunks
SOURCE_ID_DIVISOR = 2 ** (59 - 2 * HEALPIX_PARTITION_LEVEL)
N_PARTITIONS = 12 * 4 ** HEALPIX_PARTITION_LEVEL

TAP_MAX_RETRIES = 5
TAP_RETRY_BASE_SLEEP = 20.0          # seconds; exponential backoff


# --------------------------------------------------------------------------
# Sample cuts.  Applied server-side in ADQL wherever possible.
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class SampleCuts:
    # --- astrometry -------------------------------------------------------
    # parallax_over_error > 20 makes the inverse-parallax distance bias
    # (Bailer-Jones 2015, PASP 127, 994) smaller than 0.4%, i.e. < 0.01 mag in
    # distance modulus, without needing a distance prior that would imprint its
    # own spatial structure on the residuals.
    parallax_over_error_min: float = 20.0

    # RUWE < 1.4 is the standard unresolved-binary / bad-astrometric-solution
    # rejector (Lindegren 2018, GAIA-C3-TN-LU-LL-124).  Note the sign of the
    # residual bias this leaves behind: an unresolved companion makes a star
    # OVER-luminous in G, which is the *opposite* of a harvesting deficit, so
    # residual binary contamination biases the search conservative.  This claim
    # is verified numerically in tests/test_binary_bias.py rather than assumed.
    ruwe_max: float = 1.4

    # Astrometric solution quality.  ipd_frac_multi_peak flags sources whose
    # image-parameter-determination windows contained a second peak, i.e.
    # marginally resolved doubles that RUWE can miss.
    ipd_frac_multi_peak_max: int = 2      # percent
    ipd_frac_odd_win_max: int = 7         # percent
    astrometric_excess_noise_sig_max: float = 2.0
    visibility_periods_used_min: int = 10

    # --- distance ---------------------------------------------------------
    # Primary sample: 500 pc.  Extinction is the dominant systematic and its
    # map-to-map disagreement grows with distance; 500 pc keeps the bulk of the
    # sample at A_V < 0.5 while still giving ~10^7 stars.
    distance_max_pc: float = 500.0
    distance_min_pc: float = 10.0         # below this Gaia saturates (G < 4)

    # Trade-study sample for the sensitivity-vs-systematics study (step 3).
    distance_max_pc_wide: float = 1250.0  # Edenhofer+2023 map outer boundary

    # --- photometry -------------------------------------------------------
    # Faint-end floor.  At G > 19 the Gaia photometric uncertainty exceeds
    # 0.01 mag and BP/RP become unreliable for the excess-factor cut.
    g_mag_max: float = 19.0
    g_mag_min: float = 4.0                # brighter than this, Gaia saturates
    g_flux_over_error_min: float = 100.0  # -> sigma_G < 0.011 mag

    # 2MASS Ks.  ph_qual 'A' means SNR > 10 and a good profile fit
    # (Skrutskie et al. 2006, AJ 131, 1163).  cc_flg '0' = no contamination or
    # confusion.  These are the load-bearing 2MASS quality flags.
    ks_sigma_max: float = 0.05
    tmass_ph_qual_ks: str = "A"
    tmass_cc_flg_clean: str = "0"

    # BP/RP excess factor.  C* = phot_bp_rp_excess_factor corrected for colour
    # per Riello et al. 2021 (A&A 649, A3), eq. 6 and Table 2; sources are kept
    # within N sigma of the corrected locus.  Rejects blends and extended
    # objects whose BP/RP flux exceeds the G aperture.
    bp_rp_excess_nsigma_max: float = 3.0

    # --- stellar type -----------------------------------------------------
    # Lower main sequence box in (M_Ks, (BP-RP)_0).  Boundaries chosen in
    # DECISIONS.md; the goal is a locus that is single-valued in M_Ks, avoids
    # the subgiant branch and the pre-main-sequence, and has dense enough
    # sampling that the fiducial spline is not extrapolating anywhere.
    m_ks_min: float = 3.0                 # ~ G0V; brighter -> turnoff/subgiants
    m_ks_max: float = 8.0                 # ~ M4V; fainter -> 2MASS incomplete
    bp_rp0_min: float = 0.7               # ~ G2V
    bp_rp0_max: float = 3.6               # ~ M4V
    # Vertical half-width of the MS band, in mag, used for the initial
    # (pre-fit) rejection of giants / white dwarfs / obvious binaries.
    ms_band_halfwidth_mag: float = 1.0

    # --- variability / multiplicity --------------------------------------
    reject_variable: bool = True          # phot_variable_flag = 'VARIABLE'
    reject_non_single_star: bool = True   # non_single_star != 0
    reject_duplicated_source: bool = True

    # --- extinction -------------------------------------------------------
    # Above this the colour-dependent extinction law's own uncertainty
    # (~5% of A_G) exceeds the target sensitivity, so the star is dropped from
    # the primary sample and used only in the high-extinction null test.
    a_g_max_primary: float = 0.5


CUTS = SampleCuts()


# --------------------------------------------------------------------------
# Empirical constants.  Every one of these is somebody's measurement.
# --------------------------------------------------------------------------

# Band effective wavelengths (micron), Vega-referenced.
#
# Gaia DR3 passbands from the SVO Filter Profile Service (Rodrigo & Solano
# 2020), GAIA/GAIA3 set, which tabulates the Riello et al. 2021 (A&A 649, A3)
# revised response curves.  NOTE these are the DR3 values; the DR2-era numbers
# (G = 0.673, BP = 0.532, RP = 0.797) are still widely quoted and are NOT the
# same.
#
# 2MASS: Cohen, Wheaton & Megeath 2003, AJ 126, 1090, Table 2.
#
# CAVEAT that matters for pipeline/anchor.py: lambda_eff is SED-dependent, and
# these are Vega-referenced.  The Gaia G band is extraordinarily broad
# (400-950 nm), so for the K and M dwarfs that dominate this sample the
# flux-weighted effective wavelength is far redder than the Vega value -- close
# to the old DR2 number, by coincidence.  Anything that needs a wavelength for
# THESE stars must integrate the actual SED through the band rather than use
# this constant; see anchor.effective_wavelength().
LAMBDA_EFF_UM = {
    "G": 0.58224,
    "BP": 0.50358,
    "RP": 0.76200,
    "J": 1.235,
    "H": 1.662,
    "Ks": 2.159,
}

# Pivot wavelengths, same sources. Included because the pivot wavelength is
# SED-independent and is the right choice when a single number must stand in
# for a band irrespective of the source spectrum.
LAMBDA_PIVOT_UM = {
    "G": 0.62176,
    "BP": 0.51097,
    "RP": 0.77690,
}

# Extinction law A -- constant ratios A_band / A_V for R_V = 3.1.
# Wang & Chen 2019, ApJ 877, 116, Table 3 (optical/NIR extinction law from
# Gaia DR2 + 2MASS + APOGEE red clump stars).
EXTINCTION_RATIOS_WANG_CHEN_2019 = {
    "G": 0.789,
    "BP": 1.002,
    "RP": 0.589,
    "J": 0.243,
    "H": 0.131,
    "Ks": 0.078,
}

# Dust-map unit conversions.
#
# Bayestar19 (Green et al. 2019, ApJ 887, 93) returns a reddening in its own
# arbitrary unit.  Green et al. 2019 section 5 / Table 1 give per-band
# extinction coefficients for that unit directly; the value below converts
# Bayestar19 units to A_V for R_V = 3.3, which is what their coefficients
# assume.  Reference: Green et al. 2019, Table 1 caption.
BAYESTAR19_TO_AV = 2.742 * 1.0   # A_V = 2.742 * E_bayestar  (Green+2019 Tab. 1)

# Edenhofer et al. 2023 (A&A 685, A82) returns a unitless extinction density
# integrated along the line of sight, on the Zhang, Green & Rix (2023) E-scale.
# Their published conversion to A_V (Johnson V, R_V=3.1) is A_V = 2.8 * E.
# Reference: Edenhofer et al. 2023, section 2 and the dustmaps documentation.
EDENHOFER23_TO_AV = 2.8

# Solar metallicity reference for [M/H] normalisation.
# Asplund et al. 2009, ARA&A 47, 481.
SOLAR_Z = 0.0134


# --------------------------------------------------------------------------
# Fiducial-relation fit
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class FitConfig:
    # Natural cubic spline in M_Ks, linear (then quadratic, CV-selected) in
    # [M/H].  Knot count is chosen by k-fold CV over this grid -- see
    # DECISIONS.md for why a spline rather than a GP (10^7 points makes an
    # exact GP intractable and a sparse GP adds its own inducing-point
    # systematic that would be indistinguishable from the signal).
    knot_grid: tuple = (6, 8, 10, 12, 16, 20, 25, 30, 40)
    cv_folds: int = 5
    spline_degree: int = 3
    # Robust loss: the MS has a genuine binary sequence 0.75 mag above it that
    # least-squares would drag the fiducial into.  Huber with this threshold in
    # units of the running sigma.
    huber_delta: float = 2.0
    max_irls_iter: int = 20
    irls_tol: float = 1e-6
    # Metallicity term
    mh_degree_grid: tuple = (1, 2, 3)
    # Outlier flag threshold for the individual-candidate search.
    outlier_nsigma: float = 5.0


FIT = FitConfig()


# --------------------------------------------------------------------------
# Null tests (step 3) and injection (step 4)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class NullConfig:
    # Subsample sizes at which mean-residual scatter is measured, to find where
    # the 1/sqrt(N) scaling breaks.
    n_grid: tuple = (100, 300, 1000, 3000, 10_000, 30_000,
                     100_000, 300_000, 1_000_000, 3_000_000)
    n_realisations: int = 200
    rng_seed: int = 20260815


@dataclass(frozen=True)
class InjectionConfig:
    f_grid: tuple = (1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
    # Sparse case: fraction p of stars harvested at f_sparse.
    sparse_p_grid: tuple = (1e-5, 1e-4, 1e-3, 1e-2)
    sparse_f_grid: tuple = (0.1, 0.3, 0.5)
    n_realisations: int = 50
    rng_seed: int = 8675309


NULLS = NullConfig()
INJECT = InjectionConfig()
