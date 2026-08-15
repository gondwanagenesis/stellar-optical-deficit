"""ADQL construction for the primary sample pull.

Design notes
------------
*Partitioning.*  Gaia ``source_id`` has a level-12 HEALPix index in its high
bits, so a ``source_id BETWEEN lo AND hi`` clause is simultaneously (a) a
primary-key range scan the archive can index and (b) a contiguous sky patch.
That second property is what makes the hemisphere and crowding null tests in
step 3 cheap: the partition id *is* a sky position.

*Which cuts go server-side.*  Everything that does not need extinction.  The
absolute-magnitude pre-box uses the *uncorrected* Ks, widened by the largest
plausible A_Ks within 500 pc, so it can never clip a star that the final
extinction-corrected box would have kept.

*What the Gaia-hosted 2MASS mirror lacks.*  ``gaiadr1.tmass_original_valid`` is
a reduced-column mirror: it carries ``ph_qual`` but not ``cc_flg``, ``rd_flg`` or
``bl_flg``.  Contamination/confusion screening therefore relies on ph_qual plus
the Gaia-side cross-match uniqueness flags here, and the full 2MASS PSC flags
are pulled from VizieR only for the handful of step-7 candidates.
"""

from __future__ import annotations

from . import config as cfg

# --------------------------------------------------------------------------
# Column lists
# --------------------------------------------------------------------------

GAIA_COLS = [
    "source_id", "ra", "dec", "l", "b",
    "parallax", "parallax_error", "parallax_over_error",
    "pmra", "pmdec", "ruwe",
    "phot_g_mean_mag", "phot_g_mean_flux_over_error",
    "phot_bp_mean_mag", "phot_bp_mean_flux_over_error",
    "phot_rp_mean_mag", "phot_rp_mean_flux_over_error",
    "phot_bp_rp_excess_factor", "bp_rp",
    "phot_bp_n_obs", "phot_rp_n_obs",
    "ipd_frac_multi_peak", "ipd_frac_odd_win",
    "astrometric_excess_noise", "astrometric_excess_noise_sig",
    "astrometric_params_solved", "visibility_periods_used",
    "phot_variable_flag", "non_single_star", "duplicated_source",
    "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat",
    "in_qso_candidates", "in_galaxy_candidates",
    "classprob_dsc_combmod_star",
]

AP_COLS = [
    "teff_gspphot", "logg_gspphot", "mh_gspphot",
    "ag_gspphot", "azero_gspphot", "ebpminrp_gspphot",
    "teff_gspspec", "logg_gspspec", "mh_gspspec", "flags_gspspec",
]

TMASS_COLS = ["j_m", "j_msigcom", "h_m", "h_msigcom",
              "ks_m", "ks_msigcom", "ph_qual", "ext_key"]

WISE_COLS = ["w1mpro", "w1mpro_error", "w2mpro", "w2mpro_error",
             "w3mpro", "w3mpro_error", "w4mpro", "w4mpro_error",
             "cc_flags", "ext_flag", "var_flag"]


def _select_clause() -> str:
    parts = [f"  g.{c}" for c in GAIA_COLS]
    parts += [f"  ap.{c}" for c in AP_COLS]
    parts += ["  tm.designation AS tmass_designation"]
    parts += [f"  tm.{c} AS tmass_{c}" for c in TMASS_COLS]
    parts += [
        "  xm.angular_distance AS tmass_xm_dist",
        "  xm.number_of_neighbours AS tmass_xm_nnb",
        "  xm.number_of_mates AS tmass_xm_nmates",
    ]
    parts += ["  aw.designation AS allwise_designation"]
    parts += [f"  aw.{c} AS wise_{c}" for c in WISE_COLS]
    parts += ["  aw.ph_qual AS wise_ph_qual",
              "  wxm.angular_distance AS wise_xm_dist",
              "  wxm.number_of_neighbours AS wise_xm_nnb",
              "  wxm.number_of_mates AS wise_xm_nmates"]
    return ",\n".join(parts)


FROM_CLAUSE = """
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm
  ON xm.source_id = g.source_id
JOIN gaiadr3.tmass_psc_xsc_join AS xj
  ON xj.clean_tmass_psc_xsc_oid = xm.clean_tmass_psc_xsc_oid
JOIN gaiadr1.tmass_original_valid AS tm
  ON tm.designation = xj.original_psc_source_id
LEFT OUTER JOIN gaiadr3.astrophysical_parameters AS ap
  ON ap.source_id = g.source_id
LEFT OUTER JOIN gaiadr3.allwise_best_neighbour AS wxm
  ON wxm.source_id = g.source_id
LEFT OUTER JOIN gaiadr1.allwise_original_valid AS aw
  ON aw.allwise_oid = wxm.allwise_oid
"""


def _where_clause(distance_max_pc: float) -> str:
    c = cfg.CUTS
    plx_min = 1000.0 / distance_max_pc
    plx_max = 1000.0 / c.distance_min_pc

    # Widen the observed-colour and observed-M_Ks pre-boxes by the largest
    # reddening/extinction that can occur inside the distance limit, so the
    # server-side pre-filter is guaranteed not to clip anything the final
    # extinction-corrected box would keep.
    #   E(BP-RP) = (A_BP/A_V - A_RP/A_V) * A_V   with Wang & Chen (2019) ratios
    r = cfg.EXTINCTION_RATIOS_WANG_CHEN_2019
    a_v_max = 3.0 if distance_max_pc > 600 else 1.5      # generous envelope
    e_bp_rp_max = (r["BP"] - r["RP"]) * a_v_max
    a_ks_max = r["Ks"] * a_v_max

    bp_rp_lo = c.bp_rp0_min - 0.05
    bp_rp_hi = c.bp_rp0_max + e_bp_rp_max + 0.05
    mks_lo = c.m_ks_min - 0.3
    mks_hi = c.m_ks_max + a_ks_max + 0.3

    return f"""
WHERE g.source_id BETWEEN {{lo}} AND {{hi}}
  AND g.parallax_over_error > {c.parallax_over_error_min}
  AND g.parallax > {plx_min}
  AND g.parallax < {plx_max}
  AND g.ruwe < {c.ruwe_max}
  AND g.phot_g_mean_mag BETWEEN {c.g_mag_min} AND {c.g_mag_max}
  AND g.phot_g_mean_flux_over_error > {c.g_flux_over_error_min}
  AND g.phot_bp_mean_mag IS NOT NULL
  AND g.phot_rp_mean_mag IS NOT NULL
  AND g.bp_rp BETWEEN {bp_rp_lo:.3f} AND {bp_rp_hi:.3f}
  AND g.ipd_frac_multi_peak <= {c.ipd_frac_multi_peak_max}
  AND g.ipd_frac_odd_win <= {c.ipd_frac_odd_win_max}
  AND g.astrometric_excess_noise_sig < {c.astrometric_excess_noise_sig_max}
  AND g.visibility_periods_used >= {c.visibility_periods_used_min}
  AND g.duplicated_source = 'false'
  AND tm.ks_m IS NOT NULL
  AND tm.ks_msigcom < {c.ks_sigma_max}
  AND xm.number_of_mates = 0
  AND xm.number_of_neighbours = 1
  AND (tm.ks_m + 5 * LOG10(g.parallax) - 10) BETWEEN {mks_lo:.3f} AND {mks_hi:.3f}
"""


def partition_bounds(index: int) -> tuple[int, int]:
    """source_id range for HEALPix partition ``index`` at the configured level."""
    if not 0 <= index < cfg.N_PARTITIONS:
        raise ValueError(f"partition {index} out of range 0..{cfg.N_PARTITIONS - 1}")
    lo = index * cfg.SOURCE_ID_DIVISOR
    hi = (index + 1) * cfg.SOURCE_ID_DIVISOR - 1
    return lo, hi


def chunk_query(index: int, *, distance_max_pc: float | None = None,
                top: int | None = None) -> str:
    """Full ADQL for one sky partition."""
    dmax = distance_max_pc if distance_max_pc is not None else cfg.CUTS.distance_max_pc
    lo, hi = partition_bounds(index)
    top_clause = f"TOP {top} " if top else ""
    where = _where_clause(dmax).format(lo=lo, hi=hi)
    return f"SELECT {top_clause}\n{_select_clause()}\n{FROM_CLAUSE}{where}"


def count_query(index: int, *, distance_max_pc: float | None = None) -> str:
    dmax = distance_max_pc if distance_max_pc is not None else cfg.CUTS.distance_max_pc
    lo, hi = partition_bounds(index)
    where = _where_clause(dmax).format(lo=lo, hi=hi)
    return f"SELECT COUNT(*) AS n\n{FROM_CLAUSE}{where}"
