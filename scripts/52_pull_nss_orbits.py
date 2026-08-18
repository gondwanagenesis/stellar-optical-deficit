#!/usr/bin/env python
"""Pull Gaia DR3 non-single-star orbital solutions within 500 pc.

    run.sh scripts/52_pull_nss_orbits.py

WHY
---
Every channel so far measured photons. This pulls the data for a channel that
measures GRAVITY: astrometric orbits, from which the companion mass follows
without the companion ever emitting anything.

An absorber that hides a star from every photometric channel -- the f -> 1
case our optical-deficit estimator is structurally blind to, because a fully
enshrouded star simply leaves the photometric sample -- still pulls on its
companion with undiminished mass. In an astrometric binary that is a maximal
signal, not a null one.

WHAT IS PULLED
--------------
Solution types that carry a photocentre orbit (Thiele-Innes coefficients):

    Orbital                          49,399 within 500 pc
    AstroSpectroSB1                  22,011
    OrbitalTargetedSearch(+Validated)   492

joined to gaia_source for photometry/astrometry and to binary_masses for the
primary mass m1 where Gaia derived one.

The Thiele-Innes coefficients (A, B, F, G, in mas) give the photocentre
semi-major axis, which is the quantity the Search E estimator needs.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.tap import run_adql

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("pull_nss")

NSS_QUERY = """
SELECT
  nss.source_id,
  nss.nss_solution_type,
  nss.a_thiele_innes, nss.a_thiele_innes_error,
  nss.b_thiele_innes, nss.b_thiele_innes_error,
  nss.f_thiele_innes, nss.f_thiele_innes_error,
  nss.g_thiele_innes, nss.g_thiele_innes_error,
  nss.period, nss.period_error,
  nss.eccentricity, nss.eccentricity_error,
  nss.inclination, nss.inclination_error,
  nss.parallax AS nss_parallax, nss.parallax_error AS nss_parallax_error,
  nss.goodness_of_fit, nss.significance, nss.efficiency,
  nss.bit_index, nss.flags,
  nss.astrometric_n_good_obs_al,
  nss.g_luminosity_ratio,
  nss.mass_ratio, nss.mass_ratio_error,
  nss.semi_amplitude_primary,
  bm.m1, bm.m1_lower, bm.m1_upper,
  bm.m2, bm.m2_lower, bm.m2_upper,
  bm.fluxratio, bm.combination_method,
  g.ra, g.dec, g.l, g.b,
  g.parallax, g.parallax_error, g.parallax_over_error,
  g.pmra, g.pmdec, g.radial_velocity,
  g.ruwe, g.astrometric_excess_noise,
  g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
  g.bp_rp, g.phot_bp_rp_excess_factor,
  g.teff_gspphot, g.mh_gspphot, g.logg_gspphot,
  g.non_single_star
FROM gaiadr3.nss_two_body_orbit AS nss
JOIN gaiadr3.gaia_source AS g
  ON g.source_id = nss.source_id
LEFT JOIN gaiadr3.binary_masses AS bm
  ON bm.source_id = nss.source_id
WHERE nss.nss_solution_type IN (
        'Orbital', 'AstroSpectroSB1',
        'OrbitalTargetedSearch', 'OrbitalTargetedSearchValidated')
  AND g.parallax > 2
  AND g.parallax_over_error > 20
  AND nss.period IS NOT NULL
  AND nss.a_thiele_innes IS NOT NULL
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    df = run_adql(NSS_QUERY, name="nss_orbits_500pc", force=args.force)
    log.info("pulled %d orbital solutions within 500 pc", len(df))

    log.info("by solution type:")
    for t, n in df["nss_solution_type"].value_counts().items():
        log.info("  %-34s %7d", t, n)

    n_m1 = int(df["m1"].notna().sum())
    log.info("with a Gaia-derived primary mass m1: %d (%.1f%%)",
             n_m1, 100 * n_m1 / max(len(df), 1))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
