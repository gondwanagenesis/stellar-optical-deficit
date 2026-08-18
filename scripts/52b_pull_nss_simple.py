#!/usr/bin/env python
"""Pull Gaia DR3 astrometric orbits, no server-side joins.

    run.sh scripts/52b_pull_nss_simple.py

WHY THIS EXISTS
---------------
The obvious query joins nss_two_body_orbit to gaia_source and binary_masses to
collect orbits, photometry and primary masses in one shot. That query ran for
98 minutes on the ESA archive without returning: a three-table join against
gaia_source is the expensive part, and it is entirely avoidable.

nss_two_body_orbit carries its own parallax, so the 500 pc cut needs no join.
And binary_masses.parquet -- already on disk from script 16 -- carries m1 plus
the photometry and sky position for the same source_ids. So we pull the orbit
table alone and merge locally, which costs one pandas join instead of a
server-side scan of 1.8 billion rows.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from pipeline import config as cfg
from pipeline.tap import run_adql

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("pull_nss")

NSS_QUERY = """
SELECT
  source_id, nss_solution_type,
  a_thiele_innes, a_thiele_innes_error,
  b_thiele_innes, b_thiele_innes_error,
  f_thiele_innes, f_thiele_innes_error,
  g_thiele_innes, g_thiele_innes_error,
  period, period_error,
  eccentricity, eccentricity_error,
  inclination, inclination_error,
  parallax, parallax_error,
  goodness_of_fit, significance, efficiency,
  bit_index, flags,
  astrometric_n_good_obs_al,
  g_luminosity_ratio,
  mass_ratio, mass_ratio_error,
  semi_amplitude_primary
FROM gaiadr3.nss_two_body_orbit
WHERE nss_solution_type IN (
        'Orbital', 'AstroSpectroSB1',
        'OrbitalTargetedSearch', 'OrbitalTargetedSearchValidated')
  AND parallax > 2
  AND period IS NOT NULL
  AND a_thiele_innes IS NOT NULL
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    nss = run_adql(NSS_QUERY, name="nss_orbits_raw", force=args.force)
    log.info("orbital solutions within 500 pc: %d", len(nss))
    for t, n in nss["nss_solution_type"].value_counts().items():
        log.info("  %-34s %7d", t, n)

    bm = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    log.info("binary_masses on disk: %d", len(bm))

    bm_cols = [c for c in [
        "source_id", "m1", "m1_lower", "m1_upper",
        "m2", "m2_lower", "m2_upper", "fluxratio", "combination_method",
        "ra", "dec", "l", "b", "ruwe",
        "phot_g_mean_mag", "bp_rp",
        "tmass_ks_m", "tmass_ph_qual"] if c in bm.columns]

    merged = nss.merge(bm[bm_cols], on="source_id", how="left",
                       suffixes=("", "_bm"))
    n_m1 = int(merged["m1"].notna().sum())
    log.info("merged: %d rows, %d with a primary mass (%.1f%%)",
             len(merged), n_m1, 100 * n_m1 / max(len(merged), 1))

    out = cfg.RAW_DIR / "nss_orbits_500pc.parquet"
    merged.to_parquet(out, index=False)
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
