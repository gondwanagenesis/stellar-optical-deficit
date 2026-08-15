#!/usr/bin/env python
"""Validate table names, column names and join logic against the live TAP
service before committing to a multi-hour download.

Run:  ~/sd-venv/bin/python scripts/00_probe_schema.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from pipeline import config as cfg
from pipeline.tap import run_adql, tap_client

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 60)

TABLES = [
    "gaiadr3.gaia_source",
    "gaiadr3.astrophysical_parameters",
    "gaiadr3.tmass_psc_xsc_best_neighbour",
    "gaiadr3.tmass_psc_xsc_join",
    "gaiadr1.tmass_original_valid",
    "gaiadr3.allwise_best_neighbour",
    "gaiadr1.allwise_original_valid",
]

WANT = {
    "gaiadr3.gaia_source": [
        "source_id", "ra", "dec", "l", "b", "parallax", "parallax_error",
        "parallax_over_error", "pmra", "pmdec", "ruwe", "phot_g_mean_mag",
        "phot_g_mean_flux_over_error", "phot_bp_mean_mag", "phot_rp_mean_mag",
        "phot_bp_rp_excess_factor", "bp_rp", "ipd_frac_multi_peak",
        "ipd_frac_odd_win", "astrometric_excess_noise_sig",
        "visibility_periods_used", "phot_variable_flag", "non_single_star",
        "duplicated_source", "nu_eff_used_in_astrometry", "pseudocolour",
        "ecl_lat", "astrometric_params_solved", "phot_bp_n_obs", "phot_rp_n_obs",
        "in_qso_candidates", "in_galaxy_candidates", "classprob_dsc_combmod_star",
    ],
    "gaiadr3.astrophysical_parameters": [
        "source_id", "teff_gspphot", "logg_gspphot", "mh_gspphot", "ag_gspphot",
        "azero_gspphot", "ebpminrp_gspphot", "distance_gspphot",
        "teff_gspspec", "mh_gspspec", "logg_gspspec", "flags_gspspec",
        "alphafe_gspspec", "mh_gspspec_upper", "mh_gspspec_lower",
    ],
    "gaiadr3.tmass_psc_xsc_best_neighbour": [
        "source_id", "clean_tmass_psc_xsc_oid", "original_ext_source_id",
        "angular_distance", "number_of_neighbours", "number_of_mates", "xm_flag",
    ],
    "gaiadr3.tmass_psc_xsc_join": [
        "clean_tmass_psc_xsc_oid", "original_psc_source_id",
        "original_xsc_source_id",
    ],
    "gaiadr1.tmass_original_valid": [
        "designation", "j_m", "j_msigcom", "h_m", "h_msigcom", "ks_m",
        "ks_msigcom", "ph_qual", "cc_flg", "rd_flg", "bl_flg", "ext_key",
    ],
    "gaiadr3.allwise_best_neighbour": [
        "source_id", "allwise_oid", "original_ext_source_id",
        "angular_distance", "number_of_neighbours", "number_of_mates",
    ],
    "gaiadr1.allwise_original_valid": [
        "allwise_oid", "designation", "w1mpro", "w1mpro_error", "w2mpro",
        "w2mpro_error", "w3mpro", "w3mpro_error", "w4mpro", "w4mpro_error",
        "cc_flags", "ext_flag", "var_flag", "ph_qual",
    ],
}


def probe_columns() -> dict[str, set[str]]:
    tap = tap_client()
    found: dict[str, set[str]] = {}
    for tbl in TABLES:
        q = (f"SELECT column_name, datatype FROM TAP_SCHEMA.columns "
             f"WHERE table_name = '{tbl}'")
        job = tap.launch_job(q)
        cols = set(str(c) for c in job.get_results()["column_name"])
        found[tbl] = cols
        print(f"\n=== {tbl}  ({len(cols)} columns)")
        missing = [c for c in WANT[tbl] if c not in cols]
        present = [c for c in WANT[tbl] if c in cols]
        print(f"    present : {len(present)}/{len(WANT[tbl])}")
        if missing:
            print(f"    MISSING : {missing}")
            near = {m: sorted(c for c in cols if m.split('_')[0] in c)[:8]
                    for m in missing}
            for m, cand in near.items():
                print(f"        {m!r} -> candidates {cand}")
    return found


JOIN_TEST = """
SELECT TOP 1000
  g.source_id, g.ra, g.dec, g.parallax, g.parallax_over_error, g.ruwe,
  g.phot_g_mean_mag, g.bp_rp,
  ap.teff_gspphot, ap.mh_gspphot, ap.ag_gspphot,
  tm.designation AS tmass_designation, tm.ks_m, tm.ks_msigcom, tm.ph_qual, tm.cc_flg,
  aw.designation AS allwise_designation, aw.w1mpro, aw.w3mpro, aw.w4mpro, aw.ext_flag
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
WHERE g.parallax_over_error > 20
  AND g.parallax > 2.0
  AND g.ruwe < 1.4
  AND g.source_id BETWEEN 4611686018427387904 AND 4620289320301428735
"""


def main() -> int:
    found = probe_columns()

    n_missing = sum(len([c for c in WANT[t] if c not in found[t]]) for t in TABLES)
    print(f"\n\n>>> total missing columns across all tables: {n_missing}")

    print("\n=== join + 1000-row validation ===")
    df = run_adql(JOIN_TEST, name="probe_join_1000", force=True)
    print(f"rows: {len(df)}   cols: {df.shape[1]}")
    print(df.dtypes)
    print(df.head(10).to_string())
    print("\nnull fractions:")
    print((df.isna().mean().sort_values(ascending=False) * 100).round(1).to_string())
    print("\nduplicate source_id rows:", int(df["source_id"].duplicated().sum()))
    return 0 if n_missing == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
