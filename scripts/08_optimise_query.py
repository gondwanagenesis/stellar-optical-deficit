#!/usr/bin/env python
"""Find out what actually costs time in the production query.

A count over the identical FROM/WHERE returns in ~30 s, but retrieving 21k rows
of 74 columns takes ~250 s of serial-equivalent time.  So the cost is in result
materialisation and serialisation, not in the scan.  This script times variants
serially on one fixed partition to find the expensive piece.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.adql import (FROM_CLAUSE, GAIA_COLS, _select_clause,
                           _where_clause, partition_bounds)
from pipeline.tap import tap_client

PART = 11
DMAX = 1250.0


def run(label: str, select: str, from_clause: str) -> None:
    lo, hi = partition_bounds(PART)
    where = _where_clause(DMAX).format(lo=lo, hi=hi)
    q = f"SELECT\n{select}\n{from_clause}{where}"
    t0 = time.time()
    try:
        job = tap_client().launch_job_async(q, dump_to_file=False)
        r = job.get_results()
        print(f"{label:44s} {time.time()-t0:7.1f}s  rows={len(r):7d} cols={len(r.colnames)}")
    except Exception as exc:                            # noqa: BLE001
        print(f"{label:44s} {time.time()-t0:7.1f}s  FAILED {type(exc).__name__}: "
              f"{str(exc)[:70]}")


# Trimmed FROM: drop astrophysical_parameters (a 226-column, ~470M-row table)
FROM_NO_AP = """
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON xm.source_id = g.source_id
JOIN gaiadr3.tmass_psc_xsc_join AS xj ON xj.clean_tmass_psc_xsc_oid = xm.clean_tmass_psc_xsc_oid
JOIN gaiadr1.tmass_original_valid AS tm ON tm.designation = xj.original_psc_source_id
LEFT OUTER JOIN gaiadr3.allwise_best_neighbour AS wxm ON wxm.source_id = g.source_id
LEFT OUTER JOIN gaiadr1.allwise_original_valid AS aw ON aw.allwise_oid = wxm.allwise_oid
"""

FROM_NO_AP_NO_WISE = """
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON xm.source_id = g.source_id
JOIN gaiadr3.tmass_psc_xsc_join AS xj ON xj.clean_tmass_psc_xsc_oid = xm.clean_tmass_psc_xsc_oid
JOIN gaiadr1.tmass_original_valid AS tm ON tm.designation = xj.original_psc_source_id
"""

MINIMAL_SELECT = ",\n".join(
    [f"  g.{c}" for c in ["source_id", "l", "b", "parallax", "parallax_over_error",
                          "phot_g_mean_mag", "phot_g_mean_flux_over_error",
                          "bp_rp", "phot_bp_rp_excess_factor", "ruwe"]]
    + ["  tm.ks_m", "  tm.ks_msigcom", "  tm.ph_qual"])

GAIA_ONLY_SELECT = ",\n".join(f"  g.{c}" for c in GAIA_COLS)

if __name__ == "__main__":
    print(f"partition {PART}, d < {DMAX:.0f} pc, serial timings\n")
    run("A. full production query (74 cols)", _select_clause(), FROM_CLAUSE)
    run("B. no astrophysical_parameters join", GAIA_ONLY_SELECT + ",\n"
        "  tm.designation AS tmass_designation,\n  tm.ks_m,\n  tm.ks_msigcom,\n"
        "  tm.ph_qual,\n  tm.j_m,\n  tm.h_m,\n"
        "  aw.w1mpro,\n  aw.w2mpro,\n  aw.w3mpro,\n  aw.w4mpro,\n  aw.cc_flags,\n"
        "  aw.ext_flag", FROM_NO_AP)
    run("C. no AP, no WISE", GAIA_ONLY_SELECT + ",\n"
        "  tm.designation AS tmass_designation,\n  tm.ks_m,\n  tm.ks_msigcom,\n"
        "  tm.ph_qual,\n  tm.j_m,\n  tm.h_m", FROM_NO_AP_NO_WISE)
    run("D. minimal 13 columns, no AP/WISE", MINIMAL_SELECT, FROM_NO_AP_NO_WISE)
