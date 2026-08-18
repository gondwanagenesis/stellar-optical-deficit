#!/usr/bin/env python
"""Probe the gaiadr3.nss_two_body_orbit schema before building Search E.

    run.sh scripts/51_probe_nss_columns.py

The NSS tables carry solution-type-dependent columns, and selecting a column
that does not exist fails the whole ADQL job. Enumerate first, query second.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.tap import tap_client


def main() -> int:
    tap = tap_client()

    q = ("SELECT column_name, datatype, unit FROM TAP_SCHEMA.columns "
         "WHERE table_name = 'gaiadr3.nss_two_body_orbit'")
    rows = tap.launch_job(q).get_results()
    print(f"gaiadr3.nss_two_body_orbit: {len(rows)} columns\n")
    for r in rows:
        unit = str(r["unit"]) if r["unit"] is not None else ""
        print(f"  {str(r['column_name']):34s} {str(r['datatype']):10s} {unit}")

    print("\n--- nss_solution_type values and counts (parallax > 2 mas) ---")
    q2 = ("SELECT nss.nss_solution_type, COUNT(*) AS n "
          "FROM gaiadr3.nss_two_body_orbit AS nss "
          "JOIN gaiadr3.gaia_source AS g ON g.source_id = nss.source_id "
          "WHERE g.parallax > 2 "
          "GROUP BY nss.nss_solution_type ORDER BY n DESC")
    rows2 = tap.launch_job(q2).get_results()
    for r in rows2:
        print(f"  {str(r['nss_solution_type']):34s} {int(r['n']):8d}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
