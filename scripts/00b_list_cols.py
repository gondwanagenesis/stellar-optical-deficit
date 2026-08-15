#!/usr/bin/env python
"""Dump the full column list of the external-catalogue tables."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline.tap import tap_client

for tbl in ["gaiadr1.tmass_original_valid", "gaiadr1.allwise_original_valid",
            "gaiadr3.tmass_psc_xsc_best_neighbour", "gaiadr3.allwise_best_neighbour"]:
    job = tap_client().launch_job(
        f"SELECT column_name, datatype, description FROM TAP_SCHEMA.columns "
        f"WHERE table_name = '{tbl}' ORDER BY column_name")
    r = job.get_results()
    print(f"\n===== {tbl}")
    for row in r:
        desc = str(row['description'])[:80].replace("\n", " ")
        print(f"  {str(row['column_name']):28s} {str(row['datatype']):10s} {desc}")
