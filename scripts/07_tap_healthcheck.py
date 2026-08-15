#!/usr/bin/env python
"""Time a trivial and a representative query to see whether the ESA archive is
degraded or whether our query is the problem."""
from __future__ import annotations
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.adql import chunk_query, count_query
from pipeline.tap import tap_client

def timed(label, q, sync=False):
    t0 = time.time()
    try:
        tap = tap_client()
        job = tap.launch_job(q) if sync else tap.launch_job_async(q, dump_to_file=False)
        r = job.get_results()
        print(f"{label:38s} {time.time()-t0:7.1f}s  rows={len(r)}")
    except Exception as exc:
        print(f"{label:38s} {time.time()-t0:7.1f}s  FAILED {type(exc).__name__}: {str(exc)[:80]}")

timed("sync: SELECT 1 row from gaia_source",
      "SELECT TOP 1 source_id FROM gaiadr3.gaia_source", sync=True)
timed("async: SELECT 1 row from gaia_source",
      "SELECT TOP 1 source_id FROM gaiadr3.gaia_source")
timed("async: count d500 partition 11", count_query(11, distance_max_pc=500))
timed("async: full d500 partition 11 TOP 5000",
      chunk_query(11, distance_max_pc=500, top=5000))
