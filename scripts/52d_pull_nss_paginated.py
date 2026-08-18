#!/usr/bin/env python
"""Pull Gaia DR3 astrometric orbits by keyset pagination over synchronous TAP.

    run.sh scripts/52d_pull_nss_paginated.py

THE DIAGNOSIS THAT LED HERE
---------------------------
Three earlier attempts hung, and the cause was never the query:

  * ESA's archive was returning nothing at all -- during this run even
    "SELECT TOP 3" did not complete inside 120 s.
  * The ARI Heidelberg mirror answered the same statement in 1.1 s.
  * But ASYNC jobs on the mirror sat in EXECUTING for 490 s and counting,
    while SYNC queries came back in 0.6 s.

So the async queue is backed up on both services while the synchronous
endpoint is healthy. Synchronous TAP caps the response at 2000 rows, which is
why one big request could never work.

The fix is keyset pagination: repeatedly ask for the next 2000 rows ordered by
source_id, carrying the last id forward as the lower bound. source_id is
unique, so pages cannot overlap or skip, which OFFSET-based paging cannot
guarantee when the underlying scan order varies. 72,334 rows is 37 pages at
roughly a second each.

MOD() was tried as a partition key first and silently failed to filter -- the
count came back as the full table -- so it is not used.
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("pull_nss")

MIRRORS = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
    ("ESA", cfg.TAP_URL),
]

COLUMNS = ("source_id, nss_solution_type, "
           "a_thiele_innes, a_thiele_innes_error, "
           "b_thiele_innes, b_thiele_innes_error, "
           "f_thiele_innes, f_thiele_innes_error, "
           "g_thiele_innes, g_thiele_innes_error, "
           "period, period_error, "
           "eccentricity, eccentricity_error, "
           "parallax, parallax_error, "
           "goodness_of_fit, significance, bit_index")

SOLUTION_TYPES = ("'Orbital', 'AstroSpectroSB1', "
                  "'OrbitalTargetedSearch', 'OrbitalTargetedSearchValidated'")

PAGE = 2000
MAX_PAGES = 200


def pick_mirror():
    from astroquery.utils.tap.core import TapPlus
    for name, url in MIRRORS:
        t0 = time.time()
        try:
            tap = TapPlus(url=url)
            tap.launch_job("SELECT TOP 1 source_id "
                           "FROM gaiadr3.nss_two_body_orbit").get_results()
            log.info("using %s -- probe %.1fs", name, time.time() - t0)
            return tap
        except Exception as exc:
            log.warning("%s unusable (%.1fs): %s", name,
                        time.time() - t0, str(exc)[:110])
    raise RuntimeError("no Gaia TAP mirror responded")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out = cfg.RAW_DIR / "nss_orbits_500pc.parquet"
    if out.exists() and not args.force:
        log.info("cached: %s", out)
        return 0

    tap = pick_mirror()
    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)

    pages = []
    last_id = -1
    t_start = time.time()

    for page in range(MAX_PAGES):
        q = (f"SELECT TOP {PAGE} {COLUMNS} "
             f"FROM gaiadr3.nss_two_body_orbit "
             f"WHERE nss_solution_type IN ({SOLUTION_TYPES}) "
             f"AND parallax > 2 "
             f"AND period IS NOT NULL "
             f"AND a_thiele_innes IS NOT NULL "
             f"AND source_id > {last_id} "
             f"ORDER BY source_id ASC")
        t0 = time.time()
        df = tap.launch_job(q).get_results().to_pandas()
        if len(df) == 0:
            break
        pages.append(df)
        last_id = int(df["source_id"].iloc[-1])
        total = sum(len(p) for p in pages)
        log.info("page %3d: %4d rows (%6d total) in %4.1fs",
                 page, len(df), total, time.time() - t0)
        if len(df) < PAGE:
            break
    else:
        log.warning("hit MAX_PAGES; result may be incomplete")

    nss = pd.concat(pages, ignore_index=True).drop_duplicates("source_id")
    log.info("pulled %d orbital solutions in %.0fs",
             len(nss), time.time() - t_start)
    for t, n in nss["nss_solution_type"].value_counts().items():
        log.info("  %-34s %7d", t, n)

    bm = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    bm_cols = [c for c in [
        "source_id", "m1", "m1_lower", "m1_upper", "m2", "fluxratio",
        "combination_method", "ra", "dec", "l", "b", "ruwe",
        "phot_g_mean_mag", "bp_rp", "tmass_ks_m"] if c in bm.columns]

    merged = nss.merge(bm[bm_cols], on="source_id", how="left",
                       suffixes=("", "_bm"))
    n_m1 = int(merged["m1"].notna().sum())
    log.info("merged: %d rows, %d with a primary mass (%.1f%%)",
             len(merged), n_m1, 100 * n_m1 / max(len(merged), 1))

    merged.to_parquet(out, index=False)
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
