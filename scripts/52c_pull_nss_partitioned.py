#!/usr/bin/env python
"""Pull Gaia DR3 astrometric orbits in parallax partitions, with visible progress.

    run.sh scripts/52c_pull_nss_partitioned.py

WHY A THIRD ATTEMPT
-------------------
Two earlier versions hung. The first joined nss_two_body_orbit to gaia_source
and binary_masses server-side and ran 98 minutes without returning. The second
dropped the joins but still asked for all 72,334 matching rows in one request,
and blocked for over two hours inside a single blocking call with no output.

A COUNT(*) over the same table with the same WHERE clause returns in 144 s, so
the archive is reachable and the predicate is cheap. The failure is in
retrieving one large VOTable over a loaded service, not in evaluating the
query.

This is the pattern the rest of the project already uses -- the main sample is
191 partition files -- so it is applied here too:

  * split by parallax into bands of roughly equal population;
  * submit each as a true background job and POLL, so a stall is visible
    rather than silent;
  * cache each partition separately, so an interrupted run resumes instead of
    restarting;
  * merge locally with binary_masses.parquet, which already carries m1 and the
    photometry.
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

# The ESA archive was the actual fault, not the query: during this run even
# "SELECT TOP 3" did not return inside 120 s, while the same statement against
# the ARI Heidelberg mirror returned in 1.1 s. Mirrors are tried in order and
# the first responsive one is used.
MIRRORS = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
    ("ESA", cfg.TAP_URL),
]

_CLIENT = None
_CLIENT_NAME = None


def get_client():
    """First mirror that answers a trivial query wins, and is then reused."""
    global _CLIENT, _CLIENT_NAME
    if _CLIENT is not None:
        return _CLIENT

    from astroquery.utils.tap.core import TapPlus
    for name, url in MIRRORS:
        t0 = time.time()
        try:
            tap = TapPlus(url=url)
            tap.launch_job("SELECT TOP 1 source_id "
                           "FROM gaiadr3.nss_two_body_orbit").get_results()
            log.info("using %s (%s) -- probe %.1fs", name, url, time.time() - t0)
            _CLIENT, _CLIENT_NAME = tap, name
            return tap
        except Exception as exc:
            log.warning("%s unusable (%.1fs): %s", name, time.time() - t0,
                        str(exc)[:120])
    raise RuntimeError("no Gaia TAP mirror responded")

COLUMNS = """
  source_id, nss_solution_type,
  a_thiele_innes, a_thiele_innes_error,
  b_thiele_innes, b_thiele_innes_error,
  f_thiele_innes, f_thiele_innes_error,
  g_thiele_innes, g_thiele_innes_error,
  period, period_error,
  eccentricity, eccentricity_error,
  parallax, parallax_error,
  goodness_of_fit, significance, bit_index
"""

SOLUTION_TYPES = ("'Orbital', 'AstroSpectroSB1', "
                  "'OrbitalTargetedSearch', 'OrbitalTargetedSearchValidated'")

# Parallax bands chosen so no partition holds more than ~12k rows.
BANDS = [(2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, 5.0),
         (5.0, 6.5), (6.5, 9.0), (9.0, 14.0), (14.0, 1e6)]

POLL_SECONDS = 10
MAX_WAIT_SECONDS = 900


def fetch_band(lo: float, hi: float, idx: int, force: bool) -> pd.DataFrame:
    part = cfg.RAW_DIR / f"nss_part_{idx:02d}.parquet"
    if part.exists() and not force:
        d = pd.read_parquet(part)
        log.info("  [%02d] cached  plx %.1f-%.1f : %5d rows", idx, lo, hi, len(d))
        return d

    q = (f"SELECT {COLUMNS} FROM gaiadr3.nss_two_body_orbit "
         f"WHERE nss_solution_type IN ({SOLUTION_TYPES}) "
         f"AND parallax >= {lo} AND parallax < {hi} "
         f"AND period IS NOT NULL AND a_thiele_innes IS NOT NULL")

    tap = get_client()
    t0 = time.time()
    job = tap.launch_job_async(q, dump_to_file=False, background=True)

    # Poll rather than block, so a stalled job is visible and bounded.
    while True:
        phase = str(job.get_phase()).upper()
        waited = time.time() - t0
        if phase in ("COMPLETED", "ERROR", "ABORTED"):
            break
        if waited > MAX_WAIT_SECONDS:
            log.error("  [%02d] plx %.1f-%.1f still %s after %.0fs -- abandoning",
                      idx, lo, hi, phase, waited)
            raise TimeoutError(f"partition {idx} stalled in {phase}")
        log.info("  [%02d] plx %.1f-%.1f : %s (%.0fs)", idx, lo, hi, phase, waited)
        time.sleep(POLL_SECONDS)

    if phase != "COMPLETED":
        raise RuntimeError(f"partition {idx} ended in phase {phase}")

    d = job.get_results().to_pandas()
    log.info("  [%02d] plx %.1f-%.1f : %5d rows in %5.0fs",
             idx, lo, hi, len(d), time.time() - t0)
    d.to_parquet(part, index=False)
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)

    frames = []
    for i, (lo, hi) in enumerate(BANDS):
        try:
            frames.append(fetch_band(lo, hi, i, args.force))
        except Exception as exc:
            log.error("partition %d failed: %s", i, exc)

    if not frames:
        log.error("no partitions succeeded")
        return 1

    nss = pd.concat(frames, ignore_index=True).drop_duplicates("source_id")
    log.info("total orbital solutions within 500 pc: %d", len(nss))
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

    out = cfg.RAW_DIR / "nss_orbits_500pc.parquet"
    merged.to_parquet(out, index=False)
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
