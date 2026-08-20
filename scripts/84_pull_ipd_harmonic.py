#!/usr/bin/env python
"""Pull the IPD harmonic columns onto the high-RUWE population.

    run.sh scripts/84_pull_ipd_harmonic.py

WHY
---
Gaia fits a PSF to every windowed transit and records how badly that fit
degrades as a function of the position angle of the scan direction. For a
single point source the goodness-of-fit does not care which way the satellite
crossed it, so the variation is noise. For a source with unresolved structure
below the PSF -- a companion, a blend, an extended envelope -- the fit is worst
when the scan runs along the elongation and best when it runs across it, so the
GoF picks up a sinusoid at twice the scan angle. `ipd_gof_harmonic_amplitude`
is that sinusoid's amplitude and `ipd_gof_harmonic_phase` is its phase, which
is the position angle of the structure modulo 180 degrees.

That phase is the only direct handle in DR3 on the *direction* of sub-PSF
structure, and no channel in this project has used it. Neither, as far as we
can tell, has any published technosignature search. It is what channel 22
needs.

HOW
---
The same identifier grid scripts 80 and 81 used -- 1000 equal chunks of the DR3
id space, every 5th taken, indexed BETWEEN range scans rather than keyset
pagination -- with the same WHERE clause, so the row set is identical by
construction and the merge is total. Only the column list changes.

Do not "improve" the grid. Scripts 80, 81 and 84 must agree exactly or the
join silently drops rows.
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("pull_ipd")

MIRRORS = [("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
           ("ESA", "https://gea.esac.esa.int/tap-server/tap")]

# Identical to scripts 80 and 81. Any drift here silently changes the row set.
SOURCE_ID_MAX = 6917529027641081856
N_CHUNKS = 1000
TAKE_EVERY = 5

NEW_COLS = ["ipd_gof_harmonic_amplitude", "ipd_gof_harmonic_phase",
            "ipd_frac_odd_win", "matched_transits",
            "visibility_periods_used", "astrometric_n_good_obs_al"]
COLUMNS = "source_id, " + ", ".join(NEW_COLS)


def pick_mirror():
    from astroquery.utils.tap.core import TapPlus
    for name, url in MIRRORS:
        try:
            tap = TapPlus(url=url)
            tap.launch_job("SELECT TOP 1 source_id "
                           "FROM gaiadr3.gaia_source").get_results()
            log.info("using %s", name)
            return tap
        except Exception as exc:
            log.warning("%s unusable: %s", name, str(exc)[:100])
    raise RuntimeError("no mirror responded")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--take-every", type=int, default=TAKE_EVERY)
    ap.add_argument("--ruwe-min", type=float, default=1.4)
    args = ap.parse_args()

    target = cfg.RAW_DIR / "high_ruwe_500pc.parquet"
    if not target.exists():
        log.error("missing %s -- run scripts/80_pull_high_ruwe_fast.py first",
                  target)
        return 1

    d = pd.read_parquet(target)
    if all(c in d.columns for c in NEW_COLS):
        log.info("all IPD columns already present; nothing to do")
        return 0
    log.info("target has %d rows, missing %s", len(d),
             [c for c in NEW_COLS if c not in d.columns])

    tap = pick_mirror()
    edges = np.linspace(0, SOURCE_ID_MAX, N_CHUNKS + 1, dtype=np.int64)
    picks = list(range(0, N_CHUNKS, args.take_every))
    base = (f"ruwe >= {args.ruwe_min} AND parallax > 2 "
            f"AND parallax_over_error > 20")

    # Each chunk is cached on arrival. The whole-run write happens only at the
    # end and refuses to proceed on an incomplete join, so without this a
    # single permanently-failed chunk near the end discards ~50 minutes of
    # completed queries. A rerun now replays the cache in seconds.
    cache = cfg.RAW_DIR / "ipd_chunks"
    cache.mkdir(exist_ok=True)

    frames, t0, n_capped, n_failed, n_cached = [], time.time(), 0, 0, 0
    for n, c in enumerate(picks):
        lo, hi = int(edges[c]), int(edges[c + 1])
        cf = cache / f"chunk_{c:04d}.parquet"
        if cf.exists():
            frames.append(pd.read_parquet(cf))
            n_cached += 1
            continue
        q = (f"SELECT {COLUMNS} FROM gaiadr3.gaia_source "
             f"WHERE source_id BETWEEN {lo} AND {hi} AND {base}")
        df = None
        for attempt in range(4):
            try:
                df = tap.launch_job(q).get_results().to_pandas()
                break
            except Exception as exc:
                if attempt == 3:
                    log.warning("  chunk %d failed permanently: %s",
                                c, str(exc)[:90])
                    n_failed += 1
                    break
                time.sleep(5 * (attempt + 1))
        if df is None:
            continue
        if len(df) >= 2000:
            # The sync endpoint truncates at 2000 without saying so. A capped
            # chunk is a silently incomplete row set, not a small one.
            n_capped += 1
            log.warning("  chunk %d hit the 2000-row sync cap", c)
        # Written even when empty: an empty chunk is a real, reusable answer,
        # and caching it stops a rerun from re-querying a genuinely empty
        # slice of the id space.
        df.to_parquet(cf, index=False)
        if len(df):
            frames.append(df)
        if n % 20 == 0 or n == len(picks) - 1:
            tot = sum(len(f) for f in frames)
            log.info("  chunk %4d/%4d  rows=%5d  total=%7d  elapsed %.0fs",
                     n + 1, len(picks), len(df), tot, time.time() - t0)

    if not frames:
        log.error("nothing retrieved")
        return 1

    a = pd.concat(frames, ignore_index=True).drop_duplicates("source_id")
    log.info("IPD rows: %d in %.0fs (%d capped, %d failed, %d from cache)",
             len(a), time.time() - t0, n_capped, n_failed, n_cached)

    merged = d.merge(a, on="source_id", how="left")
    if len(merged) != len(d):
        log.error("merge changed the row count %d -> %d; refusing to write",
                  len(d), len(merged))
        return 1

    # matched_transits is populated for every source in the catalogue, so it
    # is the join test. The harmonic amplitude and phase are legitimately null
    # for sources whose IPD never ran, so their coverage is a property of the
    # data and is only reported.
    joined = merged["matched_transits"].notna().sum()
    amp_cov = merged["ipd_gof_harmonic_amplitude"].notna().sum()
    pha_cov = merged["ipd_gof_harmonic_phase"].notna().sum()
    log.info("join coverage: matched_transits %d/%d (%.2f%%)",
             joined, len(merged), 100 * joined / len(merged))
    log.info("data coverage: amplitude %.2f%%, phase %.2f%%",
             100 * amp_cov / len(merged), 100 * pha_cov / len(merged))

    if joined < len(merged):
        log.error("%d rows did not join at all; the chunk grid must differ "
                  "from script 80. Refusing to write.",
                  int(len(merged) - joined))
        return 1

    merged.to_parquet(target, index=False)
    log.info("wrote %s with %d columns", target, merged.shape[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
