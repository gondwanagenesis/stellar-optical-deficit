#!/usr/bin/env python
"""Add the four astrometry columns script 80 did not pull.

    run.sh scripts/81_pull_ruwe_astrometry_cols.py

WHY
---
Script 80 replaced script 65's broken pull and is correct in every way that
matters for speed and for the epoch propagation. But its column list dropped
four fields that the main pipeline needs, and Search N crashed on the first
of them:

    nu_eff_used_in_astrometry
    pseudocolour
    ecl_lat
    astrometric_params_solved

All four feed pipeline.sample.parallax_zero_point, the Lindegren et al. (2021)
correction. That is not optional bookkeeping. The correction is -20 to -40 uas,
which at parallax 2 mas is a 1-2% distance error and so 0.02-0.04 mag of
distance modulus -- an order of magnitude above the sensitivity this project is
trying to reach. Dropping it, or substituting zero, would silently change the
absolute magnitudes that every residual is built from.

It would have been possible to compute ecl_lat locally, since it is a pure
coordinate rotation of (ra, dec), but the other three exist only in the
archive, so all four are re-queried together.

HOW
---
The same identifier grid script 80 used -- 1000 equal chunks of the DR3 id
space, every 5th taken -- with the same WHERE clause, so the row set is
identical by construction. Only the narrow column list changes, which makes
this far lighter than the original pull. The result is merged onto the
existing parquet by source_id, and the script refuses to write if the join
does not cover every row.
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
log = logging.getLogger("pull_astrom")

MIRRORS = [("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
           ("ESA", "https://gea.esac.esa.int/tap-server/tap")]

# Identical to script 80. Any drift here silently changes the row set.
SOURCE_ID_MAX = 6917529027641081856
N_CHUNKS = 1000
TAKE_EVERY = 5

COLUMNS = ("source_id, nu_eff_used_in_astrometry, pseudocolour, "
           "ecl_lat, astrometric_params_solved")
NEW_COLS = ["nu_eff_used_in_astrometry", "pseudocolour",
            "ecl_lat", "astrometric_params_solved"]


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
    have = [c for c in NEW_COLS if c in d.columns]
    if len(have) == len(NEW_COLS):
        log.info("all four columns already present; nothing to do")
        return 0
    log.info("target has %d rows, missing %s", len(d),
             [c for c in NEW_COLS if c not in d.columns])

    tap = pick_mirror()
    edges = np.linspace(0, SOURCE_ID_MAX, N_CHUNKS + 1, dtype=np.int64)
    picks = list(range(0, N_CHUNKS, args.take_every))
    base = (f"ruwe >= {args.ruwe_min} AND parallax > 2 "
            f"AND parallax_over_error > 20")

    frames, t0, n_capped = [], time.time(), 0
    for n, c in enumerate(picks):
        lo, hi = int(edges[c]), int(edges[c + 1])
        q = (f"SELECT {COLUMNS} FROM gaiadr3.gaia_source "
             f"WHERE source_id BETWEEN {lo} AND {hi} AND {base}")
        for attempt in range(4):
            try:
                df = tap.launch_job(q).get_results().to_pandas()
                break
            except Exception as exc:
                if attempt == 3:
                    log.warning("  chunk %d failed permanently: %s",
                                c, str(exc)[:90])
                    df = None
                    break
                time.sleep(5 * (attempt + 1))
        if df is None:
            continue
        if len(df) >= 2000:
            n_capped += 1
            log.warning("  chunk %d hit the 2000-row sync cap", c)
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
    log.info("astrometry rows: %d in %.0fs (%d capped)",
             len(a), time.time() - t0, n_capped)

    merged = d.merge(a, on="source_id", how="left")
    if len(merged) != len(d):
        log.error("merge changed the row count %d -> %d; refusing to write",
                  len(d), len(merged))
        return 1

    cov = merged["nu_eff_used_in_astrometry"].notna().sum()
    solved_cov = merged["astrometric_params_solved"].notna().sum()
    log.info("join coverage: nu_eff %d/%d (%.2f%%), solved %d/%d (%.2f%%)",
             cov, len(merged), 100 * cov / len(merged),
             solved_cov, len(merged), 100 * solved_cov / len(merged))

    # nu_eff is null by design for 2-parameter solutions, so it is not a
    # coverage test. astrometric_params_solved is populated for every source
    # in the catalogue, so if that is missing the join itself failed.
    if solved_cov < len(merged):
        log.error("%d rows did not join at all; the chunk grid must differ "
                  "from script 80. Refusing to write.",
                  int(len(merged) - solved_cov))
        return 1

    merged.to_parquet(target, index=False)
    log.info("wrote %s with %d columns", target, merged.shape[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
