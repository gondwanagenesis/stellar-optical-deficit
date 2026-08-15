#!/usr/bin/env python
"""Full resumable pull of the Gaia x 2MASS x AllWISE sample.

Pulls the *wide* (1250 pc) superset once; the 500 pc primary sample is a
client-side subset of it.  The wide sample is only 3.2x larger than the primary
(measured, scripts/01b_count_wide.py) so a single download serves both the
science sample and the distance trade study.

Resumability: each HEALPix partition is an independent cached Parquet chunk.
Re-running skips anything already in the manifest with a matching query hash.
Kill and restart freely.

    ~/sd-venv/bin/python scripts/02_pull_sample.py --workers 3
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config as cfg
from pipeline.adql import chunk_query
from pipeline.tap import is_cached, run_adql

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("pull")


def chunk_name(idx: int, dmax: float) -> str:
    return f"sample_d{int(dmax)}_p{idx:03d}"


def pull_one(idx: int, dmax: float) -> tuple[int, int, bool]:
    """Returns (partition, n_rows, was_cached)."""
    q = chunk_query(idx, distance_max_pc=dmax)
    name = chunk_name(idx, dmax)
    cached = is_cached(name, q)
    df = run_adql(q, name=name)
    return idx, len(df), cached


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3,
                    help="concurrent TAP jobs; the ESA archive tolerates a "
                         "handful for anonymous users, do not be greedy")
    ap.add_argument("--distance-max-pc", type=float,
                    default=cfg.CUTS.distance_max_pc_wide)
    ap.add_argument("--partitions", type=str, default=None,
                    help="comma list or lo-hi range; default all")
    args = ap.parse_args()

    if args.partitions:
        if "-" in args.partitions and "," not in args.partitions:
            a, b = args.partitions.split("-")
            parts = list(range(int(a), int(b) + 1))
        else:
            parts = [int(x) for x in args.partitions.split(",")]
    else:
        parts = list(range(cfg.N_PARTITIONS))

    dmax = args.distance_max_pc
    todo = [p for p in parts if not is_cached(chunk_name(p, dmax),
                                              chunk_query(p, distance_max_pc=dmax))]
    log.info("partitions requested=%d already cached=%d to download=%d",
             len(parts), len(parts) - len(todo), len(todo))
    if not todo:
        log.info("nothing to do")
        return 0

    t0 = time.time()
    total_rows = 0
    done = 0
    failures: list[tuple[int, str]] = []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(pull_one, p, dmax): p for p in todo}
        for fut in as_completed(futures):
            p = futures[fut]
            try:
                idx, n, cached = fut.result()
                total_rows += n
                done += 1
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed else 0
                eta = (len(todo) - done) / rate if rate else float("nan")
                log.info("p%03d %7d rows | %3d/%3d done | %6d rows total | "
                         "ETA %5.1f min", idx, n, done, len(todo), total_rows,
                         eta / 60)
            except Exception as exc:                      # noqa: BLE001
                done += 1
                failures.append((p, str(exc)))
                log.error("p%03d FAILED: %s", p, exc)

    log.info("finished in %.1f min; %d rows; %d failures",
             (time.time() - t0) / 60, total_rows, len(failures))
    for p, msg in failures:
        log.error("  failed partition %d: %s", p, msg)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
