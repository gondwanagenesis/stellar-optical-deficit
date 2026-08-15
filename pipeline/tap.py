"""Resumable, cached TAP client for the ESA Gaia archive.

Rules this module enforces:
  * Nothing is ever downloaded twice.  A chunk is identified by (name, sha256 of
    the ADQL text).  If the Parquet file exists and the manifest records a
    matching query hash, the chunk is skipped.
  * Every download appends a manifest record with the full ADQL, the row count,
    the job id and a UTC timestamp, so the provenance of every row is
    recoverable after the fact.
  * Transient archive failures are retried with exponential backoff.  A chunk
    that returns exactly TAP_ROW_CAP rows is treated as a *failure*, not a
    success, because the archive truncates silently.
"""

from __future__ import annotations

import hashlib
import json
import logging
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from astroquery.utils.tap.core import TapPlus

from . import config as cfg

log = logging.getLogger(__name__)

# astroquery's TapPlus keeps per-instance connection state, so each worker
# thread gets its own client.
_local = threading.local()
_manifest_lock = threading.Lock()
_manifest_cache: dict[str, dict] | None = None


def tap_client() -> TapPlus:
    client = getattr(_local, "tap", None)
    if client is None:
        client = TapPlus(url=cfg.TAP_URL)
        _local.tap = client
    return client


def adql_hash(query: str) -> str:
    """Whitespace-insensitive hash so cosmetic reformatting doesn't bust cache."""
    normalised = " ".join(query.split())
    return hashlib.sha256(normalised.encode()).hexdigest()[:16]


# --------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------

def _read_manifest(refresh: bool = False) -> dict[str, dict]:
    """Manifest keyed by chunk name; later records win.

    Cached in memory because the resumability check reads it once per chunk and
    the file grows to thousands of lines during a full pull.
    """
    global _manifest_cache
    with _manifest_lock:
        if _manifest_cache is not None and not refresh:
            return _manifest_cache
        out: dict[str, dict] = {}
        if cfg.MANIFEST_PATH.exists():
            with cfg.MANIFEST_PATH.open() as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    rec = json.loads(line)
                    out[rec["name"]] = rec
        _manifest_cache = out
        return out


def _append_manifest(rec: dict) -> None:
    global _manifest_cache
    cfg.MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    with _manifest_lock:
        with cfg.MANIFEST_PATH.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        if _manifest_cache is not None:
            _manifest_cache[rec["name"]] = rec


# --------------------------------------------------------------------------
# Table -> DataFrame
# --------------------------------------------------------------------------

def _to_dataframe(table) -> pd.DataFrame:
    df = table.to_pandas()
    # VOTable char columns arrive as bytes under some astropy/pyarrow combos.
    for col in df.columns:
        if df[col].dtype == object:
            sample = df[col].dropna()
            if len(sample) and isinstance(sample.iloc[0], (bytes, bytearray)):
                df[col] = df[col].str.decode("utf-8", errors="replace")
    return df


# --------------------------------------------------------------------------
# Core
# --------------------------------------------------------------------------

def run_adql(query: str, *, name: str, cache_dir: Path | None = None,
             force: bool = False, expect_truncation_cap: int | None = None) -> pd.DataFrame:
    """Execute ``query``, caching the result to ``<cache_dir>/<name>.parquet``.

    Returns the DataFrame either from cache or from a fresh download.
    """
    cache_dir = cache_dir or cfg.RAW_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{name}.parquet"
    qhash = adql_hash(query)
    cap = expect_truncation_cap if expect_truncation_cap is not None else cfg.TAP_ROW_CAP

    if path.exists() and not force:
        manifest = _read_manifest()
        rec = manifest.get(name)
        if rec and rec.get("adql_hash") == qhash:
            log.info("cache hit  %-28s %8d rows", name, rec.get("n_rows", -1))
            return pd.read_parquet(path)
        log.warning("cache file %s exists but query hash changed; re-downloading", path.name)

    last_err: Exception | None = None
    for attempt in range(1, cfg.TAP_MAX_RETRIES + 1):
        try:
            t0 = time.time()
            job = tap_client().launch_job_async(query, dump_to_file=False,
                                                background=False)
            table = job.get_results()
            df = _to_dataframe(table)
            elapsed = time.time() - t0

            if len(df) >= cap:
                raise RuntimeError(
                    f"chunk {name} returned {len(df)} rows, at or above the "
                    f"archive row cap {cap}: the result is silently truncated. "
                    f"Subdivide this partition."
                )

            # Write-then-rename so an interrupted run never leaves a truncated
            # Parquet file that a later run would treat as a cache hit.
            tmp = path.with_suffix(".parquet.tmp")
            df.to_parquet(tmp, index=False, compression="zstd")
            tmp.replace(path)
            _append_manifest({
                "name": name,
                "adql_hash": qhash,
                "adql": " ".join(query.split()),
                "n_rows": int(len(df)),
                "n_cols": int(df.shape[1]),
                "job_id": getattr(job, "jobid", None),
                "path": str(path.relative_to(cfg.PROJECT_ROOT)),
                "utc": datetime.now(timezone.utc).isoformat(),
                "elapsed_s": round(elapsed, 1),
                "tap_url": cfg.TAP_URL,
                "attempt": attempt,
            })
            log.info("downloaded %-28s %8d rows in %6.1fs", name, len(df), elapsed)
            return df

        except Exception as exc:                       # noqa: BLE001
            last_err = exc
            if "silently truncated" in str(exc):
                raise                                   # not retryable
            sleep = cfg.TAP_RETRY_BASE_SLEEP * (2 ** (attempt - 1))
            log.warning("chunk %s attempt %d/%d failed: %s -- retrying in %.0fs",
                        name, attempt, cfg.TAP_MAX_RETRIES, exc, sleep)
            if attempt < cfg.TAP_MAX_RETRIES:
                time.sleep(sleep)

    raise RuntimeError(f"chunk {name} failed after {cfg.TAP_MAX_RETRIES} attempts") from last_err


def is_cached(name: str, query: str, cache_dir: Path | None = None) -> bool:
    cache_dir = cache_dir or cfg.RAW_DIR
    path = cache_dir / f"{name}.parquet"
    if not path.exists():
        return False
    rec = _read_manifest().get(name)
    return bool(rec and rec.get("adql_hash") == adql_hash(query))


def manifest_summary() -> pd.DataFrame:
    recs = list(_read_manifest().values())
    if not recs:
        return pd.DataFrame(columns=["name", "n_rows", "utc"])
    return pd.DataFrame(recs)[["name", "n_rows", "n_cols", "elapsed_s", "utc"]]
