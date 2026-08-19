#!/usr/bin/env python
"""Pull the high-RUWE population by indexed source_id ranges.

    run.sh scripts/80_pull_high_ruwe_fast.py

WHY THIS REPLACES THE KEYSET PAGINATION IN SCRIPT 65
----------------------------------------------------
Script 65 paginated with

    WHERE ruwe >= 1.4 AND parallax > 2 AND ... AND source_id > X
    ORDER BY source_id ASC

which is the pattern that worked well for the 72k-row NSS table. On
gaia_source it is catastrophic. Timed directly against the mirror:

    filter + ORDER BY, one page   >300 s   (timed out)
    source_id BETWEEN range scan    16.8 s  -> 746 rows

because source_id is the PRIMARY KEY, so a BETWEEN clause is an indexed range
scan, whereas sorting a filtered slice of 1.8 billion rows is not indexed at
all. The original job ran for 2 h 47 m without completing a single logged page;
at >300 s per page and 564 pages it would have needed roughly 47 hours.

SAMPLING RATHER THAN COMPLETENESS, AND WHY THAT IS FINE HERE
------------------------------------------------------------
Even with range scans, covering the whole identifier space costs about
1000 x 17 s = 4.7 hours. Search N does not need every star: it compares the
dim-to-bright ratio of astrometrically disturbed stars against a low-RUWE
reference, and a few hundred thousand objects saturate that comparison long
before the systematics do.

So this takes every Nth chunk of the identifier space rather than a
contiguous block. Gaia source_id encodes HEALPix level 12, so a contiguous
range is a contiguous patch of sky and would bias the sample toward one
Galactic region; an interleaved sample covers the whole sky uniformly at the
same cost. The result is a uniform random sky sample, not a magnitude- or
position-limited subset, which is exactly what the ratio test wants.

Epoch propagation for the 2MASS cross-match is inherited from script 65 and is
not optional: this sample is nearby and fast-moving, and matching at the Gaia
epoch would both lose the highest-proper-motion stars and corrupt the
tmass_xm_dist cleanliness cut that Search N uses as its central discriminant.
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
log = logging.getLogger("pull_ruwe_fast")

MIRRORS = [("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
           ("ESA", "https://gea.esac.esa.int/tap-server/tap")]

SOURCE_ID_MAX = 6917529027641081856     # 2^62.6, the top of the DR3 id space
N_CHUNKS = 1000                          # granularity of the id grid
TAKE_EVERY = 5                           # 1 chunk in 5 -> ~20% sky, uniform

COLUMNS = ("source_id, ra, dec, l, b, parallax, parallax_error, "
           "parallax_over_error, pmra, pmdec, ruwe, astrometric_excess_noise, "
           "astrometric_excess_noise_sig, ipd_frac_multi_peak, "
           "duplicated_source, phot_g_mean_mag, phot_bp_mean_mag, "
           "phot_rp_mean_mag, bp_rp, phot_bp_rp_excess_factor, "
           "phot_variable_flag, non_single_star, teff_gspphot, mh_gspphot")


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
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--take-every", type=int, default=TAKE_EVERY)
    ap.add_argument("--ruwe-min", type=float, default=1.4)
    args = ap.parse_args()

    out = cfg.RAW_DIR / "high_ruwe_500pc.parquet"
    if out.exists() and not args.force:
        log.info("cached: %s", out)
        return 0

    tap = pick_mirror()
    edges = np.linspace(0, SOURCE_ID_MAX, N_CHUNKS + 1, dtype=np.int64)
    picks = list(range(0, N_CHUNKS, args.take_every))
    log.info("sampling %d of %d id chunks (~%.0f%% of sky, uniform)",
             len(picks), N_CHUNKS, 100 * len(picks) / N_CHUNKS)

    base = (f"ruwe >= {args.ruwe_min} AND parallax > 2 "
            f"AND parallax_over_error > 20")
    frames, t_start = [], time.time()
    for n, c in enumerate(picks):
        lo, hi = int(edges[c]), int(edges[c + 1])
        q = (f"SELECT {COLUMNS} FROM gaiadr3.gaia_source "
             f"WHERE source_id BETWEEN {lo} AND {hi} AND {base}")
        try:
            df = tap.launch_job(q).get_results().to_pandas()
        except Exception as exc:
            log.warning("  chunk %d failed: %s", c, str(exc)[:90])
            continue
        if len(df):
            frames.append(df)
        if n % 20 == 0 or n == len(picks) - 1:
            tot = sum(len(f) for f in frames)
            log.info("  chunk %4d/%4d  rows=%5d  total=%7d  elapsed %.0fs",
                     n + 1, len(picks), len(df), tot, time.time() - t_start)
        if len(df) >= 2000:
            log.warning("  chunk %d hit the 2000-row sync cap: "
                        "the id grid is too coarse here and this chunk is "
                        "truncated", c)

    if not frames:
        log.error("nothing retrieved")
        return 1

    d = pd.concat(frames, ignore_index=True).drop_duplicates("source_id")
    log.info("high-RUWE sources: %d in %.0fs", len(d), time.time() - t_start)
    log.info("  ruwe median %.2f, 90th pct %.2f, max %.1f",
             float(d["ruwe"].median()), float(d["ruwe"].quantile(0.9)),
             float(d["ruwe"].max()))

    # ---- 2MASS via CDS XMatch, at the 2MASS epoch ------------------------
    from astroquery.xmatch import XMatch
    from astropy.table import Table
    import astropy.units as u

    dt = 1999.5 - 2016.0
    dec_r = np.radians(d["dec"].to_numpy(float))
    ra99 = d["ra"].to_numpy(float) + (
        d["pmra"].fillna(0).to_numpy(float)
        / np.maximum(np.cos(dec_r), 1e-6)) * dt / 3.6e6
    de99 = d["dec"].to_numpy(float) + d["pmdec"].fillna(0).to_numpy(float) * dt / 3.6e6
    shift = np.hypot((ra99 - d["ra"].to_numpy(float)) * np.cos(dec_r),
                     de99 - d["dec"].to_numpy(float)) * 3600.0
    log.info("epoch propagation to 1999.5: median %.2f\", 99th pct %.2f\", "
             "max %.1f\"", float(np.median(shift)),
             float(np.percentile(shift, 99)), float(shift.max()))

    src = pd.DataFrame({"source_id": d["source_id"].to_numpy(),
                        "ra": ra99, "dec": de99})
    chunks, step = [], 25000
    for i0 in range(0, len(src), step):
        part = src.iloc[i0:i0 + step]
        t0 = time.time()
        r = XMatch.query(cat1=Table.from_pandas(part),
                         cat2="vizier:II/246/out", max_distance=3 * u.arcsec,
                         colRA1="ra", colDec1="dec").to_pandas()
        log.info("  xmatch %d-%d -> %d (%.0fs)", i0,
                 min(i0 + step, len(src)), len(r), time.time() - t0)
        chunks.append(r)

    xm = pd.concat(chunks, ignore_index=True)
    keep = {"source_id": "source_id", "angDist": "tmass_xm_dist",
            "Kmag": "tmass_ks_m", "e_Kmag": "tmass_ks_msigcom",
            "Jmag": "tmass_j_m", "Qflg": "tmass_ph_qual"}
    keep = {k: v for k, v in keep.items() if k in xm.columns}
    xm = xm[list(keep)].rename(columns=keep)
    nnb = xm.groupby("source_id").size().rename("tmass_xm_nnb")
    xm = xm.sort_values("tmass_xm_dist").drop_duplicates("source_id")
    xm = xm.merge(nnb, on="source_id", how="left")

    d = d.merge(xm, on="source_id", how="left")
    n_ks = int(d["tmass_ks_m"].notna().sum())
    log.info("with 2MASS Ks: %d of %d (%.1f%%)", n_ks, len(d),
             100 * n_ks / max(len(d), 1))

    d.to_parquet(out, index=False)
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
