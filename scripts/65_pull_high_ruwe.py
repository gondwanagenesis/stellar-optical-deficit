#!/usr/bin/env python
"""Pull the high-RUWE population that the sample build never downloaded.

SUPERSEDED -- USE scripts/80_pull_high_ruwe_fast.py INSTEAD.
============================================================
This script does not work and will refuse to run. It is kept because the
reasoning in the notes below is still the reasoning behind Search N, and
because the way it failed is worth not repeating.

Two independent defects, both fixed in script 80:

1. KEYSET PAGINATION IS CATASTROPHIC ON gaia_source. The
   "WHERE ... AND source_id > X ORDER BY source_id ASC" pattern below was
   lifted from script 52d, where it works well -- but that table has 72k
   rows. gaia_source has 1.8e9, and an open-ended lower bound gives the
   query planner no bound on how far it must scan, so it sorts a filtered
   slice of the whole catalogue. Measured against the mirror: >300 s for a
   single page (never returned), against 16.8 s for an equivalent
   "source_id BETWEEN lo AND hi" range scan, which is indexed because
   source_id is the primary key. One run of this script burned 2 h 47 m
   without completing a single logged page; at that rate it needed ~47 h.

2. THE 2MASS CROSS-MATCH IGNORED PROPER MOTION. It matched Gaia epoch-2016.0
   positions against a ~1999.5 survey inside 3 arcsec. This sample is
   parallax > 2 mas, so it is nearby and fast-moving, and the fastest stars
   simply fall outside the radius. Worse, script 66 uses tmass_xm_dist <
   0.5 arcsec as its photometric-cleanliness cut: unpropagated, that
   separation measures proper-motion displacement rather than blending, so
   the "pristine" subsample would have been selected on kinematics and the
   comparison against the low-RUWE reference would have been confounded at
   the heart of the discriminant.

Script 80 carries indexed range scans, epoch propagation to 1999.5, and an
interleaved chunk sample so the sky coverage is uniform rather than one patch.

    run.sh scripts/80_pull_high_ruwe_fast.py

WHY THIS DATA DOES NOT EXIST YET
--------------------------------
The ADQL that built this project's sample carried a server-side RUWE cut. The
local partitions top out at ruwe = 1.387, contain zero duplicated sources and
zero multi-peak sources. So every result in this project -- all nineteen
channels -- was computed on a sample from which the astrometrically-disturbed
stars had already been removed, and that population was never on disk to audit.

That is the wrong population to have thrown away for this particular search.
RUWE measures how badly a single-star astrometric model fits: it is elevated
precisely when an unseen mass is pulling the star around. A completely
enshrouded companion is an unseen mass. The f -> 1 case that every photometric
channel is structurally blind to would announce itself in RUWE and nowhere
else -- and RUWE is what we filtered on.

This pulls it, with the same 2MASS join and the same volume, so the optical
deficit can be computed on it with the identical estimator.

WHAT IT IS NOT
--------------
High RUWE overwhelmingly means an ordinary unresolved binary; roughly a third
of nearby stars show it. Finding the population dim is expected and means
little on its own, because an unresolved companion inflates 2MASS Ks in a 4
arcsec beam far more than it inflates Gaia G, which manufactures a deficit.
The point of pulling it is that the question has never been asked, not that
the answer is likely to be interesting.
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
log = logging.getLogger("pull_ruwe")

MIRRORS = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
    ("ESA", "https://gea.esac.esa.int/tap-server/tap"),
]

COLUMNS = """
  g.source_id AS source_id, g.ra, g.dec, g.l, g.b,
  g.parallax, g.parallax_error, g.parallax_over_error,
  g.pmra, g.pmdec, g.ruwe, g.astrometric_excess_noise,
  g.astrometric_excess_noise_sig, g.astrometric_params_solved,
  g.ipd_frac_multi_peak, g.duplicated_source,
  g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
  g.bp_rp, g.phot_bp_rp_excess_factor, g.phot_variable_flag,
  g.non_single_star, g.classprob_dsc_combmod_star,
  g.nu_eff_used_in_astrometry, g.pseudocolour, g.ecl_lat,
  g.teff_gspphot, g.mh_gspphot
"""

PAGE = 2000
MAX_PAGES = 600
RUWE_MIN = 1.4


def pick_mirror():
    from astroquery.utils.tap.core import TapPlus
    for name, url in MIRRORS:
        t0 = time.time()
        try:
            tap = TapPlus(url=url)
            tap.launch_job("SELECT TOP 1 source_id "
                           "FROM gaiadr3.gaia_source").get_results()
            log.info("using %s -- probe %.1fs", name, time.time() - t0)
            return tap
        except Exception as exc:
            log.warning("%s unusable (%.1fs): %s", name,
                        time.time() - t0, str(exc)[:110])
    raise RuntimeError("no Gaia TAP mirror responded")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--ruwe-min", type=float, default=RUWE_MIN)
    args = ap.parse_args()

    log.error("scripts/65_pull_high_ruwe.py is SUPERSEDED and will not run.")
    log.error("  Its keyset pagination needs ~47 h against gaia_source, and its")
    log.error("  2MASS cross-match ignores proper motion, which corrupts the")
    log.error("  tmass_xm_dist cleanliness cut that Search N depends on.")
    log.error("  Use: run.sh scripts/80_pull_high_ruwe_fast.py")
    return 2

    out = cfg.RAW_DIR / "high_ruwe_500pc.parquet"
    if out.exists() and not args.force:
        log.info("cached: %s", out)
        return 0

    tap = pick_mirror()
    pages = []
    last_id = -1
    t_start = time.time()

    for page in range(MAX_PAGES):
        # The Heidelberg mirror carries DR3 but not the gaiadr1 2MASS
        # auxiliary tables, so the photometry is joined afterwards through
        # CDS XMatch rather than server-side.
        q = (f"SELECT TOP {PAGE} {COLUMNS} "
             f"FROM gaiadr3.gaia_source AS g "
             f"WHERE g.ruwe >= {args.ruwe_min} "
             f"AND g.parallax > 2 "
             f"AND g.parallax_over_error > 20 "
             f"AND g.source_id > {last_id} "
             # DaCHS rejects a table-qualified column in ORDER BY
             # ("Encountered '.'"), so order by the SELECT alias instead.
             f"ORDER BY source_id ASC")
        t0 = time.time()
        try:
            df = tap.launch_job(q).get_results().to_pandas()
        except Exception as exc:
            log.error("page %d failed: %s", page, str(exc)[:160])
            break
        if len(df) == 0:
            break
        pages.append(df)
        last_id = int(df["source_id"].iloc[-1])
        total = sum(len(p) for p in pages)
        if page % 10 == 0 or len(df) < PAGE:
            log.info("page %3d: %4d rows (%7d total) in %4.1fs",
                     page, len(df), total, time.time() - t0)
        if len(df) < PAGE:
            break
    else:
        log.warning("hit MAX_PAGES; result may be incomplete")

    if not pages:
        log.error("nothing retrieved")
        return 1

    d = pd.concat(pages, ignore_index=True).drop_duplicates("source_id")
    log.info("high-RUWE sources within 500 pc: %d in %.0fs",
             len(d), time.time() - t_start)
    log.info("  ruwe: median %.2f, 90th pct %.2f, max %.1f",
             float(d["ruwe"].median()), float(d["ruwe"].quantile(0.9)),
             float(d["ruwe"].max()))

    # ---- 2MASS photometry via CDS XMatch --------------------------------
    from astroquery.xmatch import XMatch
    from astropy.table import Table
    import astropy.units as u

    log.info("cross-matching %d positions against 2MASS ...", len(d))
    chunks = []
    step = 25000
    src = d[["source_id", "ra", "dec"]]
    for i0 in range(0, len(src), step):
        part = src.iloc[i0:i0 + step]
        t0 = time.time()
        r = XMatch.query(cat1=Table.from_pandas(part),
                         cat2="vizier:II/246/out",
                         max_distance=3.0 * u.arcsec,
                         colRA1="ra", colDec1="dec")
        log.info("  rows %d-%d -> %d matches (%.1fs)",
                 i0, min(i0 + step, len(src)), len(r), time.time() - t0)
        chunks.append(r.to_pandas())

    xm = pd.concat(chunks, ignore_index=True)
    keep = {"source_id": "source_id", "angDist": "tmass_xm_dist",
            "Kmag": "tmass_ks_m", "e_Kmag": "tmass_ks_msigcom",
            "Jmag": "tmass_j_m", "Qflg": "tmass_ph_qual"}
    keep = {k: v for k, v in keep.items() if k in xm.columns}
    xm = xm[list(keep)].rename(columns=keep)
    # nearest match per source, plus a crowding indicator
    nnb = xm.groupby("source_id").size().rename("tmass_xm_nnb")
    xm = xm.sort_values("tmass_xm_dist").drop_duplicates("source_id")
    xm = xm.merge(nnb, on="source_id", how="left")

    d = d.merge(xm, on="source_id", how="left")
    n_ks = int(d["tmass_ks_m"].notna().sum()) if "tmass_ks_m" in d else 0
    log.info("with 2MASS Ks: %d of %d (%.1f%%)",
             n_ks, len(d), 100 * n_ks / max(len(d), 1))

    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
    d.to_parquet(out, index=False)
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
