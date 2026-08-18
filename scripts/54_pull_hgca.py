#!/usr/bin/env python
"""Pull the Hipparcos-Gaia Catalog of Accelerations (Brandt 2021) for Search F.

    run.sh scripts/54_pull_hgca.py

WHY THIS CATALOGUE
------------------
Gaia DR3 gives position and velocity, never acceleration. But Hipparcos
observed the same stars ~25 years earlier, and the difference between

    mu_HG    the mean proper motion over the 1991->2016 baseline
    mu_Gaia  the near-instantaneous proper motion at the Gaia epoch

is a measured change in velocity: an acceleration proxy, for ~115,000 stars.

Brandt (2021, ApJS 254, 42) publishes both proper motions on a common,
cross-calibrated reference frame with inflated and validated uncertainties.
That calibration is the hard part and the reason we take his catalogue rather
than differencing the raw archives ourselves.

WHAT SEARCH F DOES WITH IT
--------------------------
The proper-motion anomaly is already a standard companion-detection tool, so
a large anomaly on a single star means "unseen companion", not "engine". The
discriminant is not magnitude, it is DIRECTION: an unseen companion points the
anomaly wherever the orbit happens to lie, uniformly at random across a
population. Coordinated thrust across a region would not.

So the search statistic is local directional coherence, and the two systematics
that also produce coherence must come out first:

  * the Hipparcos-to-Gaia frame rotation residual, which is global and
    removable as a rigid spin;
  * perspective acceleration from radial motion, which Brandt supplies
    per-star as (dpmRA, dpmDE).
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("pull_hgca")

VIZIER_CAT = "J/ApJS/254/42/catalog"

WANT = [
    "HIP", "Gaia", "RA_ICRS", "DE_ICRS", "RV", "e_RV",
    "plx", "e_plx",
    "pmRA", "pmDE", "e_pmRA", "e_pmDE", "pmRApmDEcor",
    "pmRAhg", "pmDEhg", "e_pmRAhg", "e_pmDEhg", "pmRApmDEhg",
    "pmRAhip", "pmDEhip", "e_pmRAhip", "e_pmDEhip", "pmRApmDEhip",
    "EpochRAgaia", "EpochDEgaia", "EpochRAhip", "EpochDEhip",
    "dpmRA", "dpmDE", "chi2",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out = cfg.RAW_DIR / "hgca_brandt2021.parquet"
    if out.exists() and not args.force:
        df = pd.read_parquet(out)
        log.info("cached: %s (%d rows)", out, len(df))
        return 0

    from astroquery.vizier import Vizier

    v = Vizier(columns=WANT, row_limit=-1)
    log.info("querying VizieR %s ...", VIZIER_CAT)
    tables = v.get_catalogs(VIZIER_CAT)
    t = tables[0]
    log.info("retrieved %d rows, %d columns", len(t), len(t.colnames))

    df = t.to_pandas()

    # VizieR hands back the Gaia identifier as a string column.
    if "Gaia" in df.columns:
        df["source_id"] = pd.to_numeric(df["Gaia"], errors="coerce").astype("Int64")

    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    log.info("wrote %s (%d rows)", out, len(df))

    n_plx = int((df["plx"] > 2).sum()) if "plx" in df else 0
    log.info("within 500 pc (plx > 2 mas): %d", n_plx)
    log.info("with a Gaia source_id: %d", int(df["source_id"].notna().sum()))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
