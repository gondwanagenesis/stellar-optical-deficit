#!/usr/bin/env python
"""Validate the production ADQL on a small slice before the full pull.

Pulls TOP 1000 from one partition, checks join sanity, and estimates the total
sample size by counting a few partitions and extrapolating.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline.adql import chunk_query, count_query, partition_bounds
from pipeline.tap import run_adql

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 100)

# Three partitions spread over the sky: one high-|b| north, one near the
# galactic plane, one high-|b| south.  Partition index maps to HEALPix ring
# order at nside=4 (nested), so we just take widely separated indices.
PROBE_PARTITIONS = [10, 95, 180]


def main() -> int:
    idx = PROBE_PARTITIONS[0]
    lo, hi = partition_bounds(idx)
    print(f"partition {idx}: source_id {lo} .. {hi}")

    q = chunk_query(idx, top=1000)
    print("\n--- ADQL ---")
    print(q)

    df = run_adql(q, name=f"validate_top1000_p{idx}", force=True)
    print(f"\nrows={len(df)}  cols={df.shape[1]}")
    print("\ndtypes:")
    print(df.dtypes.to_string())

    print("\nhead:")
    show = ["source_id", "l", "b", "parallax", "parallax_over_error", "ruwe",
            "phot_g_mean_mag", "bp_rp", "tmass_ks_m", "tmass_ks_msigcom",
            "tmass_ph_qual", "mh_gspphot", "ag_gspphot", "wise_w1mpro",
            "wise_w4mpro"]
    print(df[show].head(12).to_string())

    print("\nnull fraction (%) for columns with any nulls:")
    nf = (df.isna().mean() * 100).round(1)
    print(nf[nf > 0].sort_values(ascending=False).to_string())

    print("\nsanity checks")
    print(f"  duplicate source_id            : {int(df['source_id'].duplicated().sum())}")
    print(f"  parallax_over_error min        : {df['parallax_over_error'].min():.1f}")
    print(f"  ruwe max                       : {df['ruwe'].max():.3f}")
    print(f"  distance range (pc)            : "
          f"{1000/df['parallax'].max():.1f} .. {1000/df['parallax'].min():.1f}")
    m_ks = df["tmass_ks_m"] + 5 * np.log10(df["parallax"]) - 10
    print(f"  uncorrected M_Ks range         : {m_ks.min():.2f} .. {m_ks.max():.2f}")
    print(f"  2MASS ph_qual Ks == 'A' frac   : "
          f"{(df['tmass_ph_qual'].str[2] == 'A').mean():.3f}")
    print(f"  AllWISE match frac             : {df['allwise_designation'].notna().mean():.3f}")
    print(f"  GSP-Phot mh available frac     : {df['mh_gspphot'].notna().mean():.3f}")
    print(f"  GSP-Spec mh available frac     : {df['mh_gspspec'].notna().mean():.3f}")

    # --- size estimate -----------------------------------------------------
    print("\n=== partition counts (for total-size extrapolation) ===")
    counts = {}
    for p in PROBE_PARTITIONS:
        cq = count_query(p)
        n = int(run_adql(cq, name=f"count_p{p}", force=True)["n"].iloc[0])
        counts[p] = n
        print(f"  partition {p:3d}: {n:8d}")
    mean_n = float(np.mean(list(counts.values())))
    print(f"\n  mean per partition: {mean_n:.0f}")
    print(f"  extrapolated total over {cfg.N_PARTITIONS} partitions: "
          f"{mean_n * cfg.N_PARTITIONS:,.0f}")
    print("  (crude -- partitions are equal-area but stellar density is not)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
