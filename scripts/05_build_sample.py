#!/usr/bin/env python
"""Step 1 (client side): raw TAP chunks -> analysis sample.

    run.sh scripts/05_build_sample.py --pattern 'sample_d500_p*' --tag primary

Writes:
    data/cache/derived/<tag>.parquet
    results/cutflow_<tag>.csv
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import sample as smp

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("build")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="sample_d500_p*")
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--dust-map", default="edenhofer23")
    ap.add_argument("--band-law", default="fitz19")
    ap.add_argument("--all-maps", default="edenhofer23",
                    help="comma list of maps to attach A_0 for")
    ap.add_argument("--distance-max-pc", type=float, default=None)
    ap.add_argument("--max-chunks", type=int, default=None)
    args = ap.parse_args()

    files = sorted(cfg.RAW_DIR.glob(args.pattern + ".parquet"))
    if args.max_chunks:
        files = files[: args.max_chunks]
    if not files:
        log.error("no chunks match %s in %s", args.pattern, cfg.RAW_DIR)
        return 1
    log.info("loading %d chunks", len(files))
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    log.info("loaded %d rows", len(df))

    # A source can appear once per chunk only, but guard anyway.
    n0 = len(df)
    df = df.drop_duplicates(subset="source_id").reset_index(drop=True)
    if len(df) != n0:
        log.warning("dropped %d duplicate source_id rows", n0 - len(df))

    flow = smp.CutFlow()

    log.info("applying parallax zero point ...")
    df = smp.add_astrometry(df)
    frac_faint_clip = float((df["phot_g_mean_mag"] < 6.0).mean())
    log.info("  G < 6 (zero point held at G=6): %.4f%% of sample",
             100 * frac_faint_clip)
    log.info("  parallax zero point: median %.1f uas, 16-84%% %.1f..%.1f uas",
             1000 * df["parallax_zp"].median(),
             1000 * df["parallax_zp"].quantile(0.16),
             1000 * df["parallax_zp"].quantile(0.84))

    maps = [m.strip() for m in args.all_maps.split(",") if m.strip()]
    df = smp.add_extinction(df, maps=tuple(maps))
    df = smp.add_absolute_magnitudes(df, args.dust_map, args.band_law)
    df = smp.add_sky_density(df)

    df = smp.apply_quality_cuts(df, flow)
    df = smp.apply_extinction_and_ms_cuts(df, flow,
                                          distance_max_pc=args.distance_max_pc)

    out = cfg.DERIVED_DIR / f"{args.tag}.parquet"
    df.to_parquet(out, index=False, compression="zstd")
    log.info("wrote %s  (%d rows, %d cols)", out, len(df), df.shape[1])

    cf = flow.to_frame()
    cf_path = cfg.RESULT_DIR / f"cutflow_{args.tag}.csv"
    cf.to_csv(cf_path, index=False)
    print("\n=== cut flow ===")
    print(cf.to_string(index=False))

    print("\n=== sample summary ===")
    for col in ["dist_pc", "M_G", "M_Ks", "bp_rp0", "A_0", "A_G", "A_Ks",
                "mh_gspphot", "cstar_nsigma"]:
        if col in df:
            s = df[col]
            print(f"  {col:14s} n={s.notna().sum():8d}  "
                  f"median={s.median():8.3f}  "
                  f"16-84%={s.quantile(0.16):8.3f}..{s.quantile(0.84):8.3f}")
    print(f"\n  |b| < 10 deg fraction : {(df['b'].abs() < 10).mean():.3f}")
    print(f"  north (b>0) fraction  : {(df['b'] > 0).mean():.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
