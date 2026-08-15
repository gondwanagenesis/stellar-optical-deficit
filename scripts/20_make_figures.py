#!/usr/bin/env python
"""Regenerate every figure. One script, one command.

    run.sh scripts/20_make_figures.py --tag testrun_nir
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import figures as figs


def maybe(path: Path):
    return pd.read_csv(path) if path.exists() else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="testrun_nir")
    args = ap.parse_args()
    t = args.tag
    made = []

    meta_path = cfg.RESULT_DIR / f"fiducial_{t}.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    sigma = meta.get("sigma_observed_mag", 0.1)

    resid_path = cfg.DERIVED_DIR / f"{t}_resid.parquet"
    if resid_path.exists():
        d = pd.read_parquet(resid_path)
        made.append(figs.fig_sample(d))
        made.append(figs.fig_residual_distribution(
            d["residual"].to_numpy(float), sigma))
        n_total = len(d)
    else:
        d, n_total = None, 1

    scaling = maybe(cfg.RESULT_DIR / f"nulls_{t}_scaling.csv")
    random = maybe(cfg.RESULT_DIR / f"nulls_{t}_random.csv")
    if scaling is not None and random is not None:
        made.append(figs.fig_floor(scaling, random, sigma, n_total))

    splits = maybe(cfg.RESULT_DIR / f"nulls_{t}_splits.csv")
    floor_path = cfg.RESULT_DIR / f"floor_{t}.json"
    floor = json.loads(floor_path.read_text()) if floor_path.exists() else {}
    if splits is not None:
        made.append(figs.fig_splits(splits, floor.get("floor_all_axes_mag")))

    lev_path = cfg.RESULT_DIR / "spectral_leverage.json"
    lev_tab = maybe(cfg.RESULT_DIR / "spectral_leverage_table.csv")
    if lev_path.exists() and lev_tab is not None:
        lev = json.loads(lev_path.read_text())
        made.append(figs.fig_leverage(lev_tab, lev["slope_dMG_dMKs"],
                                      lev["alpha_blind_numeric"]))

    inj = maybe(cfg.RESULT_DIR / f"injection_{t}.csv")
    if inj is not None:
        made.append(figs.fig_injection(inj))

    excl = maybe(cfg.RESULT_DIR / f"exclusion_{t}.csv")
    if excl is not None:
        made.append(figs.fig_exclusion(excl))

    cf = maybe(cfg.RESULT_DIR / f"cutflow_{meta.get('source_tag', t)}.csv")
    if cf is None:
        cf = maybe(cfg.RESULT_DIR / f"cutflow_{t}.csv")
    if cf is not None:
        made.append(figs.fig_cutflow(cf))

    print(f"wrote {len(made)} figures:")
    for m in made:
        print("  " + m)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
