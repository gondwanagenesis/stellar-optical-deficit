#!/usr/bin/env python
"""Step 4: injection-recovery.

    run.sh scripts/13_injection.py --tag testrun_nir --max-n 200000
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import fiducial as fid
from pipeline import injection as inj
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("inject")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="testrun_nir")
    ap.add_argument("--max-n", type=int, default=None)
    ap.add_argument("--n-realisations", type=int, default=None)
    ap.add_argument("--k", type=float, default=cfg.FIT.outlier_nsigma)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet")
    meta = json.loads((cfg.RESULT_DIR / f"fiducial_{args.tag}.json").read_text())
    knots, mh_degree = meta["n_interior_knots"], meta["mh_degree"]
    nir = bool(meta.get("nir_control", False))
    cov_cols = ("mh_gspphot", "j_ks0") if nir else ("mh_gspphot",)

    resid0 = d["residual"].to_numpy(dtype=float)
    covs0 = [(d["mh_gspphot"].to_numpy(float), mh_degree)]
    if nir:
        covs0.append((d["j_ks0"].to_numpy(float), 1))
    fit0 = fid.fit_fiducial(d["M_Ks"].to_numpy(float), covs0,
                            d["M_G"].to_numpy(float), knots)

    log.info("campaign: N=%d knots=%d mh_degree=%d nir=%s",
             len(d), knots, mh_degree, nir)
    camp = inj.injection_campaign(
        d, fit0, resid0, knots, mh_degree, cov_columns=cov_cols,
        n_realisations=args.n_realisations, k=args.k, max_n=args.max_n)
    camp.to_parquet(cfg.DERIVED_DIR / f"injection_{args.tag}.parquet",
                    index=False, compression="zstd")

    summ = inj.summarise(camp)
    summ.to_csv(cfg.RESULT_DIR / f"injection_{args.tag}.csv", index=False)

    print("\n=== uniform injection (self-calibrated recovery must be ~0) ===")
    u = summ[summ["mode"] == "uniform"]
    print(u[["f_injected", "delta_injected", "mean_residual",
             "mean_residual_std", "n_pos_5sig", "anchored_f"]].to_string(
        index=False, float_format=lambda v: f"{v:12.6f}"))

    print("\n=== sparse injection (recovery is in the tail) ===")
    s = summ[summ["mode"] == "sparse"]
    print(s[["f_injected", "p_injected", "p_recovered", "p_recovered_std",
             "p_recovered_over_injected", "n_pos_5sig"]].to_string(
        index=False, float_format=lambda v: f"{v:12.6f}"))

    # Smallest recoverable p at each f: first p whose recovery is >3 sigma from 0
    print("\n=== smallest recoverable p per f (>3 sigma) ===")
    rows = []
    for f, sub in s.groupby("f_injected"):
        sub = sub.sort_values("p_injected")
        ok = sub[sub["p_recovered"] > 3 * sub["p_recovered_std"].replace(0, np.nan)]
        pmin = float(ok["p_injected"].min()) if len(ok) else np.nan
        rows.append({"f": f, "p_min_recoverable": pmin,
                     "mean_f_min": pmin * f if np.isfinite(pmin) else np.nan})
    rec = pd.DataFrame(rows)
    rec.to_csv(cfg.RESULT_DIR / f"injection_threshold_{args.tag}.csv", index=False)
    print(rec.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
