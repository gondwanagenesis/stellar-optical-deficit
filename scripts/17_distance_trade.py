#!/usr/bin/env python
"""Sensitivity-versus-systematics trade study across the distance limit.

    run.sh scripts/17_distance_trade.py --knots 6 --mh-degree 1

Pushing the distance limit out buys stars -- volume goes as d^3 until the
2MASS Ks depth bites -- which improves sigma/sqrt(N). It also buys extinction,
which is the dominant systematic. This script measures both sides on the SAME
sky partitions so the comparison is like for like, and answers whether the
extra stars are worth having.

All distance cuts use the same spline hyperparameters so that differences are
attributable to the sample and not to model complexity.
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
from pipeline import extinction as ext
from pipeline import fiducial as fid
from pipeline import nulls
from pipeline import sample as smp

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("trade")

DISTANCES = [200.0, 300.0, 500.0, 750.0, 1000.0, 1250.0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="sample_d1250_p*")
    ap.add_argument("--knots", type=int, default=6)
    ap.add_argument("--mh-degree", type=int, default=1)
    ap.add_argument("--nir-control", action="store_true", default=True)
    ap.add_argument("--a-g-max", type=float, default=None,
                    help="override the A_G cut; default uses the config value "
                         "at every distance so the comparison isolates distance")
    args = ap.parse_args()

    files = sorted(cfg.RAW_DIR.glob(args.pattern + ".parquet"))
    if not files:
        log.error("no chunks match %s", args.pattern)
        return 1
    log.info("loading %d wide chunks", len(files))
    raw = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    raw = raw.drop_duplicates(subset="source_id").reset_index(drop=True)
    log.info("%d rows", len(raw))

    raw = smp.add_astrometry(raw)
    raw = smp.add_extinction(raw, maps=("edenhofer23",))
    raw = smp.add_absolute_magnitudes(raw, "edenhofer23", "fitz19")
    raw = smp.add_sky_density(raw)
    a_j = ext.deredden("J", np.nan_to_num(raw["A_0"].to_numpy(float)),
                       raw["bp_rp"].to_numpy(float), law="fitz19")
    raw["j_ks0"] = ((raw["tmass_j_m"].to_numpy(float) - a_j)
                    - (raw["tmass_ks_m"].to_numpy(float)
                       - raw["A_Ks"].to_numpy(float)))

    flow0 = smp.CutFlow()
    base = smp.apply_quality_cuts(raw, flow0)

    rows = []
    for dmax in DISTANCES:
        flow = smp.CutFlow()
        d = smp.apply_extinction_and_ms_cuts(base, flow, distance_max_pc=dmax,
                                             a_g_max=args.a_g_max)
        d = d[d["mh_gspphot"].notna() & np.isfinite(d["j_ks0"])].reset_index(drop=True)
        if len(d) < 5000:
            log.warning("d<%.0f pc: only %d stars, skipping", dmax, len(d))
            continue

        covs = [(d["mh_gspphot"].to_numpy(float), args.mh_degree)]
        if args.nir_control:
            covs.append((d["j_ks0"].to_numpy(float), 1))
        fit = fid.fit_fiducial(d["M_Ks"].to_numpy(float), covs,
                               d["M_G"].to_numpy(float), args.knots)
        resid = fit.residuals(d["M_Ks"].to_numpy(float), covs,
                              d["M_G"].to_numpy(float))
        sigma = fit.sigma_robust

        splits = nulls.standard_splits(d, resid)
        worst = splits.loc[splits["n_sigma"].abs().idxmax()]
        spatial = nulls.spatial_scaling(d, resid)
        big = spatial[spatial["mean_group_n"] > 0.01 * len(d)]
        plateau = float(np.nanmedian(big["excess"])) if len(big) else np.nan

        floor = float(np.nanmax([abs(float(worst["difference"])),
                                 plateau if np.isfinite(plateau) else 0.0]))
        rows.append({
            "distance_max_pc": dmax,
            "n_stars": len(d),
            "median_A_0": float(np.nanmedian(d["A_0"])),
            "p90_A_0": float(np.nanpercentile(d["A_0"], 90)),
            "sigma_mag": float(sigma),
            "naive_sigma_over_sqrtN": float(sigma / np.sqrt(len(d))),
            "spatial_plateau_mag": plateau,
            "worst_split": str(worst["test"]),
            "worst_split_mag": float(worst["difference"]),
            "floor_mag": floor,
            "floor_over_naive": floor / (sigma / np.sqrt(len(d))),
            "floor_implied_f": float(1 - 10 ** (-floor / 2.5)),
        })
        log.info("d<%4.0f pc: N=%7d sigma=%.4f floor=%.4f (%.0fx naive)",
                 dmax, len(d), sigma, floor, rows[-1]["floor_over_naive"])

    out = pd.DataFrame(rows)
    out.to_csv(cfg.RESULT_DIR / "distance_trade.csv", index=False)
    print("\n=== distance trade study ===")
    print(out.to_string(index=False, float_format=lambda v: f"{v:12.5g}"))

    if len(out):
        best = out.loc[out["floor_mag"].idxmin()]
        print(f"\nlowest floor at d < {best['distance_max_pc']:.0f} pc: "
              f"{best['floor_mag']:.4f} mag with N = {best['n_stars']:,}")
        print("Note: N grows with distance while the floor does too, so the "
              "optimum is NOT\nthe largest sample -- which is the whole point "
              "of the trade study.")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    ax.plot(out["distance_max_pc"], out["naive_sigma_over_sqrtN"], "o--",
            color="black", label=r"naive $\sigma/\sqrt{N}$")
    ax.plot(out["distance_max_pc"], out["floor_mag"], "s-", color="#c0392b",
            label="measured floor")
    ax.plot(out["distance_max_pc"], out["spatial_plateau_mag"], "^:",
            color="#e67e22", label="spatial plateau")
    ax.set_yscale("log")
    ax.set_xlabel("distance limit (pc)")
    ax.set_ylabel("mag")
    ax.grid(alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(out["distance_max_pc"], out["n_stars"], "d-", color="#2471a3",
             alpha=0.6)
    ax2.set_ylabel("stars in sample", color="#2471a3")
    ax2.set_yscale("log"); ax2.grid(False)
    ax.legend(fontsize=8, loc="center left")
    ax.set_title("More stars, more extinction: the trade")
    fig.savefig(cfg.FIG_DIR / "F9_distance_trade.png", dpi=130, bbox_inches="tight")
    print(f"\nwrote {cfg.FIG_DIR / 'F9_distance_trade.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
