#!/usr/bin/env python
"""Figures for the multi-channel grabby-alien search.

    run.sh scripts/45_grabby_figures.py --tag primary
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

from pipeline import config as cfg       # noqa: E402

plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 140, "font.size": 9,
                     "axes.grid": True, "grid.alpha": 0.25,
                     "axes.axisbelow": True, "figure.facecolor": "white",
                     "savefig.bbox": "tight"})
POS, NEG, NEU = "#c0392b", "#2471a3", "#555555"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()
    R, F = cfg.RESULT_DIR, cfg.FIG_DIR

    sl = pd.read_parquet(cfg.DERIVED_DIR / f"spectral_slope_{args.tag}.parquet")
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["source_id", "l", "b", "A_0"])
    m = sl.merge(d, on="source_id", how="left", suffixes=("", "_d"))
    strong = m[(m["deficit_G_mag"] > 0.10) & (m["deficit_G_mag"] < 3.0)
               & (m["dchi2"] > 25)]

    # ---------- F11: spectral slope, the specificity plot -----------------
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))

    ax[0].hist(strong["alpha"], bins=np.linspace(-1, 6, 71), color=NEU,
               histtype="stepfilled", alpha=0.85)
    ax[0].axvline(2.0, color=POS, lw=2, label="interstellar dust, $\\alpha=2$")
    ax[0].axvline(0.19, color=NEG, lw=2, ls="--",
                  label="our blind spot, $\\alpha=0.19$")
    ax[0].axvline(float(strong["alpha"].median()), color="black", lw=1.2, ls=":",
                  label=f"median {strong['alpha'].median():.2f}")
    ax[0].set_xlabel(r"absorber spectral slope $\alpha$")
    ax[0].set_ylabel("stars")
    ax[0].set_title(f"Slope of {len(strong):,} significant absorbers")
    ax[0].legend(fontsize=7)

    lo = strong["alpha_upper_2sig"] < 1.5
    du = strong["alpha"].between(1.8, 2.2)
    bins = np.linspace(0, 0.6, 40)
    ax[1].hist(strong.loc[du, "A_0"], bins=bins, density=True, color=POS,
               alpha=0.55, label=r"$\alpha\approx2$ (dust-like)")
    ax[1].hist(strong.loc[lo, "A_0"], bins=bins, density=True, color=NEG,
               alpha=0.55, label=r"$\alpha<1.5$ (flat)")
    ax[1].set_xlabel("$A_0$ (mag)"); ax[1].set_ylabel("density")
    ax[1].set_title("Flat slopes track extinction\n(grain growth, not engineering)")
    ax[1].legend(fontsize=7)

    bb = np.linspace(0, 90, 46)
    ax[2].hist(np.abs(strong.loc[du, "b"]), bins=bb, density=True, color=POS,
               alpha=0.55, label=r"$\alpha\approx2$")
    ax[2].hist(np.abs(strong.loc[lo, "b"]), bins=bb, density=True, color=NEG,
               alpha=0.55, label=r"$\alpha<1.5$")
    ax[2].set_xlabel("$|b|$ (deg)"); ax[2].set_ylabel("density")
    ax[2].set_title("...and hug the Galactic plane")
    ax[2].legend(fontsize=7)
    fig.savefig(F / "F11_spectral_slope.png"); plt.close(fig)
    print("wrote F11_spectral_slope.png")

    # ---------- F12: mirror control + survivor sky map --------------------
    nul = json.loads((R / f"searchA_null_{args.tag}.json").read_text())
    surv = pd.read_csv(R / f"searchA_candidates_{args.tag}.csv")

    fig = plt.figure(figsize=(13, 4.2))
    a0 = fig.add_subplot(1, 3, 1)
    a0.bar(["dimmed\n(physical)", "brightened\n(mirror)"],
           [nul["n_dimmed"], nul["n_brightened"]],
           color=[POS, NEG], width=0.6)
    a0.set_yscale("log"); a0.set_ylabel("stars surviving identical cuts")
    a0.set_title(f"Mirror control: {nul['ratio']:.1f}:1 asymmetry\n"
                 f"(an absorber can only dim)")
    for i, v in enumerate([nul["n_dimmed"], nul["n_brightened"]]):
        a0.text(i, v * 1.15, str(v), ha="center", fontsize=9)

    a1 = fig.add_subplot(1, 3, (2, 3), projection="mollweide")
    lam = np.radians(((surv["l"] + 180) % 360) - 180)
    bet = np.radians(surv["b"])
    sc = a1.scatter(lam, bet, c=surv["deficit_G_mag"], s=9, cmap="magma_r",
                    vmin=0.1, vmax=1.0, alpha=0.85, edgecolors="none")
    a1.grid(alpha=0.3)
    a1.set_xticklabels([])
    a1.set_title(f"{len(surv):,} Search-A survivors in Galactic coordinates\n"
                 "(low extinction, $|b|>20^\\circ$, flat slope, clean)",
                 fontsize=9)
    fig.colorbar(sc, ax=a1, shrink=0.7, label="$G$-band deficit (mag)")
    fig.savefig(F / "F12_mirror_and_sky.png"); plt.close(fig)
    print("wrote F12_mirror_and_sky.png")

    # ---------- F13: channel coverage map --------------------------------
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    channels = [
        ("optical deficit vs $M_{Ks}$\n(this work, pairs)", 0.19, 6.0, 0.51, POS),
        ("IR excess\n(Suazo+22, WISE)", 0.5, 6.0, 0.50, "#16a085"),
        ("dynamical mass\n(this work)", -1.0, 6.0, 0.73, "#8e44ad"),
        ("spectral slope fit\n(this work, Search A)", 0.3, 6.0, 0.30, "#e67e22"),
    ]
    for i, (name, a_lo, a_hi, f_min, col) in enumerate(channels):
        ax.barh(i, a_hi - a_lo, left=a_lo, height=0.55, color=col, alpha=0.75)
        ax.text(a_hi + 0.1, i, f"$f\\gtrsim{f_min:.2f}$", va="center", fontsize=8)
    ax.axvline(0.0, color="black", lw=1.2, ls="-")
    ax.text(0.02, 3.6, "grey", fontsize=8, rotation=90, va="top")
    ax.axvline(2.0, color=NEU, lw=1.2, ls="--")
    ax.text(2.05, 3.6, "dust", fontsize=8, rotation=90, va="top")
    ax.axvspan(-1.0, 0.3, color="grey", alpha=0.15)
    ax.set_yticks(range(len(channels)))
    ax.set_yticklabels([c[0] for c in channels], fontsize=8)
    ax.set_xlabel(r"absorber spectral slope $\alpha$   ($\tau\propto\lambda^{-\alpha}$)")
    ax.set_xlim(-1.2, 6.6)
    ax.set_title("Which channel sees which absorber\n"
                 "grey shading: only the dynamical-mass channel reaches it",
                 fontsize=10)
    fig.savefig(F / "F13_channel_coverage.png"); plt.close(fig)
    print("wrote F13_channel_coverage.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
