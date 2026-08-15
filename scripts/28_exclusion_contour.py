#!/usr/bin/env python
"""The paper's headline figure: exclusion in the (f, p) plane.

    run.sh scripts/28_exclusion_contour.py

f = optical fraction intercepted per harvested star
p = fraction of stars that are harvested

Every constraint in this work is a curve in this plane, and so is every prior
limit worth comparing against.  Plotting them together is the only honest way
to show that the channels are complementary rather than competing: the
infrared-excess searches own the high-f, low-p corner for objects with warm
dust, and this work owns a different region and is blind in others.

The shaded blind regions are as important as the excluded ones.
"""
from __future__ import annotations

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


def main() -> int:
    R = cfg.RESULT_DIR
    excl_a = pd.read_csv(R / "exclusion_primary.csv")
    excl_b = pd.read_csv(R / "exclusion_primary_optcol.csv")
    pair = json.loads((R / "pair_limit_primary.json").read_text())
    pair_tab = pd.DataFrame(pair["table"])

    fig, ax = plt.subplots(figsize=(7.6, 5.6))

    def curve(df, fcol, pcol, **kw):
        d = df[(df[pcol] < 0.999) & (df[pcol] > 0)].sort_values(fcol)
        if len(d):
            ax.plot(d[fcol], d[pcol], **kw)
        return d

    a = curve(excl_a, "f", "p_upper_limit", color="#c0392b", lw=2.0,
              marker="o", ms=4, label="This work: single star, any selective absorber")
    b = curve(excl_b, "f", "p_upper_limit", color="#e67e22", lw=2.0,
              marker="s", ms=4, ls="--",
              label="This work: single star, grey-across-optical only")
    pt = pair_tab[pair_tab["p_UL"] < 1].sort_values("f_detectable")
    ax.plot(pt["f_detectable"], pt["p_UL"], color="#2471a3", lw=2.4,
            marker="D", ms=5, label="This work: clean wide pairs (best)")

    # Everything above a curve is excluded; shade the best one.
    if len(pt):
        ax.fill_between(pt["f_detectable"], pt["p_UL"], 1.0,
                        color="#2471a3", alpha=0.10)

    # --- blind regions ----------------------------------------------------
    ax.axvspan(1e-3, 0.10, color="grey", alpha=0.13)
    ax.text(0.0095, 1.2e-5, "BLIND\nno constraint\nat any $p$",
            fontsize=9, ha="center", va="center", color="#333333",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#888888",
                      alpha=0.85))

    # --- prior work -------------------------------------------------------
    # Zackrisson et al. 2015: <3% of local disc galaxies host galaxy-scale
    # astroengineering dimming the Tully-Fisher relation. Not the same axis
    # (galaxies, not stars) -- drawn as an annotation, not a curve.
    ax.axhline(0.03, color="#7f8c8d", ls=":", lw=1.4)
    ax.text(0.62, 0.036, "Zackrisson+15: <3% of disc galaxies\n"
                          "(galaxy-scale, different population)",
            fontsize=7.5, color="#555555")

    # Suazo et al. 2024 (Hephaistos II): 7 candidates in ~5e6 sources, i.e.
    # a rate <~1.4e-6, but ONLY for objects with detectable warm dust
    # re-emission. Drawn as a point with an arrow to mark the different regime.
    ax.plot([0.5], [1.4e-6], marker="*", ms=15, color="#16a085", ls="none",
            label="Hephaistos II (IR excess; requires warm re-emission)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-3, 1.0)
    ax.set_ylim(1e-6, 1.0)
    ax.set_xlabel("$f$  —  optical fraction intercepted per harvested star")
    ax.set_ylabel("$p$  —  fraction of stars harvested")
    ax.set_title("95% CL exclusion in the $(f, p)$ plane\n"
                 "shaded above the blue curve is excluded by this work",
                 fontsize=11)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.5, loc="upper right", bbox_to_anchor=(1.0, 0.62),
              framealpha=0.95)

    # constant mean-harvested-fraction diagonals
    for fbar, lab in [(1e-2, r"$\bar f=10^{-2}$"), (1e-3, r"$\bar f=10^{-3}$"),
                      (1e-4, r"$\bar f=10^{-4}$")]:
        ff = np.logspace(-3, 0, 50)
        ax.plot(ff, fbar / ff, color="black", lw=0.6, ls="-.", alpha=0.45)
        ax.text(0.85, fbar / 0.85 * 1.15, lab, fontsize=6.5, alpha=0.7)

    path = cfg.FIG_DIR / "F10_exclusion_contour.png"
    fig.savefig(path, dpi=140, bbox_inches="tight")
    print(f"wrote {path}")

    # machine-readable summary for the paper
    summary = {
        "single_star_variantA_best_mean_f": float(
            excl_a[excl_a["p_upper_limit"] < 0.999]["mean_f_upper_limit"].min()),
        "single_star_variantB_best_mean_f": float(
            excl_b[excl_b["p_upper_limit"] < 0.999]["mean_f_upper_limit"].min()),
        "pair_best_mean_f": float(pair["best_mean_f_UL"]),
        "pair_best_p": float(pair["best_p_UL"]),
        "pair_f_at_best": float(pair["f_detectable_at_best"]),
        "blind_below_f": 0.1,
    }
    (R / "exclusion_contour.json").write_text(json.dumps(summary, indent=2))
    for k, v in summary.items():
        print(f"  {k:36s} {v:.4g}" if isinstance(v, float) else f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
