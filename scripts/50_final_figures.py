#!/usr/bin/env python
"""Final figures for the multi-channel grabby-alien search report.

    run.sh scripts/50_final_figures.py --tag primary

Generates F14-F18 covering channels 9-12 and the grand summary.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

from pipeline import config as cfg

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 150, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25,
    "axes.axisbelow": True, "figure.facecolor": "white",
    "savefig.bbox": "tight",
})
POS, NEG, NEU, HI = "#c0392b", "#2471a3", "#555555", "#27ae60"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()
    tag = args.tag
    R, F = cfg.RESULT_DIR, cfg.FIG_DIR

    # ── F14: Search B — intercept-and-re-emit (W3/W4 excess) ────────────
    sb = json.loads((R / f"searchB_cold_{tag}.json").read_text())
    sb_cand = pd.read_csv(R / f"searchB_cold_candidates_{tag}.csv")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    ax = axes[0]
    labels = ["optical deficit\n(>3σ)", "W3 excess\n(>3σ)",
              "W4 excess\n(>3σ)", "BOTH\n(outside SFRs)", "mirror\ncontrol"]
    vals = [sb["n_optical_deficit"], sb["n_ir_excess_w3"],
            sb["n_ir_excess_w4"], sb["n_intercept_reemit_any"],
            sb["n_mirror_control"]]
    colors = [NEU, NEU, NEU, POS, NEG]
    bars = ax.bar(range(len(labels)), vals, color=colors, width=0.6)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_yscale("log")
    ax.set_ylabel("count")
    for i, v in enumerate(vals):
        ax.text(i, v * 1.3, f"{v:,}", ha="center", fontsize=7)
    ax.set_title(f"Search B filter chain\n({sb['signal_mirror_ratio']:.0f}:1 asymmetry)")

    ax = axes[1]
    w3 = sb_cand["resid_w3"].dropna()
    w4 = sb_cand["resid_w4"].dropna()
    ax.scatter(w3, w4, s=4, alpha=0.4, c=POS, edgecolors="none")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("W3 residual (mag, negative=IR-bright)")
    ax.set_ylabel("W4 residual (mag, negative=IR-bright)")
    ax.set_title(f"IR excess of {len(sb_cand):,} candidates\n(debris disc population)")

    ax = axes[2]
    w34 = sb_cand["w3w4_colour"].dropna()
    ax.hist(w34, bins=50, color=POS, alpha=0.8, histtype="stepfilled")
    ax.axvline(w34.median(), color="k", ls=":", lw=1.5,
               label=f"median={w34.median():.2f}")
    ax.set_xlabel("W3 − W4 (mag)")
    ax.set_ylabel("count")
    ax.set_title("W3−W4 colour → warm dust")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(F / "F14_searchB_ir_excess.png")
    plt.close(fig)
    print("wrote F14_searchB_ir_excess.png")

    # ── F15: Search A diagnosis — blend indicators ──────────────────────
    diag = json.loads((R / f"searchA_diagnosis_{tag}.json").read_text())

    fig, axes = plt.subplots(1, 4, figsize=(14, 3.5))

    indicators = [
        ("WISE xm_dist", "wise_xm_dist_surv", "wise_xm_dist_ctrl", "″"),
        ("W1−W2 colour", "w12_surv", "w12_ctrl", "mag"),
        ("M-dwarf frac", "frac_mdwarf_surv", "frac_mdwarf_ctrl", "%"),
        ("BP-RP excess", "bprp_excess_surv", "bprp_excess_ctrl", ""),
    ]

    for i, (label, sk, ck, ulabel) in enumerate(indicators):
        ax = axes[i]
        sv = diag.get(sk, 0)
        cv = diag.get(ck, 0)
        if "frac" in sk:
            sv *= 100; cv *= 100
        bars = ax.bar(["survivors", "control"], [sv, cv],
                      color=[POS, NEG], width=0.5, alpha=0.85)
        ax.set_title(label, fontsize=9)
        if ulabel:
            ax.set_ylabel(ulabel)
        for j, v in enumerate([sv, cv]):
            ax.text(j, v * 1.05, f"{v:.3f}" if v < 10 else f"{v:.1f}",
                    ha="center", fontsize=7)

    fig.suptitle("Search A diagnosis: all indicators point to unresolved NIR companions",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(F / "F15_searchA_diagnosis.png")
    plt.close(fig)
    print("wrote F15_searchA_diagnosis.png")

    # ── F16: Search D — domain edge diagnosis ───────────────────────────
    sd = json.loads((R / f"searchD_domain_edge_{tag}.json").read_text())

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    ax = axes[0]
    methods = ["planar\nscan", "radial\nscan", "local\ngradient"]
    sigmas = [sd["planar_scan"]["excess_over_null_sigma"],
              sd["radial_scan"]["excess_over_null_sigma"],
              sd["local_gradient"]["excess_over_null_sigma"]]
    colors_bar = [POS if s > 5 else NEG if s < -2 else NEU for s in sigmas]
    ax.bar(methods, sigmas, color=colors_bar, width=0.5, alpha=0.85)
    ax.axhline(5, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel("excess over null (σ)")
    ax.set_title("Three independent edge detectors")
    for i, v in enumerate(sigmas):
        ax.text(i, max(v, 0) + 3, f"{v:.1f}σ", ha="center", fontsize=8)

    ax = axes[1]
    ps = sd["planar_scan"]
    ax.text(0.5, 0.7, f"Planar: {ps['excess_over_null_sigma']:.1f}σ",
            transform=ax.transAxes, ha="center", fontsize=14, color=POS)
    ax.text(0.5, 0.5, f"Δ = {abs(ps['best_delta_mag']):.4f} mag",
            transform=ax.transAxes, ha="center", fontsize=11)
    lb = ps["best_normal_lbdeg"]
    ax.text(0.5, 0.3, f"normal → l={lb[0]:.0f}°, b={lb[1]:.0f}°",
            transform=ax.transAxes, ha="center", fontsize=10)
    ax.text(0.5, 0.12, "(North Galactic Cap = lowest extinction)",
            transform=ax.transAxes, ha="center", fontsize=8, color=NEU)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Strongest planar signal")

    ax = axes[2]
    verdicts = [
        ("A₀<0.05 reduces by 40%", True),
        ("local gradient null (−0.4σ)", True),
        ("points at low-extinction cap", True),
        ("Gaia scanning-law ~0.01 mag", True),
    ]
    for i, (txt, kills) in enumerate(verdicts):
        ax.text(0.05, 0.85 - i * 0.2, "✗ " + txt if kills else "? " + txt,
                transform=ax.transAxes, fontsize=9,
                color=NEG if kills else NEU)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Diagnosis: extinction + scanning law")

    fig.tight_layout()
    fig.savefig(F / "F16_searchD_domain_edge.png")
    plt.close(fig)
    print("wrote F16_searchD_domain_edge.png")

    # ── F17: Search C — radio null ──────────────────────────────────────
    sc = json.loads((R / f"searchC_radio_{tag}.json").read_text())

    fig, ax = plt.subplots(figsize=(7, 4))
    cats = []
    cand_rates = []
    ctrl_rates = []
    for cat_name, info in sc["catalogs"].items():
        if "error" in info:
            continue
        cats.append(cat_name)
        cand_rates.append(info["n_candidate_matches"])
        ctrl_rates.append(info["n_control_matches"])

    x = np.arange(len(cats))
    w = 0.35
    ax.bar(x - w/2, cand_rates, w, color=POS, alpha=0.85, label="dimmed candidates")
    ax.bar(x + w/2, ctrl_rates, w, color=NEG, alpha=0.85, label="control")
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("radio matches (out of 5,000)")
    ax.set_title("Search C: radio counterparts\n"
                 "Main-sequence stars are radio-quiet (as expected)")
    ax.legend(fontsize=8)
    for i in range(len(cats)):
        ax.text(x[i] - w/2, cand_rates[i] + 0.15, str(cand_rates[i]),
                ha="center", fontsize=9)
        ax.text(x[i] + w/2, ctrl_rates[i] + 0.15, str(ctrl_rates[i]),
                ha="center", fontsize=9)
    ax.set_ylim(0, max(max(cand_rates + ctrl_rates) + 2, 6))

    fig.tight_layout()
    fig.savefig(F / "F17_searchC_radio.png")
    plt.close(fig)
    print("wrote F17_searchC_radio.png")

    # ── F18: Grand summary — all 12 channels ────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 7))

    channels = [
        (1, "optical deficit (single)", "p < 6.1e−3", "null", NEU),
        (2, "optical deficit (wide pair)", "p < 5.3e−4", "null", NEU),
        (3, "mid-IR veto (beamed)", "p < 4.3e−4", "null", NEU),
        (4, "dynamical mass (grey-safe)", "p < 2.2e−4, 0 positives", "null", HI),
        (5, "spectral slope α", "372 asymmetric survivors", "interesting", POS),
        (6, "3D spatial front", "19.5σ → dust", "explained", NEG),
        (7, "kinematic clustering", "4301σ → young stars", "explained", NEG),
        (8, "deficit colour type", "all reddening-like", "null", NEU),
        (9, "W3/W4 intercept+re-emit", "2,632 → debris discs", "explained", NEG),
        (10, "radio (NVSS/FIRST)", "1/5000 (0.25× control)", "null", NEU),
        (11, "3D domain edge", "148.5σ → extinction+scanning", "explained", NEG),
        (12, "Search A diagnosis", "NIR companions", "explained", NEG),
    ]

    y = np.arange(len(channels))[::-1]
    cat_colors = {"null": NEU, "explained": NEG, "interesting": POS}

    for i, (ch, name, result, cat, col) in enumerate(channels):
        ax.barh(y[i], 1.0, height=0.6, color=col, alpha=0.7)
        ax.text(0.02, y[i], f"Ch {ch}: {name}", va="center", fontsize=8.5,
                fontweight="bold" if cat == "interesting" else "normal")
        ax.text(1.02, y[i], result, va="center", fontsize=8, ha="left")

    ax.set_xlim(0, 2.4)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=NEU, alpha=0.7, label="null (no signal)"),
        Patch(facecolor=NEG, alpha=0.7, label="explained (mundane population)"),
        Patch(facecolor=POS, alpha=0.7, label="interesting (real signal, likely blends)"),
        Patch(facecolor=HI, alpha=0.7, label="strongest constraint"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    ax.set_title("All 12 search channels: multi-band, multi-method sweep\n"
                 "3.3M Gaia DR3 main-sequence stars within 500 pc\n"
                 "Joint constraint: p < 6.2×10⁻⁴ (fewer than 1 in 1,614 "
                 "stars intercepts ≥51% of optical output)",
                 fontsize=10, pad=10)

    fig.tight_layout()
    fig.savefig(F / "F18_grand_summary.png")
    plt.close(fig)
    print("wrote F18_grand_summary.png")

    print("\nAll final figures written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
