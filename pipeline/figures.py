"""Figure generation. One figure per claim; all regenerable by
scripts/20_make_figures.py."""

from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

from . import config as cfg              # noqa: E402

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "figure.facecolor": "white", "savefig.bbox": "tight",
})

POS = "#c0392b"      # deficit-like / fainter
NEG = "#2471a3"      # over-luminous
NEU = "#555555"


def _save(fig, name: str) -> str:
    path = cfg.FIG_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    return str(path)


# --------------------------------------------------------------------------

def fig_sample(d: pd.DataFrame, name="F1_sample_and_fiducial") -> str:
    """Claim: the sample occupies a clean, single-valued lower main sequence."""
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))

    h = ax[0].hexbin(d["M_Ks"], d["M_G"], gridsize=110, bins="log",
                     cmap="viridis", mincnt=1)
    q = d.groupby(pd.cut(d["M_Ks"], 60), observed=True).agg(
        x=("M_Ks", "median"), y=("M_G", "median"))
    ax[0].plot(q["x"], q["y"], color="white", lw=1.6)
    ax[0].plot(q["x"], q["y"], color="black", lw=0.8, ls="--",
               label="running median")
    ax[0].invert_yaxis(); ax[0].invert_xaxis()
    ax[0].set_xlabel(r"$M_{K_s}$"); ax[0].set_ylabel(r"$M_G$")
    ax[0].set_title("Main-sequence box")
    ax[0].legend(loc="lower left", fontsize=7)
    fig.colorbar(h, ax=ax[0], label="stars/bin")

    ax[1].hexbin(d["bp_rp0"], d["M_G"], gridsize=110, bins="log",
                 cmap="magma", mincnt=1)
    ax[1].invert_yaxis()
    ax[1].set_xlabel(r"$(BP-RP)_0$"); ax[1].set_ylabel(r"$M_G$")
    ax[1].set_title("Colour-magnitude")

    ax[2].hist(d["A_0"], bins=80, color=NEU)
    ax[2].set_xlabel(r"$A_0$ (mag)"); ax[2].set_ylabel("stars")
    ax[2].set_title("Extinction distribution")
    ax[2].set_yscale("log")

    fig.suptitle(f"Sample: N = {len(d):,}", y=1.02)
    return _save(fig, name)


def fig_residual_distribution(resid: np.ndarray, sigma: float,
                              name="F2_residual_distribution") -> str:
    """Claim: the residual distribution has a strong one-sided positive tail."""
    r = resid[np.isfinite(resid)]
    med = float(np.median(r))
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))

    lim = 6 * sigma
    bins = np.linspace(med - lim, med + lim, 220)
    ax[0].hist(r, bins=bins, color=NEU, histtype="stepfilled", alpha=0.85)
    g = (len(r) * (bins[1] - bins[0]) / (sigma * np.sqrt(2 * np.pi))
         * np.exp(-0.5 * ((bins - med) / sigma) ** 2))
    ax[0].plot(bins, g, color=POS, lw=1.2, label="Gaussian, same $\\sigma$")
    ax[0].set_yscale("log")
    # Counts below one star are not meaningful; without this the Gaussian
    # reference curve trails down to 1e-5 and stretches the axis pointlessly.
    ax[0].set_ylim(bottom=0.5)
    ax[0].set_xlabel("residual $M_G - $ fiducial (mag)")
    ax[0].set_ylabel("stars")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Residuals vs a Gaussian of equal width")

    # Folded comparison: positive vs mirrored negative tail
    k = np.linspace(0.5, 7, 60)
    npos = [(r > med + kk * sigma).sum() for kk in k]
    nneg = [(r < med - kk * sigma).sum() for kk in k]
    ax[1].plot(k, npos, color=POS, label="fainter than fiducial (deficit-like)")
    ax[1].plot(k, nneg, color=NEG, label="brighter than fiducial (control)")
    ax[1].set_yscale("log")
    ax[1].set_xlabel(r"threshold $k$ ($\times\sigma$)")
    ax[1].set_ylabel("count beyond threshold")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Tail asymmetry")
    return _save(fig, name)


def fig_floor(scaling: pd.DataFrame, random: pd.DataFrame, sigma: float,
              n_total: int, name="F3_systematic_floor") -> str:
    """THE headline figure: group-mean scatter vs group size stops falling."""
    fig, ax = plt.subplots(figsize=(6.6, 4.6))

    nn = np.logspace(1, np.log10(max(n_total, 100)), 100)
    ax.plot(nn, sigma / np.sqrt(nn), color="black", ls="--", lw=1.2,
            label=r"$\sigma/\sqrt{N}$ (the naive expectation)")

    ax.plot(random["n"], random["rms_of_means"], "o-", color=NEG, ms=4,
            label="random subsamples (control)")

    markers = {"sky (HEALPix)": "s", "A_0": "^", "phot_g_mean_mag": "v",
               "dist_pc": "D", "sky_density": "P"}
    colours = {"sky (HEALPix)": "#c0392b", "A_0": "#e67e22",
               "phot_g_mean_mag": "#8e44ad", "dist_pc": "#16a085",
               "sky_density": "#7f8c8d"}
    labels = {"sky (HEALPix)": "sky patches", "A_0": "extinction bins",
              "phot_g_mean_mag": "apparent-$G$ bins",
              "dist_pc": "distance bins", "sky_density": "crowding bins"}
    for axis, sub in scaling.groupby("axis"):
        sub = sub.sort_values("mean_group_n")
        ax.plot(sub["mean_group_n"], sub["rms_group_means"],
                marker=markers.get(axis, "o"), ms=4, lw=1.1,
                color=colours.get(axis, NEU), label=labels.get(axis, axis))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("stars per group $N$")
    ax.set_ylabel("RMS of group mean residuals (mag)")
    ax.set_title("Structured subsamples stop averaging down;\n"
                 "random ones do not. The plateau is the sensitivity.")
    ax.legend(fontsize=7, loc="lower left")
    return _save(fig, name)


def fig_splits(splits: pd.DataFrame, floor: float | None = None,
               name="F4_null_splits") -> str:
    """Claim: every split that must return zero does not."""
    s = splits.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    y = np.arange(len(s))
    colours = [POS if v > 0 else NEG for v in s["difference"]]
    ax.errorbar(s["difference"], y, xerr=s["difference_err"], fmt="o",
                ms=5, lw=1.2, ecolor=NEU, capsize=3,
                markerfacecolor="none", markeredgewidth=1.4)
    for i, (v, e, c) in enumerate(zip(s["difference"], s["difference_err"], colours)):
        ax.plot([v], [i], "o", color=c, ms=5)
    ax.axvline(0, color="black", lw=1)
    if floor:
        ax.axvspan(-floor, floor, color="grey", alpha=0.15,
                   label=f"measured floor $\\pm${floor:.3f} mag")
        ax.legend(fontsize=7)
    ax.set_yticks(y); ax.set_yticklabels(s["test"], fontsize=7.5)
    ax.set_xlabel("mean residual difference between halves (mag)")
    ax.set_title("Null split tests")
    return _save(fig, name)


def fig_leverage(tab: pd.DataFrame, slope: float, alpha_blind: float,
                 name="F5_spectral_leverage") -> str:
    """Claim: the blind spot is not the flat absorber."""
    from . import anchor
    a = np.linspace(-0.5, 4.0, 300)
    lev = np.array([anchor.leverage(x, slope) for x in a])

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(a, lev, color=NEU, lw=1.8)
    ax.axhline(0, color="black", lw=1)
    ax.axvline(alpha_blind, color=POS, ls="--", lw=1.2,
               label=fr"blind spot $\alpha={alpha_blind:.2f}$")
    ax.plot([0], [anchor.leverage(0.0, slope)], "o", color=NEG, ms=7,
            label=f"grey absorber: leverage {anchor.leverage(0.0, slope):+.2f}")
    ax.plot([2.0], [anchor.leverage(2.0, slope)], "s", color="#e67e22", ms=7,
            label=f"dust-like $\\alpha=2$: {anchor.leverage(2.0, slope):+.2f}")
    ax.fill_between(a, -1, 0, color=NEG, alpha=0.07)
    ax.fill_between(a, 0, 1, color=POS, alpha=0.07)
    ax.text(3.2, 0.45, "deficit-like\n(fainter)", color=POS, fontsize=8, ha="center")
    ax.text(-0.2, -0.5, "over-luminous", color=NEG, fontsize=8, ha="left")
    ax.set_xlabel(r"absorber spectral slope $\alpha$   ($\tau \propto \lambda^{-\alpha}$)")
    ax.set_ylabel(r"leverage $= (\delta m_G - s\,\delta m_{K_s})/\delta m_G$")
    ax.set_ylim(-1.0, 1.1)
    ax.set_title(f"Spectral leverage, fiducial slope $s={slope:.2f}$")
    ax.legend(fontsize=7.5, loc="lower right")
    return _save(fig, name)


def fig_injection(summary: pd.DataFrame, name="F6_injection_recovery") -> str:
    """Claim: what the pipeline can and cannot recover."""
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))

    uni = summary[summary["mode"] == "uniform"].sort_values("f_injected")
    if len(uni):
        ax[0].errorbar(uni["f_injected"], np.abs(uni["mean_residual"]),
                       yerr=uni["mean_residual_std"], fmt="o-", color=NEG,
                       ms=5, capsize=3, label="self-calibrated (mean residual)")
        ax[0].errorbar(uni["f_injected"], uni["anchored_f"],
                       yerr=uni["anchored_f_std"], fmt="s-", color=POS,
                       ms=5, capsize=3, label="model-anchored")
        ax[0].plot(uni["f_injected"], uni["f_injected"], color="black",
                   ls="--", lw=1, label="perfect recovery")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("injected uniform $f$")
    ax[0].set_ylabel("recovered signal")
    ax[0].set_title("Uniform injection:\nself-calibrated recovery is zero")
    ax[0].legend(fontsize=7)

    sp = summary[summary["mode"] == "sparse"]
    for f, sub in sp.groupby("f_injected"):
        sub = sub.sort_values("p_injected")
        ax[1].errorbar(sub["p_injected"], sub["p_recovered"],
                       yerr=sub["p_recovered_std"], fmt="o-", ms=4, capsize=3,
                       label=f"$f={f:g}$")
    lims = [sp["p_injected"].min() * 0.5, sp["p_injected"].max() * 2] if len(sp) else [1e-6, 1e-1]
    ax[1].plot(lims, lims, color="black", ls="--", lw=1, label="perfect")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("injected harvested fraction $p$")
    ax[1].set_ylabel("recovered $p$")
    ax[1].set_title("Sparse injection:\nrecovered in the tail, not the mean")
    ax[1].legend(fontsize=7)
    return _save(fig, name)


def fig_exclusion(excl: pd.DataFrame, name="F7_exclusion") -> str:
    """Claim: the (f, p) region excluded at 95% CL."""
    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    ax.plot(excl["f"], excl["p_upper_limit"], "o-", color=POS, ms=4)
    ax.fill_between(excl["f"], excl["p_upper_limit"], 1.0, color=POS, alpha=0.12)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("per-star harvested optical fraction $f$")
    ax.set_ylabel("upper limit on harvested fraction of stars $p$")
    ax.set_title("95% CL exclusion (shaded region excluded)")
    ax2 = ax.twinx()
    ax2.plot(excl["f"], excl["mean_f_upper_limit"], color=NEG, ls=":", lw=1.4)
    ax2.set_yscale("log")
    ax2.set_ylabel(r"implied limit on mean $\bar f = p f$", color=NEG)
    ax2.grid(False)
    return _save(fig, name)


def fig_cutflow(cf: pd.DataFrame, name="F8_cutflow") -> str:
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    y = np.arange(len(cf))
    ax.barh(y, cf["n_after"], color=NEU)
    ax.set_yticks(y); ax.set_yticklabels(cf["cut"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("stars surviving")
    ax.set_title("Cut flow")
    for i, v in enumerate(cf["n_after"]):
        ax.text(v * 1.05, i, f"{v:,}", va="center", fontsize=6.5)
    return _save(fig, name)
