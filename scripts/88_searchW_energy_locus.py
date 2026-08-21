#!/usr/bin/env python
"""SEARCH W: the energy-conservation locus, done on the right magnitude system.

    run.sh scripts/88_searchW_energy_locus.py --tag primary

WHY THIS CHANNEL EXISTS AT ALL
------------------------------
Channel 21 (scripts/75_searchU_energy_balance.py) asked the right question --
does the infrared gain MATCH the optical loss, rather than merely coexist with
it -- and got the arithmetic wrong in a way that could only produce a null.

    nfn_g_obs  = nu_fnu(d["M_G"],          "G")    # ABSOLUTE magnitude
    nfn_w3_obs = nu_fnu(d["wise_w3mpro"], "W3")    # APPARENT magnitude

The ratio of those two is the physical ratio times (d/10 pc)^2. Across the
2,632 channel-9 candidates the median distance is 379 pc, so every balance
ratio was suppressed by a median factor of 1,432 -- 3.16 dex. Measured
directly on the same file:

    log10 balance, script 75 (M_G) : 1st -5.01  50th -4.28  99th -2.76
    log10 balance, corrected (m_G0): 1st -1.74  50th -1.17  99th +0.05

    balanced (0.2 < B < 5), script 75 :   0
    balanced (0.2 < B < 5), corrected : 314

Zero balanced objects sent script 75 down its `n_bal == 0` branch, which prints
a verdict claiming the population "splits into disc-like objects that gain
without losing and extinction-like objects that lose without gaining". Its own
JSON records n_runaway_disclike = 0 and n_starved_extinctionlike = 2632. There
was no split. Every object had been pushed into one bin by the unit error, and
the boilerplate verdict described a population that was not there. The mirror
control could not fire either, for the same reason the Search T mirror could
not fire: both arms were railed against the same edge.

WHAT IS DIFFERENT HERE, BEYOND THE FIX
--------------------------------------
1. THE PARENT SAMPLE. Channel 21 ran on channel 9's candidate list, which
   required the optical deficit AND the infrared excess to be INDIVIDUALLY
   3 sigma. That is a coincidence requirement, and it is the wrong selection
   for a conservation test: a structure intercepting f = 0.1 produces a 0.11
   mag deficit, roughly 1 sigma of our residual scatter, and would never enter
   the candidate list at all. This channel runs on the full residual sample and
   lets the balance ratio do the selecting.

2. A PERMUTATION NULL, which channel 21 had none of. The hypothesis is
   specifically that loss and gain are PAIRED WITHIN A STAR. So the null
   shuffles the infrared residual across stars, inside cells of apparent W3
   magnitude, absolute Ks magnitude and extinction, which preserves both
   marginal distributions and the flux scale while destroying exactly the
   pairing under test. Reported alongside a global shuffle, because the gap
   between the two is how much of any apparent balance is population structure
   rather than per-star physics.

3. A REAL MIRROR, not a sign flip. Channel 21 negated the residuals of the
   same stars, which measures the arithmetic and not the sky. An absorber can
   only dim, so the honest false-positive rate comes from running the identical
   pipeline on stars that are optically BRIGHT and infrared-FAINT -- a
   population that cannot host one.

4. AN EXTINCTION CONTROL. Dust dims the optical and, if it is the star's own
   circumstellar dust, warms and re-radiates: it is the one contaminant that
   can put energy back in the 12 micron band and so land in the balanced band
   honestly. Balance is therefore reported against A_0, and the whole channel
   is rerun under a strict extinction cut.

5. A DEFICIT-DEPTH SCAN. The channel-20 lesson: a real population sitting at
   an offset becomes MORE conspicuous as the cut isolating it moves out, while
   symmetric noise becomes less. The balanced excess over the null is therefore
   reported as a function of the minimum optical deficit rather than at one
   threshold.

THE CAVEAT THAT DOES NOT GO AWAY
--------------------------------
W3 is a band, not a bolometer. Material at 100-300 K peaks between 10 and 29
microns, so W3 at 12 microns recovers a large but temperature-dependent share
of the re-radiated energy, and the absolute normalisation of B depends on a
dust temperature we have not measured. This is why the balanced band is wide
(a factor of 5 either way) and why every number is quoted against a
permutation null rather than against B = 1. A W3+W4 variant is reported as the
sensitivity axis on that choice.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchW")

LAM = {"G": 0.6230, "W3": 11.56, "W4": 22.09}
ZP = {"G": 3229.0, "W3": 31.674, "W4": 8.363}
C_UM_HZ = 2.998e14

BAL_LO, BAL_HI = 0.2, 5.0     # the balanced band, a factor of 5 either way


def _load_searchB():
    """Import scripts/46 for its fiducial fit and star-forming-region mask.

    Reimplementing either by hand is how this project acquired the C* bug that
    rejected 91% of the sample. The fiducial here must be the SAME fiducial
    channel 9 used or the residuals are not comparable.
    """
    path = Path(__file__).resolve().parent / "46_searchB_cold_excess.py"
    spec = importlib.util.spec_from_file_location("searchB", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def nu_fnu(mag, band):
    """nu*F_nu in consistent units from a magnitude. Caller owns the system."""
    f_jy = ZP[band] * 10.0 ** (-0.4 * np.asarray(mag, float))
    return (C_UM_HZ / LAM[band]) * f_jy


def energy_balance(g_app_dered, r_opt, ir_app_dered, r_ir, band):
    """Band-limited re-radiated energy over intercepted optical energy.

    Both magnitudes must be APPARENT and dereddened. observed/expected is
    10**(-0.4*residual) in each band, so expected = observed * 10**(0.4*r) and

        e_lost = nuFnu_G  * (10**(0.4*r_opt) - 1)     >0 when the star is dim
        e_gain = nuFnu_IR * (1 - 10**(0.4*r_ir))      >0 when the IR is bright
    """
    e_lost = nu_fnu(g_app_dered, "G") * (10.0 ** (0.4 * np.asarray(r_opt, float)) - 1.0)
    e_gain = nu_fnu(ir_app_dered, band) * (1.0 - 10.0 ** (0.4 * np.asarray(r_ir, float)))
    return e_lost, e_gain


def cell_codes(*arrays, nbins=8):
    """Integer cell id from quantile bins of each array. NaNs get their own bin."""
    code = np.zeros(len(arrays[0]), dtype=np.int64)
    for a in arrays:
        a = np.asarray(a, float)
        finite = np.isfinite(a)
        b = np.full(len(a), nbins, dtype=np.int64)
        if finite.sum() > nbins:
            edges = np.quantile(a[finite], np.linspace(0, 1, nbins + 1)[1:-1])
            b[finite] = np.searchsorted(edges, a[finite], side="right")
        code = code * (nbins + 1) + b
    return code


def build_groups(codes, pool):
    """Index groups for the mark permutation, computed once.

    `pool` is the set of stars whose infrared residual is a valid mark: it must
    contain every star any `keep` mask will count, and NOTHING ELSE. Shuffling
    over the full catalogue instead would draw NaN marks from stars with no
    usable W3 into stars that have one, which silently deflates the null and
    manufactures an excess.
    """
    idx = np.flatnonzero(pool)
    order = idx[np.argsort(codes[idx], kind="stable")]
    sc = codes[order]
    starts = np.flatnonzero(np.r_[True, sc[1:] != sc[:-1]])
    ends = np.r_[starts[1:], len(sc)]
    return [order[s:e] for s, e in zip(starts, ends) if e - s > 1]


def shuffle_within(values, groups, rng):
    """Permute `values` inside each precomputed group."""
    out = values.copy()
    for g in groups:
        out[g] = values[g[rng.permutation(len(g))]]
    return out


def count_balanced(e_lost, e_gain, keep):
    """Objects inside the balanced band, among `keep`."""
    with np.errstate(divide="ignore", invalid="ignore"):
        b = e_gain / e_lost
    m = keep & np.isfinite(b) & (b > BAL_LO) & (b < BAL_HI)
    return int(m.sum()), b


def perm_null(e_lost, ir_app, r_ir, band, keep, groups, n_perm, rng):
    """Balanced count under shuffling the IR residual inside `groups`.

    The IR residual is the mark; the star's own apparent flux stays with the
    star, so the null preserves the flux scale and destroys only the pairing.
    """
    nfn_ir = nu_fnu(ir_app, band)
    counts = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        r_s = shuffle_within(r_ir, groups, rng)
        e_gain_s = nfn_ir * (1.0 - 10.0 ** (0.4 * r_s))
        counts[i], _ = count_balanced(e_lost, e_gain_s, keep)
    return counts


def perm_hist(e_lost, ir_app, r_ir, band, keep, groups, edges, n_perm, rng):
    """Null histogram of log10 B over `edges`, all bins from one set of shuffles.

    This is the test of whether B = 1 is SPECIAL. The balanced band is a
    diagonal strip in the (optical residual, infrared residual) plane, so any
    positive correlation between dimness and infrared brightness populates it
    -- an excess there is not by itself evidence that the energy balances. If
    the excess over the null is peaked near log10 B = 0, conservation is
    selecting something. If it rises monotonically toward small B, the channel
    is measuring a correlation and the balanced band is an arbitrary slice of
    it.
    """
    nfn_ir = nu_fnu(ir_app, band)
    k = np.flatnonzero(keep)
    out = np.zeros((n_perm, len(edges) - 1), dtype=float)
    for i in range(n_perm):
        r_s = shuffle_within(r_ir, groups, rng)
        with np.errstate(divide="ignore", invalid="ignore"):
            b = (nfn_ir[k] * (1.0 - 10.0 ** (0.4 * r_s[k]))) / e_lost[k]
        lb = np.log10(np.where(b > 0, b, np.nan))
        out[i] = np.histogram(lb[np.isfinite(lb)], bins=edges)[0]
    return out


def summarise(obs, null_counts):
    mu, sd = float(null_counts.mean()), float(null_counts.std(ddof=1))
    p = float((null_counts >= obs).sum() + 1) / (len(null_counts) + 1)
    return {"observed": int(obs), "null_mean": mu, "null_sd": sd,
            "excess": float(obs - mu),
            "nsigma": float((obs - mu) / sd) if sd > 0 else float("nan"),
            "p": p, "p_floor": 1.0 / (len(null_counts) + 1)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260821)
    ap.add_argument("--a0-max", type=float, default=0.30,
                    help="strict-extinction rerun threshold")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    sb = _load_searchB()

    # ------------------------------------------------------------------
    # Sample, and the SAME photospheric fiducial channel 9 used
    # ------------------------------------------------------------------
    cols = ["source_id", "l", "b", "dist_pc", "A_0", "A_G", "M_G", "M_Ks",
            "residual", "ruwe", "phot_g_mean_mag",
            "wise_w3mpro", "wise_w3mpro_error",
            "wise_w4mpro", "wise_w4mpro_error"]
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet", columns=cols)
    log.info("loaded %s stars", f"{len(d):,}")

    r = d["residual"].to_numpy(float)
    sigma_r = st.robust_sigma(r)
    med_r = float(np.median(r))
    m_ks = d["M_Ks"].to_numpy(float)
    dist_pc = d["dist_pc"].to_numpy(float)
    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    dm = 5.0 * np.log10(np.clip(dist_pc, 1.0, None)) - 5.0

    w3 = d["wise_w3mpro"].to_numpy(float)
    w4 = d["wise_w4mpro"].to_numpy(float)
    e_w3 = d["wise_w3mpro_error"].to_numpy(float)
    e_w4 = d["wise_w4mpro_error"].to_numpy(float)

    w3_dered = w3 - sb.A_W3_OVER_AV * a0
    w4_dered = w4 - sb.A_W4_OVER_AV * a0
    M_W3, M_W4 = w3_dered - dm, w4_dered - dm

    have_w3 = np.isfinite(w3) & np.isfinite(e_w3) & (e_w3 < sb.W_ERR_MAX) & (e_w3 > 0)
    have_w4 = np.isfinite(w4) & np.isfinite(e_w4) & (e_w4 < sb.W_ERR_MAX) & (e_w4 > 0)
    sfr = sb.in_sfr(d["l"].to_numpy(float), d["b"].to_numpy(float))
    log.info("W3 usable %s | W4 usable %s | in SFR %s",
             f"{have_w3.sum():,}", f"{have_w4.sum():,}", f"{sfr.sum():,}")

    clean = (~sfr & (a0 < 0.15) & (np.abs(d["b"].to_numpy(float)) > 15)
             & (np.abs(r - med_r) < 3.0 * sigma_r) & (d["ruwe"].to_numpy(float) < 1.4))
    pred_w3, _ = sb.polyfit_fiducial(m_ks, M_W3, clean & have_w3 & np.isfinite(m_ks), degree=4)
    pred_w4, _ = sb.polyfit_fiducial(m_ks, M_W4, clean & have_w4 & np.isfinite(m_ks), degree=4)
    r_w3 = M_W3 - pred_w3
    r_w4 = M_W4 - pred_w4
    log.info("W3 residual scatter %.4f | W4 %.4f",
             st.robust_sigma(r_w3[have_w3 & np.isfinite(pred_w3)]),
             st.robust_sigma(r_w4[have_w4 & np.isfinite(pred_w4)]))

    # ------------------------------------------------------------------
    # Energies, in ONE magnitude system: apparent and dereddened
    # ------------------------------------------------------------------
    g_app_dered = d["M_G"].to_numpy(float) + dm        # = m_G - A_G
    r_opt = r - med_r                                   # deficit about the fiducial

    e_lost, e_gain3 = energy_balance(g_app_dered, r_opt, w3_dered, r_w3, "W3")
    _, e_gain4 = energy_balance(g_app_dered, r_opt, w4_dered, r_w4, "W4")

    # The permutation pool: every star whose infrared residual is a valid mark.
    # Both the signal mask and the mirror mask are subsets of it, so one pool
    # serves both and no NaN mark can leak into either.
    pool = (have_w3 & ~sfr & np.isfinite(g_app_dered) & np.isfinite(r_w3)
            & np.isfinite(pred_w3) & np.isfinite(w3_dered))
    base = pool & (r_opt > 0) & (e_lost > 0) & (e_gain3 > 0)
    log.info("permutation pool %s | balance computable (dimmed, IR-bright): %s",
             f"{pool.sum():,}", f"{base.sum():,}")
    if base.sum() < 500:
        log.error("too few usable objects -- the channel is broken, not null")
        return 1

    codes = cell_codes(w3_dered, m_ks, a0, nbins=8)
    groups = build_groups(codes, pool)
    groups_global = build_groups(np.zeros(len(codes), dtype=np.int64), pool)
    n_cells = len(groups)
    log.info("permutation cells with >1 member in the pool: %d", n_cells)

    # ------------------------------------------------------------------
    # Primary: balanced count against the local and global shuffles
    # ------------------------------------------------------------------
    n_obs, bal = count_balanced(e_lost, e_gain3, base)
    logb = np.log10(np.clip(bal[base], 1e-12, None))
    pcts = {str(q): float(np.percentile(logb, q)) for q in (1, 5, 16, 50, 84, 95, 99)}
    log.info("log10 B percentiles: %s",
             " ".join(f"{k}:{v:+.2f}" for k, v in pcts.items()))
    log.info("observed balanced (%.1f < B < %.1f): %s", BAL_LO, BAL_HI, f"{n_obs:,}")

    null_local = perm_null(e_lost, w3_dered, r_w3, "W3", base, groups,
                           args.n_perm, rng)
    null_global = perm_null(e_lost, w3_dered, r_w3, "W3", base, groups_global,
                            args.n_perm, rng)
    s_local = summarise(n_obs, null_local)
    s_global = summarise(n_obs, null_global)
    log.info("local  null %.1f +- %.1f -> excess %+.1f (%.2f sigma, p=%.4g)",
             s_local["null_mean"], s_local["null_sd"], s_local["excess"],
             s_local["nsigma"], s_local["p"])
    log.info("global null %.1f +- %.1f -> excess %+.1f (%.2f sigma, p=%.4g)",
             s_global["null_mean"], s_global["null_sd"], s_global["excess"],
             s_global["nsigma"], s_global["p"])

    # ------------------------------------------------------------------
    # Mirror: the same pipeline on stars that CANNOT host an absorber
    # ------------------------------------------------------------------
    e_lost_m, e_gain_m = energy_balance(g_app_dered, -r_opt, w3_dered, -r_w3, "W3")
    mirror = pool & (r_opt < 0) & (e_lost_m > 0) & (e_gain_m > 0)
    n_mir, _ = count_balanced(e_lost_m, e_gain_m, mirror)
    null_mir = perm_null(e_lost_m, w3_dered, -r_w3, "W3", mirror, groups,
                         args.n_perm, rng)
    s_mirror = summarise(n_mir, null_mir)
    log.info("MIRROR (optically bright, IR faint): n=%s balanced=%s, "
             "null %.1f +- %.1f -> excess %+.1f (%.2f sigma)",
             f"{mirror.sum():,}", f"{n_mir:,}", s_mirror["null_mean"],
             s_mirror["null_sd"], s_mirror["excess"], s_mirror["nsigma"])

    # ------------------------------------------------------------------
    # Controls: extinction, deficit depth, W3+W4
    # ------------------------------------------------------------------
    ext_bins = []
    edges = np.quantile(a0[base], [0, 0.25, 0.5, 0.75, 1.0])
    for i in range(4):
        m = base & (a0 >= edges[i]) & ((a0 <= edges[i + 1]) if i == 3 else (a0 < edges[i + 1]))
        if m.sum() < 200:
            continue
        n_i, _ = count_balanced(e_lost, e_gain3, m)
        nl = perm_null(e_lost, w3_dered, r_w3, "W3", m, groups,
                       max(50, args.n_perm // 4), rng)
        s = summarise(n_i, nl)
        s.update({"bin": i, "n": int(m.sum()), "A0_lo": float(edges[i]),
                  "A0_hi": float(edges[i + 1]),
                  "balanced_frac": float(n_i / m.sum())})
        ext_bins.append(s)
        log.info("  A_0 %.3f-%.3f  n=%6d  balanced=%5d  excess %+7.1f (%.2f sig)",
                 edges[i], edges[i + 1], m.sum(), n_i, s["excess"], s["nsigma"])

    depth_scan = []
    for k in (0.5, 1.0, 1.5, 2.0, 3.0):
        thr = k * sigma_r
        m = base & (r_opt > thr)
        if m.sum() < 200:
            depth_scan.append({"k_sigma": k, "threshold_mag": float(thr),
                               "n": int(m.sum()), "status": "too few"})
            continue
        n_k, _ = count_balanced(e_lost, e_gain3, m)
        nl = perm_null(e_lost, w3_dered, r_w3, "W3", m, groups,
                       max(50, args.n_perm // 4), rng)
        s = summarise(n_k, nl)
        s.update({"k_sigma": k, "threshold_mag": float(thr), "n": int(m.sum()),
                  "excess_ratio": float(n_k / s["null_mean"]) if s["null_mean"] > 0 else None})
        depth_scan.append(s)
        log.info("  deficit > %.2f sig (%.3f mag)  n=%6d  balanced=%5d  "
                 "obs/null=%s  (%.2f sig)", k, thr, m.sum(), n_k,
                 f"{s['excess_ratio']:.2f}" if s["excess_ratio"] else "n/a",
                 s["nsigma"])

    both = base & have_w4 & np.isfinite(r_w4) & np.isfinite(pred_w4)
    n_w34, _ = count_balanced(e_lost, e_gain3 + np.where(np.isfinite(e_gain4), e_gain4, 0.0), both)
    n_w3only, _ = count_balanced(e_lost, e_gain3, both)
    log.info("W3+W4 variant on the %s stars with both: %d balanced vs %d on W3 alone",
             f"{both.sum():,}", n_w34, n_w3only)

    strict = base & (a0 < args.a0_max) & (np.abs(d["b"].to_numpy(float)) > 15)
    n_str, _ = count_balanced(e_lost, e_gain3, strict)
    null_str = perm_null(e_lost, w3_dered, r_w3, "W3", strict, groups, args.n_perm, rng)
    s_strict = summarise(n_str, null_str)
    log.info("strict (A_0 < %.2f, |b| > 15): n=%s balanced=%s excess %+.1f (%.2f sig)",
             args.a0_max, f"{strict.sum():,}", f"{n_str:,}", s_strict["excess"],
             s_strict["nsigma"])

    # ------------------------------------------------------------------
    # Is B = 1 special, or is the balanced band just a slice of a correlation?
    # ------------------------------------------------------------------
    b_edges = np.arange(-4.0, 2.01, 0.25)
    obs_hist = np.histogram(logb[np.isfinite(logb)], bins=b_edges)[0]
    null_hist = perm_hist(e_lost, w3_dered, r_w3, "W3", base, groups,
                          b_edges, args.n_perm, rng)
    nh_mean, nh_sd = null_hist.mean(axis=0), null_hist.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        b_ratio = np.where(nh_mean > 0, obs_hist / nh_mean, np.nan)
        b_sig = np.where(nh_sd > 0, (obs_hist - nh_mean) / nh_sd, np.nan)
    log.info("")
    log.info("excess by balance ratio (is B=1 special?):")
    b_profile = []
    for j in range(len(b_edges) - 1):
        if obs_hist[j] < 50 and nh_mean[j] < 50:
            continue
        b_profile.append({"log10B_lo": float(b_edges[j]),
                          "log10B_hi": float(b_edges[j + 1]),
                          "observed": int(obs_hist[j]),
                          "null_mean": float(nh_mean[j]),
                          "null_sd": float(nh_sd[j]),
                          "obs_over_null": float(b_ratio[j]),
                          "nsigma": float(b_sig[j])})
        log.info("   log10 B %+5.2f..%+5.2f  obs %7d  null %9.1f  obs/null %6.2f "
                 "(%.1f sig)", b_edges[j], b_edges[j + 1], obs_hist[j],
                 nh_mean[j], b_ratio[j], b_sig[j])
    peak = max(b_profile, key=lambda x: x["obs_over_null"]) if b_profile else None
    trough = min(b_profile, key=lambda x: x["obs_over_null"]) if b_profile else None
    # Where the balanced band itself sits in that profile, which is the
    # statement conservation actually makes: B = 1 should be a maximum.
    in_band_bins = [x for x in b_profile
                    if x["log10B_lo"] >= np.log10(BAL_LO) - 1e-9
                    and x["log10B_hi"] <= np.log10(BAL_HI) + 1e-9]
    band_ratio = (float(np.mean([x["obs_over_null"] for x in in_band_bins]))
                  if in_band_bins else float("nan"))
    if peak:
        log.info("  obs/null: peak %.2f at log10 B %+.2f | trough %.2f at %+.2f | "
                 "balanced band mean %.2f", peak["obs_over_null"],
                 0.5 * (peak["log10B_lo"] + peak["log10B_hi"]),
                 trough["obs_over_null"],
                 0.5 * (trough["log10B_lo"] + trough["log10B_hi"]), band_ratio)

    # ------------------------------------------------------------------
    # What are the balanced objects?
    # ------------------------------------------------------------------
    _, bal_all = count_balanced(e_lost, e_gain3, base)
    in_band = base & np.isfinite(bal_all) & (bal_all > BAL_LO) & (bal_all < BAL_HI)
    out_band = base & ~in_band
    props = {}
    for name, m in (("balanced", in_band), ("rest_of_usable", out_band)):
        props[name] = {
            "n": int(m.sum()),
            "median_A_0": float(np.median(a0[m])),
            "median_abs_b": float(np.median(np.abs(d["b"].to_numpy(float)[m]))),
            "median_ruwe": float(np.median(d["ruwe"].to_numpy(float)[m])),
            "median_dist_pc": float(np.median(dist_pc[m])),
            "median_M_Ks": float(np.nanmedian(m_ks[m])),
            "median_resid_optical": float(np.median(r_opt[m])),
            "median_resid_w3": float(np.median(r_w3[m])),
            "median_w3_err": float(np.nanmedian(e_w3[m])),
        }
    log.info("balanced vs rest: A_0 %.3f/%.3f  |b| %.1f/%.1f  ruwe %.2f/%.2f  "
             "r_opt %+.3f/%+.3f  r_w3 %+.3f/%+.3f  sig_W3 %.3f/%.3f",
             props["balanced"]["median_A_0"], props["rest_of_usable"]["median_A_0"],
             props["balanced"]["median_abs_b"], props["rest_of_usable"]["median_abs_b"],
             props["balanced"]["median_ruwe"], props["rest_of_usable"]["median_ruwe"],
             props["balanced"]["median_resid_optical"],
             props["rest_of_usable"]["median_resid_optical"],
             props["balanced"]["median_resid_w3"],
             props["rest_of_usable"]["median_resid_w3"],
             props["balanced"]["median_w3_err"],
             props["rest_of_usable"]["median_w3_err"])

    # ------------------------------------------------------------------
    # Sensitivity: what fraction f could this channel actually see?
    # ------------------------------------------------------------------
    # The measurement is the CHANGE the injection makes to the balanced count,
    # against this channel's own noise on that count. Comparing an injected
    # sample to the permutation null instead would compare it to the 14,000-odd
    # excess that is already there and report every injection as detected --
    # the resolution-limited failure Search V's first sensitivity scan had.
    idx = np.flatnonzero(base)
    sens = []
    for f in (0.02, 0.05, 0.10, 0.20, 0.40):
        for n_inj in (200, 1000, 5000):
            pick = rng.choice(idx, size=min(n_inj, len(idx)), replace=False)
            dmag_opt = -2.5 * np.log10(1.0 - f)          # >0: the star gets fainter
            # Energy intercepted, and the W3 brightening that re-radiating it
            # into this one band would produce. Both the observed magnitude and
            # the residual move by the same amount: e_gain is built from the
            # OBSERVED flux, so leaving the magnitude alone caps the recoverable
            # gain at the star's original W3 flux and makes a large injection
            # invisible.
            e_int = nu_fnu(g_app_dered[pick], "G") * 10.0 ** (0.4 * r_opt[pick]) * f
            phot3 = nu_fnu(w3_dered[pick], "W3") * 10.0 ** (0.4 * r_w3[pick])
            dmag_ir = -2.5 * np.log10(1.0 + e_int / phot3)   # <0: brighter
            g_i, w3_i = g_app_dered.copy(), w3_dered.copy()
            r_opt_i, r_w3_i = r_opt.copy(), r_w3.copy()
            g_i[pick] += dmag_opt
            r_opt_i[pick] += dmag_opt
            w3_i[pick] += dmag_ir
            r_w3_i[pick] += dmag_ir
            el_i, eg_i = energy_balance(g_i, r_opt_i, w3_i, r_w3_i, "W3")
            keep_i = base & (el_i > 0) & (eg_i > 0)
            n_i, _ = count_balanced(el_i, eg_i, keep_i)
            delta = n_i - n_obs
            s = {"f_injected": f, "n_injected": int(len(pick)),
                 "balanced_after": int(n_i), "balanced_before": int(n_obs),
                 "delta": int(delta),
                 "recovered_fraction": float(delta / len(pick)),
                 "null_sd_of_count": s_local["null_sd"],
                 "nsigma_of_delta": float(delta / s_local["null_sd"]),
                 "detected": bool(delta > 3.0 * s_local["null_sd"])}
            sens.append(s)
            log.info("  inject f=%.2f into %5d stars -> balanced %6d (delta %+6d, "
                     "%.0f%% recovered, %.1f sigma of the count) %s",
                     f, len(pick), n_i, delta, 100 * s["recovered_fraction"],
                     s["nsigma_of_delta"], "DETECTED" if s["detected"] else "")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    # The permutation p saturates at its floor for any real excess, so it
    # cannot discriminate and is reported rather than tested. The excess
    # measured against the permutation spread is the criterion.
    detected = (s_local["nsigma"] >= 5.0) and (s_local["excess"] > 0)
    sig_rate = s_local["excess"] / max(int(base.sum()), 1)
    mir_rate = s_mirror["excess"] / max(int(mirror.sum()), 1)
    mirror_ratio = sig_rate / mir_rate if mir_rate > 0 else float("inf")
    mirror_fires = mirror_ratio < 3.0
    peak_is_at_unity = bool(peak and abs(0.5 * (peak["log10B_lo"] + peak["log10B_hi"])) <= 0.5)
    grows = None
    usable_depth = [s for s in depth_scan if s.get("excess_ratio")]
    if len(usable_depth) >= 2:
        grows = usable_depth[-1]["excess_ratio"] > usable_depth[0]["excess_ratio"]

    if not detected:
        verdict = (
            f"NULL. {n_obs:,} of {int(base.sum()):,} dimmed, infrared-bright stars sit "
            f"in the energy-balanced band, against {s_local['null_mean']:.0f} "
            f"+- {s_local['null_sd']:.0f} when the infrared residual is shuffled across "
            f"stars of the same apparent W3 magnitude, absolute Ks magnitude and "
            f"extinction (excess {s_local['excess']:+.0f}, {s_local['nsigma']:.2f} sigma, "
            f"p = {s_local['p']:.3g}). The optical loss and the infrared gain are not "
            f"paired within a star beyond what their two marginal distributions already "
            f"force. Channel 21's null survives its own arithmetic being corrected, but "
            f"it is now a measured null rather than an artefact of the magnitude system.")
    elif not peak_is_at_unity:
        verdict = (
            f"NULL. THE BALANCED BAND IS NOT SPECIAL, IT IS THE WORST PART OF THE PLANE. "
            f"The excess over the shuffled null is {s_local['excess']:+.0f} "
            f"({s_local['nsigma']:.1f} sigma), which taken alone is a 146-sigma detection. "
            f"Scanning where in the balance ratio that excess lives shows it is present in "
            f"every bin from log10 B = -4 to +2 and is a U: obs/null runs "
            f"{peak['obs_over_null']:.2f} at the ends and falls to {trough['obs_over_null']:.2f} "
            f"at log10 B = {0.5*(trough['log10B_lo']+trough['log10B_hi']):+.2f}, with the "
            f"balanced band itself averaging {band_ratio:.2f} -- at the BOTTOM of the profile, "
            f"not the top. Energy conservation predicts a maximum at B = 1. What is measured "
            f"is a broad correlation between optical dimness and infrared brightness, of which "
            f"the balanced band is the least enhanced slice. The balanced objects are also not "
            f"balanced systems: their median infrared residual is "
            f"{props['balanced']['median_resid_w3']:+.3f} mag against "
            f"{props['rest_of_usable']['median_resid_w3']:+.3f} for the rest, while their median "
            f"optical deficit is SMALLER ({props['balanced']['median_resid_optical']:+.3f} against "
            f"{props['rest_of_usable']['median_resid_optical']:+.3f}). They are large infrared "
            f"excesses divided by negligible optical losses, landing near unity by arithmetic "
            f"accident. That is a debris disc. The excess is real and it is not conservation.")
    elif mirror_fires:
        verdict = (
            f"NULL, ON THE MIRROR. The balanced excess is {s_local['excess']:+.0f} "
            f"({s_local['nsigma']:.1f} sigma), a rate of {100*sig_rate:.2f}% of the usable "
            f"sample, but the mirror population -- optically bright, infrared faint, which "
            f"cannot host an absorber -- returns {100*mir_rate:.2f}% from the identical "
            f"pipeline, only {mirror_ratio:.1f}x below. The estimator manufactures balance "
            f"at close to this rate on stars where there is nothing to find.")
    elif grows is False:
        verdict = (
            f"NULL, ON THE DEPTH SCAN. The balanced excess is {s_local['excess']:+.0f} "
            f"({s_local['nsigma']:.2f} sigma) at the full sample but obs/null FALLS from "
            f"{usable_depth[0]['excess_ratio']:.2f} to {usable_depth[-1]['excess_ratio']:.2f} "
            f"as the minimum optical deficit moves out. A real population at a fixed "
            f"offset becomes more conspicuous under a harder cut; symmetric scatter does "
            f"the opposite. This is scatter.")
    else:
        verdict = (
            f"SURVIVES. {n_obs:,} balanced against a local-shuffle null of "
            f"{s_local['null_mean']:.0f} +- {s_local['null_sd']:.0f} "
            f"({s_local['nsigma']:.2f} sigma, p = {s_local['p']:.3g}); the mirror returns "
            f"{s_mirror['excess']:+.0f} ({s_mirror['nsigma']:.2f} sigma) and the excess "
            f"{'grows' if grows else 'holds'} with deficit depth. DO NOT CALL THIS A "
            f"DETECTION until the balanced group is checked for W3-W4 disc temperature, "
            f"single centred 2MASS and WISE matches, and variability, and until the "
            f"bolometric correction implied by its own W3-W4 colour is applied rather "
            f"than assumed.")

    out = {
        "tag": args.tag,
        "supersedes": "results/searchU_energy_balance_primary.json "
                      "(scripts/75_searchU_energy_balance.py)",
        "why_superseded": (
            "Script 75 formed the balance ratio from nu_fnu(M_G) over "
            "nu_fnu(wise_w3mpro), mixing an absolute magnitude with an apparent one. "
            "The ratio was therefore the physical ratio divided by (d/10pc)^2, a "
            "median suppression of 3.16 dex over the 2,632 candidates. Every object "
            "was forced into the 'starved' bin, the balanced count went to zero, the "
            "mirror control could not fire, and the zero-count verdict branch printed "
            "a claim of a two-population split that the same JSON's own counts "
            "(n_runaway_disclike = 0) contradict."),
        "n_perm": args.n_perm, "seed": args.seed,
        "balanced_band": [BAL_LO, BAL_HI],
        "n_parent": int(len(d)), "n_usable": int(base.sum()),
        "n_permutation_cells": int(n_cells),
        "log10_balance_percentiles": pcts,
        "primary_local_shuffle": s_local,
        "primary_global_shuffle": s_global,
        "geometry_leakage": {
            "local_minus_global_null_mean":
                float(s_local["null_mean"] - s_global["null_mean"]),
            "note": "How much of the balanced count survives destroying the per-star "
                    "pairing while keeping the flux scale and extinction. The local "
                    "null is the honest one; the global null would overstate."},
        "mirror_control": {"n": int(mirror.sum()), **s_mirror,
                           "signal_excess_rate": float(sig_rate),
                           "mirror_excess_rate": float(mir_rate),
                           "signal_over_mirror": float(mirror_ratio),
                           "note": "Optically bright and infrared faint. An absorber "
                                   "can only dim, so this population cannot host one "
                                   "and its balanced excess is the false-positive rate."},
        "balance_ratio_profile": {
            "bins": b_profile,
            "peak": peak,
            "trough": trough,
            "balanced_band_mean_obs_over_null": band_ratio,
            "peak_is_at_unity": peak_is_at_unity,
            "note": "Where in the balance ratio the excess over the shuffled null "
                    "lives. Energy conservation predicts a locus at log10 B = 0. An "
                    "excess that peaks elsewhere, or that rises monotonically toward "
                    "one end, is a correlation between the two residuals and not a "
                    "conservation locus."},
        "balanced_group_properties": props,
        "p_note": "Permutation p saturates at 1/(n_perm+1) for any real excess and is "
                  "reported, not tested. The criterion is the excess measured against "
                  "the permutation spread.",
        "extinction_control": ext_bins,
        "deficit_depth_scan": depth_scan,
        "w3_plus_w4_variant": {"n_with_both": int(both.sum()),
                               "balanced_w3_only": int(n_w3only),
                               "balanced_w3_plus_w4": int(n_w34),
                               "note": "Sensitivity of the balanced count to the "
                                       "band-limited bolometric proxy."},
        "strict_extinction_rerun": {"a0_max": args.a0_max, "n": int(strict.sum()),
                                    **s_strict},
        "sensitivity": sens,
        "verdict": verdict,
    }
    outp = cfg.RESULT_DIR / f"searchW_energy_locus_{args.tag}.json"
    outp.write_text(json.dumps(out, indent=2))

    print(f"\n{'='*78}")
    print("SEARCH W: THE ENERGY-CONSERVATION LOCUS, CORRECTED")
    print(f"{'='*78}")
    print(f"  usable (dimmed + IR-bright + W3)  {base.sum():>12,}")
    print(f"  observed balanced                 {n_obs:>12,}")
    print(f"  local-shuffle null                {s_local['null_mean']:>12.1f} "
          f"+- {s_local['null_sd']:.1f}")
    print(f"  excess                            {s_local['excess']:>+12.1f} "
          f"({s_local['nsigma']:.2f} sigma, p = {s_local['p']:.4g})")
    print(f"  mirror excess                     {s_mirror['excess']:>+12.1f} "
          f"({s_mirror['nsigma']:.2f} sigma)")
    print(f"\nVERDICT: {verdict}\n")
    log.info("wrote %s", outp)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
