#!/usr/bin/env python
"""Why are variable stars one-sidedly dim? Test the mean-vs-epoch bias.

    run.sh scripts/68_variable_deficit_test.py --tag primary

WHAT SURVIVED
-------------
Of the four discarded classes Search M flagged, three collapsed toward the
reference once the 2MASS cross-match was required to be single, uncontested
and centred better than 0.3 arcsec: C* failures 98 -> 8.8, poor-Ks 70 -> 14.6,
extragalactic 60 -> 33, against a reference of 6.0.

The variables did not. They went from 180:1 to 218:1 -- a factor of 36 above
the reference, in stars whose photometry is otherwise pristine. That class was
the single largest cut in the whole project, 382,319 sources removed before
any channel ran, and Wright et al. (2016) name aperodic variability as the
primary predicted signature of a transiting megastructure.

So it has to be explained properly rather than waved at.

THE MECHANISM TO EXCLUDE FIRST
------------------------------
The estimator compares two magnitudes measured in different ways:

    Gaia G   is a MEAN over ~34 months of scanning.
    2MASS Ks is a SINGLE EPOCH, taken once between 1997 and 2001.

For a constant star that difference is harmless. For a variable it is not, and
crucially the bias is ONE-SIGNED whenever the light curve is asymmetric.
A star that spends most of its time near maximum and occasionally dips --
eclipsing binaries, transiting systems, dippers, anything occulted -- has a
mean fainter than a randomly-sampled epoch. Averaging G over the dips while
catching Ks out of them produces exactly a G-band deficit, with no absorber
involved.

That mechanism makes a sharp, falsifiable prediction: the deficit must scale
with variability amplitude. A star varying at the millimagnitude level cannot
produce a tenth-magnitude offset no matter what shape its light curve has.

THE AMPLITUDE PROXY
-------------------
Gaia does not publish an amplitude for every source, but it publishes enough
to reconstruct one. phot_g_mean_flux_over_error is the mean flux divided by
its standard error, and the standard error is the epoch-to-epoch scatter over
sqrt(N). So

    amplitude ~ sqrt(N_obs) / (flux / sigma_flux)

is the fractional epoch scatter, which is the standard Gaia variability proxy
(used this way by Belokurov+2017, Deason+2017). G transits are counted per CCD
and BP/RP per field-of-view, nine CCDs to a field, so N_G is taken as
9 x phot_bp_n_obs where the G count is absent.
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
from pipeline import sample as smp
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("vartest")

COLS = [
    "source_id", "l", "b", "parallax", "parallax_error",
    "phot_g_mean_mag", "phot_g_mean_flux_over_error",
    "phot_bp_n_obs", "phot_rp_n_obs",
    "bp_rp", "phot_bp_rp_excess_factor", "phot_variable_flag",
    "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat",
    "astrometric_params_solved",
    "tmass_ks_m", "tmass_ph_qual",
    "tmass_xm_dist", "tmass_xm_nnb", "tmass_xm_nmates",
]

K_SIGMA = 3.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=K_SIGMA)
    args = ap.parse_args()

    import pyarrow.parquet as pq
    files = sorted(cfg.RAW_DIR.glob("sample_d500_p*.parquet"))
    frames = []
    for f in files:
        have = set(pq.read_schema(f).names)
        frames.append(pd.read_parquet(f, columns=[c for c in COLS if c in have]))
    d = pd.concat(frames, ignore_index=True)
    log.info("raw rows: %d", len(d))

    d = smp.add_astrometry(d)
    a0 = ext.query_a0("edenhofer23", d["l"].to_numpy(float),
                      d["b"].to_numpy(float), d["dist_pc"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = d["dist_mod"].to_numpy(float)
    M_G = d["phot_g_mean_mag"].to_numpy(float) - mu - a_g
    M_Ks = d["tmass_ks_m"].to_numpy(float) - mu - a_ks

    cstar_n = smp.corrected_excess_factor(
        bp_rp, d["phot_bp_rp_excess_factor"].to_numpy(float)) / np.maximum(
        smp.excess_factor_sigma(d["phot_g_mean_mag"].to_numpy(float)), 1e-9)

    is_var = (d["phot_variable_flag"].astype(str).str.upper() == "VARIABLE").to_numpy()
    ks_q = d["tmass_ph_qual"].astype(str)

    base = (np.isfinite(M_G) & np.isfinite(M_Ks) & np.isfinite(a0)
            & (M_Ks > 3.0) & (M_Ks < 8.0)
            & (bp_rp > 0.7) & (bp_rp < 3.6)
            & (d["dist_pc"].to_numpy() > 10) & (d["dist_pc"].to_numpy() < 500)
            & (a_g < 0.5)
            # keep only pristine cross-matches, so the beam plays no part
            & (d["tmass_xm_nnb"].fillna(1).to_numpy() <= 1)
            & (d["tmass_xm_nmates"].fillna(1).to_numpy() <= 1)
            & (d["tmass_xm_dist"].fillna(0).to_numpy() < 0.3))

    clean = base & ~is_var & (np.abs(cstar_n) < 3.0) & (ks_q.str[2] == "A").to_numpy()

    x, y = M_Ks[clean], M_G[clean]
    best = None
    for deg in (3, 4, 5):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, sigma, deg = best
    resid = M_G - np.polyval(coef, M_Ks)
    thr = args.k * sigma
    log.info("fiducial degree %d sigma %.4f; threshold %.3f mag",
             deg, sigma, thr)

    # ---- the amplitude proxy ---------------------------------------------
    fove = d["phot_g_mean_flux_over_error"].to_numpy(float)
    n_bp = d["phot_bp_n_obs"].fillna(0).to_numpy(float)
    n_g = np.where(n_bp > 0, 9.0 * n_bp, np.nan)
    amp = np.sqrt(n_g) / np.maximum(fove, 1e-9)      # fractional epoch scatter
    amp_mag = 2.5 / np.log(10) * amp                  # to magnitudes

    ok = base & np.isfinite(resid) & np.isfinite(amp_mag) & (amp_mag > 0)
    log.info("usable with an amplitude proxy: %d", int(ok.sum()))
    log.info("amplitude proxy: variables median %.4f mag, "
             "non-variables median %.4f mag",
             float(np.median(amp_mag[ok & is_var])),
             float(np.median(amp_mag[ok & ~is_var])))

    # ---- does the deficit scale with amplitude? --------------------------
    log.info("")
    log.info("=== deficit vs variability amplitude ===")
    edges = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 10.0]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        for label, sel in (("variable", is_var), ("non-variable", ~is_var)):
            m = ok & sel & (amp_mag >= lo) & (amp_mag < hi)
            if m.sum() < 200:
                continue
            r = resid[m]
            med = float(np.median(r))
            nd = int(np.count_nonzero(r > med + thr))
            nb = int(np.count_nonzero(r < med - thr))
            rows.append({"amp_lo": lo, "amp_hi": hi, "group": label,
                         "n": int(m.sum()), "median_resid": med,
                         "dim": nd, "bright": nb,
                         "ratio": (nd / nb) if nb else None})
            log.info("  amp %.3f-%-6.3f %-13s n=%7d  med=%+.4f  "
                     "dim=%6d bright=%5d  ratio=%s",
                     lo, hi, label, int(m.sum()), med, nd, nb,
                     f"{nd/nb:7.2f}" if nb else "    inf")

    # ---- the decisive correlation ----------------------------------------
    mv = ok & is_var
    corr = float(pd.Series(resid[mv]).corr(pd.Series(amp_mag[mv]),
                                           method="spearman"))
    log.info("")
    log.info("Spearman corr(residual, amplitude) among variables = %+.3f", corr)

    # Compare the LOWEST-amplitude variables to non-variables: if the deficit
    # is the mean-vs-epoch bias it must vanish where there is nothing to
    # average over.
    low_amp = ok & (amp_mag < 0.01)
    rv = resid[low_amp & is_var]
    rn = resid[low_amp & ~is_var]
    if len(rv) > 200 and len(rn) > 200:
        mv_, mn_ = float(np.median(rv)), float(np.median(rn))
        ndv = int(np.count_nonzero(rv > mv_ + thr))
        nbv = int(np.count_nonzero(rv < mv_ - thr))
        ndn = int(np.count_nonzero(rn > mn_ + thr))
        nbn = int(np.count_nonzero(rn < mn_ - thr))
        rat_v = (ndv / nbv) if nbv else np.inf
        rat_n = (ndn / nbn) if nbn else np.inf
        log.info("low-amplitude (<0.01 mag) variables : n=%d ratio=%s",
                 len(rv), f"{rat_v:.2f}" if np.isfinite(rat_v) else "inf")
        log.info("low-amplitude non-variables         : n=%d ratio=%.2f",
                 len(rn), rat_n)
    else:
        rat_v = rat_n = float("nan")

    # ---- verdict ----------------------------------------------------------
    if corr > 0.15 and np.isfinite(rat_v) and np.isfinite(rat_n) \
            and rat_v < 3.0 * rat_n:
        verdict = (
            f"MEAN-VERSUS-EPOCH BIAS. The deficit scales with variability "
            f"amplitude (Spearman {corr:+.3f}), and at amplitudes below "
            f"0.01 mag the variables are no more one-sidedly dim than "
            f"constant stars ({rat_v:.2f} vs {rat_n:.2f}). The excess is the "
            f"comparison of a 34-month Gaia mean against a single 2MASS epoch "
            f"from ~1999: any light curve that dips has a mean fainter than a "
            f"random sample of it. No absorber is required, and the effect "
            f"disappears exactly where there is nothing to average over.")
    elif corr > 0.15:
        verdict = (
            f"PARTLY the mean-versus-epoch bias: the deficit scales with "
            f"amplitude (Spearman {corr:+.3f}), but low-amplitude variables "
            f"remain more one-sidedly dim than constant stars "
            f"({rat_v:.2f} vs {rat_n:.2f}). Amplitude does not account for "
            f"all of it and the residual needs a separate explanation.")
    else:
        verdict = (
            f"NOT the mean-versus-epoch bias. The deficit does not scale with "
            f"variability amplitude (Spearman {corr:+.3f}), so averaging over "
            f"dips is not what makes these stars faint. The one-sided excess "
            f"in the variable class survives both the beam test and the "
            f"amplitude test and has no explanation yet.")

    print(f"\n{'='*72}")
    print("WHY ARE VARIABLE STARS ONE-SIDEDLY DIM?")
    print(f"{'='*72}")
    print(f"  Spearman corr(residual, amplitude) among variables : {corr:+.3f}")
    print(f"  low-amplitude (<0.01 mag) variables  dim:bright    : "
          f"{rat_v:.2f}" if np.isfinite(rat_v) else "  low-amp variables: inf")
    print(f"  low-amplitude non-variables          dim:bright    : "
          f"{rat_n:.2f}" if np.isfinite(rat_n) else "  low-amp constants: inf")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"variable_deficit_test_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "threshold_mag": float(thr),
        "fiducial_sigma": float(sigma),
        "spearman_resid_amplitude": corr,
        "amplitude_bins": rows,
        "low_amplitude_variable_ratio": (
            float(rat_v) if np.isfinite(rat_v) else None),
        "low_amplitude_constant_ratio": (
            float(rat_n) if np.isfinite(rat_n) else None),
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
