#!/usr/bin/env python
"""Search N: mass without light, in the population no channel ever saw.

    run.sh scripts/66_searchN_mass_without_light.py --tag primary

THE POPULATION
--------------
Every one of the nineteen channels above ran on a sample whose ADQL carried a
server-side RUWE cut; the local partitions top out at ruwe = 1.387. So the
astrometrically-disturbed stars were never downloaded, never examined, and are
excluded from every limit this project has published.

For this particular search that is the worst possible thing to have discarded.
RUWE measures how badly a single-star astrometric model fits, and it rises
precisely when an unseen mass pulls the star around. The f -> 1 case that every
photometric channel is structurally blind to -- a companion carrying full mass
and emitting nothing -- would announce itself in RUWE and essentially nowhere
else.

THE DISCRIMINANT
----------------
High RUWE by itself means "unresolved binary" and is unremarkable: roughly a
third of nearby stars have it. The question is not whether these stars are
disturbed, it is whether the disturbance is accompanied by light.

An ordinary companion CONTRIBUTES FLUX. It inflates the BP/RP excess factor,
it shows up as a second 2MASS entry or a displaced cross-match, and it lifts
the system above the main sequence. So for real binaries, astrometric
disturbance and photometric blend evidence rise together.

An enshrouded companion carries the same mass and emits nothing. It breaks
that relation: full RUWE, zero photometric trace.

So the search is for stars that are astrometrically disturbed and
photometrically pristine, and then significantly dim on top.

WHY THIS RESCUES THE MIRROR CONTROL
-----------------------------------
Search M could not interpret its own one-sidedness, because this project's
dominant contaminant -- 2MASS aperture mismatch, where a 4 arcsec beam catches
a neighbour that Gaia resolves out -- makes Ks brighter and never fainter, so
the inferred deficit is one-signed by construction.

Requiring pristine blend indicators suppresses exactly that contaminant. In a
subsample with a single, well-centred 2MASS match and a normal BP/RP excess
factor, there is no unmodelled flux to inflate Ks, so photometric scatter goes
both ways again and the bright tail becomes a usable false-positive rate.
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
log = logging.getLogger("searchN")

K_SIGMA = 3.0


def absolute_mags(d: pd.DataFrame):
    """Distance, extinction and absolute magnitudes, as the main pipeline does."""
    d = smp.add_astrometry(d)
    a0 = ext.query_a0("edenhofer23", d["l"].to_numpy(float),
                      d["b"].to_numpy(float), d["dist_pc"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = d["dist_mod"].to_numpy(float)
    d = d.assign(
        A_0=a0,
        A_G=a_g,
        M_G=d["phot_g_mean_mag"].to_numpy(float) - mu - a_g,
        M_Ks=d["tmass_ks_m"].to_numpy(float) - mu - a_ks,
    )
    return d


def fit_reference_fiducial(tag: str):
    """Refit M_G(M_Ks) on the CLEAN retained sample, so the comparison is fair."""
    r = pd.read_parquet(cfg.DERIVED_DIR / f"{tag}_resid.parquet",
                        columns=["M_G", "M_Ks", "ruwe", "cstar_nsigma"])
    ok = (r["M_G"].notna() & r["M_Ks"].notna()
          & (r["ruwe"] < 1.4) & (r["cstar_nsigma"].abs() < 3))
    x = r.loc[ok, "M_Ks"].to_numpy(float)
    y = r.loc[ok, "M_G"].to_numpy(float)
    best = None
    for deg in (3, 4, 5):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, sigma, deg = best
    log.info("reference fiducial from %d clean stars: degree %d, sigma %.4f",
             int(ok.sum()), deg, sigma)
    return coef, sigma


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=K_SIGMA)
    args = ap.parse_args()

    src = cfg.RAW_DIR / "high_ruwe_500pc.parquet"
    if not src.exists():
        log.error("missing %s -- run scripts/65_pull_high_ruwe.py first", src)
        return 1

    d = pd.read_parquet(src)
    log.info("high-RUWE sources pulled: %d", len(d))

    d = d[d["tmass_ks_m"].notna()].reset_index(drop=True)
    log.info("with 2MASS Ks: %d", len(d))

    d = absolute_mags(d)

    coef, sigma = fit_reference_fiducial(args.tag)
    resid = d["M_G"].to_numpy(float) - np.polyval(coef, d["M_Ks"].to_numpy(float))
    d = d.assign(residual=resid)

    # ---- same main-sequence box as the rest of the project ---------------
    bp_rp = d["bp_rp"].to_numpy(float)
    base = (np.isfinite(resid) & np.isfinite(d["A_0"].to_numpy(float))
            & (d["M_Ks"] > 3.0) & (d["M_Ks"] < 8.0)
            & (bp_rp > 0.7) & (bp_rp < 3.6)
            & (d["dist_pc"] > 10) & (d["dist_pc"] < 500)
            & (d["A_G"] < 0.5))
    d = d[base].reset_index(drop=True)
    resid = d["residual"].to_numpy(float)
    log.info("in the lower-main-sequence box: %d", len(d))

    # ---- photometric cleanliness -----------------------------------------
    cstar = smp.corrected_excess_factor(
        d["bp_rp"].to_numpy(float),
        d["phot_bp_rp_excess_factor"].to_numpy(float))
    csig = smp.excess_factor_sigma(d["phot_g_mean_mag"].to_numpy(float))
    cstar_n = cstar / np.maximum(csig, 1e-9)
    d = d.assign(cstar_nsigma=cstar_n)

    nnb = d.get("tmass_xm_nnb", pd.Series(1, index=d.index)).fillna(1).to_numpy()
    xmd = d.get("tmass_xm_dist", pd.Series(0.0, index=d.index)).fillna(0.0).to_numpy()
    ipd = d.get("ipd_frac_multi_peak", pd.Series(0, index=d.index)).fillna(0).to_numpy()
    dup = d.get("duplicated_source", pd.Series(False, index=d.index)).fillna(False).to_numpy(bool)

    pristine = ((np.abs(cstar_n) < 3.0) & (nnb <= 1) & (xmd < 0.5)
                & (ipd <= 2) & ~dup)
    log.info("photometrically pristine despite high RUWE: %d (%.1f%%)",
             int(pristine.sum()), 100 * pristine.mean())

    thr = args.k * sigma
    f_det = float(st.fraction_from_delta(thr))
    log.info("threshold %.0f sigma = %.3f mag (f >= %.3f)",
             args.k, thr, f_det)

    # ---- does astrometric disturbance track photometric blending? --------
    # For real unresolved binaries it must. Breaking that relation is the
    # signature of mass carrying no light.
    log.info("")
    log.info("=== does RUWE track photometric blend evidence? ===")
    ruwe = d["ruwe"].to_numpy(float)
    rows = []
    for lo, hi in [(1.4, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 10.0), (10.0, 1e9)]:
        m = (ruwe >= lo) & (ruwe < hi)
        if m.sum() < 50:
            continue
        rows.append({
            "ruwe_lo": lo, "ruwe_hi": hi, "n": int(m.sum()),
            "frac_pristine": float(pristine[m].mean()),
            "median_cstar": float(np.median(cstar_n[m])),
            "median_resid": float(np.median(resid[m])),
        })
        log.info("  RUWE %5.1f-%-6.1f n=%6d  pristine=%.3f  "
                 "med C*=%+.2f  med resid=%+.3f",
                 lo, hi, int(m.sum()), float(pristine[m].mean()),
                 float(np.median(cstar_n[m])), float(np.median(resid[m])))

    # ---- the search, inside the pristine subsample -----------------------
    log.info("")
    log.info("=== deficit within the pristine high-RUWE subsample ===")
    r_p = resid[pristine]
    med = float(np.median(r_p))
    n_dim = int(np.count_nonzero(r_p > med + thr))
    n_bright = int(np.count_nonzero(r_p < med - thr))
    ratio = n_dim / n_bright if n_bright else np.inf
    log.info("  n=%d  median resid=%+.4f  robust sigma=%.4f",
             len(r_p), med, float(st.robust_sigma(r_p)))
    log.info("  dim=%d  bright=%d  ratio=%s", n_dim, n_bright,
             f"{ratio:.2f}" if np.isfinite(ratio) else "inf")

    # blended comparison group: same RUWE range, photometry NOT pristine
    r_b = resid[~pristine]
    med_b = float(np.median(r_b))
    nd_b = int(np.count_nonzero(r_b > med_b + thr))
    nb_b = int(np.count_nonzero(r_b < med_b - thr))
    ratio_b = nd_b / nb_b if nb_b else np.inf
    log.info("  blended comparison: n=%d median=%+.4f dim=%d bright=%d "
             "ratio=%s", len(r_b), med_b, nd_b, nb_b,
             f"{ratio_b:.2f}" if np.isfinite(ratio_b) else "inf")

    # ---- reference: the clean low-RUWE sample ----------------------------
    ref = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                          columns=["residual", "ruwe", "cstar_nsigma"])
    rok = (ref["ruwe"] < 1.4) & (ref["cstar_nsigma"].abs() < 3) \
        & ref["residual"].notna()
    r_ref = ref.loc[rok, "residual"].to_numpy(float)
    med_r = float(np.median(r_ref))
    nd_r = int(np.count_nonzero(r_ref > med_r + thr))
    nb_r = int(np.count_nonzero(r_ref < med_r - thr))
    ratio_r = nd_r / nb_r if nb_r else np.inf
    log.info("  low-RUWE reference : n=%d dim=%d bright=%d ratio=%.2f",
             len(r_ref), nd_r, nb_r, ratio_r)

    from scipy import stats
    p_ref = ratio_r / (1.0 + ratio_r)
    p_val = stats.binomtest(n_dim, max(n_dim + n_bright, 1), p_ref,
                            alternative="greater").pvalue if (n_dim + n_bright) else 1.0
    log.info("  p(pristine high-RUWE more one-sided than reference) = %.3g",
             p_val)

    # ---- verdict ----------------------------------------------------------
    excess_over_ref = (ratio / ratio_r) if (np.isfinite(ratio) and ratio_r) else np.inf

    if n_dim == 0:
        verdict = (
            f"NULL. Among {int(pristine.sum()):,} astrometrically disturbed "
            f"but photometrically pristine stars, none is significantly dim. "
            f"The population no other channel examined contains no case of "
            f"mass without light at f >= {f_det:.2f}.")
    elif p_val > 1e-3:
        verdict = (
            f"NULL. The pristine high-RUWE subsample is no more one-sidedly "
            f"dim than the clean low-RUWE reference "
            f"({ratio:.2f} vs {ratio_r:.2f}, p = {p_val:.2g}). Astrometric "
            f"disturbance alone does not produce an optical deficit once "
            f"blended photometry is excluded.")
    else:
        verdict = (
            f"EXCESS: pristine high-RUWE stars are {excess_over_ref:.1f}x more "
            f"one-sidedly dim than the clean low-RUWE reference "
            f"({ratio:.2f} vs {ratio_r:.2f}, p = {p_val:.2g}), on "
            f"{n_dim} objects. These are astrometrically disturbed, "
            f"photometrically clean, and faint. That is the mass-without-light "
            f"signature. Required next: confirm RUWE is not itself driven by a "
            f"resolved neighbour below the 2MASS matching radius, and obtain "
            f"spectroscopy -- a real absorber leaves the photosphere untouched "
            f"while a faint companion does not.")

    print(f"\n{'='*74}")
    print("SEARCH N: MASS WITHOUT LIGHT (high-RUWE, never before examined)")
    print(f"{'='*74}")
    print(f"  high-RUWE sources in the MS box       : {len(d):,}")
    print(f"  photometrically pristine among them   : {int(pristine.sum()):,}")
    print(f"  threshold {args.k:.0f} sigma                    : "
          f"{thr:.3f} mag (f >= {f_det:.3f})")
    print()
    print(f"  pristine high-RUWE  dim/bright        : "
          f"{n_dim}/{n_bright} = "
          f"{ratio:.2f}" if np.isfinite(ratio) else "  pristine: inf")
    print(f"  blended  high-RUWE  dim/bright        : {nd_b}/{nb_b} = "
          f"{ratio_b:.2f}" if np.isfinite(ratio_b) else "  blended: inf")
    print(f"  clean low-RUWE reference              : {nd_r}/{nb_r} = "
          f"{ratio_r:.2f}")
    print(f"  p vs reference                        : {p_val:.3g}")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "n_high_ruwe_pulled": int(len(d)),
        "n_pristine": int(pristine.sum()),
        "k_sigma": args.k,
        "threshold_mag": float(thr),
        "f_detectable": f_det,
        "ruwe_vs_blend": rows,
        "pristine": {"n": len(r_p), "median": med, "dim": n_dim,
                     "bright": n_bright,
                     "ratio": float(ratio) if np.isfinite(ratio) else None},
        "blended": {"n": len(r_b), "median": med_b, "dim": nd_b,
                    "bright": nb_b,
                    "ratio": float(ratio_b) if np.isfinite(ratio_b) else None},
        "reference_low_ruwe": {"n": len(r_ref), "dim": nd_r, "bright": nb_r,
                               "ratio": float(ratio_r)},
        "p_vs_reference": float(p_val),
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchN_mass_without_light_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)

    if n_dim:
        cand = d[pristine & (resid > med + thr)].copy()
        csv = cfg.RESULT_DIR / f"searchN_candidates_{args.tag}.csv"
        cand.sort_values("residual", ascending=False).to_csv(csv, index=False)
        log.info("wrote %d candidates to %s", len(cand), csv)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
