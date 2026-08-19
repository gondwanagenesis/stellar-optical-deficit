#!/usr/bin/env python
"""Search M: audit the 20% of the sample we threw away before looking.

    run.sh scripts/64_searchM_discard_pile.py --tag primary

THE PROBLEM WITH EVERY RESULT ABOVE
-----------------------------------
Nineteen channels all ran on a sample that had already been cleaned. The
cutflow removed 995,789 of 4,879,956 sources -- 20.4% -- before any search
began, and every one of those cuts was chosen to remove things that look like
bad data.

But the target of this project does not look like good data. Each rejection
criterion is also a plausible signature of the thing being hunted:

  phot_variable_flag (382,319 removed)
      The largest single cut. Wright et al. (2016) list aperiodic variability
      of varied depths as the PRIMARY expected signature of a transiting
      megastructure -- "many, aperiodic transit signatures of varied shapes,
      varied depths, no wavelength dependence". We deleted that population
      first, as noise.

  RUWE > 1.4 (applied per channel)
      RUWE measures astrometric excess scatter, which means an unseen
      companion is pulling on the star. A completely enshrouded companion is
      an unseen companion. We removed exactly the population where one would
      have to live.

  |C*| > 3 (33,612 removed)
      BP/RP flux inconsistent with G. That is what a source whose energy
      distribution does not match a bare photosphere looks like.

  classprob_dsc_combmod_star < 0.5 (173,856 removed)
      Rejected as extragalactic. A star dimmed enough, or with a sufficiently
      odd colour, is misclassified by a classifier trained on normal stars.

  non_single_star != 0, duplicated_source, ipd_frac_multi_peak
      Partially resolved or blended sources.

THE TEST
--------
Recompute the optical deficit for every discarded class against the SAME
fiducial fitted on the retained sample, then ask a question that noise cannot
fake: is the excess ONE-SIDED?

An absorber can only dim. Photometric pathology, blending, misclassification
and variability all scatter a residual in BOTH directions. So for each discard
class we count significantly-dim and significantly-bright sources and compare
the two. A class that is symmetrically broadened is explained by its rejection
reason. A class that is dim-heavy is not, and has to be explained some other
way.

This does not assume the discards are clean. It assumes only that whatever
made them dirty is not systematically one-signed -- and where that assumption
fails (extinction, which only ever dims) the extinction is measured per source
and controlled explicitly.
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
log = logging.getLogger("searchM")

NEEDED = [
    "source_id", "l", "b", "parallax", "parallax_error", "parallax_over_error",
    "ruwe", "phot_g_mean_mag", "bp_rp", "phot_bp_rp_excess_factor",
    "phot_variable_flag", "non_single_star", "duplicated_source",
    "in_qso_candidates", "in_galaxy_candidates", "classprob_dsc_combmod_star",
    "ipd_frac_multi_peak", "astrometric_excess_noise_sig",
    "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat",
    "astrometric_params_solved",
    "tmass_ks_m", "tmass_ks_msigcom", "tmass_ph_qual",
    "mh_gspphot",
]

K_SIGMA = 3.0


def load_all(pattern="sample_d500_p*.parquet") -> pd.DataFrame:
    files = sorted(cfg.RAW_DIR.glob(pattern))
    log.info("loading %d raw partitions", len(files))
    frames = []
    for i, f in enumerate(files):
        import pyarrow.parquet as pq
        have = set(pq.read_schema(f).names)
        cols = [c for c in NEEDED if c in have]
        frames.append(pd.read_parquet(f, columns=cols))
        if (i + 1) % 50 == 0:
            log.info("  %d/%d", i + 1, len(files))
    d = pd.concat(frames, ignore_index=True)
    log.info("raw joined rows: %d", len(d))
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=K_SIGMA)
    args = ap.parse_args()

    d = load_all()

    # ---- label WHY each source would be discarded ------------------------
    var = d.get("phot_variable_flag")
    is_var = (var.astype(str).str.upper() == "VARIABLE") if var is not None \
        else pd.Series(False, index=d.index)

    nss = d.get("non_single_star", pd.Series(0, index=d.index)).fillna(0)
    dup = d.get("duplicated_source", pd.Series(False, index=d.index)).fillna(False)
    qso = d.get("in_qso_candidates", pd.Series(False, index=d.index)).fillna(False)
    gal = d.get("in_galaxy_candidates", pd.Series(False, index=d.index)).fillna(False)
    dsc = d.get("classprob_dsc_combmod_star", pd.Series(1.0, index=d.index)).fillna(1.0)
    ruwe = d.get("ruwe", pd.Series(1.0, index=d.index)).fillna(1.0)
    ipd = d.get("ipd_frac_multi_peak", pd.Series(0, index=d.index)).fillna(0)
    ks_q = d.get("tmass_ph_qual", pd.Series("AAA", index=d.index)).astype(str)

    # C* : use the pipeline's own Riello+2021 implementation. The locus is
    # piecewise in colour; approximating it with a single polynomial rejects
    # ~91% of the sample instead of 0.7%.
    bp_rp = d["bp_rp"].to_numpy(float)
    cstar = smp.corrected_excess_factor(
        bp_rp, d["phot_bp_rp_excess_factor"].to_numpy(float))
    sigma_cstar = smp.excess_factor_sigma(
        d["phot_g_mean_mag"].to_numpy(float))
    cstar_n = cstar / np.maximum(sigma_cstar, 1e-9)

    # ---- geometry and photometry for EVERYTHING --------------------------
    d = smp.add_astrometry(d)
    a0 = ext.query_a0("edenhofer23", d["l"].to_numpy(float),
                      d["b"].to_numpy(float), d["dist_pc"].to_numpy(float))
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = d["dist_mod"].to_numpy(float)

    M_G = d["phot_g_mean_mag"].to_numpy(float) - mu - a_g
    M_Ks = d["tmass_ks_m"].to_numpy(float) - mu - a_ks
    d = d.assign(A_0=a0, M_G=M_G, M_Ks=M_Ks, cstar_nsigma=cstar_n)

    # ---- the retained sample, reproducing the published cutflow ----------
    base = (np.isfinite(M_G) & np.isfinite(M_Ks) & np.isfinite(a0)
            & (M_Ks > 3.0) & (M_Ks < 8.0)
            & (bp_rp > 0.7) & (bp_rp < 3.6)
            & (d["dist_pc"] > 10) & (d["dist_pc"] < 500)
            & (a_g < 0.5))

    clean = (base
             & ~is_var.to_numpy()
             & (nss.to_numpy() == 0)
             & ~qso.to_numpy(bool) & ~gal.to_numpy(bool)
             & (dsc.to_numpy() > 0.5)
             & (np.abs(cstar_n) < 3.0)
             & (ks_q.str[2] == "A").to_numpy())

    log.info("reproduced retained sample: %d", int(clean.sum()))

    # ---- fiducial fitted on the RETAINED sample only ---------------------
    x = M_Ks[clean]
    y = M_G[clean]
    best = None
    for deg in (3, 4, 5):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, sigma, deg = best
    log.info("fiducial: degree %d, robust scatter %.4f mag", deg, sigma)

    resid = M_G - np.polyval(coef, M_Ks)
    d = d.assign(residual=resid)

    thr = args.k * sigma
    log.info("threshold at %.0f sigma = %.3f mag "
             "(intercepted fraction f >= %.3f)",
             args.k, thr, float(st.fraction_from_delta(thr)))

    # ---- the discard classes ---------------------------------------------
    usable = base & np.isfinite(resid)

    # Some rejection classes are NOT auditable from this data at all, because
    # they were excluded server-side in the ADQL that built the sample: the
    # local partitions top out at ruwe = 1.326 and contain zero duplicated or
    # multi-peak sources. Those populations were never downloaded, so no
    # channel in this project has ever seen them. They are reported as an
    # explicit hole rather than silently scored as empty.
    ruwe_arr = ruwe.to_numpy()
    unauditable = {}
    if float(np.nanmax(ruwe_arr)) < 1.4:
        unauditable["RUWE > 1.4"] = (
            f"excluded server-side; local max ruwe = {np.nanmax(ruwe_arr):.3f}")
    if not dup.to_numpy(bool).any():
        unauditable["duplicated_source"] = "excluded server-side; none present"
    if not (ipd.to_numpy() > 2).any():
        unauditable["IPD multi-peak > 2"] = "excluded server-side; none present"
    for k, v in unauditable.items():
        log.warning("NOT AUDITABLE  %-24s %s", k, v)

    classes = {
        "RETAINED (reference)": clean & usable,
        "variable (phot_variable_flag)": usable & is_var.to_numpy(),
        "non_single_star != 0": usable & (nss.to_numpy() != 0),
        "C* beyond 3 sigma": usable & (np.abs(cstar_n) >= 3.0),
        "classified extragalactic": usable & ((dsc.to_numpy() <= 0.5)
                                              | qso.to_numpy(bool)
                                              | gal.to_numpy(bool)),
        "2MASS Ks not 'A'": usable & (ks_q.str[2] != "A").to_numpy(),
    }
    for k in list(unauditable):
        classes.pop(k, None)

    rows = []
    for name, m in classes.items():
        n = int(m.sum())
        if n < 200:
            log.info("%-32s n=%6d  (too few)", name, n)
            continue
        r = resid[m]
        med = float(np.median(r))
        n_dim = int(np.count_nonzero(r > thr))
        n_bright = int(np.count_nonzero(r < -thr))
        ratio = n_dim / n_bright if n_bright else np.inf
        rows.append({
            "class": name, "n": n,
            "median_resid": med,
            "robust_sigma": float(st.robust_sigma(r)),
            "n_dim": n_dim, "n_bright": n_bright,
            "frac_dim": n_dim / n,
            "dim_bright_ratio": float(ratio) if np.isfinite(ratio) else None,
        })
        log.info("%-32s n=%7d  med=%+.3f  sig=%.3f  dim=%6d  bright=%6d  "
                 "ratio=%s", name, n, med, float(st.robust_sigma(r)),
                 n_dim, n_bright,
                 f"{ratio:6.2f}" if np.isfinite(ratio) else "   inf")

    tab = pd.DataFrame(rows)

    # ---- which classes are ONE-SIDED beyond the reference? ---------------
    ref = tab[tab["class"].str.startswith("RETAINED")].iloc[0]
    ref_ratio = ref["dim_bright_ratio"] or 1.0
    log.info("")
    log.info("reference dim:bright ratio in the retained sample = %.2f",
             ref_ratio)

    from scipy import stats
    flagged = []
    for _, r in tab.iterrows():
        if r["class"].startswith("RETAINED"):
            continue
        nd, nb = int(r["n_dim"]), int(r["n_bright"])
        if nd + nb < 20:
            continue
        # Is this class more one-sided than the retained sample?
        p_ref = ref_ratio / (1.0 + ref_ratio)
        p = stats.binomtest(nd, nd + nb, p_ref, alternative="greater").pvalue
        r = r.to_dict()
        r["p_vs_reference"] = float(p)
        flagged.append(r)
        log.info("  %-32s dim:bright %6.2f   p(more one-sided than "
                 "retained) = %.3g", r["class"],
                 r["dim_bright_ratio"] or np.inf, p)

    sig = [f for f in flagged if f["p_vs_reference"] < 1e-3]

    # ---- extinction control ----------------------------------------------
    # Extinction is the one contaminant that IS one-signed, so any flagged
    # class has to survive a low-extinction cut.
    log.info("")
    log.info("extinction control (A_0 < 0.05):")
    low = usable & (a0 < 0.05)
    for f in sig:
        m = classes[f["class"]] & low
        if m.sum() < 200:
            f["low_ext_ratio"] = None
            continue
        r = resid[m]
        nd = int(np.count_nonzero(r > thr))
        nb = int(np.count_nonzero(r < -thr))
        f["low_ext_n"] = int(m.sum())
        f["low_ext_dim"] = nd
        f["low_ext_bright"] = nb
        f["low_ext_ratio"] = float(nd / nb) if nb else None
        log.info("  %-32s n=%6d  dim=%5d bright=%5d  ratio=%s",
                 f["class"], int(m.sum()), nd, nb,
                 f"{nd/nb:.2f}" if nb else "inf")

    # ---- verdict ----------------------------------------------------------
    survivors = [f for f in sig
                 if f.get("low_ext_ratio") is not None
                 and f["low_ext_ratio"] > ref_ratio * 1.5]

    if not sig:
        verdict = (
            "No discarded class is more one-sidedly dim than the retained "
            "sample. The 995,789 sources removed by the cutflow do not hide a "
            "population of anomalously faint stars: whatever made them fail "
            "quality control scatters their photometry in both directions, "
            "which is what pathology does and what an absorber cannot do.")
    elif not survivors:
        verdict = (
            f"{len(sig)} discarded classes are one-sidedly dim, but none "
            f"survives the low-extinction control at more than 1.5x the "
            f"retained sample's own ratio. The asymmetry tracks dust, which "
            f"is the one contaminant that only ever dims.")
    else:
        names = ", ".join(f["class"] for f in survivors)
        verdict = (
            f"ONE-SIDED DIM EXCESS survives in: {names}. These classes were "
            f"discarded before any channel ran, are dimmer than the retained "
            f"population in a way their rejection reason does not explain, "
            f"and the asymmetry persists at low extinction. This is the "
            f"population every result in this project excluded by "
            f"construction and it needs to be worked through individually.")

    print(f"\n{'='*78}")
    print("SEARCH M: THE DISCARD PILE")
    print(f"{'='*78}")
    print(tab.to_string(index=False, float_format=lambda v: f"{v:9.3f}"))
    print(f"\nthreshold {args.k:.0f} sigma = {thr:.3f} mag "
          f"(f >= {float(st.fraction_from_delta(thr)):.3f})")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "n_raw": int(len(d)),
        "n_retained": int(clean.sum()),
        "n_discarded": int(len(d) - clean.sum()),
        "fiducial_degree": int(deg),
        "fiducial_sigma": float(sigma),
        "k_sigma": args.k,
        "threshold_mag": float(thr),
        "reference_dim_bright_ratio": float(ref_ratio),
        "classes": rows,
        "flagged": flagged,
        "survivors_after_extinction_control": survivors,
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchM_discard_pile_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    log.info("wrote %s", out)

    tab.to_csv(cfg.RESULT_DIR / f"searchM_classes_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
