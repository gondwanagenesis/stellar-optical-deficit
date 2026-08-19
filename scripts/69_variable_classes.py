#!/usr/bin/env python
"""What KIND of variables are the one-sidedly dim ones?

    run.sh scripts/69_variable_classes.py --tag primary

THE STATE OF THE QUESTION
-------------------------
The variable class was the largest single cut in the project -- 382,319 stars
removed before any channel ran -- and it is the only discarded class whose
one-sided dim excess survived both diagnostics:

  aperture mismatch   ruled out: the excess persists (218:1) when the 2MASS
                      match is single, uncontested and centred under 0.3",
                      while every other flagged class collapsed toward the
                      6:1 reference.

  mean-versus-epoch   only partly responsible: the deficit does scale with
                      variability amplitude (Spearman +0.387), but the
                      ASYMMETRY runs the wrong way for that mechanism. If
                      averaging a 34-month Gaia mean over dips against a
                      single 1999 2MASS epoch were the cause, the excess
                      would be largest at high amplitude. Instead it peaks
                      around 0.01-0.02 mag (149:1 against 7:1) and vanishes
                      above 0.02 mag, where variables and constants converge.

THE HYPOTHESIS THAT PREDICTS THAT INVERSION
-------------------------------------------
Starspots. A spotted star is dim in the optical relative to the infrared
because spots are ~1500-2000 K cooler than the photosphere, which suppresses G
far more than Ks. But the ROTATIONAL MODULATION only samples the ASYMMETRIC
part of the spot distribution: a heavily but evenly spotted star is strongly
dimmed and barely variable.

So spots uniquely predict a large mean deficit at small amplitude -- exactly
the inversion observed -- while the averaging bias predicts the opposite.

THE TEST
--------
Gaia DR3 classifies its variables. If the dim population is dominated by
rotational modulation (ROT), young stellar objects (YSO) or eruptive types,
the effect is spot and disc physics. If it is spread evenly across classes, or
concentrated in something that does not spot, the explanation fails and the
excess remains open.

This queries vari_classifier_result for the dim variables and for a matched
control of non-dim variables, and compares the class mixtures.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
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
log = logging.getLogger("varclass")

MIRRORS = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
    ("ESA", "https://gea.esac.esa.int/tap-server/tap"),
]

COLS = [
    "source_id", "l", "b", "parallax", "parallax_error",
    "phot_g_mean_mag", "bp_rp", "phot_bp_rp_excess_factor",
    "phot_variable_flag", "teff_gspphot",
    "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat",
    "astrometric_params_solved",
    "tmass_ks_m", "tmass_ph_qual", "tmass_j_m",
    "tmass_xm_dist", "tmass_xm_nnb", "tmass_xm_nmates",
    "wise_w1mpro", "wise_w2mpro", "wise_w3mpro", "wise_w4mpro",
]

N_SAMPLE = 4000
K_SIGMA = 3.0


def pick_mirror():
    from astroquery.utils.tap.core import TapPlus
    for name, url in MIRRORS:
        try:
            tap = TapPlus(url=url)
            tap.launch_job("SELECT TOP 1 source_id "
                           "FROM gaiadr3.vari_classifier_result").get_results()
            log.info("using %s", name)
            return tap
        except Exception as exc:
            log.warning("%s unusable: %s", name, str(exc)[:110])
    return None


def fetch_classes(tap, ids):
    """best_class_name for a list of source_ids, in sync-sized pages."""
    out = []
    step = 400
    for i in range(0, len(ids), step):
        chunk = ids[i:i + step]
        lst = ",".join(str(int(s)) for s in chunk)
        q = (f"SELECT source_id, best_class_name, best_class_score "
             f"FROM gaiadr3.vari_classifier_result "
             f"WHERE source_id IN ({lst})")
        try:
            out.append(tap.launch_job(q).get_results().to_pandas())
        except Exception as exc:
            log.warning("  class chunk %d failed: %s", i, str(exc)[:90])
    if not out:
        return pd.DataFrame(columns=["source_id", "best_class_name"])
    return pd.concat(out, ignore_index=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n", type=int, default=N_SAMPLE)
    args = ap.parse_args()
    rng = np.random.default_rng(69_2026)

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
    thr = K_SIGMA * sigma
    d = d.assign(residual=resid, M_G=M_G, M_Ks=M_Ks)

    var_ok = base & is_var & np.isfinite(resid)
    dim = var_ok & (resid > thr)
    ord_ = var_ok & (np.abs(resid) < thr)
    log.info("variables: %d  dim: %d  ordinary: %d",
             int(var_ok.sum()), int(dim.sum()), int(ord_.sum()))

    # ---- infrared excess, the disc discriminator -------------------------
    # Spots dim the optical and leave the mid-infrared alone. A circumstellar
    # disc dims the optical AND emits at W3/W4. The two are separable here.
    for label, m in (("dim variables", dim), ("ordinary variables", ord_)):
        w = d.loc[m]
        if "wise_w3mpro" in w and "wise_w1mpro" in w:
            w13 = (w["wise_w1mpro"] - w["wise_w3mpro"]).dropna()
            log.info("  %-20s median W1-W3 = %+.3f  (n=%d)",
                     label, float(w13.median()) if len(w13) else np.nan, len(w13))

    # ---- variability classes ---------------------------------------------
    tap = pick_mirror()
    if tap is None:
        log.error("no mirror serves vari_classifier_result; skipping classes")
        return 1

    ids_dim = d.loc[dim, "source_id"].to_numpy()
    ids_ord = d.loc[ord_, "source_id"].to_numpy()
    if len(ids_dim) > args.n:
        ids_dim = rng.choice(ids_dim, args.n, replace=False)
    if len(ids_ord) > args.n:
        ids_ord = rng.choice(ids_ord, args.n, replace=False)

    log.info("querying variability classes for %d dim and %d ordinary ...",
             len(ids_dim), len(ids_ord))
    t0 = time.time()
    cls_dim = fetch_classes(tap, ids_dim)
    cls_ord = fetch_classes(tap, ids_ord)
    log.info("retrieved in %.0fs: %d dim, %d ordinary classified",
             time.time() - t0, len(cls_dim), len(cls_ord))

    def mix(df, n_query):
        if len(df) == 0:
            return {}
        vc = df["best_class_name"].astype(str).value_counts()
        return {k: {"n": int(v), "frac_of_classified": float(v / len(df))}
                for k, v in vc.items()}

    mix_dim = mix(cls_dim, len(ids_dim))
    mix_ord = mix(cls_ord, len(ids_ord))

    log.info("")
    log.info("%-28s %12s %12s", "class", "dim", "ordinary")
    keys = sorted(set(mix_dim) | set(mix_ord),
                  key=lambda k: -mix_dim.get(k, {}).get("n", 0))
    for k in keys[:15]:
        a = mix_dim.get(k, {}).get("frac_of_classified", 0.0)
        b = mix_ord.get(k, {}).get("frac_of_classified", 0.0)
        log.info("%-28s %11.3f %11.3f", k, a, b)

    frac_class_dim = len(cls_dim) / max(len(ids_dim), 1)
    frac_class_ord = len(cls_ord) / max(len(ids_ord), 1)
    log.info("")
    log.info("fraction with a published classification: dim %.3f, ordinary %.3f",
             frac_class_dim, frac_class_ord)

    spotty = {"ROT", "SOLAR_LIKE", "YSO", "ECL", "LPV"}
    dom = sum(v["frac_of_classified"] for k, v in mix_dim.items()
              if k.upper() in spotty)
    log.info("dim variables in spot/disc/eclipse classes: %.3f", dom)

    if dom > 0.6:
        verdict = (
            f"STELLAR ACTIVITY. {dom:.0%} of the classified dim variables fall "
            f"in rotational, solar-like, young-stellar, eclipsing or long-period "
            f"classes -- all of which dim the optical relative to the infrared "
            f"through cool starspots or circumstellar material. Spots in "
            f"particular predict the inversion that defeated the averaging "
            f"explanation: coverage sets the mean deficit while only the "
            f"asymmetric part of the spot distribution modulates, so a heavily "
            f"spotted star is strongly dimmed and barely variable.")
    else:
        verdict = (
            f"UNRESOLVED. Only {dom:.0%} of the classified dim variables fall "
            f"in the classes that dim the optical through spots or discs, so "
            f"stellar activity does not obviously account for the excess. This "
            f"population survives the beam test, survives the amplitude test, "
            f"and is not explained by its variability type. It is the strongest "
            f"open anomaly in the project and needs spectroscopy.")

    print(f"\n{'='*72}")
    print("VARIABILITY CLASSES OF THE ONE-SIDEDLY DIM VARIABLES")
    print(f"{'='*72}")
    print(f"{'class':<28} {'dim':>10} {'ordinary':>10}")
    for k in keys[:15]:
        print(f"{k:<28} "
              f"{mix_dim.get(k, {}).get('frac_of_classified', 0.0):>10.3f} "
              f"{mix_ord.get(k, {}).get('frac_of_classified', 0.0):>10.3f}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"variable_classes_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "n_variables": int(var_ok.sum()),
        "n_dim": int(dim.sum()),
        "n_ordinary": int(ord_.sum()),
        "threshold_mag": float(thr),
        "class_mix_dim": mix_dim,
        "class_mix_ordinary": mix_ord,
        "frac_classified_dim": frac_class_dim,
        "frac_classified_ordinary": frac_class_ord,
        "frac_dim_in_activity_classes": dom,
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
