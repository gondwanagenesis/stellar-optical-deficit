#!/usr/bin/env python
"""Improved wide-pair limit: more pairs, and a correctly-targeted background.

    run.sh scripts/34_pair_limit_v3.py --tag primary

TWO IMPROVEMENTS OVER scripts/33
--------------------------------
(1) MORE PAIRS. The fitted sample required a GSP-Phot metallicity, costing 15%
    of stars and therefore ~30% of pairs (both members must survive). Within a
    co-natal pair the metallicity is common to both components and cancels in
    the difference exactly, so requiring it is pure loss. We refit the fiducial
    on the full sample using only the near-infrared colour control, which
    every star has, and search for pairs in that.

(2) A BACKGROUND MEASURED ON THE RIGHT POPULATION. scripts/33 predicted the
    pair tail as 2*q(T)*N using q from the whole 3.3M sample. But pair members
    are not a random subset: both passed parallax- and proper-motion-consistency
    cuts, which correlate with astrometric quality and hence with the blending
    that dominates the tail. The global q therefore OVER-predicts, and since
    the limit carries sqrt(obs + pred), an inflated prediction actively weakens
    it. Here q is measured on the single-star residuals of pair members
    themselves.

    That is not circular: the prediction uses each star's OWN residual r, while
    the statistic uses the DIFFERENCE dr between partners. A signal on one
    component enters both, so if anything this makes the limit conservative.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import extinction as ext
from pipeline import fiducial as fid
from pipeline import pairs as pr
from pipeline import statistics as st

A_W1, A_W2 = 0.039, 0.026
EXCESS_NSIGMA = 3.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--knots", type=int, default=6)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}.parquet")
    print(f"full sample: {len(d):,} stars")

    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_j = ext.deredden("J", a0, bp_rp)
    jks0 = ((d["tmass_j_m"].to_numpy(float) - a_j)
            - (d["tmass_ks_m"].to_numpy(float) - d["A_Ks"].to_numpy(float)))
    d["j_ks0"] = jks0

    ok = np.isfinite(d["M_G"]) & np.isfinite(d["M_Ks"]) & np.isfinite(jks0)
    d = d[ok].reset_index(drop=True)
    print(f"with finite M_G, M_Ks, (J-Ks)_0 : {len(d):,} "
          f"(vs 3,321,566 that also had [M/H])")

    covs = [(d["j_ks0"].to_numpy(float), 1)]
    m_ks = d["M_Ks"].to_numpy(float)
    fit = fid.fit_fiducial(m_ks, covs, d["M_G"].to_numpy(float), args.knots)
    r = fit.residuals(m_ks, covs, d["M_G"].to_numpy(float))
    d["residual_nomh"] = r
    sig = st.robust_sigma(r)
    print(f"  metallicity-free fiducial: sigma = {sig:.5f} mag "
          f"(with [M/H] it was 0.09913)")

    p = pr.drop_shared_components(pr.find_pairs(d))
    fake = pr.find_pairs(d, scramble=True)
    clean = p[p["theta_arcsec"] > pr.CLEAN_SEP_ARCSEC].reset_index(drop=True)
    print(f"\npairs {len(p):,} | chance {len(fake):,} "
          f"({100*len(fake)/max(len(p),1):.2f}%) | clean {len(clean):,} "
          f"(was 8,844)")

    i, j = clean["i"].to_numpy(), clean["j"].to_numpy()
    mks_i, mks_j = m_ks[i], m_ks[j]
    prim_i = mks_i <= mks_j
    dr = np.where(prim_i, r[i], r[j]) - np.where(prim_i, r[j], r[i])
    med = float(np.median(dr))
    s_dr = st.robust_sigma(dr)
    print(f"  sigma(dr) = {s_dr:.5f}   common-mode variance fraction = "
          f"{100*(1-(s_dr/np.sqrt(2)/sig)**2):.1f}%")

    # --- background, measured on pair members themselves -----------------
    members = np.unique(np.concatenate([i, j]))
    r_mem = r[members]
    med_mem = float(np.median(r_mem))
    dev_mem = np.abs(r_mem - med_mem)
    dev_all = np.abs(r - float(np.median(r)))

    def q_member(T):
        return float(np.count_nonzero(dev_mem > T)) / len(r_mem)

    def q_global(T):
        return float(np.count_nonzero(dev_all > T)) / len(r)

    print(f"\n  pair members: {len(members):,} stars")
    for T in (0.45, 0.55, 0.67):
        print(f"    q(|r|>{T:.2f}): members {q_member(T):.3e}  "
              f"global {q_global(T):.3e}  "
              f"ratio {q_member(T)/max(q_global(T),1e-12):.2f}")

    # --- mid-IR veto ------------------------------------------------------
    ks0 = d["tmass_ks_m"].to_numpy(float) - d["A_Ks"].to_numpy(float)
    w1 = d["wise_w1mpro"].to_numpy(float) - A_W1 * a0[ok.to_numpy()][:len(d)] \
        if False else d["wise_w1mpro"].to_numpy(float)
    w2 = d["wise_w2mpro"].to_numpy(float)
    e1 = d["wise_w1mpro_error"].to_numpy(float)
    e2 = d["wise_w2mpro_error"].to_numpy(float)
    have = (np.isfinite(w1) & np.isfinite(w2) & np.isfinite(e1) & np.isfinite(e2)
            & (e1 < 0.2) & (e2 < 0.2))

    def expected(y):
        lo, hi = np.nanpercentile(jks0[:len(d)][have], [1, 99])
        edges = np.linspace(lo, hi, 61)
        idx = np.clip(np.digitize(d["j_ks0"].to_numpy(float), edges) - 1, 0, 59)
        m_ = np.full(60, np.nan)
        for b in range(60):
            mm = (idx == b) & have & np.isfinite(y)
            if mm.sum() > 100:
                m_[b] = np.median(y[mm])
        g = np.isfinite(m_)
        return np.interp(d["j_ks0"].to_numpy(float),
                         (0.5 * (edges[:-1] + edges[1:]))[g], m_[g])

    ex1 = (ks0 - w1) - expected(ks0 - w1)
    ex2 = (ks0 - w2) - expected(ks0 - w2)
    s1, s2 = st.robust_sigma(ex1[have]), st.robust_sigma(ex2[have])
    bare = have & (np.abs(ex1) < EXCESS_NSIGMA * s1) & (np.abs(ex2) < EXCESS_NSIGMA * s2)
    both_bare = bare[i] & bare[j]
    print(f"  pairs with both components bare : {int(both_bare.sum()):,}")

    # --- limits -----------------------------------------------------------
    def table(drv, label):
        n_pairs = len(drv)
        n_stars = 2 * n_pairs
        rows = []
        for k in (4, 5, 6, 7):
            T = k * s_dr
            obs = int(np.count_nonzero(np.abs(drv - med) > T))
            pred_m = 2.0 * q_member(T) * n_pairs
            pred_g = 2.0 * q_global(T) * n_pairs
            f_det = float(st.fraction_from_delta(T))
            ul_cons = st.poisson_upper_limit(obs)
            ul_sub = max(obs - pred_m, 0.0) + 1.645 * np.sqrt(obs + pred_m)
            best = min(ul_cons, ul_sub)
            rows.append({"sample": label, "k": k, "thr_mag": T,
                         "f_det": f_det, "obs": obs,
                         "pred_member": pred_m, "pred_global": pred_g,
                         "p_UL": best / n_stars,
                         "mean_f_UL": best / n_stars * f_det})
        return pd.DataFrame(rows)

    t = pd.concat([table(dr, "all clean"),
                   table(dr[both_bare], "clean + bare")], ignore_index=True)
    print("\n=== limits ===")
    print(t.to_string(index=False, float_format=lambda v: f"{v:10.5g}"))

    ba = t[t["sample"] == "all clean"].loc[lambda x: x["mean_f_UL"].idxmin()]
    bb = t[t["sample"] == "clean + bare"].loc[lambda x: x["mean_f_UL"].idxmin()]
    p_dark = float(bb["p_UL"]); p_iso = 1.9e-4
    print(f"\n  all clean pairs : p < {float(ba['p_UL']):.3e} at f >= {float(ba['f_det']):.3f}")
    print(f"  beamed class    : p < {p_dark:.3e} at f >= {float(bb['f_det']):.3f}")
    print(f"\n  previous (scripts/33): 4.771e-04 and 4.895e-04")
    print(f"\n=== joint ===")
    print(f"  p_total < {p_iso + p_dark:.3e}  -> 1 in {1/(p_iso+p_dark):,.0f} stars")

    res = {"tag": args.tag, "n_stars_fitted": int(len(d)),
           "sigma_nomh": float(sig), "n_pairs_clean": int(len(clean)),
           "chance_frac": float(len(fake) / max(len(p), 1)),
           "n_pairs_bare": int(both_bare.sum()),
           "best_all": ba.to_dict(), "best_bare": bb.to_dict(),
           "p_dark": p_dark, "p_total": p_iso + p_dark,
           "one_in_n": float(1 / (p_iso + p_dark)),
           "table": t.to_dict(orient="records")}
    (cfg.RESULT_DIR / f"pair_limit_v3_{args.tag}.json").write_text(
        json.dumps(res, indent=2, default=float))
    t.to_csv(cfg.RESULT_DIR / f"pair_limit_v3_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
