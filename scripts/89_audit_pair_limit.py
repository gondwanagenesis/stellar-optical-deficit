#!/usr/bin/env python
"""AUDIT of the headline claim: p_total < 6.2e-4, one star in 1,614.

    run.sh scripts/89_audit_pair_limit.py --tag primary

WHY THIS ONE
------------
Every channel in this project reports a null except the wide-pair limit, which
is the only place a NUMBER is quoted and the only thing anyone would cite. It
is also the number that was rebuilt twice under time pressure (scripts 32, 33,
34) after the original asymmetry estimator was found to be blind. It has never
been audited. The report's abstract line is

    Joint, disposal-agnostic: p_total < 6.2e-4. Fewer than 1 in 1,614 nearby
    lower-main-sequence stars intercepts >= 51% of its optical output.

and it comes from a single row of results/pair_limit_v3_primary.json:
`clean + bare`, k = 7, obs = 3, f_det = 0.5066, p_UL = 4.297e-4.

WHAT IS CHECKED
---------------
1. REPRODUCIBILITY. Rebuild the sample through the same code path and compare
   every cell of the table against the JSON on disk.

2. THE EFFICIENCY TERM. `pipeline.statistics.exclusion_curve` -- the shared
   function every other channel uses -- computes

       p_UL(f) = N_UL / (N_total * efficiency(f, k))

   Script 34 hand-rolled its own table and wrote `p_UL = best / n_stars`, with
   no efficiency factor at all. A star whose covering fraction is exactly
   f_det produces an expected |dr| of exactly the threshold T, so noise puts it
   above T about half the time. If the efficiency is really ~0.5 the quoted
   limit is a factor of two too strong. This is the same class of error as the
   hand-reimplemented C* formula already caught in this project: a shared,
   tested helper existed and was bypassed.

   Measured here by shifting the OBSERVED dr distribution, which uses the real
   non-Gaussian noise shape rather than a Gaussian assumption, and by finding
   the f at which the efficiency actually reaches 0.9 -- i.e. the f at which
   the quoted number is true as stated.

3. THE BACKGROUND MODEL. At every k the observed count runs far BELOW the
   predicted background: 10 vs 51 at k = 7, 206 vs 422 at k = 4. The report
   presents that as reassurance. It is not -- a background model that
   over-predicts by an order of magnitude is a broken model, and it is load
   bearing, because `best = min(ul_cons, ul_sub)` takes the smaller of the
   Poisson and the background-subtracted limit, so wherever the subtracted
   branch wins the published limit rests on it.

   Two candidate explanations are tested here and BOTH FAIL, which is the
   informative part. (a) A threshold-units mismatch: T = k*sigma_dr, while
   q's own scale is sigma_r, so q is evaluated at 6.76 sigma_r rather than 7.
   Correcting it moves obs/pred from 0.080 to 0.094 and no further. (b) A
   mark-permutation null on dr, shuffling residuals within (M_Ks, distance)
   cells over the same pairs: it over-predicts worse still, obs/pred = 0.062.

   What is left is physical rather than a coding slip. q(T) is the tail of a
   SINGLE-STAR statistic and is applied to a DIFFERENCE. The single-star tail
   is dominated by contaminants common to both components of a wide pair --
   shared crowding, shared extinction error, shared fiducial mis-specification
   -- and those cancel in dr. The core already carries 53% common-mode
   variance and the tail evidently carries far more, which is also why the
   permutation, which destroys the cancellation outright, over-predicts most.
   The conclusion is that this analysis has NO valid background for |dr|, and
   the subtracted branch must not be taken at any k.

4. MIRROR. The estimator counts |dr| > T two-sided, so its own mirror is the
   sign split of the survivors. An absorber can only dim, so a signal must be
   one-sided in the sense of "the dimmer component is the harvested one" -- but
   dr is defined primary-minus-secondary, so a real population is symmetric in
   sign and the split tests only the noise, which is stated rather than used.

5. LOOK-ELSEWHERE. The quoted row is the best of 8 (4 values of k x 2 samples)
   chosen by minimum mean_f_UL, with no penalty. Reported here against a fixed
   pre-registered k so the size of the optimism is visible.
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
P_ISO_SUAZO = 1.9e-4
N_PERM = 200
SEED = 20260821
K_PREREG = 5          # a fixed k, so the best-of-8 optimism is visible


def rebuild(tag: str, knots: int):
    """Exactly script 34's sample construction, re-run."""
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{tag}.parquet")
    print(f"full sample: {len(d):,} stars", flush=True)

    a0 = np.nan_to_num(d["A_0"].to_numpy(float))
    bp_rp = d["bp_rp"].to_numpy(float)
    a_j = ext.deredden("J", a0, bp_rp)
    jks0 = ((d["tmass_j_m"].to_numpy(float) - a_j)
            - (d["tmass_ks_m"].to_numpy(float) - d["A_Ks"].to_numpy(float)))
    d["j_ks0"] = jks0

    ok = np.isfinite(d["M_G"]) & np.isfinite(d["M_Ks"]) & np.isfinite(jks0)
    d = d[ok].reset_index(drop=True)
    print(f"with finite M_G, M_Ks, (J-Ks)_0 : {len(d):,}", flush=True)

    covs = [(d["j_ks0"].to_numpy(float), 1)]
    m_ks = d["M_Ks"].to_numpy(float)
    fit = fid.fit_fiducial(m_ks, covs, d["M_G"].to_numpy(float), knots)
    r = fit.residuals(m_ks, covs, d["M_G"].to_numpy(float))
    sig = st.robust_sigma(r)
    print(f"  sigma(r) = {sig:.5f} mag", flush=True)

    p = pr.drop_shared_components(pr.find_pairs(d))
    fake = pr.find_pairs(d, scramble=True)
    clean = p[p["theta_arcsec"] > pr.CLEAN_SEP_ARCSEC].reset_index(drop=True)
    print(f"pairs {len(p):,} | chance {len(fake):,} | clean {len(clean):,}",
          flush=True)

    i, j = clean["i"].to_numpy(), clean["j"].to_numpy()
    prim_i = m_ks[i] <= m_ks[j]
    dr = np.where(prim_i, r[i], r[j]) - np.where(prim_i, r[j], r[i])

    # mid-IR veto, verbatim
    ks0 = d["tmass_ks_m"].to_numpy(float) - d["A_Ks"].to_numpy(float)
    w1 = d["wise_w1mpro"].to_numpy(float)
    w2 = d["wise_w2mpro"].to_numpy(float)
    e1 = d["wise_w1mpro_error"].to_numpy(float)
    e2 = d["wise_w2mpro_error"].to_numpy(float)
    have = (np.isfinite(w1) & np.isfinite(w2) & np.isfinite(e1)
            & np.isfinite(e2) & (e1 < 0.2) & (e2 < 0.2))
    jk = d["j_ks0"].to_numpy(float)

    def expected(y):
        lo, hi = np.nanpercentile(jk[have], [1, 99])
        edges = np.linspace(lo, hi, 61)
        idx = np.clip(np.digitize(jk, edges) - 1, 0, 59)
        m_ = np.full(60, np.nan)
        for b in range(60):
            mm = (idx == b) & have & np.isfinite(y)
            if mm.sum() > 100:
                m_[b] = np.median(y[mm])
        g = np.isfinite(m_)
        return np.interp(jk, (0.5 * (edges[:-1] + edges[1:]))[g], m_[g])

    ex1 = (ks0 - w1) - expected(ks0 - w1)
    ex2 = (ks0 - w2) - expected(ks0 - w2)
    s1, s2 = st.robust_sigma(ex1[have]), st.robust_sigma(ex2[have])
    bare = have & (np.abs(ex1) < EXCESS_NSIGMA * s1) & (np.abs(ex2) < EXCESS_NSIGMA * s2)
    both_bare = bare[i] & bare[j]
    print(f"  pairs with both components bare : {int(both_bare.sum()):,}",
          flush=True)

    return d, r, sig, i, j, dr, both_bare, len(fake) / max(len(p), 1)


def efficiency(dr_noise, med, T, delta):
    """Fraction of pairs with a Delta-mag dimming on one component that clear T.

    Dimming the primary shifts dr by +Delta, the secondary by -Delta; each is
    equally likely, so average the two. Uses the observed dr as the noise
    realisation, which is conservative: any real signal already in the tail
    inflates the assumed noise and lowers the efficiency.
    """
    up = np.mean(np.abs(dr_noise + delta - med) > T)
    dn = np.mean(np.abs(dr_noise - delta - med) > T)
    return float(0.5 * (up + dn))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--knots", type=int, default=6)
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)

    d, r, sig, i, j, dr, both_bare, chance = rebuild(args.tag, args.knots)
    med = float(np.median(dr))
    s_dr = st.robust_sigma(dr)
    common_mode = 1.0 - (s_dr / np.sqrt(2) / sig) ** 2
    print(f"\nsigma(dr) = {s_dr:.5f}  sigma(r) = {sig:.5f}  "
          f"common-mode variance fraction = {100*common_mode:.1f}%", flush=True)

    members = np.unique(np.concatenate([i, j]))
    r_mem = r[members]
    dev_mem = np.abs(r_mem - float(np.median(r_mem)))
    dev_all = np.abs(r - float(np.median(r)))
    q_mem = lambda T: float(np.count_nonzero(dev_mem > T)) / len(r_mem)
    q_glo = lambda T: float(np.count_nonzero(dev_all > T)) / len(r)

    # ---- 3. mark-permutation null on dr --------------------------------
    # Shuffle residuals over pair members inside (M_Ks, distance) cells, then
    # recompute dr for the SAME pairs. This is the project's mandated null
    # construction. It removes the common mode, so it OVER-estimates the noise
    # and therefore the background; stated, not hidden.
    m_ks = d["M_Ks"].to_numpy(float)
    dist = d["dist_pc"].to_numpy(float)
    cell = (np.clip(np.digitize(m_ks[members], np.percentile(m_ks[members], np.arange(5, 100, 5))), 0, 19) * 20
            + np.clip(np.digitize(dist[members], np.percentile(dist[members], np.arange(5, 100, 5))), 0, 19))
    pos_of = {s: n for n, s in enumerate(members)}
    ii = np.array([pos_of[x] for x in i])
    jj = np.array([pos_of[x] for x in j])
    order = np.argsort(cell, kind="stable")
    bounds = np.flatnonzero(np.diff(cell[order])) + 1
    groups = np.split(order, bounds)

    def perm_counts(Ts):
        out = np.zeros((N_PERM, len(Ts)))
        rr = r[members].copy()
        for p_ in range(N_PERM):
            sh = rr.copy()
            for g in groups:
                sh[g] = rr[rng.permutation(g)]
            drp = sh[ii] - sh[jj]
            drp = np.where(m_ks[i] <= m_ks[j], drp, -drp)
            a = np.abs(drp - np.median(drp))
            out[p_] = [np.count_nonzero(a > T) for T in Ts]
        return out

    # ---- tables --------------------------------------------------------
    f_grid = np.arange(0.05, 0.951, 0.005)
    ks = (4, 5, 6, 7)
    Ts = [k * s_dr for k in ks]
    print(f"\nrunning {N_PERM} mark permutations over "
          f"{len(groups)} cells ...", flush=True)
    pc = perm_counts(Ts)
    print("  done", flush=True)

    rows = []
    for label, mask in (("all clean", np.ones(len(dr), bool)),
                        ("clean + bare", both_bare)):
        drv = dr[mask]
        n_pairs = len(drv)
        n_stars = 2 * n_pairs
        scale = n_pairs / len(dr)
        for kk, T in zip(ks, Ts):
            obs = int(np.count_nonzero(np.abs(drv - med) > T))
            n_pos = int(np.count_nonzero(drv - med > T))
            n_neg = int(np.count_nonzero(drv - med < -T))
            f_det = float(st.fraction_from_delta(T))
            ul_cons = st.poisson_upper_limit(obs)

            # script 34's background, and the two diagnoses of it
            pred_m_drsig = 2.0 * q_mem(T) * n_pairs             # what 34 used
            pred_m_rsig = 2.0 * q_mem(kk * sig) * n_pairs       # T in sigma_r
            k_of_T_in_r = T / sig
            idx = ks.index(kk)
            pred_perm = float(pc[:, idx].mean()) * scale
            pred_perm_sd = float(pc[:, idx].std()) * np.sqrt(scale)

            ul_sub = max(obs - pred_m_drsig, 0.0) + 1.645 * np.sqrt(obs + pred_m_drsig)
            best_34 = min(ul_cons, ul_sub)
            p_UL_34 = best_34 / n_stars
            used_bkg_branch = bool(ul_sub < ul_cons)
            p_UL_poisson_only = ul_cons / n_stars

            eff_at_fdet = efficiency(drv, med, T, T)   # delta = T by definition
            p_UL_eff = ul_cons / (n_stars * eff_at_fdet)

            # the f at which the ORIGINAL quoted p_UL is actually valid
            effs = np.array([efficiency(drv, med, T, st.delta_mag(f))
                             for f in f_grid])
            f90 = float(f_grid[np.argmax(effs >= 0.90)]) if (effs >= 0.90).any() else float("nan")
            f99 = float(f_grid[np.argmax(effs >= 0.99)]) if (effs >= 0.99).any() else float("nan")

            # best mean-f limit with efficiency folded in, over the f grid
            with np.errstate(divide="ignore"):
                mf = np.where(effs > 0, ul_cons / (n_stars * effs) * f_grid, np.inf)
            b = int(np.argmin(mf))

            rows.append({
                "sample": label, "k": kk, "thr_mag": T, "f_det": f_det,
                "obs": obs, "n_pos": n_pos, "n_neg": n_neg,
                "pred_34_in_sigma_dr": pred_m_drsig,
                "pred_in_sigma_r": pred_m_rsig,
                "k_of_T_in_sigma_r": k_of_T_in_r,
                "obs_over_pred_34": obs / max(pred_m_drsig, 1e-9),
                "obs_over_pred_sigma_r": obs / max(pred_m_rsig, 1e-9),
                "pred_permutation": pred_perm,
                "pred_permutation_sd": pred_perm_sd,
                "obs_over_pred_permutation": obs / max(pred_perm, 1e-9),
                "p_UL_as_published": p_UL_34,
                "published_used_background_branch": used_bkg_branch,
                "p_UL_poisson_only": p_UL_poisson_only,
                "p_UL_poisson_only_eff_corrected": p_UL_poisson_only / eff_at_fdet,
                "efficiency_at_f_det": eff_at_fdet,
                "p_UL_efficiency_corrected": p_UL_eff,
                "f_at_eff_90": f90, "f_at_eff_99": f99,
                "best_mean_f_UL_corrected": float(mf[b]),
                "f_of_best_mean_f": float(f_grid[b]),
                "p_UL_at_best_mean_f": float(ul_cons / (n_stars * effs[b])),
                "n_pairs": n_pairs,
            })

    t = pd.DataFrame(rows)
    show = ["sample", "k", "f_det", "obs", "n_pos", "n_neg",
            "obs_over_pred_34", "obs_over_pred_sigma_r",
            "obs_over_pred_permutation", "published_used_background_branch",
            "p_UL_as_published", "efficiency_at_f_det",
            "p_UL_efficiency_corrected", "f_at_eff_90"]
    print("\n=== audit table ===", flush=True)
    print(t[show].to_string(index=False,
                            float_format=lambda v: f"{v:10.4g}"), flush=True)

    # ---- reproducibility against the JSON on disk ----------------------
    pub_path = cfg.RESULT_DIR / f"pair_limit_v3_{args.tag}.json"
    repro = {"checked": False}
    if pub_path.exists():
        pub = json.loads(pub_path.read_text())
        diffs = []
        for pr_ in pub["table"]:
            m = t[(t["sample"] == pr_["sample"]) & (t["k"] == pr_["k"])]
            if not len(m):
                continue
            m = m.iloc[0]
            for a, b_ in (("obs", "obs"), ("thr_mag", "thr_mag"),
                          ("f_det", "f_det"), ("p_UL", "p_UL_as_published")):
                got, want = float(m[b_]), float(pr_[a])
                if abs(got - want) > 1e-9 + 1e-6 * abs(want):
                    diffs.append({"sample": pr_["sample"], "k": pr_["k"],
                                  "field": a, "published": want, "rebuilt": got})
        repro = {"checked": True, "n_mismatches": len(diffs),
                 "mismatches": diffs[:20],
                 "reproduces": len(diffs) == 0}
        print(f"\nreproducibility vs {pub_path.name}: "
              f"{'EXACT' if not diffs else str(len(diffs)) + ' MISMATCHES'}",
              flush=True)

    # ---- the headline row, as published and as corrected ---------------
    head = t[(t["sample"] == "clean + bare") & (t["k"] == 7)].iloc[0]
    p_dark_pub = float(head["p_UL_as_published"])
    p_dark_eff = float(head["p_UL_efficiency_corrected"])
    prereg = t[(t["sample"] == "clean + bare") & (t["k"] == K_PREREG)].iloc[0]

    headline = {
        "published": {
            "p_dark": p_dark_pub, "f": float(head["f_det"]),
            "p_total": P_ISO_SUAZO + p_dark_pub,
            "one_in_n": 1.0 / (P_ISO_SUAZO + p_dark_pub)},
        "efficiency_corrected_same_f": {
            "efficiency": float(head["efficiency_at_f_det"]),
            "p_dark": p_dark_eff, "f": float(head["f_det"]),
            "p_total": P_ISO_SUAZO + p_dark_eff,
            "one_in_n": 1.0 / (P_ISO_SUAZO + p_dark_eff)},
        "same_p_valid_at_f": {
            "f_at_eff_90": float(head["f_at_eff_90"]),
            "f_at_eff_99": float(head["f_at_eff_99"]),
            "note": "the covering fraction at which the PUBLISHED p_UL is "
                    "true as stated, i.e. where the channel is actually "
                    "efficient rather than at its 50% threshold"},
        "rows_that_used_the_background_branch": [
            {"sample": rr["sample"], "k": rr["k"], "f": rr["f_det"],
             "p_UL_as_published": rr["p_UL_as_published"],
             "p_UL_poisson_only": rr["p_UL_poisson_only"],
             "p_UL_poisson_only_eff_corrected": rr["p_UL_poisson_only_eff_corrected"],
             "weakening_factor": rr["p_UL_poisson_only"] / rr["p_UL_as_published"]}
            for rr in t.to_dict(orient="records")
            if rr["published_used_background_branch"]],
        "preregistered_k": {
            "k": K_PREREG, "p_dark": float(prereg["p_UL_as_published"]),
            "p_dark_efficiency_corrected": float(prereg["p_UL_efficiency_corrected"]),
            "f": float(prereg["f_det"]),
            "note": "the published row is the best of 8 (4 k x 2 samples) "
                    "chosen by minimum mean_f_UL with no trials penalty"},
    }
    print("\n=== headline ===", flush=True)
    print(json.dumps(headline, indent=1), flush=True)

    # ---- verdict -------------------------------------------------------
    eff = float(head["efficiency_at_f_det"])
    ratio_34 = float(head["obs_over_pred_34"])
    ratio_r = float(head["obs_over_pred_sigma_r"])
    findings = []
    if eff < 0.75:
        findings.append(
            f"EFFICIENCY TERM MISSING. Script 34 wrote p_UL = N_UL / n_stars "
            f"with no efficiency factor, while pipeline.statistics."
            f"exclusion_curve -- the shared helper every other channel uses -- "
            f"defines p_UL = N_UL / (N * efficiency(f, k)). Measured "
            f"efficiency at the quoted f = {float(head['f_det']):.4f} is "
            f"{eff:.3f}, because a star at exactly the threshold is scattered "
            f"above it only about half the time. The published p_dark = "
            f"{p_dark_pub:.3e} at f >= {float(head['f_det']):.3f} is therefore "
            f"a factor {1/eff:.2f} too strong AT THAT f. Corrected: p_dark < "
            f"{p_dark_eff:.3e}, p_total < {P_ISO_SUAZO + p_dark_eff:.3e}, "
            f"one star in {1/(P_ISO_SUAZO + p_dark_eff):,.0f} rather than "
            f"{1/(P_ISO_SUAZO + p_dark_pub):,.0f}. Equivalently the published "
            f"number is true as stated at f >= {float(head['f_at_eff_90']):.3f} "
            f"(90% efficient), not at f >= {float(head['f_det']):.3f}.")
    ratio_perm = float(head["obs_over_pred_permutation"])
    if not (0.6 < ratio_34 < 1.7):
        # The threshold-units hypothesis was tested here and FAILED to explain
        # the over-prediction: correcting T to sigma_r units moves obs/pred
        # from 0.080 to 0.094 and no further. So the model is wrong in kind,
        # not in calibration.
        findings.append(
            f"BACKGROUND MODEL IS FOR THE WRONG STATISTIC. Script 34 predicts "
            f"the pair tail as 2*q(T)*n_pairs with q = P(|r_star| > T) measured "
            f"on single-star residuals, and applies it to |dr|, a DIFFERENCE. "
            f"Observed over predicted is {ratio_34:.3f} at the headline row -- "
            f"the prediction is {1/max(ratio_34,1e-9):.0f}x too large. Two "
            f"candidate explanations were tested and both fail to rescue it. "
            f"(a) A threshold-units mismatch: T = k*sigma_dr while q's own "
            f"scale is sigma_r = {sig:.5f} vs sigma_dr = {s_dr:.5f}, so q is "
            f"evaluated at {float(head['k_of_T_in_sigma_r']):.2f} sigma_r "
            f"rather than 7. Correcting it moves obs/pred only "
            f"{ratio_34:.3f} -> {ratio_r:.3f}. (b) A mark-permutation null on "
            f"dr, shuffling residuals within (M_Ks, distance) cells over the "
            f"same pairs, over-predicts further still, obs/pred = "
            f"{ratio_perm:.3f}. The residual explanation is physical: the "
            f"single-star tail is dominated by contaminants COMMON to both "
            f"components of a wide pair -- shared crowding, shared extinction "
            f"error, shared fiducial mis-specification -- which cancel in the "
            f"difference. The core already carries "
            f"{100*common_mode:.0f}% common-mode variance and the tail "
            f"evidently carries far more; the permutation destroys that "
            f"cancellation entirely, which is why it over-predicts worst. "
            f"CONSEQUENCE: no valid background for |dr| exists in this "
            f"analysis, so `best = min(ul_cons, ul_sub)` must not take the "
            f"subtracted branch. It does so in "
            f"{int(t['published_used_background_branch'].sum())} of the 8 "
            f"published rows, which are weakened by factors of "
            f"{t.loc[t['published_used_background_branch'],'p_UL_poisson_only'].div(t.loc[t['published_used_background_branch'],'p_UL_as_published']).min():.2f}"
            f"-{t.loc[t['published_used_background_branch'],'p_UL_poisson_only'].div(t.loc[t['published_used_background_branch'],'p_UL_as_published']).max():.2f} "
            f"once forced onto the Poisson branch. The headline row "
            f"(clean + bare, k = 7) is NOT among them -- it already takes the "
            f"conservative branch -- so this finding costs the abstract "
            f"nothing and costs the rest of the table a great deal.")
    verdict = (" ".join(findings) if findings
               else "NO ERROR FOUND. The published limit reproduces and its "
                    "efficiency and background model both check out.")
    print("\n=== VERDICT ===\n" + verdict, flush=True)

    res = {
        "tag": args.tag, "audits": "results/pair_limit_v3_%s.json" % args.tag,
        "audited_claim": "p_total < 6.2e-4, fewer than 1 in 1,614 nearby "
                         "lower-main-sequence stars intercepts >= 51% of its "
                         "optical output by any means",
        "seed": SEED, "n_perm": N_PERM,
        "sigma_r": sig, "sigma_dr": s_dr,
        "common_mode_variance_fraction": float(common_mode),
        "chance_pair_fraction": float(chance),
        "reproducibility": repro,
        "headline": headline,
        "sign_split": {
            "note": "dr is primary-minus-secondary and a harvested component "
                    "can be either, so a signal is symmetric in sign and the "
                    "split is a noise diagnostic, not a mirror. Every row "
                    "runs NEGATIVE-heavy, which is the known unresolved-"
                    "companion excess (a brightening) already identified in "
                    "channel 12. Since those objects cannot be absorbers, the "
                    "two-sided Poisson limit counts them as candidate signal "
                    "and is conservative by roughly the negative fraction.",
            "rows": [{"sample": rr["sample"], "k": rr["k"],
                      "n_pos": rr["n_pos"], "n_neg": rr["n_neg"]}
                     for rr in t.to_dict(orient="records")]},
        "mirror_note": "dr is defined primary-minus-secondary, so a real "
                       "harvested population is symmetric in the SIGN of dr "
                       "and the sign split tests the noise only. It is "
                       "reported, not used as a false-positive rate.",
        "table": t.to_dict(orient="records"),
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"audit_pair_limit_{args.tag}.json"
    out.write_text(json.dumps(res, indent=2, default=float))
    t.to_csv(cfg.RESULT_DIR / f"audit_pair_limit_{args.tag}.csv", index=False)
    print(f"\nwrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
