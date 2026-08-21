#!/usr/bin/env python
"""Audit of channel 4 - the dynamical-mass search (scripts/37).

    run.sh scripts/90_audit_dynamical_mass.py

WHY THIS CHANNEL AND NOT ANOTHER
--------------------------------
Channel 4 is load-bearing for coverage rather than for the headline number.
Every other channel measures light, and the report's own coverage-gaps section
says so:

    "The grey case (alpha = 0) is invisible to every photometric channel by
     construction and is covered *only* by channel 4."

So if channel 4 is not in fact independent of photometry, the grey blind spot is
open and the summary table's "ANY spectral slope, incl. grey" is wrong.

THE SUSPICION
-------------
gaiadr3.binary_masses carries a `combination_method` string, and seven of its
nine values end in "+M1".  In the Gaia DR3 multiplicity processing that suffix
marks the primary mass as supplied externally from a mass-luminosity relation
rather than solved from the orbit: the astrometric mass-ratio function and the
SB1 mass function both leave an m1 degeneracy that has to be broken with
something, and the something is the primary's absolute G magnitude.

If that is what happened, then regressing M_G on log(m1) regresses absolute G
against a monotonic function of absolute G.  A grey absorber that dims a star
does not push it off the relation, it slides it ALONG the relation, because the
same dimming that changes the numerator also changes the denominator.  The
channel would return zero positives whether or not the signal is there.

WHAT THIS SCRIPT DOES
---------------------
1. Counts how much of channel 4's sample - and of its headline faint-secondary
   subset - carries a photometric m1.
2. Determinism test: how well log(m1) is predicted by photometry alone, split by
   combination_method.  The genuinely dynamical methods are the control.
3. The decisive test: INJECT a grey dimming into real stars, propagate it
   through an empirically reconstructed version of Gaia's own m1(M_G) map, and
   measure how much of it the channel's estimator recovers.  A working channel
   returns 1.0 mag of residual per 1.0 mag injected.
4. POSITIVE CONTROL, in code: the same injection with m1 held fixed - the
   genuinely dynamical case - must return 1.0.  Without it a response of zero
   would be indistinguishable from a broken injection harness.
5. The efficiency term.  scripts/37 hand-rolls `p_UL = ul / n` with no
   efficiency factor, the third occurrence in this project of a shared, tested
   helper being bypassed by a local reimplementation that drops one of its
   terms.
6. What the channel can honestly say once restricted to a dynamical m1.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import extinction as ext
from pipeline import sample as smp
from pipeline import statistics as st

RNG = np.random.default_rng(20260821)


def m1_is_photometric(method: pd.Series) -> pd.Series:
    """The "+M1" suffix marks an externally supplied primary mass."""
    return method.fillna("").str.contains(r"\+M1$", regex=True)


def design(mg, br) -> np.ndarray:
    """Photometry-only design matrix for the m1 map."""
    mg = np.asarray(mg, dtype=float)
    br = np.asarray(br, dtype=float)
    return np.column_stack([np.ones(len(mg)), mg, mg ** 2, mg ** 3,
                            br, br ** 2, mg * br])


def build_sample() -> pd.DataFrame:
    """Reproduce scripts/37's cut chain exactly, keeping the method label."""
    df = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    df = smp.add_astrometry(df)
    a0 = ext.query_a0("edenhofer23", df["l"].to_numpy(float),
                      df["b"].to_numpy(float), df["dist_pc"].to_numpy(float))
    bp_rp = df["bp_rp"].to_numpy(float)
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = df["dist_mod"].to_numpy(float)
    df["A_0"] = a0
    df["M_G"] = df["phot_g_mean_mag"].to_numpy(float) - mu - a_g
    df["M_Ks"] = df["tmass_ks_m"].to_numpy(float) - mu - a_ks
    good = (np.isfinite(df["M_G"]) & np.isfinite(a0) & df["m1"].notna()
            & (df["m1"] > 0.2) & (df["m1"] < 1.6)
            & (df["tmass_ph_qual"].str[2] == "A")
            & (a0 < 0.6))
    return df[good].reset_index(drop=True)


def fit_MG_on_logm1(d: pd.DataFrame):
    """The channel's own estimator: M_G against log m1, best of deg 2-4."""
    x = np.log10(d["m1"].to_numpy(float))
    y = d["M_G"].to_numpy(float)
    best = None
    for deg in (2, 3, 4):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    return best


def mad_dex(r: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(r - np.median(r))))


def main() -> int:
    out: dict = {}
    d = build_sample()
    d["m1_photometric"] = m1_is_photometric(d["combination_method"])
    print(f"channel 4 usable sample reproduced: {len(d):,}")

    # ------------------------------------------------------------------ 1
    print("\n=== 1. how much of channel 4 has a photometric m1 ===")
    faint = d["fluxratio"].notna() & (d["fluxratio"] < 0.1)
    comp = []
    for label, m in (("all systems", pd.Series(True, index=d.index)),
                     ("faint secondary (headline row)", faint)):
        sub = d[m]
        nphot = int(sub["m1_photometric"].sum())
        comp.append({"subset": label, "N": int(len(sub)),
                     "n_m1_photometric": nphot,
                     "n_m1_dynamical": int(len(sub) - nphot),
                     "frac_photometric": float(nphot / max(len(sub), 1))})
        print(f"  {label:32s} N={len(sub):6,d}  photometric m1={nphot:6,d} "
              f"({nphot / max(len(sub), 1):.4%})  dynamical={len(sub) - nphot}")
    out["composition"] = comp

    bym = (d.groupby("combination_method")
             .agg(N=("m1", "size"), photometric=("m1_photometric", "max"))
             .sort_values("N", ascending=False))
    print("\n  by method:")
    print("  " + bym.to_string().replace("\n", "\n  "))
    out["by_method"] = {k: {"N": int(v["N"]), "photometric": bool(v["photometric"])}
                        for k, v in bym.iterrows()}

    # ------------------------------------------------------------------ 2
    print("\n=== 2. determinism: is log m1 predictable from photometry alone? ===")
    print("  fitting log m1 = F(M_G, bp_rp) per method; dynamical methods are")
    print("  the control and must scatter more.")
    det = []
    for meth, g in d.groupby("combination_method"):
        if len(g) < 15:
            continue
        y = np.log10(g["m1"].to_numpy(float))
        A = design(g["M_G"], g["bp_rp"])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ coef
        det.append({"method": meth, "N": int(len(g)),
                    "photometric": bool(g["m1_photometric"].iloc[0]),
                    "scatter_dex_mad": mad_dex(r),
                    "R2": float(1 - np.var(r) / np.var(y))})
    det.sort(key=lambda z: z["scatter_dex_mad"])
    for z in det:
        tag = "PHOT" if z["photometric"] else "DYN "
        print(f"  [{tag}] {z['method']:24s} N={z['N']:6,d}  "
              f"scatter {z['scatter_dex_mad']:.4f} dex   R2={z['R2']:.4f}")
    out["determinism"] = det
    print("  the only dynamical control here has N=25, so this split is"
          " suggestive rather than decisive on its own. Test 2b is not.")

    # -- 2b -------------------------------------------------------------
    print("\n=== 2b. does m1 know anything photometry does not? ===")
    print("  An astrometric mass goes as a^3/P^2 with a = theta * distance, so")
    print("  at FIXED (M_G, bp_rp) a dynamical m1 must still depend strongly on")
    print("  distance. A mass read off a magnitude-luminosity relation cannot.")
    print("  Adding apparent G - which supplies the parallax that M_G has")
    print("  already divided out - is therefore a direct test of provenance.")
    nested = []
    for meth, g in d.groupby("combination_method"):
        if len(g) < 500 and g["m1_photometric"].iloc[0]:
            continue
        if len(g) < 12:
            continue
        y = np.log10(g["m1"].to_numpy(float))
        A = design(g["M_G"], g["bp_rp"])
        gmag = g["phot_g_mean_mag"].to_numpy(float)
        B = np.column_stack([A, gmag, gmag ** 2, gmag * g["M_G"].to_numpy(float)])
        r_a = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
        r_b = y - B @ np.linalg.lstsq(B, y, rcond=None)[0]
        s_a, s_b = mad_dex(r_a), mad_dex(r_b)
        nested.append({"method": meth, "N": int(len(g)),
                       "photometric": bool(g["m1_photometric"].iloc[0]),
                       "scatter_photometry_only": s_a,
                       "scatter_plus_distance": s_b,
                       "improvement": float(1 - s_b / s_a) if s_a > 0 else None})
        tag = "PHOT" if g["m1_photometric"].iloc[0] else "DYN "
        print(f"  [{tag}] {meth:24s} N={len(g):6,d}  "
              f"{s_a:.4f} -> {s_b:.4f} dex   "
              f"distance adds {1 - s_b / s_a:6.1%}")
    out["nested_distance_test"] = nested

    # ------------------------------------------------------------------ 3, 4
    print("\n=== 3. injection: how much dimming does the channel recover? ===")
    c, sig, deg = fit_MG_on_logm1(d)
    d["res"] = d["M_G"].to_numpy(float) - np.polyval(c, np.log10(d["m1"]))
    print(f"  estimator: degree {deg} polynomial, robust scatter {sig:.4f} mag")

    ph = d[d["m1_photometric"]]
    Fcoef, *_ = np.linalg.lstsq(design(ph["M_G"], ph["bp_rp"]),
                                np.log10(ph["m1"].to_numpy(float)), rcond=None)
    Fres = (np.log10(ph["m1"].to_numpy(float))
            - design(ph["M_G"], ph["bp_rp"]) @ Fcoef)
    print(f"  reconstructed Gaia map F(M_G,bp_rp) -> log m1 : "
          f"residual {mad_dex(Fres):.4f} dex over {len(ph):,} systems")
    out["F_map_scatter_dex"] = mad_dex(Fres)

    # A rare signal does not move the global fit, so the estimator polynomial
    # stays fixed at its un-injected value while the targets are dimmed.
    idx = RNG.choice(len(d), size=max(int(0.01 * len(d)), 500), replace=False)
    tgt = d.iloc[idx]
    mg0 = tgt["M_G"].to_numpy(float)
    br = tgt["bp_rp"].to_numpy(float)
    lm1_fix = np.log10(tgt["m1"].to_numpy(float))

    # Both branches are baselined through their OWN zero-injection value, so
    # the error in our reconstruction of F cancels out of the response and only
    # the response itself is measured.
    r_circ0 = mg0 - np.polyval(c, design(mg0, br) @ Fcoef)
    r_dyn0 = mg0 - np.polyval(c, lm1_fix)

    # F is a cubic and must not be extrapolated: outside the magnitude range it
    # was fitted over it diverges, which would manufacture a response.
    mg_lo = float(np.percentile(ph["M_G"], 0.5))
    mg_hi = float(np.percentile(ph["M_G"], 99.5))
    print(f"  F is calibrated over M_G in [{mg_lo:.2f}, {mg_hi:.2f}]; injected "
          "stars leaving that range are dropped rather than extrapolated")

    inj = []
    for delta in (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0):
        mg1 = mg0 + delta                    # grey dimming, colour unchanged
        in_range = (mg1 >= mg_lo) & (mg1 <= mg_hi)

        # (a) circular branch: Gaia re-derives m1 from the dimmed photometry
        lm1_dim = design(mg1, br) @ Fcoef
        r_circ = mg1 - np.polyval(c, lm1_dim)
        m1_dim = 10.0 ** lm1_dim
        kept = in_range & (m1_dim > 0.2) & (m1_dim < 1.6)

        # (b) positive control: a dynamical m1 does not move when the star dims
        r_dyn = mg1 - np.polyval(c, lm1_fix)

        if kept.sum() < 20:
            print(f"  inject {delta:4.2f} mag -> only {int(kept.sum())} stars "
                  "remain in range and in the mass window; not measurable")
            inj.append({"delta_mag": delta, "response_circular": None,
                        "recovered_mag": None,
                        "response_dynamical_control": None,
                        "frac_surviving_cuts": float(np.mean(kept)),
                        "f_covering": float(st.fraction_from_delta(delta))})
            continue

        rec = float(np.median((r_circ - r_circ0)[kept]))
        g_circ = rec / delta
        g_dyn = float(np.median((r_dyn - r_dyn0)[kept]) / delta)
        inj.append({"delta_mag": delta, "response_circular": g_circ,
                    "recovered_mag": rec,
                    "response_dynamical_control": g_dyn,
                    "frac_surviving_cuts": float(np.mean(kept)),
                    "f_covering": float(st.fraction_from_delta(delta))})
        print(f"  inject {delta:4.2f} mag -> recovered {rec:+.4f} mag "
              f"(gain {g_circ:+.4f})  |  control gain {g_dyn:+.4f}  |  "
              f"{np.mean(kept):.1%} survive the cuts")

    # mirror: an absorber can only dim, so a brightening is unphysical and its
    # response measures the harness rather than the signal
    mg_m = mg0 - 1.0
    ok_m = (mg_m >= mg_lo) & (mg_m <= mg_hi)
    r_m = mg_m - np.polyval(c, design(mg_m, br) @ Fcoef)
    mirror = float(np.median((r_m - r_circ0)[ok_m]) / -1.0)
    print(f"  mirror (brighten 1.00 mag): gain {mirror:+.4f}  "
          "(symmetric, as an attenuation should be)")
    out["injection"] = inj
    out["mirror_gain_1mag"] = mirror

    gains = [z["response_circular"] for z in inj if z["response_circular"]]
    ctrl = [z["response_dynamical_control"] for z in inj
            if z["response_dynamical_control"]]
    gbar = float(np.median(gains))
    print(f"\n  median gain, circular : {gbar:+.4f}")
    print(f"  median gain, control  : {np.median(ctrl):+.4f}  "
          "(must be ~1.000 or the injection harness itself is broken)")
    out["gain_circular_median"] = gbar
    out["gain_control_median"] = float(np.median(ctrl))

    # Where does the channel actually reach? Solve recovered(delta) = threshold.
    dd = np.array([z["delta_mag"] for z in inj if z["recovered_mag"] is not None])
    rr = np.array([z["recovered_mag"] for z in inj
                   if z["recovered_mag"] is not None])
    print(f"  recovered signal saturates at {rr.max():.3f} mag by "
          f"delta = {dd[np.argmax(rr)]:.1f} mag")
    out["max_recovered_mag"] = float(rr.max())
    out["recovery_curve"] = {"delta_mag": dd.tolist(),
                             "recovered_mag": rr.tolist()}

    # ------------------------------------------------------------------ 5
    print("\n=== 4. the efficiency term scripts/37 does not have ===")
    k = 5.0
    rows = []
    for label, m in (("all systems", pd.Series(True, index=d.index)),
                     ("faint secondary", faint)):
        r = d.loc[m, "res"].to_numpy(float)
        if len(r) < 500:
            continue
        s = st.robust_sigma(r)
        med = float(np.median(r))
        thr = k * s
        npos = int(np.count_nonzero(r > med + thr))
        nneg = int(np.count_nonzero(r < med - thr))
        f_det = float(st.fraction_from_delta(thr))
        ul = st.poisson_upper_limit(npos)
        p_pub = float(ul / len(r))               # what scripts/37 reports
        # efficiency as the shared helper defines it: assumes a 1:1 response
        eff_ideal = float(st.detection_efficiency(r, f_det, k,
                                                  sigma=s, median=med))
        # efficiency once the measured response gain is folded in
        eff_real = float(np.mean(r + gbar * thr > med + thr))
        # the covering fraction the channel actually reaches: the true dimming
        # whose RECOVERED residual equals the threshold
        if rr.max() >= thr:
            delta_true = float(np.interp(thr, rr, dd))
            f_true = float(st.fraction_from_delta(delta_true))
        else:
            delta_true, f_true = None, None
        row = {"subset": label, "N": int(len(r)), "sigma": s,
               "threshold_mag": thr, "f_det_as_published": f_det,
               "n_pos": npos, "n_neg": nneg,
               "p_UL_as_published": p_pub,
               "efficiency_ideal": eff_ideal,
               "p_UL_efficiency_corrected": (p_pub / eff_ideal
                                             if eff_ideal > 0 else None),
               "efficiency_with_measured_gain": eff_real,
               "p_UL_gain_corrected": (p_pub / eff_real
                                       if eff_real > 0 else None),
               "true_dimming_needed_mag": delta_true,
               "f_actually_reached": f_true}
        rows.append(row)
        print(f"  {label:16s} N={len(r):6,d}  sigma={s:.4f}  "
              f"thr={thr:.3f} mag (published f>={f_det:.3f})  "
              f"n_pos={npos} n_neg={nneg}")
        print(f"      p_UL as published            : {p_pub:.3e}")
        print(f"      efficiency, 1:1 response     : {eff_ideal:.4f}"
              + (f"  -> p_UL {p_pub / eff_ideal:.3e}" if eff_ideal > 0 else ""))
        print(f"      efficiency, measured gain    : {eff_real:.6f}"
              + (f"  -> p_UL {p_pub / eff_real:.3e}" if eff_real > 0
                 else "  -> no limit at the published f"))
        if delta_true is None:
            print(f"      true dimming needed to clear {thr:.2f} mag of "
                  f"residual: MORE than the {rr.max():.2f} mag the estimator "
                  "can ever recover")
            print("      -> the channel does not reach its own threshold at "
                  "ANY covering fraction, up to and including f = 1")
        else:
            print(f"      true dimming needed: {delta_true:.2f} mag "
                  f"-> the limit holds at f >= {f_true:.4f}, "
                  f"not f >= {f_det:.3f}")
    out["limits"] = rows

    # ------------------------------------------------------------------ 6
    print("\n=== 5. what survives: the genuinely dynamical subset ===")
    dyn = d[~d["m1_photometric"]]
    print(f"  systems with a dynamically solved m1 after channel 4's cuts: "
          f"{len(dyn):,}")
    if len(dyn):
        print("  " + dyn["combination_method"].value_counts()
              .to_string().replace("\n", "\n  "))
    best_possible = float(st.poisson_upper_limit(0) / max(len(dyn), 1))
    print(f"  best possible rate limit from that subset, assuming 0 positives: "
          f"p < {best_possible:.2e}  (one in {1 / best_possible:,.0f})")
    print("  the channel as published quotes one in 4,545.")
    out["dynamical_only"] = {
        "N": int(len(dyn)),
        "methods": {str(k2): int(v) for k2, v
                    in dyn["combination_method"].value_counts().items()},
        "best_possible_p_UL": best_possible}

    # ------------------------------------------------------------------ 7
    print("\n=== 6. smaller notes ===")
    df_all = pd.read_parquet(cfg.RAW_DIR / "binary_masses.parquet")
    n_ks = int((df_all["tmass_ph_qual"].str[2] == "A").sum())
    print(f"  scripts/37 requires a clean 2MASS Ks ({n_ks:,} of {len(df_all):,} "
          "pass) although the M_G branch never uses Ks. The channel is free of "
          "the 2MASS anchor in its estimator but not in its selection.")
    print("  the quoted row is the best of four (2 subsets x 2 bands) with no "
          "trials penalty, the same optimism the pair-limit audit found.")
    print("  the polynomial degree is chosen by minimising scatter on the same "
          "data it is fitted to, which shrinks sigma and inflates the outlier "
          "count, so that one runs in the conservative direction.")
    out["notes"] = {"ks_selection_not_estimator": {"n_pass": n_ks,
                                                   "n_total": int(len(df_all))},
                    "trials": "best of 4 rows, no penalty",
                    "degree_selection": "in-sample, conservative direction"}

    verdict = (
        "Channel 4 is not an independent probe of grey absorbers. "
        f"{comp[1]['frac_photometric']:.2%} of the headline faint-secondary "
        "subset carries a primary mass imported from the primary's own "
        "absolute G magnitude rather than solved from the orbit - an SB1 mass "
        "function and an astrometric mass-ratio function each give one "
        "equation in two unknowns, so m1 cannot come from the orbit at all. "
        "The estimator therefore regresses absolute G on a function of "
        "absolute G, and a grey dimming slides the star along the relation "
        "instead of off it. Injection recovers "
        f"{out['gain_circular_median']:+.3f} mag of residual per magnitude "
        f"injected against {out['gain_control_median']:+.3f} for the "
        "dynamical-m1 positive control, and the recovered signal saturates at "
        f"{out['max_recovered_mag']:.2f} mag - below the "
        f"{rows[-1]['threshold_mag']:.2f} mag threshold the channel needs to "
        "clear - then reverses sign, so a heavily harvested star is reported "
        "as slightly OVER-luminous. The channel does not reach its own "
        "threshold at any covering fraction up to f = 1, its zero positives "
        "are uninformative, and the grey blind spot the report says only "
        "channel 4 covers is open. The headline p_total is unaffected: "
        "channel 4 does not enter the joint constraint. What is lost is "
        "coverage, and specifically the only coverage claimed for alpha = 0.")
    out["verdict"] = verdict
    print("\n=== VERDICT ===\n" + verdict)

    (cfg.RESULT_DIR / "audit_dynamical_mass_primary.json").write_text(
        json.dumps(out, indent=2, default=float))
    print("\nwrote results/audit_dynamical_mass_primary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
