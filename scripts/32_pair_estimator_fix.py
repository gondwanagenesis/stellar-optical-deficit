#!/usr/bin/env python
"""Is the wide-pair ASYMMETRY estimator sensitive to the signal at all?

    run.sh scripts/32_pair_estimator_fix.py --tag primary

THE WORRY
---------
scripts/25 and 30 set limits from the asymmetry n_pos - n_neg of
dr = r_primary - r_secondary, justified by the measurement that the background
is symmetric beyond 10 arcsec.  A symmetric background is necessary for that
estimator but not sufficient: the SIGNAL must be asymmetric too.

If harvesting strikes either component with equal probability -- which is the
default assumption absent a reason to prefer one -- then a fraction p/2 of
pairs shift by +Delta and p/2 by -Delta.  The asymmetry is then zero FOR THE
SIGNAL AS WELL, and the estimator has no sensitivity whatever.

This script injects a symmetric signal and measures the response of both
estimators.  If the asymmetry estimator does not respond, the limits in
sections 5.6-5.8 of the paper are wrong and must be rebuilt on the two-sided
count.

THE REPLACEMENT
---------------
Two-sided counting has no free background subtraction, so the honest version is
conservative: every pair in either tail is allowed to be signal.  The pair
sample still wins over single stars, but for a different reason than claimed --
common-mode cancellation removes much of the BACKGROUND, not just the scatter.
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
from pipeline import statistics as st

CLEAN_SEP_ARCSEC = 10.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-real", type=int, default=200)
    args = ap.parse_args()
    rng = np.random.default_rng(31415)

    p = pd.read_parquet(cfg.DERIVED_DIR / f"wide_binaries_{args.tag}.parquet")
    clean = p[p["theta_arcsec"] > CLEAN_SEP_ARCSEC].reset_index(drop=True)
    dr0 = clean["dr"].to_numpy(float)
    n_pairs = len(dr0)
    sigma = st.robust_sigma(dr0)
    med = float(np.median(dr0))
    print(f"clean pairs: {n_pairs:,}   sigma(dr) = {sigma:.5f} mag\n")

    def stats(dr, k):
        npos = int(np.count_nonzero(dr > med + k * sigma))
        nneg = int(np.count_nonzero(dr < med - k * sigma))
        return npos, nneg, npos - nneg, npos + nneg

    print("=== injection test: does each estimator respond? ===")
    print("injecting a fraction p of pairs with ONE component dimmed by "
          "Delta(f), component chosen at random\n")
    rows = []
    for f in (0.3, 0.5):
        delta = float(st.delta_mag(f))
        for p_inj in (0.0, 1e-3, 5e-3, 2e-2):
            asym, tot = [], []
            for _ in range(args.n_real):
                dr = dr0.copy()
                hit = rng.random(n_pairs) < p_inj
                # which component is harvested is a coin flip -> symmetric
                sign = rng.choice([-1.0, 1.0], size=n_pairs)
                dr[hit] += delta * sign[hit]
                _, _, a, t = stats(dr, 5.0)
                asym.append(a); tot.append(t)
            rows.append({"f": f, "p_injected": p_inj,
                         "asymmetry_mean": float(np.mean(asym)),
                         "asymmetry_std": float(np.std(asym)),
                         "twosided_mean": float(np.mean(tot)),
                         "twosided_std": float(np.std(tot))})
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda v: f"{v:11.4g}"))

    base_a = t[t["p_injected"] == 0]["asymmetry_mean"].mean()
    base_t = t[t["p_injected"] == 0]["twosided_mean"].mean()
    big = t[(t["p_injected"] == 2e-2) & (t["f"] == 0.5)].iloc[0]
    d_asym = abs(big["asymmetry_mean"] - base_a) / max(big["asymmetry_std"], 1e-9)
    d_two = abs(big["twosided_mean"] - base_t) / max(big["twosided_std"], 1e-9)
    print(f"\nresponse to p=2e-2 at f=0.5, in units of its own scatter:")
    print(f"  asymmetry estimator : {d_asym:6.2f} sigma")
    print(f"  two-sided estimator : {d_two:6.2f} sigma")
    verdict = ("ASYMMETRY ESTIMATOR IS BLIND -- limits must be rebuilt"
               if d_asym < 3 * 1.0 and d_two > 3 * d_asym else
               "asymmetry estimator responds; original limits stand")
    print(f"\nVERDICT: {verdict}")

    # --- corrected, conservative two-sided limits -----------------------
    print("\n=== corrected two-sided limits (no background subtraction) ===")
    n_stars = 2 * n_pairs
    out_rows = []
    for k in (3, 4, 5, 6):
        npos, nneg, a, tot = stats(dr0, k)
        delta = k * sigma
        f_det = float(st.fraction_from_delta(delta))
        ul = st.poisson_upper_limit(tot)
        out_rows.append({"k": k, "f_detectable": f_det,
                         "n_pos": npos, "n_neg": nneg,
                         "two_sided_total": tot, "poisson_UL": ul,
                         "p_UL": ul / n_stars,
                         "mean_f_UL": ul / n_stars * f_det})
    o = pd.DataFrame(out_rows)
    print(o.to_string(index=False, float_format=lambda v: f"{v:11.5g}"))
    best = o.loc[o["mean_f_UL"].idxmin()]
    print(f"\nCORRECTED best: p < {best['p_UL']:.3e} at f >= "
          f"{best['f_detectable']:.3f}  (mean f_bar < {best['mean_f_UL']:.3e})")

    # what the single-star test gives at the same threshold, for contrast
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["residual"])
    r = d["residual"].to_numpy(float)
    s_s = st.robust_sigma(r); m_s = float(np.median(r))
    for k in (6,):
        n_s = int(np.count_nonzero(np.abs(r - m_s) > k * s_s))
        print(f"\n  single stars at k={k}: {n_s:,} of {len(r):,} "
              f"-> rate {n_s/len(r):.2e}")
        print(f"  clean pairs  at k={k}: {stats(dr0,k)[3]:,} of {n_stars:,} "
              f"-> rate {stats(dr0,k)[3]/n_stars:.2e}")
        print(f"  pair sample background is "
              f"{(n_s/len(r))/(stats(dr0,k)[3]/n_stars):.1f}x cleaner")

    res = {"tag": args.tag, "n_pairs": int(n_pairs),
           "sigma_dr": float(sigma),
           "injection_response_asymmetry_sigma": float(d_asym),
           "injection_response_twosided_sigma": float(d_two),
           "verdict": verdict,
           "corrected_best_p_UL": float(best["p_UL"]),
           "corrected_best_f": float(best["f_detectable"]),
           "corrected_best_mean_f_UL": float(best["mean_f_UL"]),
           "table": o.to_dict(orient="records")}
    (cfg.RESULT_DIR / f"pair_estimator_fix_{args.tag}.json").write_text(
        json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
