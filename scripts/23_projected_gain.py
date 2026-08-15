#!/usr/bin/env python
"""How much does the limit improve if we cut on the contamination proxies?

    run.sh scripts/23_projected_gain.py --tag primary

scripts/18_binary_hypothesis.py showed the 5-sigma positive tail rate varies by
20x across quartiles of the BP/RP excess factor C*, and 2-3x across RUWE and
astrometric_excess_noise.  The model-free upper limit is p_UL ~ n_pos / N, so
throwing away stars in the contaminated quartiles costs N linearly but can buy
back the background rate by a much larger factor.

This scans cuts and reports the projected limit.  It is an ESTIMATE: it reuses
the existing residuals rather than refitting the fiducial inside each
subsample, so the real gain will differ (probably improve, since a cleaner
sample also fits better).

IMPORTANT CAVEAT ON WHETHER THE CUT IS FREE
-------------------------------------------
C* compares BP+RP flux to G flux.  An absorber that is GREY ACROSS THE OPTICAL
attenuates all three equally and leaves C* unchanged, so cutting on C* costs no
signal.  An absorber that is spectrally selective *within* the optical does
change C*, and a hard C* cut could then remove real signal.  The injection
framework cannot test this, because injection adds a deficit to M_G only and
by construction leaves C* alone.  Any limit derived after a C* cut is therefore
conditional on optical greyness and must be labelled as such.
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


def limit_for(resid: np.ndarray, k: float = 5.0, f: float = 0.5) -> dict:
    """Model-free 95% CL limit on p from a set of residuals."""
    r = resid[np.isfinite(resid)]
    n = len(r)
    if n < 1000:
        return {"n": n, "p_ul": np.nan, "n_pos": np.nan, "sigma": np.nan}
    s = st.robust_sigma(r)
    med = float(np.median(r))
    n_pos = int(np.count_nonzero(r > med + k * s))
    eff = st.detection_efficiency(r, f, k, sigma=s, median=med)
    n_ul = st.poisson_upper_limit(n_pos)
    p_ul = n_ul / (n * eff) if eff > 0 else np.inf
    return {"n": n, "sigma": s, "n_pos": n_pos, "pos_rate": n_pos / n,
            "efficiency": eff, "p_ul": min(p_ul, 1.0),
            "mean_f_ul": min(p_ul, 1.0) * f}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["residual", "cstar_nsigma", "ruwe",
                                 "astrometric_excess_noise", "sky_density"])
    r = d["residual"].to_numpy(float)
    base = limit_for(r, args.k)
    print(f"BASELINE  N={base['n']:,}  sigma={base['sigma']:.5f}  "
          f"n_pos={base['n_pos']:,}  rate={base['pos_rate']:.5f}")
    print(f"          p_UL={base['p_ul']:.3e}   mean_f_UL={base['mean_f_ul']:.3e}\n")

    rows = [{"cut": "none (baseline)", **base, "gain": 1.0}]

    cstar = d["cstar_nsigma"].to_numpy(float)
    ruwe = d["ruwe"].to_numpy(float)
    aen = d["astrometric_excess_noise"].to_numpy(float)

    scans = []
    for q in [0.75, 0.50, 0.25, 0.10, 0.05]:
        thr = np.nanquantile(cstar, q)
        scans.append((f"C* below {int(q*100)}th pct ({thr:+.3f})", cstar <= thr))
    for q in [0.50, 0.25, 0.10]:
        thr = np.nanquantile(ruwe, q)
        scans.append((f"RUWE below {int(q*100)}th pct ({thr:.3f})", ruwe <= thr))
    # combinations
    for qc, qr in [(0.50, 0.50), (0.25, 0.50), (0.25, 0.25), (0.10, 0.25),
                   (0.10, 0.10), (0.05, 0.25)]:
        m = (cstar <= np.nanquantile(cstar, qc)) & (ruwe <= np.nanquantile(ruwe, qr))
        scans.append((f"C*<{int(qc*100)}pct AND RUWE<{int(qr*100)}pct", m))
    # add astrometric excess noise to the best-looking combination
    m = ((cstar <= np.nanquantile(cstar, 0.10))
         & (ruwe <= np.nanquantile(ruwe, 0.25))
         & (aen <= np.nanquantile(aen, 0.50)))
    scans.append(("C*<10pct AND RUWE<25pct AND aen<50pct", m))

    for label, mask in scans:
        res = limit_for(r[mask], args.k)
        res["gain"] = base["p_ul"] / res["p_ul"] if res["p_ul"] > 0 else np.nan
        rows.append({"cut": label, **res})

    out = pd.DataFrame(rows)
    out.to_csv(cfg.RESULT_DIR / f"projected_gain_{args.tag}.csv", index=False)
    print(out[["cut", "n", "sigma", "n_pos", "pos_rate", "p_ul",
               "mean_f_ul", "gain"]].to_string(
        index=False, float_format=lambda v: f"{v:11.5g}"))

    best = out.loc[out["p_ul"].idxmin()]
    print(f"\nBEST: {best['cut']}")
    print(f"  N = {int(best['n']):,}  (from {base['n']:,}, "
          f"keeping {100*best['n']/base['n']:.0f}%)")
    print(f"  background rate {base['pos_rate']:.5f} -> {best['pos_rate']:.5f}"
          f"  ({base['pos_rate']/best['pos_rate']:.1f}x cleaner)")
    print(f"  p_UL  {base['p_ul']:.3e} -> {best['p_ul']:.3e}"
          f"   ({best['gain']:.1f}x better)")
    print(f"  mean f_UL {base['mean_f_ul']:.3e} -> {best['mean_f_ul']:.3e}")
    print("\n  CONDITIONAL on the absorber being grey across the optical "
          "(see module docstring).")

    (cfg.RESULT_DIR / f"projected_gain_{args.tag}.json").write_text(json.dumps({
        "baseline": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                     for k, v in base.items()},
        "best_cut": str(best["cut"]),
        "best_p_ul": float(best["p_ul"]),
        "best_mean_f_ul": float(best["mean_f_ul"]),
        "gain": float(best["gain"]),
        "caveat": ("conditional on optical greyness; a C* cut can remove "
                   "signal from an absorber that is selective within the "
                   "optical, and injection cannot test this"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
