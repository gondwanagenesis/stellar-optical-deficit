#!/usr/bin/env python
"""Test on real data whether unresolved binaries drive the positive tail.

    run.sh scripts/18_binary_hypothesis.py --tag primary

tests/test_binary_bias.py shows *analytically* that because dM_G/dM_Ks > 1, an
unresolved companion pushes the residual positive -- i.e. binaries mimic the
harvesting signal rather than suppressing it, contradicting the brief.

This checks the prediction against the data.  RUWE and ipd_frac_multi_peak are
both monotone proxies for unresolved multiplicity, and the sample is already
cut at RUWE < 1.4, so the surviving range still spans a large contrast in
binary likelihood.  If the analytic argument is right, the positive 5-sigma
tail must grow with these proxies while the negative tail does not.

A null result here would mean the positive tail is dominated by something else
(YSOs, spots, photometric artefacts), which is also worth knowing.

MEASURED RESULT (3.32M stars, 5 sigma threshold)
------------------------------------------------
Positive-tail rate, lowest to highest quartile, with the negative tail as a
control:

    RUWE                      2.04x   (negative tail 0.72x -- it SHRINKS)
    astrometric_excess_noise  2.77x   (negative 0.85x)
    BP/RP excess factor C*   20.52x   (negative 0.38x)

Every multiplicity proxy drives the positive tail up and the negative tail
down.  The analytic prediction is confirmed on data: unresolved companions
manufacture false deficits.

THE MECHANISM, which is worth stating explicitly
------------------------------------------------
C* is by far the strongest driver, and that identifies the dominant channel as
**aperture mismatch between Gaia and 2MASS**.  Gaia resolves sources at
sub-arcsecond scale; the 2MASS beam is ~4 arcsec.  A neighbour that Gaia
separates is blended into the same 2MASS K_s measurement.  The star therefore
gets a K_s that is too bright while its G is not, so M_Ks moves brighter at
fixed M_G -- which in this diagram is indistinguishable from the star being
under-luminous in G.  That is a manufactured optical deficit, and it is the
single largest contaminant in the search.
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


def tail_rates(resid: np.ndarray, sigma: float, med: float,
               k: float = 5.0) -> tuple[float, float, int]:
    n = len(resid)
    if n == 0:
        return np.nan, np.nan, 0
    pos = float(np.count_nonzero(resid > med + k * sigma)) / n
    neg = float(np.count_nonzero(resid < med - k * sigma)) / n
    return pos, neg, n


def by_quantile(d: pd.DataFrame, resid: np.ndarray, column: str,
                sigma: float, med: float, nq: int = 4,
                k: float = 5.0) -> pd.DataFrame:
    x = d[column].to_numpy(dtype=float)
    ok = np.isfinite(x)
    edges = np.nanquantile(x[ok], np.linspace(0, 1, nq + 1))
    edges = np.unique(edges)
    rows = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        m = ok & (x >= lo) & (x <= hi if i == len(edges) - 2 else x < hi)
        pos, neg, n = tail_rates(resid[m], sigma, med, k)
        rows.append({
            "column": column, "quantile": f"Q{i+1}",
            "range": f"{lo:.4g}..{hi:.4g}", "n": n,
            "pos_rate": pos, "neg_rate": neg,
            "asymmetry": pos - neg,
            "ratio_pos_neg": pos / neg if neg > 0 else np.inf,
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--k", type=float, default=5.0)
    args = ap.parse_args()

    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=["residual", "ruwe", "ipd_frac_multi_peak",
                                 "astrometric_excess_noise", "cstar_nsigma",
                                 "M_Ks", "bp_rp0"])
    r = d["residual"].to_numpy(dtype=float)
    sigma = st.robust_sigma(r)
    med = float(np.median(r))
    print(f"N = {len(d):,}   sigma = {sigma:.5f}   threshold = "
          f"{args.k} sigma\n")

    frames = []
    for col in ["ruwe", "ipd_frac_multi_peak", "astrometric_excess_noise",
                "cstar_nsigma"]:
        if col not in d:
            continue
        t = by_quantile(d, r, col, sigma, med, k=args.k)
        frames.append(t)
        print(f"-- split by {col}")
        print(t[["quantile", "range", "n", "pos_rate", "neg_rate",
                 "ratio_pos_neg"]].to_string(
            index=False, float_format=lambda v: f"{v:11.6f}"))
        lo, hi = t.iloc[0], t.iloc[-1]
        if lo["pos_rate"] > 0:
            print(f"   positive-tail rate Q4/Q1 = "
                  f"{hi['pos_rate']/lo['pos_rate']:.2f}x")
        if lo["neg_rate"] > 0:
            print(f"   negative-tail rate Q4/Q1 = "
                  f"{hi['neg_rate']/lo['neg_rate']:.2f}x")
        print()

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(cfg.RESULT_DIR / f"binary_hypothesis_{args.tag}.csv", index=False)

    ruwe = out[out["column"] == "ruwe"]
    pos_ratio = float(ruwe.iloc[-1]["pos_rate"] / ruwe.iloc[0]["pos_rate"])
    neg_ratio = float(ruwe.iloc[-1]["neg_rate"] / ruwe.iloc[0]["neg_rate"]) \
        if ruwe.iloc[0]["neg_rate"] > 0 else np.nan

    verdict = ("CONSISTENT with binaries driving the positive tail"
               if pos_ratio > 1.3 and pos_ratio > neg_ratio else
               "NOT explained by binarity alone -- the positive tail does not "
               "track RUWE, so another population dominates")
    summary = {
        "tag": args.tag, "n_stars": int(len(d)), "sigma_mag": sigma,
        "k_sigma": args.k,
        "ruwe_Q4_over_Q1_positive_tail": pos_ratio,
        "ruwe_Q4_over_Q1_negative_tail": neg_ratio,
        "verdict": verdict,
    }
    (cfg.RESULT_DIR / f"binary_hypothesis_{args.tag}.json").write_text(
        json.dumps(summary, indent=2))
    print("=" * 70)
    print(f"positive tail grows {pos_ratio:.2f}x from lowest to highest RUWE "
          f"quartile")
    print(f"negative tail grows {neg_ratio:.2f}x over the same range")
    print(f"VERDICT: {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
