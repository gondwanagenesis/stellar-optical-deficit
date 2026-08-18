#!/usr/bin/env python
"""Search J: cold blackbody radiators in the 3-100 K band nobody has searched.

    run.sh scripts/59_searchJ_cold_blackbody.py

THE COVERAGE GAP
----------------
Every published waste-heat search covers 100-1000 K: Carrigan (2009) on IRAS,
Suazo et al. (2022, 2024) on WISE, Wright et al.'s G-hat on galaxies. Below
100 K there are no searches at all -- not because the sensitivity is missing
(Planck, AKARI and IRAS reach ~1e-3 Lsun at 100 pc for temperatures down to a
few K) but because everything cold looks like Galactic cirrus and gets
catalogued as a molecular cloud core.

That gap is exactly where thermodynamics points. Erasing a bit costs kT ln2,
so computation per joule scales as 1/T: a civilisation optimising for total
computation rather than raw power is driven toward cold radiators, and would
be invisible to every search ever run.

THE DISCRIMINANT
----------------
Interstellar dust is not a blackbody. Its emissivity rises with frequency,
    kappa(nu) ~ nu^beta,   beta = 1.5-2.0 for silicate/carbonaceous grains,
because the grains are far smaller than the wavelength. That is grain physics,
not a fitting convention.

An engineered thermal radiator has no such constraint. A surface built to dump
heat is designed to be a near-ideal emitter, which means beta ~ 0.

So the search is: among cold sources, find the ones whose spectrum is a
BLACKBODY rather than a dust spectrum.

WHY THE PLANCK COLD CLUMPS
--------------------------
The PGCC (Planck 2015 results XXVIII) is all-sky, contains 13,188 sources at
6-20 K, and -- crucially -- publishes a fit with beta as a FREE parameter
(betaC) alongside the conventional fixed-beta=2 fit. So the measurement we
need already exists in a public catalogue and has simply never been read this
way. We do not refit; we select the tail.

THE CONTROL THAT MAKES IT A MEASUREMENT
---------------------------------------
beta has a physical range. Dust cannot have beta < ~1, and no known material
has beta > ~3. Both tails are therefore unphysical, and a noisy SED fit
scatters into them symmetrically. Counting the beta > 3.5 tail measures the
rate at which this catalogue manufactures unphysical emissivities from noise
alone, which is the false-positive rate for the beta < 0.5 tail we care about.
An absorber can only dim; a radiator can only be at most a blackbody -- the
asymmetry is the signal.
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchJ")

PGCC_CAT = "J/A+A/594/A28/pgcc"

BETA_BLACKBODY_MAX = 0.5     # engineered radiator: near-ideal emitter
BETA_UNPHYSICAL_MIN = 3.5    # mirror control: no material does this either
BETA_DUST_LO, BETA_DUST_HI = 1.0, 2.5

# VizieR caps the returned VOTable by total volume, not row count: asking for
# ~27 columns silently truncates at 2000 rows even with ROW_LIMIT = -1. Keep
# the column list to what the search actually needs.
WANT = ["Name", "GLON", "GLAT", "SNR",
        "q_Flux", "FBlend",
        "TempC", "e_TempC", "betaC", "e_betaC"]


def load_pgcc(force=False):
    out = cfg.RAW_DIR / "pgcc_planck2015.parquet"
    if out.exists() and not force:
        d = pd.read_parquet(out)
        log.info("cached PGCC: %d sources", len(d))
        return d

    from astroquery.vizier import Vizier
    # The row limit must be set on the CLASS before instantiation; passing
    # row_limit=-1 to the constructor is silently ignored and caps at 2000.
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=WANT)
    v.ROW_LIMIT = -1
    log.info("querying VizieR %s ...", PGCC_CAT)
    t = v.get_catalogs(PGCC_CAT)[0]
    log.info("retrieved %d rows", len(t))
    d = t.to_pandas()
    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
    d.to_parquet(out, index=False)
    log.info("wrote %s (%d sources)", out, len(d))
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    d = load_pgcc(force=args.force)
    n_all = len(d)
    log.info("PGCC sources: %d", n_all)

    # ---- require a usable free-beta fit ----------------------------------
    have = d["betaC"].notna() & d["e_betaC"].notna() & d["TempC"].notna()
    d = d[have].copy()
    log.info("with a free-beta fit: %d", len(d))

    beta = d["betaC"].to_numpy(float)
    beta_err = d["e_betaC"].to_numpy(float)
    temp = d["TempC"].to_numpy(float)

    log.info("beta distribution: median %.2f, 16-84%% [%.2f, %.2f]",
             float(np.median(beta)),
             float(np.percentile(beta, 16)), float(np.percentile(beta, 84)))
    log.info("temperature: median %.1f K, range %.1f - %.1f K",
             float(np.median(temp)), float(np.nanmin(temp)), float(np.nanmax(temp)))

    # ---- quality: the fit must actually constrain beta -------------------
    # A beta consistent with zero AND with two is not evidence of anything.
    good = (beta_err > 0) & (beta_err < 0.5)
    if "q_Flux" in d.columns:
        good &= d["q_Flux"].fillna(3).to_numpy(float) <= 2
    if "FBlend" in d.columns:
        good &= d["FBlend"].fillna(1).to_numpy(float) == 0
    good &= d["SNR"].fillna(0).to_numpy(float) > 4.0

    dq = d[good].copy()
    beta_q = beta[good]
    beta_err_q = beta_err[good]
    temp_q = temp[good]
    log.info("well-constrained beta, clean flux, unblended, SNR>4: %d", len(dq))

    if len(dq) < 50:
        log.error("too few sources survive quality cuts")
        return 1

    # ---- the two tails ----------------------------------------------------
    # Require the tail membership to be significant, not just a point estimate.
    blackbody = (beta_q < BETA_BLACKBODY_MAX) & \
                (beta_q + 3 * beta_err_q < BETA_DUST_LO)
    unphysical = (beta_q > BETA_UNPHYSICAL_MIN) & \
                 (beta_q - 3 * beta_err_q > BETA_DUST_HI)
    dustlike = (beta_q >= BETA_DUST_LO) & (beta_q <= BETA_DUST_HI)

    n_bb = int(blackbody.sum())
    n_un = int(unphysical.sum())
    n_dust = int(dustlike.sum())

    # Point-estimate counts too, so the significance requirement is visible
    # rather than silently doing all the work.
    n_bb_raw = int((beta_q < BETA_BLACKBODY_MAX).sum())
    n_un_raw = int((beta_q > BETA_UNPHYSICAL_MIN).sum())
    log.info("")
    log.info("point estimate only: beta < %.1f -> %d, beta > %.1f -> %d",
             BETA_BLACKBODY_MAX, n_bb_raw, BETA_UNPHYSICAL_MIN, n_un_raw)
    log.info("beta < %.1f  (blackbody-like, 3-sigma)  : %d",
             BETA_BLACKBODY_MAX, n_bb)
    log.info("beta > %.1f  (unphysical mirror control): %d",
             BETA_UNPHYSICAL_MIN, n_un)
    log.info("beta in [%.1f, %.1f] (dust-like)        : %d",
             BETA_DUST_LO, BETA_DUST_HI, n_dust)

    ratio = n_bb / n_un if n_un > 0 else float("inf") if n_bb > 0 else 0.0
    log.info("asymmetry ratio (blackbody : unphysical) = %s",
             f"{ratio:.2f}" if np.isfinite(ratio) else "inf")

    # ---- Poisson significance of the asymmetry ---------------------------
    # Under the null the two tails are populated at the same rate by fit
    # noise, so n_bb ~ Poisson(n_un). Use the exact two-sided binomial split.
    from scipy import stats
    if n_bb + n_un > 0:
        p_binom = stats.binomtest(n_bb, n_bb + n_un, 0.5,
                                  alternative="greater").pvalue
    else:
        p_binom = 1.0
    log.info("P(n_blackbody >= observed | symmetric noise) = %.3g", p_binom)

    # ---- where do the blackbody-like sources sit? ------------------------
    dq = dq.assign(beta=beta_q, beta_err=beta_err_q, temp_K=temp_q,
                   tail=np.where(blackbody, "blackbody",
                                 np.where(unphysical, "unphysical",
                                          np.where(dustlike, "dust", "other"))))

    if n_bb > 0:
        bb = dq[dq["tail"] == "blackbody"]
        med_absb_bb = float(np.abs(bb["GLAT"]).median())
        med_absb_dust = float(np.abs(dq[dq["tail"] == "dust"]["GLAT"]).median())
        med_T_bb = float(bb["temp_K"].median())
        med_T_dust = float(dq[dq["tail"] == "dust"]["temp_K"].median())

        log.info("")
        log.info("blackbody-tail sources: median |b| = %.1f deg "
                 "(dust-like: %.1f)", med_absb_bb, med_absb_dust)
        log.info("blackbody-tail sources: median T = %.1f K "
                 "(dust-like: %.1f)", med_T_bb, med_T_dust)

        # A genuine engineered radiator should NOT be preferentially in the
        # Galactic plane, where cirrus confusion and beta degeneracy live.
        high_lat = bb[np.abs(bb["GLAT"]) > 15.0]
        log.info("blackbody-tail at |b| > 15 deg: %d", len(high_lat))
    else:
        med_absb_bb = med_absb_dust = med_T_bb = med_T_dust = float("nan")
        high_lat = dq.iloc[0:0]

    # ---- verdict ----------------------------------------------------------
    if n_bb == 0:
        verdict = (
            "NULL. No Planck cold clump has an emissivity index significantly "
            "below the dust floor. Every cold source in the all-sky catalogue "
            "is spectrally consistent with grains, not with an engineered "
            "blackbody radiator.")
    elif p_binom > 0.05:
        verdict = (
            f"NULL. The {n_bb} blackbody-like sources are matched by {n_un} in "
            f"the equally-unphysical high-beta tail (p = {p_binom:.2g} against "
            f"symmetric fit noise). The low-beta tail is SED-fitting scatter, "
            f"not a distinct population.")
    else:
        verdict = (
            f"ASYMMETRIC: {n_bb} sources with beta below the dust floor "
            f"against {n_un} in the mirror tail (p = {p_binom:.2g}). "
            f"Before this means anything, note that low beta is degenerate "
            f"with temperature in a two-band fit, is known to be produced by "
            f"line-of-sight temperature mixing along the beam, and that "
            f"{len(high_lat)} of these sit at |b| > 15 deg where cirrus "
            f"confusion is lowest. Those are the ones worth following up.")

    print(f"\n{'='*68}")
    print("SEARCH J: COLD BLACKBODY RADIATORS (3-100 K)")
    print(f"{'='*68}")
    print(f"  PGCC sources                        : {n_all:,}")
    print(f"  with a well-constrained free beta   : {len(dq):,}")
    print(f"  temperature range                   : "
          f"{float(np.nanmin(temp_q)):.1f} - {float(np.nanmax(temp_q)):.1f} K")
    print()
    print(f"  beta < {BETA_BLACKBODY_MAX} (blackbody-like)          : {n_bb}")
    print(f"  beta > {BETA_UNPHYSICAL_MIN} (unphysical, mirror)     : {n_un}")
    print(f"  beta in [{BETA_DUST_LO}, {BETA_DUST_HI}] (dust)            : {n_dust}")
    print(f"  asymmetry                           : "
          f"{ratio:.2f} : 1" if np.isfinite(ratio) else "  asymmetry: inf")
    print(f"  p(symmetric noise)                  : {p_binom:.3g}")
    if n_bb:
        print(f"  blackbody tail median |b|           : {med_absb_bb:.1f} deg "
              f"(dust {med_absb_dust:.1f})")
        print(f"  blackbody tail median T            : {med_T_bb:.1f} K "
              f"(dust {med_T_dust:.1f})")
        print(f"  blackbody tail at |b| > 15 deg     : {len(high_lat)}")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "catalogue": PGCC_CAT,
        "n_pgcc": n_all,
        "n_with_free_beta": int(len(dq)),
        "beta_thresholds": {
            "blackbody_max": BETA_BLACKBODY_MAX,
            "unphysical_min": BETA_UNPHYSICAL_MIN,
            "dust_range": [BETA_DUST_LO, BETA_DUST_HI],
        },
        "n_blackbody_tail": n_bb,
        "n_unphysical_mirror": n_un,
        "n_blackbody_point_estimate": n_bb_raw,
        "n_unphysical_point_estimate": n_un_raw,
        "n_dust_like": n_dust,
        "asymmetry_ratio": float(ratio) if np.isfinite(ratio) else None,
        "p_symmetric_noise": float(p_binom),
        "median_absb_blackbody": med_absb_bb,
        "median_absb_dust": med_absb_dust,
        "median_temp_blackbody": med_T_bb,
        "median_temp_dust": med_T_dust,
        "n_blackbody_high_latitude": int(len(high_lat)),
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchJ_cold_blackbody_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)

    flagged = dq[dq["tail"].isin(["blackbody", "unphysical"])]
    if len(flagged):
        csv = cfg.RESULT_DIR / f"searchJ_candidates_{args.tag}.csv"
        flagged.sort_values("beta").to_csv(csv, index=False)
        log.info("wrote %d tail sources to %s", len(flagged), csv)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
