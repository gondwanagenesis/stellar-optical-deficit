#!/usr/bin/env python
"""Search E: dark stellar-mass companions that violate the white-dwarf channel.

    run.sh scripts/56_searchE_dark_companions.py

THE HYPOTHESIS
--------------
Every photometric channel in this project measures a star's light. A star that
is COMPLETELY enshrouded (f -> 1) does not produce a small deficit -- it leaves
the photometric catalogue altogether, which is why channels 1-12 are blind to
exactly the case they were built to find.

Gravity is not. In an astrometric binary the companion's mass is measured from
the photocentre wobble whether or not the companion emits anything at all.
A fully enshrouded star is therefore a MAXIMAL astrometric signal.

WHAT IS ALREADY KNOWN, SO THAT THIS IS NOT A REDISCOVERY
--------------------------------------------------------
Shahaf et al. (2019, 2023, 2024) built exactly this machinery -- the
astrometric mass-ratio function -- and published the classification of 101,380
Gaia DR3 orbits into main-sequence companions (class I), hierarchical triples
(class II) and compact objects (class III), plus ~3,145 white-dwarf
candidates. "There is a dark companion here" is their result, not ours.

The AMRF measures mass and darkness. It cannot say WHY something is dark, so
every class-III system is already "a companion too dark for its dynamical
mass". Our contribution has to be a discriminant applied ON TOP of that, and
the two obvious ones are spent: Shahaf used GALEX ultraviolet excess to
confirm hot white dwarfs and optical colour excess to split off M-dwarf pairs.

THE DISCRIMINANT WE USE INSTEAD: THE MASS-PERIOD RELATION
---------------------------------------------------------
Below ~0.50 Msun no white dwarf can form from a single star in a Hubble time
(Cummings+2018: the lowest cluster-calibrated white dwarf is 0.507 Msun from a
0.85 Msun progenitor in a 12 Gyr cluster). Anything lighter is a HELIUM white
dwarf and must have been stripped by a companion.

That stripping welds mass to orbit. The donor's degenerate core mass sets both
the remnant mass and the donor's radius at the end of Roche-lobe overflow, so

    P_orb = b * (M_WD - c)^a        (Tauris & Savonije 1999, Pop I:
                                     a = 4.50, b = 1.2e5 d, c = 0.120 Msun)

is a one-parameter physics-locked relation, and it TERMINATES near 0.47 Msun
where helium ignites. A 0.25 Msun helium white dwarf belongs at ~10 days. Find
one at 2000 days and it is off the relation by two orders of magnitude.

Stable Roche-lobe overflow also circularises. Phinney's (1992) fluctuation-
dissipation floor predicts residual e ~ 1e-3 at P ~ 1000 d, so an eccentric
wide orbit is a second, independent violation.

WHAT WE DO NOT CLAIM
--------------------
KIC 8145411 is a real, natural 0.20 Msun white dwarf at P = 450 d, e = 0.14 --
sitting inside the box this search calls anomalous. It took two dedicated
theory papers to explain. So the honest claim is never "forbidden", it is
"requires a channel whose Galactic rate is ~1e3-1e4 per Gyr": dynamical
exchange in a cluster, or a supernova in a hierarchical triple. Both of those
peak at P ~ 450-770 d, which is why the genuinely clean region is P > 1500 d.

AND THE STRUCTURAL PROBLEM, STATED UP FRONT
-------------------------------------------
DR3's astrometric baseline is ~1000 days. The physically clean region begins
at 1500 days. The region where this search is most powerful is therefore
exactly the region where the orbits are extrapolated and least trustworthy.
That tension is not resolvable with DR3 and is reported rather than hidden;
Gaia DR4 (66-month baseline, epoch astrometry) is the release that fixes it.
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
log = logging.getLogger("searchE")

# Tauris & Savonije (1999) He-WD mass-period relation, Population I
TS99_A, TS99_B, TS99_C = 4.50, 1.2e5, 0.120
HE_WD_MAX_MASS = 0.47      # helium ignites above this -> CO WD instead
CO_WD_MIN_MASS = 0.50      # Cummings+2018 single-star floor
CHANDRASEKHAR = 1.40

DAYS_PER_YEAR = 365.25
RNG_SEED = 56_2026


# ---------------------------------------------------------------------------
# AMRF machinery
# ---------------------------------------------------------------------------

def photocentre_semimajor(a_ti, b_ti, f_ti, g_ti):
    """Photocentre semi-major axis a0 (mas) from the Thiele-Innes constants.

    Gaia DR3 documentation, Halbwachs et al. 2023 (A&A 674, A9):
        u = (A^2 + B^2 + F^2 + G^2) / 2
        v = A*G - B*F
        a0 = sqrt(u + sqrt(u^2 - v^2))
    """
    u = (a_ti ** 2 + b_ti ** 2 + f_ti ** 2 + g_ti ** 2) / 2.0
    v = a_ti * g_ti - b_ti * f_ti
    disc = np.maximum(u ** 2 - v ** 2, 0.0)
    return np.sqrt(np.maximum(u + np.sqrt(disc), 0.0))


def amrf(a0_mas, parallax_mas, m1, period_days):
    """Astrometric mass-ratio function (Shahaf et al. 2019, eq. 1)."""
    p_yr = period_days / DAYS_PER_YEAR
    return (a0_mas / parallax_mas) * m1 ** (-1.0 / 3.0) * p_yr ** (-2.0 / 3.0)


def q_from_amrf_dark(script_a, n_iter=80):
    """Mass ratio q for a COMPLETELY DARK companion, i.e. flux ratio S = 0.

    Then A = q / (1+q)^(2/3), which is monotonic in q, so bisect. Any light
    from the companion reduces the photocentre wobble, so the dark solution is
    a strict LOWER bound on q -- and therefore on the companion mass.
    """
    lo = np.zeros_like(script_a)
    hi = np.full_like(script_a, 50.0)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        val = mid / (1.0 + mid) ** (2.0 / 3.0)
        too_small = val < script_a
        lo = np.where(too_small, mid, lo)
        hi = np.where(too_small, hi, mid)
    return 0.5 * (lo + hi)


def ts99_period(m_wd):
    """Period a helium white dwarf of this mass MUST have (days)."""
    return TS99_B * np.maximum(m_wd - TS99_C, 1e-6) ** TS99_A


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--max-dist-pc", type=float, default=300.0)
    args = ap.parse_args()

    src = cfg.RAW_DIR / "nss_orbits_500pc.parquet"
    if not src.exists():
        log.error("missing %s -- run scripts/52_pull_nss_orbits.py first", src)
        return 1

    d = pd.read_parquet(src)
    log.info("NSS orbital solutions within 500 pc: %d", len(d))

    n0 = len(d)
    cut_log = [("pulled", n0)]

    # ---- primary mass is required ----------------------------------------
    d = d[d["m1"].notna() & (d["m1"] > 0.1) & (d["m1"] < 3.0)]
    cut_log.append(("has Gaia primary mass m1", len(d)))

    # Shahaf restricts to the combinations where m1 comes from an isochrone
    # plus luminosity, i.e. a genuine main-sequence primary.
    if "combination_method" in d.columns:
        d = d[d["combination_method"].isin(["Orbital+M1", "AstroSpectroSB1+M1",
                                            "Orbital+SB1+M1"])]
        cut_log.append(("m1 from isochrone+luminosity", len(d)))

    # ---- Thiele-Innes quality (Shahaf et al. 2023, eq. 4) -----------------
    rel = np.zeros(len(d))
    for c in ("a_thiele_innes", "b_thiele_innes",
              "f_thiele_innes", "g_thiele_innes"):
        val = d[c].to_numpy(float)
        err = d[f"{c}_error"].to_numpy(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel += np.where(np.abs(val) > 0, (err / val) ** 2, np.inf)
    d = d[np.isfinite(rel) & (rel <= 36.0)]
    cut_log.append(("Thiele-Innes relative error", len(d)))

    # ---- DR3 pipeline validity cuts (Halbwachs+2023) ----------------------
    P = d["period"].to_numpy(float)
    plx = d["parallax"].to_numpy(float)
    plx_err = d["parallax_error"].to_numpy(float)
    ecc_err = d["eccentricity_error"].to_numpy(float)

    a0 = photocentre_semimajor(
        d["a_thiele_innes"].to_numpy(float), d["b_thiele_innes"].to_numpy(float),
        d["f_thiele_innes"].to_numpy(float), d["g_thiele_innes"].to_numpy(float))
    sig = d["significance"].to_numpy(float)          # = a0 / sigma_a0

    keep = (
        (ecc_err < 0.079 * np.log(np.maximum(P, 1.1)) - 0.244)
        & (plx / plx_err > 20000.0 / np.maximum(P, 1.0))
        & (sig > 158.0 / np.sqrt(np.maximum(P, 1.0)))
    )
    d, a0, P, sig = d[keep], a0[keep], P[keep], sig[keep]
    cut_log.append(("DR3 pipeline validity", len(d)))

    # ---- magnitude-conditioned goodness of fit ---------------------------
    # F2 is not N(0,1) in practice: for G < 13 it peaks near 3.5 and 13% of
    # solutions exceed 10, so a flat cut would preferentially delete bright
    # stars (El-Badry et al. 2023, appendix E).
    gof = d["goodness_of_fit"].to_numpy(float)
    gmag = d["phot_g_mean_mag"].to_numpy(float)
    gof_ok = np.where(gmag > 13.0, gof < 5.0, gof < 15.0)
    d, a0, P, sig = d[gof_ok], a0[gof_ok], P[gof_ok], sig[gof_ok]
    cut_log.append(("magnitude-conditioned F2", len(d)))

    # ---- scanning-law period rejection -----------------------------------
    # Spurious solutions pile up at the satellite's 63 d precession period,
    # at 1 yr, and at their aliases (Holl et al. 2023).
    bad = np.zeros(len(d), dtype=bool)
    for p_bad in (63.0, 182.6, 365.25, 730.5):
        bad |= np.abs(P - p_bad) / p_bad < 0.15
    d, a0, P, sig = d[~bad], a0[~bad], P[~bad], sig[~bad]
    cut_log.append(("scanning-law periods removed", len(d)))

    # ---- distance ---------------------------------------------------------
    plx = d["parallax"].to_numpy(float)
    dist_pc = 1000.0 / plx
    near = dist_pc < args.max_dist_pc
    d, a0, P, sig, dist_pc = d[near], a0[near], P[near], sig[near], dist_pc[near]
    cut_log.append((f"within {args.max_dist_pc:.0f} pc", len(d)))

    for label, n in cut_log:
        log.info("  %-34s %7d", label, n)

    if len(d) < 100:
        log.error("too few systems survive; aborting")
        return 1

    # ---- the estimator ----------------------------------------------------
    m1 = d["m1"].to_numpy(float)
    plx = d["parallax"].to_numpy(float)
    ecc = d["eccentricity"].to_numpy(float)

    script_a = amrf(a0, plx, m1, P)
    q_dark = q_from_amrf_dark(script_a)
    m2_min = q_dark * m1

    d = d.assign(a0_mas=a0, amrf=script_a, q_dark=q_dark,
                 m2_min=m2_min, dist_pc=dist_pc)

    log.info("AMRF: median %.4f, 90th pct %.4f",
             float(np.median(script_a)), float(np.percentile(script_a, 90)))
    log.info("minimum companion mass: median %.3f Msun", float(np.median(m2_min)))

    # ---- the mass-period violation ---------------------------------------
    p_required = ts99_period(m2_min)
    violation = P / np.maximum(p_required, 1e-6)
    d = d.assign(p_ts99_required=p_required, mp_violation=violation)

    sub_co = m2_min < CO_WD_MIN_MASS          # cannot be a single-star WD
    he_range = (m2_min >= 0.20) & (m2_min <= 0.45)

    log.info("companion below the single-star WD floor (%.2f Msun): %d",
             CO_WD_MIN_MASS, int(sub_co.sum()))
    log.info("in the helium-WD mass window 0.20-0.45 Msun: %d",
             int(he_range.sum()))

    # ---- the anomaly box --------------------------------------------------
    # Every natural channel is excluded by at least two of these at once.
    anomalous = he_range & (P > 1500.0) & (ecc > 0.15)
    marginal = he_range & (P > 800.0) & (ecc > 0.15) & ~anomalous

    log.info("ANOMALY BOX (0.20-0.45 Msun, P > 1500 d, e > 0.15): %d",
             int(anomalous.sum()))
    log.info("marginal box    (0.20-0.45 Msun, P >  800 d, e > 0.15): %d",
             int(marginal.sum()))

    # ---- the control that makes this a measurement ------------------------
    # The AMRF is a magnitude, so it has no sign to mirror. The honest control
    # is the eccentricity axis: stable Roche-lobe overflow circularises, so
    # genuine helium white dwarfs must sit at LOW eccentricity. Counting the
    # same mass-and-period box at low e gives the rate at which this pipeline
    # populates the box for reasons that have nothing to do with the signal.
    circular_box = he_range & (P > 1500.0) & (ecc < 0.05)
    log.info("circular control (same M2 and P, e < 0.05): %d",
             int(circular_box.sum()))

    # Distance-eccentricity systematic: Shahaf documents fitted eccentricity
    # rising with distance, which would manufacture exactly this signal.
    if int(anomalous.sum()) > 0:
        med_d_anom = float(d.loc[anomalous, "dist_pc"].median())
    else:
        med_d_anom = float("nan")
    med_d_all = float(d["dist_pc"].median())
    ecc_dist_corr = float(pd.Series(ecc).corr(pd.Series(dist_pc)))
    log.info("eccentricity-distance correlation: r = %+.3f "
             "(positive => fitted-e inflation with distance)", ecc_dist_corr)

    # ---- the baseline caveat, quantified ---------------------------------
    beyond_baseline = P > 1000.0
    log.info("periods beyond the DR3 ~1000 d baseline (extrapolated): "
             "%d (%.1f%%)", int(beyond_baseline.sum()),
             100 * beyond_baseline.mean())

    # ---- verdict ----------------------------------------------------------
    n_anom = int(anomalous.sum())
    n_ctrl = int(circular_box.sum())

    if n_anom == 0:
        verdict = (
            "NULL. No companion in the helium-white-dwarf mass window sits at "
            "a period and eccentricity that the stable-Roche-lobe channel "
            "cannot produce. Every dark companion recovered here is consistent "
            "with a white dwarf, a hierarchical triple, or a compact object -- "
            "the populations Shahaf et al. already catalogued.")
    elif n_anom <= n_ctrl:
        verdict = (
            f"NULL. The {n_anom} systems in the anomaly box are matched by "
            f"{n_ctrl} in the circular control box, so the eccentric excess is "
            f"not significant: the box is being populated by orbit-fitting "
            f"scatter rather than by a distinct population.")
    else:
        verdict = (
            f"{n_anom} systems in the anomaly box against {n_ctrl} in the "
            f"circular control. These violate the helium-white-dwarf "
            f"mass-period relation by a median factor of "
            f"{float(d.loc[anomalous, 'mp_violation'].median()):.0f}. "
            f"BEFORE this means anything: {int((d.loc[anomalous, 'period'] > 1000).sum())} "
            f"have periods beyond the DR3 baseline and are extrapolated; each "
            f"needs radial-velocity confirmation, a tertiary search to exclude "
            f"Kozai-Lidov pumping, and a colour-excess check to exclude an "
            f"M-dwarf inner pair.")

    print(f"\n{'='*66}")
    print("SEARCH E: DARK COMPANIONS VS THE WHITE-DWARF CHANNEL")
    print(f"{'='*66}")
    print(f"  clean astrometric orbits within {args.max_dist_pc:.0f} pc : {len(d):,}")
    print(f"  companion below single-star WD floor        : {int(sub_co.sum()):,}")
    print(f"  in helium-WD mass window (0.20-0.45)        : {int(he_range.sum()):,}")
    print(f"  ANOMALY BOX (P>1500 d, e>0.15)              : {n_anom}")
    print(f"  circular control (P>1500 d, e<0.05)         : {n_ctrl}")
    print(f"  marginal (P>800 d, e>0.15)                  : {int(marginal.sum())}")
    print(f"  eccentricity-distance correlation           : r = {ecc_dist_corr:+.3f}")
    print(f"  periods beyond DR3 baseline                 : "
          f"{int(beyond_baseline.sum()):,} ({100*beyond_baseline.mean():.1f}%)")
    print(f"\nVERDICT: {verdict}")

    summary = {
        "tag": args.tag,
        "max_dist_pc": args.max_dist_pc,
        "cut_flow": [{"cut": c, "n": int(n)} for c, n in cut_log],
        "n_clean": int(len(d)),
        "n_below_co_wd_floor": int(sub_co.sum()),
        "n_he_wd_mass_window": int(he_range.sum()),
        "n_anomaly_box": n_anom,
        "n_circular_control": n_ctrl,
        "n_marginal": int(marginal.sum()),
        "median_amrf": float(np.median(script_a)),
        "median_m2_min": float(np.median(m2_min)),
        "ecc_distance_correlation": ecc_dist_corr,
        "median_dist_anomalous_pc": med_d_anom,
        "median_dist_all_pc": med_d_all,
        "n_beyond_dr3_baseline": int(beyond_baseline.sum()),
        "ts99_params": {"a": TS99_A, "b": TS99_B, "c": TS99_C},
        "verdict": verdict,
    }
    out = cfg.RESULT_DIR / f"searchE_dark_companions_{args.tag}.json"
    out.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", out)

    keep_cols = [c for c in [
        "source_id", "nss_solution_type", "ra", "dec", "l", "b", "dist_pc",
        "period", "period_error", "eccentricity", "eccentricity_error",
        "a0_mas", "significance", "goodness_of_fit",
        "m1", "m2", "m2_min", "q_dark", "amrf",
        "p_ts99_required", "mp_violation",
        "phot_g_mean_mag", "bp_rp", "ruwe", "teff_gspphot",
        "combination_method"] if c in d.columns]

    flagged = d[anomalous | marginal].copy()
    flagged["box"] = np.where(anomalous[anomalous | marginal],
                              "anomalous", "marginal")
    if len(flagged):
        csv = cfg.RESULT_DIR / f"searchE_candidates_{args.tag}.csv"
        flagged[keep_cols + ["box"]].sort_values(
            "mp_violation", ascending=False).to_csv(csv, index=False)
        log.info("wrote %d flagged systems to %s", len(flagged), csv)

    allcsv = cfg.RESULT_DIR / f"searchE_all_amrf_{args.tag}.csv"
    d[keep_cols].to_csv(allcsv, index=False)
    log.info("wrote full AMRF table to %s", allcsv)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
