#!/usr/bin/env python
"""Audit every channel for the failure modes that have already bitten us.

    run.sh scripts/79_audit_all_channels.py

WHY
---
Three self-inflicted errors have been found in this project so far, and all
three produced a plausible-looking number rather than an obvious crash:

  a cross-match that ignored proper motion   -> a clean-looking FALSE NULL
    (Search O returned zero rankable stars because 2MASS observed in 1999 and
     Gaia in 2016, and Barnard's Star moves 165 arcsec in that gap)

  a permutation null over a constrained quantity -> a FALSE DEFICIT
    (Search F shuffled transverse vectors between different lines of sight,
     comparing real 2D neighbourhoods against synthetic 3D ones)

  a fit running to its grid boundary          -> a near-miss FALSE POSITIVE
    (Search T returned 1023 candidates all sitting exactly on the beta grid
     edge, with a mirror control that could not fire by construction)

None of those would have been caught by reading the headline number. So this
audits the stored results of every channel against the classes of failure that
have actually occurred, plus the one dimension never applied uniformly:

STATISTICAL POWER. A null is only informative if the test could have seen the
thing. Several channels report "no candidates" on samples small enough that
the strongest possible signal would also have produced no candidates. Those
are not nulls, they are non-measurements, and they should not be quoted
alongside the channels that genuinely constrain something.

WHAT IT CHECKS
--------------
  mirror validity   a control that reports zero may be a clean null OR a
                    control that cannot fire. The two are distinguished by
                    whether the mirror had any objects to draw from.

  power             for a counting channel, the Poisson 95% upper limit on a
                    rate given zero detections is 3/N. If 3/N is larger than
                    the effect being sought, the channel is uninformative.

  marginal nulls    p between 0.01 and 0.1 reported as if decisive.

  edge effects      any recorded best-fit parameter sitting on a stated bound.

  sample provenance whether the channel measured its own quantity or inherited
                    a literature value (Search H's stellar term).
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("audit")

# Channel -> (result file, the count fields that define signal and mirror)
CHANNELS = {
    "5  spectral slope":      ("searchA_null_primary.json",
                               "n_dimmed", "n_brightened"),
    "9  W3/W4 re-emit":       ("searchB_cold_primary.json",
                               "n_intercept_reemit_any", "n_mirror_control"),
    "10 radio":               ("searchC_radio_primary.json", None, None),
    "11 domain edge":         ("searchD_domain_edge_primary.json", None, None),
    "13 coherent accel":      ("searchF_corrected_null.json", None, None),
    "14 mass ledger":         ("searchH_mass_budget_primary.json", None, None),
    "15 cold blackbody PGCC": ("searchJ_cold_blackbody_primary.json",
                               "n_blackbody_tail", "n_unphysical_mirror"),
    "16 dark companions":     ("searchE_dark_companions_primary.json",
                               "n_anomaly_box", "n_circular_control"),
    "17 mass function":       ("searchK_mass_function_primary_blt20.json",
                               None, None),
    "18 bolometric closure":  ("searchL_bolometric_primary.json",
                               "n_deficient", "n_excessive"),
    "19 discard pile":        ("searchM_discard_pile_primary.json", None, None),
    "20 nearest 10 pc":       ("searchO_nearest10pc_primary.json",
                               "n_deficit_3sig", "n_excess_3sig"),
    "21 quiescence":          ("searchP_quiescence_primary.json", None, None),
    "22 SED vs spectrum":     ("searchQ_sed_inconsistency_primary.json",
                               "n_veiled", "n_inverse_mirror"),
    "23 joint distribution":  ("searchR_joint_anomaly_primary.json", None, None),
    "24 energy balance":      ("searchU_energy_balance_primary.json",
                               "n_balanced", "n_mirror_balanced"),
    "25 Planck excluded bin": ("searchT_pccs2e_primary.json",
                               "n_cold_blackbody", "n_unphysical_mirror"),
}

# Sample size actually tested, for the power calculation
N_TESTED = {
    "5  spectral slope": 2_991_398,
    "9  W3/W4 re-emit": 104_299,
    "10 radio": 5_000,
    "13 coherent accel": 31_989,
    "15 cold blackbody PGCC": 3_030,
    "16 dark companions": 17_717,
    "17 mass function": 11_225,
    "18 bolometric closure": 94_123,
    "20 nearest 10 pc": 91,
    "22 SED vs spectrum": 32_016,
    "23 joint distribution": 2_565_377,
    "24 energy balance": 2_632,
    "25 Planck excluded bin": 20_932,
}


def main() -> int:
    findings = []
    log.info("=" * 78)
    log.info("CHANNEL AUDIT")
    log.info("=" * 78)

    for name, (fname, sig_k, mir_k) in CHANNELS.items():
        p = cfg.RESULT_DIR / fname
        if not p.exists():
            log.info("\n%-24s  MISSING RESULT FILE (%s)", name, fname)
            findings.append({"channel": name, "issue": "missing result file",
                             "severity": "blocking"})
            continue
        d = json.loads(p.read_text())
        log.info("\n%-24s  %s", name, fname)

        n = N_TESTED.get(name)
        issues = []

        # ---- mirror validity ---------------------------------------------
        if sig_k and mir_k and sig_k in d and mir_k in d:
            s, m = d[sig_k], d[mir_k]
            log.info("   signal %-8s = %-8s   mirror %-8s = %s",
                     sig_k, s, mir_k, m)
            if m == 0 and s == 0:
                issues.append("BOTH tails empty: the control never had "
                              "anything to draw from, so this is a "
                              "non-measurement rather than a null")
            elif m == 0 and s > 0:
                issues.append("mirror is ZERO while signal is not: verify the "
                              "control can fire at all before reading the "
                              "asymmetry (this is exactly how Search T's "
                              "grid-edge false positive presented)")
            elif m > 0:
                ratio = s / m
                log.info("   signal:mirror = %.2f", ratio)
                if ratio <= 1.0:
                    log.info("   -> null, and conservative: the mirror "
                             "exceeds the signal")

        # ---- power -------------------------------------------------------
        if n:
            ul = 3.0 / n           # Poisson 95% UL on a rate given 0 events
            log.info("   n tested = %-10s  ->  95%% rate UL if zero found = %.2e",
                     f"{n:,}", ul)
            if ul > 1e-2:
                issues.append(f"LOW POWER: with n={n:,} the best possible "
                              f"limit is a rate of {ul:.1e}, so this channel "
                              f"cannot constrain anything rarer than about "
                              f"1 in {int(1/ul):,}")

        # ---- marginal significance ---------------------------------------
        for key in ("p_empirical", "p_symmetric_noise", "p_vs_reference"):
            if key in d and d[key] is not None:
                pv = float(d[key])
                log.info("   %s = %.4f", key, pv)
                if 0.01 < pv < 0.10:
                    issues.append(f"MARGINAL: {key} = {pv:.3f} is neither a "
                                  f"detection nor a clean null and should not "
                                  f"be quoted as decisive")

        # ---- explicit self-reported caveats -------------------------------
        v = str(d.get("verdict", ""))
        for phrase, note in (
                ("too few", "the test reported insufficient objects"),
                ("could not", "the test reported an inability to proceed"),
                ("unreliable", "the test flagged its own inputs as unreliable"),
                ("grid", "grid/boundary effects mentioned in the verdict")):
            if phrase in v.lower():
                issues.append(f"verdict self-flags: {note}")

        if issues:
            for i in issues:
                log.info("   [!] %s", i)
                findings.append({"channel": name, "issue": i})
        else:
            log.info("   no flags")

    # ---- channels whose result is known to rest on inherited values ------
    log.info("\n" + "=" * 78)
    log.info("PROVENANCE AND KNOWN LIMITATIONS, from the record")
    log.info("=" * 78)
    known = [
        ("14 mass ledger",
         "our own star count reached only 6% of the literature stellar "
         "density, so the luminous term is the PUBLISHED value, not ours. "
         "The channel is a literature synthesis with our own mass scale, not "
         "an independent measurement."),
        ("24 energy balance",
         "the sign-reversed mirror returned ZERO usable objects, because the "
         "sample was pre-selected for positive optical deficits and reversing "
         "the sign makes every fractional loss negative. The null rests on "
         "the four-order-of-magnitude energy gap alone."),
        ("22 SED vs spectrum",
         "GSP-Spec needs RVS spectra, so the 32k sample is restricted to "
         "bright stars and is not the same population the photometric "
         "channels constrain."),
        ("20 nearest 10 pc",
         "91 clean stars is a dossier, not a rate limit. Reported as such."),
        ("16 dark companions",
         "every candidate orbit has a period beyond the ~1000 d DR3 baseline "
         "and is extrapolated. Unresolvable until DR4."),
        ("13 coherent accel",
         "the original mark-permutation null was geometrically invalid and "
         "gave -3.9 sigma; corrected to -0.03 by rotating each vector about "
         "its own line of sight. Earlier sensitivity claims are superseded."),
        ("25 Planck excluded bin",
         "survivors still sit one grid step from the boundary, and PCCS2E "
         "photometry is unreliable by construction. Null is safe only "
         "because the mirror exceeds the signal."),
    ]
    for ch, note in known:
        log.info("\n  %-24s %s", ch, note)
        findings.append({"channel": ch, "issue": note, "source": "record"})

    # ---- the summary that matters ----------------------------------------
    log.info("\n" + "=" * 78)
    log.info("WHICH CHANNELS ACTUALLY CONSTRAIN SOMETHING?")
    log.info("=" * 78)
    strong, weak = [], []
    for name, n in sorted(N_TESTED.items(), key=lambda kv: -kv[1]):
        ul = 3.0 / n
        (strong if ul < 1e-3 else weak).append((name, n, ul))
    log.info("\n  informative (can reach rarer than 1 in 1,000):")
    for nm, n, ul in strong:
        log.info("    %-24s n=%-10s  rate UL %.1e", nm, f"{n:,}", ul)
    log.info("\n  low power (cannot):")
    for nm, n, ul in weak:
        log.info("    %-24s n=%-10s  rate UL %.1e  (1 in %s)",
                 nm, f"{n:,}", ul, f"{int(1/ul):,}")

    out = cfg.RESULT_DIR / "channel_audit.json"
    out.write_text(json.dumps({
        "n_findings": len(findings),
        "findings": findings,
        "informative_channels": [
            {"channel": nm, "n": n, "rate_upper_limit": ul}
            for nm, n, ul in strong],
        "low_power_channels": [
            {"channel": nm, "n": n, "rate_upper_limit": ul}
            for nm, n, ul in weak],
    }, indent=2))
    log.info("\nwrote %s  (%d findings)", out, len(findings))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
