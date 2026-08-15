#!/usr/bin/env python
"""Independent mass anchor: how well does M_Ks stand in for stellar mass?

    run.sh scripts/16_mass_anchor.py

Why this exists
---------------
The whole method treats M_Ks as "how big is this star really".  If an absorber
dims K_s as well as G, that anchor moves and the inferred deficit is wrong.  A
mass measured *dynamically* does not move, so the scatter of M_Ks about a
mass-M_Ks relation bounds how much unmodelled K_s dimming the population can
carry.

THE UNAVOIDABLE TENSION
-----------------------
Dynamical masses come from binaries.  The science sample deliberately removes
binaries (RUWE < 1.4, non_single_star = 0) because they are the dominant
contaminant.  The control population is therefore *not* the science population,
and the bound it gives is an argument by analogy, not a direct measurement on
the same stars.  This is stated as a limitation rather than papered over.
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
from pipeline import extinction as ext
from pipeline.tap import run_adql, tap_client

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("anchor")

MASS_QUERY = """
SELECT
  bm.source_id, bm.m1, bm.m1_lower, bm.m1_upper,
  bm.m2, bm.m2_lower, bm.m2_upper, bm.fluxratio, bm.combination_method,
  g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_over_error, g.ruwe,
  g.phot_g_mean_mag, g.bp_rp, g.nu_eff_used_in_astrometry, g.pseudocolour,
  g.ecl_lat, g.astrometric_params_solved,
  tm.ks_m AS tmass_ks_m, tm.ks_msigcom AS tmass_ks_msigcom,
  tm.j_m AS tmass_j_m, tm.ph_qual AS tmass_ph_qual
FROM gaiadr3.binary_masses AS bm
JOIN gaiadr3.gaia_source AS g ON g.source_id = bm.source_id
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON xm.source_id = g.source_id
JOIN gaiadr3.tmass_psc_xsc_join AS xj
  ON xj.clean_tmass_psc_xsc_oid = xm.clean_tmass_psc_xsc_oid
JOIN gaiadr1.tmass_original_valid AS tm
  ON tm.designation = xj.original_psc_source_id
WHERE bm.m1 IS NOT NULL
  AND g.parallax_over_error > 20
  AND g.parallax > 2
  AND tm.ks_m IS NOT NULL
  AND tm.ks_msigcom < 0.05
"""


def table_exists(name: str) -> bool:
    """Check by listing every table, not with a WHERE clause.

    A WHERE on TAP_SCHEMA.tables triggers a query-rewriting bug in the ESA
    service ("missing FROM-clause entry for table cte_authorisation_NNN") that
    surfaces as a 500 and looks exactly like "table not present". Selecting
    `description` is also avoided: some rows carry non-ASCII bytes that the
    VOTable BINARY parser decodes as ASCII and dies on.
    """
    try:
        rows = tap_client().launch_job(
            "SELECT table_name FROM TAP_SCHEMA.tables").get_results()
        return name in {str(r["table_name"]) for r in rows}
    except Exception as exc:                            # noqa: BLE001
        log.warning("table listing failed (%s); attempting the query anyway", exc)
        return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not table_exists("gaiadr3.binary_masses"):
        log.error("gaiadr3.binary_masses not present on this TAP service")
        return 1

    df = run_adql(MASS_QUERY, name="binary_masses", force=args.force)
    log.info("%d systems with dynamical masses and clean 2MASS Ks", len(df))
    if len(df) < 50:
        log.error("too few systems for a useful anchor")
        return 1

    # Distances and absolute Ks. Extinction is small for this bright,
    # nearby set, but apply it consistently with the main analysis.
    from pipeline import sample as smp
    df = smp.add_astrometry(df)
    a0 = ext.query_a0("edenhofer23", df["l"].to_numpy(float),
                      df["b"].to_numpy(float), df["dist_pc"].to_numpy(float))
    a_ks = ext.deredden("Ks", np.where(np.isfinite(a0), a0, 0.0),
                        df["bp_rp"].to_numpy(float), law="fitz19")
    df["A_0"] = a0
    df["M_Ks"] = (df["tmass_ks_m"].to_numpy(float) - df["dist_mod"].to_numpy(float)
                  - np.where(np.isfinite(a0), a_ks, np.nan))

    good = (np.isfinite(df["M_Ks"]) & df["m1"].notna()
            & (df["m1"] > 0.15) & (df["m1"] < 2.0)
            & (df["tmass_ph_qual"].str[2] == "A"))
    d = df[good].reset_index(drop=True)
    log.info("%d systems in the usable mass range", len(d))
    if len(d) < 30:
        log.error("too few after quality cuts")
        return 1

    def analyse(dd: pd.DataFrame, label: str) -> dict:
        """Fit M_Ks(log mass) and split the scatter into mass error vs the rest."""
        x = np.log10(dd["m1"].to_numpy(float))
        y = dd["M_Ks"].to_numpy(float)
        results = {}
        for deg in (1, 2, 3):
            cc = np.polyfit(x, y, deg)
            rr = y - np.polyval(cc, x)
            results[deg] = (cc, 1.4826 * np.median(np.abs(rr - np.median(rr))))
        deg = min(results, key=lambda k: results[k][1])
        cc, s_mks = results[deg]

        m_err = ((dd["m1_upper"] - dd["m1_lower"]) / 2.0).to_numpy(float)
        rel_m_err = float(np.nanmedian(m_err / dd["m1"].to_numpy(float)))
        dMks_dlogm = float(np.polyval(np.polyder(cc), np.median(x)))
        s_from_mass = abs(dMks_dlogm) * rel_m_err / np.log(10)
        s_ks = float(np.sqrt(max(s_mks ** 2 - s_from_mass ** 2, 0.0)))
        return {
            "subset": label,
            "n_systems": int(len(dd)),
            "poly_degree": int(deg),
            "scatter_M_Ks_about_mass_relation_mag": float(s_mks),
            "median_relative_mass_error": rel_m_err,
            "dMKs_dlog10M": dMks_dlogm,
            "scatter_attributable_to_mass_error_mag": float(s_from_mass),
            "residual_M_Ks_scatter_mag": s_ks,
            "flat_absorber_bound_frac": float(1 - 10 ** (-s_ks / 2.5)),
        }

    subsets = {"all systems": d}
    # The dominant nuisance is light from the SECONDARY, which inflates M_Ks at
    # fixed primary mass and has nothing to do with any absorber.  Where the
    # flux ratio is measured, restricting to faint secondaries removes most of
    # it and tightens the bound substantially.
    if "fluxratio" in d and d["fluxratio"].notna().any():
        for thr in (0.30, 0.10):
            sub = d[d["fluxratio"].notna() & (d["fluxratio"] < thr)]
            if len(sub) >= 50:
                subsets[f"fluxratio < {thr:.2f}"] = sub

    rows = [analyse(dd, lab) for lab, dd in subsets.items()]
    best = min(rows, key=lambda r: r["flat_absorber_bound_frac"])

    out = {
        "mass_range_msun": [float(d["m1"].min()), float(d["m1"].max())],
        "subsets": rows,
        "best_subset": best["subset"],
        "flat_absorber_bound_frac": best["flat_absorber_bound_frac"],
        "residual_M_Ks_scatter_mag": best["residual_M_Ks_scatter_mag"],
        "n_systems": best["n_systems"],
        "caveats": [
            "control population is BINARIES, which the science sample "
            "explicitly excludes; this is an argument by analogy",
            "unmodelled secondary light inflates M_Ks at fixed primary mass, "
            "so the bound is conservative",
            "a flat absorber affecting every star IDENTICALLY is absorbed into "
            "the fitted relation and is unbounded by any self-calibrated method",
        ],
    }
    (cfg.RESULT_DIR / "mass_anchor.json").write_text(json.dumps(out, indent=2))
    d.to_parquet(cfg.DERIVED_DIR / "mass_anchor.parquet", index=False)

    print("\n=== independent mass anchor ===")
    print(pd.DataFrame(rows).to_string(index=False,
                                       float_format=lambda v: f"{v:10.4f}"))
    print(f"\n  mass range: {out['mass_range_msun'][0]:.2f} - "
          f"{out['mass_range_msun'][1]:.2f} Msun")
    for c in out["caveats"]:
        print(f"  caveat: {c}")
    print(f"\n  A population-wide flat (grey) absorber larger than "
          f"{out['flat_absorber_bound_frac']:.3f}")
    print("  would show up as excess scatter in M_Ks at fixed dynamical mass, "
          "unless it\n  affected every star identically -- in which case it is "
          "absorbed into the\n  fitted relation and remains unbounded by ANY "
          "self-calibrated method.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
