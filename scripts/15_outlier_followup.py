#!/usr/bin/env python
"""Step 7: cross-check optical-deficit candidates against WISE and SIMBAD.

    run.sh scripts/15_outlier_followup.py --tag primary --max-objects 300

Two things are being asked of the mid-infrared, and they point opposite ways:

  * mid-IR EXCESS  -> a partial sphere radiating its waste heat isotropically,
    but ALSO the signature of every mundane contaminant that produces an
    optical deficit: debris discs, YSOs, edge-on discs, and reddening pockets.
    Excess is therefore evidence *against* a clean interpretation far more
    often than for one.

  * NO mid-IR excess -> consistent with the beamed-waste-heat case that
    motivates the optical channel in the first place, but equally consistent
    with the candidate being a photometric or extinction artefact.

Neither outcome is decisive alone, which is why the SIMBAD classification and
the negative-residual control list matter: any mechanism that produces a
symmetric tail should populate both sides equally.
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("followup")

# Expected main-sequence W1-W2 and W1-W3 colours are ~0 for Rayleigh-Jeans
# photospheres. Excess thresholds below are deliberately loose so that
# borderline cases are flagged rather than silently dropped.
# W1-W2 > 0.2 and W1-W3 > 0.5 are conventional debris-disc / YSO screens
# (e.g. Koenig & Leisawitz 2014, ApJ 791, 131, Sect. 3).
W1W2_EXCESS = 0.2
W1W3_EXCESS = 0.5
W1W4_EXCESS = 1.0


def wise_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for a, b, thr, name in [("wise_w1mpro", "wise_w2mpro", W1W2_EXCESS, "w1w2"),
                            ("wise_w1mpro", "wise_w3mpro", W1W3_EXCESS, "w1w3"),
                            ("wise_w1mpro", "wise_w4mpro", W1W4_EXCESS, "w1w4")]:
        if a in out and b in out:
            col = out[a] - out[b]
            out[f"{name}_colour"] = col
            out[f"{name}_excess"] = col > thr
        else:
            out[f"{name}_colour"] = np.nan
            out[f"{name}_excess"] = False
    out["any_ir_excess"] = (out[["w1w2_excess", "w1w3_excess",
                                 "w1w4_excess"]].any(axis=1))
    out["ir_measured"] = out[["w1w2_colour", "w1w3_colour",
                              "w1w4_colour"]].notna().any(axis=1)
    # "No excess" and "no data" are different statements and must not be
    # merged: an undetected W3/W4 makes a candidate LOOK beaming-consistent
    # when in fact nothing was measured.  The beamed-waste-heat case requires a
    # positive statement that the mid-IR was observed and found normal.
    out["ir_no_excess_measured"] = out["ir_measured"] & ~out["any_ir_excess"]
    out["ir_unknown"] = ~out["ir_measured"]
    return out


def query_simbad(ra, dec, radius_arcsec=3.0):
    """Classify candidates with a single SIMBAD TAP cone-join.

    A per-batch `Simbad.query_region` with a vector of coordinates does not
    return a reliable row->input mapping in current astroquery (the
    script_number_id field is gone), which silently produced an all-None
    classification.  A TAP query with an explicit input table keeps the
    association exact.
    """
    from astroquery.utils.tap.core import TapPlus

    n = len(ra)
    main = np.array([None] * n, dtype=object)
    otype = np.array([None] * n, dtype=object)
    sep = np.full(n, np.nan)

    tap = TapPlus(url="https://simbad.cds.unistra.fr/simbad/sim-tap")
    r_deg = radius_arcsec / 3600.0

    # SIMBAD's ADQL parser rejects UNION inside a subquery, so the input list
    # is expressed as an OR of cone predicates and matched back to the inputs
    # client-side by angular separation.
    chunk = 25
    for i0 in range(0, n, chunk):
        i1 = min(i0 + chunk, n)
        cones = " OR ".join(
            f"CONTAINS(POINT('ICRS', b.ra, b.dec), "
            f"CIRCLE('ICRS', {ra[j]:.8f}, {dec[j]:.8f}, {r_deg:.8f})) = 1"
            for j in range(i0, i1))
        q = (f"SELECT b.main_id, b.otype_txt, b.sp_type, b.ra, b.dec "
             f"FROM basic AS b WHERE b.ra IS NOT NULL AND ({cones})")
        for attempt in range(3):
            try:
                res = tap.launch_job(q).get_results()
                break
            except Exception as exc:                   # noqa: BLE001
                log.warning("SIMBAD TAP %d-%d attempt %d: %s",
                            i0, i1, attempt + 1, str(exc)[:120])
                res = None
                time.sleep(4 * (attempt + 1))
        if res is None or len(res) == 0:
            log.info("SIMBAD %d/%d (no matches in chunk)", i1, n)
            continue
        rdf = res.to_pandas()
        rr = rdf["ra"].to_numpy(dtype=float)
        rd = rdf["dec"].to_numpy(dtype=float)
        for j in range(i0, i1):
            # small-angle separation is fine at 3 arcsec
            cosd = np.cos(np.radians(dec[j]))
            ds = np.hypot((rr - ra[j]) * cosd, rd - dec[j])
            if len(ds) == 0:
                continue
            m = int(np.argmin(ds))
            if ds[m] <= r_deg:
                sep[j] = ds[m]
                mid, oty = rdf.iloc[m].get("main_id"), rdf.iloc[m].get("otype_txt")
                main[j] = mid.decode() if isinstance(mid, bytes) else mid
                otype[j] = oty.decode() if isinstance(oty, bytes) else oty
        log.info("SIMBAD %d/%d (%d rows returned)", i1, n, len(rdf))
    return main, otype


CONTAMINANT_OTYPES = {
    # young stellar objects and their discs -- the classic optical-deficit
    # false positive, since circumstellar material genuinely dims the star
    "YSO", "Y*O", "Y*?", "TT*", "TT?", "Or*", "pr*", "out", "HH",
    # emission-line and shell stars
    "Em*", "Be*", "WR*", "s*b", "sg*",
    # evolved stars that leave the dwarf sequence
    "PN", "PN?", "RG*", "AGB", "C*", "S*", "HB*", "pA*", "WD*", "WD?",
    # pulsators and irregulars
    "LP?", "LP*", "Mi*", "Ce*", "RR*", "dS*", "V*", "Ir*", "cC*", "gam",
    # SPOTTED / ACTIVE ROTATORS. These matter more than they look: starspots
    # dim the optical continuum at roughly constant K_s, which is precisely
    # the signal being searched for. BY Dra and RS CVn stars are therefore an
    # irreducible astrophysical background, not a data-quality issue.
    "BY*", "RS*", "Ro*", "Fl*", "UV*", "Er*", "Ir*",
    # multiples
    "EB*", "SB*", "Al*", "bL*", "WU*", "**", "SB?", "EB?",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--max-objects", type=int, default=300)
    ap.add_argument("--no-simbad", action="store_true")
    args = ap.parse_args()

    path = cfg.RESULT_DIR / f"outliers_{args.tag}.csv"
    if not path.exists():
        log.error("%s not found -- run scripts/14_analyse.py first", path)
        return 1
    o = pd.read_csv(path)
    log.info("%d outliers (%d positive)", len(o),
             int(o["side"].str.startswith("positive").sum()))

    o = wise_flags(o)

    # Rank by |residual| and follow up the strongest, both sides.
    o["abs_sigma"] = o["residual_sigma"].abs()
    pos = o[o["side"].str.startswith("positive")].nlargest(
        args.max_objects, "abs_sigma")
    neg = o[o["side"].str.startswith("negative")].nlargest(
        args.max_objects, "abs_sigma")
    sub = pd.concat([pos, neg], ignore_index=True)

    if not args.no_simbad:
        log.info("querying SIMBAD for %d objects ...", len(sub))
        main_id, otype = query_simbad(sub["ra"].to_numpy(),
                                      sub["dec"].to_numpy())
        sub["simbad_main_id"] = main_id
        sub["simbad_otype"] = otype
    else:
        sub["simbad_main_id"] = None
        sub["simbad_otype"] = None

    sub["known_contaminant"] = sub["simbad_otype"].isin(CONTAMINANT_OTYPES)
    # A "clean" candidate must have the mid-IR actually MEASURED and normal.
    # Candidates with no WISE photometry are counted separately, not promoted.
    sub["clean_candidate"] = (sub["side"].str.startswith("positive")
                              & sub["ir_no_excess_measured"]
                              & ~sub["known_contaminant"])
    sub["unverifiable_no_ir_data"] = (sub["side"].str.startswith("positive")
                                      & sub["ir_unknown"])

    sub.to_csv(cfg.RESULT_DIR / f"outlier_followup_{args.tag}.csv", index=False)

    print("\n=== WISE mid-infrared screening ===")
    for side, g in sub.groupby("side"):
        n = len(g)
        print(f"  {side}: n={n}")
        print(f"    with any WISE colour measured : {int(g['ir_measured'].sum())}")
        for c in ["w1w2_excess", "w1w3_excess", "w1w4_excess", "any_ir_excess"]:
            print(f"    {c:16s} : {int(g[c].sum())}")

    print("\n=== SIMBAD classification (top types) ===")
    if sub["simbad_otype"].notna().any():
        for side, g in sub.groupby("side"):
            print(f"  -- {side}")
            vc = g["simbad_otype"].value_counts(dropna=True).head(12)
            for k, v in vc.items():
                mark = "  <-- contaminant" if k in CONTAMINANT_OTYPES else ""
                print(f"     {str(k):12s} {v:5d}{mark}")
            print(f"     unmatched in SIMBAD: {int(g['simbad_otype'].isna().sum())}")
    else:
        print("  (SIMBAD not queried)")

    n_clean = int(sub["clean_candidate"].sum())
    n_pos_f = int(sub["side"].str.startswith("positive").sum())
    n_unver = int(sub["unverifiable_no_ir_data"].sum())
    print(f"\n=== surviving candidates ===")
    print(f"  positive outliers followed up            : {n_pos_f}")
    print(f"  mid-IR MEASURED and normal, not a known contaminant : {n_clean}")
    print(f"  no WISE photometry at all (unverifiable) : {n_unver}")
    print(f"  rejected by IR excess or SIMBAD class    : "
          f"{n_pos_f - n_clean - n_unver}")
    if n_clean:
        cols = [c for c in ["source_id", "ra", "dec", "residual_sigma",
                            "implied_f", "w1w2_colour", "w1w3_colour",
                            "simbad_otype", "ruwe", "A_0"] if c in sub.columns]
        print(sub[sub["clean_candidate"]].nlargest(20, "abs_sigma")[cols]
              .to_string(index=False))
    print(f"\n  full table: results/outlier_followup_{args.tag}.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
