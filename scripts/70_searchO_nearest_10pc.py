#!/usr/bin/env python
"""Search O: the nearest 10 parsecs, examined one star at a time.

    run.sh scripts/70_searchO_nearest_10pc.py

THE HOLE
--------
The sample that every channel in this project runs on carries the cut

    10 < d < 500 pc

so all 315 Gaia sources inside 10 parsecs were excluded from all of it. The
lower bound exists for good instrumental reasons -- nearby stars are bright,
and bright means Gaia saturation below G ~ 6 and 2MASS saturation below
Ks ~ 4 -- but the effect is that the volume closest to us has never been
looked at by any of the nineteen channels.

That is the wrong volume to have skipped. Every limit published here is a
statement about a shell from 10 to 500 pc, and says nothing whatsoever about
the Solar neighbourhood.

WHY THIS RUNS DIFFERENTLY FROM EVERY OTHER CHANNEL
--------------------------------------------------
315 objects is far too few for the population statistics the rest of the
project relies on -- a single anomalous system would be one count in a bin,
indistinguishable from noise, and the mirror-control logic has nothing to work
with. So this does not compute a limit. It builds a dossier: every measurement
we can assemble for every star, ranked, so that individual objects can be
inspected rather than averaged away.

This is the right method when the question is "is there one here", and the
wrong method for "how common are they". It complements the other channels
rather than adding to them.

THE ONE ADVANTAGE OF THIS VOLUME
--------------------------------
Extinction inside 10 pc is essentially zero. The dominant systematic of the
entire project -- getting the dust right -- simply does not apply, so a
deficit measured here cannot be blamed on the extinction map.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchO")

MIRRORS = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
    ("ESA", "https://gea.esac.esa.int/tap-server/tap"),
]

# Saturation limits. Below these the photometry is not trustworthy and the
# star is reported but excluded from ranking.
G_SAT = 6.0
KS_SAT = 4.5
W1_SAT = 8.0
W2_SAT = 7.0

QUERY = """
SELECT source_id, ra, dec, l, b,
       parallax, parallax_error, parallax_over_error,
       pmra, pmdec, radial_velocity, ruwe,
       astrometric_excess_noise, astrometric_excess_noise_sig,
       ipd_frac_multi_peak, duplicated_source,
       phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp,
       phot_bp_rp_excess_factor, phot_variable_flag,
       non_single_star, teff_gspphot, logg_gspphot, mh_gspphot
FROM gaiadr3.gaia_source
WHERE parallax > 100
"""


def pick_mirror():
    from astroquery.utils.tap.core import TapPlus
    for name, url in MIRRORS:
        try:
            tap = TapPlus(url=url)
            tap.launch_job("SELECT TOP 1 source_id "
                           "FROM gaiadr3.gaia_source").get_results()
            log.info("using %s", name)
            return tap
        except Exception as exc:
            log.warning("%s unusable: %s", name, str(exc)[:110])
    raise RuntimeError("no mirror responded")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cache = cfg.RAW_DIR / "nearest_10pc.parquet"
    if cache.exists() and not args.force:
        d = pd.read_parquet(cache)
        log.info("cached: %d sources", len(d))
    else:
        tap = pick_mirror()
        d = tap.launch_job(QUERY).get_results().to_pandas()
        log.info("pulled %d sources within 10 pc", len(d))

        # 2MASS and WISE by position, since the mirror lacks the DR1 aux tables
        from astroquery.xmatch import XMatch
        from astropy.table import Table
        import astropy.units as u

        # These are the nearest stars, so they are also the fastest moving.
        # 2MASS observed ~1999 and Gaia's reference epoch is 2016.0, and
        # Barnard's Star covers 165 arcsec in that gap. Cross-matching at the
        # Gaia epoch silently loses exactly the nearest objects, so propagate
        # each position back to the survey's own epoch first.
        def at_epoch(df, epoch):
            dt = epoch - 2016.0
            dec_rad = np.radians(df["dec"].to_numpy(float))
            pmra = df["pmra"].fillna(0).to_numpy(float)     # already *cos(dec)
            pmdec = df["pmdec"].fillna(0).to_numpy(float)
            ra = df["ra"].to_numpy(float) + \
                (pmra / np.maximum(np.cos(dec_rad), 1e-6)) * dt / 3.6e6
            dec = df["dec"].to_numpy(float) + pmdec * dt / 3.6e6
            return pd.DataFrame({"source_id": df["source_id"].to_numpy(),
                                 "ra": ra, "dec": dec})

        for cat, cols, pre, epoch in (
                ("vizier:II/246/out",
                 {"Kmag": "tmass_ks_m", "Jmag": "tmass_j_m",
                  "Hmag": "tmass_h_m", "Qflg": "tmass_ph_qual"}, "2MASS", 1999.5),
                ("vizier:II/328/allwise",
                 {"W1mag": "W1", "W2mag": "W2", "W3mag": "W3",
                  "W4mag": "W4"}, "WISE", 2010.5)):
            src = Table.from_pandas(at_epoch(d, epoch))
            t0 = time.time()
            r = XMatch.query(cat1=src, cat2=cat, max_distance=10 * u.arcsec,
                             colRA1="ra", colDec1="dec").to_pandas()
            log.info("  %s: %d matches (%.1fs)", pre, len(r), time.time() - t0)
            keep = {"source_id": "source_id", "angDist": f"{pre}_xm_dist"}
            keep.update({k: v for k, v in cols.items() if k in r.columns})
            r = r[list(keep)].rename(columns=keep)
            r = r.sort_values(f"{pre}_xm_dist").drop_duplicates("source_id")
            d = d.merge(r, on="source_id", how="left")

        cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
        d.to_parquet(cache, index=False)
        log.info("wrote %s", cache)

    # ---- geometry. Extinction inside 10 pc is negligible. ----------------
    dist_pc = 1000.0 / d["parallax"].to_numpy(float)
    mu = 5.0 * np.log10(dist_pc) - 5.0
    d = d.assign(dist_pc=dist_pc)
    M_G = d["phot_g_mean_mag"].to_numpy(float) - mu
    M_Ks = d["tmass_ks_m"].to_numpy(float) - mu if "tmass_ks_m" in d else np.full(len(d), np.nan)
    d = d.assign(M_G=M_G, M_Ks=M_Ks)

    # ---- fiducial from the main sample, same estimator -------------------
    ref = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                          columns=["M_G", "M_Ks", "ruwe", "cstar_nsigma"])
    ok = (ref["M_G"].notna() & ref["M_Ks"].notna()
          & (ref["ruwe"] < 1.4) & (ref["cstar_nsigma"].abs() < 3))
    x, y = ref.loc[ok, "M_Ks"].to_numpy(float), ref.loc[ok, "M_G"].to_numpy(float)
    best = None
    for deg in (3, 4, 5):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, sigma, deg = best
    log.info("fiducial from %d clean stars: degree %d, sigma %.4f mag",
             int(ok.sum()), deg, sigma)

    in_range = np.isfinite(M_Ks) & (M_Ks > 3.0) & (M_Ks < 8.0)
    resid = np.where(in_range, M_G - np.polyval(coef, np.clip(M_Ks, 3.0, 8.0)),
                     np.nan)
    d = d.assign(residual=resid, resid_nsigma=resid / sigma)

    # ---- saturation and quality flags ------------------------------------
    g = d["phot_g_mean_mag"].to_numpy(float)
    ks = d["tmass_ks_m"].to_numpy(float) if "tmass_ks_m" in d else np.full(len(d), np.nan)
    w1 = d["W1"].to_numpy(float) if "W1" in d else np.full(len(d), np.nan)
    w2 = d["W2"].to_numpy(float) if "W2" in d else np.full(len(d), np.nan)

    # Saturation is per band. WISE saturates on almost every star this close,
    # which must not remove them from the OPTICAL test -- that needs only G
    # and Ks. Blanket-flagging on any band emptied the sample entirely.
    sat_opt = (g < G_SAT) | (ks < KS_SAT)
    sat_ir = (w1 < W1_SAT) | (w2 < W2_SAT)
    d = d.assign(saturated_optical=sat_opt, saturated_ir=sat_ir)
    sat = sat_opt
    log.info("within 10 pc: %d sources, %d with usable M_Ks", len(d),
             int(in_range.sum()))
    log.info("  saturated in G or Ks (optical test): %d", int(sat_opt.sum()))
    log.info("  saturated in W1 or W2 (IR test)    : %d", int(sat_ir.sum()))

    # ---- mid-infrared excess: does anything re-radiate? ------------------
    # W1-W2 for a bare photosphere is ~0 for these types; a warm shell is red.
    w12 = w1 - w2
    d = d.assign(W1_W2=w12)
    good_ir = np.isfinite(w12) & ~sat_ir
    if good_ir.sum() > 20:
        med12 = float(np.median(w12[good_ir]))
        sig12 = float(st.robust_sigma(w12[good_ir]))
        d = d.assign(w12_nsigma=(w12 - med12) / max(sig12, 1e-6))
        log.info("W1-W2 among unsaturated: median %+.3f, robust sigma %.3f",
                 med12, sig12)
    else:
        d = d.assign(w12_nsigma=np.nan)

    # ---- shared infrared counterparts -------------------------------------
    # Inside 10 pc this is the dominant failure mode, and it is worst exactly
    # here: Gaia resolves a binary at sub-arcsecond scale while 2MASS (4") and
    # WISE (6") do not, so BOTH components inherit the same combined infrared
    # magnitude. That Ks is brighter than either star alone, which the
    # estimator reads as a large G deficit in both.
    from scipy.spatial import cKDTree
    ra_r = np.radians(d["ra"].to_numpy(float))
    dec_r = np.radians(d["dec"].to_numpy(float))
    xyz = np.column_stack([np.cos(dec_r) * np.cos(ra_r),
                           np.cos(dec_r) * np.sin(ra_r), np.sin(dec_r)])
    tree = cKDTree(xyz)
    # 6 arcsec chord on the unit sphere
    pairs = tree.query_ball_point(xyz, r=2 * np.sin(np.radians(6.0 / 3600) / 2))
    n_gaia_within_6as = np.array([len(p) - 1 for p in pairs])

    # and directly: do two Gaia sources carry identical infrared photometry?
    key = (d["W1"].round(3).astype(str) + "_" + d["W2"].round(3).astype(str))
    shared_ir = key.map(key.value_counts()) > 1
    shared_ir &= d["W1"].notna()

    d = d.assign(n_gaia_within_6as=n_gaia_within_6as,
                 shares_ir_counterpart=shared_ir.to_numpy())
    log.info("sources with another Gaia star inside the 6\" WISE beam: %d",
             int((n_gaia_within_6as > 0).sum()))
    log.info("sources sharing identical W1/W2 with another Gaia source: %d",
             int(shared_ir.sum()))

    # ---- the dossier ------------------------------------------------------
    rank_ok = in_range & ~sat & np.isfinite(resid)
    log.info("rankable (unsaturated, in the MS box): %d", int(rank_ok.sum()))

    dd = d[rank_ok].copy()
    dd["abs_resid_nsigma"] = dd["resid_nsigma"].abs()

    # An absorber DIMS, so the deficit tail is the one that matters; the
    # bright tail is reported alongside as this volume's own false-positive
    # indicator.
    dim = dd[dd["resid_nsigma"] > 3.0].sort_values("resid_nsigma",
                                                   ascending=False)
    bright = dd[dd["resid_nsigma"] < -3.0].sort_values("resid_nsigma")
    ir_red = dd[dd["w12_nsigma"] > 3.0].sort_values("w12_nsigma",
                                                    ascending=False)
    hi_ruwe = dd[dd["ruwe"] > 1.4].sort_values("ruwe", ascending=False)

    log.info("")
    log.info("optical deficit > 3 sigma : %d", len(dim))
    log.info("optical excess  > 3 sigma : %d  (this volume's own noise level)",
             len(bright))
    log.info("W1-W2 red       > 3 sigma : %d", len(ir_red))
    log.info("RUWE > 1.4                : %d", len(hi_ruwe))

    if len(dim):
        log.info("")
        log.info("most optically deficient systems:")
        for _, r in dim.head(10).iterrows():
            log.info("  %d  d=%.2f pc  M_Ks=%.2f  deficit=%+.3f mag (%.1f sig)"
                     "  RUWE=%.2f  W1-W2=%s  %s",
                     int(r["source_id"]), r["dist_pc"], r["M_Ks"],
                     r["residual"], r["resid_nsigma"], r["ruwe"],
                     f"{r['W1_W2']:+.3f}" if np.isfinite(r["W1_W2"]) else "  n/a",
                     str(r.get("phot_variable_flag", "")))

    # ---- is the deficit tail just blended multiples? ---------------------
    blended = (dd["shares_ir_counterpart"].to_numpy()
               | (dd["n_gaia_within_6as"].to_numpy() > 0)
               | (dd["ruwe"].to_numpy() > 1.4))
    clean_dd = dd[~blended]
    n_dim_clean = int((clean_dd["resid_nsigma"] > 3.0).sum())
    n_bright_clean = int((clean_dd["resid_nsigma"] < -3.0).sum())
    log.info("")
    log.info("after removing shared-IR / close-neighbour / high-RUWE systems:")
    log.info("  clean rankable : %d", len(clean_dd))
    log.info("  deficit > 3sig : %d", n_dim_clean)
    log.info("  excess  > 3sig : %d", n_bright_clean)
    if n_dim_clean:
        log.info("  surviving deficient systems:")
        for _, r in clean_dd[clean_dd["resid_nsigma"] > 3.0].sort_values(
                "resid_nsigma", ascending=False).iterrows():
            log.info("    %d  d=%.2f pc  M_Ks=%.2f  deficit=%+.3f (%.1f sig)"
                     "  RUWE=%.2f", int(r["source_id"]), r["dist_pc"],
                     r["M_Ks"], r["residual"], r["resid_nsigma"], r["ruwe"])

    # ---- identify what survives ------------------------------------------
    survivors = clean_dd[clean_dd["resid_nsigma"] > 3.0]
    if len(survivors):
        try:
            from astroquery.simbad import Simbad
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            sb = Simbad()
            sb.add_votable_fields("otype", "sptype")
            log.info("")
            log.info("SIMBAD identification of the survivors:")
            for _, r in survivors.iterrows():
                c = SkyCoord(r["ra"] * u.deg, r["dec"] * u.deg)
                try:
                    res = sb.query_region(c, radius=20 * u.arcsec)
                except Exception:
                    res = None
                if res is None or len(res) == 0:
                    log.info("  %d : no SIMBAD match", int(r["source_id"]))
                    continue
                row = res[0]
                name = str(row[res.colnames[0]])
                otype = str(row["otype"]) if "otype" in res.colnames else "?"
                sptype = str(row["sptype"]) if "sptype" in res.colnames else "?"
                log.info("  %d : %s  type=%s  sp=%s",
                         int(r["source_id"]), name, otype, sptype)
        except Exception as exc:
            log.warning("SIMBAD lookup unavailable: %s", str(exc)[:100])

    # ---- verdict ----------------------------------------------------------
    n_dim, n_bright = len(dim), len(bright)
    if n_dim_clean <= n_bright_clean:
        verdict = (
            f"NULL after the blend cut. Of {n_dim} apparently deficient "
            f"systems inside 10 pc, {n_dim - n_dim_clean} carry a shared "
            f"infrared counterpart, a Gaia neighbour inside the 6 arcsec WISE "
            f"beam, or RUWE > 1.4. That is aperture mismatch in its most "
            f"extreme form, and this volume is where it is worst: Gaia "
            f"resolves nearby binaries that 2MASS and WISE blend, so both "
            f"components inherit one combined Ks brighter than either star. "
            f"Among the {len(clean_dd)} clean stars that remain, "
            f"{n_dim_clean} are deficient against {n_bright_clean} "
            f"over-luminous -- symmetric at this sample size, i.e. noise. The "
            f"Solar neighbourhood contains no optical-deficit candidate.")
    elif n_dim_clean == 0 and n_dim > 0:
        verdict = (
            f"RESOLVED MULTIPLES. The {n_dim} deficient systems inside 10 pc "
            f"all carry a shared infrared counterpart, a Gaia neighbour inside "
            f"the 6 arcsec WISE beam, or RUWE > 1.4. None survives. This is "
            f"aperture mismatch in its most extreme form: Gaia resolves nearby "
            f"binaries that 2MASS and WISE blend, so both components inherit "
            f"one combined Ks that is brighter than either star, and the "
            f"estimator reads that as a deficit in both. The effect is worst "
            f"in this volume precisely because nearby binaries are the widest "
            f"on the sky. After the cut, {len(clean_dd)} clean stars inside "
            f"10 pc show no optical deficit at all.")
    elif n_dim == 0:
        verdict = (
            f"NULL. Of {int(rank_ok.sum())} unsaturated main-sequence stars "
            f"inside 10 pc, none shows an optical deficit beyond 3 sigma. The "
            f"volume excluded from every other channel contains no candidate. "
            f"Note this is a dossier, not a limit: 315 objects cannot bound a "
            f"rate, and the value of the result is that each was inspected "
            f"rather than averaged.")
    elif n_dim <= n_bright:
        verdict = (
            f"NULL. {n_dim} systems are deficient and {n_bright} are equally "
            f"over-luminous, so the tail is symmetric scatter at this sample "
            f"size rather than a population that only dims.")
    else:
        verdict = (
            f"{n_dim} deficient against {n_bright} over-luminous inside 10 pc. "
            f"With this few objects that asymmetry is not statistically "
            f"meaningful on its own, but these are the nearest, best-measured, "
            f"least-extinguished stars available and each is individually "
            f"checkable against decades of literature. Listed in the CSV for "
            f"identification.")

    print(f"\n{'='*72}")
    print("SEARCH O: THE NEAREST 10 PARSECS")
    print(f"{'='*72}")
    print(f"  Gaia sources within 10 pc          : {len(d)}")
    print(f"  with usable M_Ks                   : {int(in_range.sum())}")
    print(f"  saturated (excluded from ranking)  : {int(sat.sum())}")
    print(f"  rankable                           : {int(rank_ok.sum())}")
    print()
    print(f"  optical deficit > 3 sigma          : {n_dim}")
    print(f"  optical excess  > 3 sigma (mirror) : {n_bright}")
    print(f"  W1-W2 red > 3 sigma                : {len(ir_red)}")
    print(f"  RUWE > 1.4                         : {len(hi_ruwe)}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchO_nearest10pc_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "n_within_10pc": int(len(d)),
        "n_usable_mks": int(in_range.sum()),
        "n_saturated": int(sat.sum()),
        "n_rankable": int(rank_ok.sum()),
        "fiducial_sigma": float(sigma),
        "n_deficit_3sig": n_dim,
        "n_excess_3sig": n_bright,
        "n_ir_red_3sig": len(ir_red),
        "n_high_ruwe": len(hi_ruwe),
        "saturation_limits": {"G": G_SAT, "Ks": KS_SAT,
                              "W1": W1_SAT, "W2": W2_SAT},
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)

    cols = [c for c in ["source_id", "ra", "dec", "dist_pc", "phot_g_mean_mag",
                        "bp_rp", "M_G", "M_Ks", "residual", "resid_nsigma",
                        "W1", "W2", "W1_W2", "w12_nsigma", "ruwe",
                        "astrometric_excess_noise", "phot_variable_flag",
                        "non_single_star", "teff_gspphot", "saturated"]
            if c in d.columns]
    d[cols].sort_values("resid_nsigma", ascending=False).to_csv(
        cfg.RESULT_DIR / f"searchO_nearest10pc_dossier_{args.tag}.csv",
        index=False)
    log.info("wrote the full 10 pc dossier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
