#!/usr/bin/env python
"""Search C: radio counterparts to optically-dimmed stars.

    run.sh scripts/47_searchC_radio.py --tag primary

THE HYPOTHESIS
--------------
A main-sequence star should not be a radio source.  The radio luminosity of
quiescent FGK/M dwarfs is far below the detection threshold of any existing
wide-field survey (Gudel 1994, ApJ 429, L29; Gudel 2002, ARA&A 40, 217);
detections at 1.4 GHz require active or close-binary stars, which RUWE and
variability cuts have already removed.

A radio detection coincident with a star that is also anomalously dim in the
optical is therefore interesting for several reasons:

  * Non-thermal synchrotron emission from engineered structures, analogous to
    waste-heat searches at mid-IR but in a completely different band where the
    expected background is zero for a clean main-sequence star.

  * Communication leakage, beamed or isotropic, at GHz frequencies -- the
    canonical SETI band.

  * An active magnetosphere hidden by a RUWE < 1.4 solution.  This is a
    mundane explanation, but it would still be a useful contaminant flag for
    the optical-deficit candidate list.

  * A chance alignment with a background AGN or radio galaxy.  The expected
    rate of such alignments is quantified via a control sample below.

SURVEYS QUERIED
---------------
Three wide-area 1--3 GHz surveys available through VizieR TAP:

    NVSS   (1.4 GHz, Condon+1998, AJ 115, 1693)    -- VIII/65/nvss
    VLASS  (3 GHz, Lacy+2020, PASP 132, 035001)     -- VIII/110/vlass2
    FIRST  (1.4 GHz, Helfand+2015, ApJ 801, 26)     -- VIII/92/first14

Cross-match radius is 5 arcsec, applied as a cone search on VizieR TAP.

CONTROL SAMPLE
--------------
To estimate the false-positive rate from chance alignments, we query the same
surveys with a control sample of stars that are photometrically similar but
show no optical deficit (residual within 1 sigma of zero).  The ratio of
match rates gives the excess, if any.
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchC_radio")

# --------------------------------------------------------------------------
# VizieR TAP configuration
# --------------------------------------------------------------------------

VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"

RADIO_CATALOGS = {
    "NVSS": {
        "table": "\"VIII/65/nvss\"",
        "ra_col": "RAJ2000",
        "dec_col": "DEJ2000",
        "flux_col": "S1_4",       # VizieR normalises "S1.4" -> "S1_4" in TAP
        "freq_ghz": 1.4,
        "cite": "Condon+1998",
    },
    # VLASS is not available in VizieR TAP as of 2026-08.  VIII/110 is GLEAM-X,
    # not VLASS.  We proceed with NVSS + FIRST (both 1.4 GHz) which together
    # cover the northern sky to Dec > -40 degrees.
    "FIRST": {
        "table": "\"VIII/92/first14\"",
        "ra_col": "RAJ2000",
        "dec_col": "DEJ2000",
        "flux_col": "Fpeak",      # mJy at 1.4 GHz
        "freq_ghz": 1.4,
        "cite": "Helfand+2015",
    },
}

MATCH_RADIUS_ARCSEC = 5.0
MATCH_RADIUS_DEG = MATCH_RADIUS_ARCSEC / 3600.0

MAX_CANDIDATES = 5000        # top optically-dimmed stars to cross-match
CONTROL_N = 5000             # size of the non-dimmed control sample
CHUNK_SIZE = 50              # stars per TAP query (keep ADQL short for VizieR)
MAX_RETRIES = 3
RETRY_BASE_SLEEP = 10.0     # seconds

RNG_SEED = 47_2026


# --------------------------------------------------------------------------
# TAP query helpers
# --------------------------------------------------------------------------

def _vizier_tap():
    """Get a VizieR TapPlus client."""
    from astroquery.utils.tap.core import TapPlus
    return TapPlus(url=VIZIER_TAP)


def _cone_union_adql(cat: dict, ra_arr, dec_arr) -> str:
    """Build an ADQL query that ORs together cone searches for a batch.

    VizieR column names like ``S1.4`` contain dots, which break qualified
    references (``t."S1.4"`` is mis-parsed by the TAP engine).  Avoid table
    aliases entirely — use the bare table name in FROM and unqualified column
    names in SELECT/WHERE.
    """
    ra_col, dec_col = cat["ra_col"], cat["dec_col"]
    cones = " OR ".join(
        f"CONTAINS(POINT('ICRS',{ra_col},{dec_col}),"
        f"CIRCLE('ICRS',{ra:.7f},{dec:.7f},{MATCH_RADIUS_DEG:.8f}))=1"
        for ra, dec in zip(ra_arr, dec_arr)
    )
    flux_col = cat["flux_col"]
    return (
        f"SELECT {ra_col} AS ra_radio, "
        f"{dec_col} AS dec_radio, "
        f"{flux_col} AS flux_mJy "
        f"FROM {cat['table']} "
        f"WHERE {cones}"
    )


def _query_with_retry(tap, adql: str, label: str) -> pd.DataFrame | None:
    """Execute a TAP query with retry logic. Returns DataFrame or None."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            job = tap.launch_job(adql)
            result = job.get_results()
            df = result.to_pandas()
            # Decode bytes columns (VOTable char type)
            for col in df.columns:
                if df[col].dtype == object:
                    sample = df[col].dropna()
                    if len(sample) and isinstance(sample.iloc[0], (bytes, bytearray)):
                        df[col] = df[col].str.decode("utf-8", errors="replace")
            return df
        except Exception as exc:
            sleep = RETRY_BASE_SLEEP * (2 ** (attempt - 1))
            log.warning("%s attempt %d/%d failed: %s -- retrying in %.0fs",
                        label, attempt, MAX_RETRIES, exc, sleep)
            if attempt < MAX_RETRIES:
                time.sleep(sleep)
    log.error("%s failed after %d attempts", label, MAX_RETRIES)
    return None


def _match_nearest(star_ra, star_dec, star_ids, radio_df):
    """For each star, find the nearest radio source within the match radius.

    Returns a list of dicts for stars that have at least one match.
    """
    if radio_df is None or len(radio_df) == 0:
        return []
    r_ra = radio_df["ra_radio"].to_numpy(float)
    r_dec = radio_df["dec_radio"].to_numpy(float)
    r_flux = radio_df["flux_mJy"].to_numpy(float)

    matches = []
    for i in range(len(star_ra)):
        cos_dec = np.cos(np.radians(star_dec[i]))
        sep = np.hypot((r_ra - star_ra[i]) * cos_dec, r_dec - star_dec[i])
        sep_arcsec = sep * 3600.0
        within = sep_arcsec <= MATCH_RADIUS_ARCSEC
        if within.any():
            k = int(np.argmin(sep_arcsec))
            matches.append({
                "source_id": star_ids[i],
                "ra_star": star_ra[i],
                "dec_star": star_dec[i],
                "ra_radio": float(r_ra[k]),
                "dec_radio": float(r_dec[k]),
                "sep_arcsec": float(sep_arcsec[k]),
                "flux_mJy": float(r_flux[k]),
            })
    return matches


def cross_match_catalog(stars: pd.DataFrame, cat_name: str, cat: dict):
    """Cross-match a set of stars against one radio catalog via VizieR TAP.

    Parameters
    ----------
    stars : DataFrame with columns ra, dec, source_id
    cat_name : human label (NVSS, VLASS, FIRST)
    cat : entry from RADIO_CATALOGS

    Returns
    -------
    list of match dicts, or empty list on total failure.
    """
    try:
        tap = _vizier_tap()
    except Exception as exc:
        log.error("cannot connect to VizieR TAP: %s", exc)
        return []

    ra = stars["ra"].to_numpy(float)
    dec = stars["dec"].to_numpy(float)
    ids = stars["source_id"].to_numpy()
    n = len(stars)

    all_matches = []
    n_failed_chunks = 0

    for i0 in range(0, n, CHUNK_SIZE):
        i1 = min(i0 + CHUNK_SIZE, n)
        label = f"{cat_name} chunk {i0}-{i1}"

        adql = _cone_union_adql(cat, ra[i0:i1], dec[i0:i1])
        radio_df = _query_with_retry(tap, adql, label)

        if radio_df is None:
            n_failed_chunks += 1
            continue

        chunk_matches = _match_nearest(
            ra[i0:i1], dec[i0:i1], ids[i0:i1], radio_df)
        all_matches.extend(chunk_matches)

        log.info("  %s: %d matches in %d stars",
                 label, len(chunk_matches), i1 - i0)

    if n_failed_chunks > 0:
        log.warning("%s: %d/%d chunks failed",
                    cat_name, n_failed_chunks,
                    (n + CHUNK_SIZE - 1) // CHUNK_SIZE)

    return all_matches


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Search C: radio counterparts to optically-dimmed stars")
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--max-candidates", type=int, default=MAX_CANDIDATES)
    ap.add_argument("--control-n", type=int, default=CONTROL_N)
    args = ap.parse_args()

    tag = args.tag

    # --- load residuals ---------------------------------------------------
    log.info("loading %s_resid.parquet", tag)
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{tag}_resid.parquet")
    n_all = len(d)
    log.info("loaded %d stars", n_all)

    residual = d["residual"].to_numpy(float)
    sigma_meas = d["sigma_meas"].to_numpy(float)
    ruwe = d["ruwe"].to_numpy(float)
    cstar = d["cstar_nsigma"].to_numpy(float)

    # --- select optically-dimmed candidates -------------------------------
    # Positive residual = star is fainter in G than expected from Ks, i.e.
    # something is intercepting optical light.
    # Require: residual > 3 * sigma_meas AND clean photometry.
    significance = np.where(sigma_meas > 0, residual / sigma_meas, 0.0)
    clean = (ruwe < 1.4) & (np.abs(cstar) < 3.0)
    dimmed = (residual > 0) & (significance > 3.0) & clean

    n_dimmed = int(dimmed.sum())
    log.info("optically dimmed (resid > 3 sigma_meas, clean): %d (%.3f%%)",
             n_dimmed, 100 * n_dimmed / n_all)

    # Take the top candidates sorted by significance
    d_dimmed = d[dimmed].copy()
    d_dimmed["significance"] = significance[dimmed]
    d_dimmed = (d_dimmed
                .nlargest(args.max_candidates, "significance")
                .reset_index(drop=True))
    n_cand = len(d_dimmed)
    log.info("top %d candidates selected for cross-match "
             "(significance range %.1f -- %.1f)",
             n_cand,
             d_dimmed["significance"].min(),
             d_dimmed["significance"].max())

    # --- build control sample ---------------------------------------------
    # Stars with similar brightness (G mag) but residual within 1 sigma of
    # zero -- no optical deficit.
    control_mask = (np.abs(significance) < 1.0) & clean
    d_control_pool = d[control_mask].copy()

    # Match the G-mag distribution of the candidates roughly
    g_min = d_dimmed["phot_g_mean_mag"].quantile(0.01)
    g_max = d_dimmed["phot_g_mean_mag"].quantile(0.99)
    d_control_pool = d_control_pool[
        (d_control_pool["phot_g_mean_mag"] >= g_min)
        & (d_control_pool["phot_g_mean_mag"] <= g_max)
    ]

    rng = np.random.default_rng(RNG_SEED)
    n_ctrl = min(args.control_n, len(d_control_pool))
    if n_ctrl < 100:
        log.warning("control pool too small (%d); results will be unreliable",
                    n_ctrl)
    ctrl_idx = rng.choice(len(d_control_pool), size=n_ctrl, replace=False)
    d_control = d_control_pool.iloc[ctrl_idx].reset_index(drop=True)
    log.info("control sample: %d stars (G in [%.1f, %.1f])",
             n_ctrl, g_min, g_max)

    # --- cross-match against radio catalogs -------------------------------
    results = {}
    all_candidate_matches = []
    vizier_reachable = True

    for cat_name, cat in RADIO_CATALOGS.items():
        log.info("=== %s (%s, %.1f GHz) ===", cat_name, cat["cite"],
                 cat["freq_ghz"])

        # Candidates
        log.info("querying %d candidates...", n_cand)
        cand_matches = cross_match_catalog(d_dimmed, cat_name, cat)

        if not cand_matches and n_cand > 0:
            # Could be genuinely zero, or could be VizieR unreachable.
            # Test with a tiny query.
            try:
                tap = _vizier_tap()
                test_q = (f"SELECT TOP 1 {cat['ra_col']} "
                          f"FROM {cat['table']}")
                test_r = _query_with_retry(tap, test_q, f"{cat_name}_health")
                if test_r is None:
                    log.error("%s: VizieR unreachable, skipping", cat_name)
                    vizier_reachable = False
                    results[cat_name] = {"error": "VizieR unreachable"}
                    continue
            except Exception:
                vizier_reachable = False
                results[cat_name] = {"error": "VizieR unreachable"}
                continue

        # Control
        log.info("querying %d control stars...", n_ctrl)
        ctrl_matches = cross_match_catalog(d_control, cat_name, cat)

        cand_rate = len(cand_matches) / max(n_cand, 1)
        ctrl_rate = len(ctrl_matches) / max(n_ctrl, 1)

        # Annotate candidate matches with catalog name
        for m in cand_matches:
            m["catalog"] = cat_name
            m["freq_ghz"] = cat["freq_ghz"]
        all_candidate_matches.extend(cand_matches)

        excess_factor = (cand_rate / ctrl_rate) if ctrl_rate > 0 else np.inf
        # Poisson uncertainty on the control rate
        ctrl_rate_err = np.sqrt(max(len(ctrl_matches), 1)) / max(n_ctrl, 1)

        results[cat_name] = {
            "freq_ghz": cat["freq_ghz"],
            "cite": cat["cite"],
            "n_candidates": n_cand,
            "n_candidate_matches": len(cand_matches),
            "candidate_match_rate": float(cand_rate),
            "n_control": n_ctrl,
            "n_control_matches": len(ctrl_matches),
            "control_match_rate": float(ctrl_rate),
            "control_match_rate_err": float(ctrl_rate_err),
            "excess_factor": float(excess_factor) if np.isfinite(excess_factor) else None,
        }

        log.info("  candidates: %d/%d = %.4f",
                 len(cand_matches), n_cand, cand_rate)
        log.info("  control   : %d/%d = %.4f",
                 len(ctrl_matches), n_ctrl, ctrl_rate)
        log.info("  excess factor: %s",
                 f"{excess_factor:.2f}" if np.isfinite(excess_factor) else "inf")

    # --- aggregate across surveys ----------------------------------------
    cand_ids_with_radio = set()
    for m in all_candidate_matches:
        cand_ids_with_radio.add(m["source_id"])

    n_unique_radio = len(cand_ids_with_radio)
    n_multi_survey = sum(
        1 for sid in cand_ids_with_radio
        if sum(1 for m in all_candidate_matches if m["source_id"] == sid) > 1
    )

    log.info("\n=== summary ===")
    log.info("unique candidates with any radio counterpart: %d / %d",
             n_unique_radio, n_cand)
    log.info("candidates detected in >1 survey: %d", n_multi_survey)

    # --- interpretation ---------------------------------------------------
    total_cand_matches = sum(
        r.get("n_candidate_matches", 0) for r in results.values()
        if isinstance(r, dict) and "n_candidate_matches" in r
    )
    total_ctrl_matches = sum(
        r.get("n_control_matches", 0) for r in results.values()
        if isinstance(r, dict) and "n_control_matches" in r
    )

    if total_ctrl_matches > 0:
        overall_excess = (total_cand_matches / max(n_cand, 1)) / \
                         (total_ctrl_matches / max(n_ctrl, 1))
    else:
        overall_excess = float("inf") if total_cand_matches > 0 else 0.0

    if not vizier_reachable:
        interpretation = (
            "VizieR TAP was unreachable for at least one catalog. Results "
            "are partial. Re-run when the service is available.")
    elif total_cand_matches == 0:
        interpretation = (
            "No radio counterparts found for any optically-dimmed candidate. "
            "This is the expected result if the dimming is caused by a "
            "spectrally selective absorber that does not emit at GHz "
            "frequencies, or if the dimming is a photometric/extinction "
            "artefact. The non-detection places no constraint on the "
            "harvesting hypothesis but rules out strong GHz emission "
            "correlated with the optical deficit.")
    elif overall_excess > 3.0:
        interpretation = (
            f"Radio detection rate among dimmed stars exceeds the control "
            f"rate by {overall_excess:.1f}x. This excess warrants follow-up: "
            f"individual matches should be checked against SIMBAD for known "
            f"active stars, AGN, or radio galaxies before any anomalous "
            f"interpretation is entertained.")
    elif overall_excess > 1.5:
        interpretation = (
            f"Mild excess ({overall_excess:.1f}x) in radio detections among "
            f"dimmed candidates relative to control. Likely explained by "
            f"residual active-star contamination or chance; not compelling.")
    else:
        interpretation = (
            f"Radio detection rate among dimmed stars ({overall_excess:.2f}x "
            f"the control rate) is consistent with chance alignments. No "
            f"evidence of anomalous radio emission correlated with the "
            f"optical deficit.")

    print(f"\n{'='*60}")
    print(f"SEARCH C: RADIO COUNTERPARTS")
    print(f"{'='*60}")
    for cat_name, r in results.items():
        if "error" in r:
            print(f"  {cat_name:8s}: {r['error']}")
        else:
            print(f"  {cat_name:8s}: {r['n_candidate_matches']:3d}/{n_cand} "
                  f"dimmed   {r['n_control_matches']:3d}/{n_ctrl} control   "
                  f"excess={r.get('excess_factor', 'N/A')}")
    print(f"\n  unique dimmed stars with radio: {n_unique_radio}")
    print(f"  detected in >1 survey:         {n_multi_survey}")
    print(f"\nINTERPRETATION: {interpretation}")

    # --- save outputs -----------------------------------------------------
    # JSON summary
    summary = {
        "tag": tag,
        "match_radius_arcsec": MATCH_RADIUS_ARCSEC,
        "n_all_stars": n_all,
        "n_optically_dimmed": n_dimmed,
        "n_candidates_queried": n_cand,
        "n_control": n_ctrl,
        "catalogs": results,
        "n_unique_radio_matches": n_unique_radio,
        "n_multi_survey_detections": n_multi_survey,
        "overall_excess_factor": (
            float(overall_excess) if np.isfinite(overall_excess) else None),
        "vizier_reachable": vizier_reachable,
        "interpretation": interpretation,
    }
    json_path = cfg.RESULT_DIR / f"searchC_radio_{tag}.json"
    json_path.write_text(json.dumps(summary, indent=2))
    log.info("wrote %s", json_path)

    # Candidates CSV
    if all_candidate_matches:
        match_df = pd.DataFrame(all_candidate_matches)
        # Merge back candidate properties for context
        merge_cols = ["source_id", "ra", "dec", "l", "b",
                      "phot_g_mean_mag", "M_Ks", "residual", "sigma_meas"]
        merge_cols = [c for c in merge_cols if c in d_dimmed.columns]
        match_df = match_df.merge(
            d_dimmed[merge_cols].drop_duplicates(subset="source_id"),
            on="source_id", how="left", suffixes=("", "_star"))
        csv_path = cfg.RESULT_DIR / f"searchC_radio_candidates_{tag}.csv"
        match_df.to_csv(csv_path, index=False)
        log.info("wrote %d radio-matched candidates to %s",
                 len(match_df), csv_path)
    else:
        log.info("no radio matches; no candidates CSV written")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
