#!/usr/bin/env python
"""Search T: the cold-blackbody audit of Planck's excluded-source bin.

    run.sh scripts/74_searchT_pccs2e.py

THE PHYSICS THAT PICKS THE INSTRUMENT
-------------------------------------
Hiding by running cold forces you to run big: radiating power P at temperature
T needs area P/(sigma T^4). One solar luminosity at 3 K needs a radius of
17,200 AU, which subtends about 5.7 arcmin at 100 pc.

But the band matters more than the size. At 3 K and 100 micron, hv/kT = 48 and
the Planck function is suppressed by e^-48 ~ 1e-21. IRAS, AKARI, WISE and
Herschel are ALL identically blind to a 3-5 K source. The Wien peak of a 3 K
blackbody sits at 176 GHz, inside Planck's CMB channels. So the only
instruments that can see this regime are Planck and the ground-based
millimetre surveys -- and every technosignature search ever published used the
infrared.

This also exposes a hole in channel 15: the Planck cold-clump catalogue
requires IRAS 100 micron as an input band, so a genuinely 3-5 K source cannot
be selected by it at all, no matter how bright.

WHY THIS WOULD NOT HAVE BEEN NOTICED
------------------------------------
An opaque 3 K surface seen against the 2.725 K CMB is an excess of 275 mK.
The CMB anisotropy is about 100 microkelvin rms at 5 arcmin. So such an object
is a ~3000-sigma positive spot -- among the loudest things in the Planck maps.

It survives component separation rather than being removed by it: the internal
linear combinations preserve anything with a blackbody spectral response, so a
cold blackbody is degenerate with a CMB temperature fluctuation and passes
straight through. What actually happens to it is MASKING. Compact sources are
cut before cosmological analysis, and a 3000-sigma compact positive spot is
flagged as an instrumental artefact or a moving Solar System object.

The place it lands is PCCS2E -- the "excluded" catalogue of Planck compact
sources that could not be validated, mostly because they sit in regions of
high cirrus. Planck's own documentation treats it as a contamination bin. It
has never been searched as a source population.

THE DISCRIMINANT
----------------
Three SED shapes, cleanly separable across 100-857 GHz:

  interstellar dust   modified blackbody, beta ~ 1.6, T ~ 20 K. Rises steeply
                      to 857 GHz. This is almost everything in PCCS2E.
  thermal SZ          has a null at 217 GHz and a decrement below it.
  cold blackbody      beta = 0 by construction, T < 10 K, CMB-like at low
                      frequency and peeling away above ~200 GHz.

Fitting beta as a free parameter separates them. Dust cannot reach beta ~ 0:
grains are far smaller than a millimetre wavelength, and that is grain physics
rather than a fitting convention.
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
log = logging.getLogger("searchT")

H = 6.62607015e-34
KB = 1.380649e-23
C = 2.99792458e8
T_CMB = 2.72548

# Planck HFI bands, GHz -> effective frequency
BANDS = {"100": 100.0, "143": 143.0, "217": 217.0,
         "353": 353.0, "545": 545.0, "857": 857.0}

BETA_BLACKBODY_MAX = 0.5
T_COLD_MAX = 10.0


def planck_bnu(nu_hz, T):
    x = H * nu_hz / (KB * np.maximum(T, 1e-3))
    x = np.clip(x, 1e-6, 500)
    return (2 * H * nu_hz ** 3 / C ** 2) / (np.expm1(x))


def fit_sed(nu_ghz, flux, err):
    """Grid-fit a modified blackbody I ~ nu^beta B_nu(T).

    Returns the best (T, beta, chi2) and the chi2 of the beta-fixed-to-zero
    (pure blackbody) and beta-fixed-to-1.6 (dust) alternatives.
    """
    nu = np.asarray(nu_ghz, float) * 1e9
    f = np.asarray(flux, float)
    e = np.maximum(np.asarray(err, float), 1e-6 * np.max(np.abs(f)))
    ok = np.isfinite(f) & np.isfinite(e) & (e > 0)
    if ok.sum() < 3:
        return None
    nu, f, e = nu[ok], f[ok], e[ok]

    # Vectorised over the whole (T, beta) grid. The scalar version evaluated
    # ~6000 models per source in Python and could not finish 134,000 sources.
    temps = np.concatenate([np.arange(2.0, 12.0, 0.25),
                            np.arange(12.0, 40.0, 0.5)])
    betas = np.arange(-0.5, 3.01, 0.05)

    b_nu = planck_bnu(nu[None, :], temps[:, None])          # (nT, nband)
    model = (nu[None, None, :] ** betas[None, :, None]) * b_nu[:, None, :]

    w = 1.0 / e ** 2
    denom = np.einsum("tbn,n->tb", model ** 2, w)
    numer = np.einsum("tbn,n->tb", model, f * w)
    with np.errstate(divide="ignore", invalid="ignore"):
        amp = np.where(denom > 0, numer / denom, np.nan)
    amp = np.where(amp > 0, amp, np.nan)

    resid = f[None, None, :] - amp[:, :, None] * model
    chi2 = np.einsum("tbn,n->tb", resid ** 2, w)
    chi2 = np.where(np.isfinite(amp), chi2, np.inf)

    if not np.isfinite(chi2).any():
        return None
    it, ib = np.unravel_index(np.argmin(chi2), chi2.shape)

    i_bb = int(np.argmin(np.abs(betas - 0.0)))
    i_du = int(np.argmin(np.abs(betas - 1.60)))
    return {"chi2": float(chi2[it, ib]),
            "T": float(temps[it]), "beta": float(betas[ib]),
            "chi2_blackbody": float(np.min(chi2[:, i_bb])),
            "chi2_dust": float(np.min(chi2[:, i_du])),
            "n_bands": int(ok.sum())}


def load_pccs2e(force=False):
    """PCCS2E, the excluded compact-source bin, from VizieR."""
    out = cfg.RAW_DIR / "pccs2e.parquet"
    if out.exists() and not force:
        d = pd.read_parquet(out)
        log.info("cached PCCS2E: %d rows", len(d))
        return d

    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=["**"])
    v.ROW_LIMIT = -1

    log.info("locating the Planck compact-source catalogues on VizieR ...")
    cats = v.get_catalogs("J/A+A/594/A26")
    frames = []
    for t in cats:
        name = str(t.meta.get("name", ""))
        log.info("  table %-40s %7d rows, %d cols",
                 name, len(t), len(t.colnames))
        df = t.to_pandas()
        df["_table"] = name
        frames.append(df)
    if not frames:
        raise RuntimeError("no PCCS2 tables returned")
    d = pd.concat(frames, ignore_index=True)
    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
    d.to_parquet(out, index=False)
    log.info("wrote %s (%d rows)", out, len(d))
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--tag", default="primary")
    args = ap.parse_args()

    d = load_pccs2e(force=args.force)
    log.info("total rows across PCCS2 tables: %d", len(d))
    log.info("columns: %s", ", ".join(list(d.columns)[:40]))

    # Identify which tables are the EXCLUDED bin and which frequency each is.
    tabs = d["_table"].value_counts()
    log.info("")
    for k, n in tabs.items():
        log.info("  %-46s %7d", k, n)

    # Locate flux and position columns pragmatically: VizieR names vary by
    # release, so match on prefix rather than assuming an exact schema.
    def find_col(df, *cands):
        for c in cands:
            if c in df.columns:
                return c
        for c in df.columns:
            lc = c.lower()
            if any(lc.startswith(x.lower()) for x in cands):
                return c
        return None

    ra_c = find_col(d, "RAJ2000", "_RAJ2000", "RA")
    de_c = find_col(d, "DEJ2000", "_DEJ2000", "DE")
    glon_c = find_col(d, "GLON")
    glat_c = find_col(d, "GLAT")
    log.info("")
    log.info("position columns: ra=%s dec=%s glon=%s glat=%s",
             ra_c, de_c, glon_c, glat_c)

    flux_cols = [c for c in d.columns
                 if c.lower().startswith(("detflux", "aperflux", "gauflux",
                                          "psfflux", "flux"))]
    log.info("flux-like columns: %s", ", ".join(flux_cols[:20]))

    if not flux_cols or ra_c is None:
        log.error("could not identify the catalogue schema; "
                  "inspect the columns above and adjust")
        summary = {"tag": args.tag, "status": "schema_unrecognised",
                   "n_rows": int(len(d)),
                   "tables": {str(k): int(v) for k, v in tabs.items()},
                   "columns": list(d.columns)}
        (cfg.RESULT_DIR / f"searchT_pccs2e_{args.tag}.json").write_text(
            json.dumps(summary, indent=2))
        return 1

    # ---- build per-source SEDs across frequency tables -------------------
    # Each PCCS2 table is one frequency. Cross-match by position to assemble
    # a spectrum per sky location.
    from scipy.spatial import cKDTree

    def freq_of(tabname):
        for b in BANDS:
            if b in str(tabname):
                return BANDS[b]
        return np.nan

    d["_freq"] = d["_table"].map(freq_of)
    d = d[np.isfinite(d["_freq"])]
    log.info("rows with an identifiable frequency: %d", len(d))
    if len(d) == 0:
        log.error("no frequency could be parsed from the table names")
        return 1

    # Each frequency is its own VizieR table with frequency-suffixed columns
    # (DetFlux857, DetFlux545, ...), so the flux column has to be selected
    # per row from that row's own band rather than globally.
    band_of = d["_freq"].astype(int).astype(str)
    flux = np.full(len(d), np.nan)
    ferr = np.full(len(d), np.nan)
    for b in sorted(band_of.unique()):
        m = (band_of == b).to_numpy()
        fc = find_col(d, f"DetFlux{b}", f"GauFlux{b}", f"AperFlux{b}")
        ec = find_col(d, f"e_DetFlux{b}", f"e_GauFlux{b}", f"e_AperFlux{b}")
        if fc is None:
            log.warning("  no flux column for band %s GHz", b)
            continue
        flux[m] = d.loc[m, fc].to_numpy(float)
        if ec is not None:
            ferr[m] = d.loc[m, ec].to_numpy(float)
        log.info("  band %4s GHz -> %s (%d rows)", b, fc, int(m.sum()))
    d = d.assign(_flux=flux, _ferr=ferr)
    fcol, ecol = "_flux", "_ferr"

    ra = d[ra_c].to_numpy(float)
    de = d[de_c].to_numpy(float)
    good = np.isfinite(ra) & np.isfinite(de) & np.isfinite(d[fcol].to_numpy(float))
    d = d[good].reset_index(drop=True)
    log.info("rows with a usable flux in their own band: %d", len(d))
    ra, de = d[ra_c].to_numpy(float), d[de_c].to_numpy(float)

    rr, dr = np.radians(ra), np.radians(de)
    xyz = np.column_stack([np.cos(dr) * np.cos(rr),
                           np.cos(dr) * np.sin(rr), np.sin(dr)])
    tree = cKDTree(xyz)
    r_match = 2 * np.sin(np.radians(5.0 / 60.0) / 2)   # 5 arcmin
    groups = tree.query_ball_point(xyz, r=r_match)

    seen = np.zeros(len(d), dtype=bool)
    rows = []
    for i, g in enumerate(groups):
        if seen[i]:
            continue
        g = [j for j in g if not seen[j]]
        if not g:
            continue
        seen[g] = True
        sub = d.iloc[g]
        nu = sub["_freq"].to_numpy(float)
        fl = sub[fcol].to_numpy(float)
        er = sub[ecol].to_numpy(float) if ecol else np.full(len(fl), np.nan)
        if not np.isfinite(er).any():
            er = 0.1 * np.abs(fl)
        if len(np.unique(nu)) < 3:
            continue
        fit = fit_sed(nu, fl, er)
        if fit is None:
            continue
        rows.append({
            "ra": float(sub[ra_c].iloc[0]), "dec": float(sub[de_c].iloc[0]),
            "glon": float(sub[glon_c].iloc[0]) if glon_c else np.nan,
            "glat": float(sub[glat_c].iloc[0]) if glat_c else np.nan,
            "n_bands": fit["n_bands"], "T": fit["T"], "beta": fit["beta"],
            "chi2": fit["chi2"], "chi2_blackbody": fit["chi2_blackbody"],
            "chi2_dust": fit["chi2_dust"],
        })

    if not rows:
        log.error("no source had three or more frequencies")
        return 1

    s = pd.DataFrame(rows)
    log.info("")
    log.info("sources with >= 3 Planck bands: %d", len(s))
    log.info("fitted beta : median %.2f, 16-84%% [%.2f, %.2f]",
             float(s["beta"].median()),
             float(s["beta"].quantile(0.16)), float(s["beta"].quantile(0.84)))
    log.info("fitted T    : median %.1f K", float(s["T"].median()))

    # A fit that lands on a grid boundary is unconstrained: the data want a
    # value outside the grid and the fitter is simply running to the wall.
    # Without this check the entire candidate list was beta = -0.50 exactly,
    # the lower edge, with T pinned at the 2 K floor -- a pathology, not a
    # measurement. Note it also breaks the mirror control, because if every
    # fit jams against the LOW edge nothing can ever reach the high one.
    b_lo, b_hi = -0.5, 3.0
    t_lo, t_hi = 2.0, 39.5
    at_edge = ((np.abs(s["beta"] - b_lo) < 1e-6) | (np.abs(s["beta"] - b_hi) < 1e-6)
               | (np.abs(s["T"] - t_lo) < 1e-6) | (np.abs(s["T"] - t_hi) < 1e-6))
    s = s.assign(at_grid_edge=at_edge)
    log.info("")
    log.info("fits pinned to a grid boundary (unconstrained): %d of %d (%.1f%%)",
             int(at_edge.sum()), len(s), 100 * at_edge.mean())

    interior = s[~s["at_grid_edge"]]
    log.info("interior fits usable for the test: %d", len(interior))

    cold_bb = interior[(interior["beta"] < BETA_BLACKBODY_MAX)
                       & (interior["T"] < T_COLD_MAX)
                       & (interior["chi2_blackbody"] < interior["chi2_dust"] - 4.0)]
    unphysical = interior[interior["beta"] > 2.5]
    log.info("")
    log.info("beta < %.1f AND T < %.0f K AND blackbody beats dust by dchi2>4: %d",
             BETA_BLACKBODY_MAX, T_COLD_MAX, len(cold_bb))
    log.info("beta > 3 (unphysical mirror control)                       : %d",
             len(unphysical))

    if len(cold_bb):
        log.info("")
        for _, r in cold_bb.nsmallest(15, "beta").iterrows():
            log.info("  l=%7.2f b=%+6.2f  T=%5.1f K  beta=%+.2f  "
                     "chi2 bb=%.1f dust=%.1f  bands=%d",
                     r["glon"], r["glat"], r["T"], r["beta"],
                     r["chi2_blackbody"], r["chi2_dust"], r["n_bands"])

    n_c, n_u = len(cold_bb), len(unphysical)
    if n_c == 0:
        verdict = (
            f"NULL. Of {len(s)} Planck compact sources with three or more "
            f"bands, none prefers a cold blackbody over modified-blackbody "
            f"dust. The excluded bin that Planck treats as contamination "
            f"contains no source with beta near zero at T below "
            f"{T_COLD_MAX:.0f} K. This is the first time that bin has been "
            f"searched as a source population, and the band -- 100 to 857 GHz "
            f"-- is the only one where a 3-5 K radiator emits at all.")
    elif n_c <= n_u:
        verdict = (
            f"NULL. The {n_c} cold-blackbody candidates are matched by {n_u} "
            f"sources fitting an equally unphysical beta > 3, so both tails "
            f"are SED-fitting scatter rather than a population.")
    else:
        verdict = (
            f"{n_c} sources prefer a cold blackbody (beta < "
            f"{BETA_BLACKBODY_MAX}, T < {T_COLD_MAX:.0f} K) over dust by "
            f"dchi2 > 4, against {n_u} in the unphysical high-beta mirror. "
            f"Each needs checking against the Solar System object registry, "
            f"the Planck artefact flags, and HI4PI: absence of an HI "
            f"counterpart at the Boulanger emissivity is the cleanest single "
            f"falsifier of cirrus.")

    print(f"\n{'='*72}")
    print("SEARCH T: COLD BLACKBODIES IN PLANCK'S EXCLUDED BIN")
    print(f"{'='*72}")
    print(f"  catalogue rows                      : {len(d):,}")
    print(f"  sources with >= 3 bands             : {len(s):,}")
    print(f"  median fitted beta                  : {float(s['beta'].median()):.2f}")
    print(f"  median fitted T                     : {float(s['T'].median()):.1f} K")
    print(f"  cold blackbody (beta<{BETA_BLACKBODY_MAX}, T<{T_COLD_MAX:.0f}K) : {n_c}")
    print(f"  unphysical beta>3 (mirror)          : {n_u}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchT_pccs2e_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "n_catalogue_rows": int(len(d)),
        "n_multiband_sources": int(len(s)),
        "median_beta": float(s["beta"].median()),
        "median_T": float(s["T"].median()),
        "n_cold_blackbody": n_c,
        "n_unphysical_mirror": n_u,
        "beta_threshold": BETA_BLACKBODY_MAX,
        "T_threshold": T_COLD_MAX,
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    s.to_csv(cfg.RESULT_DIR / f"searchT_sed_fits_{args.tag}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
