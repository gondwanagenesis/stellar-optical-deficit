#!/usr/bin/env python
"""Search T, corrected: cold blackbodies in Planck's excluded compact-source bin.

    run.sh scripts/78_searchT_pccs2e_fixed.py --tag primary

THE PHYSICS THAT PICKS THE INSTRUMENT
-------------------------------------
Hiding by running cold forces you to run big: radiating power P at temperature
T needs area P/(sigma T^4). One solar luminosity at 3 K needs a radius of
17,200 AU, about 5.7 arcmin at 100 pc.

The band matters more than the size. At 3 K and 100 micron, hv/kT = 48 and the
Planck function is suppressed by e^-48 ~ 1e-21. IRAS, AKARI, WISE and Herschel
are ALL identically blind to a 3-5 K source. The Wien peak of a 3 K blackbody
sits at 176 GHz, inside Planck's CMB channels. So the only instruments that can
see this regime are Planck and the ground-based millimetre surveys, and every
technosignature search ever published used the infrared.

An opaque 3 K surface seen against the 2.725 K CMB is an excess of 275 mK
against ~100 microkelvin rms anisotropy: a ~3000-sigma compact positive spot.
It survives component separation rather than being removed by it, because an
internal linear combination preserves anything with a blackbody spectral
response. What happens to it is MASKING, and the place it lands is PCCS2E,
the "excluded" bin of Planck compact sources that could not be validated.
Planck's own documentation treats that bin as contamination. It has never been
searched as a source population.

WHY THIS SCRIPT EXISTS: THE MIRROR CONTROL IN SCRIPT 74 WAS DEAD
-----------------------------------------------------------------
Script 74 fitted beta on the grid np.arange(-0.5, 3.01, 0.05), whose largest
element is 2.9999999999999991, and then defined its unphysical mirror control
as beta > 3.0. That condition can never be satisfied. The mirror was
identically zero for any input whatsoever.

That is the most dangerous available failure mode for this project: the verdict
logic compares n_cold against n_mirror, so a candidate list of any size would
have been reported against a false-positive rate measured as exactly zero.

Three changes, all of which make the test harder to pass rather than easier:

1. The beta grid runs to 5.0, so the high-beta tail is measured rather than
   railed against the grid edge.

2. The mirror is made SYMMETRIC ABOUT DUST rather than pinned to an arbitrary
   number. Interstellar dust sits at beta_dust = 1.6. A cold blackbody is
   beta = 0, an offset of -1.6. The equal-and-opposite absurdity is therefore
   beta = 3.2, not "beta > 3". Both tails are cut at the same
   |beta - 1.6| >= 1.1 and carry the same T < 10 K requirement and the same
   dchi2 > 4 preference over the dust fit, so the mirror measures this
   pipeline's own rate of manufacturing an equally unphysical SED from the
   same sources, the same bands and the same noise.

   The mirror is meaningful here. Unlike 2MASS aperture mismatch, which can
   only make Ks brighter and so gives a one-sided contaminant where
   one-sidedness proves nothing, SED-fit scatter over three to six noisy bands
   has no reason to prefer one sign of beta. Both tails should populate
   equally under the null, and the dominant real contaminant (cirrus) lives at
   beta ~ 1.5-2.0, between them.

3. The fit is vectorised over the (T, beta) grid, which is what makes the
   extended grid affordable: script 74's pure-Python double loop takes tens of
   ms per source, this takes under 1 ms.

Two further controls that script 74 did not have:

  Galactic-latitude split. Cirrus is a disc population and concentrates at low
  |b|. An extrasolar population would not. If candidates track |b| the way the
  bulk catalogue does, they are cirrus.

  Permuted-flux recovery. The same fitter is run on the same sources with
  fluxes shuffled across position WITHIN each band. That is mark permutation
  on a scalar mark, which is legitimate (unlike the vector case that broke
  Search F), and it preserves the flux distribution and the sky geometry
  exactly while destroying any real SED.
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
log = logging.getLogger("searchT2")

H = 6.62607015e-34
KB = 1.380649e-23
C = 2.99792458e8

BANDS = {"100": 100.0, "143": 143.0, "217": 217.0,
         "353": 353.0, "545": 545.0, "857": 857.0}

BETA_DUST = 1.6          # interstellar dust emissivity index
BETA_OFFSET = 1.1        # how far from dust a fit must sit to count
T_COLD_MAX = 10.0        # applied to BOTH tails, so the mirror stays symmetric
DCHI2_MIN = 4.0

BETA_LO_MAX = BETA_DUST - BETA_OFFSET     # 0.5  -> cold-blackbody side
BETA_HI_MIN = BETA_DUST + BETA_OFFSET     # 2.7  -> unphysical mirror side

TEMPS = np.concatenate([np.arange(2.0, 12.0, 0.25),
                        np.arange(12.0, 40.0, 0.5)])
BETAS = np.arange(-0.5, 5.001, 0.05)      # reaches 5.0, mirror not railed


def planck_bnu(nu_hz, T):
    x = H * nu_hz / (KB * np.maximum(T, 1e-3))
    x = np.clip(x, 1e-6, 500.0)
    return (2 * H * nu_hz ** 3 / C ** 2) / np.expm1(x)


def fit_sed_vec(nu_ghz, flux, err):
    """Grid-fit I ~ nu^beta B_nu(T), vectorised over the whole (T, beta) grid.

    Amplitude is solved analytically at every grid point, so the search is only
    two-dimensional. Returns the global best fit plus the profile chi2 of the
    two fixed-beta hypotheses (pure blackbody, dust).
    """
    nu = np.asarray(nu_ghz, float) * 1e9
    f = np.asarray(flux, float)
    e = np.asarray(err, float)
    ok = np.isfinite(f) & np.isfinite(e) & (e > 0)
    if ok.sum() < 3:
        return None
    nu, f, e = nu[ok], f[ok], e[ok]

    b_nu = planck_bnu(nu[None, :], TEMPS[:, None])            # (nT, nband)
    nu_beta = nu[None, :] ** BETAS[:, None]                   # (nbeta, nband)
    model = b_nu[:, None, :] * nu_beta[None, :, :]            # (nT, nbeta, nband)

    w = 1.0 / e ** 2
    denom = np.einsum("tbk,k->tb", model ** 2, w)
    numer = np.einsum("tbk,k->tb", model, f * w)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(denom > 0, numer / denom, np.nan)
    a = np.where(a > 0, a, np.nan)                            # emission only

    resid = f[None, None, :] - a[:, :, None] * model
    with np.errstate(invalid="ignore"):
        chi2 = np.einsum("tbk,k->tb", resid ** 2, w)
    chi2 = np.where(np.isfinite(a), chi2, np.inf)
    if not np.isfinite(chi2).any():
        return None

    it, ib = np.unravel_index(np.argmin(chi2), chi2.shape)
    j_bb = int(np.argmin(np.abs(BETAS - 0.0)))
    j_dust = int(np.argmin(np.abs(BETAS - BETA_DUST)))
    return {"chi2": float(chi2[it, ib]),
            "T": float(TEMPS[it]), "beta": float(BETAS[ib]),
            "T_blackbody": float(TEMPS[int(np.argmin(chi2[:, j_bb]))]),
            "chi2_blackbody": float(np.min(chi2[:, j_bb])),
            "chi2_dust": float(np.min(chi2[:, j_dust])),
            "n_bands": int(ok.sum())}


def find_col(df, *cands):
    for c in cands:
        if c in df.columns:
            return c
    for c in df.columns:
        if any(c.lower().startswith(x.lower()) for x in cands):
            return c
    return None


def load_pccs2e(force=False):
    out = cfg.RAW_DIR / "pccs2e.parquet"
    if out.exists() and not force:
        d = pd.read_parquet(out)
        log.info("cached PCCS2/PCCS2E: %d rows", len(d))
        return d
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=["**"])
    v.ROW_LIMIT = -1
    frames = []
    for t in v.get_catalogs("J/A+A/594/A26"):
        df = t.to_pandas()
        df["_table"] = str(t.meta.get("name", ""))
        frames.append(df)
    d = pd.concat(frames, ignore_index=True)
    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
    d.to_parquet(out, index=False)
    return d


def assemble_seds(d):
    """Attach a frequency and a flux to every single-band Planck detection."""
    def freq_of(t):
        for b in BANDS:
            if b in str(t):
                return BANDS[b]
        return np.nan

    d = d.assign(_freq=d["_table"].map(freq_of))
    d = d[np.isfinite(d["_freq"])].reset_index(drop=True)
    d["_excluded"] = d["_table"].astype(str).str.endswith("e")

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
        log.info("  band %4s GHz -> %s (%d detections)", b, fc, int(m.sum()))
    d = d.assign(_flux=flux, _ferr=ferr)

    ra = d["RAJ2000"].to_numpy(float)
    de = d["DEJ2000"].to_numpy(float)
    d = d[np.isfinite(ra) & np.isfinite(de) & np.isfinite(d["_flux"])]
    return d.reset_index(drop=True)


def group_by_position(d, radius_arcmin=5.0):
    from scipy.spatial import cKDTree
    rr = np.radians(d["RAJ2000"].to_numpy(float))
    dr = np.radians(d["DEJ2000"].to_numpy(float))
    xyz = np.column_stack([np.cos(dr) * np.cos(rr),
                           np.cos(dr) * np.sin(rr), np.sin(dr)])
    tree = cKDTree(xyz)
    r = 2 * np.sin(np.radians(radius_arcmin / 60.0) / 2)
    nbr = tree.query_ball_point(xyz, r=r)
    seen = np.zeros(len(d), bool)
    groups = []
    for i, g in enumerate(nbr):
        if seen[i]:
            continue
        g = [j for j in g if not seen[j]]
        if not g:
            continue
        seen[g] = True
        groups.append(g)
    return groups


def run_fits(d, groups, flux_override=None):
    fl_all = (d["_flux"].to_numpy(float) if flux_override is None
              else flux_override)
    er_all = d["_ferr"].to_numpy(float)
    nu_all = d["_freq"].to_numpy(float)
    glon = d["GLON"].to_numpy(float)
    glat = d["GLAT"].to_numpy(float)
    ra = d["RAJ2000"].to_numpy(float)
    de = d["DEJ2000"].to_numpy(float)
    exc = d["_excluded"].to_numpy(bool)

    rows = []
    for g in groups:
        nu = nu_all[g]
        if len(np.unique(nu)) < 3:
            continue
        fl = fl_all[g]
        er = er_all[g]
        er = np.where(np.isfinite(er) & (er > 0), er, 0.1 * np.abs(fl))
        fit = fit_sed_vec(nu, fl, er)
        if fit is None:
            continue
        fit.update(ra=float(ra[g[0]]), dec=float(de[g[0]]),
                   glon=float(glon[g[0]]), glat=float(glat[g[0]]),
                   excluded=bool(exc[g].any()))
        rows.append(fit)
    return pd.DataFrame(rows)


def select(s):
    """Symmetric selection: same T cut, same dchi2 margin, mirrored beta.

    On the cold side the hypothesis is beta = 0 exactly, so the profile chi2 at
    beta = 0 is the right statistic. On the mirror side there is no privileged
    beta, so the profile minimum over the high-beta region is used, which is
    the more permissive of the two and therefore biases the mirror UP -- the
    conservative direction for a claim of excess.
    """
    if len(s) == 0:
        return s, s
    cold = s[(s["beta"] < BETA_LO_MAX) & (s["T"] < T_COLD_MAX)
             & (s["chi2_blackbody"] < s["chi2_dust"] - DCHI2_MIN)]
    mirror = s[(s["beta"] > BETA_HI_MIN) & (s["T"] < T_COLD_MAX)
               & (s["chi2"] < s["chi2_dust"] - DCHI2_MIN)]
    return cold, mirror


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--n-perm", type=int, default=3)
    args = ap.parse_args()
    rng = np.random.default_rng(78_2026)

    d = load_pccs2e(force=args.force)
    d = assemble_seds(d)
    log.info("single-band detections with a usable flux: %d", len(d))
    log.info("  of which in the EXCLUDED bin: %d", int(d["_excluded"].sum()))

    groups = group_by_position(d)
    log.info("position groups within 5 arcmin: %d", len(groups))

    s = run_fits(d, groups)
    log.info("sources with >= 3 Planck bands: %d", len(s))
    if len(s) == 0:
        log.error("no source had three or more bands: broken, not null")
        return 1

    log.info("fitted beta : median %.2f, 16-84%% [%.2f, %.2f], max %.2f",
             float(s["beta"].median()), float(s["beta"].quantile(0.16)),
             float(s["beta"].quantile(0.84)), float(s["beta"].max()))
    log.info("fitted T    : median %.1f K", float(s["T"].median()))

    cold, mirror = select(s)
    n_c, n_m = len(cold), len(mirror)
    n_dead = int((s["beta"] > 3.0).sum())
    log.info("")
    log.info("cold-blackbody tail  (beta < %.1f, T < %.0f K, dchi2 > %.0f): %d",
             BETA_LO_MAX, T_COLD_MAX, DCHI2_MIN, n_c)
    log.info("unphysical mirror    (beta > %.1f, same T and dchi2)       : %d",
             BETA_HI_MIN, n_m)
    log.info("sources with beta > 3.0 (unreachable on script 74's grid)  : %d",
             n_dead)

    # --- control 1: Galactic latitude. cirrus is a disc population ---------
    def lowb(x):
        return (float(np.mean(np.abs(x["glat"]) < 10.0))
                if len(x) else float("nan"))
    lat = {"all_sources": lowb(s), "cold_tail": lowb(cold),
           "mirror_tail": lowb(mirror)}
    log.info("")
    log.info("fraction at |b| < 10 deg -- all %.3f, cold %.3f, mirror %.3f",
             lat["all_sources"], lat["cold_tail"], lat["mirror_tail"])

    # --- control 2: permuted-flux recovery --------------------------------
    perm_cold, perm_mirror = [], []
    band = d["_freq"].to_numpy(float)
    base = d["_flux"].to_numpy(float)
    for i in range(args.n_perm):
        fp = base.copy()
        for b in np.unique(band):
            m = band == b
            fp[m] = rng.permutation(fp[m])
        cp, mp = select(run_fits(d, groups, flux_override=fp))
        perm_cold.append(int(len(cp)))
        perm_mirror.append(int(len(mp)))
        log.info("  permutation %d/%d -> cold %d, mirror %d",
                 i + 1, args.n_perm, len(cp), len(mp))
    pc = float(np.mean(perm_cold)) if perm_cold else float("nan")
    pm = float(np.mean(perm_mirror)) if perm_mirror else float("nan")

    if n_c == 0:
        verdict = (
            f"NULL. Of {len(s)} Planck compact sources with three or more "
            f"bands ({int(s['excluded'].sum())} of them from the excluded "
            f"bin), none prefers a cold blackbody over modified-blackbody dust "
            f"by dchi2 > {DCHI2_MIN:.0f} at T < {T_COLD_MAX:.0f} K. The "
            f"symmetric unphysical mirror at beta > {BETA_HI_MIN:.1f} returns "
            f"{n_m}, and permuted flux returns {pc:.1f} cold / {pm:.1f} "
            f"mirror, so the fitter's two-sided false-positive rate is "
            f"measured rather than assumed and the absence is not a selection "
            f"artefact. This is the first search of Planck's excluded "
            f"compact-source bin as a source population, in the only band "
            f"where a 3-5 K radiator emits at all.")
    elif n_c <= max(n_m, pc):
        verdict = (
            f"NULL. The {n_c} cold-blackbody candidates are matched by {n_m} "
            f"in the symmetric unphysical mirror and {pc:.1f} under flux "
            f"permutation, so both tails are SED-fitting scatter rather than "
            f"a population.")
    else:
        verdict = (
            f"{n_c} sources prefer a cold blackbody against {n_m} in the "
            f"symmetric mirror and {pc:.1f} expected from permuted flux; "
            f"excess over mirror {n_c - n_m}. Every survivor needs checking "
            f"against the Solar System object registry, the Planck artefact "
            f"flags, and HI4PI, where absence of an HI counterpart at the "
            f"Boulanger emissivity is the cleanest single falsifier of "
            f"cirrus. Low-|b| fraction is {lat['cold_tail']:.2f} against "
            f"{lat['all_sources']:.2f} for the full catalogue; if those agree "
            f"the candidates are a disc population and therefore cirrus.")

    print(f"\n{'=' * 74}")
    print("SEARCH T (CORRECTED): COLD BLACKBODIES IN PLANCK'S EXCLUDED BIN")
    print(f"{'=' * 74}")
    print(f"  single-band detections              : {len(d):,}")
    print(f"  sources with >= 3 bands             : {len(s):,}")
    print(f"  median fitted beta / T              : "
          f"{float(s['beta'].median()):.2f} / {float(s['T'].median()):.1f} K")
    print(f"  cold blackbody  (beta < {BETA_LO_MAX:.1f})       : {n_c}")
    print(f"  unphysical mirror (beta > {BETA_HI_MIN:.1f})     : {n_m}")
    print(f"  script 74's dead mirror (beta > 3.0): {n_dead} "
          f"(unreachable on its grid, capped at 2.99999)")
    print(f"  permuted-flux recovery cold/mirror  : {pc:.1f} / {pm:.1f}")
    print(f"  |b|<10 fraction all/cold/mirror     : "
          f"{lat['all_sources']:.2f} / {lat['cold_tail']:.2f} / "
          f"{lat['mirror_tail']:.2f}")
    print(f"\nVERDICT: {verdict}")

    out = cfg.RESULT_DIR / f"searchT_pccs2e_{args.tag}.json"
    out.write_text(json.dumps({
        "tag": args.tag,
        "n_detections": int(len(d)),
        "n_excluded_bin_detections": int(d["_excluded"].sum()),
        "n_multiband_sources": int(len(s)),
        "n_multiband_from_excluded_bin": int(s["excluded"].sum()),
        "median_beta": float(s["beta"].median()),
        "median_T": float(s["T"].median()),
        "max_beta": float(s["beta"].max()),
        "beta_dust": BETA_DUST, "beta_offset": BETA_OFFSET,
        "beta_lo_max": BETA_LO_MAX, "beta_hi_min": BETA_HI_MIN,
        "T_threshold": T_COLD_MAX, "dchi2_min": DCHI2_MIN,
        "n_cold_blackbody": int(n_c),
        "n_unphysical_mirror": int(n_m),
        "n_beta_gt_3_unreachable_on_script74_grid": n_dead,
        "permuted_flux_cold": perm_cold,
        "permuted_flux_mirror": perm_mirror,
        "lowb_fraction": lat,
        "verdict": verdict,
    }, indent=2))
    log.info("wrote %s", out)
    s.to_csv(cfg.RESULT_DIR / f"searchT_sed_fits_{args.tag}.csv", index=False)
    if n_c:
        cold.to_csv(cfg.RESULT_DIR / f"searchT_cold_candidates_{args.tag}.csv",
                    index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
