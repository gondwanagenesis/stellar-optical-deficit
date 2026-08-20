#!/usr/bin/env python
"""Search V (channel 22): does sub-PSF structure point along the proper motion?

    run.sh scripts/85_searchV_harmonic_phase.py --tag primary

THE IDEA
--------
Every technosignature search that works from photometric excess or deficit
has the same dominant false positive: a chance-aligned background source
inside the point spread function, adding light the star did not emit or
distorting the astrometry so the star looks disturbed. It is the class that
destroyed the Project Hephaistos Dyson-sphere candidates. The usual defences
are statistical -- crowding priors, blend flags, image-quality cuts -- and
they are all indirect, because none of them measures the *direction* of the
contaminating structure.

Gaia DR3 does. For every source the Image Parameter Determination fits a PSF
to each windowed transit, and the goodness of that fit varies with the
position angle of the scan direction whenever there is structure below the
PSF: worst when the scan runs along the elongation, best when it runs across.
Gaia fits that variation as a sinusoid in twice the scan angle and publishes
its amplitude (`ipd_gof_harmonic_amplitude`) and its phase
(`ipd_gof_harmonic_phase`). The phase is the position angle of the sub-PSF
structure, modulo 180 degrees. It is the only direct handle in DR3 on the
orientation of whatever is sitting inside the PSF.

THE DISCRIMINANT
----------------
A background blend is stationary; the star is not. The separation vector is

    s(t) = s0 - mu * (t - t0)

so as the mission baseline lengthens the drift term takes over and the
separation swings onto the proper-motion axis. In the limit
|mu| * dt >> |s0| the structure's position angle IS the proper motion's
position angle, PA_PM = atan2(pmra, pmdec), modulo 180.

A comoving structure -- a bound companion, a disc, an enshrouding shell -- has
no such preference. Its orientation is set by an orbit or an inclination and
is uncorrelated with where the star happens to be heading through the Galaxy.

So define

    dphi = (ipd_gof_harmonic_phase - PA_PM) mod 180

and ask whether it is concentrated. Concentration means background blends.
Uniformity means the structure moves with the star.

This is a per-object discriminant, not a population prior, and as far as we
can find it has not been used in any technosignature search.

WHY THE STATISTIC IS ROTATION-INVARIANT, AND WHY THE MIRROR CONTROL IS NOT
--------------------------------------------------------------------------
`ipd_gof_harmonic_phase` is a scan-angle phase and this project has no
independent way to verify that Gaia's zero point for it is the same
north-through-east zero point as PA_PM. A constant convention offset would
move the concentration to some other dphi without destroying it. The test is
therefore the Rayleigh statistic on the doubled angle 2*dphi, which detects a
concentration at *any* dphi, and the location of the concentration is reported
separately as a diagnostic rather than assumed.

That choice has a price, and it is stated here rather than buried: because the
statistic is invariant under rotation, the perpendicular mirror control this
project uses elsewhere is **degenerate**. Measuring dphi against PA_PM + 90
subtracts 90 degrees from every dphi and leaves the Rayleigh R numerically
unchanged. The perpendicular tail carries exactly zero information here, so it
cannot be quoted as a false-positive rate. The false-positive rate has to come
from permutation, below.

THE CONFOUND THAT WOULD FAKE THIS, AND THE NULL THAT KILLS IT
-------------------------------------------------------------
Both angles are smooth functions of sky position:

* PA_PM is dominated by solar reflex motion, so proper motions converge on the
  antapex and PA_PM is a smooth field on the sky.
* Gaia's scanning law samples only certain scan angles at a given position,
  and that sampling pattern varies with ecliptic latitude. A phase fitted from
  an incompletely sampled sinusoid is pulled toward the sampled directions, so
  the fitted phase also carries a smooth positional field.

Two smooth fields correlate with each other for reasons that have nothing to
do with any individual star. A global mark permutation -- shuffling PA_PM over
all positions -- destroys that correlation along with the signal, so it would
report the geometry as if it were discovery. This is the survey-geometry trap
in its usual costume.

The honest null is a **local** mark permutation: shuffle PA_PM only among
stars sharing a HEALPix cell (nside 8, ~7.3 deg). Within a cell both smooth
fields are near constant, so the shuffle preserves them and destroys only the
per-star pairing, which is the thing under test. Both nulls are run and the
gap between them is reported, because that gap is a direct measurement of how
much of the raw alignment is survey geometry.

POSITIVE CONTROLS
-----------------
A null here would be worthless without evidence that the instrument is
capable of showing the effect at all, and this project has already been bitten
once by a null that turned out to be a broken test. Three internal gradients
are predicted by the mechanism and checked:

1. AMPLITUDE. The phase of a sinusoid that was not detected is noise. So the
   alignment must strengthen with `ipd_gof_harmonic_amplitude`. If it does not,
   the phase column is not carrying orientation information and the channel is
   broken, not null.
2. PROPER MOTION. The drift term only wins when |mu| * dt exceeds |s0|. So the
   alignment must strengthen with total proper motion.
3. GAIA'S OWN BLEND FLAG. Sources with `ipd_frac_multi_peak` > 0 are ones where
   Gaia itself saw a second peak in the window. They must align more strongly
   than sources with a clean single peak.

If all three gradients are present, the discriminant works and a null on the
comoving population means something. If they are absent, the channel reports
itself broken.

THE APPLICATION
---------------
Channel 20 found a dim tail in the high-RUWE population and attributed it to
symmetric astrometric noise on the strength of how the excess moved with the
cut. That argument was indirect. This channel tests the same population
directly: if the dim tail is background blends, it must align with PA_PM more
than the bright tail does. The bright tail is the natural comparison here and
it is a genuine one, because the alignment statistic has no reason to prefer
either sign of the photometric residual.
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
from pipeline import sample as smp
from pipeline import statistics as st

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("searchV")

DIST_MIN, DIST_MAX = 69.0, 500.0
NSIDE = 8            # ~7.3 deg cells: small enough that PA_PM and the scanning
                     # law are near constant, large enough to hold ~200 stars.
N_PERM = 500
RNG_SEED = 20260820


# ---------------------------------------------------------------------------
# axial statistics
# ---------------------------------------------------------------------------

def rayleigh(dphi_deg: np.ndarray) -> dict:
    """Mean resultant length of axial data in [0, 180), on the doubled angle.

    The accompanying p is the textbook Rayleigh p = exp(-n R^2), and it is
    DESCRIPTIVE ONLY. It assumes both marginals are uniform, and in this data
    neither is: the raw harmonic phase piles up on the 0/180 axis (the scanning
    law does not sample scan angles evenly), and PA_PM is shaped by solar
    reflex motion. Two independent axial variables with non-uniform marginals
    have an expected resultant of roughly the product of their marginal
    resultants, not zero, so the Rayleigh p would report that product as
    significance. Every p this channel actually reasons from is a permutation
    p, which preserves both marginals by construction. The marginals are
    measured and stored alongside so the size of the effect is on the record.
    """
    th = np.deg2rad(2.0 * np.asarray(dphi_deg, float))
    n = th.size
    if n < 8:
        return {"n": int(n), "R": None, "Z": None, "p": None,
                "preferred_dphi_deg": None}
    c, s = float(np.mean(np.cos(th))), float(np.mean(np.sin(th)))
    r = float(np.hypot(c, s))
    z = n * r * r
    return {"n": int(n), "R": r, "Z": float(z),
            "p": float(np.exp(-z)),
            "preferred_dphi_deg": float(np.rad2deg(0.5 * np.arctan2(s, c)) % 180.0)}


def resultant(phase_deg: np.ndarray, pa_deg: np.ndarray) -> float:
    """Mean resultant length R of (phase - pa) as axial data. The test statistic."""
    th = np.deg2rad(2.0 * ((phase_deg - pa_deg) % 180.0))
    return float(np.hypot(np.mean(np.cos(th)), np.mean(np.sin(th))))


def permute_global(phase, pa, rng, n_perm=N_PERM):
    """Shuffle the marks over all positions. Destroys sky geometry too."""
    return np.array([resultant(phase, rng.permutation(pa)) for _ in range(n_perm)])


def permute_local(phase, pa, cell, rng, n_perm=N_PERM):
    """Shuffle the marks only within HEALPix cells.

    Preserves any smooth positional field in either angle, so what is left to
    destroy is the per-star pairing. This is the null the verdict rests on.
    """
    order = np.argsort(cell, kind="stable")
    cs = cell[order]
    bounds = np.flatnonzero(np.diff(cs)) + 1
    groups = np.split(order, bounds)
    groups = [g for g in groups if g.size > 1]
    out = np.empty(n_perm)
    for i in range(n_perm):
        pa_s = pa.copy()
        for g in groups:
            pa_s[g] = pa[rng.permutation(g)]
        out[i] = resultant(phase, pa_s)
    return out


def perm_p(obs: float, null: np.ndarray) -> float:
    """One-sided permutation p with the standard +1 so it can never be zero."""
    return float((np.count_nonzero(null >= obs) + 1) / (null.size + 1))


def within_cell_dispersion(pa_deg, cell) -> float:
    """Mean circular sd of PA_PM inside a cell, in degrees, on the doubled angle.

    This is what gives the local null its power. If PA_PM were nearly constant
    within a cell, shuffling inside the cell would barely change anything and
    the null would swallow a real per-star alignment along with the geometry.
    """
    th = np.deg2rad(2.0 * np.asarray(pa_deg, float))
    df = pd.DataFrame({"c": cell, "x": np.cos(th), "y": np.sin(th)})
    g = df.groupby("c")[["x", "y"]].agg(["mean", "count"])
    r = np.hypot(g[("x", "mean")], g[("y", "mean")]).to_numpy()
    n = g[("x", "count")].to_numpy()
    r = np.clip(r[n >= 20], 1e-6, 1 - 1e-9)
    # circular sd = sqrt(-2 ln R), halved to undo the angle doubling
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(r))).mean() / 2.0)


def sensitivity(phase, pa, cell, rng, n_perm,
                fractions=(0.0025, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20)):
    """Smallest injected aligned fraction this sample can still show.

    Sources are chosen at random and their phase replaced by their own PA_PM
    plus 15 deg of scatter, which is a deliberately imperfect alignment. The
    detection threshold is the same local-shuffle null the verdict uses.
    """
    rows, n = [], phase.size
    for f in fractions:
        ph = phase.copy()
        m = rng.random(n) < f
        ph[m] = (pa[m] + rng.normal(0, 15, int(m.sum()))) % 180.0
        r_obs = resultant(ph, pa)
        nb = permute_local(ph, pa, cell, rng, n_perm)
        p = perm_p(r_obs, nb)
        rows.append({"injected_fraction": f, "R": r_obs,
                     "null_R": float(nb.mean()), "null_sd": float(nb.std()),
                     "p_local": p, "detected": bool(p < 0.01)})
    det = [r["injected_fraction"] for r in rows if r["detected"]]
    return {"scan": rows,
            "min_detectable_fraction": (min(det) if det else None),
            "permutation_p_floor": float(1.0 / (n_perm + 1)),
            "note": ("Aligned fraction recoverable at p<0.01 against the local "
                     "null, injected at 15 deg scatter on the real positions. "
                     "A null from this channel excludes aligned populations at "
                     "or above this fraction and says nothing below it. The "
                     "permutation p cannot fall below the floor quoted here, "
                     "so a scan that detects every injected fraction is "
                     "resolution-limited and its lowest entry is an upper "
                     "bound on the true threshold, not the threshold.")}


def summarise(name, phase, pa, out_list):
    r = rayleigh((phase - pa) % 180.0)
    r["subset"] = name
    out_list.append(r)
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    # Exists so the analysis can be exercised end-to-end against a partial
    # merge while the full pull is still running. A dry run on a subset finds
    # the crash in the extinction path; it does not license a verdict, so the
    # tag must be changed too and the JSON lands under that tag.
    ap.add_argument("--source", default=None,
                    help="override the input parquet (dry runs only)")
    # The reference sample that defines the fiducial is a property of the
    # main pipeline run, not of this channel's output label, so a dry run can
    # relabel its output without losing the reference.
    ap.add_argument("--ref-tag", default=None,
                    help="tag of the reference residual parquet (default: --tag)")
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    src = Path(args.source) if args.source else cfg.RAW_DIR / "high_ruwe_500pc.parquet"
    d = pd.read_parquet(src)
    need = ["ipd_gof_harmonic_amplitude", "ipd_gof_harmonic_phase"]
    missing = [c for c in need if c not in d.columns]
    if missing:
        log.error("%s lacks %s -- run scripts/84_pull_ipd_harmonic.py first",
                  src.name, missing)
        return 1

    n_raw = len(d)
    # The phase is only defined where the harmonic was fitted at all.
    ok = (d["ipd_gof_harmonic_phase"].notna()
          & d["ipd_gof_harmonic_amplitude"].notna()
          & d["pmra"].notna() & d["pmdec"].notna())
    d = d[ok].reset_index(drop=True)
    log.info("IPD phase available for %d/%d (%.1f%%)",
             len(d), n_raw, 100 * len(d) / n_raw)
    if len(d) < 1000:
        log.error("only %d usable sources; this test cannot run", len(d))
        return 1

    pmra = d["pmra"].to_numpy(float)
    pmdec = d["pmdec"].to_numpy(float)
    pa_pm = np.rad2deg(np.arctan2(pmra, pmdec)) % 180.0   # north through east
    pm_tot = np.hypot(pmra, pmdec)
    phase = d["ipd_gof_harmonic_phase"].to_numpy(float) % 180.0
    amp = d["ipd_gof_harmonic_amplitude"].to_numpy(float)
    dphi = (phase - pa_pm) % 180.0

    import healpy as hp
    cell = hp.ang2pix(NSIDE, np.radians(90.0 - d["dec"].to_numpy(float)),
                      np.radians(d["ra"].to_numpy(float)))

    out = {"tag": args.tag, "n_input": int(n_raw), "n_usable": int(len(d)),
           "nside": NSIDE, "n_perm": int(args.n_perm), "seed": RNG_SEED,
           "mirror_control": (
               "DEGENERATE BY CONSTRUCTION and therefore not quoted. The "
               "statistic is the Rayleigh R of the doubled angle, which is "
               "invariant under a constant rotation of dphi, so measuring "
               "against PA_PM+90 returns the identical R. The false-positive "
               "rate here comes from the local mark permutation instead.")}

    # ---- how non-uniform each marginal is, on its own --------------------
    # This is why the analytic Rayleigh p cannot be the test. Two independent
    # axial variables with these marginals already produce a resultant near
    # the product of the two, and the permutation nulls should reproduce that
    # number. If they do not, the permutation is not preserving what it claims.
    marg = {"phase": rayleigh(phase), "pa_pm": rayleigh(pa_pm)}
    out["marginals"] = {
        "phase_R": marg["phase"]["R"], "pa_pm_R": marg["pa_pm"]["R"],
        "product_R": float(marg["phase"]["R"] * marg["pa_pm"]["R"]),
        "note": ("Expected resultant for two INDEPENDENT axial variables with "
                 "these marginals, against which the permutation null means "
                 "are the empirical check.")}
    log.info("marginal R: phase %.5f, PA_PM %.5f -> independent product %.5f",
             marg["phase"]["R"], marg["pa_pm"]["R"],
             out["marginals"]["product_R"])

    # ---- global result and the two nulls ---------------------------------
    obs = rayleigh(dphi)
    obs["p_is_descriptive_only"] = True
    out["global"] = obs
    log.info("global: n=%d R=%.5f preferred dphi=%.1f deg Rayleigh p=%.3g",
             obs["n"], obs["R"], obs["preferred_dphi_deg"], obs["p"])

    log.info("global mark permutation (%d draws) ...", args.n_perm)
    ng = permute_global(phase, pa_pm, rng, args.n_perm)
    log.info("local mark permutation within nside=%d cells (%d draws) ...",
             NSIDE, args.n_perm)
    nl = permute_local(phase, pa_pm, cell, rng, args.n_perm)

    out["null_global_shuffle"] = {
        "mean_R": float(ng.mean()), "sd_R": float(ng.std()),
        "p_obs": perm_p(obs["R"], ng)}
    out["null_local_shuffle"] = {
        "n_cells_used": int(np.unique(cell).size),
        "mean_R": float(nl.mean()), "sd_R": float(nl.std()),
        "p_obs": perm_p(obs["R"], nl)}
    out["geometry_leakage"] = {
        "excess_R_of_local_null_over_global_null":
            float(nl.mean() - ng.mean()),
        "note": ("How much apparent alignment survives destroying the "
                 "per-star pairing while keeping sky position. Any of the "
                 "observed R below this level is survey geometry, not "
                 "structure.")}
    log.info("null R: global %.5f+-%.5f, local %.5f+-%.5f; observed %.5f",
             ng.mean(), ng.std(), nl.mean(), nl.std(), obs["R"])
    log.info("p(observed | global null) = %.4g ; p(observed | local null) = %.4g",
             out["null_global_shuffle"]["p_obs"],
             out["null_local_shuffle"]["p_obs"])

    # ---- positive controls: three predicted gradients ---------------------
    n_bin_perm = max(200, args.n_perm // 2)

    def gradient(values, label, edges=None):
        """Alignment excess in quartiles of `values`, each with its own null.

        Strict monotonicity across four bins is not used as the criterion. It
        is a knife-edge on noise: one bin pair inverting by less than its own
        error kills a gradient that is plainly present. The criterion is the
        top-minus-bottom difference measured against the two bins' own
        permutation spreads, with the Spearman rank correlation across bins
        reported alongside as the shape check.
        """
        q = np.nanpercentile(values, [25, 50, 75]) if edges is None else edges
        bins, rows = np.digitize(values, q), []
        for b in range(len(q) + 1):
            m = bins == b
            if m.sum() < 200:
                continue
            r_obs = resultant(phase[m], pa_pm[m])
            nb = permute_local(phase[m], pa_pm[m], cell[m], rng, n_bin_perm)
            rows.append({"bin": int(b), "n": int(m.sum()),
                         "lo": float(values[m].min()),
                         "hi": float(values[m].max()),
                         "R": r_obs, "null_R": float(nb.mean()),
                         "null_sd": float(nb.std()),
                         "excess_R": float(r_obs - nb.mean()),
                         "p_local": perm_p(r_obs, nb)})
        log.info("gradient %s: excess_R by bin = %s", label,
                 ", ".join(f"{r['excess_R']:+.5f}" for r in rows))
        return {"bins": rows, "trend": trend(rows)}

    def trend(rows):
        """Top-minus-bottom excess in units of its own permutation error."""
        if len(rows) < 2:
            return {"status": "too few bins"}
        lo, hi = rows[0], rows[-1]
        diff = hi["excess_R"] - lo["excess_R"]
        sd = float(np.hypot(hi["null_sd"], lo["null_sd"]))
        e = [r["excess_R"] for r in rows]
        rho = float(np.corrcoef(np.argsort(np.argsort(e)),
                                np.arange(len(e)))[0, 1]) if len(e) > 2 else None
        return {"top_minus_bottom": float(diff), "sigma": sd,
                "nsigma": float(diff / sd) if sd > 0 else None,
                "spearman_rho_over_bins": rho,
                # Two sigma and the right sign. Below that the gradient is not
                # established, which is a statement about this test's power on
                # this sample and not about the physics.
                "passes": bool(sd > 0 and diff / sd > 2.0)}

    out["control_amplitude"] = gradient(amp, "ipd_gof_harmonic_amplitude")
    out["control_proper_motion"] = gradient(pm_tot, "total proper motion")

    mp = d["ipd_frac_multi_peak"].fillna(0).to_numpy(float)
    out["control_multi_peak"] = gradient(mp, "ipd_frac_multi_peak",
                                         edges=np.array([0.5, 2.5, 10.0]))

    grads = {k: bool(out[f"control_{k}"]["trend"].get("passes"))
             for k in ("amplitude", "proper_motion", "multi_peak")}
    out["controls_pass"] = grads
    n_pass = sum(grads.values())
    for k in ("amplitude", "proper_motion", "multi_peak"):
        t = out[f"control_{k}"]["trend"]
        log.info("control %-13s top-bottom %+.5f = %+.2f sigma  rho=%s  %s",
                 k, t.get("top_minus_bottom", float("nan")),
                 t.get("nsigma") or float("nan"), t.get("spearman_rho_over_bins"),
                 "PASS" if grads[k] else "fail")

    # ---- sensitivity: inject alignment into THIS sky and re-measure -------
    # The generic injection test lives in tests/. This one runs on the real
    # positions, the real PA_PM field and the real cell occupancies, because
    # the local null's power depends on how much PA_PM actually varies inside
    # a cell. If it varies little the null would absorb a genuine signal along
    # with the geometry, and the channel's null would mean nothing. This
    # measures the smallest injected aligned fraction the channel can still
    # see, and that number is the honest statement of what a null excludes.
    out["within_cell_pa_dispersion_deg"] = float(within_cell_dispersion(pa_pm, cell))
    out["sensitivity"] = sensitivity(phase, pa_pm, cell, rng, n_bin_perm)
    log.info("within-cell PA_PM dispersion %.1f deg; detectable aligned "
             "fraction >= %s", out["within_cell_pa_dispersion_deg"],
             out["sensitivity"].get("min_detectable_fraction"))

    # ---- application: channel 20's dim tail against its bright tail -------
    app = apply_to_dim_tail(d, phase, pa_pm, cell, rng, args)
    out["dim_vs_bright_tail"] = app

    # ---- verdict ----------------------------------------------------------
    p_local = out["null_local_shuffle"]["p_obs"]
    excess_R = obs["R"] - nl.mean()
    if n_pass == 0:
        verdict = (
            f"BROKEN, NOT NULL. None of the three predicted gradients is "
            f"present: alignment does not strengthen with harmonic amplitude, "
            f"with proper motion, or with Gaia's own multi-peak flag. The "
            f"phase column is therefore not carrying usable orientation "
            f"information for this population, and no statement about "
            f"comoving structure can be made from it. Observed R = "
            f"{obs['R']:.5f} against a local-shuffle null of {nl.mean():.5f}.")
    elif p_local > 0.01:
        verdict = (
            f"NULL, with the discriminant demonstrated to work. "
            f"{n_pass}/3 positive controls show the predicted gradient at "
            f">2 sigma, so the phase column does carry orientation. Yet the "
            f"population as a whole shows no alignment with the proper-motion "
            f"axis beyond what survives destroying the per-star pairing: "
            f"R = {obs['R']:.5f} against a local-shuffle null of "
            f"{nl.mean():.5f} +- {nl.std():.5f}, p = {p_local:.3g}. The null "
            f"excludes an aligned fraction at or above "
            f"{out['sensitivity'].get('min_detectable_fraction')} and says "
            f"nothing below it. Note that the global shuffle would have "
            f"reported p = {out['null_global_shuffle']['p_obs']:.3g}; the "
            f"difference is survey geometry, and taking the global number "
            f"would have been the mistake this channel was built to avoid.")
    else:
        verdict = (
            f"ALIGNMENT DETECTED. R = {obs['R']:.5f} against a local-shuffle "
            f"null of {nl.mean():.5f} +- {nl.std():.5f} (p = {p_local:.3g}), "
            f"an excess of {excess_R:+.5f} concentrated near dphi = "
            f"{obs['preferred_dphi_deg']:.1f} deg. {n_pass}/3 positive "
            f"controls pass. This is a measurement of the background-"
            f"blend contamination of the high-RUWE population, not a "
            f"technosignature: alignment with the proper-motion axis is the "
            f"signature of a stationary contaminant the star is drifting away "
            f"from. The comoving -- that is, non-aligned -- remainder is the "
            f"population any future channel should search.")
    out["verdict"] = verdict

    print("\n" + "=" * 78)
    print("SEARCH V: IPD harmonic phase against the proper-motion axis")
    print("=" * 78)
    print(f"  usable sources                    {obs['n']:,}")
    print(f"  observed R                        {obs['R']:.5f}")
    print(f"  local-shuffle null R              {nl.mean():.5f} +- {nl.std():.5f}")
    print(f"  global-shuffle null R             {ng.mean():.5f} +- {ng.std():.5f}")
    print(f"  p (local null, the honest one)    {p_local:.4g}")
    print(f"  p (global null, would overstate)  {out['null_global_shuffle']['p_obs']:.4g}")
    print(f"  positive controls passed          {n_pass}/3  {grads}")
    if app.get("dim") and app.get("bright"):
        print(f"  dim tail   excess R               {app['dim']['excess_R']:+.5f}"
              f"  (n={app['dim']['n']:,})")
        print(f"  bright tail excess R              {app['bright']['excess_R']:+.5f}"
              f"  (n={app['bright']['n']:,})")
    print(f"\nVERDICT: {verdict}\n")

    p = cfg.RESULT_DIR / f"searchV_harmonic_phase_{args.tag}.json"
    p.write_text(json.dumps(out, indent=2))
    log.info("wrote %s", p)
    return 0


def apply_to_dim_tail(d, phase, pa_pm, cell, rng, args) -> dict:
    """Channel 20's dim tail against its bright tail, in alignment.

    Rebuilds channel 20's residual with channel 20's estimator -- the same
    polynomial refitted on the clean low-RUWE reference -- so that "dim" here
    means exactly what it meant there. Using a different estimator for the two
    is the error that manufactured channel 20's first null, and it is not
    going to be repeated in the channel that audits it.
    """
    try:
        ref_tag = args.ref_tag or args.tag
        ref = pd.read_parquet(cfg.DERIVED_DIR / f"{ref_tag}_resid.parquet",
                              columns=["M_G", "M_Ks", "ruwe", "cstar_nsigma"])
    except Exception as exc:                       # noqa: BLE001
        log.warning("reference residuals unavailable (%s); skipping the "
                    "dim/bright application", str(exc)[:80])
        return {"status": f"unavailable: {str(exc)[:80]}"}

    clean = (ref["M_G"].notna() & ref["M_Ks"].notna()
             & (ref["ruwe"] < 1.4) & (ref["cstar_nsigma"].abs() < 3))
    x = ref.loc[clean, "M_Ks"].to_numpy(float)
    y = ref.loc[clean, "M_G"].to_numpy(float)
    best = None
    for deg in (3, 4, 5):
        c = np.polyfit(x, y, deg)
        s = st.robust_sigma(y - np.polyval(c, x))
        if best is None or s < best[1]:
            best = (c, s, deg)
    coef, _, deg = best

    e = d.copy()
    keep = e["tmass_ks_m"].notna().to_numpy()
    e = smp.add_astrometry(e[keep].reset_index(drop=True))
    ph, pa, ce = phase[keep], pa_pm[keep], cell[keep]

    a0 = ext.query_a0("edenhofer23", e["l"].to_numpy(float),
                      e["b"].to_numpy(float), e["dist_pc"].to_numpy(float))
    bp_rp = e["bp_rp"].to_numpy(float)
    a_g = ext.deredden("G", np.nan_to_num(a0), bp_rp)
    a_ks = ext.deredden("Ks", np.nan_to_num(a0), bp_rp)
    mu = e["dist_mod"].to_numpy(float)
    m_g = e["phot_g_mean_mag"].to_numpy(float) - mu - a_g
    m_ks = e["tmass_ks_m"].to_numpy(float) - mu - a_ks
    resid = m_g - np.polyval(coef, m_ks)

    box = (np.isfinite(resid) & np.isfinite(a0)
           & (m_ks > 3.0) & (m_ks < 8.0)
           & (bp_rp > 0.7) & (bp_rp < 3.6)
           & (e["dist_pc"].to_numpy(float) > DIST_MIN)
           & (e["dist_pc"].to_numpy(float) < DIST_MAX)
           & (a_g < 0.5))
    resid, ph, pa, ce = resid[box], ph[box], pa[box], ce[box]
    if resid.size < 500:
        return {"status": f"only {int(resid.size)} sources in the box"}

    sig, med = float(st.robust_sigma(resid)), float(np.median(resid))
    res = {"status": "ok", "fit_degree": int(deg), "n_box": int(resid.size),
           "sigma": sig, "median": med, "k": 2.0}
    for name, m in (("dim", resid > med + 2.0 * sig),
                    ("bright", resid < med - 2.0 * sig)):
        if m.sum() < 200:
            res[name] = {"n": int(m.sum()), "status": "too few"}
            continue
        r_obs = resultant(ph[m], pa[m])
        nb = permute_local(ph[m], pa[m], ce[m], rng, max(200, args.n_perm // 2))
        res[name] = {"n": int(m.sum()), "R": r_obs,
                     "null_R": float(nb.mean()), "null_sd": float(nb.std()),
                     "excess_R": float(r_obs - nb.mean()),
                     "p_local": perm_p(r_obs, nb)}
    if isinstance(res.get("dim"), dict) and "excess_R" in res["dim"] \
            and isinstance(res.get("bright"), dict) and "excess_R" in res["bright"]:
        res["dim_minus_bright_excess_R"] = (res["dim"]["excess_R"]
                                            - res["bright"]["excess_R"])
        res["interpretation"] = (
            "A positive dim-minus-bright difference means the dim tail is "
            "more contaminated by stationary background structure than the "
            "bright tail, which is the direct version of the argument "
            "channel 20 had to make indirectly from how its excess moved "
            "with the cut. A difference consistent with zero means the two "
            "tails are equally blended and the dim tail's origin is "
            "elsewhere.")
    return res


if __name__ == "__main__":
    raise SystemExit(main())
