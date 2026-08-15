#!/usr/bin/env python
"""Differential test inside co-natal wide binaries.

    run.sh scripts/24_wide_binaries.py --tag primary

WHY THIS SHOULD BEAT THE SINGLE-STAR TEST
------------------------------------------
The systematic budget is dominated by a 0.068 mag term that is *not*
instrumental: main-sequence structure at fixed (M_Ks, [M/H], J-Ks) driven by
age, alpha-enhancement, rotation and activity.  Two components of a wide binary
were born together, so they share age, metallicity, distance and (to first
order) the dust column.  Take the difference of their residuals and every
common-mode term cancels:

    dr = r_A - r_B

What survives is the genuinely per-star part: photometric noise, blending, and
whatever activity/rotation difference the two stars have.  If the common-mode
fraction is large, sigma(dr) will come out well below sqrt(2) * sigma_single,
and the effective per-star scatter sigma_individual = sigma(dr)/sqrt(2) is the
number that sets the new per-star reach.

If sigma(dr) ~ sqrt(2) * sigma_single instead, there is no common-mode
cancellation, the 0.068 mag term is per-star after all, and this whole avenue
is dead.  That is a real possible outcome and the script reports it plainly.

PAIRS ARE FOUND IN OUR OWN SAMPLE, NOT DOWNLOADED
--------------------------------------------------
Wide binaries are *resolved*, so both components appear as separate Gaia
sources and both already survived every quality cut, including RUWE.  A
self-join on the analysis sample is therefore the right move: no new catalogue,
no new cross-match, and both components already have residuals from the same
fitted fiducial.

Selection follows the logic of El-Badry, Rix & Heintz (2021, MNRAS 506, 2269):
small projected separation, consistent parallax, and a proper-motion difference
small enough to be a bound orbit.  Chance alignments are estimated by a
scramble test rather than assumed negligible.

KNOWN APPROXIMATION: proper-motion uncertainties were not pulled in the
original ADQL.  Gaia DR3 astrometric errors scale together, so sigma_mu is
approximated as 1.2 * parallax_error; the factor is crude but the criterion is
dominated by the orbital-velocity term at these separations, and the scramble
test measures whatever contamination it lets through.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from pipeline import config as cfg
from pipeline import statistics as st

# --- selection constants, all justified in the docstring ------------------
MAX_ANGSEP_ARCSEC = 120.0     # generous; physical separation cut does the work
MAX_PHYS_SEP_AU = 20_000.0    # El-Badry+2021 contamination stays low below this
MIN_PHYS_SEP_AU = 100.0       # below this Gaia's own completeness drops
PARALLAX_CONSISTENCY_NSIGMA = 3.0
MU_ERR_PER_PARALLAX_ERR = 1.2  # crude DR3 scaling, see docstring
MU_MARGIN_NSIGMA = 2.0
# Total mass assumed when computing the maximum bound orbital velocity.
# Deliberately generous (both components could be ~1 Msun) so the criterion
# errs toward keeping pairs; the scramble test then measures contamination.
M_TOT_MSUN = 2.0


def orbital_mu_max(sep_au: np.ndarray, dist_pc: np.ndarray) -> np.ndarray:
    """Max proper-motion difference for a bound orbit, mas/yr.

    v_esc [km/s] = 29.78 * sqrt(2 * M_tot / a_AU)   (29.78 km/s = Earth's orbit)
    mu [mas/yr]  = v [km/s] / (4.74 * d_kpc)
    """
    v_esc = 29.78 * np.sqrt(2.0 * M_TOT_MSUN / np.maximum(sep_au, 1.0))
    return v_esc / (4.74 * np.maximum(dist_pc, 1.0) / 1000.0)


def find_pairs(d: pd.DataFrame, scramble: bool = False,
               seed: int = 7) -> pd.DataFrame:
    """Self-join the sample for co-moving, co-distant, close pairs."""
    ra = np.radians(d["ra"].to_numpy(float))
    dec = np.radians(d["dec"].to_numpy(float))
    xyz = np.column_stack([np.cos(dec) * np.cos(ra),
                           np.cos(dec) * np.sin(ra),
                           np.sin(dec)])

    plx = d["parallax_corr"].to_numpy(float)
    plx_err = d["parallax_error"].to_numpy(float)
    pmra = d["pmra"].to_numpy(float)
    pmdec = d["pmdec"].to_numpy(float)
    dist = d["dist_pc"].to_numpy(float)

    if scramble:
        # Destroy real pairs while preserving every marginal distribution:
        # keep positions, permute the astrometry.
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(d))
        plx, plx_err = plx[perm], plx_err[perm]
        pmra, pmdec, dist = pmra[perm], pmdec[perm], dist[perm]

    r_chord = 2.0 * np.sin(np.radians(MAX_ANGSEP_ARCSEC / 3600.0) / 2.0)
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(r_chord, output_type="ndarray")
    if len(pairs) == 0:
        return pd.DataFrame()

    i, j = pairs[:, 0], pairs[:, 1]
    dot = np.clip(np.einsum("ij,ij->i", xyz[i], xyz[j]), -1.0, 1.0)
    theta_arcsec = np.degrees(np.arccos(dot)) * 3600.0

    dist_mean = 0.5 * (dist[i] + dist[j])
    sep_au = theta_arcsec * dist_mean

    dplx = np.abs(plx[i] - plx[j])
    sig_dplx = np.hypot(plx_err[i], plx_err[j])
    dmu = np.hypot(pmra[i] - pmra[j], pmdec[i] - pmdec[j])
    sig_dmu = MU_ERR_PER_PARALLAX_ERR * sig_dplx
    mu_max = orbital_mu_max(sep_au, dist_mean)

    keep = ((sep_au > MIN_PHYS_SEP_AU) & (sep_au < MAX_PHYS_SEP_AU)
            & (dplx < PARALLAX_CONSISTENCY_NSIGMA * sig_dplx)
            & (dmu < mu_max + MU_MARGIN_NSIGMA * sig_dmu))

    return pd.DataFrame({
        "i": i[keep], "j": j[keep],
        "theta_arcsec": theta_arcsec[keep], "sep_au": sep_au[keep],
        "dist_pc": dist_mean[keep], "dplx_nsigma": (dplx / sig_dplx)[keep],
        "dmu": dmu[keep], "mu_max": mu_max[keep],
    })


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--resid-col", default="residual")
    args = ap.parse_args()

    cols = ["source_id", "ra", "dec", "parallax_corr", "parallax_error",
            "pmra", "pmdec", "dist_pc", "M_Ks", "M_G", "bp_rp0", "A_0",
            "ruwe", "cstar_nsigma", args.resid_col]
    d = pd.read_parquet(cfg.DERIVED_DIR / f"{args.tag}_resid.parquet",
                        columns=cols).reset_index(drop=True)
    r = d[args.resid_col].to_numpy(float)
    sigma_single = st.robust_sigma(r)
    print(f"sample: {len(d):,} stars, single-star sigma = {sigma_single:.5f} mag\n")

    pairs = find_pairs(d)
    print(f"pairs found                 : {len(pairs):,}")
    if len(pairs) < 200:
        print("too few pairs for a useful measurement")
        return 1

    # --- chance alignment ------------------------------------------------
    fake = find_pairs(d, scramble=True)
    contam = len(fake) / max(len(pairs), 1)
    print(f"scrambled-sky false pairs   : {len(fake):,}")
    print(f"estimated chance-alignment contamination : {100*contam:.2f}%\n")

    i, j = pairs["i"].to_numpy(), pairs["j"].to_numpy()
    # A star can appear in several candidate pairs (triples, or chance
    # overlaps). Keep only stars used once, so every pair is independent.
    used = np.concatenate([i, j])
    counts = pd.Series(used).value_counts()
    multi = set(counts[counts > 1].index)
    solo = ~(pd.Series(i).isin(multi).to_numpy()
             | pd.Series(j).isin(multi).to_numpy())
    print(f"pairs after removing shared components : {int(solo.sum()):,} "
          f"(dropped {int((~solo).sum()):,})")
    pairs = pairs[solo].reset_index(drop=True)
    i, j = pairs["i"].to_numpy(), pairs["j"].to_numpy()

    # --- the measurement -------------------------------------------------
    r_i, r_j = r[i], r[j]
    # Order within each pair by M_Ks so "primary" is the more luminous star;
    # this makes dr a defined quantity rather than a sign-ambiguous one.
    mks_i, mks_j = d["M_Ks"].to_numpy(float)[i], d["M_Ks"].to_numpy(float)[j]
    prim_is_i = mks_i <= mks_j
    r_prim = np.where(prim_is_i, r_i, r_j)
    r_sec = np.where(prim_is_i, r_j, r_i)
    dr = r_prim - r_sec

    sigma_dr = st.robust_sigma(dr)
    sigma_indiv = sigma_dr / np.sqrt(2.0)
    no_cancel = np.sqrt(2.0) * sigma_single
    common_frac = 1.0 - (sigma_indiv / sigma_single) ** 2

    print(f"\n=== differential scatter ===")
    print(f"  sigma(single star)                 : {sigma_single:.5f} mag")
    print(f"  sqrt(2)*sigma  [no cancellation]   : {no_cancel:.5f} mag")
    print(f"  sigma(dr) MEASURED                 : {sigma_dr:.5f} mag")
    print(f"  -> sigma per star after cancelling : {sigma_indiv:.5f} mag")
    print(f"  -> common-mode variance fraction   : {100*common_frac:.1f}%")

    thr_single = 5.0 * sigma_single
    thr_pair = 5.0 * sigma_indiv
    print(f"\n  5-sigma deficit threshold, single  : {thr_single:.4f} mag "
          f"-> f = {st.fraction_from_delta(thr_single):.3f}")
    print(f"  5-sigma deficit threshold, paired  : {thr_pair:.4f} mag "
          f"-> f = {st.fraction_from_delta(thr_pair):.3f}")

    # --- tail of the differential ---------------------------------------
    med = float(np.median(dr))
    tail = []
    for k in (3, 4, 5, 6):
        npos = int(np.count_nonzero(dr > med + k * sigma_dr))
        nneg = int(np.count_nonzero(dr < med - k * sigma_dr))
        tail.append({"k": k, "n_pos": npos, "n_neg": nneg,
                     "frac_pos": npos / len(dr), "frac_neg": nneg / len(dr)})
    tdf = pd.DataFrame(tail)
    print("\n=== differential tail (sign convention: primary fainter = +) ===")
    print(tdf.to_string(index=False, float_format=lambda v: f"{v:10.6f}"))
    print("  NOTE: unlike the single-star test this statistic is intrinsically")
    print("  two-sided -- either component can be the harvested one -- so the")
    print("  two tails are each other's control.")

    # --- why is the tail asymmetric? -------------------------------------
    # dr < 0 means the SECONDARY looks too faint in G at its M_Ks.  If that is
    # the 2MASS aperture-blending mechanism of RESULTS.md 7b, it must be
    # driven by the brighter primary leaking into the secondary's 4-arcsec Ks
    # beam -- and it must therefore get worse at small ANGULAR separation and
    # at large magnitude contrast.  If instead it is intrinsic M-dwarf
    # activity, it should be flat in angular separation.
    print("\n=== is the asymmetry angular-separation dependent? ===")
    theta = pairs["theta_arcsec"].to_numpy()
    dmag = np.abs(d["M_G"].to_numpy(float)[i] - d["M_G"].to_numpy(float)[j])
    rows = []
    edges = np.nanquantile(theta, [0, 0.25, 0.5, 0.75, 1.0])
    for a in range(len(edges) - 1):
        m = (theta >= edges[a]) & (theta <= edges[a + 1] if a == len(edges) - 2
                                   else theta < edges[a + 1])
        if m.sum() < 50:
            continue
        s_loc = st.robust_sigma(dr[m])
        npos = int(np.count_nonzero(dr[m] > med + 5 * sigma_dr))
        nneg = int(np.count_nonzero(dr[m] < med - 5 * sigma_dr))
        rows.append({"theta_range": f"{edges[a]:.1f}-{edges[a+1]:.1f}\"",
                     "n_pairs": int(m.sum()), "sigma_dr": s_loc,
                     "n_pos_5sig": npos, "n_neg_5sig": nneg,
                     "neg_over_pos": nneg / max(npos, 1),
                     "median_dmag_G": float(np.median(dmag[m]))})
    sep_tab = pd.DataFrame(rows)
    print(sep_tab.to_string(index=False, float_format=lambda v: f"{v:9.4f}"))
    sep_tab.to_csv(cfg.RESULT_DIR / f"wide_binaries_sepsplit_{args.tag}.csv",
                   index=False)

    out = pairs.copy()
    out["source_id_prim"] = np.where(prim_is_i, d["source_id"].to_numpy()[i],
                                     d["source_id"].to_numpy()[j])
    out["source_id_sec"] = np.where(prim_is_i, d["source_id"].to_numpy()[j],
                                    d["source_id"].to_numpy()[i])
    out["r_prim"], out["r_sec"], out["dr"] = r_prim, r_sec, dr
    out["dr_nsigma"] = (dr - med) / sigma_dr
    out.to_parquet(cfg.DERIVED_DIR / f"wide_binaries_{args.tag}.parquet",
                   index=False)

    summary = {
        "tag": args.tag, "n_stars": int(len(d)),
        "n_pairs": int(len(pairs)),
        "chance_contamination_frac": float(contam),
        "median_sep_au": float(np.median(pairs["sep_au"])),
        "sigma_single_mag": float(sigma_single),
        "sigma_dr_mag": float(sigma_dr),
        "sqrt2_sigma_single_mag": float(no_cancel),
        "sigma_individual_mag": float(sigma_indiv),
        "common_mode_variance_fraction": float(common_frac),
        "f_threshold_single": float(st.fraction_from_delta(thr_single)),
        "f_threshold_paired": float(st.fraction_from_delta(thr_pair)),
        "verdict": ("common-mode cancellation WORKS; paired test reaches lower f"
                    if sigma_dr < 0.9 * no_cancel else
                    "NO useful cancellation; the scatter is per-star, this "
                    "avenue is dead"),
    }
    (cfg.RESULT_DIR / f"wide_binaries_{args.tag}.json").write_text(
        json.dumps(summary, indent=2))
    tdf.to_csv(cfg.RESULT_DIR / f"wide_binaries_tail_{args.tag}.csv", index=False)
    print(f"\nVERDICT: {summary['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
