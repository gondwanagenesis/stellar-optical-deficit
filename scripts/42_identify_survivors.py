#!/usr/bin/env python
"""Identify the 513 stars that survived every cut in Search A.

    run.sh scripts/42_identify_survivors.py --tag primary

These have: a significant dimming absorber, a spectral slope significantly
flatter than interstellar dust, an acceptable fit, clean photometry and
astrometry, low extinction, high Galactic latitude, and no membership in a known
star-forming complex. They are what Search A produces as candidates.

SIMBAD is queried by an OR-of-cones TAP join -- astroquery's vector
query_region gives no reliable row-to-input mapping in current versions, which
silently returns all-None classifications.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from pipeline import config as cfg

CONTAMINANT = {
    "YSO", "Y*O", "Y*?", "TT*", "TT?", "Or*", "pr*", "out", "HH",
    "Em*", "Be*", "WR*", "s*b", "sg*", "PN", "PN?", "RG*", "AGB", "C*", "S*",
    "HB*", "pA*", "WD*", "WD?", "LP?", "LP*", "Mi*", "Ce*", "RR*", "dS*",
    "V*", "Ir*", "cC*", "gam", "BY*", "RS*", "Ro*", "Fl*", "UV*", "Er*",
    "EB*", "SB*", "Al*", "bL*", "WU*", "**", "SB?", "EB?", "No*", "LM*",
}


def simbad(ra, dec, radius_arcsec=3.0, chunk=25):
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url="https://simbad.cds.unistra.fr/simbad/sim-tap")
    r_deg = radius_arcsec / 3600.0
    n = len(ra)
    main = np.array([None] * n, dtype=object)
    oty = np.array([None] * n, dtype=object)
    sep = np.full(n, np.nan)
    for i0 in range(0, n, chunk):
        i1 = min(i0 + chunk, n)
        cones = " OR ".join(
            f"CONTAINS(POINT('ICRS',b.ra,b.dec),"
            f"CIRCLE('ICRS',{ra[j]:.7f},{dec[j]:.7f},{r_deg:.7f}))=1"
            for j in range(i0, i1))
        q = (f"SELECT b.main_id,b.otype_txt,b.sp_type,b.ra,b.dec FROM basic AS b "
             f"WHERE b.ra IS NOT NULL AND ({cones})")
        res = None
        for a in range(3):
            try:
                res = tap.launch_job(q).get_results()
                break
            except Exception as e:
                time.sleep(4 * (a + 1))
        if res is None or len(res) == 0:
            continue
        t = res.to_pandas()
        rr, rd = t["ra"].to_numpy(float), t["dec"].to_numpy(float)
        for j in range(i0, i1):
            cd = np.cos(np.radians(dec[j]))
            ds = np.hypot((rr - ra[j]) * cd, rd - dec[j])
            k = int(np.argmin(ds))
            if ds[k] <= r_deg:
                sep[j] = ds[k] * 3600
                a_, b_ = t.iloc[k].get("main_id"), t.iloc[k].get("otype_txt")
                main[j] = a_.decode() if isinstance(a_, bytes) else a_
                oty[j] = b_.decode() if isinstance(b_, bytes) else b_
        print(f"  SIMBAD {i1}/{n}", flush=True)
    return main, oty, sep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="primary")
    ap.add_argument("--max-objects", type=int, default=520)
    args = ap.parse_args()

    s = pd.read_csv(cfg.RESULT_DIR / f"slope_survivors_{args.tag}.csv")
    s = s.nlargest(args.max_objects, "deficit_G_mag").reset_index(drop=True)
    print(f"identifying {len(s):,} Search-A survivors")

    m, o, sp = simbad(s["ra"].to_numpy(float), s["dec"].to_numpy(float))
    s["simbad_id"], s["otype"], s["sep_arcsec"] = m, o, sp
    s["known_contaminant"] = s["otype"].isin(CONTAMINANT)

    matched = s["otype"].notna()
    print(f"\nmatched in SIMBAD : {int(matched.sum()):,} / {len(s):,}")
    print("\n=== object types ===")
    vc = s.loc[matched, "otype"].value_counts()
    for k, v in vc.head(20).items():
        flag = "  <-- known contaminant" if k in CONTAMINANT else ""
        print(f"  {str(k):12s} {v:5d}{flag}")

    clean = s[~s["known_contaminant"]]
    print(f"\n=== after removing known contaminant classes ===")
    print(f"  {len(clean):,} remain ({int(matched.sum() - s['known_contaminant'].sum()):,} "
          f"of them SIMBAD-classified as something benign, "
          f"{int((~matched).sum()):,} unmatched)")

    cols = ["source_id", "ra", "dec", "l", "b", "alpha", "alpha_upper_2sig",
            "deficit_G_mag", "A_0", "chi2_fit", "simbad_id", "otype"]
    top = clean.nlargest(20, "deficit_G_mag")[cols]
    print()
    print(top.to_string(index=False, float_format=lambda v: f"{v:9.4g}"))

    s.to_csv(cfg.RESULT_DIR / f"searchA_candidates_{args.tag}.csv", index=False)
    out = {"n_survivors": int(len(s)),
           "n_simbad_matched": int(matched.sum()),
           "n_known_contaminant": int(s["known_contaminant"].sum()),
           "n_remaining": int(len(clean)),
           "types": {str(k): int(v) for k, v in vc.head(25).items()}}
    (cfg.RESULT_DIR / f"searchA_candidates_{args.tag}.json").write_text(
        json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
