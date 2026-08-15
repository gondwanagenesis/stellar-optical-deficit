#!/usr/bin/env python
"""Is the Lindegren+2021 parallax zero point actually being applied?"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import pandas as pd
from pipeline import config as cfg
from pipeline import sample as smp

f = sorted(cfg.RAW_DIR.glob("sample_d500_p0*.parquet"))[0]
df = pd.read_parquet(f)
print(f"{f.name}: {len(df)} rows")

print("\ninputs:")
print(f"  astrometric_params_solved value counts:\n"
      f"{df['astrometric_params_solved'].value_counts().to_string()}")
for c in ["phot_g_mean_mag", "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat"]:
    s = df[c]
    print(f"  {c:28s} n_finite={s.notna().sum():6d}  "
          f"min={s.min():8.3f} med={s.median():8.3f} max={s.max():8.3f}")

zp = smp.parallax_zero_point(df)
print(f"\nzero point (mas): n_finite={np.isfinite(zp).sum()}  "
      f"n_zero={(zp == 0).sum()}")
print(f"  median = {np.median(zp)*1000:.2f} uas")
print(f"  16-84% = {np.percentile(zp,16)*1000:.2f} .. {np.percentile(zp,84)*1000:.2f} uas")
print(f"  min/max = {zp.min()*1000:.2f} / {zp.max()*1000:.2f} uas")

# direct call, no wrapper
from zero_point import zpt
zpt.load_tables()
g = np.clip(df["phot_g_mean_mag"].to_numpy(float), 6.0, 21.0)
nu = df["nu_eff_used_in_astrometry"].to_numpy(float)
ps = df["pseudocolour"].to_numpy(float)
el = df["ecl_lat"].to_numpy(float)
sol = df["astrometric_params_solved"].to_numpy(int)
raw = zpt.get_zpt(g, np.where(np.isfinite(nu), nu, 1.45),
                  np.where(np.isfinite(ps), ps, 1.45), el, sol, _warnings=False)
raw = np.asarray(raw, float)
print(f"\nraw zpt.get_zpt output: n_finite={np.isfinite(raw).sum()}  "
      f"median={np.nanmedian(raw):.3f}  units as returned")
print(f"  16-84% = {np.nanpercentile(raw,16):.3f} .. {np.nanpercentile(raw,84):.3f}")

plx = df["parallax"].to_numpy(float)
print(f"\nparallax median {np.median(plx):.4f} mas")
print(f"implied distance-modulus shift from the zero point: "
      f"{np.median(5*np.log10(plx/(plx-zp)))*1:.5f} mag")
