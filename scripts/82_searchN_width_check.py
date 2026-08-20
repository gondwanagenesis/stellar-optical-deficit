#!/usr/bin/env python
"""Is Search N's null an artefact of a fixed threshold on unequal widths?

Search N thresholds both populations at 3 x sigma of the LOW-RUWE reference
fit (0.530 mag). But the pristine high-RUWE residuals are 1.65x wider. At a
fixed absolute cut, a wider symmetric distribution necessarily populates both
tails and so looks LESS one-sided -- which is the direction that produced the
null. This re-runs the same counting with each population thresholded at
3 x its OWN robust sigma, so one-sidedness is measured at matched
significance rather than matched magnitude.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from pipeline import config as cfg
from pipeline import statistics as st

ref = pd.read_parquet(cfg.DERIVED_DIR / "primary_resid.parquet",
                      columns=["residual", "ruwe", "cstar_nsigma"])
ok = (ref["ruwe"] < 1.4) & (ref["cstar_nsigma"].abs() < 3) & ref["residual"].notna()
r_ref = ref.loc[ok, "residual"].to_numpy(float)

# The 6,529-row candidate CSV that script 66 emitted alongside these numbers
# has been deleted rather than kept: it is the dim tail at what turned out to
# be a 1.81-sigma cut, so it is not a candidate list in any useful sense and
# leaving it on disk invites a future run to cite it as one.

s_ref = st.robust_sigma(r_ref)
# pristine high-RUWE numbers as reported by script 66
s_pri = 0.2921
print(f"robust sigma  low-RUWE reference : {s_ref:.4f}")
print(f"robust sigma  pristine high-RUWE : {s_pri:.4f}  "
      f"({s_pri/s_ref:.2f}x wider)")
print(f"fixed threshold used by Search N : 0.530 mag")
print(f"  = {0.530/s_ref:.2f} sigma for the reference")
print(f"  = {0.530/s_pri:.2f} sigma for the pristine high-RUWE sample\n")

med = np.median(r_ref)
for k in (1.8, 3.0):
    thr = k * s_ref
    nd = int((r_ref > med + thr).sum()); nb = int((r_ref < med - thr).sum())
    print(f"reference at {k:.1f} x its own sigma ({thr:.3f} mag): "
          f"dim={nd} bright={nb} ratio={nd/max(nb,1):.2f}")
print()
print("Search N compared the pristine sample at 1.81 sigma against the "
      "reference at 3.00 sigma. The line above gives the reference at the "
      "SAME significance the pristine sample was actually cut at.")
