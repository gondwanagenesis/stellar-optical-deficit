#!/usr/bin/env python
"""Print the injection-recovery summary."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pandas as pd
from pipeline import config as cfg

tag = sys.argv[1] if len(sys.argv) > 1 else "primary"
s = pd.read_csv(cfg.RESULT_DIR / f"injection_{tag}.csv")

print("UNIFORM  (self-calibrated recovery must be identically zero)")
u = s[s["mode"] == "uniform"].sort_values("f_injected")
print(u[["f_injected", "delta_injected", "n_stars", "mean_residual",
         "n_pos_5sig", "anchored_f"]].to_string(
    index=False, float_format=lambda v: f"{v:14.8f}"))

print("\nSPARSE  (recovery lives in the tail)")
sp = s[s["mode"] == "sparse"].sort_values(["f_injected", "p_injected"])
print(sp[["f_injected", "p_injected", "p_recovered", "p_recovered_std",
          "p_recovered_over_injected", "n_pos_5sig"]].to_string(
    index=False, float_format=lambda v: f"{v:14.7f}"))

thr = cfg.RESULT_DIR / f"injection_threshold_{tag}.csv"
if thr.exists():
    print("\nSmallest recoverable p (>3 sigma)")
    print(pd.read_csv(thr).to_string(index=False))
