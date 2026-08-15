#!/usr/bin/env python
"""How exactly does a uniform offset get absorbed? Quantify the leakage."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from pipeline import fiducial as fid
from pipeline import statistics as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests"))
from test_fiducial import synth                     # noqa: E402

m_ks, mh, m_g, _ = synth()
covs = [(mh, 1)]
base = fid.fit_fiducial(m_ks, covs, m_g, 10)
r_base = base.residuals(m_ks, covs, m_g)
print(f"baseline: mean={np.mean(r_base):+.3e} sigma={st.robust_sigma(r_base):.5f}")

print("\nuniform injection leakage into residuals:")
for f in [1e-5, 1e-3, 1e-2, 0.1, 0.5]:
    d = st.delta_mag(f)
    fit = fid.fit_fiducial(m_ks, covs, m_g + d, 10)
    r = fit.residuals(m_ks, covs, m_g + d)
    print(f"  f={f:<8g} delta={d:8.5f}  max|r-r_base|={np.abs(r-r_base).max():.3e}"
          f"  mean shift={np.mean(r)-np.mean(r_base):+.3e}"
          f"  leak fraction={(np.mean(r)-np.mean(r_base))/d:.3e}")

print("\nsparse injection p=0.01 f=0.5:")
rng = np.random.default_rng(3)
m_ks2, mh2, m_g2, _ = synth(n=60_000, seed=7)
covs2 = [(mh2, 1)]
b2 = fid.fit_fiducial(m_ks2, covs2, m_g2, 10)
rb = b2.residuals(m_ks2, covs2, m_g2)
sb = st.robust_sigma(rb)
npb = int((rb > np.median(rb) + 5 * sb).sum())
mask = rng.random(len(m_g2)) < 0.01
inj = m_g2.copy(); inj[mask] += st.delta_mag(0.5)
f2 = fid.fit_fiducial(m_ks2, covs2, inj, 10)
r2 = f2.residuals(m_ks2, covs2, inj)
s2 = st.robust_sigma(r2)
np2 = int((r2 > np.median(r2) + 5 * s2).sum())
print(f"  n_pos 5sig: {npb} -> {np2}")
print(f"  mean shift: {np.mean(r2)-np.mean(rb):+.5f}   p*delta = {0.01*st.delta_mag(0.5):.5f}")
print(f"  ratio mean_shift/(p*delta) = {(np.mean(r2)-np.mean(rb))/(0.01*st.delta_mag(0.5)):.3f}")
print(f"  sigma: {sb:.5f} -> {s2:.5f}")
