# DECISIONS

Judgement calls not specified in the brief, with reasoning. Anything here that
could plausibly move a headline result by more than a factor of two is flagged
**[MATERIAL]**.

---

## D1. Run the pipeline under WSL, not native Windows

`healpy` has no Windows wheel (it needs `cfitsio`, and the pkg-config fallback
fails), and `dustmaps` requires `healpy` for every 3D map. Rather than shim
`healpy` onto `astropy-healpix`, the pipeline runs under the machine's
`kali-linux` WSL distribution (Python 3.13, 8 cores, 15 GB RAM). Code and
results live on `/mnt/c/...` so they stay visible from Windows; the multi-GB
dust maps live on native ext4 (`~/dustmaps_data`) because DrvFs is roughly an
order of magnitude slower for the random HDF5/FITS reads the maps do.

`run.sh` wraps the venv interpreter so Windows-side callers never fight quoting.

---

## D1a. Heavy stages must run serially — 15 GB is the binding constraint

The 3.9M-star sample is ~1.3 GB on disk and several GB in memory once derived
columns exist, and the Edenhofer map adds ~3 GB more when its integrated form
is built. Running the analysis chain and the distance trade study concurrently
drove available memory to 2 GB and risked the OOM killer taking whichever
process happened to allocate next — including a 25-minute cross-validation run.

Practical rule: one memory-heavy stage at a time. TAP pulls (network-bound) can
overlap with local compute; two dust-map consumers cannot.

---

## D2. Sky partitioning on `source_id` ranges

Gaia `source_id` encodes a level-12 HEALPix index in its high bits, so
`source_id BETWEEN lo AND hi` is simultaneously a primary-key range scan and a
contiguous sky patch. Partitioning this way gets indexed queries and free sky
locality — which the hemisphere and crowding null tests need — from the same
clause. Level 2 (nside 4, 192 partitions) puts ~28k rows (500 pc) or ~90k rows
(1250 pc) in each chunk, far below the 3e6-row archive cap.

A chunk that returns exactly the row cap is treated as a **failure**, not a
success, because the ESA archive truncates silently.

---

## D3. Pull the wide (1250 pc) sample as a superset — revised

Measured counts: 5.42M stars at 500 pc, 17.3M at 1250 pc (only 3.2x). The
original plan was to pull the wide superset once and subset client-side.

**Revised after measuring throughput.** The ESA archive throttles concurrent
anonymous async jobs: 2 workers completed a 120k-row partition every ~13 min,
but 5 workers completed *zero* partitions in 45 minutes. At the achievable
2–3 worker rate the wide pull is ~20 h versus ~7 h for the 500 pc primary.
The 500 pc primary is therefore the critical path, and the wide sample is
pulled for a subset of partitions for the distance trade study only.

---

## D4. Which cuts run server-side

Everything that does not require extinction. The observed-colour and
observed-`M_Ks` pre-boxes are widened by the largest reddening possible inside
the distance limit, so the server-side filter provably cannot clip a star the
final extinction-corrected box would have kept.

---

## D5. 2MASS quality flags: only `ph_qual` is available

`gaiadr1.tmass_original_valid` is a reduced-column mirror carrying `ph_qual`
but **not** `cc_flg`, `rd_flg` or `bl_flg`. Contamination screening therefore
uses `ph_qual[Ks] == 'A'`, `ks_msigcom < 0.05`, and the Gaia-side cross-match
uniqueness flags (`number_of_mates = 0`, `number_of_neighbours = 1`). The full
2MASS PSC flags are pulled from VizieR only for step-7 candidates, where the
per-object cost is affordable.

---

## D6. Extinction: two maps × two band laws **[MATERIAL]**

Band law is *not* a nuisance to average over. The Gaia/Fitzpatrick-2019
coefficient `k_G = A_G/A_0` runs from 1.00 at `(BP-RP)_0 = 0` to 0.65 at
`(BP-RP)_0 = 3`. A constant `A_G/A_0` would imprint a colour-dependent
residual across the main sequence — which is precisely the shape a
mass-dependent harvesting signal would have. The colour-dependent law is
primary; Wang & Chen (2019) constant ratios are the systematic variant.

The two laws also disagree by a factor 2.5 in the near-infrared
(`A_Ks/A_0 = 0.194` Fitz19 vs `A_Ks/A_V = 0.078` Wang & Chen). At `A_0 = 0.3`
that is a 0.035 mag shift in `M_Ks`, far above the target sensitivity. This is
a real published disagreement about the NIR slope of the extinction curve and
it is propagated, not averaged away.

`A_G < 0.5` for the primary sample keeps the ~5% band-law uncertainty on `A_G`
below 0.025 mag.

---

## D7. Bayestar19 is unavailable; map substitutions **[MATERIAL]**

Harvard Dataverse returns HTTP 202 with an empty body for *every* path
(`/`, `/api/info/version`, dataset pages, DOI redirects), from both a plain
client and a browser user-agent. Bayestar19 and SFD98 are therefore
unreachable and `dustmaps.fetch()` dies with a `JSONDecodeError` when it tries
to parse the gate page as JSON. Zenodo-hosted maps are unaffected.

Consequence for the "two independent maps" requirement:

- **Primary:** Edenhofer et al. 2024 (A&A 685, A82) — full sky within 1.25 kpc,
  Gaia DR3 XP spectra, nside 256. Downloaded (3.1 GB).
- **Second:** Leike et al. 2020 (Zenodo) if it fetches. Caveat stated
  explicitly in LIMITATIONS.md: Leike2020 and Edenhofer2023 are *not* fully
  independent — same group, same variational-inference lineage, overlapping
  input data. Their difference bounds reconstruction/resolution systematics,
  not the systematics common to the whole Gaia-based approach.
- **Third, genuinely independent:** Gaia's own `ag_gspphot` — a *per-star*
  extinction from BP/RP spectra, not a map at all. It has completely different
  systematic structure (no spatial smoothing; model-dependent instead). This
  substitutes for the lost photometric-map independence and is arguably a
  better independence test than Bayestar19 would have been.

Bayestar19 will be retried; if Dataverse recovers it is added as a fourth.

**`Edenhofer2023Query` must be constructed with `integrated=True`.** The
default returns extinction *density per parsec* (~1e-4), which silently
produces an almost-zero `A_0`. This bug was caught in the first end-to-end
test only because the median `A_0` printed as 0.000 — it is invisible in any
summary that merely looks "small but plausible".

---

## D8. Parallax zero-point correction is mandatory **[MATERIAL]**

Lindegren et al. (2021, A&A 649, A4). The zero point is −20 to −40 µas, i.e.
1–2% in distance at parallax 2 mas, i.e. 0.02–0.04 mag of distance modulus —
an order of magnitude above the sensitivity being chased.

The published correction is calibrated for 6 < G < 21; the sample floor is
G = 4, so for G < 6 the correction is held at its G = 6 value rather than
extrapolating a fitted spline outside its support. The affected fraction is
reported in the cut flow.

Note the partial cancellation: a distance-modulus error moves `M_G` and `M_Ks`
together, so it enters the residual only through `(1 − dM_G/dM_Ks)`. On the
lower main sequence `dM_G/dM_Ks ≈ 1.1–1.4`, so the cancellation is strong but
incomplete, and the sign of the leakage flips across the colour range.

### D8a. The zero point was applied 1000× too small (found and fixed)

`gaiadr3-zeropoint`'s `zpt.get_zpt()` returns the offset in **milliarcseconds**
— the same units as `parallax` — not microarcseconds. The first implementation
divided by 1000, reducing the correction to −0.04 µas: numerically present,
physically absent.

The failure mode is the dangerous kind. Nothing raised, no NaNs appeared, and
the diagnostic line printed `median -0.0 uas`, which reads as "the zero point
is negligible here" rather than as "this is broken". It was caught only because
that printed value contradicted the literature figure quoted two lines above it
in the same docstring.

Impact: a 1.4% distance-scale error at the sample's median parallax of 3.05 mas,
i.e. **−0.031 mag of distance modulus**, entering the residual through
`(1 − s) ≈ −0.25` as ~0.008 mag — the same order as the measured systematic
floor, and distance-dependent, so it would have leaked directly into the
near/far null split.

Guarded now: `parallax_zero_point()` raises if the median falls outside
(−150, −1) µas, and `tests/test_zeropoint.py` checks the magnitude, the units
ratio against the parallax, the sign of the distance-modulus shift, and that a
deliberately rescaled package is rejected.

**All interim (validation-sample) numbers in this repository's history were
computed with the buggy value.** The full-sample results use the fix.

---

## D9. Fiducial: cubic B-spline at quantile knots, not a GP

An exact GP on ~10^6–10^7 points is intractable, and a sparse GP introduces an
inducing-point systematic that would be indistinguishable from the signal.
A cubic B-spline with knots at *quantiles* of `M_Ks` (density-adaptive, so no
knot interval is starved) tensored with a metallicity polynomial gives a
flexible fit whose complexity is a single integer chosen by 5-fold CV.

CV loss is **Huber, not squared**: the main sequence has a genuine equal-mass
binary sequence 0.75 mag above it, and a squared loss would select the knot
count that best chases it.

### D9a. The CV grid was truncated on monotonicity evidence **[MATERIAL]**

Full-sample 5-fold CV on a 200k subsample (`results/cv_primary.csv`):

| knots | deg 1 | deg 2 | deg 3 |
|---|---|---|---|
| **6** | **0.777804** | 0.796967 | 0.802826 |
| 8 | 0.779414 | 0.799195 | 0.803106 |
| 10 | 0.780631 | 0.799219 | 0.804972 |
| 12 | 0.780968 | 0.799725 | 0.805281 |
| 16 | 0.781087 | — | — |

The loss rises monotonically in *both* directions across 13 evaluated points;
the minimum sits at the grid corner (6 interior knots, metallicity degree 1).
Evaluating the remaining grid (20, 25, 30, 40 knots) would have cost ~3 hours
per variant, because IRLS cost grows steeply with knot count, in order to
confirm a trend that is already unambiguous over a factor of 2.7 in knots.

The grid was therefore truncated at 16 and the corner selection adopted. This
is a decision about computational scope taken **before unblinding** and with
the evidence recorded above; it is flagged MATERIAL because a genuinely
non-monotonic loss surface would invalidate it. The same complexity was
adopted for variant B (which adds the optical-colour covariate) rather than
re-running its own grid; the independent CV on the validation sample also
selected (6, 1) for both variants, which is the check that makes this safe.

That 6 interior knots wins is itself informative: the lower main sequence in
`M_G` versus `M_Ks` is very nearly a smooth low-order curve, and the extra
freedom of 20+ knots buys nothing but variance.

Fitting is Huber IRLS for the same reason. Because binaries are *over*-luminous,
any residual pull is toward brighter residuals — it suppresses a deficit
signal rather than creating one.

---

## D10. Metallicity control uses GSP-Phot, with variants **[MATERIAL]**

GSP-Phot `[M/H]` is available for ~85% of the sample (GSP-Spec for only 0.7%).
But GSP-Phot metallicity is derived from *the same BP/RP photometry* whose
residuals we are measuring, so it can in principle absorb a real deficit by
reassigning metallicity. Three variants are run and the spread reported:

- **Primary:** GSP-Phot `[M/H]`.
- **Variant A:** no metallicity term at all.
- **Variant B:** GSP-Spec `[M/H]` (spectroscopic, independent of BP/RP
  photometry) on the 0.7% subsample.

---

## D11. The self-calibrated test is blind to a uniform offset **[MATERIAL]**

This corrects the framing in the brief and is the single most important
methodological point in the project.

The fiducial relation is fit *to the data*. A harvesting fraction `f` applied
uniformly to every star shifts `M_G` by a constant `−2.5 log10(1−f)` and leaves
`M_Ks` untouched; the fit absorbs that constant into its intercept and the
residuals return to zero. **Sensitivity of the self-calibrated test to a
uniform population offset is exactly zero, not `sigma/sqrt(N)`.**

What the empirical test *can* see:
- a sparse subpopulation offset from the rest (individual outliers, and a
  shifted mean if `p·f` is large enough);
- a deficit whose dependence on `M_Ks`, colour or metallicity the spline basis
  cannot absorb.

Bounding the *uniform* case requires an external absolute anchor — stellar
models, or dynamical/asteroseismic masses — and is then limited by that
anchor's systematics rather than by counting statistics. Both analyses are
run and reported separately. The step-4 injection tests include a uniform
injection precisely to demonstrate the zero recovery rather than assert it.

---

## D12. Structured, not random, subsamples measure the floor **[MATERIAL]**

Random subsamples of a single catalogue cannot manufacture a systematic, so
their mean-residual scatter tracks `sigma/sqrt(N)` all the way to `N_total`.
Running that test and reporting the (inevitable) agreement would be a way of
fooling oneself.

The floor is measured with **structured** groups that share a systematic:
contiguous HEALPix patches over a range of nside, and quantile bins in
extinction, apparent magnitude, distance and crowding. Group-mean scatter
follows `sigma/sqrt(N)` while statistics dominate and flattens at the coherent
level. That plateau is the reported sensitivity. The random-subsample curve is
reported alongside as an explicit control.

---

## D13. Boundaries chosen by inspection

- Main-sequence box `3.0 < M_Ks < 8.0`: brighter than 3.0 the turnoff and
  subgiants contaminate; fainter than 8.0 2MASS becomes incomplete at 500 pc.
- `0.7 < (BP−RP)_0 < 3.6`: roughly G2V to M4V, the range over which the
  relation is single-valued and densely sampled.
- Outlier threshold 5 sigma for individual candidates — chosen before
  unblinding, and with ~5e6 stars a 5-sigma Gaussian cut yields ~3 expected
  false positives, few enough to follow up individually.
- Crowding proxy is sample stars per nside-64 pixel (~0.84 deg²). This is the
  density of the *selected* sample, not true stellar density; a monotone proxy
  is all the split test requires.
