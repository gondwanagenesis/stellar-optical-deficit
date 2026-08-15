# LIMITATIONS

Ordered roughly by how much they cost. The first two are structural: they
change what the measurement *is*, not merely how precise it is.

---

## 1. The self-calibrated test cannot see a uniform deficit at all

**Not "is imprecise". Cannot see it. Sensitivity exactly zero.**

The fiducial relation `M_G = f(M_Ks, [M/H], …)` is fit to the same stars being
searched. A harvesting fraction applied uniformly to the whole population
shifts `M_G` by the constant `−2.5 log10(1−f)` and leaves `M_Ks` alone; a spline
with a free intercept absorbs that constant exactly, and the residuals return
to zero.

This is not a theoretical worry — it is measured. In the injection tests
(`results/injection_*.csv`), uniform injections spanning `f = 1e-5` to `0.1`
all returned an identical mean residual to six decimal places. Four orders of
magnitude of injected signal, zero response.

The brief's proposed statistic — "a nonzero MEAN `Delta_M_G` across the
population" — is therefore not a measurable quantity in a self-calibrated
analysis, and no amount of `N` helps. What *is* measurable is the **shape** of
the residual distribution at fixed `M_Ks`: a one-sided fainter tail, which
sparse harvesting produces and a uniform offset does not.

Bounding the uniform case requires an external absolute anchor (stellar models,
or dynamical masses). That analysis is then limited by the anchor's systematic
error — of order 0.05 mag for isochrones, i.e. `f ≈ 5%` — not by counting
statistics.

---

## 2. The flat-absorber caveat is real, but it is not where the brief puts it

The brief states that a spectrally flat absorber "suppresses numerator and
denominator together and is invisible to this test". Worked through, that is
not correct, and the actual blind spot is somewhere else.

The residual for an attenuated star is

```
r = dm_G − s · dm_Ks        where s = dM_G/dM_Ks along the fiducial
```

and **`s` is not 1**. Measured on our own fit: `s = 1.25` (16–84%: 1.04–1.45)
with the NIR colour control, `s = 1.52` (1.23–1.94) without. `G` is a bluer and
far more temperature-sensitive band than `Ks`, so `M_G` runs faster along the
main sequence. Consequences:

| absorber | `dm_G/dm_Ks` | leverage | what it looks like |
|---|---|---|---|
| grey (`α = 0`) | 1.00 | **−0.25** | *over*-luminous, not invisible |
| `α = 0.19` | 1.25 | **0.00** | genuinely invisible — the true blind spot |
| `α = 1` | 3.3 | +0.63 | deficit-like |
| interstellar dust (`α ≈ 2`) | 11.6 | +0.89 | deficit-like |

So:

- A **grey absorber is not invisible** — it produces a residual of the opposite
  sign, about a quarter the size of the naive expectation. A population of grey
  absorbers would show up as an *over-luminous* tail.
- The **true blind spot is at `α ≈ 0.19`**, an absorber that is already
  moderately wavelength-selective. Answering the brief's question directly:
  the test has leverage above about `α ≈ 0.2`, reaching 60% of the naive
  sensitivity by `α ≈ 1` and 90% by `α ≈ 2`.
- The blind `α` moves with `s`, so it varies across the sample: `α_blind` runs
  from 0.04 to 0.32 over the 16–84% range of the fitted slope. There is no
  single blind wavelength dependence for the whole population.
- **Interstellar dust sits at `α ≈ 2`, deep in the deficit-like regime.** Under-
  corrected reddening does not merely add noise; it mimics the signal with 89%
  efficiency. This is why extinction dominates the systematic budget.

See `figures/F5_spectral_leverage.png` and `results/spectral_leverage.json`.

---

## 3. Unresolved binaries mimic the signal — the brief has the sign backwards

The brief says an unresolved companion "makes a star OVERluminous, which is the
opposite of the signal, so contamination biases you conservative", and asked
for this to be verified rather than assumed. Verified, and it is false.

Because `s > 1`, and because a cool companion adds far more `Ks` flux than `G`
flux, both terms in `r = dm_G − s·dm_Ks` push the residual **positive** —
fainter than fiducial, which is the harvesting signature.

- Equal-mass binary: `r = 0.753 (s − 1) > 0`.
- Cool companion: worse, reaching several tenths of a magnitude.
- A primary at `M_Ks = 4.0` with a companion 1–2 mag fainter in `Ks` mimics
  `f > 0.1`.

Tests in `tests/test_binary_bias.py`. The consequence is that the RUWE cut is
load-bearing in the *opposite* sense to the one stated: it is not a
conservative safety margin, it is the primary defence against a contaminant
that manufactures false positives. Residual undetected binarity — companions
too close, too faint, or at unlucky orbital phase for RUWE to flag — is an
irreducible background.

---

## 4. Extinction is the dominant measured systematic

The null tests put numbers on this. Splitting the sample by extinction
quartile gives a mean-residual difference of **−0.020 mag** at 12 sigma
(`results/nulls_*_splits.csv`); after the NIR control, the colour split is
worse still. For comparison the naive `sigma/sqrt(N)` on the same sample is
`4.4e-4` mag.

Contributing pieces:

- **Map choice.** Bayestar19 was unobtainable — Harvard Dataverse returned HTTP
  202 with an empty body for every path throughout this work, so `dustmaps`
  could not fetch it (or SFD98). The primary map is Edenhofer et al. 2024.
  The intended independent cross-check is therefore weakened; see §8.
- **Band law.** Fitzpatrick 2019 (the official Gaia coefficients) versus
  Wang & Chen 2019 differ by a **factor 2.5** in `A_Ks/A_0` (0.194 vs 0.078).
  That is a real, published disagreement about the near-infrared slope of the
  extinction curve, not a rounding difference. The paired test gives an RMS
  per-star difference of 0.022 mag.
- **Per-star vs map.** Swapping the map for Gaia's own GSP-Phot `A_G` changes
  individual residuals with an RMS of **0.11 mag** — comparable to the entire
  intrinsic scatter.
- **`E → A_0` normalisation.** The Edenhofer/ZGR23 unitless extinction is
  converted with `A_V = 2.8 E`. An error in that single scalar is a global
  multiplicative extinction error and appears directly in the extinction-split
  null test.

---

## 5. Metallicity control is partly circular

GSP-Phot `[M/H]` is available for ~87% of the sample (GSP-Spec for 0.7%), but
it is derived from **the same BP/RP photometry** whose residuals are being
measured. In principle it can absorb a real deficit by reassigning metallicity.
Three variants are run (GSP-Phot, none, GSP-Spec) and the spread reported, but
none of them is both high-coverage and independent. The 0.7% GSP-Spec
subsample is the only clean control and is too small to constrain the full
sample.

Dropping the metallicity term entirely inflates the scatter from 0.096 to
0.168 mag, so the term is doing real work — which is exactly why its
circularity matters.

---

## 6. Selection function

- `parallax_over_error > 20` and `ruwe < 1.4` select against exactly the stars
  most likely to be interesting (distant, and astrometrically perturbed).
- The 2MASS `Ks` depth, not Gaia, sets the faint limit: `ks_msigcom < 0.05`
  binds well before `G < 19`. The sample is effectively `Ks`-limited, so at
  fixed `M_Ks` the selection is *approximately* independent of `M_G` — which is
  fortunate, since a `G`-limited sample would preferentially drop exactly the
  deficit stars being searched for. "Approximately" because the `G`-side quality
  cuts (flux SNR, `C*` locus) do select on `G`.
- The near/far and bright/faint null splits both fail at >10 sigma, so some
  distance- and magnitude-dependent selection or extinction residual is
  present and is not being fully modelled.

---

## 7. The mass anchor is a different population from the science sample

Dynamical masses come from binaries. The science sample removes binaries by
construction (RUWE, `non_single_star`). The flat-absorber bound from the mass
anchor is therefore an argument by analogy between two disjoint populations,
not a direct measurement on the stars being searched. Asteroseismic masses
would be independent but exist for only a handful of main-sequence dwarfs.

---

## 8. Only one genuinely independent 3D dust map was obtainable

The brief asked for at least two independent maps with their difference
treated as a systematic. Harvard Dataverse was unreachable for the duration
(HTTP 202, empty body, every path, plain client and browser user-agent), which
removed Bayestar19 and SFD98. What was used instead:

- **Edenhofer et al. 2024** (primary, from Zenodo).
- **Gaia GSP-Phot per-star `A_G`** — not a map at all, and therefore
  genuinely independent in its systematics (no spatial smoothing; model-
  dependent instead). This is a real independence test and is reported.
- Leike et al. 2020 was attempted as a second map, with the caveat that it
  shares methodology and input data with Edenhofer2023 and so bounds
  reconstruction systematics only.

The map-to-map term in the systematic budget is consequently less well
constrained than intended, and is probably **under**-estimated.

---

## 9. Astrophysical backgrounds that genuinely dim optical light

These are not errors; they are real stars that look like the signal.

- **Starspots.** BY Dra and RS CVn rotators dim the optical continuum at
  roughly constant `Ks`. This is the harvesting signature, produced by
  magnetic activity. Spot coverage of 10–30% is common on active M dwarfs.
- **Young stellar objects and debris discs.** Circumstellar material dims the
  star and reddens it. The SIMBAD follow-up on the test sample found the
  classified positive outliers dominated by `Y*?` (candidate YSO) and `Em*`
  (emission-line) types.
- **Edge-on discs**, **background reddening pockets** below the resolution of
  any 3D dust map (the maps have ~1 pc voxels at best; a dense clump in front
  of one star is unresolved).

---

## 10. What "no mid-IR excess" can and cannot mean

The beaming argument that motivates the optical channel predicts candidates
*without* mid-infrared excess. But absence of excess is only informative if the
mid-infrared was actually measured. In the test-sample follow-up, 34 of 60
positive outliers had **no AllWISE photometry at all** — and outliers
preferentially lack clean WISE matches precisely because they sit in confused
regions, which is also why they are outliers. Counting "no excess detected" as
"beaming-consistent" would have promoted 35 candidates instead of 1.

---

## 11. Scope

- Lower main sequence only (`3 < M_Ks < 8`, `0.7 < (BP−RP)_0 < 3.6`), roughly
  G2V–M4V. Nothing here constrains harvesting around hotter, more luminous
  stars, which are the more attractive engineering targets.
- Within 500 pc. This is `~10^-7` of the Galaxy by volume.
- Gaia DR3 photometry, single epoch-averaged. No time-domain information, so a
  variable deficit is not searched for.
- The `A_G < 0.5` cut removes the highest-extinction lines of sight, which are
  also where a distant harvested population would most plausibly hide.
