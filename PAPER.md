# A search for the optical deficit signature of stellar energy harvesting in 3.3 million Gaia DR3 main-sequence stars

**Status:** working paper, self-published. Not peer reviewed. All code, data
provenance and null tests are in this repository; every number below is
regenerable by a numbered script.

---

## Abstract

Essentially every survey for Dyson-type stellar energy harvesting has searched
for mid-infrared **excess** from reprocessed waste heat. That channel is
evadable: a civilisation that beams its waste heat into a solid angle Ω is
detectable by only Ω/4π of observers. The optical **deficit** — starlight that
is intercepted and never departs — is visible from every direction and cannot
be evaded by beaming. Its cost is sensitivity, and that cost has not previously
been measured.

We measure it. From 4,809,840 Gaia DR3 × 2MASS sources within 500 pc we build a
sample of 3,884,167 lower-main-sequence stars (3,321,566 with GSP-Phot
metallicities) and search for a deficit in absolute *G* at fixed absolute
*K*<sub>s</sub>, after 3D-extinction and metallicity control.

The intrinsic main-sequence scatter is **0.0886 mag**, giving a naive
sensitivity of σ/√N = 5.4 × 10⁻⁵ mag. **That number is wrong by three orders of
magnitude.** Structured null tests — sky patches, extinction bins, magnitude
bins — plateau at 1.6–2.3 × 10⁻² mag while random subsamples track σ/√N to
within 1% out to N = 3 × 10⁶. The measured systematic floor is **957× the naive
expectation**, and it is dominated by an *astrophysical* term (0.068 mag) that
survives where there is no dust.

We report **zero candidates**. Of the 300 strongest optical-deficit outliers,
none has measured-and-normal mid-infrared photometry. Our best constraint,
from 8,844 co-natal wide pairs where the background is measurably symmetric,
is **p < 4.8 × 10⁻⁴** of stars intercepting ≥52% of their optical output.
(An earlier version used an asymmetry estimator that injection shows is blind
to this signal; see §5.6.)

This is *not* the strongest published limit — [Suazo et al.
(2022)](https://doi.org/10.1093/mnras/stac1789) reach 1.9 × 10⁻⁴ at γ ≥ 0.5 by
using the infrared re-emission as well. What we can claim uniquely is a limit
on the class their estimator is not designed to constrain: requiring both pair
components to have a **measured bare photosphere** in W1 and W2 gives
**p < 4.9 × 10⁻⁴ at f ≥ 0.45** for stars with a large optical deficit *and no
warm re-emission* — the beaming-consistent, cold, or non-thermally-exporting
case that motivates the deficit channel in the first place.

Combining that with the infrared limit gives what we believe is the **first
bound on stellar energy harvesting that does not assume how the waste heat is
disposed of**: writing p_total = p_iso + p_dark for the isotropic-warm and
beamed/cold/non-thermal populations, Suazo et al. constrain the former and §5.7
constrains the latter, so **p_total < 6.8 × 10⁻⁴** — fewer than 1 in 1,473
nearby lower-main-sequence stars intercepts ≥45% of its optical output by any
means. We further show that a 10× larger wide-pair sample, well within the
existing El-Badry et al. catalogue, would push the deficit limit below the
infrared limit and make the *beamed fraction* β = p_dark/p_total a measurable
quantity for the first time.

Two further channels close blind spots the main estimator cannot reach. A
**dynamical-mass** search — regressing M_G on log M_dyn for 88,149 Gaia binary
systems — is indifferent to spectral slope and does not use 2MASS as its
anchor, giving **p < 2.2 × 10⁻⁴ at f ≥ 0.73 for an absorber of any spectral
slope**, with **zero positive outliers** among the 13,482 systems with faint
secondaries. This is the only constraint here covering the grey and α ≈ 0.19
cases. A **kinematic** search in velocity space rather than position space
returned a 4301σ excess that dissolves to ~1.4× once bound pairs, |v| and pair
separation are controlled: anomalies are kinematically cold (median |v| = 22.8
against 35.1 km/s) because the tail is enriched in young stars, which cluster
in comoving associations for mundane reasons.

Three results correct assumptions we began with, and we consider them the more
durable contribution:

1. **A uniform deficit is unmeasurable by any self-calibrated version of this
   method** — sensitivity exactly zero, not σ/√N, demonstrated by injection
   across four decades.
2. **Unresolved companions mimic the signal rather than suppressing it**,
   because d*M<sub>G</sub>*/d*M*<sub>Ks</sub> = 1.26 > 1. Confirmed on data.
3. **The spectral blind spot is not a grey absorber.** It lies at α ≈ 0.19;
   a grey absorber produces a residual of the *opposite* sign.

The dominant contaminant is identified mechanically: **aperture mismatch**
between Gaia (sub-arcsecond) and 2MASS (~4″). We localise it to θ ≲ 6.5″ in
wide pairs, where the contaminated-to-clean tail ratio is 50:1 and falls to
unity beyond 10″. This limits *any* optical-deficit search anchored on 2MASS.

---

## 1. Introduction

### 1.1 The infrared-excess paradigm and its blind spot

[Dyson (1960)](https://doi.org/10.1126/science.131.3414.1667) pointed out that a
civilisation intercepting a substantial fraction of its star's output must
re-radiate that energy at low temperature, producing an infrared excess. Every
major search since has followed that logic:
[Carrigan (2009)](https://doi.org/10.1088/0004-637X/698/2/2075) searched ~250,000
IRAS sources for 100–600 K blackbodies and found that the best candidates were
reddened dusty objects — heavily extinguished stars, protostars, Mira
variables, AGB stars and planetary nebulae. The Ĝ programme
([Wright et al. 2014a](https://arxiv.org/abs/1408.1133),
[2014b](https://doi.org/10.1088/0004-637X/792/1/27);
[Griffith et al. 2015](https://doi.org/10.1088/0067-0049/217/2/25)) extended
this to WISE with ~1000× better sensitivity and 5× better resolution, over both
Galactic sources and ~10⁵ galaxies. Project Hephaistos
([Suazo et al. 2022](https://doi.org/10.1093/mnras/stac1789);
[2024](https://doi.org/10.1093/mnras/stae1186)) applied the modern
Gaia DR3 × 2MASS × AllWISE combination to ~5 × 10⁶ sources, reporting seven
candidates after a filtering cascade.

The blind spot is geometric. Infrared excess is a statement about where the
waste heat *goes*. [Wright (2023)](https://doi.org/10.3847/1538-4357/acf44f)
showed thermodynamically that optimal use of mass favours small, hot spheres,
and that a sphere accepting radiation from a limited solid angle behaves
differently from the isotropic idealisation. A radiator that beams into solid
angle Ω is seen by a fraction Ω/4π of observers. Nothing in physics forbids
this, and a civilisation with any preference for not being conspicuous — or
merely one that finds it useful to aim its waste heat somewhere — is
under-counted by every infrared survey by that factor.

### 1.2 The optical deficit channel

The complementary observable is the light that never leaves.
[Annis (1999)](https://ui.adsabs.harvard.edu/abs/1999JBIS...52...33A)
originated the method at galaxy scale, looking for outliers in the
Tully–Fisher relation and the fundamental plane, adopting a criterion of ≥1.5
mag of dimming (a factor of 4) on 57 disks and 106 ellipticals, and finding
none. [Zackrisson et al. (2015)](https://doi.org/10.1088/0004-637X/810/1/23)
made this quantitative on 1,359 disk galaxies, limiting to **<3%** (0.3%
tentative) the fraction of local disk galaxies subject to galaxy-scale Dysonian
astroengineering.

A deficit cannot be beamed away. It is visible to every observer in every
direction, and it is agnostic about what happens to the energy afterwards — it
does not care whether the waste heat is radiated, beamed, stored, or used for
computation at the Landauer limit. That is the entire appeal.

The obvious extension is from galaxies to individual stars, where Gaia provides
10⁹ objects with parallaxes. The back-of-envelope is seductive: with a
main-sequence relation of width σ ≈ 0.1 mag, a single star gives a 3σ limit of
f ≈ 26% (useless), but N ~ 10⁸ stars would give σ/√N ~ 10⁻⁵ mag, i.e.
f̄ ~ 10⁻⁵ (decisive).

**The purpose of this paper is to show that this scaling fails, to measure
where and why, and to report what the channel actually delivers.**

### 1.3 This is not the first stellar optical-deficit search, and not the best limit

We state this plainly because we initially framed the work as more novel than
it is, and a reader deserves the correction up front.

**[Zackrisson et al. (2018)](https://doi.org/10.3847/1538-4357/aac386)** is the
direct methodological ancestor. They search for optically underluminous stars
by comparing spectrophotometric to trigonometric distance, since a dimmed star
has its spectrophotometric distance overestimated while its parallax does not
move: f_cov = 1 − (D_trig/D_spec)². On ~230,000 Gaia DR1 × RAVE DR5 stars
(8,441 after cuts) they find 75 stars at f_cov > 0.9, reducing to six
unexplained outliers after vetting.

Their method is sensitive **only to nearly complete spheres** — it takes
f_cov > 0.75 to produce the factor-of-two distance discrepancy the method
needs, and a 10% covering fraction gives a ~5% distance disagreement that is
undetectable. That is the same sensitivity regime we reach, from a different
observable.

Their prime candidate, TYC 6111-1162-1 (f_cov ≈ 0.77, and notably *without*
mid-infrared excess), resolved to a **single-lined spectroscopic binary** whose
unseen ~1 M☉ white-dwarf companion perturbed the Gaia parallax. We regard this
as important independent corroboration: their one candidate died of binarity,
and we find binarity to be the dominant contaminant by a completely different
route (§5.5). Two methods, two different observables, the same false-positive
population.

**[Suazo et al. (2022)](https://doi.org/10.1093/mnras/stac1789)** (Hephaistos I)
already set upper limits on partial Dyson spheres using Gaia DR2 × AllWISE with
the AGENT formalism, modelling the system as
L_ν = (1−γ)L_ν,⋆ + γ·BB_ν(T_DS) — that is, **using the optical obscuration and
the infrared re-emission together**. Within 100 pc at T_DS = 300 K they obtain

| covering fraction | upper limit on fraction of stars |
|---|---|
| γ ≥ 0.1 | 6.6 × 10⁻³ |
| γ ≥ 0.5 | **1.9 × 10⁻⁴** |
| γ ≥ 0.9 | 1.8 × 10⁻⁵ |

**Their limit at γ ≥ 0.5 is roughly 2.4× stronger than our best
(4.6 × 10⁻⁴ at f ≥ 0.46).** We do not set a record and we do not claim one.

### 1.4 What is actually new here

Given the above, the contribution of this work is not the limit. It is:

1. **A measurement of the systematic floor of the deficit channel.** No prior
   work quantifies where σ/√N stops applying. We find it stops **957–1376×**
   above the naive value, and we identify what sets it.
2. **Identification and localisation of the dominant contaminant** as
   Gaia/2MASS aperture mismatch, pinned to θ ≲ 6.5″ (§5.5). This is structural
   and limits any 2MASS-anchored search, including future ones.
3. **Three methodological corrections** that we believe are not stated
   correctly elsewhere: a uniform deficit is unmeasurable by a self-calibrated
   relation (§3.2); unresolved companions mimic rather than suppress the signal
   because dM_G/dM_Ks > 1 (§5.5); and the spectral blind spot lies at α ≈ 0.19
   rather than at a grey absorber (§3.3). The last bears directly on Hephaistos'
   choice to model obscuration as grey.
4. **A differential wide-pair technique** with a *measured* symmetric
   background, which is the only estimator here that permits background
   subtraction (§5.6).

A note on method complementarity: our estimator regresses against an empirical
main-sequence relation, which is why a uniform or grey deficit is invisible to
it. Zackrisson et al. and Suazo et al. use absolute SED fitting against a
parallax, which does **not** share that blind spot — a grey absorber makes a
star underluminous in every band and is directly detectable that way. Our blind
spots and theirs are different, which is an argument for running both rather
than for preferring either.

### 1.3 Related approaches

Other lines of attack share our motivation of looking for absence rather than
emission. The VASCO project
([Villarroel et al. 2016](https://doi.org/10.3847/1538-3881/ab570f), 2022)
compares century-old photographic plates against modern surveys for sources
that have *vanished* — a limiting case of a deficit, and one with the enormous
advantage of a pre-Sputnik baseline free of satellite contamination.
[Boyajian et al. (2016)](https://doi.org/10.1093/mnras/stv1233) established
that deep, aperiodic, otherwise-unexplained dimming of a single main-sequence
star is detectable in practice, which is the time-domain analogue of what we
search for photometrically.

---

## 2. Data

### 2.1 Query and sample construction

We query the ESA Gaia archive TAP service directly, joining
`gaiadr3.gaia_source` to 2MASS via `tmass_psc_xsc_best_neighbour` →
`tmass_psc_xsc_join` → `gaiadr1.tmass_original_valid`, to AllWISE via
`allwise_best_neighbour` → `gaiadr1.allwise_original_valid`, and to
`gaiadr3.astrophysical_parameters`. The sky is partitioned by `source_id`
range, which encodes a level-12 HEALPix index in its high bits, so each chunk
is simultaneously an indexed primary-key scan and a contiguous sky patch.

192 partitions, **4,809,840 rows, zero failures**, 252 minutes at three
concurrent workers.

### 2.2 Cuts

| stage | stars | removed |
|---|---|---|
| raw joined rows (server-side ADQL cuts) | 4,879,956 | — |
| 2MASS Ks `ph_qual` = 'A' | 4,865,655 | 0.29% |
| not photometrically variable | 4,497,637 | 7.83% |
| `non_single_star` = 0 | 4,497,475 | 7.84% |
| not QSO/galaxy candidate | 4,497,459 | 7.84% |
| DSC star probability > 0.5 | 4,323,603 | 11.40% |
| \|C*\| < 3σ (Riello et al. 2021) | 4,289,991 | 12.09% |
| dust-map coverage | 4,275,099 | 0.35% |
| 10 < d < 500 pc | 4,275,099 | 0.35% |
| A_G < 0.5 | 3,925,746 | 8.49% |
| 3.0 < M_Ks < 8.0 | 3,885,617 | 9.43% |
| 0.7 < (BP−RP)₀ < 3.6 | **3,884,167** | 9.46% |
| with GSP-Phot [M/H] | **3,321,566** | |

Server-side cuts include `parallax_over_error > 20`, `ruwe < 1.4`,
`ipd_frac_multi_peak ≤ 2`, `visibility_periods_used ≥ 10`, and 1:1 cross-match
uniqueness (`number_of_mates = 0`, `number_of_neighbours = 1`).

### 2.3 Astrometry and extinction

Parallaxes are corrected using the
[Lindegren et al. (2021)](https://doi.org/10.1051/0004-6361/202039653)
zero point; the median correction is **−39.0 µas** (16–84%: −43.8 to −33.3),
i.e. a 1.3% distance-scale effect and 0.031 mag of distance modulus at our
median parallax. This is an order of magnitude above the target sensitivity and
cannot be neglected.

Extinction uses the [Edenhofer et al. (2024)](https://doi.org/10.1051/0004-6361/202347628)
3D map (full sky within 1.25 kpc, Gaia XP-based), converted to A₀ via
A_V = 2.8 E, with the official Gaia colour-dependent coefficients from the
[Fitzpatrick et al. (2019)](https://doi.org/10.3847/1538-4357/ab4c3a)
extinction law. The colour dependence is not optional: k_G = A_G/A₀ runs from
1.00 at (BP−RP)₀ = 0 to 0.65 at 3, and a constant ratio would imprint a
colour-dependent residual degenerate with a mass-dependent signal.

**Bayestar19 ([Green et al. 2019](https://doi.org/10.3847/1538-4357/ab5362))
was unobtainable.** Harvard Dataverse returned HTTP 202 with an empty body on
every path for the duration of this work. We substitute Gaia's own per-star
`ag_gspphot`, which is not a map at all and therefore has genuinely different
systematics.

---

## 3. Method

### 3.1 The observable and the signal model

For a star at distance *d*,

- M_G = G − 5 log₁₀ d + 5 − A_G
- M_Ks = Ks − 5 log₁₀ d + 5 − A_Ks

If a fraction *f* of the optical output is intercepted,
**ΔM_G = −2.5 log₁₀(1 − f)** (positive = fainter).

We fit a fiducial relation M_G = 𝔣(M_Ks, [M/H], (J−Ks)₀) as a cubic B-spline in
M_Ks at quantile knots, tensored with covariate polynomials, by Huber IRLS.
Complexity is chosen by 5-fold cross-validation; the measured loss surface is
monotone in both directions with its minimum at 6 interior knots and
metallicity degree 1 (`results/cv_primary.csv`).

Controlling on the *near-infrared* colour (J−Ks)₀ is safe and reduces the
scatter from 0.175 to 0.099 mag: an optically selective absorber leaves a NIR
colour alone, so the control cannot launder the signal. We run a second variant
that also controls on the *optical* (BP−RP)₀, which reduces the scatter to
0.028 mag but restricts sensitivity to absorbers grey across the optical.

### 3.2 What is measurable, and what is not

**A uniform deficit is not measurable.** The fiducial is fit to the same stars
being searched. A constant offset lies exactly in the span of the basis and is
absorbed. This is not a bound, it is an identity, and we verify it: injections
from f = 10⁻⁵ to 10⁻¹ return a mean residual identical to eight decimal places
(0.00629139 in every case) with an identical 5σ tail count.

What survives is the **shape** of the residual distribution at fixed M_Ks.
Harvesting adds a one-sided fainter tail.

A *sparse* deficit does move the mean, contrary to the naive expectation: the
robust fit down-weights injected outliers rather than following them, so ~77%
of the p·Δ shift survives. We nonetheless use the tail, because the mean is
where the coherent systematics live: at the measured floor a mean-based limit
at f = 0.5 gives p < 0.09, against p < 5 × 10⁻³ from the tail.

### 3.3 Spectral leverage: where the method is actually blind

The residual for an attenuated star is **r = δm_G − s · δm_Ks** with
s = dM_G/dM_Ks. The critical point is that **s ≠ 1**: we measure
s = 1.2566 (16–84%: 1.0597–1.4567), so the slope exceeds unity essentially
everywhere in the sample.

For a power-law absorber τ ∝ λ^−α:

| absorber | δm_G/δm_Ks | leverage |
|---|---|---|
| grey (α = 0) | 1.00 | **−0.257** |
| **α = 0.193** | 1.257 | **0.000 — blind spot** |
| α = 1 | 3.33 | +0.623 |
| α = 2 (dust-like) | 11.66 | +0.892 |

**A grey absorber is not invisible; it produces a residual of the opposite
sign.** The true blind spot is a moderately selective absorber at α ≈ 0.19.
Interstellar dust sits at α ≈ 2, mimicking the signal at 89% efficiency, which
is why extinction carries so much of the systematic budget.

λ_eff is SED-dependent and this matters: Gaia G spans 400–950 nm, so for a
4500 K photosphere the flux-weighted effective wavelength is 0.679 µm, not the
Vega-referenced 0.582 µm — a 17% difference propagating to ~10% in α_blind.

---

## 4. The systematic floor

**This section was completed before any signal was examined.**

### 4.1 Null splits

Every two-sided split that must return zero fails, the mildest at 7σ:

| split | difference (mag) | significance |
|---|---|---|
| **colour, blue vs red** | **−0.05205** | **−353σ** |
| extinction quartile | −0.04456 | −212σ |
| apparent G, bright vs faint | −0.02818 | −189σ |
| crowding | −0.02612 | −174σ |
| galactic latitude | +0.02730 | +173σ |
| hemisphere N vs S | −0.00550 | −37σ |
| distance, near vs far | −0.00106 | −7σ |

### 4.2 Random subsamples are a control, not a measurement

| N | RMS of subsample means | σ/√N |
|---|---|---|
| 10,000 | 0.001361 | 0.001367 |
| 1,000,000 | 0.000132 | 0.000137 |
| **3,000,000** | **0.000080** | **0.000079** |

Agreement to 1% at N = 3 × 10⁶. Random subsampling cannot manufacture a
systematic, so it reproduces the naive scaling however wrong that scaling is as
a sensitivity. **Anyone quoting σ/√N would find this reassuring and would be
wrong by three orders of magnitude.**

Structured groups on identical residuals plateau instead: apparent-G bins at
0.0229, extinction bins 0.0183, sky patches 0.0169, crowding 0.0153, distance
0.0034 mag.

### 4.3 The floor is astrophysical

Repeating the worst split *inside the lowest-extinction quartile*:

| split | full sample | lowest-A₀ quartile |
|---|---|---|
| colour | −0.05205 (353σ) | **−0.06779 (259σ) — grows** |
| distance | −0.00106 (7σ) | +0.01031 (38σ) |

The colour systematic is **larger** where there is no dust. It is real
main-sequence structure that M_Ks, [M/H] and (J−Ks) do not capture — age,
α-enhancement, rotation, activity. No better dust map moves it.

### 4.4 Budget and the number

| term | mag |
|---|---|
| **astrophysical: MS structure vs colour** | **0.06779** |
| extinction: low vs high A₀ | 0.04456 |
| photometric: bright vs faint G | 0.02818 |
| crowding | 0.02612 |
| spatial coherence plateau | 0.01694 |
| extinction fractional error × median A_G | 0.01314 |
| map vs Gaia per-star A_G | 0.00702 |
| parallax zero-point residual | 0.00195 |
| metallicity calibration | 0.00099 |
| band law (Fitz19 vs Wang & Chen 19) | 0.00005 |
| **quadrature sum** | **0.09258** |
| **largest single term** | **0.06779** |
| naive σ/√N | 0.0000544 |
| **ratio** | **1246×** |

Regressing residual on A₀ gives 0.1253 ± 0.0005 mag per unit A₀ (247σ),
implying the extinction correction is ≤19% short.

### 4.5 The distance trade

Same partitions, same complexity, varying only the distance limit:

| d (pc) | N | median A₀ | σ/√N | spatial plateau | floor/naive |
|---|---|---|---|---|---|
| 200 | 67,061 | 0.033 | 3.28e−4 | 0.0154 | 149× |
| 500 | 480,700 | 0.110 | 1.42e−4 | 0.0227 | 345× |
| 1250 | 1,383,542 | 0.166 | 7.51e−5 | 0.0235 | 536× |

21× more stars, 4.4× better σ/√N, and the ratio degrades 3.6×. **The extra
stars are bought at a worse price than they are worth.**

---

## 5. Results

### 5.1 Injection–recovery

Uniform injections return identically zero (§3.2). Sparse recovery at N = 8×10⁵:

| f | p injected | p recovered | ratio |
|---|---|---|---|
| 0.5 | 1e−2 | 9.48e−3 | 0.95 |
| 0.5 | 1e−4 | 9.84e−5 | 0.98 |
| 0.3 | 1e−3 | 6.23e−4 | 0.62 |
| 0.1 | any | consistent with zero | — |

**f = 0.1 is unrecoverable at any p**: Δ = 0.114 mag is 1.15σ and never reaches
the tail.

### 5.2 Blinding, and its failure

A secret offset (uniform ±0.05 mag, SHA-256 committed) was added to all
residuals; the analysis was frozen in git before unblinding. Revealed offset
+0.036403 mag, commitment verified.

**The blind did not work for the primary statistic**, and we report this rather
than omit it. Tail counts are defined relative to the sample median, so a
constant offset shifts residuals and median together and the counts are
identical blinded and unblinded. Only the mean and model-anchored numbers were
genuinely blinded. What protects this result is **pre-registration**: all
thresholds were fixed in `config.py` and `DECISIONS.md` before the full sample
existed. That is weaker than a working blind and it is what we have.

### 5.3 Single-star limits

| f | Δ (mag) | efficiency | p upper limit | mean f̄ |
|---|---|---|---|---|
| 0.03 | 0.033 | 0.008 | 0.770 | 0.0231 |
| 0.10 | 0.114 | 0.016 | 0.388 | 0.0388 |
| 0.30 | 0.387 | 0.169 | 0.0358 | 0.0107 |
| **0.50** | 0.753 | 0.986 | **6.13e−3** | **3.07e−3** |

The validation subset of 47,927 stars gave p_UL = 5.1 × 10⁻³; the full
3,321,566-star sample gives 6.1 × 10⁻³. **Sixty-nine times more stars, and the
limit did not improve** — the signature of a background that scales with N.

Variant B (grey-across-optical) reaches f̄ < 9.7 × 10⁻⁴ with per-star reach
f ≳ 0.2.

### 5.4 Candidates

At 5σ: **19,844 positive, 1,126 negative** — a 17.6:1 asymmetry at 129σ. Of the
300 strongest positive outliers:

| outcome | count |
|---|---|
| rejected by mid-IR excess or SIMBAD class | 22 |
| **no AllWISE photometry at all** | **278** |
| **measured-and-normal mid-IR, not a known contaminant** | **0** |

SIMBAD classes on the positive side are contaminant types (Em* ×5, Y*? ×1,
Y*O ×1); the negative control side draws ordinary stars (PM* ×17). That 278/300
lack WISE is itself diagnostic: extreme outliers sit in fields too confused for
a clean match, which is also why they are outliers.

**Zero surviving candidates.**

### 5.5 The contaminant, identified

The positive tail tracks every multiplicity proxy while the negative control
tail *shrinks*:

| proxy | positive tail | negative tail |
|---|---|---|
| RUWE | 2.04× | 0.72× |
| astrometric_excess_noise | 2.77× | 0.85× |
| **BP/RP excess factor C*** | **20.5×** | 0.38× |

C* dominating by 10× over RUWE identifies the mechanism as **aperture
mismatch**: Gaia resolves sub-arcsecond, the 2MASS beam is ~4″, so a neighbour
Gaia separates is blended into the same Ks measurement — inflating Ks and
making the star look under-luminous in G.

Wide pairs catch this in the act and localise it:

| separation | n_pos | n_neg | ratio |
|---|---|---|---|
| 3.0–6.5″ | 3 | 152 | **50.7** |
| 6.5–10.3″ | 3 | 8 | 2.7 |
| 10.3–19.3″ | 8 | 11 | 1.4 |
| 19.3–120″ | 7 | 4 | 0.57 |

**This is structural, not a data-quality problem.** It limits any
optical-deficit search anchored on 2MASS, and removing it requires a NIR anchor
of comparable angular resolution.

### 5.6 The best limit: clean wide pairs

> **CORRECTION (self-caught, post-unblinding).** An earlier version of this
> section built the limit on the *asymmetry* n_pos − n_neg of
> Δr = r_prim − r_sec, justified by the measured symmetry of the background.
> That justification is necessary but not sufficient: **the signal is symmetric
> too.** Harvesting strikes either component with equal probability, so a
> fraction p/2 of pairs shift by +Δ and p/2 by −Δ, and the asymmetry vanishes
> for the signal as well as the background.
>
> Injection confirms it (`scripts/32_pair_estimator_fix.py`): against an
> injected p = 2 × 10⁻² at f = 0.5, the asymmetry estimator responds at
> **0.02σ** while a two-sided count responds at **13.2σ**. The asymmetry
> estimator has essentially no sensitivity and the limits derived from it were
> void. Everything below uses two-sided counting.
>
> This is the kind of error that survives peer review when the estimator is
> validated only against the background and never against an injected signal.
> Our own injection framework existed and we had not pointed it at this
> estimator.

17,268 co-natal pairs from a self-join of the analysis sample (chance-alignment
contamination 0.02% by scramble test). Common-mode cancellation is real but
partial: σ(Δr) = 0.1138 against 0.1402 for none at all, so **only 34% of the
main-sequence variance is common-mode**. Two-thirds is genuinely per-star.

Beyond 10″ the tails balance, so the background is *measurably* symmetric and
can be subtracted — which the single-star test can never do.

| | single star | **clean pairs (θ > 10″)** |
|---|---|---|
| stars | 3,321,566 | **17,688** |
| 5σ counts | 19,844 vs 1,126 | **15 vs 18 (−3 ± 5.7)** |
| p_UL | 5.9e−3 | **4.6e−4** |
| **mean f̄** | 2.95e−3 | **2.14e−4** |

**13.8× better from 188× fewer stars.** The gain is not statistical.

### 5.7 A limit on the beaming-consistent class

The comparison in §1.3 points at what this dataset can uniquely contribute.
Suazo et al. (2022) constrain spheres that re-radiate isotropically at an
assumed temperature. Nobody has put a number on the class that does *not*.

Of 3,321,566 stars, 2,993,277 (90.1%) have usable W1 and W2, and 2,897,066
(87.2%) have a **measured** bare photosphere: |excess| < 3σ in both bands
against a (J−Ks)₀-calibrated photospheric colour. Absence of WISE data is never
counted as absence of excess — §5.4 showed that would fabricate candidates.

The veto removes very little of our tail:

| subsample | N | 5σ positive | rate | pos/neg |
|---|---|---|---|---|
| all stars | 3,321,566 | 19,844 | 0.00597 | 17.6 |
| bare photosphere | 2,897,070 | 13,611 | 0.00470 | 15.2 |
| IR excess present | 73,123 | 687 | 0.00940 | 21.5 |

The single-star limit improves by only **1.13×**. That is itself a result:
**our optical-deficit tail is not made of dusty objects.** Stars with genuine
mid-IR excess carry a rate only 2× the bare-photosphere rate and are 2.4% of
the sample. The tail is blends, as §5.5 established by a different route.

Applying the veto where our best estimator lives — clean wide pairs with
**both** components bare, 6,431 pairs / 12,862 stars:

| k | threshold (mag) | f detectable | n_pos | n_neg | asymmetry | p_UL | f̄ limit |
|---|---|---|---|---|---|---|---|
| 3 | 0.327 | 0.260 | 143 | 174 | −31 ± 17.8 | 2.5e−3 | 6.5e−4 |
| 5 | 0.545 | 0.395 | 8 | 12 | −4 ± 4.5 | 8.0e−4 | 3.2e−4 |
| **6** | 0.655 | **0.453** | 1 | 3 | −2 ± 2.0 | **4.9e−4** | **2.2e−4** |

**p < 4.9 × 10⁻⁴ at f ≥ 0.45 for stars with an optical deficit and a measured
absence of mid-infrared excess.** Every asymmetry is consistent with zero.

> **What "bare photosphere" actually means, quantified.** W1 and W2 sit at 3.4
> and 4.6 µm, so they only respond to *hot* shells.
> `scripts/38_temperature_coverage.py` computes the detectability of a shell
> re-radiating a solar luminosity at 500 pc: SNR is 1.8 × 10³ at 800 K, 454 at
> 400 K, **2.9 at 200 K**, and 7.6 × 10⁻⁶ at 100 K. Our veto is therefore only
> informative above roughly **300 K**.
>
> This makes the class broader than "beamed" — it is **beamed, OR cooler than
> ~300 K, OR non-thermal**. That is the correct reading of the limit and it is
> more useful than the narrow one, since a cold radiator is what a
> computation-limited civilisation is driven toward: the Landauer cost of a
> bit erasure is kT ln 2, so the value of energy is set by the sink
> temperature, and colder is strictly better.

This is the only limit in this paper not superseded by prior work, and it is
complementary rather than competitive: it covers the corner of parameter space
that an IR-excess estimator weights to zero.

### 5.8 Joint constraint: the first bound that does not assume a disposal mode

Every published limit constrains one way of getting rid of the waste heat.
Infrared searches constrain spheres that re-radiate isotropically in the WISE
bands. Deficit searches constrain interception regardless of disposal, but
weakly. Partition the stars intercepting a fraction ≥ f:

> **p_total = p_iso + p_dark**
>
> **p_iso** — re-radiates isotropically and warm; infrared searches see it.
> **p_dark** — beams, radiates too cold for WISE, or exports non-thermally;
> infrared searches do not see it, and only a deficit search can.

§5.7 measures **p_dark directly** by requiring a measured bare photosphere
alongside the deficit. Suazo et al. (2022) measure **p_iso**. They combine:

| | constraint | source |
|---|---|---|
| p_iso (isotropic, warm, 100 pc) | < 1.9 × 10⁻⁴ | Suazo et al. 2022 |
| **p_dark** (beamed / cold / non-thermal, 500 pc) | **< 4.9 × 10⁻⁴** | **this work, §5.7** |
| **p_total** | **< 6.8 × 10⁻⁴** | **combined** |

**Fewer than 1 in 1,473 nearby lower-main-sequence stars intercepts ≥45% of its
optical output by any means, warm or dark.** To our knowledge this is the first
bound on stellar energy harvesting that does not assume how the energy is
disposed of, and it is the natural product of running the two channels against
each other rather than in isolation.

#### When the beamed fraction becomes measurable

The quantity the whole "infrared searches are evadable" argument turns on is
β = p_dark/p_total. Both inputs are non-detections, so β is bounded only
trivially today. It becomes informative the moment the deficit limit drops
*below* the infrared limit, because then a detection in one channel and not the
other is decisive rather than ambiguous.

That requires **2.6×** improvement. Because the pair estimator is
background-subtracted, its limit scales as √N/N rather than saturating:

| pair sample | stars | projected p_UL | below IR limit? |
|---|---|---|---|
| current | 12,862 | 4.9 × 10⁻⁴ | no |
| 5× | 64,310 | 2.2 × 10⁻⁴ | no |
| **10×** | **128,620** | **1.5 × 10⁻⁴** | **yes** |
| 50× | 643,100 | 6.9 × 10⁻⁵ | yes |

El-Badry, Rix & Heintz (2021) catalogue ~1.3 × 10⁶ wide pairs. We use 6,431.
**A 10× increase — well within that catalogue — makes β a measured quantity
rather than a rhetorical one.** That is the single most valuable follow-up
this work identifies, and it requires no new observations.

### 5.8b Corrected limits, and a background prediction that tightens them

Two-sided counting is correct but conservative: every pair in either tail is
allowed to be signal. We can do better without circularity by *predicting* the
pair background from an independent measurement.

The pair tails are populated by per-star events — blends, activity, unresolved
tertiaries — which do **not** cancel in the difference. If a single star
exceeds |r| > T with probability q(T), a pair exceeds |Δr| > T when either
component does, so to first order

> **E[pair tail] = 2 · q(T) · N_pairs**

with q(T) measured on the 3.3M single-star sample at the same *absolute*
threshold in magnitudes. That prediction never touches the pair data, so
subtracting it is legitimate.

| threshold | f detectable | observed | predicted | excess |
|---|---|---|---|---|
| 0.448 mag | 0.338 | 146 | 169.3 | −23.3 |
| 0.560 mag | 0.403 | 33 | 65.8 | −32.8 |
| 0.673 mag | 0.462 | 10 | 32.0 | −22.0 |
| 0.785 mag | 0.515 | 8 | 18.3 | −10.3 |

**Observed counts sit consistently and substantially below prediction** — 10
against 32 at the 0.67 mag threshold. That is common-mode cancellation
removing roughly two thirds of the tail, measured rather than assumed, and it
is the quantitative justification for the pair estimator. Every excess is
negative: no signal anywhere.

**Corrected limits:**

| | limit | at |
|---|---|---|
| all clean pairs | **p < 4.8 × 10⁻⁴** | f ≥ 0.515 |
| **beamed class** (both components bare) | **p < 4.9 × 10⁻⁴** | f ≥ 0.505 |
| single stars, same threshold | p ≈ 2.9 × 10⁻³ | — |

The pair sample's background rate is **5.2× cleaner** than the single-star
sample at matched threshold. The joint constraint of §5.8 is unchanged at
p_total < 6.8 × 10⁻⁴ because the corrected p_dark is numerically almost
identical to the (invalid) asymmetry value — a coincidence, not a vindication.

### 5.8c Improved estimator: dropping the metallicity control *helps*

Two changes, one of which was counter-intuitive.

**(1) Requiring a metallicity was pure loss.** The fitted sample demanded a
GSP-Phot [M/H], costing 15% of stars and ~30% of pairs since both members must
survive. But within a co-natal pair the metallicity is *common to both
components and cancels exactly in the difference*. Refitting the fiducial with
only the near-infrared colour control, on all 3,884,167 stars:

| | with [M/H] | without [M/H] |
|---|---|---|
| σ (single star) | 0.09913 mag | 0.11338 mag |
| clean pairs | 8,844 | **12,214** (+38%) |
| pairs, both bare | 6,431 | **9,022** (+40%) |
| σ(Δr) | 0.11209 mag | **0.10957 mag** |
| common-mode variance | 34% | **53%** |

The single-star scatter gets *worse* (0.099 → 0.113) and the pair scatter gets
*better* (0.112 → 0.110). Removing a control for a quantity that is shared
within pairs moves that variance from the "per-star" bucket into the
"common-mode" bucket, where the difference cancels it for free. Common-mode
cancellation rises from 34% to 53%. **For a differential estimator, controlling
for a common-mode variable is actively counterproductive.**

Chance-alignment contamination stays at 0.03%.

**(2) The background prediction was targeted at the wrong population — and our
reasoning about which way it would go was wrong.** We expected pair members to
be *cleaner* than average, since both passed astrometric-consistency cuts.
Measured, they are **dirtier**: q(|r| > 0.67) is 3.40 × 10⁻³ for pair members
against 2.41 × 10⁻³ globally, a ratio of 1.41. Wide binarity correlates with
higher-order multiplicity, so a wide-pair member is *more* likely to have a
close companion as well.

The prediction therefore over-predicts by far more than we assumed —
63.5 expected against 7 observed at the 0.66 mag threshold — because
common-mode cancellation suppresses the tail so strongly. With observed counts
this small the conservative Poisson limit beats the background-subtracted one,
so **subtraction contributes nothing and the honest estimator is the plain
two-sided count.**

**Improved limits:**

| | limit | at | previous |
|---|---|---|---|
| all clean pairs | **p < 5.3 × 10⁻⁴** | f ≥ 0.507 | 4.8 × 10⁻⁴ |
| **beamed class** | **p < 4.3 × 10⁻⁴** | f ≥ 0.507 | 4.9 × 10⁻⁴ |
| **p_total (joint)** | **< 6.2 × 10⁻⁴** | | 6.8 × 10⁻⁴ |

**Fewer than 1 in 1,614 nearby lower-main-sequence stars intercepts ≥51% of its
optical output by any disposal mode.** A 12% improvement from 40% more pairs —
Poisson-limited statistics scale weakly when the counts are single digits,
which is exactly why the deeper sample matters more than any estimator tweak.

### 5.8d Dynamical mass: closing the spectral blind spot

Every limit above regresses M_G on M_Ks, which inherits two structural blind
spots: a **grey** absorber dims both bands and produces a residual of the
*wrong sign* (§3.3), and the anchor is 2MASS, whose beam is our dominant
contaminant (§5.5). A **dynamical** mass has neither problem. Gravity is
indifferent to what an absorber does to photons, and the mass comes from Gaia
astrometry rather than from a low-resolution infrared survey.

Regressing M_G directly on log M_dyn for 88,149 `gaiadr3.binary_masses`
systems (0.2–1.6 M☉, A₀ < 0.6, clean Ks):

| subset | N | σ | threshold | f detectable | n_pos | n_neg | p_UL |
|---|---|---|---|---|---|---|---|
| all systems | 88,149 | 0.277 mag | 1.386 mag | 0.721 | 14 | 41 | 2.5 × 10⁻⁴ |
| **faint secondary** (fluxratio < 0.1) | 13,482 | 0.288 mag | 1.440 mag | **0.734** | **0** | 8 | **2.2 × 10⁻⁴** |

**Zero positive outliers in 13,482 systems with clean dynamical masses.**

The negative excess (41 against 14, and 8 against 0 in the clean subset) is
the expected signature of unresolved secondary light: in a binary the observed
G contains both components while m1 is the primary alone, inflating the
luminosity. That biases toward *over*-luminosity — against the signal — so the
search is conservative by construction.

**p < 2.2 × 10⁻⁴ at f ≥ 0.73, for an absorber of ANY spectral slope.** This is
the first constraint here that covers the grey case and the α ≈ 0.19 blind
spot, which no version of the M_G–M_Ks estimator can reach at any sample size.
It costs a higher covering-fraction threshold (0.73 against 0.51) because
dynamical masses are noisy, but it constrains a class that was previously
entirely unconstrained by this work.

### 5.8e Kinematic clustering: a 4301σ excess that is young stars

§5.10 searched for a spreading front in 3D *position* and found dust.
But interstellar transfer cost is set by **Δv**, not distance: two stars 50 pc
apart co-moving at 2 km/s are adjacent in transfer cost, two stars 5 pc apart
differing by 60 km/s are not. Position-space clustering is arguably the wrong
metric, so we repeated the search in velocity space.

The first pass reported a **4301σ** excess of low-Δv anomaly pairs, 178× the
control rate. It is not a detection. Three controls were missing, and the
sequence of removing them is the result:

| control added | excess ratio at Δv < 2 km/s |
|---|---|
| none | 178× |
| exclude gravitationally bound pairs (Δv = 0 by construction) | 160× |
| match controls on \|v\| | 43× |
| match on pair separation, and restrict to \|v\| > 40 km/s | **~1.4×** |

The decisive diagnostic: **anomalies are kinematically cold**, median |v| =
22.8 km/s against 35.1 for the field, with 43% below 20 km/s against 23%.
Kinematically cold means young, young stars sit in comoving associations, and
associations cluster in velocity for entirely mundane reasons. Restricting to
the dynamically heated disc (|v| > 40 km/s), where coherent young associations
are absent, collapses ratios of 8–29× to ~1.4×.

No kinematic front. But the test has real teeth — it detected a genuine
population at enormous significance — and the failure mode is instructive: in
velocity space the young-star contaminant is *more* concentrated than in
position space, so this metric needs stronger controls, not weaker.

### 5.9 Deficit spectral type

Classifying each outlier by whether its colour excess matches the reddening
vector: positive outliers have median δC_obs/δC_dust = **1.09**. They redden.
There is no grey population and no exotic-spectrum population — nothing near
the α ≈ 0.19 blind spot where a deliberate evader would sit.

### 5.10 Search for a spreading front

A spreading civilisation produces *spatially correlated* anomalies, and
correlation statistics detect clustered weak signals far below the per-object
threshold. We regress out all non-positional systematics with a model never
shown (l, b, d) — so it cannot absorb a front — and scan 31,184 cells of 25 pc.

**Real structure is found**: cell-mean scatter 1.554σ against a shuffled null of
1.001, 137 cells above 5σ where 40 are expected, and a two-point correlation of
the cleaned residual field significant at **19.5σ** on 25 pc scales.

It is dust. Anomalous cells sit at median |b| = 15.5° against 30.7° overall,
with median A₀ = 0.176 against 0.090. Restricting to |b| > 30° and A₀ < 0.15
the structure weakens but persists and *still* points at the residual dust.
The morphology is wrong for a front: 157 separate clumps, largest 22 cells
against a null of 5.4, with no edge.

**No front.** But the test detected a real 19σ correlated signal that every
single-star statistic here was blind to, and unlike the tail statistic it is
extinction-limited rather than background-limited — its sensitivity still
improves with N.

---

## 6. Discussion

### 6.1 Answer to the motivating estimate

The σ/√N argument fails for three independent reasons:

1. It does not apply to the mean, which is not measurable at all (§3.2).
2. Where it could apply it stops at 2–5 × 10⁻² mag, **957–1376× above naive**,
   on coherent systematics that do not average down (§4).
3. The tail-based limit is background-limited by companions, YSOs,
   emission-line stars and spotted rotators, and that background scales with N
   (§5.3).

**The optical-deficit channel saturates at f̄ ~ 2 × 10⁻⁴ and cannot reach
10⁻⁵.**

### 6.2 Comparison with prior work

Against [Zackrisson et al. (2015)](https://doi.org/10.1088/0004-637X/810/1/23):
their <3% is a limit on the fraction of *galaxies*, ours on the fraction of
*stars* within 500 pc. Not directly comparable, but both are limited by the
intrinsic width of the scaling relation regressed against — Tully–Fisher for
them, the main sequence for us — rather than by photon noise. Our relation is
several times tighter per object; our limiting factor is contamination.

Against [Suazo et al. (2022)](https://doi.org/10.1093/mnras/stac1789),
Hephaistos I — the direct numerical comparison:

| | covering fraction | limit on p |
|---|---|---|
| Suazo+22, 100 pc, 300 K | γ ≥ 0.5 | **1.9 × 10⁻⁴** |
| **This work, clean pairs, 500 pc** | f ≥ 0.46 | 4.6 × 10⁻⁴ |
| Suazo+22, 100 pc, 300 K | γ ≥ 0.9 | 1.8 × 10⁻⁵ |

**Their limit is ~2.4× stronger at comparable covering fraction, from ~12×
fewer stars.** The reason is instructive: their estimator uses the infrared
re-emission as well as the obscuration, and infrared excess is measured against
a near-zero photospheric background whereas a deficit is measured against a
0.089 mag main-sequence width. Adding a second, near-background-free channel is
worth more than adding an order of magnitude in N — which is the same lesson
§4.2 teaches from the other direction.

The compensating argument for the deficit channel is coverage, not depth. A
limit derived from a model that includes γ·BB(T_DS) constrains spheres that
re-radiate into 4π at the assumed temperature. A beamed radiator, a sphere cold
enough to fall outside WISE, or one that exports its energy non-thermally is
weighted differently by that estimator and not at all by ours. We cannot
quantify how much of that class Hephaistos retains without re-running their
pipeline with the re-emission term suppressed, and we do not claim a number
for it.

Against [Zackrisson et al. (2018)](https://doi.org/10.3847/1538-4357/aac386):
comparable per-object sensitivity (their f_cov > 0.75, our f ≳ 0.46) from a
different observable, on 400× fewer stars. Their single candidate was killed by
an unseen binary companion; we find binarity to be the dominant contaminant
from an independent direction. The agreement on the *failure mode* is stronger
evidence than either limit.

Against [Zackrisson et al. (2015)](https://doi.org/10.1088/0004-637X/810/1/23):
their <3% is a limit on the fraction of *galaxies*, ours on the fraction of
*stars*. Not directly comparable, but both are limited by the intrinsic width
of the scaling relation regressed against — Tully–Fisher for them, the main
sequence for us — rather than by photon noise.

### 6.3 What this does and does not constrain

The honest accounting. Within 500 pc, among lower-main-sequence stars, fewer
than 1 in 2,200 intercepts ≥46% of its optical output via a spectrally
selective absorber. Three structural blind spots remain: a **uniform** deficit
(unbounded by any self-calibrated method), an absorber at **α ≈ 0.19**, and
**f ≲ 0.1** at any abundance.

That last is the one that matters for interpretation. A civilisation taking 1%
of its star — still ~10²⁴ W, vastly more than any plausible need — is entirely
invisible to this method. **We constrain profligate near-total conversion, not
sipping.**

In terms of hypothesis space: the searchable (f, p) plane spans roughly 21
decade², and this work excludes roughly 0.8 of them, i.e. **~4–8%**. A null
result at this sensitivity is genuine but weak evidence, because the hypothesis
space extends decades below any version of this method's reach.

One thing *is* meaningfully constrained: were we inside an active conversion
front, we would see it locally regardless of volume fraction. We do not.

### 6.4 The programme this points to

The comparison with Suazo et al. (2022) contains its own instruction. They beat
us by 2.4× using ~12× fewer stars because they combine two channels, one of
which is measured against a near-zero background. The route to a competitive
deficit-channel limit is therefore *not* more stars — §4.5 shows that fails
outright — but a better estimator on the stars we already have.

In the order we would run them:

**(i) Scale the wide-pair differential.**
[El-Badry, Rix & Heintz (2021)](https://doi.org/10.1093/mnras/stab323) catalogue
~1.3 × 10⁶ Gaia pairs to 1 kpc against our 8,844. The pair estimator is the
only one here whose background is *measurably* symmetric and therefore
subtractable. Reaching even 10⁵ clean pairs would put p ~ 10⁻⁵–10⁻⁴ at
f ≈ 0.45 — **competitive with the isotropic-case limits, in a class those
limits do not cover.** Highest-value next step, and it needs no new
observations.

**(ii) Extend the stellar mass range.** We restricted to 3 < M_Ks < 8 (roughly
G2V–M4V). More luminous stars are both more attractive engineering targets and
better represented in the comparison samples. The turnoff and subgiant
contamination that motivated the cut is tractable with Gaia parallaxes and
log g.

**(iii) Then the joint constraint.** With a deficit limit below the IR-excess
limit at matched covering fraction, the *difference* between the two bounds the
fraction of harvesters whose waste heat is beamed, cold, or non-thermal —
Ω/4π in the language of §1.1. Nobody has measured that, it is the quantity the
entire beaming argument turns on, and it requires exactly one channel to be
pushed below the other. Step (i) is what makes it possible.

**(iv) A NIR anchor at Gaia-like resolution** (VISTA VHS/VVV, UKIDSS) attacks
the dominant contaminant at its cause (§5.5) rather than filtering on a proxy
for it, which §6.5 shows saturates at 1.9×.

### 6.5 Where the remaining leverage is

In measured order of value:

1. **A NIR anchor at Gaia-like resolution** (VISTA VHS/VVV, UKIDSS). §5.5 shows
   this is the physical cause of the dominant contaminant, and it is the only
   lever that plausibly moves the *tail fraction* rather than filtering on a
   proxy for it. Hard cuts on the proxies saturate at 1.9× because the residual
   distribution is **scale-free** — cleaning shrinks σ and the tail together.
2. **More clean wide pairs.** [El-Badry, Rix & Heintz (2021)](https://doi.org/10.1093/mnras/stab323)
   have ~1.3 × 10⁶ pairs to 1 kpc against our 8,844.
3. **Gaia XP spectra.** Stop looking for missing flux and look for the wrong
   spectral *shape*. Every background here changes the SED smoothly, within the
   space stellar atmospheres already span; an engineered absorber plausibly has
   a band edge, which no atmosphere produces. This also *measures* α instead of
   being blind at 0.19.
4. **Eclipsing binaries and asteroseismology** — the only route to bounding the
   uniform case, since R and T_eff give an absolute luminosity prediction with
   no main-sequence-width systematic.

### 6.6 Errors we made

Recorded because they are instructive and because a reader should be able to
judge the pipeline's reliability.

- The **parallax zero point was applied 1000× too small** (`gaiadr3-zeropoint`
  returns mas, not µas). Silent — the diagnostic printed "−0.0 µas", which
  reads as negligible rather than broken. Caught only because it contradicted
  the literature value quoted two lines above in the same docstring.
- **Gaia band wavelengths were DR2-era values**, mis-cited as DR3.
- We assumed hard cuts on contamination proxies would buy ~10×, from a real
  20× dependence of the tail rate on C*. Measured: **1.9×**, because the tail is
  scale-free. Trusting the 20× would have been wrong by an order of magnitude.
- We expected wide-binary common-mode cancellation to be large. It is **34%**.
- The blinding scheme did not blind the primary statistic (§5.2).
- **The wide-pair limit was built on an estimator with no sensitivity.** We
  validated the asymmetry statistic against the *background* (measuring it to be
  symmetric) and never against an injected *signal*, which is also symmetric.
  Injection response: 0.02σ. Caught only by asking, late, whether the estimator
  responds at all. The injection machinery to catch it existed from §5.1 and we
  had not pointed it here. Corrected in §5.6 and §5.8b.
- **We initially framed this as the first stellar optical-deficit search and as
  setting the best limit. Both were wrong** — Zackrisson et al. (2018) precedes
  it methodologically and Suazo et al. (2022) already had a stronger limit at
  comparable covering fraction. The literature review that found this was done
  after the analysis rather than before it, which is the wrong order and is
  exactly how a novelty claim gets overstated. §1.3 and §6.2 carry the
  correction.

---

## 7. Conclusions

1. The optical-deficit channel is **background-limited, not statistics-limited**.
   The naive σ/√N sensitivity overstates the truth by **957–1376×**. This is
   the main result and we believe it is new.
2. Our constraint is **f̄ < 2.1 × 10⁻⁴**, **p < 4.6 × 10⁻⁴** at f ≥ 0.46, from
   clean wide pairs, with **zero candidates**. This is *not* the strongest
   published limit: Suazo et al. (2022) reach 1.9 × 10⁻⁴ at γ ≥ 0.5 using the
   infrared re-emission as well.
3. Requiring a **measured bare photosphere** gives **p < 4.9 × 10⁻⁴ at
   f ≥ 0.45** for the beamed / cold / non-thermal class — the one limit here
   not superseded, because infrared estimators weight that class to zero.
4. Combining the two channels gives **p_total < 6.8 × 10⁻⁴**, the first
   disposal-agnostic bound: fewer than 1 in 1,473 nearby lower-main-sequence
   stars intercepts ≥45% of its output by any means. A 10× larger wide-pair
   sample makes the beamed fraction β measurable.
5. **A uniform deficit is unmeasurable** by any self-calibrated version of this
   method.
6. **Unresolved companions mimic the signal**, because dM_G/dM_Ks = 1.26 > 1.
7. **The blind spot is α ≈ 0.19, not a grey absorber**; grey gives the opposite
   sign.
8. The dominant contaminant is **Gaia/2MASS aperture mismatch**, localised to
   θ ≲ 6.5″, and it limits any 2MASS-anchored version of this search.

---

## Reproducing

```bash
bash setup_env.sh
wsl -d kali-linux bash run.sh scripts/02_pull_sample.py --workers 3 --distance-max-pc 500
wsl -d kali-linux bash run.sh scripts/05_build_sample.py --pattern 'sample_d500_p*' --tag primary
wsl -d kali-linux bash run_primary.sh
wsl -d kali-linux bash run_after_unblind.sh
```

40 tests: `python -m pytest tests/ -q`. See `RESULTS.md`, `LIMITATIONS.md`,
`DECISIONS.md`.

## References

- Annis J., 1999, JBIS, 52, 33
- Boyajian T. S. et al., 2016, MNRAS, 457, 3988
- Carrigan R. A., 2009, ApJ, 698, 2075
- Dyson F. J., 1960, Science, 131, 1667
- Edenhofer G. et al., 2024, A&A, 685, A82
- El-Badry K., Rix H.-W., Heintz T. M., 2021, MNRAS, 506, 2269
- Fitzpatrick E. L. et al., 2019, ApJ, 886, 108
- Green G. M. et al., 2019, ApJ, 887, 93
- Griffith R. L. et al., 2015, ApJS, 217, 25
- Lindegren L. et al., 2021, A&A, 649, A4
- Riello M. et al., 2021, A&A, 649, A3
- Skrutskie M. F. et al., 2006, AJ, 131, 1163
- Suazo M., Zackrisson E., Wright J. T., Korn A. J., Huston M., 2022, MNRAS, 512, 2988 — *Project Hephaistos I: upper limits on partial Dyson spheres*
- Suazo M. et al., 2024, MNRAS, 531, 695 — *Project Hephaistos II*
- Villarroel B. et al., 2020, AJ, 159, 8 — *VASCO I*
- Zackrisson E., Korn A. J., Wehrhahn A., Reiter J., 2018, ApJ, 862, 21 — *SETI with Gaia: nearly complete Dyson spheres*
- Wang S., Chen X., 2019, ApJ, 877, 116
- Wright J. T. et al., 2014a, ApJ, 792, 26
- Wright J. T. et al., 2014b, ApJ, 792, 27
- Wright J. T., 2023, ApJ, 956, 34
- Zackrisson E. et al., 2015, ApJ, 810, 23
