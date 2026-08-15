# A search for the optical deficit signature of stellar energy harvesting in 3.3 million Gaia DR3 main-sequence stars

**Status:** working paper, self-published. Not peer reviewed. All code, data
provenance and null tests are in this repository; every number below is
regenerable by a numbered script.

---

## Abstract

Essentially every survey for Dyson-type stellar energy harvesting has searched
for mid-infrared **excess** from reprocessed waste heat. That channel is
evadable: a civilisation that beams its waste heat into a solid angle Î© is
detectable by only Î©/4Ï€ of observers. The optical **deficit** â€” starlight that
is intercepted and never departs â€” is visible from every direction and cannot
be evaded by beaming. Its cost is sensitivity, and that cost has not previously
been measured.

We measure it. From 4,809,840 Gaia DR3 Ã— 2MASS sources within 500 pc we build a
sample of 3,884,167 lower-main-sequence stars (3,321,566 with GSP-Phot
metallicities) and search for a deficit in absolute *G* at fixed absolute
*K*<sub>s</sub>, after 3D-extinction and metallicity control.

The intrinsic main-sequence scatter is **0.0886 mag**, giving a naive
sensitivity of Ïƒ/âˆšN = 5.4 Ã— 10â»âµ mag. **That number is wrong by three orders of
magnitude.** Structured null tests â€” sky patches, extinction bins, magnitude
bins â€” plateau at 1.6â€“2.3 Ã— 10â»Â² mag while random subsamples track Ïƒ/âˆšN to
within 1% out to N = 3 Ã— 10â¶. The measured systematic floor is **957Ã— the naive
expectation**, and it is dominated by an *astrophysical* term (0.068 mag) that
survives where there is no dust.

We report **zero candidates**. Of the 300 strongest optical-deficit outliers,
none has measured-and-normal mid-infrared photometry. From 8,844 co-natal wide
pairs, where the background is measurably symmetric and can therefore be
subtracted, we obtain **fÌ„ < 2.1 Ã— 10â»â´** and **p < 4.6 Ã— 10â»â´** of stars
intercepting â‰¥46% of their optical output.

This is *not* the strongest published limit â€” [Suazo et al.
(2022)](https://doi.org/10.1093/mnras/stac1789) reach 1.9 Ã— 10â»â´ at Î³ â‰¥ 0.5 by
using the infrared re-emission as well. What we can claim uniquely is a limit
on the class their estimator is not designed to constrain: requiring both pair
components to have a **measured bare photosphere** in W1 and W2, we obtain
**p < 4.9 Ã— 10â»â´ at f â‰¥ 0.45** for stars with a large optical deficit *and no
warm re-emission* â€” the beaming-consistent, cold, or non-thermally-exporting
case that motivates the deficit channel in the first place.

Three results correct assumptions we began with, and we consider them the more
durable contribution:

1. **A uniform deficit is unmeasurable by any self-calibrated version of this
   method** â€” sensitivity exactly zero, not Ïƒ/âˆšN, demonstrated by injection
   across four decades.
2. **Unresolved companions mimic the signal rather than suppressing it**,
   because d*M<sub>G</sub>*/d*M*<sub>Ks</sub> = 1.26 > 1. Confirmed on data.
3. **The spectral blind spot is not a grey absorber.** It lies at Î± â‰ˆ 0.19;
   a grey absorber produces a residual of the *opposite* sign.

The dominant contaminant is identified mechanically: **aperture mismatch**
between Gaia (sub-arcsecond) and 2MASS (~4â€³). We localise it to Î¸ â‰² 6.5â€³ in
wide pairs, where the contaminated-to-clean tail ratio is 50:1 and falls to
unity beyond 10â€³. This limits *any* optical-deficit search anchored on 2MASS.

---

## 1. Introduction

### 1.1 The infrared-excess paradigm and its blind spot

[Dyson (1960)](https://doi.org/10.1126/science.131.3414.1667) pointed out that a
civilisation intercepting a substantial fraction of its star's output must
re-radiate that energy at low temperature, producing an infrared excess. Every
major search since has followed that logic:
[Carrigan (2009)](https://doi.org/10.1088/0004-637X/698/2/2075) searched ~250,000
IRAS sources for 100â€“600 K blackbodies and found that the best candidates were
reddened dusty objects â€” heavily extinguished stars, protostars, Mira
variables, AGB stars and planetary nebulae. The Äœ programme
([Wright et al. 2014a](https://arxiv.org/abs/1408.1133),
[2014b](https://doi.org/10.1088/0004-637X/792/1/27);
[Griffith et al. 2015](https://doi.org/10.1088/0067-0049/217/2/25)) extended
this to WISE with ~1000Ã— better sensitivity and 5Ã— better resolution, over both
Galactic sources and ~10âµ galaxies. Project Hephaistos
([Suazo et al. 2022](https://doi.org/10.1093/mnras/stac1789);
[2024](https://doi.org/10.1093/mnras/stae1186)) applied the modern
Gaia DR3 Ã— 2MASS Ã— AllWISE combination to ~5 Ã— 10â¶ sources, reporting seven
candidates after a filtering cascade.

The blind spot is geometric. Infrared excess is a statement about where the
waste heat *goes*. [Wright (2023)](https://doi.org/10.3847/1538-4357/acf44f)
showed thermodynamically that optimal use of mass favours small, hot spheres,
and that a sphere accepting radiation from a limited solid angle behaves
differently from the isotropic idealisation. A radiator that beams into solid
angle Î© is seen by a fraction Î©/4Ï€ of observers. Nothing in physics forbids
this, and a civilisation with any preference for not being conspicuous â€” or
merely one that finds it useful to aim its waste heat somewhere â€” is
under-counted by every infrared survey by that factor.

### 1.2 The optical deficit channel

The complementary observable is the light that never leaves.
[Annis (1999)](https://ui.adsabs.harvard.edu/abs/1999JBIS...52...33A)
originated the method at galaxy scale, looking for outliers in the
Tullyâ€“Fisher relation and the fundamental plane, adopting a criterion of â‰¥1.5
mag of dimming (a factor of 4) on 57 disks and 106 ellipticals, and finding
none. [Zackrisson et al. (2015)](https://doi.org/10.1088/0004-637X/810/1/23)
made this quantitative on 1,359 disk galaxies, limiting to **<3%** (0.3%
tentative) the fraction of local disk galaxies subject to galaxy-scale Dysonian
astroengineering.

A deficit cannot be beamed away. It is visible to every observer in every
direction, and it is agnostic about what happens to the energy afterwards â€” it
does not care whether the waste heat is radiated, beamed, stored, or used for
computation at the Landauer limit. That is the entire appeal.

The obvious extension is from galaxies to individual stars, where Gaia provides
10â¹ objects with parallaxes. The back-of-envelope is seductive: with a
main-sequence relation of width Ïƒ â‰ˆ 0.1 mag, a single star gives a 3Ïƒ limit of
f â‰ˆ 26% (useless), but N ~ 10â¸ stars would give Ïƒ/âˆšN ~ 10â»âµ mag, i.e.
fÌ„ ~ 10â»âµ (decisive).

**The purpose of this paper is to show that this scaling fails, to measure
where and why, and to report what the channel actually delivers.**

### 1.3 This is not the first stellar optical-deficit search, and not the best limit

We state this plainly because we initially framed the work as more novel than
it is, and a reader deserves the correction up front.

**[Zackrisson et al. (2018)](https://doi.org/10.3847/1538-4357/aac386)** is the
direct methodological ancestor. They search for optically underluminous stars
by comparing spectrophotometric to trigonometric distance, since a dimmed star
has its spectrophotometric distance overestimated while its parallax does not
move: f_cov = 1 âˆ’ (D_trig/D_spec)Â². On ~230,000 Gaia DR1 Ã— RAVE DR5 stars
(8,441 after cuts) they find 75 stars at f_cov > 0.9, reducing to six
unexplained outliers after vetting.

Their method is sensitive **only to nearly complete spheres** â€” it takes
f_cov > 0.75 to produce the factor-of-two distance discrepancy the method
needs, and a 10% covering fraction gives a ~5% distance disagreement that is
undetectable. That is the same sensitivity regime we reach, from a different
observable.

Their prime candidate, TYC 6111-1162-1 (f_cov â‰ˆ 0.77, and notably *without*
mid-infrared excess), resolved to a **single-lined spectroscopic binary** whose
unseen ~1 Mâ˜‰ white-dwarf companion perturbed the Gaia parallax. We regard this
as important independent corroboration: their one candidate died of binarity,
and we find binarity to be the dominant contaminant by a completely different
route (Â§5.5). Two methods, two different observables, the same false-positive
population.

**[Suazo et al. (2022)](https://doi.org/10.1093/mnras/stac1789)** (Hephaistos I)
already set upper limits on partial Dyson spheres using Gaia DR2 Ã— AllWISE with
the AGENT formalism, modelling the system as
L_Î½ = (1âˆ’Î³)L_Î½,â‹† + Î³Â·BB_Î½(T_DS) â€” that is, **using the optical obscuration and
the infrared re-emission together**. Within 100 pc at T_DS = 300 K they obtain

| covering fraction | upper limit on fraction of stars |
|---|---|
| Î³ â‰¥ 0.1 | 6.6 Ã— 10â»Â³ |
| Î³ â‰¥ 0.5 | **1.9 Ã— 10â»â´** |
| Î³ â‰¥ 0.9 | 1.8 Ã— 10â»âµ |

**Their limit at Î³ â‰¥ 0.5 is roughly 2.4Ã— stronger than our best
(4.6 Ã— 10â»â´ at f â‰¥ 0.46).** We do not set a record and we do not claim one.

### 1.4 What is actually new here

Given the above, the contribution of this work is not the limit. It is:

1. **A measurement of the systematic floor of the deficit channel.** No prior
   work quantifies where Ïƒ/âˆšN stops applying. We find it stops **957â€“1376Ã—**
   above the naive value, and we identify what sets it.
2. **Identification and localisation of the dominant contaminant** as
   Gaia/2MASS aperture mismatch, pinned to Î¸ â‰² 6.5â€³ (Â§5.5). This is structural
   and limits any 2MASS-anchored search, including future ones.
3. **Three methodological corrections** that we believe are not stated
   correctly elsewhere: a uniform deficit is unmeasurable by a self-calibrated
   relation (Â§3.2); unresolved companions mimic rather than suppress the signal
   because dM_G/dM_Ks > 1 (Â§5.5); and the spectral blind spot lies at Î± â‰ˆ 0.19
   rather than at a grey absorber (Â§3.3). The last bears directly on Hephaistos'
   choice to model obscuration as grey.
4. **A differential wide-pair technique** with a *measured* symmetric
   background, which is the only estimator here that permits background
   subtraction (Â§5.6).

A note on method complementarity: our estimator regresses against an empirical
main-sequence relation, which is why a uniform or grey deficit is invisible to
it. Zackrisson et al. and Suazo et al. use absolute SED fitting against a
parallax, which does **not** share that blind spot â€” a grey absorber makes a
star underluminous in every band and is directly detectable that way. Our blind
spots and theirs are different, which is an argument for running both rather
than for preferring either.

### 1.3 Related approaches

Other lines of attack share our motivation of looking for absence rather than
emission. The VASCO project
([Villarroel et al. 2016](https://doi.org/10.3847/1538-3881/ab570f), 2022)
compares century-old photographic plates against modern surveys for sources
that have *vanished* â€” a limiting case of a deficit, and one with the enormous
advantage of a pre-Sputnik baseline free of satellite contamination.
[Boyajian et al. (2016)](https://doi.org/10.1093/mnras/stv1233) established
that deep, aperiodic, otherwise-unexplained dimming of a single main-sequence
star is detectable in practice, which is the time-domain analogue of what we
search for photometrically.

---

## 2. Data

### 2.1 Query and sample construction

We query the ESA Gaia archive TAP service directly, joining
`gaiadr3.gaia_source` to 2MASS via `tmass_psc_xsc_best_neighbour` â†’
`tmass_psc_xsc_join` â†’ `gaiadr1.tmass_original_valid`, to AllWISE via
`allwise_best_neighbour` â†’ `gaiadr1.allwise_original_valid`, and to
`gaiadr3.astrophysical_parameters`. The sky is partitioned by `source_id`
range, which encodes a level-12 HEALPix index in its high bits, so each chunk
is simultaneously an indexed primary-key scan and a contiguous sky patch.

192 partitions, **4,809,840 rows, zero failures**, 252 minutes at three
concurrent workers.

### 2.2 Cuts

| stage | stars | removed |
|---|---|---|
| raw joined rows (server-side ADQL cuts) | 4,879,956 | â€” |
| 2MASS Ks `ph_qual` = 'A' | 4,865,655 | 0.29% |
| not photometrically variable | 4,497,637 | 7.83% |
| `non_single_star` = 0 | 4,497,475 | 7.84% |
| not QSO/galaxy candidate | 4,497,459 | 7.84% |
| DSC star probability > 0.5 | 4,323,603 | 11.40% |
| \|C*\| < 3Ïƒ (Riello et al. 2021) | 4,289,991 | 12.09% |
| dust-map coverage | 4,275,099 | 0.35% |
| 10 < d < 500 pc | 4,275,099 | 0.35% |
| A_G < 0.5 | 3,925,746 | 8.49% |
| 3.0 < M_Ks < 8.0 | 3,885,617 | 9.43% |
| 0.7 < (BPâˆ’RP)â‚€ < 3.6 | **3,884,167** | 9.46% |
| with GSP-Phot [M/H] | **3,321,566** | |

Server-side cuts include `parallax_over_error > 20`, `ruwe < 1.4`,
`ipd_frac_multi_peak â‰¤ 2`, `visibility_periods_used â‰¥ 10`, and 1:1 cross-match
uniqueness (`number_of_mates = 0`, `number_of_neighbours = 1`).

### 2.3 Astrometry and extinction

Parallaxes are corrected using the
[Lindegren et al. (2021)](https://doi.org/10.1051/0004-6361/202039653)
zero point; the median correction is **âˆ’39.0 Âµas** (16â€“84%: âˆ’43.8 to âˆ’33.3),
i.e. a 1.3% distance-scale effect and 0.031 mag of distance modulus at our
median parallax. This is an order of magnitude above the target sensitivity and
cannot be neglected.

Extinction uses the [Edenhofer et al. (2024)](https://doi.org/10.1051/0004-6361/202347628)
3D map (full sky within 1.25 kpc, Gaia XP-based), converted to Aâ‚€ via
A_V = 2.8 E, with the official Gaia colour-dependent coefficients from the
[Fitzpatrick et al. (2019)](https://doi.org/10.3847/1538-4357/ab4c3a)
extinction law. The colour dependence is not optional: k_G = A_G/Aâ‚€ runs from
1.00 at (BPâˆ’RP)â‚€ = 0 to 0.65 at 3, and a constant ratio would imprint a
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

- M_G = G âˆ’ 5 logâ‚â‚€ d + 5 âˆ’ A_G
- M_Ks = Ks âˆ’ 5 logâ‚â‚€ d + 5 âˆ’ A_Ks

If a fraction *f* of the optical output is intercepted,
**Î”M_G = âˆ’2.5 logâ‚â‚€(1 âˆ’ f)** (positive = fainter).

We fit a fiducial relation M_G = ð”£(M_Ks, [M/H], (Jâˆ’Ks)â‚€) as a cubic B-spline in
M_Ks at quantile knots, tensored with covariate polynomials, by Huber IRLS.
Complexity is chosen by 5-fold cross-validation; the measured loss surface is
monotone in both directions with its minimum at 6 interior knots and
metallicity degree 1 (`results/cv_primary.csv`).

Controlling on the *near-infrared* colour (Jâˆ’Ks)â‚€ is safe and reduces the
scatter from 0.175 to 0.099 mag: an optically selective absorber leaves a NIR
colour alone, so the control cannot launder the signal. We run a second variant
that also controls on the *optical* (BPâˆ’RP)â‚€, which reduces the scatter to
0.028 mag but restricts sensitivity to absorbers grey across the optical.

### 3.2 What is measurable, and what is not

**A uniform deficit is not measurable.** The fiducial is fit to the same stars
being searched. A constant offset lies exactly in the span of the basis and is
absorbed. This is not a bound, it is an identity, and we verify it: injections
from f = 10â»âµ to 10â»Â¹ return a mean residual identical to eight decimal places
(0.00629139 in every case) with an identical 5Ïƒ tail count.

What survives is the **shape** of the residual distribution at fixed M_Ks.
Harvesting adds a one-sided fainter tail.

A *sparse* deficit does move the mean, contrary to the naive expectation: the
robust fit down-weights injected outliers rather than following them, so ~77%
of the pÂ·Î” shift survives. We nonetheless use the tail, because the mean is
where the coherent systematics live: at the measured floor a mean-based limit
at f = 0.5 gives p < 0.09, against p < 5 Ã— 10â»Â³ from the tail.

### 3.3 Spectral leverage: where the method is actually blind

The residual for an attenuated star is **r = Î´m_G âˆ’ s Â· Î´m_Ks** with
s = dM_G/dM_Ks. The critical point is that **s â‰  1**: we measure
s = 1.2566 (16â€“84%: 1.0597â€“1.4567), so the slope exceeds unity essentially
everywhere in the sample.

For a power-law absorber Ï„ âˆ Î»^âˆ’Î±:

| absorber | Î´m_G/Î´m_Ks | leverage |
|---|---|---|
| grey (Î± = 0) | 1.00 | **âˆ’0.257** |
| **Î± = 0.193** | 1.257 | **0.000 â€” blind spot** |
| Î± = 1 | 3.33 | +0.623 |
| Î± = 2 (dust-like) | 11.66 | +0.892 |

**A grey absorber is not invisible; it produces a residual of the opposite
sign.** The true blind spot is a moderately selective absorber at Î± â‰ˆ 0.19.
Interstellar dust sits at Î± â‰ˆ 2, mimicking the signal at 89% efficiency, which
is why extinction carries so much of the systematic budget.

Î»_eff is SED-dependent and this matters: Gaia G spans 400â€“950 nm, so for a
4500 K photosphere the flux-weighted effective wavelength is 0.679 Âµm, not the
Vega-referenced 0.582 Âµm â€” a 17% difference propagating to ~10% in Î±_blind.

---

## 4. The systematic floor

**This section was completed before any signal was examined.**

### 4.1 Null splits

Every two-sided split that must return zero fails, the mildest at 7Ïƒ:

| split | difference (mag) | significance |
|---|---|---|
| **colour, blue vs red** | **âˆ’0.05205** | **âˆ’353Ïƒ** |
| extinction quartile | âˆ’0.04456 | âˆ’212Ïƒ |
| apparent G, bright vs faint | âˆ’0.02818 | âˆ’189Ïƒ |
| crowding | âˆ’0.02612 | âˆ’174Ïƒ |
| galactic latitude | +0.02730 | +173Ïƒ |
| hemisphere N vs S | âˆ’0.00550 | âˆ’37Ïƒ |
| distance, near vs far | âˆ’0.00106 | âˆ’7Ïƒ |

### 4.2 Random subsamples are a control, not a measurement

| N | RMS of subsample means | Ïƒ/âˆšN |
|---|---|---|
| 10,000 | 0.001361 | 0.001367 |
| 1,000,000 | 0.000132 | 0.000137 |
| **3,000,000** | **0.000080** | **0.000079** |

Agreement to 1% at N = 3 Ã— 10â¶. Random subsampling cannot manufacture a
systematic, so it reproduces the naive scaling however wrong that scaling is as
a sensitivity. **Anyone quoting Ïƒ/âˆšN would find this reassuring and would be
wrong by three orders of magnitude.**

Structured groups on identical residuals plateau instead: apparent-G bins at
0.0229, extinction bins 0.0183, sky patches 0.0169, crowding 0.0153, distance
0.0034 mag.

### 4.3 The floor is astrophysical

Repeating the worst split *inside the lowest-extinction quartile*:

| split | full sample | lowest-Aâ‚€ quartile |
|---|---|---|
| colour | âˆ’0.05205 (353Ïƒ) | **âˆ’0.06779 (259Ïƒ) â€” grows** |
| distance | âˆ’0.00106 (7Ïƒ) | +0.01031 (38Ïƒ) |

The colour systematic is **larger** where there is no dust. It is real
main-sequence structure that M_Ks, [M/H] and (Jâˆ’Ks) do not capture â€” age,
Î±-enhancement, rotation, activity. No better dust map moves it.

### 4.4 Budget and the number

| term | mag |
|---|---|
| **astrophysical: MS structure vs colour** | **0.06779** |
| extinction: low vs high Aâ‚€ | 0.04456 |
| photometric: bright vs faint G | 0.02818 |
| crowding | 0.02612 |
| spatial coherence plateau | 0.01694 |
| extinction fractional error Ã— median A_G | 0.01314 |
| map vs Gaia per-star A_G | 0.00702 |
| parallax zero-point residual | 0.00195 |
| metallicity calibration | 0.00099 |
| band law (Fitz19 vs Wang & Chen 19) | 0.00005 |
| **quadrature sum** | **0.09258** |
| **largest single term** | **0.06779** |
| naive Ïƒ/âˆšN | 0.0000544 |
| **ratio** | **1246Ã—** |

Regressing residual on Aâ‚€ gives 0.1253 Â± 0.0005 mag per unit Aâ‚€ (247Ïƒ),
implying the extinction correction is â‰¤19% short.

### 4.5 The distance trade

Same partitions, same complexity, varying only the distance limit:

| d (pc) | N | median Aâ‚€ | Ïƒ/âˆšN | spatial plateau | floor/naive |
|---|---|---|---|---|---|
| 200 | 67,061 | 0.033 | 3.28eâˆ’4 | 0.0154 | 149Ã— |
| 500 | 480,700 | 0.110 | 1.42eâˆ’4 | 0.0227 | 345Ã— |
| 1250 | 1,383,542 | 0.166 | 7.51eâˆ’5 | 0.0235 | 536Ã— |

21Ã— more stars, 4.4Ã— better Ïƒ/âˆšN, and the ratio degrades 3.6Ã—. **The extra
stars are bought at a worse price than they are worth.**

---

## 5. Results

### 5.1 Injectionâ€“recovery

Uniform injections return identically zero (Â§3.2). Sparse recovery at N = 8Ã—10âµ:

| f | p injected | p recovered | ratio |
|---|---|---|---|
| 0.5 | 1eâˆ’2 | 9.48eâˆ’3 | 0.95 |
| 0.5 | 1eâˆ’4 | 9.84eâˆ’5 | 0.98 |
| 0.3 | 1eâˆ’3 | 6.23eâˆ’4 | 0.62 |
| 0.1 | any | consistent with zero | â€” |

**f = 0.1 is unrecoverable at any p**: Î” = 0.114 mag is 1.15Ïƒ and never reaches
the tail.

### 5.2 Blinding, and its failure

A secret offset (uniform Â±0.05 mag, SHA-256 committed) was added to all
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

| f | Î” (mag) | efficiency | p upper limit | mean fÌ„ |
|---|---|---|---|---|
| 0.03 | 0.033 | 0.008 | 0.770 | 0.0231 |
| 0.10 | 0.114 | 0.016 | 0.388 | 0.0388 |
| 0.30 | 0.387 | 0.169 | 0.0358 | 0.0107 |
| **0.50** | 0.753 | 0.986 | **6.13eâˆ’3** | **3.07eâˆ’3** |

The validation subset of 47,927 stars gave p_UL = 5.1 Ã— 10â»Â³; the full
3,321,566-star sample gives 6.1 Ã— 10â»Â³. **Sixty-nine times more stars, and the
limit did not improve** â€” the signature of a background that scales with N.

Variant B (grey-across-optical) reaches fÌ„ < 9.7 Ã— 10â»â´ with per-star reach
f â‰³ 0.2.

### 5.4 Candidates

At 5Ïƒ: **19,844 positive, 1,126 negative** â€” a 17.6:1 asymmetry at 129Ïƒ. Of the
300 strongest positive outliers:

| outcome | count |
|---|---|
| rejected by mid-IR excess or SIMBAD class | 22 |
| **no AllWISE photometry at all** | **278** |
| **measured-and-normal mid-IR, not a known contaminant** | **0** |

SIMBAD classes on the positive side are contaminant types (Em* Ã—5, Y*? Ã—1,
Y*O Ã—1); the negative control side draws ordinary stars (PM* Ã—17). That 278/300
lack WISE is itself diagnostic: extreme outliers sit in fields too confused for
a clean match, which is also why they are outliers.

**Zero surviving candidates.**

### 5.5 The contaminant, identified

The positive tail tracks every multiplicity proxy while the negative control
tail *shrinks*:

| proxy | positive tail | negative tail |
|---|---|---|
| RUWE | 2.04Ã— | 0.72Ã— |
| astrometric_excess_noise | 2.77Ã— | 0.85Ã— |
| **BP/RP excess factor C*** | **20.5Ã—** | 0.38Ã— |

C* dominating by 10Ã— over RUWE identifies the mechanism as **aperture
mismatch**: Gaia resolves sub-arcsecond, the 2MASS beam is ~4â€³, so a neighbour
Gaia separates is blended into the same Ks measurement â€” inflating Ks and
making the star look under-luminous in G.

Wide pairs catch this in the act and localise it:

| separation | n_pos | n_neg | ratio |
|---|---|---|---|
| 3.0â€“6.5â€³ | 3 | 152 | **50.7** |
| 6.5â€“10.3â€³ | 3 | 8 | 2.7 |
| 10.3â€“19.3â€³ | 8 | 11 | 1.4 |
| 19.3â€“120â€³ | 7 | 4 | 0.57 |

**This is structural, not a data-quality problem.** It limits any
optical-deficit search anchored on 2MASS, and removing it requires a NIR anchor
of comparable angular resolution.

### 5.6 The best limit: clean wide pairs

17,268 co-natal pairs from a self-join of the analysis sample (chance-alignment
contamination 0.02% by scramble test). Common-mode cancellation is real but
partial: Ïƒ(Î”r) = 0.1138 against 0.1402 for none at all, so **only 34% of the
main-sequence variance is common-mode**. Two-thirds is genuinely per-star.

Beyond 10â€³ the tails balance, so the background is *measurably* symmetric and
can be subtracted â€” which the single-star test can never do.

| | single star | **clean pairs (Î¸ > 10â€³)** |
|---|---|---|
| stars | 3,321,566 | **17,688** |
| 5Ïƒ counts | 19,844 vs 1,126 | **15 vs 18 (âˆ’3 Â± 5.7)** |
| p_UL | 5.9eâˆ’3 | **4.6eâˆ’4** |
| **mean fÌ„** | 2.95eâˆ’3 | **2.14eâˆ’4** |

**13.8Ã— better from 188Ã— fewer stars.** The gain is not statistical.

### 5.7 A limit on the beaming-consistent class

The comparison in Â§1.3 points at what this dataset can uniquely contribute.
Suazo et al. (2022) constrain spheres that re-radiate isotropically at an
assumed temperature. Nobody has put a number on the class that does *not*.

Of 3,321,566 stars, 2,993,277 (90.1%) have usable W1 and W2, and 2,897,066
(87.2%) have a **measured** bare photosphere: |excess| < 3Ïƒ in both bands
against a (Jâˆ’Ks)â‚€-calibrated photospheric colour. Absence of WISE data is never
counted as absence of excess â€” Â§5.4 showed that would fabricate candidates.

The veto turns out to remove very little of our tail:

| subsample | N | 5Ïƒ positive | rate | pos/neg |
|---|---|---|---|---|
| all stars | 3,321,566 | 19,844 | 0.00597 | 17.6 |
| bare photosphere | 2,897,070 | 13,611 | 0.00470 | 15.2 |
| IR excess present | 73,123 | 687 | 0.00940 | 21.5 |

The single-star limit improves by only **1.13Ã—**. That is itself a result:
**our optical-deficit tail is not made of dusty objects.** Stars with genuine
mid-IR excess carry a rate only 2Ã— the bare-photosphere rate and are 2.4% of
the sample. The tail is blends, as Â§5.5 established by a different route.

Applying the veto where our best estimator lives â€” clean wide pairs with
**both** components bare, 6,431 pairs / 12,862 stars:

| k | threshold (mag) | f detectable | n_pos | n_neg | asymmetry | p_UL | fÌ„ limit |
|---|---|---|---|---|---|---|---|
| 3 | 0.327 | 0.260 | 143 | 174 | âˆ’31 Â± 17.8 | 2.5eâˆ’3 | 6.5eâˆ’4 |
| 5 | 0.545 | 0.395 | 8 | 12 | âˆ’4 Â± 4.5 | 8.0eâˆ’4 | 3.2eâˆ’4 |
| **6** | 0.655 | **0.453** | 1 | 3 | âˆ’2 Â± 2.0 | **4.9eâˆ’4** | **2.2eâˆ’4** |

**p < 4.9 Ã— 10â»â´ at f â‰¥ 0.45 for stars with an optical deficit and a measured
absence of mid-infrared excess.** Every asymmetry is consistent with zero.

This is the only limit in this paper not superseded by prior work, and it is
complementary rather than competitive: it covers the corner of parameter space
that an IR-excess estimator weights to zero.

### 5.8 Deficit spectral type

Classifying each outlier by whether its colour excess matches the reddening
vector: positive outliers have median Î´C_obs/Î´C_dust = **1.09**. They redden.
There is no grey population and no exotic-spectrum population â€” nothing near
the Î± â‰ˆ 0.19 blind spot where a deliberate evader would sit.

### 5.8 Search for a spreading front

A spreading civilisation produces *spatially correlated* anomalies, and
correlation statistics detect clustered weak signals far below the per-object
threshold. We regress out all non-positional systematics with a model never
shown (l, b, d) â€” so it cannot absorb a front â€” and scan 31,184 cells of 25 pc.

**Real structure is found**: cell-mean scatter 1.554Ïƒ against a shuffled null of
1.001, 137 cells above 5Ïƒ where 40 are expected, and a two-point correlation of
the cleaned residual field significant at **19.5Ïƒ** on 25 pc scales.

It is dust. Anomalous cells sit at median |b| = 15.5Â° against 30.7Â° overall,
with median Aâ‚€ = 0.176 against 0.090. Restricting to |b| > 30Â° and Aâ‚€ < 0.15
the structure weakens but persists and *still* points at the residual dust.
The morphology is wrong for a front: 157 separate clumps, largest 22 cells
against a null of 5.4, with no edge.

**No front.** But the test detected a real 19Ïƒ correlated signal that every
single-star statistic here was blind to, and unlike the tail statistic it is
extinction-limited rather than background-limited â€” its sensitivity still
improves with N.

---

## 6. Discussion

### 6.1 Answer to the motivating estimate

The Ïƒ/âˆšN argument fails for three independent reasons:

1. It does not apply to the mean, which is not measurable at all (Â§3.2).
2. Where it could apply it stops at 2â€“5 Ã— 10â»Â² mag, **957â€“1376Ã— above naive**,
   on coherent systematics that do not average down (Â§4).
3. The tail-based limit is background-limited by companions, YSOs,
   emission-line stars and spotted rotators, and that background scales with N
   (Â§5.3).

**The optical-deficit channel saturates at fÌ„ ~ 2 Ã— 10â»â´ and cannot reach
10â»âµ.**

### 6.2 Comparison with prior work

Against [Zackrisson et al. (2015)](https://doi.org/10.1088/0004-637X/810/1/23):
their <3% is a limit on the fraction of *galaxies*, ours on the fraction of
*stars* within 500 pc. Not directly comparable, but both are limited by the
intrinsic width of the scaling relation regressed against â€” Tullyâ€“Fisher for
them, the main sequence for us â€” rather than by photon noise. Our relation is
several times tighter per object; our limiting factor is contamination.

Against [Suazo et al. (2022)](https://doi.org/10.1093/mnras/stac1789),
Hephaistos I â€” the direct numerical comparison:

| | covering fraction | limit on p |
|---|---|---|
| Suazo+22, 100 pc, 300 K | Î³ â‰¥ 0.5 | **1.9 Ã— 10â»â´** |
| **This work, clean pairs, 500 pc** | f â‰¥ 0.46 | 4.6 Ã— 10â»â´ |
| Suazo+22, 100 pc, 300 K | Î³ â‰¥ 0.9 | 1.8 Ã— 10â»âµ |

**Their limit is ~2.4Ã— stronger at comparable covering fraction, from ~12Ã—
fewer stars.** The reason is instructive: their estimator uses the infrared
re-emission as well as the obscuration, and infrared excess is measured against
a near-zero photospheric background whereas a deficit is measured against a
0.089 mag main-sequence width. Adding a second, near-background-free channel is
worth more than adding an order of magnitude in N â€” which is the same lesson
Â§4.2 teaches from the other direction.

The compensating argument for the deficit channel is coverage, not depth. A
limit derived from a model that includes Î³Â·BB(T_DS) constrains spheres that
re-radiate into 4Ï€ at the assumed temperature. A beamed radiator, a sphere cold
enough to fall outside WISE, or one that exports its energy non-thermally is
weighted differently by that estimator and not at all by ours. We cannot
quantify how much of that class Hephaistos retains without re-running their
pipeline with the re-emission term suppressed, and we do not claim a number
for it.

Against [Zackrisson et al. (2018)](https://doi.org/10.3847/1538-4357/aac386):
comparable per-object sensitivity (their f_cov > 0.75, our f â‰³ 0.46) from a
different observable, on 400Ã— fewer stars. Their single candidate was killed by
an unseen binary companion; we find binarity to be the dominant contaminant
from an independent direction. The agreement on the *failure mode* is stronger
evidence than either limit.

Against [Zackrisson et al. (2015)](https://doi.org/10.1088/0004-637X/810/1/23):
their <3% is a limit on the fraction of *galaxies*, ours on the fraction of
*stars*. Not directly comparable, but both are limited by the intrinsic width
of the scaling relation regressed against â€” Tullyâ€“Fisher for them, the main
sequence for us â€” rather than by photon noise.

### 6.3 What this does and does not constrain

The honest accounting. Within 500 pc, among lower-main-sequence stars, fewer
than 1 in 2,200 intercepts â‰¥46% of its optical output via a spectrally
selective absorber. Three structural blind spots remain: a **uniform** deficit
(unbounded by any self-calibrated method), an absorber at **Î± â‰ˆ 0.19**, and
**f â‰² 0.1** at any abundance.

That last is the one that matters for interpretation. A civilisation taking 1%
of its star â€” still ~10Â²â´ W, vastly more than any plausible need â€” is entirely
invisible to this method. **We constrain profligate near-total conversion, not
sipping.**

In terms of hypothesis space: the searchable (f, p) plane spans roughly 21
decadeÂ², and this work excludes roughly 0.8 of them, i.e. **~4â€“8%**. A null
result at this sensitivity is genuine but weak evidence, because the hypothesis
space extends decades below any version of this method's reach.

One thing *is* meaningfully constrained: were we inside an active conversion
front, we would see it locally regardless of volume fraction. We do not.

### 6.4 The programme this points to

The comparison with Suazo et al. (2022) contains its own instruction. They beat
us by 2.4Ã— using ~12Ã— fewer stars because they combine two channels, one of
which is measured against a near-zero background. The route to a competitive
deficit-channel limit is therefore *not* more stars â€” Â§4.5 shows that fails
outright â€” but a better estimator on the stars we have.

Concretely, in the order we would run them:

**(i) Scale the wide-pair differential.**
[El-Badry, Rix & Heintz (2021)](https://doi.org/10.1093/mnras/stab323) catalogue
~1.3 Ã— 10â¶ Gaia pairs to 1 kpc against our 8,844. The pair estimator is the
only one here whose background is *measurably* symmetric and therefore
subtractable, and its limit scales roughly as 1/N rather than 1/âˆšN in the
background-subtracted regime. Reaching even 10âµ clean pairs would put
p ~ 10â»âµâ€“10â»â´ at f â‰ˆ 0.45, i.e. **competitive with or better than the
isotropic-case limits, in a class those limits do not cover.** This is the
single highest-value next step and it needs no new observations.

**(ii) Extend the stellar mass range.** We deliberately restricted to
3 < M_Ks < 8 (roughly G2Vâ€“M4V). More luminous stars are both more attractive
engineering targets and better represented in the comparison samples. The
turnoff and subgiant contamination that motivated our cut is tractable with
Gaia parallaxes and log g.

**(iii) Then, and only then, the joint constraint.** With a deficit limit below
the IR-excess limit at matched covering fraction, the *difference* between the
two bounds the fraction of harvesters whose waste heat is beamed, cold, or
non-thermal â€” Î©/4Ï€ in the language of Â§1.1. No one has measured that, it is the
quantity the whole beaming argument turns on, and it requires exactly one
channel to be pushed below the other. Step (i) is what makes it possible.

**(iv) A NIR anchor at Gaia-like resolution** (VISTA VHS/VVV, UKIDSS) attacks
the dominant contaminant at its cause (Â§5.5) rather than filtering on a proxy
for it, which Â§6.4 below shows saturates at 1.9Ã—.

### 6.5 Where the remaining leverage is

In measured order of value:

1. **A NIR anchor at Gaia-like resolution** (VISTA VHS/VVV, UKIDSS). Â§5.5 shows
   this is the physical cause of the dominant contaminant, and it is the only
   lever that plausibly moves the *tail fraction* rather than filtering on a
   proxy for it. Hard cuts on the proxies saturate at 1.9Ã— because the residual
   distribution is **scale-free** â€” cleaning shrinks Ïƒ and the tail together.
2. **More clean wide pairs.** [El-Badry, Rix & Heintz (2021)](https://doi.org/10.1093/mnras/stab323)
   have ~1.3 Ã— 10â¶ pairs to 1 kpc against our 8,844.
3. **Gaia XP spectra.** Stop looking for missing flux and look for the wrong
   spectral *shape*. Every background here changes the SED smoothly, within the
   space stellar atmospheres already span; an engineered absorber plausibly has
   a band edge, which no atmosphere produces. This also *measures* Î± instead of
   being blind at 0.19.
4. **Eclipsing binaries and asteroseismology** â€” the only route to bounding the
   uniform case, since R and T_eff give an absolute luminosity prediction with
   no main-sequence-width systematic.

### 6.6 Errors we made

Recorded because they are instructive and because a reader should be able to
judge the pipeline's reliability.

- The **parallax zero point was applied 1000Ã— too small** (`gaiadr3-zeropoint`
  returns mas, not Âµas). Silent â€” the diagnostic printed "âˆ’0.0 Âµas", which
  reads as negligible rather than broken. Caught only because it contradicted
  the literature value quoted two lines above in the same docstring.
- **Gaia band wavelengths were DR2-era values**, mis-cited as DR3.
- We assumed hard cuts on contamination proxies would buy ~10Ã—, from a real
  20Ã— dependence of the tail rate on C*. Measured: **1.9Ã—**, because the tail is
  scale-free. Trusting the 20Ã— would have been wrong by an order of magnitude.
- We expected wide-binary common-mode cancellation to be large. It is **34%**.
- The blinding scheme did not blind the primary statistic (Â§5.2).
- **We initially framed this as the first stellar optical-deficit search and as
  setting the best limit. Both were wrong** â€” Zackrisson et al. (2018) precedes
  it methodologically and Suazo et al. (2022) already had a stronger limit at
  comparable covering fraction. The literature review that found this was done
  after the analysis rather than before it, which is the wrong order and is
  exactly how a novelty claim gets overstated. Â§1.3 and Â§6.2 carry the
  correction.

---

## 7. Conclusions

1. The optical-deficit channel is **background-limited, not statistics-limited**.
   The naive Ïƒ/âˆšN sensitivity overstates the truth by **957â€“1376Ã—**. This is
   the main result and we believe it is new.
2. Our constraint is **fÌ„ < 2.1 Ã— 10â»â´**, **p < 4.6 Ã— 10â»â´** at f â‰¥ 0.46, from
   clean wide pairs, with **zero candidates**. This is *not* the strongest
   published limit: Suazo et al. (2022) reach 1.9 Ã— 10â»â´ at Î³ â‰¥ 0.5 using the
   infrared re-emission as well.
3. **A uniform deficit is unmeasurable** by any self-calibrated version of this
   method.
4. **Unresolved companions mimic the signal**, because dM_G/dM_Ks = 1.26 > 1.
5. **The blind spot is Î± â‰ˆ 0.19, not a grey absorber**; grey gives the opposite
   sign.
6. The dominant contaminant is **Gaia/2MASS aperture mismatch**, localised to
   Î¸ â‰² 6.5â€³, and it limits any 2MASS-anchored version of this search.

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
- Suazo M., Zackrisson E., Wright J. T., Korn A. J., Huston M., 2022, MNRAS, 512, 2988 â€” *Project Hephaistos I: upper limits on partial Dyson spheres*
- Suazo M. et al., 2024, MNRAS, 531, 695 â€” *Project Hephaistos II*
- Villarroel B. et al., 2020, AJ, 159, 8 â€” *VASCO I*
- Zackrisson E., Korn A. J., Wehrhahn A., Reiter J., 2018, ApJ, 862, 21 â€” *SETI with Gaia: nearly complete Dyson spheres*
- Wang S., Chen X., 2019, ApJ, 877, 116
- Wright J. T. et al., 2014a, ApJ, 792, 26
- Wright J. T. et al., 2014b, ApJ, 792, 27
- Wright J. T., 2023, ApJ, 956, 34
- Zackrisson E. et al., 2015, ApJ, 810, 23

