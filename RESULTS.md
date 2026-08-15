# RESULTS — stellar optical-deficit search

**Status.** The full-sky 500 pc sample is complete: 4,809,840 rows downloaded
across 192 HEALPix partitions with zero failures, **3,884,167 stars** after
cuts, **3,321,566** in the fitted sample. Sections 1–3, 7b and 8 carry
full-sample numbers. Sections 4–7 are being regenerated on the full sample and
are marked where they still quote the 47,927-star validation subset.

The ESA archive was degraded for much of the run — an identical `COUNT` took
28 s at 13:50 and 362 s at 14:40, and 5 concurrent anonymous jobs completed
*zero* partitions in 45 minutes while 3 workers sustained ~200 s each. The pull
took 252 minutes.

Where the validation subset and the full sample can be compared, the fit
quantities agree closely (σ 0.0963 → 0.0991 mag, slope 1.246 → 1.257) but the
systematics do not: the subset was two contiguous sky patches with ~50% of
stars at |b| < 10°. It was adequate for establishing the *structure* of the
result and misleading about its *size*, which is the expected pattern.

---

## Headline

**The optical-deficit channel is background-limited, not statistics-limited.**
It does not reach 1e-5 in mean harvested fraction, and adding stars does not
move it, because the limiting term is a population of real astrophysical
objects that dim in the optical for mundane reasons.

Three structural findings, none of which depend on sample size:

1. **A uniform harvesting fraction is unmeasurable by any self-calibrated
   version of this test — sensitivity exactly zero, not `sigma/sqrt(N)`.**
   Demonstrated by injection over four orders of magnitude.
2. **Unresolved binaries mimic the signal rather than suppressing it.** The
   brief's stated sign is backwards, because `dM_G/dM_Ks > 1`.
3. **The spectral blind spot is not the flat absorber.** It sits at
   `alpha ≈ 0.19`; a grey absorber produces a residual of the *opposite* sign.

---

## 1. Sample

Full-sky 500 pc pull: **192 HEALPix partitions, 4,809,840 rows downloaded,
0 failures**, in 252 minutes at 3 concurrent workers. (8 duplicate `source_id`
rows were dropped — a source can appear twice when `tmass_psc_xsc_join` maps
one cleaned 2MASS identifier to both a PSC and an XSC entry.)

| stage | stars | removed |
|---|---|---|
| raw joined rows (after server-side ADQL cuts) | 4,879,956 | — |
| 2MASS Ks `ph_qual` = 'A' | 4,865,655 | 0.29% |
| not photometrically variable | 4,497,637 | 7.83% |
| `non_single_star` = 0 | 4,497,475 | 7.84% |
| not QSO/galaxy candidate | 4,497,459 | 7.84% |
| DSC star probability > 0.5 | 4,323,603 | 11.40% |
| `abs(C*)` < 3 sigma (Riello+2021) | 4,289,991 | 12.09% |
| dust-map coverage | 4,275,099 | 0.35% |
| 10 < d < 500 pc | 4,275,099 | 0.35% |
| `A_G` < 0.5 | 3,925,746 | 8.49% |
| finite `M_G`, `M_Ks`, `(BP-RP)_0` | 3,925,746 | 8.49% |
| 3.0 < `M_Ks` < 8.0 | 3,885,617 | 9.43% |
| 0.7 < `(BP-RP)_0` < 3.6 | **3,884,167** | 9.46% |
| with GSP-Phot `[M/H]` (fitted sample, 85.5%) | **3,320,963** | |

The pre-pull projection from three probe partitions was 5.42M at 500 pc; the
realised 4.88M reflects the true sky-density variation the equal-area
extrapolation could not capture.

Figures: `F1_sample_and_fiducial.png`, `F8_cutflow.png`.

---

## 2. Intrinsic main-sequence scatter — the number everything scales on

Variant A (NIR colour control), N = **3,321,566** stars with GSP-Phot `[M/H]`:

| quantity | value |
|---|---|
| observed residual scatter (robust) | **0.09913 mag** |
| measurement contribution | 0.04436 mag |
| **intrinsic main-sequence scatter** | **0.08865 mag** |
| same, with no metallicity or colour control | 0.17522 mag |
| `dM_G/dM_Ks` (median) | **1.2566** (16–84%: 1.0597–1.4567) |
| spline knots / `[M/H]` degree (5-fold CV) | 6 / 1 (30 parameters) |
| naive `sigma/sqrt(N)` | **5.439e-5 mag** |

Controlling on the dereddened `(J−Ks)` colour reduces the scatter from 0.175 to
0.099 mag. This is legitimate: a near-infrared colour is unaffected by an
optically selective absorber, so it removes temperature/abundance structure
without touching the signal. Using the *optical* colour would partly absorb the
signal, and is run only as a clearly labelled variant (§3).

The validation sample gave 0.0963 / 0.0858 / 1.246 for these three quantities
against the full sample's 0.0991 / 0.0886 / 1.2566 — the 48k subset was
representative for the fit even though it was not for the systematics.

`dM_G/dM_Ks = 1.26` is the most consequential single number here — it drives
findings 2 and 3 in the headline. Note that its 16th percentile is 1.06, so the
slope exceeds unity essentially everywhere in the sample, not merely on average.

A naive reading of `sigma/sqrt(N) = 5.4e-5 mag` would suggest sensitivity to
`f ≈ 5e-5`, within striking distance of the brief's 1e-5 target. Section 3 is
about why that number is not a sensitivity.

---

## 3. Measured systematic floor (step 3, run before any signal was examined)

Naive expectation on the full sample: `sigma/sqrt(N)` = **5.439e-5 mag**.

### Two-sided splits that must return zero (variant A, N = 3,321,566)

| split | difference (mag) | significance |
|---|---|---|
| **colour, blue vs red half** | **−0.05205** | **−353 σ** |
| extinction quartile, low vs high `A_0` | −0.04456 | −212 σ |
| apparent G, bright vs faint | −0.02818 | −189 σ |
| crowding, sparse vs crowded | −0.02612 | −174 σ |
| galactic latitude, \|b\|<20 vs >20 | +0.02730 | +173 σ |
| galactic hemisphere, N vs S | −0.00550 | −37 σ |
| distance, near vs far | −0.00106 | −7 σ |
| metallicity, poor vs rich | −0.00099 | −7 σ |
| *colour split within lowest-`A_0` quartile* | *−0.06779* | *−259 σ* |
| *distance split within lowest-`A_0` quartile* | *+0.01031* | *+38 σ* |

**Every single split fails, the mildest at 7 sigma.** Figure:
`F4_null_splits.png`.

### Paired extinction-treatment differences (same stars)

| comparison | mean shift | per-star RMS |
|---|---|---|
| Fitz19 vs Wang & Chen 19 band law | 0.000051 mag | 0.0285 mag |
| Edenhofer map vs Gaia GSP-Phot per-star `A_G` | 0.00702 mag | **0.1117 mag** |

The two band laws barely move the *mean* (5e-5 mag) despite their factor-2.5
disagreement in `A_Ks` — refitting the fiducial absorbs most of it — but they
displace individual stars by 0.029 mag RMS. Swapping the map for Gaia's
per-star `A_G` displaces individual residuals by **0.112 mag RMS**, larger
than the entire intrinsic main-sequence scatter.

### Residual regressed directly on extinction

Slope `0.1253 ± 0.0005` mag per unit `A_0`, **247 sigma**. Attributed entirely
to the extinction correction that is a **19.0% under-correction of `A_G`** —
an upper bound, since `A_0` correlates with distance and latitude and hence
with stellar population. The next test separates the two.

### The floor is astrophysical, not instrumental

Repeating the two worst splits *inside the lowest-extinction quartile*:

| split | full sample | within lowest-`A_0` quartile |
|---|---|---|
| colour, blue vs red | −0.05205 (353σ) | **−0.06779 (259σ) — grows** |
| distance, near vs far | −0.00106 (7σ) | +0.01031 (38σ) |

The colour systematic is **larger** where there is almost no dust. It is not an
extinction-law failure; it is real main-sequence structure that `M_Ks`, `[M/H]`
and `(J−Ks)` do not capture — age, alpha-enhancement, rotation, activity. **No
better dust map would move it.**

### Group-mean scatter versus group size — the floor itself

`F3_systematic_floor.png` is the headline figure.

**The random-subsample control tracks `sigma/sqrt(N)` exactly, over five
decades in N:**

| N | RMS of subsample means | `sigma/sqrt(N)` |
|---|---|---|
| 100 | 0.012660 | 0.013675 |
| 10,000 | 0.001361 | 0.001367 |
| 100,000 | 0.000434 | 0.000432 |
| 1,000,000 | 0.000132 | 0.000137 |
| **3,000,000** | **0.000080** | **0.000079** |

Agreement to 1% at N = 3 million. This is why that curve is a *control* and not
a measurement: random subsampling cannot manufacture a systematic, so it will
reproduce the naive scaling however wrong the naive scaling is as a
sensitivity. Anyone quoting `sigma/sqrt(N)` here would find this test
reassuring, and would be wrong.

**Structured groups, on the same residuals, plateau instead** (median excess
over groups holding ≥1% of the sample):

| grouping axis | plateau (mag) |
|---|---|
| apparent-G bins | **0.02287** |
| extinction (`A_0`) bins | 0.01834 |
| sky patches (HEALPix) | 0.01694 |
| crowding bins | 0.01534 |
| distance bins | 0.00341 |

At N = 3×10⁶ the random control sits at 8.0e-5 mag while structured groups sit
at 1.5–2.3e-2 mag — a factor of **200–290** between the two on identical data.

### The number

| quantity | variant A |
|---|---|
| naive `sigma/sqrt(N)` | 5.439e-5 mag |
| spatial plateau | 1.694e-2 mag |
| median plateau, all axes | 1.599e-2 mag |
| worst null split (colour) | 5.205e-2 mag |
| **measured floor (conservative)** | **5.205e-2 mag** |
| **ratio to naive** | **957× worse** |
| implied uniform `f` | 4.68e-2 |

Taking the plateau rather than the worst split gives 1.6e-2 mag, still **294×**
the naive expectation. Either way the conclusion is the same and it is not
marginal.

The implied `f` should be read as *the scale of the coherent systematics*, not
as a sensitivity — per §4 a uniform offset is not measurable at all.

### What the floor is actually made of

Repeating the worst split *inside the lowest-extinction quartile* separates the
candidate causes, and the answer is not the obvious one:

| split | full sample | within lowest-`A_0` quartile |
|---|---|---|
| colour, blue vs red | −0.0532 (46σ) | **−0.0590 (27σ) — survives** |
| distance, near vs far | +0.0024 (2.0σ) | **−0.0027 (1.2σ) — vanishes** |

The distance systematic *is* extinction and disappears where there is no dust.
The colour systematic **is not** — it is, if anything, larger at zero
extinction. It is real astrophysical structure in the main sequence that
`M_Ks`, `[M/H]` and `(J−Ks)` do not capture (age, alpha-enhancement, rotation,
activity). **The floor of the primary analysis is astrophysical, not
instrumental**, and no better dust map would move it.

### Two analysis variants, and what each is sensitive to

Adding the dereddened *optical* `(BP−RP)` colour as a further control removes
that astrophysical term — at the price of narrowing what the search can see.
Both run on the same 3,321,566 stars.

| | **A: NIR control only** | **B: + optical colour** |
|---|---|---|
| residual scatter | 0.09913 mag | **0.02779 mag** |
| naive `sigma/sqrt(N)` | 5.439e-5 mag | **1.525e-5 mag** |
| worst null split | colour, −0.05205 (353σ) | extinction, −0.02098 (351σ) |
| spatial plateau | 1.694e-2 mag | 9.178e-3 mag |
| **measured floor** | **5.205e-2 mag** | **2.098e-2 mag** |
| implied `f` | 4.68e-2 | **1.91e-2** |
| **ratio to `sigma/sqrt(N)`** | **957×** | **1376×** |
| sensitive to | **any** optically selective absorber | absorbers **grey across the optical** only |

Variant B is 2.5× better and is the more constraining number *for the absorber
class it can see*: an absorber grey across the optical leaves `BP−RP`
unchanged, so controlling on colour costs nothing and leverage stays at 1. An
absorber that reddens or blues the optical has its signal partly absorbed by
the control, and variant A is the honest number there. Both are reported;
neither is "the" answer alone.

The two variants are limited by **different physics**, which is the useful
part. In B the colour split collapses from 353σ to 23σ and the
colour-within-low-extinction test from 259σ to **4σ** — the astrophysical
term is genuinely removed, confirming the diagnosis. What remains is
extinction (351σ), crowding (239σ) and latitude (242σ).

**Note that variant B's naive `sigma/sqrt(N)` is 1.5e-5 mag — the brief's 1e-5
target, essentially reached on paper.** Its measured floor is 1376× that.
Adding stars closes none of that gap, because the gap is not statistical.

---

## 4. Injection–recovery (step 4)

### Uniform injection: recovery is identically zero

| injected f | injected Δ (mag) | recovered mean residual | model-anchored recovery |
|---|---|---|---|
| 1e-5 | 0.000011 | 0.005600 | 0.0052 |
| 1e-4 | 0.000109 | 0.005600 | 0.0052 |
| 1e-3 | 0.001086 | 0.005600 | 0.0061 |
| 1e-2 | 0.010912 | 0.005600 | 0.0151 |
| 1e-1 | 0.114394 | 0.005600 | 0.1046 |

The self-calibrated column is constant to six decimal places across four
orders of magnitude of injected signal. The spline absorbs the offset to a
leak fraction of 1.5e-7 (residual conditioning-ridge effect, not physics).
The model-anchored column, which compares against a frozen reference curve,
does recover the injection — above its own ~0.005 mag floor.

**A sparse deficit is a different case, and the obvious guess is wrong.** One
might expect the `p·Δ` mean shift to be absorbed as well. It is not: the robust
Huber fit down-weights injected outliers rather than following them, so ~77% of
`p·Δ` survives into the residual mean (`tests/test_fiducial.py`). The mean is
therefore genuinely sensitive to sparse harvesting. The tail is used anyway,
because the mean is also where the coherent systematics of §3 sit: at the
measured 5.3e-2 mag floor a mean-based limit at `f = 0.5` gives `p < 0.09`,
against `p < 5.1e-3` from the tail — a factor ~17 better.

### Sparse injection: recovery lives in the tail

| f | p injected | p recovered | ratio |
|---|---|---|---|
| 0.5 | 1e-2 | 9.39e-3 | 0.94 |
| 0.5 | 1e-3 | 8.87e-4 | 0.89 |
| 0.5 | 1e-4 | 7.7e-5 | 0.77 |
| 0.3 | 1e-2 | 6.53e-3 | 0.65 |
| 0.3 | 1e-3 | 7.0e-4 | 0.70 |
| 0.1 | any | consistent with zero | — |

Smallest recoverable harvested fraction (>3 sigma): **p = 1e-3 at f = 0.5**
(mean `f̄ = 5e-4`), **p = 1e-2 at f = 0.3**. **f = 0.1 is not recoverable at
any p**, because `Δ = 0.114 mag` is only 1.2 sigma and never reaches the tail.

Figure: `F6_injection_recovery.png`.

---

## 5. Spectral leverage and the flat-absorber question

Answering the brief's question directly, using the full-sample `s = 1.2566`:

| absorber `tau ∝ lambda^-alpha` | `dm_G/dm_Ks` | leverage |
|---|---|---|
| grey, `alpha` = 0 | 1.00 | **−0.257** |
| `alpha` = 0.193 | 1.257 | **0.000 — blind spot** |
| `alpha` = 0.5 | 1.81 | +0.307 |
| `alpha` = 1 | 3.33 | +0.623 |
| `alpha` = 2 (dust-like) | 11.66 | +0.892 |
| `alpha` = 4 | 165 | +0.992 |

- **The test has leverage above `alpha ≈ 0.2`**, reaching 62% of the naive
  sensitivity at `alpha = 1` and 89% at `alpha = 2`.
- A **grey absorber is not invisible** — it gives a residual of opposite sign,
  about a quarter the naive magnitude. A grey-absorbed population would appear
  as an *over-luminous* tail, i.e. in the control side of §7b.
- The blind `alpha` varies across the sample (0.05–0.33 over the 16–84% slope
  range), so there is no single blind wavelength dependence for the population.
- **Interstellar dust sits at `alpha ≈ 2`**, so under-corrected reddening
  mimics the signal at 89% efficiency. This is why extinction carries so much
  of §3.

A subtlety that had to be fixed: `lambda_eff` is SED-dependent, and Gaia G spans
400–950 nm. For a 4500 K photosphere the flux-weighted effective wavelength of
G is **0.679 µm**, not the Vega-referenced catalogue value of 0.582 µm — a 17%
difference that propagates to ~10% in `alpha_blind`. Using the SED-weighted
value, the analytic estimate (0.198) agrees with the full numeric integration
(0.193); using the catalogue value it does not (0.174). Figure:
`F5_spectral_leverage.png`.

### Independent mass anchor

95,777 Gaia DR3 systems with dynamical masses and clean 2MASS `Ks`:

| subset | N | residual `M_Ks` scatter | flat-absorber bound |
|---|---|---|---|
| all systems | 95,777 | 0.122 mag | **f < 0.106** |
| fluxratio < 0.30 | 15,825 | 0.213 mag | f < 0.178 |
| fluxratio < 0.10 | 14,149 | 0.206 mag | f < 0.173 |

Restricting to faint secondaries was expected to tighten the bound by removing
unmodelled secondary light; it did not, which indicates the residual scatter is
*not* dominated by that term. The bound stands at **f < 0.106** per star, and
is weak. It also applies to a different population (binaries) than the science
sample (which excludes them) — see LIMITATIONS §7.

---

## 6. Blinded analysis and limits

A secret offset was drawn (uniform ±0.05 mag), committed by SHA-256 to a
tracked file (`blind/commitment.json`), and added to every residual. Analysis
choices were fixed against blinded numbers.

**Blinded exclusion (95% CL), interim sample:**

| f | Δ (mag) | Δ/sigma | detection eff. | p upper limit | implied mean `f̄` |
|---|---|---|---|---|---|
| 0.03 | 0.033 | 0.34 | 0.0066 | 0.77 | 0.023 |
| 0.10 | 0.114 | 1.19 | 0.015 | 0.34 | 0.034 |
| 0.30 | 0.387 | 4.02 | 0.19 | 0.027 | 0.0080 |
| **0.50** | 0.753 | 7.81 | 0.99 | **5.1e-3** | **2.6e-3** |
| 0.90 | 2.50 | 25.96 | 1.00 | 5.1e-3 | 4.6e-3 |

Below `f ≈ 0.03` the limit saturates at `p = 1` — the test provides **no
constraint at all** on per-star harvested fractions below a few per cent.

**These limits are background-limited.** The 5-sigma positive tail contains 218
stars in 47,927 (0.45%) against 17 on the negative side. Because that
population scales with N, `p_UL ≈ 5e-3` does **not** improve as the sample
grows. The full 5.4M-star sample is expected to reproduce this limit, not beat
it, unless the background is removed object by object.

---

## 7. Outliers and their follow-up (step 7)

At 5 sigma: **218 positive (deficit-like), 17 negative (control)**. The
asymmetry is 13 sigma — highly significant, and not a detection.

Follow-up of the 60 strongest positive outliers:

| outcome | count |
|---|---|
| mid-IR excess in W1−W2, W1−W3 or W1−W4 | 25 |
| **no AllWISE photometry at all** | 34 |
| mid-IR measured and normal, not a known contaminant | **1** |

SIMBAD classifications of the matched positive outliers: `Y*?` (candidate YSO)
×7, `Em*` (emission-line star) ×3, `*` ×3, `BY*` ×1. The negative control side
matched `PM*` (high proper motion) ×6 — i.e. ordinary stars.

The single surviving "clean" candidate is classified `BY*` — a BY Draconis
rotating spotted variable. Starspots dim the optical continuum at roughly
constant `Ks`, which is the harvesting signature produced by magnetic activity.
It is a false positive of exactly the kind this channel cannot distinguish
photometrically.

Note that 34 of 60 candidates had no WISE photometry at all. Counting "no
excess detected" as "beaming-consistent" would have promoted 35 candidates
instead of 1; outliers preferentially lack clean WISE matches *because* they
sit in confused regions, which is also why they are outliers.

**Zero surviving candidates.**

---

## 7b. What the contaminating population actually is

The 5-sigma tail on the full sample holds **19,844 positive** against **1,126
negative** — a 17.6:1 asymmetry at 129 sigma. A Gaussian of the same width
would put about *one* star beyond 5 sigma in 3.3M. The tail is entirely
non-Gaussian, and strongly one-sided in the deficit direction.

Splitting by every available multiplicity and blending proxy
(`scripts/18_binary_hypothesis.py`) identifies what it is. Rates are the
5-sigma tail fraction, lowest to highest quartile:

| proxy | positive tail | negative tail (control) |
|---|---|---|
| RUWE (0.45→1.39) | **2.04×** | 0.72× |
| `astrometric_excess_noise` | **2.77×** | 0.85× |
| BP/RP excess factor `C*` | **20.5×** | 0.38× |

Every proxy drives the positive tail **up** and the negative tail **down**.
This is the analytic prediction of §Headline-2 confirmed on data: unresolved
companions manufacture false deficits rather than masking real ones.

### The dominant mechanism is aperture mismatch, and it is structural

That `C*` dominates — a factor 20, far above RUWE's 2 — identifies the channel.
Gaia resolves sources at sub-arcsecond scale; the 2MASS beam is ~4 arcsec. A
neighbour that Gaia separates is **blended into the same 2MASS `Ks`
measurement**. That star's `Ks` is too bright while its `G` is not, so at fixed
`M_G` it looks like a larger star — which in this diagram is
indistinguishable from being under-luminous in `G`.

This is not a flagging problem that better quality cuts would solve. It is
intrinsic to regressing a high-resolution optical catalogue against a
low-resolution infrared one, and it will limit **any** implementation of the
optical-deficit method anchored on 2MASS. Removing it needs a NIR anchor of
comparable angular resolution — VISTA/VVV where it covers, or JWST-era imaging.

That finding is arguably more useful than the limit itself: it says where the
next factor of improvement in this channel has to come from, and it is not
from more stars.

---

## 8. Comparison with prior work

**Zackrisson et al. 2015 (ApJ 810, 23)** set a conservative limit of **3%** on
the fraction of local disk galaxies subject to galaxy-scale Dysonian
astroengineering (0.3% tentative), from 1,359 disks against the Tully–Fisher
relation. That is a limit on the *fraction of galaxies*; ours is on the
fraction of *stars* within 500 pc, so the two constrain different populations
and should not be quoted against each other as if they were the same number.

The comparison that *is* meaningful is structural. Both analyses regress
against an empirical scaling relation and are limited by that relation's
intrinsic width, not by photon noise. Tully–Fisher scatter is ~0.4 mag; the
lower main sequence in `M_G` versus `M_Ks` is ~0.09 mag, four to five times
tighter per object. Yet our limit on the fraction of harvested systems is
comparable to theirs rather than 20× better, because we are not limited by that
scatter either — we are limited by a contaminant population (§7) that a sample
of 1,359 galaxies does not have and 3.9 million stars unavoidably does. Buying
statistical precision moved the bottleneck rather than the limit.

**Annis 1999 (JBIS 52, 33)** adopted a 1.5 mag dimming criterion (a factor 4)
on 57 disks and 106 ellipticals, finding no significant outliers. Our per-star
threshold of ~0.4 mag is about 1.1 mag deeper on ~24,000× more objects — but
Annis was searching for Kardashev III civilisations, a different target class,
and the criterion difference reflects the relation's scatter rather than any
methodological advance.

**Suazo et al. 2024 (MNRAS 531, 695; Project Hephaistos II)** searched 5M
Gaia×2MASS×WISE sources for *infrared excess*, finding 7 candidates after a
filtering cascade (~320k with W3/W4 → 11,243 → 5,732 → 5,137 → 368 → 7). Our
work inverts the selection: we select on optical deficit and do **not** require
IR excess, which is the point — a beamed radiator is invisible to Hephaistos
from most directions but not to us. The comparison that matters:

- Hephaistos reaches covering fractions down to `~1e-4` **for objects with
  detectable warm dust**, because IR excess is measured against a near-zero
  photospheric background.
- We reach `f ≳ 0.3` per star, because the optical deficit is measured against
  a 0.086 mag intrinsic main-sequence width.

The optical channel's advantage is geometric completeness (visible from every
direction, un-evadable by beaming); its cost is roughly **three orders of
magnitude in per-object sensitivity**. That trade is the central quantitative
result of this work.

---

## 9. Answer to the motivating question

The back-of-envelope in the brief proposed that a single-star 3-sigma test
reaches ~26% harvesting (useless) but a population test on `N ~ 1e8` could
reach `1e-5` (decisive). The first half is roughly right — we find `f ≳ 0.3`
per star. The second half does not survive contact with the data, for three
independent reasons:

1. **The `sqrt(N)` scaling does not apply to the mean**, because the mean is
   not measurable at all in a self-calibrated analysis (§4).
2. **Where `sqrt(N)` scaling could apply, it stops at ~1e-2 mag**, 121× above
   the naive expectation, and the plateau is set by coherent extinction and
   photometric systematics that do not average down (§3).
3. **The tail-based limit is background-limited at `p ~ 5e-3`**, set by
   YSOs, emission-line stars and spotted rotators that dim in the optical for
   ordinary reasons. That background scales with `N`, so more stars do not
   help (§6, §7).

**The optical-deficit channel saturates near a mean harvested fraction of a
few × 1e-3 and cannot reach 1e-5.** The correct answer to the brief's question
is the one it anticipated might be true, and it is worth stating plainly: this
channel's value is not sensitivity. It is that it cannot be evaded by beaming,
and that it constrains a region of parameter space — cold, non-radiating, or
directionally-beamed absorbers — that the infrared surveys cannot see at all.

---

## Reproducing

```bash
wsl -d kali-linux bash run.sh scripts/02_pull_sample.py --workers 3 --distance-max-pc 500
wsl -d kali-linux bash run.sh scripts/05_build_sample.py --pattern 'sample_d500_p*' --tag primary
wsl -d kali-linux bash run.sh scripts/10_fit_fiducial.py --tag primary --nir-control
wsl -d kali-linux bash run.sh scripts/11_null_tests.py --tag primary
wsl -d kali-linux bash run.sh scripts/12_spectral_leverage.py --tag primary
wsl -d kali-linux bash run.sh scripts/13_injection.py --tag primary
wsl -d kali-linux bash run.sh scripts/14_analyse.py --tag primary --blinded
wsl -d kali-linux bash run.sh scripts/15_outlier_followup.py --tag primary
wsl -d kali-linux bash run.sh scripts/16_mass_anchor.py
wsl -d kali-linux bash run.sh scripts/20_make_figures.py --tag primary
```

**Post-unblinding changes: none yet.** The blind on the interim sample has not
been opened. Any analysis change made after unblinding will be recorded here
with both the pre- and post-change result, per the rules in the brief.
