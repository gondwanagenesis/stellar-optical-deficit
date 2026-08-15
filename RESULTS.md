# RESULTS — stellar optical-deficit search

**Status: interim.** All numbers below come from a 2-HEALPix-partition
validation sample (47,927 stars after cuts). The full-sky 500 pc pull is in
progress; the ESA archive has been degraded throughout (a `COUNT` that took
28 s at 13:50 took 362 s at 14:40), and 5 concurrent anonymous jobs completed
zero partitions in 45 minutes while 3 workers sustained ~200 s/partition. Every
number is regenerable by re-running the numbered scripts against the full
sample; the section headers mark which conclusions are sample-size dependent
and which are not.

The validation sample is *not* representative: it is two contiguous sky
patches with ~50% of stars at |b| < 10°, so its extinction systematics are
worse than the full-sky sample's will be. It is adequate for establishing the
structure of the result, which is what the interim conclusions rest on.

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

| stage | stars |
|---|---|
| raw joined rows (after server-side ADQL cuts) | 224,662 |
| 2MASS Ks `ph_qual` = 'A' | 224,335 |
| not photometrically variable | 209,953 |
| `non_single_star` = 0 | 209,951 |
| not QSO/galaxy candidate | 209,951 |
| DSC star probability > 0.5 | 207,475 |
| `abs(C*)` < 3 sigma (Riello+2021) | 206,620 |
| dust-map coverage | 205,958 |
| 10 < d < 500 pc | 57,336 |
| `A_G` < 0.5 | 55,538 |
| finite `M_G`, `M_Ks`, `(BP-RP)_0` | 55,538 |
| 3.0 < `M_Ks` < 8.0 | 54,911 |
| 0.7 < `(BP-RP)_0` < 3.6 | **54,885** |
| with GSP-Phot `[M/H]` (fitted sample) | **47,927** |

Full-sky projection from measured partition counts: **5.42M** stars at 500 pc,
17.3M at 1250 pc. Figure: `F1_sample_and_fiducial.png`, `F8_cutflow.png`.

---

## 2. Intrinsic main-sequence scatter — the number everything scales on

| quantity | value |
|---|---|
| observed residual scatter (robust) | **0.0963 mag** |
| measurement contribution | 0.0437 mag |
| **intrinsic main-sequence scatter** | **0.0858 mag** |
| same, without the NIR colour control | 0.168 mag |
| `dM_G/dM_Ks` (median) | **1.246** (16–84%: 1.043–1.448) |
| spline knots / `[M/H]` degree (5-fold CV) | 6 / 1 |

Controlling on the dereddened `(J−Ks)` colour reduces the scatter from 0.139 to
0.096 mag. This is legitimate: a near-infrared colour is unaffected by an
optically selective absorber, so it removes temperature/abundance structure
without touching the signal. Using the *optical* colour would have absorbed
the signal and was not done.

`dM_G/dM_Ks = 1.25` is the most consequential single number here — it drives
findings 2 and 3 in the headline.

---

## 3. Measured systematic floor (step 3, run before any signal was examined)

Naive expectation on this sample: `sigma/sqrt(N)` = **4.40e-4 mag**.

### Two-sided splits that must return zero

| split | difference (mag) | significance |
|---|---|---|
| colour, blue vs red half | −0.0532 | −39 sigma |
| extinction quartile, low vs high | −0.0199 | −12 sigma |
| apparent G, bright vs faint | −0.031 (pre-NIR) | −19 sigma |
| distance, near vs far | −0.019 (pre-NIR) | −12 sigma |
| galactic latitude, |b|<20 vs >20 | +0.020 (pre-NIR) | +11 sigma |
| crowding, sparse vs crowded | −0.013 (pre-NIR) | −8 sigma |
| galactic hemisphere, N vs S | +0.009 (pre-NIR) | +5 sigma |

**Every single split fails.** Figure: `F4_null_splits.png`.

### Paired extinction-treatment differences (same stars)

| comparison | mean shift | per-star RMS |
|---|---|---|
| Fitz19 vs Wang & Chen 19 band law | 0.00016 mag | 0.022 mag |
| Edenhofer map vs Gaia GSP-Phot per-star `A_G` | 0.0082 mag | **0.112 mag** |

The per-star RMS of the map-vs-GSP-Phot comparison is comparable to the entire
intrinsic scatter.

### Group-mean scatter versus group size — the floor itself

`F3_systematic_floor.png` is the headline figure. Random subsamples track
`sigma/sqrt(N)` down to 7e-4 mag, exactly as they must — which is why that
curve is a control and not a measurement. Structured groups plateau:

| grouping axis | plateau (mag) |
|---|---|
| apparent-G bins | ~0.020 |
| sky patches (HEALPix) | ~0.019 |
| extinction bins | ~0.008 |
| crowding bins | ~0.008 |
| distance bins | ~0.005 |

**Measured floor: 5.3e-2 mag** (conservative, set by the worst split), or
~1.1e-2 mag taking the spatial plateau alone. Against a naive `sigma/sqrt(N)`
of 4.4e-4 mag, that is **121x worse**.

In harvested-fraction units the floor corresponds to **f ≈ 4.8e-2** on a
uniform offset — except that, per §4, a uniform offset is not measurable at
all, so this number is best read as the scale of the coherent systematics
rather than as a sensitivity.

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

| | **A: NIR control only** | **B: + optical colour** |
|---|---|---|
| residual scatter | 0.0963 mag | **0.0267 mag** |
| worst null split | colour, −0.0532 | extinction, −0.0159 |
| **measured floor** | **5.3e-2 mag** | **1.6e-2 mag** |
| implied `f` | 4.8e-2 | **1.5e-2** |
| ratio to `sigma/sqrt(N)` | 121x | 131x |
| sensitive to | **any** optically selective absorber | absorbers **grey across the optical** only |

Variant B is 3.3x better and is the more constraining number *for the absorber
class it can see*: an absorber that is grey across the optical leaves `BP−RP`
unchanged, so controlling on colour costs nothing and leverage stays at 1. An
absorber that reddens or blues the optical has its signal partly absorbed by
the control, and variant A is the honest number there. Both are reported;
neither is "the" answer alone.

Note that in variant B the colour split passes (2.4σ) and **extinction becomes
the limiting term** (32σ) — the two variants are limited by different physics,
which is why the ratio to `sigma/sqrt(N)` is ~125x in both cases despite the
3.6x difference in scatter.

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

Answering the brief's question directly, using the measured `s = 1.246`:

| absorber `tau ∝ lambda^-alpha` | `dm_G/dm_Ks` | leverage |
|---|---|---|
| grey, `alpha` = 0 | 1.00 | **−0.246** |
| `alpha` = 0.19 | 1.246 | **0.000 — blind spot** |
| `alpha` = 0.5 | 1.81 | +0.313 |
| `alpha` = 1 | 3.33 | +0.626 |
| `alpha` = 2 (dust-like) | 11.6 | +0.893 |

- **The test has leverage above `alpha ≈ 0.2`**, reaching 63% of the naive
  sensitivity at `alpha = 1` and 89% at `alpha = 2`.
- A **grey absorber is not invisible** — it gives a residual of opposite sign,
  about a quarter the naive magnitude. A grey-absorbed population would appear
  as an *over-luminous* tail.
- The blind `alpha` varies across the sample (0.04–0.32 over the 16–84% slope
  range), so there is no single blind wavelength dependence.
- **Interstellar dust sits at `alpha ≈ 2`**, so under-corrected reddening
  mimics the signal at 89% efficiency. This is why extinction dominates §3.

Analytic and SED-weighted numeric agree (0.189 vs 0.186). Figure:
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

## 8. Comparison with prior work

**Zackrisson et al. 2015 (ApJ 810, 23)** set a conservative limit of **3%** on
the fraction of local disk galaxies subject to galaxy-scale Dysonian
astroengineering (0.3% tentative), from 1,359 disks against the Tully–Fisher
relation. That is a limit on the *fraction of galaxies*, whereas ours is on the
fraction of *stars* within 500 pc. The two are not directly comparable, but the
per-system dimming sensitivities are: Zackrisson resolves galaxy-scale dimming
at the few-tenths-of-a-magnitude level, and so do we
(`Δ ≳ 0.4 mag`, i.e. `f ≳ 0.3`, before a star enters our tail). Twelve extra
years, a 3,500× larger sample and per-star rather than per-galaxy photometry
have **not** bought a better per-object dimming threshold, because both
analyses are limited by the intrinsic width of the relation they regress
against, not by photon noise.

**Annis 1999 (JBIS 52, 33)** adopted a 1.5 mag dimming criterion (a factor 4)
on 57 disks and 106 ellipticals. Our per-star threshold of ~0.4 mag is about
1.1 mag deeper, on ~300× more objects — but Annis was searching for
Kardashev III civilisations, a different target class.

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
