# RESULTS — stellar optical-deficit search

**Status: complete and unblinded.** The full-sky 500 pc sample is
4,809,840 rows downloaded across 192 HEALPix partitions with zero failures,
**3,884,167 stars** after cuts, **3,321,566** in the fitted sample. Every
number below is from that sample unless explicitly labelled otherwise. The
analysis was frozen at commit `738942b` and unblinded at `0240acd`; the one
post-unblinding change is declared at the end.

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

### Distance trade study: pushing further out makes things worse

Same 26 sky partitions at every distance cut, same spline complexity, so the
comparison isolates distance (`scripts/17_distance_trade.py`,
`F9_distance_trade.png`):

| d limit (pc) | N | median `A_0` | `sigma/sqrt(N)` | spatial plateau | floor/naive |
|---|---|---|---|---|---|
| 200 | 67,061 | 0.033 | 3.28e-4 | 0.0154 | 149× |
| 300 | 172,126 | 0.061 | 2.35e-4 | 0.0157 | 287× |
| 500 | 480,700 | 0.110 | 1.42e-4 | 0.0227 | 345× |
| 750 | 945,934 | 0.144 | 9.41e-5 | 0.0230 | 339× |
| 1000 | 1,265,254 | 0.160 | 7.92e-5 | 0.0231 | 352× |
| 1250 | 1,383,542 | 0.166 | 7.51e-5 | 0.0235 | 536× |

Going from 200 pc to 1250 pc buys **21× more stars** and improves
`sigma/sqrt(N)` by 4.4×. Over the same range the median extinction rises 5×
and the spatial systematic plateau rises from 0.0154 to 0.0235 mag. The ratio
of floor to naive degrades by 3.6×.

**The extra stars are bought at a worse price than they are worth.** There is
no distance at which this reverses within the Edenhofer map's 1.25 kpc
validity, so 500 pc is not a compromise that a bolder cut would improve on.

The 500 pc row here (floor 0.0489 mag on 26 partitions) is consistent with the
full-sky 500 pc measurement (0.0520 mag on 192 partitions), which is a useful
check that the trade study's subset is not pathological.

---

## 4. Injection–recovery (step 4)

Run on an 800,000-star random subsample of the fitted sample, 6 realisations
per setting, with the spline hyperparameters held fixed. The subsample size is
recorded in the output; recovery significance scales as `p·sqrt(N)`, so these
thresholds are **not** rescaled to the full sample.

### Uniform injection: recovery is identically zero

| injected f | injected Δ (mag) | recovered mean residual | model-anchored recovery |
|---|---|---|---|
| 1e-5 | 0.0000109 | **0.00629139** | 0.00580 |
| 1e-4 | 0.0001086 | **0.00629139** | 0.00589 |
| 1e-3 | 0.0010863 | **0.00629139** | 0.00678 |
| 1e-2 | 0.0109120 | **0.00629139** | 0.01573 |
| 1e-1 | 0.1143937 | **0.00629139** | 0.10521 |

The self-calibrated column is **identical to eight decimal places** across four
orders of magnitude of injected signal, and the 5-sigma tail count is
identical too (4,846 in every row). The spline absorbs the offset completely.

The model-anchored column, comparing against a frozen reference curve, does
track the injection — recovering 0.105 for 0.1 injected and 0.0157 for 0.01 —
down to its own floor of ~0.006, set by the baseline offset of the reference.

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

| f | p injected | p recovered | ratio | 5σ tail count |
|---|---|---|---|---|
| 0.5 | 1e-2 | 9.476e-3 | 0.95 | 12,321 |
| 0.5 | 1e-3 | 9.681e-4 | 0.97 | 5,608 |
| 0.5 | 1e-4 | 9.84e-5 | 0.98 | 4,922 |
| 0.5 | 1e-5 | 1.14e-5 | 1.14 | 4,853 |
| 0.3 | 1e-2 | 6.381e-3 | 0.64 | 5,711 |
| 0.3 | 1e-3 | 6.230e-4 | 0.62 | 4,929 |
| 0.3 | 1e-4 | 5.89e-5 | 0.59 | 4,852 |
| 0.1 | any | consistent with zero | — | ~4,846 |

At `f = 0.5` recovery is essentially unbiased down to `p = 1e-5`. At `f = 0.3`
about 60% of the injected signal is recovered — the deficit `Δ = 0.387 mag` is
3.9 sigma, so a substantial fraction of harvested stars fall short of the
5-sigma threshold. **At `f = 0.1` nothing is recoverable at any `p`**, because
`Δ = 0.114 mag` is 1.15 sigma and never reaches the tail.

### The gap between injection sensitivity and the achievable limit

This is the single most informative pair of numbers in the report, and they
must not be confused:

| | value at `f = 0.5` |
|---|---|
| smallest **injected** `p` recovered at >3σ | **1e-5** (grid edge; true threshold may be lower) |
| model-free **upper limit** on `p` from the data | **6.13e-3** |
| **ratio** | **~600×** |

The injection test measures sensitivity to an *added* excess on top of a
baseline it takes as perfectly known — and it takes that baseline from the same
data. The upper limit cannot make that assumption: the 19,844 real 5-sigma
outliers are indistinguishable, star by star, from signal, so all of them must
be allowed to be signal.

**That factor of ~600 is the background limitation, quantified.** It is not a
statistical shortfall, and no amount of `N` reduces it. Closing it requires
removing the background — which §7b shows means a NIR anchor with angular
resolution comparable to Gaia's, not more stars.

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
unmodelled secondary light; it did not — it made the bound worse — which
indicates the residual scatter is *not* dominated by that term. Reporting the
failed hypothesis rather than only the subset that happened to win: the bound
stands at **f < 0.106** per star, and it is weak.

It also applies to a different population (binaries) than the science sample
(which excludes them by construction) — see LIMITATIONS §7. And note it is
bounded by the *scatter*, so a flat absorber affecting every star identically
is absorbed into the fitted relation and remains unbounded by any
self-calibrated method, mass anchor included.

---

## 6. Blinded analysis and limits

A secret offset was drawn (uniform ±0.05 mag), committed by SHA-256 to a
tracked file (`blind/commitment.json`), and added to every residual. The
analysis was frozen and committed at `738942b` with a clean working tree; the
unblinding script refuses to run otherwise and verifies the commitment hash
before revealing anything.

**Revealed offset: +0.036403 mag.** Commitment verified.

### The blind did not do its job, and that must be said

The tail statistics are defined relative to the sample median, so a constant
offset shifts the residuals and the median together and **the tail counts are
identical blinded and unblinded**. The exclusion curve computed under the blind
is bit-for-bit the one below. There was no moment at which the primary answer
was hidden. Only the mean-residual and model-anchored numbers were genuinely
blinded (mean +0.04270 blinded → +0.00630 unblinded).

What protects this result is **pre-registration, not blinding**: the 5-sigma
threshold, the conservative no-background-subtraction limit, the colour box and
the null-test battery were all fixed in `config.py` and `DECISIONS.md` before
the full sample existed — most before any data was downloaded. That is weaker
than a working blind. See DECISIONS.md D12a.

### Exclusion (95% CL), N = 3,321,566

| f | Δ (mag) | Δ/σ | detection eff. | p upper limit | implied mean `f̄` |
|---|---|---|---|---|---|
| 0.001 | 0.0011 | 0.011 | 0.0060 | 1.0 — unconstrained | — |
| 0.01 | 0.0109 | 0.110 | 0.0065 | 0.926 | — |
| 0.03 | 0.0331 | 0.334 | 0.0078 | 0.770 | 0.0231 |
| 0.10 | 0.1144 | 1.154 | 0.0156 | 0.388 | 0.0388 |
| 0.30 | 0.3873 | 3.907 | 0.169 | 0.0358 | 0.0107 |
| **0.50** | 0.7526 | 7.592 | 0.986 | **6.13e-3** | **3.07e-3** |
| 0.70 | 1.3072 | 13.19 | 1.000 | 6.04e-3 | 4.23e-3 |
| 0.90 | 2.5000 | 25.22 | 1.000 | 6.04e-3 | 5.44e-3 |

**Best constraint: `p < 6.1e-3` at `f = 0.5`, i.e. mean harvested fraction
`f̄ < 3.1e-3`.**

Below `f ≈ 0.03` the limit saturates — the test provides **no constraint at
all** on per-star harvested fractions below a few per cent, at any `p`.

**The limit is background-limited, and the full sample proves it.** The
validation subset gave `p_UL = 5.1e-3` on 47,927 stars; the full sample gives
`6.1e-3` on 3,321,566 stars — **69× more stars, and the limit got slightly
worse.** That is the signature of a background that scales with `N`. The
positive 5-sigma tail is 0.60% of the sample in both.

### Systematic budget, itemised

| term | mag |
|---|---|
| **astrophysical: main-sequence structure vs colour** | **0.06779** |
| extinction: low vs high `A_0` split | 0.04456 |
| photometric: bright vs faint G split | 0.02818 |
| crowding: sparse vs crowded split | 0.02612 |
| spatial coherence plateau | 0.01694 |
| extinction: implied fractional error × median `A_G` | 0.01314 |
| extinction: map vs Gaia per-star `A_G` | 0.00702 |
| parallax zero-point residual (10 µas, ×\|1−s\|) | 0.00195 |
| metallicity calibration | 0.00099 |
| extinction: band law (Fitz19 vs Wang & Chen 19) | 0.00005 |
| **quadrature sum** | **0.09258** |
| **largest single term** | **0.06779** |
| naive `sigma/sqrt(N)` | 0.0000544 |
| **budget / naive** | **1246×** |

The largest term is astrophysical, not instrumental. Extinction is second.
Everything Gaia-instrumental — zero point, metallicity, band law — is
negligible by comparison.

---

## 7. Outliers and their follow-up (step 7)

At 5 sigma: **19,844 positive (deficit-like), 1,126 negative (control)** — a
17.6:1 asymmetry at 129 sigma. Highly significant, and not a detection: §7b
identifies what it is.

Follow-up of the **300 strongest** positive outliers, with the 300 strongest
negative outliers as a control:

| outcome | count |
|---|---|
| rejected by mid-IR excess or SIMBAD contaminant class | 22 |
| **no AllWISE photometry at all — unverifiable** | **278** |
| **mid-IR measured and normal, not a known contaminant** | **0** |

SIMBAD classifications that matched, positive side: `Em*` (emission-line) ×5,
`Y*?` (candidate YSO) ×1, `Y*O` (YSO) ×1, `PM*` ×1, 279 unmatched. Negative
control side: `PM*` (ordinary high-proper-motion stars) ×17, `Em*` ×1,
`Y*?` ×1, 272 unmatched. The positive side draws contaminant classes; the
control side draws ordinary stars.

**278 of 300 have no WISE photometry at all.** That is itself diagnostic: the
extreme outliers sit preferentially in fields too confused for AllWISE to
produce a clean match — which is also why they are outliers. Counting "no
excess detected" as "beaming-consistent" would have promoted 278 candidates
instead of 0.

## **Zero surviving candidates.**

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
2. **Where `sqrt(N)` scaling could apply, it stops at 2–5e-2 mag** —
   **957× above the naive expectation for variant A and 1376× for variant B**
   — and the plateau is set by coherent astrophysical and extinction
   systematics that do not average down (§3). Variant B's naive `sigma/sqrt(N)`
   is 1.5e-5 mag, so the 1e-5 target *is* reachable statistically; it is three
   orders of magnitude away systematically.
3. **The tail-based limit is background-limited**, set by unresolved
   companions blended into the 2MASS beam, YSOs, emission-line stars and
   spotted rotators — all of which dim in the optical for ordinary reasons.
   That background scales with `N`, so more stars do not help (§6, §7, §7b).

**The optical-deficit channel saturates at a mean harvested fraction of
`f̄ < 3.1e-3` and cannot reach 1e-5.** It is two and a half orders of magnitude
short, and the shortfall is not statistical.

The cleanest single piece of evidence: the validation subset of 47,927 stars
gave `p_UL = 5.1e-3`; the full sample of 3,321,566 stars gave `6.1e-3`.
**Sixty-nine times more stars, and the limit did not improve.** Under
`sqrt(N)` scaling it should have improved eightfold.

The correct answer to the brief's question is the one it anticipated might be
true, and it is worth stating plainly: this channel's value is not sensitivity.
It is that it cannot be evaded by beaming, and that it constrains a region of
parameter space — cold, non-radiating, or directionally-beamed absorbers —
that the infrared surveys cannot see at all.

### What would actually improve it

Not more stars, and not a better dust map. In order of leverage:

1. **A near-infrared anchor at Gaia-like angular resolution.** §7b shows the
   dominant contaminant is neighbours blended into the 4-arcsec 2MASS `Ks`
   beam. This is worth a large factor and nothing else comes close.
2. **A mass anchor that is not a luminosity.** Asteroseismic or dynamical
   masses for single main-sequence stars would remove the flat-absorber
   degeneracy (§5) and the `M_Ks`-blending degeneracy at once.
3. **A metallicity that is not derived from the same photometry** (§5 of
   LIMITATIONS), which would let the astrophysical colour term — the single
   largest budget entry at 0.068 mag — be modelled rather than absorbed.

Improving the extinction treatment, which is where the effort naturally goes,
is fourth: it is worth 0.045 mag against the astrophysical term's 0.068 mag,
and the trade study shows the sample cannot be grown to compensate.

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

## Post-unblinding changes

The brief requires these to be declared with both versions. There is one.

**1. Substring-collision bug in the systematic-budget table
(`scripts/22_systematic_budget.py`).** The budget looked splits up by substring,
and `"extinction quartile"` also matches `"colour split WITHIN
lowest-extinction quartile"`. The extinction row therefore reported the
(larger) colour value. Fixed to exact-name matching, which now raises rather
than silently mismatching.

| | before fix | after fix |
|---|---|---|
| extinction: low vs high `A_0` split | 0.06779 mag *(wrong — colour value)* | **0.04456 mag** |
| quadrature sum | 0.10574 mag | **0.09258 mag** |
| largest single term | 0.06779 mag | 0.06779 mag (unchanged) |
| budget / naive | 1246× | 1246× (unchanged) |

This is a reporting bug in a summary table, not an analysis change: the
underlying null-test values were computed before unblinding and are unaltered,
and the headline numbers (largest term, ratio) are identical either way. The
before-fix numbers are recorded above so the correction is auditable.

**No analysis choice, cut, threshold or model was changed after unblinding.**
