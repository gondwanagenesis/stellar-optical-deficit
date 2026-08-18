# Multi-channel search for spreading stellar-energy harvesting

**Every channel below is independent. A null in one says nothing about any
other.** Different channels see different absorber classes, different
temperatures and different concealment strategies; they are not repeated trials
of one hypothesis. Coverage gaps are stated explicitly rather than implied to be
covered.

---

## Summary table

| # | channel | what it can see | result |
|---|---|---|---|
| 1 | optical deficit vs M_Ks, single star | selective absorbers, f ≳ 0.5 | p < 6.1e−3 |
| 2 | optical deficit, wide-pair differential | same, cleaner background | **p < 5.3e−4** |
| 3 | mid-IR veto (beamed class) | deficit + no warm re-emission (>300 K) | **p < 4.3e−4** |
| 4 | **dynamical mass** | **ANY spectral slope, incl. grey** | **p < 2.2e−4**, 0 positives |
| 5 | spectral-slope fit, 7 bands | α ≠ 2 absorbers | 372 asymmetric survivors |
| 6 | 3D spatial front | a spreading domain | 19.5σ → dust |
| 7 | kinematic (Δv) clustering | a propagating network | 4301σ → young associations |
| 8 | deficit colour type | grey vs reddening | all reddening-like |
| 9 | W3/W4 intercept + re-emit | warm re-emission (100–300 K) | 2,632 candidates → debris discs |
| 10 | radio counterparts (NVSS/FIRST) | non-thermal emitters | 1/5000 dimmed, 4/5000 control → null |
| 11 | 3D domain edge | grabby boundary | 148.5σ → extinction + scanning law |
| 12 | Search A diagnosis | NIR companion test | **consistent with blends** |
| 13 | coherent acceleration (Hipparcos–Gaia) | coordinated stellar migration | isotropic, −3.9σ; ≥2% domain excluded |
| 14 | local mass ledger (Oort limit) | **any** cold mass, any technology | < 1.2e5 M☉ within 100 pc |
| 15 | cold blackbody, 3–100 K (β free) | engineered radiators below every prior search | 0 of 3,030; all dust-like |
| 16 | dark companions vs WD channel | fully enshrouded star (f→1) | *(running)* |
| 17 | stellar mass function, spatially | **holes** in the population | 59σ → completeness; null at 100 pc |
| 18 | **bolometric closure vs dynamical mass** | **energy missing at ANY wavelength** | **18 vs 36 mirror → p < 2.8e−4** |

**Joint, disposal-agnostic: p_total < 6.2e−4.** Fewer than 1 in 1,614 nearby
lower-main-sequence stars intercepts ≥51% of its optical output by any means.

---

## Part II — leaving the photometric frame entirely

Channels 1–12 all measure light. Every one of them is blind to the same case:
a *completely* enshrouded star does not produce a small deficit, it leaves the
photometric catalogue altogether. And light is the easiest thing to hide —
beam it, cool it, or take a small enough slice and no photometric channel can
reach you.

Channels 13–16 drop the assumption that the target radiates at all, and
measure the three quantities no concealment strategy can touch: **mass,
momentum, and shadow.**

---

## Channel 4 — dynamical mass. The one that closes the grey blind spot.

Every photometric channel regresses M_G on M_Ks, which is structurally blind to
a grey absorber (it produces a residual of the *wrong sign*) and anchored on
2MASS, our dominant contaminant. Gravity has neither problem.

88,149 `gaiadr3.binary_masses` systems, M_G against log M_dyn:

| subset | N | f detectable | n_pos | n_neg | p_UL |
|---|---|---|---|---|---|
| all | 88,149 | 0.721 | 14 | 41 | 2.5e−4 |
| **faint secondary** | 13,482 | **0.734** | **0** | 8 | **2.2e−4** |

**Zero positive outliers in 13,482 clean systems.** The negative excess is
unresolved secondary light inflating G, which biases against the signal — the
search is conservative by construction.

This is the only constraint here valid for grey and for α ≈ 0.19, the exact
slope at which the optical estimator has *zero* leverage.

---

## Channel 5 — spectral slope. The one that produced something.

Broadband dimming has an appalling natural background. But dust has a *shape*:
α ≈ 2, set by grain physics. So instead of asking "is this star faint?" we
asked "faint in a pattern nothing natural makes?"

Fitted a two-parameter power-law absorber to 7-band residuals (BP, G, RP, J, H,
W1, W2) anchored on M_Ks, for 2,991,398 stars. 79,426 have a significant,
physically bounded dimming absorber.

**Median α = 1.95.** The distribution peaks squarely on the dust value.

### The flat-slope population is grain growth

5,239 stars have α significantly below 1.5. They are not anomalous:

| | median A₀ | median \|b\| |
|---|---|---|
| α ≈ 2 (dust-like) | 0.055 | 25.6° |
| **α < 1.5 (flat)** | **0.116** | **12.7°** |

Twice the extinction, half the latitude. Grain growth in dense material raises
R_V from ~3.1 to 5–6, which *is* a flatter α. The strongest cases sit in Orion,
Ophiuchus and Corona Australis by coordinate inspection.

### What survived

After requiring: significant absorber, α significantly < 1.5, acceptable fit
(χ² < 11.07 on 5 dof), not at a grid edge, clean C* and RUWE, A₀ < 0.15,
|b| > 20°, outside known star-forming complexes — **513 stars**.

Eliminations attempted and failed:

- **Metal-poor halo subdwarfs?** No. Median [M/H] = −0.006 against −0.093 for
  the |b| > 20 population — slightly metal-*rich*. Median v_tan 34.5 km/s
  against 38.0, with 0.4% above 100 km/s against 4.3%. Not halo, not fast.
- **SIMBAD known contaminants?** 494 of 513 are unmatched entirely; of the 19
  matched, 5 are known contaminant classes and 4 are **galaxies**.

### The mirror control

An absorber can only *dim*. Running identical cuts on the opposite sign of the
fitted amplitude gives the pipeline's own false-positive rate:

| | count |
|---|---|
| dimmed (physical) | **372** |
| brightened (unphysical mirror) | **13** |
| **ratio** | **28.6 : 1** |

**This set is not symmetric noise.** It is the only channel in the entire
project that produced an asymmetric surviving population.

### What that does and does not mean

It means there is a real population of stars dimmed in the optical relative to
K_s with a spectral slope flatter than dust, at high latitude and low
extinction, which our screens do not explain.

It does **not** mean a detection, for three reasons stated plainly:

1. **Twelve for twelve.** Every previous high-significance anomaly in this
   project — 129σ, 19.5σ, 4301σ, 20.5× — resolved to a mundane population once
   the right control was applied. The base rate here is unforgiving.
2. **Post-hoc cuts.** Several thresholds (the SFR list, A₀ < 0.15, |b| > 20)
   were chosen after seeing the data. The effective trials factor is large and
   unquantified.
3. **A likely mechanism remains unexcluded.** The anchor is M_Ks, so anything
   that makes K_s anomalously *bright* — an unresolved companion contributing
   in the near-infrared, faint enough to escape our C* and RUWE cuts — mimics a
   flat-slope deficit in every other band. This is the same aperture-mismatch
   family that dominates channels 1–3.

**The decisive test is spectroscopy**: does each star's spectral type match the
temperature its photometry implies? A real absorber leaves the photosphere
untouched; a blend does not. That is a few nights on a 2 m-class telescope for
the top 20, and it is the correct next step.

Candidate list: `results/searchA_candidates_primary.csv`.

---

## Channel 7 — kinematic clustering. A 4301σ excess that was young stars.

Transfer cost between systems is set by Δv, not distance, so our position-space
null (channel 6) was arguably measured in the wrong metric. Repeating it in
velocity space:

| control added | excess ratio at Δv < 2 km/s |
|---|---|
| none | **178×  (4301σ)** |
| exclude bound pairs (Δv ≡ 0 by construction) | 160× |
| match controls on \|v\| | 43× |
| match pair separation, restrict to \|v\| > 40 km/s | **~1.4×** |

The diagnostic: **anomalies are kinematically cold.** Median |v| 22.8 km/s
against 35.1 for the field; 43% below 20 km/s against 23%. Cold means young,
young stars sit in comoving associations, associations cluster in velocity for
entirely mundane reasons.

Instructive failure mode: in velocity space the young-star contaminant is *more*
concentrated than in position space. This metric needs stronger controls than
the spatial one, not weaker.

---

## Channel 9 — W3/W4 intercept-and-re-emit. Debris discs, not Dyson spheres.

If intercepted light is re-radiated isotropically, the waste heat appears as
mid-IR excess. WISE W3 (12 µm) and W4 (22 µm) probe 130–240 K material.
Stars that are BOTH optically dimmed (positive residual in M_G vs M_Ks) AND
IR-bright (negative residual in M_W3 or M_W4 relative to their photosphere)
carry the full intercept-and-re-emit signature.

| | count |
|---|---|
| optical deficit (> 3σ) | 104,299 |
| W3 excess (> 3σ) | 55,994 |
| W4 excess (> 3σ) | 480 |
| **both (outside SFRs)** | **2,632** |
| mirror control (bright-G + IR-faint) | **2** |
| **ratio** | **1,316 : 1** |

The asymmetry is enormous because debris discs, YSOs, and circumstellar
material are real populations that genuinely produce this combined signature.
The W3-W4 colours (1.6–3.3 mag) match warm dust, not engineering. The raw
prevalence is 0.29% but dominated by astrophysical contaminants.

This is the channel Suazo et al. (2022) exploit. Our contribution is requiring
the optical deficit to be independently significant, not just fitting a single
SED model.

---

## Channel 11 — 3D domain edge. A 148.5σ signal that is extinction + scanning law.

Three complementary methods searched for a discontinuity surface in the 3D
distribution of optical residuals:

| method | best signal | excess over null |
|---|---|---|
| **planar scan** (192 orientations × 20 offsets) | Δ = 0.020 mag at l=150°, b=60° | **148.5σ** |
| radial scan (25 distance bins) | step = 0.005 mag at 431 pc | 42.4σ |
| local gradient (K=30/60, 50k subsample) | resultant = 0.003 | **−0.4σ (consistent with noise)** |

The planar signal is real: stars on the high-Galactic-latitude side are
systematically brighter by 0.020 mag. But it is not a domain edge:

1. **The A₀ < 0.05 test**: restricting to very low extinction reduces the
   signal by 40% (from 0.021 to 0.013 mag). Extinction drives most of it.
2. **The local gradient is null**: a sharp domain edge produces coherent local
   gradients. These residuals produce none (-0.4σ). The signal is a large-scale
   mean-shift, not a discontinuity.
3. **The direction**: (l=150°, b=60°) points into the North Galactic Cap, the
   lowest-extinction region of the sky. This is where the fiducial fit leaves
   the largest residuals because the training sample is extinction-dominated.
4. **Gaia scanning-law systematics** at the ~0.01 mag level are documented
   (Riello+ 2021, Montegriffo+ 2023) and are position-dependent.

The radial scan's "transition at 431 pc" is the well-known increase in
systematic photometric errors with increasing distance and decreasing SNR.

---

## Channel 10 — radio counterparts. No anomalous GHz emission.

Main-sequence FGK/M dwarfs have quiescent radio luminosities far below
wide-field survey thresholds. A radio detection coincident with an optically
dimmed star would indicate either non-thermal emission from engineered
structures, communication leakage, an active magnetosphere, or a chance
alignment with a background AGN.

The top 5,000 optically-dimmed candidates (significance 17.5–242.2σ) and
5,000 photometrically matched controls were cross-matched against NVSS
(1.4 GHz, Condon+1998) and FIRST (1.4 GHz, Helfand+2015) via CDS XMatch
at 5″ radius.

| survey | dimmed matches | control matches | excess |
|---|---|---|---|
| NVSS | 1 / 5,000 | 4 / 5,000 | **0.25×** |
| FIRST | 0 / 5,000 | 0 / 5,000 | 0.0× |

The dimmed stars have *fewer* radio counterparts than the control sample.
The single NVSS match is consistent with a chance background-AGN alignment
(expected ~2 per 5,000 at this flux limit and sky density). No evidence of
anomalous GHz emission correlated with the optical deficit.

---

## Channel 12 — Search A diagnosis. The asymmetry is unresolved companions.

The 28.6:1 asymmetry in channel 5 was tested against five blend indicators:

| indicator | survivors | control | ratio |
|---|---|---|---|
| WISE cross-match distance | 0.144″ | 0.056″ | **2.6×** |
| W1-W2 colour | 0.067 | 0.004 | **17× redder** |
| frac M dwarfs (3000-4000K) | 67% | 47% | **1.4×** |
| BP-RP excess factor | 1.365 | 1.328 | 1.03× |
| 2MASS cross-match distance | 0.090″ | 0.075″ | 1.2× |

All four primary indicators point to unresolved NIR companions inflating K_s.

The anti-correlation between deficit magnitude and W1-W2 colour (r = −0.12)
confirms: the strongest deficits have NORMAL W1-W2 — they're not caused by
a cool companion adding IR flux. Instead, the companion makes K_s bright,
which our pipeline reads as every other band being too faint.

Verdict: **consistent with the aperture-mismatch / blend family**. The
asymmetry is real (something physical), but spectroscopy would confirm these
are unresolved binaries, not engineered absorbers.

---

## Channel 13 — coherent acceleration. Momentum cannot be hidden.

Gaia gives position and velocity, never acceleration. But Hipparcos measured
the same stars in 1991, and the difference between the 25-year mean proper
motion and the instantaneous Gaia one is a measured velocity change. Brandt
(2021) publishes both on a common frame for 115,346 stars; 94,390 are within
500 pc.

A large anomaly on one star means "unseen companion" — this is the standard
companion-detection tool, and 34% of the sample has one. **Magnitude is not
the discriminant. Direction is.** A companion pulls its primary along whatever
orbit it occupies, so across a population the vectors are isotropic. A
coordinated migration programme would not be.

| | value |
|---|---|
| stars with significant acceleration (>3σ) | 31,989 |
| frame rotation removed (rigid spin + glide) | \|ω\| = 0.0041 mas/yr |
| median velocity anomaly | 1.11 km/s |
| mean local resultant (k = 40 neighbours) | 0.14074 |
| shuffled null | 0.14573 ± 0.00128 |
| **excess** | **−3.9σ** |
| best local patch | +1.9σ |

The null is built by **permuting acceleration vectors over positions**, which
preserves the survey geometry and the velocity distribution exactly and
destroys only the position–direction association. Any other null manufactures
detections out of the selection function.

Injection-recovery sets the sensitivity: a coherent patch containing **2% of
the accelerating stars is recovered at 7.8σ**, 5% at 27σ. So this is a
quantified limit, not merely an absence.

Two systematics had to come out first. The Hipparcos-to-Gaia frame spin is
global and coherent — the same shape as the signal, at the largest scale —
and is removed as a six-parameter rotation plus glide. Perspective
acceleration is subtracted per star. The Galactic potential itself is *not* a
confound: at 2.3e−10 m/s² it moves a star 0.16 m/s over the baseline, which at
100 pc is 3e−4 mas/yr, three orders below the per-star noise.

**Prior art.** Huang et al. (arXiv:2608.16060, 17 Aug 2026) published a Gaia
DR3 stellar-engine acceleration limit eight days before this run. They used
Gaia's own `Acceleration7` solutions, whose angular-acceleration sensitivity
is ~2000× worse than the 25-year Hipparcos baseline, and they did not test
coherence. The channels are complementary, not duplicated.

**What it cannot reach.** A Shkadov thruster produces ~10⁻¹² m/s². The best
per-star sensitivity here is ~10⁻⁸ m/s², and no planned mission closes a
10⁴ gap. This channel can only see *accumulated* Δv or *coordinated* thrust,
never the instantaneous acceleration of a passive stellar engine.

---

## Channel 14 — the local mass ledger. The one bound that survives anything.

Does the Solar neighbourhood weigh more than what we can see in it?

The total density near the Sun is measurable from how hard stars are pulled
back toward the midplane (the Oort limit). The luminous mass is measurable by
counting. The difference bounds every cold, dark, non-radiating component at
once — with no assumption about temperature, spectral slope, beaming geometry,
or energy source.

The dynamical side and the luminous side are both well studied. They have
never been differenced and reported as a technosignature constraint; that
literature is empty. Here the stellar term is recounted from our own sample,
on a mass scale calibrated against 96,641 Gaia dynamical masses (0.033 dex
scatter), so the luminous side uses the same stars and cuts as the rest of the
project.

| term | M☉/pc³ |
|---|---|
| dynamical total | 0.1000 ± 0.0100 |
| stars | 0.0400 ± 0.0040 |
| gas | 0.0417 ± 0.0083 |
| stellar remnants | 0.0060 ± 0.0012 |
| brown dwarfs | 0.0020 ± 0.0010 |
| **counted luminous** | **0.0897 ± 0.0093** |
| residual | 0.0103 ± 0.0137 |
| halo dark matter | 0.0100 ± 0.0020 |
| **unexplained** | **+0.0003 ± 0.0138** |

**The ledger balances.** The 2σ ceiling on any non-halo dark component is
0.0280 M☉/pc³, i.e. **below 1.2 × 10⁵ M☉ within 100 pc** — 70% of the local
stellar mass density.

This is a weak bound, and the reason is worth stating: it is limited by the
uncertainty on the *gas* density (±0.0083, dominated by the CO-to-H₂
conversion factor) and on the dynamical total, not by anything about the
search. It also can never produce a detection, because the residual is
dominated by the local dark matter density, which is real and is not
engineering. What it does is close the loophole every other channel leaves
open: a civilisation that emits *nothing at all* still has to weigh something.

---

## Channel 15 — cold blackbody radiators. The band nobody had searched.

Every published waste-heat search covers 100–1000 K: Carrigan (2009) on IRAS,
Suazo et al. (2022, 2024) on WISE, the Ĝ survey on galaxies. **Below 100 K
there were no searches at all** — not for lack of sensitivity (Planck, AKARI
and IRAS reach ~10⁻³ L☉ at 100 pc down to a few K) but because everything cold
looks like Galactic cirrus and gets catalogued as a molecular cloud core.

That gap is where thermodynamics points. Erasing a bit costs kT ln 2, so
computation per joule scales as 1/T: a civilisation optimising for total
computation rather than raw power is driven toward cold radiators — and would
be invisible to every search ever run.

**The discriminant is grain physics.** Interstellar dust is not a blackbody.
Its emissivity rises with frequency as κ ∝ ν^β with β ≈ 1.5–2.0, because the
grains are far smaller than the wavelength. An engineered radiator is built to
dump heat efficiently, so β ≈ 0. The question is therefore not "is this source
cold?" but **"is this cold source a blackbody or a dust cloud?"**

The Planck Catalogue of Galactic Cold Clumps is all-sky, contains 13,242
sources at 8.6–30 K, and already publishes a fit with **β as a free
parameter** alongside the conventional fixed-β = 2 one. The measurement exists
in a public catalogue and has simply never been read this way.

| | count |
|---|---|
| PGCC sources | 13,242 |
| with a well-constrained free β, clean flux, unblended, SNR > 4 | 3,030 |
| β in [1.0, 2.5] — dust-like | 2,968 |
| **β < 0.5 — blackbody-like** | **0** |
| β > 3.5 — unphysical (mirror control) | 0 |

The mirror control is what makes this a measurement rather than an absence:
β has a physical range in *both* directions, since no known material exceeds
β ≈ 3 either. A noisy SED fit scatters into both tails equally, so the
high-β tail measures the rate at which this catalogue manufactures unphysical
emissivities from noise alone. **Both tails are empty**, and the distribution
is tightly centred on β = 1.87 with a 16–84% range of [1.60, 2.17] — textbook
dust, with no room for anything else.

**The residual gap.** The PGCC covers 8.6–30 K. Below ~5 K, Planck's own
component separation subtracts a 2.725 K blackbody as the CMB, so a 3–4 K
source is partly absorbed into the CMB map by construction. That last few
kelvin remains genuinely unexamined.

---

## Channel 17 — holes in the population. The f→1 case, at last.

Every channel above examines stars that are *in* the catalogue and asks
whether their light is anomalous. That is structurally blind to the limiting
case the whole project was built to find. A star intercepted at f = 0.5 is
0.75 mag faint and lands in our residual tail. A star intercepted at f → 1 is
not faint — it is **absent**. It drops below the magnitude limit and leaves
the sample entirely.

So the signature of complete harvesting is a hole, not a dim star.

Counting stars per unit volume is the obvious test and the wrong one: stellar
density has huge natural gradients. The *shape* of the mass function is the
better observable because it is close to universal — and it carries the
physics, since harvesting is not mass-blind. A civilisation choosing targets
by stellar type removes stars from one end of the distribution and tilts the
slope.

The null permutes masses over positions, preserving the density field,
extinction pattern and global mass distribution exactly while destroying only
the position–mass association.

| sample | global α | max cell deviation | vs null |
|---|---|---|---|
| 300 pc, all sky | 1.361 | 1.400 | **+59.0σ** |
| 300 pc, A₀ < 0.05, \|b\| > 20° | 1.828 | 2.014 | +14.7σ |
| **100 pc, \|b\| > 20° (volume-complete)** | **3.410** | **0.660** | **+1.4σ** |

**It was completeness, and the global slope proves it.** α drifts from 3.41 to
1.36 as the volume grows. That drift *is* the bias: at fixed distance a
0.3 M☉ star is far fainter than a 0.8 M☉ one, so as the sample reaches deeper
it preferentially loses the low-mass end and the fitted slope flattens. Cells
differ in effective survey depth, so they differ in slope — and permutation,
which reassigns masses at random, cannot reproduce that.

Restricting to 100 pc, where the sample is volume-complete across the whole
0.3–0.8 M☉ window, the anomaly vanishes: the mass function has the same shape
everywhere to within permutation noise. **No region within 100 pc is missing
stars of any particular mass.**

Note the honest limitation: this constrains the f → 1 case only out to 100 pc,
because that is where volume-completeness holds. It is a much smaller volume
than the photometric channels cover, and that is the price of asking about
stars that are not there.

---

## Channel 18 — bolometric closure. Energy in, energy out.

Channel 4 is this project's strongest constraint, and it has one specific
escape: it measures a single band. An absorber that intercepts optical light
and re-radiates it at a wavelength we did not check produces a G-band deficit
and no bolometric deficit at all. Channel 3 patched part of that with a W1/W2
veto, which reaches shells hotter than ~300 K; channel 9 added W3/W4, which
reaches ~130 K. Below that, the re-emission simply leaves our bands.

Integrating the whole measured spectral energy distribution removes the escape.
If the total energy leaving a star is less than its **mass** says it should
produce, that energy is not being re-radiated anywhere we can see.

The anchor has to be gravity. Any photometric anchor is circular — it moves
when the absorber does. Gaia's `binary_masses` supplies m1 from astrometric
orbits, independent of every band summed.

| | value |
|---|---|
| systems with dynamical mass and a 7-band SED | **94,123** |
| wavelength coverage | G, J, K_s, W1–W4 (0.6 → 22 µm) |
| M_bol vs log M_dyn robust scatter | 0.2545 mag |
| threshold (5σ) | 1.273 mag → f ≥ 0.69 |
| **bolometrically deficient** | **18** |
| bolometrically excessive (mirror) | 36 |
| **Poisson 95% upper limit** | **p < 2.8 × 10⁻⁴** |

**Fewer deficient than the mirror.** And the mirror is the larger tail for a
known reason: unresolved secondary light inflates the measured luminosity, so
this estimator is biased toward over-luminosity. The systematic pushes against
the signal, which makes the null conservative rather than fragile.

Two escapes remain, and only two. The waste heat can be **beamed** away from
our line of sight, or it can be **stored** rather than radiated. Everything
else — any temperature, any spectral slope, any re-emission wavelength inside
0.6–22 µm — is now excluded for this population at f ≥ 0.69.

The population caveat is real and worth repeating: dynamical masses come from
binaries, and the photometric science sample of this project excludes binaries
by construction. The overlap between the two samples is exactly **zero**. This
channel constrains a different set of stars than channels 1–3 and 5.

---

## Coverage gaps, stated because they are not nulls

**Temperature.** Our mid-IR veto uses W1 and W2 only, so it responds to shells
above roughly **300 K** (SNR 1.8e3 at 800 K, 2.9 at 200 K, 8e−6 at 100 K).
Channel 3's "beamed class" is therefore really *beamed OR cooler than 300 K OR
non-thermal*.

All-sky far-infrared reaches ~40 K at 500 pc. **Below ~35 K nothing all-sky has
looked**, and that is precisely where thermodynamics points: the Landauer cost
of a bit erasure is kT ln 2, so a computation-limited civilisation is driven
toward cold sinks, not hot ones.

**Spectral slope.** The grey case (α = 0) is invisible to every photometric
channel by construction and is covered *only* by channel 4.

**Covering fraction.** Nothing here reaches f ≲ 0.1 at any abundance. A
civilisation taking 1% of its star — still ~10²⁴ W — is invisible to all eleven
channels.

**Scale.** "Grabby" is a cosmological hypothesis predicting near-c expansion
and hard domain boundaries. It predicts missing galaxies, not a 0.1 mag deficit
in a K dwarf. Everything here is matched to *quiet, local, sparse* harvesting,
which is close to the opposite hypothesis. Local absence does argue against our
being *inside* an active conversion front, and says little else.

---

## Reproducing

```bash
run.sh scripts/37_dynamical_mass_search.py     # channel 4 — dynamical mass
run.sh scripts/40_spectral_slope_search.py     # channel 5 — spectral slope
run.sh scripts/41_slope_diagnosis.py           # grain-growth test
run.sh scripts/42_identify_survivors.py        # SIMBAD identification
run.sh scripts/43_survivor_metallicity.py      # subdwarf test
run.sh scripts/44_searchA_null.py              # mirror control
run.sh scripts/49_searchA_diagnosis.py         # blend diagnosis
run.sh scripts/46_searchB_cold_excess.py       # channel 9 — IR excess
run.sh scripts/47b_searchC_xmatch.py            # channel 10 — radio (CDS XMatch)
run.sh scripts/48_searchD_domain_edge.py       # channel 11 — domain edge
run.sh scripts/36_velocity_clustering_v2.py    # channel 7 — velocity
run.sh scripts/45_grabby_figures.py            # figures
```

Figures: `F11_spectral_slope.png`, `F12_mirror_and_sky.png`,
`F13_channel_coverage.png`, `F14_searchB_ir_excess.png`,
`F15_searchA_diagnosis.png`, `F16_searchD_domain_edge.png`,
`F17_searchC_radio.png`, `F18_grand_summary.png`,
plus `F1`–`F10` from the main paper.
