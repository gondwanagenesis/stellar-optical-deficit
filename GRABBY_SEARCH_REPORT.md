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
| 2 | optical deficit, wide-pair differential | same, cleaner background | **p < 1.4e−3** (was 5.3e−4; efficiency term) |
| 3 | mid-IR veto (beamed class) | deficit + no warm re-emission (>300 K) | **p < 8.6e−4** (was 4.3e−4; efficiency term) |
| 4 | **dynamical mass** | **ANY spectral slope, incl. grey** | **p < 2.2e−4**, 0 positives |
| 5 | spectral-slope fit, 7 bands | α ≠ 2 absorbers | 372 asymmetric survivors |
| 6 | 3D spatial front | a spreading domain | 19.5σ → dust |
| 7 | kinematic (Δv) clustering | a propagating network | 4301σ → young associations |
| 8 | deficit colour type | grey vs reddening | all reddening-like |
| 9 | W3/W4 intercept + re-emit | warm re-emission (100–300 K) | 2,632 candidates → debris discs |
| 10 | radio counterparts (NVSS/FIRST) | non-thermal emitters | 1/5000 dimmed, 4/5000 control → null |
| 11 | 3D domain edge | grabby boundary | 148.5σ → extinction + scanning law |
| 12 | Search A diagnosis | NIR companion test | **consistent with blends** |
| 13 | coherent acceleration (Hipparcos–Gaia) | coordinated stellar migration | isotropic, **z = −0.03** (null corrected); ≥2% domain excluded |
| 14 | local mass ledger (Oort limit) | **any** cold mass, any technology | < 1.2e5 M☉ within 100 pc |
| 15 | cold blackbody, 3–100 K (β free) | engineered radiators below every prior search | 0 of 3,030; all dust-like |
| 16 | dark companions vs WD channel | fully enshrouded star (f→1) | 486:0 → extrapolated orbits |
| 17 | stellar mass function, spatially | **holes** in the population | 59σ → completeness; null at 100 pc |
| 18 | **bolometric closure vs dynamical mass** | **energy missing at ANY wavelength** | **18 vs 36 mirror → p < 2.8e−4** |
| 19 | **Planck excluded compact-source bin** | cold radiators at 100–857 GHz, invisible to every IR survey | 36 vs 325 mirror → **null**; two grid-edge bugs caught |
| 20 | **high-RUWE mass without light** | unseen mass with no photometric trace (f→1), a population every other channel excluded server-side | **null**; excess shrinks with cut depth; two comparison bugs caught |
| 21 | **energy-conservation locus (ΔL_opt vs ΔL_IR)** | a partial absorber that re-emits what it removes, at ANY covering fraction | **null**; 146σ excess is a broad dim/IR correlation, balanced band is its trough |
| 22 | **sub-PSF orientation vs proper-motion axis** | comoving vs background blends, the Hephaistos false-positive class | **null**; alignment decays 58% as the null cell shrinks |

**Joint, disposal-agnostic: p_total < 1.05e−3.** Fewer than 1 in 953 nearby
lower-main-sequence stars intercepts ≥51% of its optical output by any means.
*Corrected 2026-08-21.* The previously quoted p_total < 6.2e−4, 1 in 1,614, omitted
the detection-efficiency term and was a factor of two too strong at that covering
fraction; it survives unchanged at f ≥ 0.58. See **Audit** below.

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
| mark-permutation null (**invalid**, see below) | 0.14574 ± 0.00123 → z = −4.06 |
| line-of-sight-spin null (**correct**) | 0.14079 ± 0.00152 → **z = −0.03** |
| best local patch | +1.9σ |

Injection-recovery sets the sensitivity against the corrected null: a coherent
patch containing **2% of the accelerating stars is recovered at 9.6σ**, 5% at
26σ, 1% at 4.9σ. So this is a quantified limit, not merely an absence.

### The null was wrong, and the −3.9σ was the symptom

This channel originally reported −3.9σ — the observed statistic sitting
*below* its own permutation null — and called it a null result. It was a bug,
caught by re-auditing the channel (scripts `76_` and `77_`).

The first hypothesis was physical: resolved wide binaries, each component
pulled toward the other, would anti-align. Script `76_` tested it and it
failed. Only six pairs lie within 0.02 pc, and removing them made the deficit
*worse* (−3.91 → −4.26σ). The nearest separation bin showed no significant
anti-alignment at all (mean cos = −0.071, z = −1.2).

The cause is geometric. **The proper-motion anomaly is a transverse quantity**
— it lives in the plane perpendicular to the line of sight, always. Stars that
are neighbours in 3D share very nearly the same line of sight, so their anomaly
vectors are confined to nearly the same plane. Mark permutation destroys
exactly that: it hands a star a vector computed for a *different* line of
sight, which is not transverse to the recipient's. Permuted neighbourhoods
therefore span the full sphere while real ones span a plane, and for *k* random
unit vectors the expected resultant is larger in three dimensions than in two:

| | E\|Σu\| / √k |
|---|---|
| 3D (permuted) | √(8/3π) = 0.921 |
| 2D (real) | √π/2 = 0.886 |

The null was inflated by construction, by 3.9%. Measured: 0.14574 / 0.14079 =
**1.035**, against the predicted 1.039. The channel was being compared against
a baseline no real data could ever reach.

The diagnostic is direct. Mean |û · LOS| is 0.00000 for the observed vectors
and 0.50065 for permuted ones — exactly the 1/2 expected for random directions
in 3D.

**The corrected null** randomises only the quantity the signal lives in:
rotate each star's anomaly vector by a random angle *about its own line of
sight*. That preserves transversality, sky positions, the density field and the
magnitude distribution, and destroys only the position-angle correlation
between neighbours — which is precisely what a coordinated thrust produces.
Against that null, z = −0.03. Still null, now calibrated.

**Standards amendment.** Mark permutation is the project's default null and is
correct for scalar marks. It is *not* automatically correct for vector marks:
if the vector is constrained by the geometry of the observation, permutation
breaks the constraint and the null stops being the same kind of object as the
data. The general rule is to randomise *within the constraint manifold* — here,
the transverse plane.

Auditing the other channels against this: channels 6 and 11 shuffle scalar
residuals and are unaffected. Channel 13 was the only permutation-null channel
carrying a vector mark. **Channel 7 is a follow-up item**, for a related but
distinct reason — its marks are tangential velocities, also 2D and also
line-of-sight-constrained, so genuinely neighbouring stars have a
systematically *smaller* Δv than randomly paired ones purely from shared
geometry, which pushes in the direction of manufacturing a small-Δv excess.
Its null is control-matched rather than permutation-based and its
separation-matched control should absorb most of that, but the size of the
residual effect has not been measured.

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

## Channel 16 — dark companions. A 486:0 asymmetry, and a broken control.

Gravity measures a companion's mass whether or not it emits anything, so a
fully enshrouded star in a binary is a *maximal* astrometric signal rather
than a null one. From 72,152 Gaia DR3 astrometric orbits we computed the
astrometric mass-ratio function from the Thiele-Innes coefficients, solved for
the minimum companion mass assuming zero companion light, and asked which
companions violate the white-dwarf formation channel.

The physics is sharp. Below ~0.50 M☉ no white dwarf forms from a single star
in a Hubble time, so anything lighter is a *helium* white dwarf and must have
been stripped by a companion. Stripping welds mass to orbit: the donor's core
mass sets both the remnant mass and the donor's radius at the end of mass
transfer, giving P_orb = 1.2×10⁵ (M_WD − 0.12)^4.5 days (Tauris & Savonije
1999), a relation that terminates near 0.47 M☉ where helium ignites. A
0.25 M☉ helium white dwarf belongs at ~10 days. And stable Roche-lobe overflow
circularises, so wide orbits should have e ~ 10⁻³.

After the full quality cascade — Thiele-Innes relative error, the DR3 pipeline
validity cuts, magnitude-conditioned goodness-of-fit, scanning-law period
rejection, within 300 pc — **17,717 clean orbits** remained, of which 8,831
had a companion in the 0.20–0.45 M☉ helium window.

**486 sat in the anomaly box (P > 1500 d, e > 0.15) against 0 in the circular
control.** Median violation of the mass-period relation: a factor of 23.

### Why that is not a detection

The zero is what gives it away. If fitted eccentricity carried real scatter,
some long-period systems would land at low e by chance. Finding *none* meant
either a real population or a broken control — and it was the control.

**88% of every astrometric orbit in this sample has e > 0.15**, including
inside the baseline. Gaia DR3 astrometric solutions are overwhelmingly
eccentric, so "e < 0.05" selected essentially nothing and the 486:0 ratio
carried no information. The control was designed on an assumption about the
eccentricity distribution that the data does not satisfy.

Three tests then settled it:

| test | result |
|---|---|
| frac(e > 0.15) within the 1000 d baseline | 0.877 |
| frac(e > 0.15) beyond 1500 d | 0.987 |
| **spread in frac(e > 0.15) across companion mass 0.2 → 1.4 M☉** | **0.083** |

The mass-independence test is decisive. The physical claim was specific to
*low-mass* companions, because only they require stripping. But wide eccentric
systems are just as common at 0.9–1.4 M☉ (0.909) as at 0.20–0.45 M☉ (0.992).
The selection was never testing the helium-white-dwarf relation — it was
selecting long-period Gaia solutions, which are eccentric regardless of what
orbits them.

And the solution quality is damning:

| | anomaly box | rest | ratio |
|---|---|---|---|
| relative period error | 0.126 | 0.004 | **31×** |
| period error (days) | 240 | 1.7 | 138× |
| significance (a₀/σ) | 10.4 | 36.1 | 0.29× |
| RUWE | 9.7 | 4.5 | 2.2× |

Every box member has a period beyond the ~1000 d DR3 baseline, so every one is
extrapolated — and they are precisely the badly-constrained tail.

**This is the structural problem stated when the channel was designed**: the
physically clean region begins at P > 1500 d, and DR3's baseline ends at
~1000 d. The region where this test has power is exactly the region where the
orbits cannot be trusted. That tension is not resolvable with DR3. Gaia DR4
(66-month baseline, epoch astrometry published for the first time) is the
release that fixes it, and this search should be re-run then.

Prior art also constrains the claim: Shahaf et al. (2019, 2023, 2024) already
classified 101,380 DR3 orbits into main-sequence, triple and compact-object
companions and published ~3,145 white-dwarf candidates. "There is a dark
companion here" is their result. Any contribution has to be a discriminant
applied on top, and the mass-period one attempted here needs DR4 to work.

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

## Channel 19 — Planck's excluded bin. Null, after two grid-edge bugs.

Hiding by running cold forces you to run big, but the band matters more than
the size. At 3 K and 100 µm, hν/kT = 48 and the Planck function is suppressed
by e⁻⁴⁸ ≈ 10⁻²¹. **IRAS, AKARI, WISE and Herschel are all identically blind to
a 3–5 K source.** The Wien peak of a 3 K blackbody is 176 GHz, inside Planck's
CMB channels — and every technosignature search ever published used the
infrared. This also closes a structural hole in channel 15: the Planck
cold-clump catalogue requires IRAS 100 µm as an input band, so a genuinely
3–5 K source is unselectable by it at any brightness.

An opaque 3 K surface against the 2.725 K CMB is a 275 mK excess on ~100 µK rms
anisotropy: a ~3000σ compact positive spot. Component separation does not
remove it — an internal linear combination preserves anything with a blackbody
spectral response. What happens to it is *masking*, and the place it lands is
**PCCS2E**, the bin of Planck compact sources that failed validation. Planck's
documentation treats it as contamination. It has never been searched as a
source population.

120,491 excluded detections, grouped by position into 20,932 sources with three
or more of the six HFI bands, each fitted for *I* ∝ ν^β *B*ν(*T*) with β and
*T* free. Dust is β ≈ 1.6; a cold blackbody is β = 0.

| | cold tail | mirror (symmetric about dust) | ratio |
|---|---|---|---|
| script `74_`, β ∈ [−0.5, 3.0] | 1,023 | **0** (control dead) | — |
| script `78_`, β ∈ [−0.5, 5.0] | 1,023 | 233 | 4.39 |
| script `80_`, β ∈ [−4.0, 7.2], interior fits only | **36** | **325** | **0.11** |

**Null.** The cold tail sits nine times *below* its own measured
false-positive rate. Getting there took two separate grid-edge bugs, and both
are worth recording because the failure mode is general.

### Bug 1 — the mirror control could not fire

Script `74_` fitted β on `np.arange(-0.5, 3.01, 0.05)`, whose largest element
is 2.9999999999999991, and then defined its unphysical mirror as β > 3.0.
**That condition can never be satisfied.** The mirror was identically zero for
any input whatsoever — and the verdict logic compares the candidate count
against it, so a candidate list of any size would have been reported against a
false-positive rate measured as exactly zero.

The rebuild (`78_`) extended the grid to β = 5, made the mirror **symmetric
about dust** rather than pinned to a round number (dust at 1.6, cold blackbody
at 0, so the equal-and-opposite absurdity is 3.2, and both tails cut at
|β − 1.6| ≥ 1.1 with the same T < 10 K and Δχ² > 4), and vectorised the fit
over the (T, β) grid — verified bit-identical to the original scalar fitter on
200 real sources before the grid was touched.

### Bug 2 — extending one edge does not validate the other

`78_` then returned 1,023 cold against 233 mirror, a 4.4× asymmetry, and the
fitted temperatures clustered at 2.75 K against T_CMB = 2.72548 K with a 56×
plane-to-pole rate gradient. That was written up as CMB confusion. **It was
wrong.** 88.6% of those candidates sat at β = −0.50 *exactly* — the lower grid
edge, inherited from `74_` and never questioned because the bug being hunted
was at the other end. A fit pinned to a boundary is not a measurement; it is
the model reporting that it cannot describe the data. The temperature
clustering the CMB argument rested on was a nuisance parameter absorbing a
slope β was not allowed to reach.

**β below −0.5 is not exotic.** In this parameterisation the flux density goes
as ν^(2+β) in the Rayleigh–Jeans limit:

| β | S ∝ | what it is |
|---|---|---|
| +1.6 | ν³·⁶ | interstellar dust |
| 0.0 | ν² | cold blackbody — the target |
| −2.0 | ν⁰ | **flat spectrum** |
| −2.7 | ν⁻⁰·⁷ | optically-thin synchrotron |

Widening the floor to β = −4 lets the pinned sources find their real slope, and
**1,071 of them land below −0.5 with a median β of −2.15**, 40% within 0.3 of
β = −2 exactly. Their |b| < 10° fraction is 0.16 against 0.62 for the
catalogue, and their rate rises 42× from plane to pole. That is the
extragalactic flat-spectrum radio population — blazars and FSRQs — isotropic on
the sky and suppressed in the plane by Galactic confusion. It reproduces the
entire gradient that had been attributed to the CMB.

On the widened grid, boundary occupancy is 0.01% at the β floor and 0.03% at
the ceiling, so the fit is no longer edge-limited where it matters. (36% rail
at the T = 39.5 K ceiling — ordinary warm dust wanting a hotter fit — which
touches neither tail, since both require T < 10 K.)

**The reusable lesson: check boundary occupancy at every edge of every fit, as
routine, and treat a railed fit as missing data.** Extending a grid at one end
proves nothing about the other, and a boundary-pinned population will
manufacture whatever correlation its nuisance parameters need.

---

## Channel 20 — mass without light. The population that was never downloaded.

Every one of the nineteen channels above ran on a sample whose ADQL carried a
server-side RUWE cut; the local partitions top out at `ruwe = 1.387`. The
astrometrically-disturbed stars were never on disk, so they were never
examined, and they are excluded from every limit this project has published.

That is the wrong population to have discarded for this particular question.
RUWE measures how badly a single-star astrometric model fits, and it rises
precisely when an unseen mass pulls a star around. A fully enshrouded
companion is an unseen mass. The f → 1 case that every photometric channel is
structurally blind to would announce itself in RUWE and essentially nowhere
else — and RUWE is what the sample build filtered on.

**224,081 sources** with `ruwe ≥ 1.4`, `parallax > 2`, `parallax_over_error >
20` were pulled across a uniform 20% all-sky sample of the identifier space;
222,757 (99.4%) carry 2MASS Ks. 149,536 fall in the same lower-main-sequence
box as the main sample, of which **92,976 (62.2%) are photometrically
pristine**.

### The discriminant, and what it says

High RUWE on its own means "unresolved binary" and is unremarkable. The
question is whether the disturbance is accompanied by *light*. An ordinary
companion contributes flux: it inflates the BP/RP excess factor, displaces the
2MASS cross-match, and lifts the system off the main sequence. An enshrouded
companion carries the same mass and emits nothing, breaking that relation.

The relation holds, cleanly and monotonically:

| RUWE | n | fraction pristine | median C\* | median residual |
|---|---|---|---|---|
| 1.4–2.0 | 63,045 | 0.674 | +0.87 | +0.191 |
| 2.0–3.0 | 38,189 | 0.613 | +1.05 | +0.212 |
| 3.0–5.0 | 27,533 | 0.595 | +1.24 | +0.226 |
| 5.0–10.0 | 16,454 | 0.545 | +1.87 | +0.256 |
| > 10 | 4,315 | 0.404 | +4.26 | +0.344 |

Astrometric disturbance and photometric blend evidence rise together across a
factor of seven in RUWE, which is what unresolved binaries do and what mass
without light cannot do.

### Two bugs in the comparison, both pushing the same way

The first run of this channel reported the pristine sample as *less*
one-sidedly dim than the low-RUWE reference — 3.96 against 22.60, p = 1 — and
called it a null. **That comparison was invalid**, for two independent reasons:

1. **One threshold, two widths.** Both populations were cut at a fixed
   0.530 mag. Against each population's own robust scatter that is **5.35 σ for
   the reference but only 1.81 σ for the pristine high-RUWE sample**. A
   dim/bright ratio always climbs as the cut moves out, so the two numbers were
   read off entirely different points of their respective tails.

2. **Two different estimators.** The high-RUWE residuals came from a degree-5
   polynomial refit inside the script (σ = 0.1767); the reference residuals
   were read from the pipeline's `residual` column, built with a local-slope
   model (σ = 0.0991). The refit exists *specifically* so the comparison would
   be like-for-like, and then the comparison was made against the other
   estimator anyway.

Both errors widen the high-RUWE distribution relative to the reference and so
suppress its apparent one-sidedness. The null was manufactured.

### The corrected channel

One polynomial applied to both populations, one selection box, one distance
support (69–500 pc, the reference's own span), counted at *k* σ of each
population's own scatter:

| k | population | σ | threshold | dim | bright | ratio |
|---|---|---|---|---|---|---|
| 2.0 | reference, low RUWE | 0.1767 | 0.353 | 182,537 | 68,104 | 2.68 |
| 2.0 | high-RUWE pristine | 0.2921 | 0.584 | 4,873 | 1,081 | **4.51** |
| 2.0 | high-RUWE blended | 0.3724 | 0.745 | 5,625 | 517 | 10.88 |
| 3.0 | reference, low RUWE | 0.1767 | 0.530 | 74,676 | 12,363 | 6.04 |
| 3.0 | high-RUWE pristine | 0.2921 | 0.876 | 976 | 214 | **4.56** |
| 3.0 | high-RUWE blended | 0.3724 | 1.117 | 2,426 | 107 | 22.67 |

The pristine sample is **1.68× the reference ratio at k = 2** (p = 6.3e−60) but
**0.76× at k = 3**. The excess *shrinks* as the cut moves out.

That direction is the whole result. A genuine population of dim objects sits at
some offset, so isolating it harder makes it more conspicuous and its excess
must **grow** with k. Extra symmetric scatter does the opposite: it inflates
both tails and drags the ratio back toward the reference's shape. The 2 σ
excess is therefore the degraded astrometry of high-RUWE stars — noisier
parallaxes, hence noisier distance moduli — and not a one-sided dim tail. This
is a discriminant a single fixed cut could not have seen at all.

The pristine cut also demonstrably works, which the first run could not show:
the blended subsample runs at 10.88 and 22.67 against 4.51 and 4.56 pristine, a
factor of five of contaminant suppression at k = 3.

**No numerical limit on f is quoted from this channel.** The dominant
contaminant here is one-signed — 2MASS aperture mismatch in a 4″ beam brightens
Ks and so manufactures a *dim* residual — so the bright tail is a conservative
mirror rather than a clean false-positive rate, and one-sidedness on its own
proves nothing.

### Infrastructure note

The original puller for this population could not work: it used the keyset
pagination (`source_id > X ORDER BY source_id`) that is correct on the 72k-row
NSS table but catastrophic against 1.8e9-row `gaia_source`, where an open-ended
lower bound leaves the planner no bound on the scan. Measured: **>300 s for a
single page** against **16.8 s for an equivalent indexed `BETWEEN` range
scan**. One run burned 2 h 47 m without completing a page.

It also cross-matched Gaia epoch-2016.0 positions against a ~1999.5 survey
inside 3″ without propagating proper motion. Measured on the delivered sample,
the 2016→1999.5 displacement has **median 0.37″, 99th percentile 3.02″, maximum
51.5″** — so over 1% of the sample falls outside the match radius outright. The
subtler damage is downstream: the pristine cut uses `tmass_xm_dist < 0.5″` as a
*blending* indicator, and at a median displacement of 0.37″ an unpropagated
separation would have been substantially a proper-motion cut, selecting on
kinematics at the exact point where this channel's discriminant lives. This is
the same class of error as the 5″ match that produced a fake null earlier in
this project, which is twice now.

---

## Channel 21 — energy conservation as a locus. A 146σ excess that is a debris disc.

Channel 9 required a star to be *both* optically dim and infrared bright and
found 2,632 candidates, all of which looked like debris discs. It never checked
the thing that would have separated a disc from an intercepting structure: that
the two effects **balance**. A partial absorber conserves bolometric
luminosity. It moves flux out of the optical and re-emits it in the infrared,
so at fixed parallax and effective temperature the two residuals must lie on a
locus of slope −1 in luminosity units — the energy taken out equals the energy
put back. Searching either side alone throws away half the signal and all of
the specificity.

The discovery statistic is therefore the **balance ratio**

  B = ΔL_IR(re-emitted) / ΔL_opt(removed),

counted in a band around B = 1, against a mark-permutation null.

### The first attempt was arithmetically broken

`scripts/75_searchU_energy_balance.py` formed B from νf_ν(M_G) over
νf_ν(W3_apparent) — an **absolute** optical magnitude divided into an
**apparent** infrared one. Every ratio was therefore the physical ratio divided
by (d/10 pc)², a median suppression of **3.16 dex** across candidates whose
median distance is 379 pc. All 2,632 objects railed into the "starved" bin,
`n_balanced` and `n_mirror_balanced` both went to zero, the mirror control could
not fire because both arms were pinned against the same edge, and the
zero-count branch printed a claim of a two-population split that the same
file's own counts contradict (`n_runaway_disclike = 0`). Correcting only the
magnitude system moves the median log₁₀B from −4.28 to −1.17 and puts 314
candidates inside the band. That JSON is now stamped `SUPERSEDED` on disk.

### The rebuilt channel

`scripts/88_searchW_energy_locus.py` rebuilds it on the full sample rather than
on channel 9's already-selected candidate list — selecting on both residuals
first and then asking whether they balance is circular — using channel 9's own
photospheric fiducial, 234,377 usable stars out of a 3,321,566 parent, a local
mark-permutation null over 511 cells in (distance, M_Ks, A₀), 200 permutations.

Taken alone it is a **146σ detection**: 30,881 objects in the balanced band
against a null mean of 16,828 ± 96. It survives every crude control. Extinction
quartiles all fire (56–94σ). Restricting to A₀ < 0.3 leaves 103σ. The excess
*grows* with deficit depth, exactly as a real population should — observed over
null runs 2.22, 2.74, 2.98, 3.02, 3.30 as the optical cut moves out to 3σ.
Injection recovers 68–87% of an injected balanced population, so the estimator
works.

### The profile kills it

The test that matters is *where in B* the excess lives. Conservation predicts a
maximum at B = 1. What the data show is a **U**:

| log₁₀B | −4.0 | −3.0 | −2.0 | **−0.9** | 0.0 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| obs/null | 2.30 | 2.22 | 1.84 | **1.66** | 1.93 | 2.18 | 2.57 |

The excess is present in **every bin from log₁₀B = −4 to +2**, and the balanced
band averages 1.91 — at the *bottom* of the profile. Energy conservation
predicts a peak at unity; the balanced band is the least enhanced slice of the
plane. What is being measured is a broad correlation between optical dimness
and infrared brightness, which any reddening-plus-disc population produces, and
the shuffle destroys that correlation wherever it sits.

The balanced objects are also not balanced systems. Their median infrared
residual is **−0.614 mag** against −0.158 for the rest of the sample, while
their median optical deficit is *smaller* — +0.023 against +0.066. They are
large infrared excesses divided by negligible optical losses, landing near unity
by arithmetic accident. That is the signature of a disc that intercepts a
geometrically small fraction of the starlight and reprocesses it into a large
mid-infrared excess, not of a structure that removes what it re-emits.

The mirror control — optically bright and infrared faint, a population that
cannot host an absorber because an absorber can only dim — fires at 76σ with an
excess rate of 1.87% against the signal's 6.00%. Signal over mirror is 3.2,
which would have looked like a detection had the profile not already explained
it.

**Channel 21 is null.** The 146σ excess is real and it is not conservation.

## Channel 22 — sub-PSF orientation vs the proper-motion axis. Null at every resolution.

`ipd_gof_harmonic_phase` is the position angle, mod 180°, of whatever
unresolved structure Gaia's image-parameter fit cannot describe with a single
PSF. That gives a discriminant nothing published uses. A **stationary
background blend** has a separation vector from the star that rotates at the
proper-motion rate, so its harmonic phase tracks — or drifts toward —
PA_PM = atan2(pmra, pmdec). A **comoving** structure keeps a fixed geometry, so
its phase is uncorrelated with PA_PM. This is exactly the false-positive class
that destroyed the Project Hephaistos Dyson-sphere candidates.

Fresh pull of `ipd_gof_harmonic_amplitude` and `ipd_gof_harmonic_phase`
(`scripts/84_pull_ipd_harmonic.py`; these are not in the local partitions),
224,081 usable sources, statistic = the Rayleigh resultant R of the doubled
angle 2(φ − PA_PM).

**The mirror control is degenerate here and is not quoted.** R of a doubled
angle is invariant under a constant rotation of Δφ, so testing against
PA_PM + 90° returns the identical R by construction. The false-positive rate has
to come from the local mark permutation instead, which is why the resolution
scan below is the whole channel rather than a footnote to it.

The first run reported ALIGNMENT DETECTED: R = 0.0409 observed against a local
nside = 8 shuffle null of 0.0375 ± 0.0013, p = 0.002. Note already how little of
the observed R is signal — the null alone sits at 92% of it. That is survey
geometry: the Gaia scanning law makes both φ and PA_PM smooth functions of
position on the sky, and two smooth fields correlate inside any cell large
enough to contain their gradients.

`scripts/87_searchV_resolution_scan.py` shrinks the null cell and watches both
the excess and the null's own power:

| nside | cell | excess R | σ | recovers injected 0.5% |
|---|---|---|---|---|
| 2 | 29.3° | +0.01185 | 8.5 | yes |
| 4 | 14.7° | +0.00589 | 4.8 | yes |
| 8 | 7.3° | +0.00349 | 2.7 | yes |
| 16 | 3.7° | +0.00216 | 1.8 | yes |
| 32 | 1.8° | +0.00146 | 1.2 | yes |
| 64 | 0.9° | +0.00142 | 1.2 | yes |
| 128 | 0.46° | +0.00148 | 1.5 | yes |

The excess falls by **58% from nside 8 to nside 128** and by 88% from nside 2,
while the null keeps recovering an injected aligned population of 0.5% at every
resolution — the shrunken null has not lost its power, the signal has lost its
support. An alignment carried by real per-star structure cannot care how the
sky is diced. One carried by two smooth fields must decay exactly like this.

**Channel 22 is null.** The first run's verdict string is stamped
`VERDICT_SUPERSEDED` on disk; its suggestion that a "comoving non-aligned
remainder" had been isolated must not be cited, because no aligned component was
ever measured and so no remainder was defined. The RUWE control is flat
(−0.71σ across quartiles, Spearman ρ = −0.4), and a group-loop versus
vectorised cross-check of the estimator agrees to within 3 standard errors.

This channel does not constrain harvesting. It establishes that Gaia's sub-PSF
orientation carries **no measurable blend signature at the population level**,
which means the blend-versus-comoving discrimination that channel 20 had to
make indirectly cannot be made directly this way either.

---

## Audit — the headline limit was a factor of two too strong

Every channel above reports a null. The wide-pair limit is the only place this
project quotes a **number**, so it is the only thing anyone would cite, and it
had never been audited. It is also the number that was rebuilt twice under
pressure (scripts 32 → 33 → 34) after the original asymmetry estimator turned
out to be blind to a symmetric signal. `scripts/89_audit_pair_limit.py` rebuilds
it through the same code path and checks it.

**It reproduces exactly.** Every cell of `pair_limit_v3_primary.json` —
thresholds, counts, f_det, p_UL — comes back bit-for-bit from a fresh fiducial
fit, pair search and mid-IR veto. The pipeline is not the problem. The formula
is.

### The missing efficiency term

`pipeline.statistics.exclusion_curve` — the shared helper channels 1, 3 and the
exclusion contour all use — defines the limit as

    p_UL(f) = N_UL / (N_total × efficiency(f, k))

`scripts/34_pair_limit_v3.py` hand-rolled its own table and wrote
`p_UL = best / n_stars`. There is no efficiency factor anywhere in it.

The consequence is exactly the factor the omission implies. f_det is defined as
the covering fraction whose magnitude deficit *equals* the threshold
T = k·σ_dr. A star sitting exactly at threshold is scattered above it half the
time. Measured by shifting the observed Δr distribution — so the real
non-Gaussian noise shape is used, not a Gaussian assumption — the efficiency at
the quoted f_det is **0.5001**, at every k, in both samples.

So the headline

> p_total < 6.2e−4, fewer than 1 in 1,614 stars intercepts ≥ 51%

is a factor of two too strong **at the covering fraction it names**. Corrected:

| | p_dark | p_total | one star in | at f ≥ |
|---|---|---|---|---|
| as published | 4.30e−4 | 6.20e−4 | 1,614 | 0.507 |
| **corrected, same f** | **8.59e−4** | **1.05e−3** | **953** | **0.507** |
| published number, corrected f | 4.30e−4 | 6.20e−4 | 1,614 | **0.580** |

Both bottom rows are true; they are the same result quoted two ways. The
published *number* is not wrong so much as attached to the wrong f — it holds
where the channel is 90% efficient (f ≥ 0.58), not at its 50% threshold. At
99% efficiency it holds at f ≥ 0.66.

This is the same failure mode as the hand-reimplemented C\* formula caught
earlier in this project: a shared, tested helper existed and was bypassed by a
local reimplementation that dropped one of its terms. That is now twice.

### The background model is for the wrong statistic

At every k the observed count runs far below the prediction — 3 against 37.7 at
the headline row, 206 against 422 at k = 4. The report treated that as
reassurance. A background model that over-predicts by 13× is a broken model,
and it is load-bearing, because `best = min(ul_cons, ul_sub)` takes the smaller
of the Poisson and the background-subtracted limit.

Two explanations were tested and **both failed**, which is the useful part.

1. *A threshold-units mismatch.* T = k·σ_dr, but q's own scale is σ_r
   (0.11338 against σ_dr = 0.10957), so q is evaluated at 6.76 σ_r rather than
   7. On a Gaussian tail that is a factor ~5, about the size of the
   discrepancy. Correcting it moves obs/pred from 0.080 to **0.094** and no
   further.
2. *A mark-permutation null on Δr* — residuals shuffled within
   (M_Ks, distance) cells over the same pairs, 200 permutations. It
   over-predicts **worse**: obs/pred = 0.062.

What is left is physical. q(T) = P(|r_star| > T) is the tail of a *single-star*
statistic, applied to a *difference*. The single-star tail is dominated by
contaminants common to both components of a wide pair — shared crowding, shared
extinction error, shared fiducial mis-specification — and those cancel in Δr.
The core already carries **53% common-mode variance**, and the tail evidently
carries far more. The permutation over-predicts most precisely because it
destroys the cancellation outright.

So this analysis has **no valid background for |Δr|**, and the subtracted branch
must not be taken at any k. It is taken in **6 of the 8 published rows**, which
weaken by factors of 1.32 to 5.61 once forced onto the Poisson branch. The
headline row is not one of them — clean + bare at k = 7 already takes the
conservative branch — so this finding costs the abstract nothing and costs the
rest of the table a great deal.

### Two smaller notes

**The sign split is negative-heavy everywhere** (4:6, 10:16, 90:116, …). Δr is
primary-minus-secondary and a harvested component can be either, so a signal is
symmetric in sign and this split is a noise diagnostic rather than a mirror.
The negative excess is the unresolved-companion brightening already identified
in channel 12. Those objects cannot host an absorber, so counting them in a
two-sided Poisson limit makes it conservative by roughly the negative fraction.

**The quoted row is the best of eight** — four k × two samples — selected by
minimum mean-f limit with no trials penalty. At a fixed pre-registered k = 5 the
same sample gives p_dark < 4.00e−3 at f ≥ 0.396 after efficiency correction.
The optimism is real but modest next to the factor of two above.

### What replaces the abstract line

**p_total < 1.05e−3 — fewer than 1 in 953 nearby lower-main-sequence stars
intercepts ≥ 51% of its optical output by any means.** The stronger 1-in-1,614
statement survives at f ≥ 0.58.

---

## Coverage gaps, stated because they are not nulls

**Temperature.** Our mid-IR veto uses W1 and W2 only, so it responds to shells
above roughly **300 K** (SNR 1.8e3 at 800 K, 2.9 at 200 K, 8e−6 at 100 K).
Channel 3's "beamed class" is therefore really *beamed OR cooler than 300 K OR
non-thermal*.

All-sky far-infrared reaches ~40 K at 500 pc. **Below ~35 K nothing all-sky had
looked**, and that is precisely where thermodynamics points: the Landauer cost
of a bit erasure is kT ln 2, so a computation-limited civilisation is driven
toward cold sinks, not hot ones. Channel 19 closes part of that band from the
millimetre side and reports a null for T ≥ 4 K.

**The CMB floor, T → 2.725 K.** Channel 19's null covers 3–10 K, but the limit
weakens to nothing as T approaches the CMB temperature, and that is physics
rather than a pipeline shortcoming. A surface in thermal equilibrium with the
CMB radiates no *net* excess and is invisible in principle, not merely below
threshold; approaching it, the contrast against a 2.725 K background that is
itself fluctuating at ~100 µK falls away, and Planck's compact-source
extraction masks exactly such objects. 2.725 K is also the coldest any passive
radiator in the present universe can be, so this is the floor of the entire
temperature axis — the one place a waste-heat argument genuinely cannot follow.
Reaching into it needs angular profile, ACT/SPT-resolution follow-up, or the
non-Gaussianity of a discrete population against a Gaussian field.

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
run.sh scripts/80_pull_high_ruwe_fast.py       # channel 20 - pull high-RUWE
run.sh scripts/81_pull_ruwe_astrometry_cols.py # channel 20 - zero-point columns
run.sh scripts/83_searchN_corrected.py         # channel 20 - corrected Search N
run.sh scripts/84_pull_ipd_harmonic.py       # channel 22 - pull IPD harmonic cols
run.sh scripts/85_searchV_harmonic_phase.py  # channel 22 - alignment (superseded verdict)
run.sh scripts/87_searchV_resolution_scan.py # channel 22 - resolution scan, the verdict
run.sh scripts/88_searchW_energy_locus.py    # channel 21 - energy-conservation locus
run.sh scripts/89_audit_pair_limit.py      # audit of the headline limit
run.sh scripts/45_grabby_figures.py            # figures
```

Figures: `F11_spectral_slope.png`, `F12_mirror_and_sky.png`,
`F13_channel_coverage.png`, `F14_searchB_ir_excess.png`,
`F15_searchA_diagnosis.png`, `F16_searchD_domain_edge.png`,
`F17_searchC_radio.png`, `F18_grand_summary.png`,
plus `F1`–`F10` from the main paper.
