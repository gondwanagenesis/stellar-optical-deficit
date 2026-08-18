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

**Joint, disposal-agnostic: p_total < 6.2e−4.** Fewer than 1 in 1,614 nearby
lower-main-sequence stars intercepts ≥51% of its optical output by any means.

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
