# Stellar optical-deficit search

A search for starlight that never leaves its system — the optical *deficit*
channel for Dyson-type harvesting, as opposed to the mid-infrared *excess*
channel that essentially all prior surveys have used.

The motivation is geometric: an optimiser that beams its waste heat into a
narrow solid angle `Ω` is detectable in the infrared by only `Ω/4π` of
observers, whereas photons that were absorbed and never departed are missing
from every direction. The cost is sensitivity, and quantifying that cost is
what this repository does.

**Read `RESULTS.md` first, then `LIMITATIONS.md`.** `DECISIONS.md` records the
judgement calls.

## The measurement

Absolute `G` magnitude at fixed absolute `Ks` magnitude, for main-sequence
stars, after extinction and metallicity control. If a fraction `f` of a star's
optical output is intercepted,

```
Δ M_G = −2.5 log10(1 − f)          (positive = fainter)
```

`M_Ks` is the "how big is this star really" anchor. Whether that anchor holds
is the central question, and `pipeline/anchor.py` answers it quantitatively.

## Layout

```
pipeline/          modular, resumable, config-driven library
  config.py        every threshold and empirical constant, with citations
  adql.py          query construction and HEALPix sky partitioning
  tap.py           cached, resumable, thread-safe TAP client
  extinction.py    3D dust maps x band extinction laws
  sample.py        raw chunks -> analysis sample, with cut-flow bookkeeping
  fiducial.py      the M_G(M_Ks, ...) relation; robust spline fit
  statistics.py    test statistics and exclusion limits
  nulls.py         step 3: the systematic-floor measurement
  injection.py     step 4: injection-recovery
  anchor.py        spectral leverage and the flat-absorber question
  blind.py         step 5: blinding with a hash commitment
  figures.py       one figure per claim
scripts/           numbered, runnable stages
tests/             30 tests, including the ones that caught two errors
data/cache/        Parquet caches (gitignored) + manifest.jsonl (tracked)
results/           JSON and CSV outputs
figures/           PNGs, all regenerable by scripts/20_make_figures.py
```

## Running it

The pipeline runs under WSL because `healpy` has no Windows wheel and
`dustmaps` needs it. `run.sh` wraps the venv interpreter.

```bash
bash setup_env.sh                                   # one-time
wsl -d kali-linux bash run.sh scripts/03_fetch_dustmaps.py
wsl -d kali-linux bash run.sh scripts/02_pull_sample.py --workers 3 --distance-max-pc 500
wsl -d kali-linux bash run.sh scripts/05_build_sample.py --pattern 'sample_d500_p*' --tag primary
wsl -d kali-linux bash run.sh scripts/10_fit_fiducial.py --tag primary --nir-control
wsl -d kali-linux bash run.sh scripts/11_null_tests.py --tag primary
wsl -d kali-linux bash run.sh scripts/20_make_figures.py --tag primary
```

Full sequence in `RESULTS.md`. Tests:

```bash
wsl -d kali-linux bash -c 'cd /mnt/c/Users/neogo/Documents/StellarDeficit && ~/sd-venv/bin/python -m pytest tests/ -q'
```

## Three things worth knowing before reading the numbers

1. **A uniform harvesting fraction is not measurable by this method.** The
   fiducial is fit to the data, so a constant offset is absorbed exactly.
   Sensitivity is zero, not `sigma/sqrt(N)`. Demonstrated by injection across
   four orders of magnitude.
2. **Unresolved binaries mimic the signal**, they do not suppress it, because
   `dM_G/dM_Ks > 1`. The RUWE cut is a defence, not a safety margin.
3. **The spectral blind spot is not the grey absorber.** It sits at
   `α ≈ 0.19`; a grey absorber gives a residual of the *opposite* sign.

## Prior art

- Annis 1999, JBIS 52, 33 — originated the optical-deficit method
- Zackrisson et al. 2015, ApJ 810, 23 — Tully-Fisher, 3% galaxy-scale limit
- Wright 2023, ApJ 956, 34 — thermodynamics; small hot spheres favoured
- Suazo et al. 2024, MNRAS 531, 695 — Project Hephaistos II, the IR-excess
  pipeline this work inverts
