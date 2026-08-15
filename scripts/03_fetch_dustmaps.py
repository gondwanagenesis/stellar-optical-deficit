#!/usr/bin/env python
"""Download the 3D dust maps.

We deliberately use two *independent* 3D maps and treat their disagreement as a
systematic error term rather than averaging it away:

  Bayestar19   Green, Schlafly, Finkbeiner et al. 2019, ApJ 887, 93.
               Pan-STARRS 1 + 2MASS + Gaia parallaxes, dec > -30 deg
               (~three quarters of the sky), nside 64-1024 adaptive.

  Edenhofer23  Edenhofer, Zucker, Frank et al. 2024, A&A 685, A82.
               Gaia DR3 XP-spectra based, full sky within 1.25 kpc,
               nside 256, 516 log-spaced distance bins.

SFD98 (Schlegel, Finkbeiner & Davis 1998, ApJ 500, 525) is fetched as a cheap
2D sanity check only; it integrates through the whole Galaxy and therefore
*overestimates* extinction for nearby stars, so it is never used for the
science measurement.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config as cfg

# Point dustmaps at native ext4 before importing any of its submodules.
cfg.DUSTMAPS_DATA_DIR.mkdir(parents=True, exist_ok=True)
from dustmaps.config import config as dm_config      # noqa: E402
dm_config["data_dir"] = str(cfg.DUSTMAPS_DATA_DIR)
print(f"dustmaps data_dir = {dm_config['data_dir']}")


def try_fetch(label: str, fn, **kw) -> None:
    print(f"\n=== fetching {label} ...", flush=True)
    try:
        fn(**kw)
        print(f"=== {label} OK", flush=True)
    except Exception as exc:                            # noqa: BLE001
        print(f"=== {label} FAILED: {type(exc).__name__}: {exc}", flush=True)


if __name__ == "__main__":
    import dustmaps.sfd
    import dustmaps.bayestar
    import dustmaps.edenhofer2023

    try_fetch("SFD98 (2D, sanity check only)", dustmaps.sfd.fetch)
    try_fetch("Bayestar19", dustmaps.bayestar.fetch, version="bayestar2019")
    # fetch_samples=False downloads mean+std only (~3 GB instead of ~30 GB);
    # we do not need the posterior samples because the map-to-map difference
    # dominates the within-map posterior width at these distances.
    try_fetch("Edenhofer2023", dustmaps.edenhofer2023.fetch, fetch_samples=False)

    print("\n--- files on disk ---")
    for p in sorted(cfg.DUSTMAPS_DATA_DIR.rglob("*")):
        if p.is_file():
            print(f"  {p.relative_to(cfg.DUSTMAPS_DATA_DIR)}  "
                  f"{p.stat().st_size / 1e9:.2f} GB")
