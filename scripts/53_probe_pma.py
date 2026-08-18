#!/usr/bin/env python
"""Probe the Hipparcos-Gaia proper-motion-anomaly catalogues for Search F.

    run.sh scripts/53_probe_pma.py

WHY
---
Gaia DR3 publishes positions and velocities, not accelerations. But the
difference between the long-baseline Hipparcos-to-Gaia mean proper motion and
the near-instantaneous Gaia proper motion IS an acceleration proxy, measured
over a ~25 yr baseline. Kervella et al. (2022) and Brandt (2021) both publish
it for the Hipparcos sample.

This matters because Search F needs a vector per star, not a scalar. Unseen
binaries produce proper-motion anomalies pointing in RANDOM directions;
coordinated thrust would produce ALIGNED ones over a region. The direction is
the whole discriminant, so we need the catalogue that keeps the vector.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

CANDIDATES = [
    ("Kervella+2022 EDR3 PMa", "J/A+A/657/A7"),
    ("Brandt 2021 HGCA EDR3", "J/ApJS/254/42"),
]


def main() -> int:
    from astroquery.vizier import Vizier

    for label, cat in CANDIDATES:
        print(f"\n{'='*70}\n{label}   ({cat})\n{'='*70}")
        try:
            cats = Vizier.find_catalogs(cat)
            if not cats:
                print("  NOT FOUND")
                continue
            for key, entry in cats.items():
                print(f"  catalog key: {key}")
                print(f"  description: {entry.description}")

            v = Vizier(columns=["**"], row_limit=3)
            tables = v.get_catalogs(cat)
            print(f"\n  {len(tables)} table(s):")
            for t in tables:
                name = t.meta.get("name", "?")
                print(f"\n  --- {name}  ({len(t)} rows shown) ---")
                print(f"      columns ({len(t.colnames)}):")
                for c in t.colnames:
                    unit = str(t[c].unit) if t[c].unit else ""
                    desc = (t[c].description or "")[:60]
                    print(f"        {c:18s} {unit:12s} {desc}")
        except Exception as exc:
            print(f"  ERROR: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
