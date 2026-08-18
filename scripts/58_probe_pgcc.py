#!/usr/bin/env python
"""Probe the Planck cold-clump catalogue for Search J.

    run.sh scripts/58_probe_pgcc.py

Search J needs per-source fluxes in at least three submillimetre bands so the
dust temperature and the emissivity index can be fitted independently. The
PGCC's own fit holds beta fixed at 2, so the published temperature column is
not usable for our purpose -- we need the raw band fluxes to refit.

This checks what is actually in the catalogue before the search is written.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

CANDIDATES = [
    ("Planck cold clumps (PGCC)", "J/A+A/594/A28"),
    ("Planck PCCS2", "J/A+A/594/A26"),
]


def main() -> int:
    from astroquery.vizier import Vizier

    for label, cat in CANDIDATES:
        print(f"\n{'='*72}\n{label}   ({cat})\n{'='*72}")
        try:
            found = Vizier.find_catalogs(cat)
            for key, entry in found.items():
                print(f"  {key}: {entry.description}")

            v = Vizier(columns=["**"], row_limit=2)
            tables = v.get_catalogs(cat)
            print(f"\n  {len(tables)} table(s)")
            for t in tables:
                name = t.meta.get("name", "?")
                print(f"\n  --- {name} ---")
                print(f"      {len(t.colnames)} columns")
                for c in t.colnames:
                    unit = str(t[c].unit) if t[c].unit else ""
                    desc = (t[c].description or "")[:58]
                    print(f"        {c:20s} {unit:10s} {desc}")
        except Exception as exc:
            print(f"  ERROR: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
