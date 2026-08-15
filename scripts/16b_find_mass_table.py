#!/usr/bin/env python
"""Find the Gaia table that carries dynamical masses / NSS orbits."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline.tap import tap_client

# A WHERE ... LIKE on TAP_SCHEMA.tables triggers a query-rewriting bug in the
# ESA service ("missing FROM-clause entry for table cte_authorisation_NNN"), so
# fetch the whole list and filter locally.
# Do not select `description`: some rows contain non-ASCII bytes that the
# VOTable BINARY parser decodes as ascii and dies on.
rows = tap_client().launch_job(
    "SELECT table_name FROM TAP_SCHEMA.tables ORDER BY table_name").get_results()
names = sorted(str(r["table_name"]) for r in rows)
print(f"{len(names)} tables visible; matches:")
for n in names:
    if any(k in n.lower() for k in ("mass", "nss", "binar", "orbit", "twobody")):
        print("  ", n)

print("\ncolumns of gaiadr3.nss_two_body_orbit (if present):")
q2 = ("SELECT column_name, datatype FROM TAP_SCHEMA.columns "
      "WHERE table_name = 'gaiadr3.nss_two_body_orbit' ORDER BY column_name")
try:
    r = tap_client().launch_job(q2).get_results()
    print("  " + ", ".join(str(c["column_name"]) for c in r))
except Exception as e:
    print("  ", e)
