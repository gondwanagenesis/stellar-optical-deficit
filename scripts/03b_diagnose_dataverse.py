#!/usr/bin/env python
"""Diagnose why dustmaps' Dataverse fetch fails for Bayestar19 / SFD."""
from __future__ import annotations
import json
import sys
import requests

DOIS = {
    "bayestar2019": "doi:10.7910/DVN/2EJ9TX",
    "sfd": "doi:10.7910/DVN/EWCNL5",
}
API = "https://dataverse.harvard.edu/api/datasets/:persistentId"

for name, doi in DOIS.items():
    print(f"\n=== {name}  {doi}")
    try:
        r = requests.get(API, params={"persistentId": doi}, timeout=60)
        print(f"  HTTP {r.status_code}  content-type={r.headers.get('content-type')}")
        body = r.text[:400].replace("\n", " ")
        print(f"  body[:400]: {body}")
        if r.headers.get("content-type", "").startswith("application/json"):
            d = r.json()
            files = d["data"]["latestVersion"]["files"]
            print(f"  {len(files)} files:")
            for f in files:
                df = f["dataFile"]
                print(f"    id={df['id']}  {df['filename']}  "
                      f"{df.get('filesize', 0)/1e9:.2f} GB  md5={df.get('md5')}")
    except Exception as exc:                            # noqa: BLE001
        print(f"  {type(exc).__name__}: {exc}")
