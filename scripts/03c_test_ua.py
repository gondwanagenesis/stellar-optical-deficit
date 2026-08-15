#!/usr/bin/env python
"""Test whether a browser User-Agent gets past the Dataverse 202 gate."""
import requests

UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
      "Chrome/126.0.0.0 Safari/537.36")
API = "https://dataverse.harvard.edu/api/datasets/:persistentId"

for label, headers in [("no UA", {}),
                       ("browser UA", {"User-Agent": UA, "Accept": "application/json"})]:
    r = requests.get(API, params={"persistentId": "doi:10.7910/DVN/2EJ9TX"},
                     headers=headers, timeout=60)
    print(f"{label:12s} HTTP {r.status_code}  ct={r.headers.get('content-type')}  "
          f"len={len(r.content)}")
    if r.status_code == 200 and "json" in (r.headers.get("content-type") or ""):
        files = r.json()["data"]["latestVersion"]["files"]
        for f in files:
            df = f["dataFile"]
            print(f"    id={df['id']:>10}  {df['filename']:<40} "
                  f"{df.get('filesize',0)/1e9:6.2f} GB")
