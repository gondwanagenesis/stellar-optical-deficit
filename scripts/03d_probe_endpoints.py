#!/usr/bin/env python
"""Probe alternative hosts for the Bayestar19 map and for backup 3D dust maps."""
import requests

UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
      "Chrome/126.0.0.0 Safari/537.36")
H = {"User-Agent": UA}

TARGETS = [
    ("dataverse root",        "https://dataverse.harvard.edu/"),
    ("dataverse api info",    "https://dataverse.harvard.edu/api/info/version"),
    ("dataverse dataset pg",  "https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/2EJ9TX"),
    ("bayestar direct doi",   "https://doi.org/10.7910/DVN/2EJ9TX"),
    ("zenodo leike2020",      "https://zenodo.org/api/records/3993082"),
    ("zenodo edenhofer",      "https://zenodo.org/api/records/8187943"),
    ("argonaut",              "http://argonaut.skymaps.info/"),
]

for label, url in TARGETS:
    try:
        r = requests.get(url, headers=H, timeout=45, allow_redirects=True,
                         stream=True)
        body = ""
        if "json" not in (r.headers.get("content-type") or ""):
            body = r.raw.read(120, decode_content=True).decode("utf8", "replace")
            body = " ".join(body.split())[:100]
        print(f"{label:22s} HTTP {r.status_code:3d}  "
              f"ct={(r.headers.get('content-type') or '')[:30]:30s} {body}")
        r.close()
    except Exception as exc:                            # noqa: BLE001
        print(f"{label:22s} {type(exc).__name__}: {str(exc)[:90]}")
