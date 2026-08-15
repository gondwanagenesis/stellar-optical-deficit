#!/usr/bin/env python
"""Fetch the official Gaia EDR3 colour-dependent extinction-law coefficients.

Source: ESA Gaia "EDR3 extinction law" page, file
Fitz19_EDR3_extinctionlawcoefficients.zip, which tabulates

    k_m = A_m / A_0
        = a1 + a2*X + a3*X^2 + a4*X^3
             + a5*A0 + a6*A0^2 + a7*A0^3
             + a8*A0*X + a9*A0*X^2 + a10*A0^2*X

for m in {G, BP, RP, J, H, Ks, G_RVS}, fitted on Fitzpatrick et al. (2019,
ApJ 886, 108) extinction curves with R_V = 3.1, over 3500-10000 K and
0.01 < A0 < 20 mag.  X is one of (BP-RP)_0, (G-K)_0 or Teff/5040 K.

We need this because A_G varies by ~10 per cent across the lower main sequence
at fixed A_0.  Using a single constant A_G/A_V would imprint a colour-dependent
residual that is *exactly degenerate* with a mass-dependent harvesting signal.
"""
from __future__ import annotations

import io
import sys
import zipfile
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline import config as cfg  # noqa: E402

OUT = cfg.DATA_DIR / "extinction_coeffs"
OUT.mkdir(parents=True, exist_ok=True)

CANDIDATE_URLS = [
    "https://www.cosmos.esa.int/documents/29201/1770596/Fitz19_EDR3_extinctionlawcoefficients.zip",
    "https://www.cosmos.esa.int/documents/29201/1769576/Fitz19_EDR3_extinctionlawcoefficients.zip",
    "https://www.cosmos.esa.int/documents/29201/1769576/Fitz19_EDR3_extinctionlawcoefficients.zip/",
]

HEADERS = {"User-Agent": "Mozilla/5.0 (research pipeline; contact via local use)"}


def find_zip_url() -> str | None:
    """Scrape the extinction-law page for the coefficient zip link."""
    try:
        r = requests.get("https://www.cosmos.esa.int/web/gaia/edr3-extinction-law",
                         headers=HEADERS, timeout=60)
        r.raise_for_status()
    except Exception as exc:                            # noqa: BLE001
        print(f"could not load page: {exc}")
        return None
    import re
    hits = re.findall(r'href="([^"]*extinctionlawcoefficients[^"]*\.zip[^"]*)"',
                      r.text, flags=re.I)
    if not hits:
        hits = re.findall(r'href="([^"]*\.zip[^"]*)"', r.text, flags=re.I)
    for h in hits:
        url = h if h.startswith("http") else "https://www.cosmos.esa.int" + h
        print(f"  found candidate link: {url}")
        return url
    return None


def main() -> int:
    urls = []
    scraped = find_zip_url()
    if scraped:
        urls.append(scraped)
    urls += CANDIDATE_URLS

    for url in urls:
        print(f"trying {url}")
        try:
            r = requests.get(url, headers=HEADERS, timeout=180)
            if r.status_code != 200 or len(r.content) < 1000:
                print(f"  -> HTTP {r.status_code}, {len(r.content)} bytes")
                continue
            zf = zipfile.ZipFile(io.BytesIO(r.content))
        except Exception as exc:                        # noqa: BLE001
            print(f"  -> {type(exc).__name__}: {exc}")
            continue

        print(f"  -> OK, {len(r.content)/1e3:.0f} kB")
        zf.extractall(OUT)
        for n in zf.namelist():
            print(f"     {n}")
        (OUT / "SOURCE_URL.txt").write_text(url + "\n")
        return 0

    print("\nFAILED to obtain the official coefficient file.")
    print("The pipeline will fall back to the Wang & Chen (2019) constant "
          "ratios plus an explicit colour-term systematic; see LIMITATIONS.md.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
