#!/usr/bin/env python
"""Watch for Gaia DR4 and re-run the channels it unblocks.

    run.sh scripts/63_check_dr4_and_rerun.py            # probe only
    run.sh scripts/63_check_dr4_and_rerun.py --rerun    # probe, then run

WHY A PROBE RATHER THAN A DATE
------------------------------
DR4 is announced for 2 December 2026, and Gaia release dates have slipped
before. Firing on the announced date would either run against data that is not
there or miss a quiet early release. So this asks the archive what tables it
actually has, and acts on the answer.

It also checks mirrors, because ESA's own service was returning nothing at all
on 2026-08-18 -- not even SELECT TOP 3 inside 120 s -- while ARI Heidelberg
answered in 1.1 s. On release day ESA will be hammered, so the mirror is the
more likely route in.

WHAT DR4 ACTUALLY UNBLOCKS, AND WHY EACH CHANNEL CARES
------------------------------------------------------
Channel 16, dark companions. This is the one that was left structurally
unfinished. The helium-white-dwarf mass-period relation only becomes a clean
test at P > 1500 d, and DR3's astrometric baseline stops near 1000 d, so every
system in the anomaly box had an extrapolated orbit with 31x worse relative
period error than the parent sample. DR4's 66-month baseline covers the region
where the test has power, and epoch astrometry allows fitting the orbit
ourselves rather than trusting the pipeline's model choice.

Channel 19, achromatic occultation. Never run, because DR3 publishes epoch
photometry for 11.75M variability-selected sources -- a biased sample far too
small for an all-sky dip search. DR4 publishes G/BP/RP time series for ~2.8
BILLION sources. That is the release that makes a colour-first transit search
possible at all, and the literature says no such search has ever been done.

Channel 13, coherent acceleration. The Hipparcos-Gaia baseline grows from
24.75 to 26.25 yr and DR4 proper motions improve by ~sqrt(2), which is only a
30-40% gain. The real prize is epoch astrometry: it allows separating a
constant acceleration from an orbital one per star, which attacks the
binary/engine degeneracy directly instead of statistically.

Channels 1-12, 17. More stars, better astrometry, deeper completeness. Channel
17 in particular is currently limited to 100 pc by volume-completeness, and
that radius grows with DR4.
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config as cfg

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("dr4watch")

MIRRORS = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap"),
    ("ESA", "https://gea.esac.esa.int/tap-server/tap"),
]

# Tables whose appearance means the corresponding channel can be re-run.
WANTED = {
    "gaiadr4.gaia_source":          "core catalogue (channels 1-12, 17)",
    "gaiadr4.nss_two_body_orbit":   "astrometric orbits (channel 16)",
    "gaiadr4.epoch_photometry":     "G/BP/RP time series (channel 19)",
    "gaiadr4.epoch_astrometry":     "per-transit astrometry (channels 13, 16)",
    "gaiadr4.binary_masses":        "dynamical masses (channels 4, 18)",
}

RERUN = [
    ("52d_pull_nss_paginated.py", ["--force"], "re-pull orbits"),
    ("56_searchE_dark_companions.py", ["--tag", "primary"],
     "channel 16 - dark companions"),
    ("62_searchE_diagnosis.py", ["--tag", "primary"],
     "channel 16 - diagnosis"),
]

STATE = cfg.RESULT_DIR / "dr4_watch_state.json"


def probe():
    """Ask each mirror which of the wanted tables exist."""
    from astroquery.utils.tap.core import TapPlus

    for name, url in MIRRORS:
        t0 = time.time()
        try:
            tap = TapPlus(url=url)
            rows = tap.launch_job(
                "SELECT table_name FROM TAP_SCHEMA.tables").get_results()
            names = {str(r["table_name"]).lower() for r in rows}
            log.info("%s reachable in %.1fs, %d tables",
                     name, time.time() - t0, len(names))
            found = {t: (t.lower() in names) for t in WANTED}
            any_dr4 = any(n.startswith("gaiadr4.") for n in names)
            return {"mirror": name, "url": url, "found": found,
                    "any_dr4_table": any_dr4,
                    "n_tables": len(names)}
        except Exception as exc:
            log.warning("%s unusable (%.1fs): %s", name,
                        time.time() - t0, str(exc)[:120])
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun", action="store_true",
                    help="run the dependent channels if DR4 orbits are live")
    args = ap.parse_args()

    now = datetime.now(timezone.utc).isoformat()
    res = probe()

    if res is None:
        log.error("no Gaia TAP mirror responded; will retry on the next tick")
        STATE.write_text(json.dumps(
            {"checked_utc": now, "status": "no_mirror"}, indent=2))
        return 1

    log.info("")
    log.info("DR4 table availability via %s:", res["mirror"])
    for t, present in res["found"].items():
        log.info("  [%s] %-32s %s", "x" if present else " ", t, WANTED[t])

    ready = res["found"].get("gaiadr4.nss_two_body_orbit", False)
    core = res["found"].get("gaiadr4.gaia_source", False)

    state = {
        "checked_utc": now,
        "mirror": res["mirror"],
        "any_dr4_table": res["any_dr4_table"],
        "tables": res["found"],
        "orbits_ready": ready,
        "core_ready": core,
    }
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(state, indent=2))

    if not (ready or core):
        print("\nGaia DR4 is not published yet. Nothing to do.")
        print("Channel 16 stays unfinished until it is: the helium-white-dwarf")
        print("test only has power at P > 1500 d, and DR3's baseline ends near")
        print("1000 d, so every candidate orbit is extrapolated.")
        return 0

    print(f"\n{'='*66}")
    print("GAIA DR4 IS LIVE")
    print(f"{'='*66}")
    for t, present in res["found"].items():
        if present:
            print(f"  available: {t:32s} -> {WANTED[t]}")

    if not args.rerun:
        print("\nRe-run with --rerun to execute the dependent channels.")
        return 0

    if not ready:
        print("\nCore catalogue is up but the orbit table is not; "
              "holding channel 16 until it appears.")
        return 0

    print("\nRe-running the channels DR4 unblocks ...\n")
    root = Path(__file__).resolve().parent
    py = sys.executable
    failures = []
    for script, extra, label in RERUN:
        log.info("--- %s (%s) ---", label, script)
        try:
            subprocess.run([py, "-u", str(root / script)] + extra,
                           cwd=str(cfg.PROJECT_ROOT), check=True)
        except subprocess.CalledProcessError as exc:
            log.error("%s failed: %s", script, exc)
            failures.append(script)

    state["reran_utc"] = datetime.now(timezone.utc).isoformat()
    state["rerun_failures"] = failures
    STATE.write_text(json.dumps(state, indent=2))

    print(f"\nDone. {len(RERUN) - len(failures)}/{len(RERUN)} steps succeeded.")
    if failures:
        print("failed:", ", ".join(failures))
    print("\nStill to build by hand, now that the data exists:")
    print("  channel 19 - colour-first achromatic transit search over the")
    print("               ~2.8e9 sources with DR4 epoch photometry. The")
    print("               literature has no systematic colour-selected dip")
    print("               search; DR3 was too small and too biased to try.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
