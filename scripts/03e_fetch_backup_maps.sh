#!/usr/bin/env bash
# Fetch backup 3D dust maps from Zenodo (Harvard Dataverse is returning HTTP 202
# with an empty body for every path, so Bayestar19 is unreachable).
# Progress bars are dumped to a log file, not stdout.
set -uo pipefail
cd "$(dirname "$0")/.."
LOG="data/cache/dustmaps_fetch.log"
mkdir -p data/cache
"$HOME/sd-venv/bin/python" - >>"$LOG" 2>&1 <<'PY'
import sys
sys.path.insert(0, ".")
from pipeline import config as cfg
from dustmaps.config import config as dm
dm["data_dir"] = str(cfg.DUSTMAPS_DATA_DIR)

import dustmaps.leike2020
try:
    dustmaps.leike2020.fetch()
    print("LEIKE2020 OK")
except Exception as e:
    print(f"LEIKE2020 FAILED: {type(e).__name__}: {e}")
PY
echo "--- tail of $LOG ---"
grep -E "OK|FAILED|Error" "$LOG" | tail -20
du -sh "$HOME/dustmaps_data"/* 2>/dev/null
