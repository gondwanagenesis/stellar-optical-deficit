#!/usr/bin/env bash
set -euo pipefail
PY="$HOME/sd-venv/bin/python"
"$HOME/sd-venv/bin/pip" install -q gaiadr3-zeropoint
"$PY" - <<'PY'
from zero_point import zpt
zpt.load_tables()
print("gaiadr3-zeropoint OK (Lindegren et al. 2021, A&A 649, A4)")
PY
