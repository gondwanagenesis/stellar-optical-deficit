#!/usr/bin/env bash
# One-time WSL environment build for the stellar optical-deficit search.
# healpy has no Windows wheel (needs cfitsio), so the pipeline runs under WSL.
set -euo pipefail

VENV="$HOME/sd-venv"

if [ ! -d "$VENV" ]; then
  python3 -m venv "$VENV"
fi

"$VENV/bin/pip" install -q --upgrade pip setuptools wheel

"$VENV/bin/pip" install -q \
  numpy scipy pandas pyarrow \
  astropy astroquery \
  scikit-learn matplotlib \
  healpy dustmaps h5py \
  requests tqdm pytest

"$VENV/bin/python" - <<'PY'
import importlib
for m in ["numpy","scipy","pandas","pyarrow","astropy","astroquery",
          "sklearn","matplotlib","healpy","dustmaps","h5py"]:
    mod = importlib.import_module(m)
    print(f"{m:12s} {getattr(mod,'__version__','?')}")
PY
