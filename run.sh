#!/usr/bin/env bash
# Thin wrapper so Windows-side callers never have to fight quoting.
#   wsl -d kali-linux bash /mnt/c/Users/neogo/Documents/StellarDeficit/run.sh scripts/foo.py --arg
set -euo pipefail
cd "$(dirname "$0")"
exec "$HOME/sd-venv/bin/python" "$@"
