#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
sed -i 's/\r$//' scripts/fix_mojibake.py
bash run.sh scripts/fix_mojibake.py PAPER.md README.md RESULTS.md LIMITATIONS.md DECISIONS.md
echo "--- section renumber (UTF-8 safe) ---"
sed -i 's/^### 6\.5 Errors we made/### 6.6 Errors we made/' PAPER.md
grep -n '^### 6\.' PAPER.md || true
echo "--- residual mojibake markers ---"
grep -c 'Ã\|Î\|â€' PAPER.md || echo "0 (clean)"
