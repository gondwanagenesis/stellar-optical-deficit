#!/usr/bin/env bash
# UTF-8 safe housekeeping. NEVER do this from PowerShell.
set -euo pipefail
cd "$(dirname "$0")"
sed -i 's/^### 5\.8 Search for a spreading front/### 5.10 Search for a spreading front/' PAPER.md
echo "--- section 5 list ---"
grep -n '^### 5\.' PAPER.md
echo "--- encoding ---"
if grep -q 'Ã\|Î\|â€' PAPER.md; then echo "MOJIBAKE PRESENT"; else echo "clean"; fi
echo "--- word count ---"
wc -w PAPER.md
