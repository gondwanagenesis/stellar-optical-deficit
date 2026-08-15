#!/usr/bin/env bash
# Everything AFTER the blind is opened. Run only once run_primary.sh has
# completed and its analysis state has been committed to git.
#
#   bash run_after_unblind.sh
#
# The unblinding step refuses to proceed without the confirmation string and
# verifies the SHA-256 commitment first, so a redrawn or edited offset aborts
# rather than silently producing a result.
set -euo pipefail
cd "$(dirname "$0")"
R="bash run.sh"

if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "ERROR: working tree is dirty. Commit the frozen analysis first --" >&2
  echo "       the whole point of the blind is that the analysis is fixed" >&2
  echo "       in git BEFORE the offset is revealed." >&2
  exit 1
fi

echo "=============== step 6: unblind ========================================"
$R scripts/14_analyse.py --tag primary --unblind

echo "=============== step 7: outlier follow-up (WISE + SIMBAD) =============="
$R scripts/15_outlier_followup.py --tag primary --max-objects 300

echo "=============== systematic budget ======================================"
$R scripts/22_systematic_budget.py --tag primary

echo "=============== distance trade study ==================================="
# Memory-heavy (loads the integrated dust map); do not run this alongside a
# fit on the full sample. See DECISIONS.md D1a.
$R scripts/17_distance_trade.py --knots 6 --mh-degree 1

echo "=============== figures ================================================"
$R scripts/20_make_figures.py --tag primary

echo "=============== collected numbers ======================================"
$R scripts/21_collect_numbers.py --tag primary | tee results/SUMMARY.txt

echo
echo "Done. RESULTS.md must record any analysis change made after this point"
echo "as post-unblinding, with both the pre- and post-change result."
