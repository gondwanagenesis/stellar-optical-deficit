#!/usr/bin/env bash
# Everything on the full primary sample UP TO the blinded analysis.
# Unblinding is deliberately NOT in here: the analysis must be frozen and
# committed to git between this script and scripts/14_analyse.py --unblind.
set -euo pipefail
cd "$(dirname "$0")"
R="bash run.sh"

# Cross-validation only selects two integers (knot count, metallicity degree),
# and those are not sensitive to N once the subsample is large compared with the
# number of parameters. 200k stars against ~30 parameters is ample, and it keeps
# the 27-point grid to ~20 min instead of ~40. The chosen complexity is then
# refit on the full 3.9M sample.
CV_N="${CV_MAX_N:-200000}"

echo "=============== step 2: fiducial (variant A, NIR control) ==============="
$R scripts/10_fit_fiducial.py --tag primary --out-tag primary --nir-control \
    --cv-max-n "$CV_N"

echo "=============== step 2b: fiducial (variant B, + optical colour) ========"
$R scripts/10_fit_fiducial.py --tag primary --out-tag primary_optcol \
    --nir-control --optical-colour-control --cv-max-n "$CV_N"

echo "=============== step 3: null tests / systematic floor =================="
$R scripts/11_null_tests.py --tag primary
$R scripts/11_null_tests.py --tag primary_optcol

echo "=============== spectral leverage ======================================"
$R scripts/12_spectral_leverage.py --tag primary

echo "=============== step 4: injection-recovery ============================="
# INJECTION_MAX_N is set from the outside so the recovery curve is always
# quoted at the N it was actually measured at. Recovery significance scales as
# p*sqrt(N), so a threshold measured on a subsample must NOT be rescaled to the
# full sample and claimed -- see the rule in the brief.
$R scripts/13_injection.py --tag primary \
    ${INJECTION_MAX_N:+--max-n "$INJECTION_MAX_N"} \
    --n-realisations "${INJECTION_REALISATIONS:-6}"

echo "=============== step 5: BLINDED analysis ==============================="
$R scripts/14_analyse.py --tag primary --blinded

echo
echo "STOP. Freeze and commit the analysis before unblinding:"
echo "    git add -A && git commit -m 'freeze analysis before unblinding'"
echo "    bash run.sh scripts/14_analyse.py --tag primary --unblind"
