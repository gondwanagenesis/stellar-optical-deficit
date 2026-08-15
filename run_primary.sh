#!/usr/bin/env bash
# Everything on the full primary sample UP TO the blinded analysis.
# Unblinding is deliberately NOT in here: the analysis must be frozen and
# committed to git between this script and scripts/14_analyse.py --unblind.
set -euo pipefail
cd "$(dirname "$0")"
R="bash run.sh"

# Spline complexity was selected by 5-fold CV on a 200k subsample; the measured
# loss surface is in results/cv_primary.csv and is monotone in both knots and
# metallicity degree, with its minimum at the grid corner (6, 1). See
# DECISIONS.md D9a for why the grid was truncated there rather than run to 40
# knots. Passing the selection explicitly keeps the expensive search out of
# every re-run while leaving the evidence on disk.
KNOTS="${FIT_KNOTS:-6}"
MHDEG="${FIT_MH_DEGREE:-1}"

echo "=============== step 2: fiducial (variant A, NIR control) ==============="
$R scripts/10_fit_fiducial.py --tag primary --out-tag primary --nir-control \
    --knots "$KNOTS" --mh-degree "$MHDEG"

echo "=============== step 2b: fiducial (variant B, + optical colour) ========"
$R scripts/10_fit_fiducial.py --tag primary --out-tag primary_optcol \
    --nir-control --optical-colour-control --knots "$KNOTS" --mh-degree "$MHDEG"

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
