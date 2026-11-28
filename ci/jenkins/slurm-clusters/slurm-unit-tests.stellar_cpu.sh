#!/bin/bash -l
# Stellar CPU unit-test payload for the shared Slurm Pipeline.
#
# This script deliberately has no #SBATCH directives. The trusted Jenkins
# pipeline owns the resource request, while this script owns only the runtime
# environment and the test command.
set -euo pipefail

# Jenkins supplies the shared output directory and candidate checkout.
: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"
: "${CI_CANDIDATE_DIR:?CI_CANDIDATE_DIR must name the candidate checkout}"

cd "$CI_CANDIDATE_DIR"

# Load the build environment explicitly in this batch shell.
. machines/module_load.stellar-intel.sh
# Add site-specific module unloads or runtime exports here, before launching tests.

# Run the unit suite and write per-test results for the Jenkins report.
started="$(date +%s)"
GKYL_UNIT_RESULTS="$CI_WORKSPACE/candidate-unit-results.txt" make unit-run
# Save elapsed time alongside the test results.
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/unit-test-seconds.txt"
echo "Unit-test runtime: $elapsed seconds"
