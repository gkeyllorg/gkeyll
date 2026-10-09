#!/bin/bash -l
# Unit-test payload for ci/jenkins/jenkinsfile.stellar_cpu.
#
# This script deliberately has no #SBATCH directives. The trusted Jenkins
# pipeline owns the resource request, while this script owns only the runtime
# environment and the test command.
set -euo pipefail

: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"
: "${CI_CANDIDATE_DIR:?CI_CANDIDATE_DIR must name the candidate checkout}"

cd "$CI_CANDIDATE_DIR"

# A batch shell does not inherit a user's interactive module selection in a
# reliable or reproducible way. Source the same environment used to configure
# and build this checkout.
. machines/module_load.stellar-intel.sh

started="$(date +%s)"
GKYL_UNIT_RESULTS="$CI_WORKSPACE/candidate-unit-results.txt" make unit-run
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/unit-test-seconds.txt"
echo "Unit-test runtime: $elapsed seconds"
