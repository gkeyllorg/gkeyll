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
if [ -n "${GKEYLL_CI_DEPENDENCY_ENV:-}" ]; then . "$GKEYLL_CI_DEPENDENCY_ENV"; fi

if [ "${GKEYLL_USE_VALGRIND:-0}" = 1 ]; then
    started="$(date +%s)"
    result=0
    python3 -I "$CI_VALGRIND_SCRIPT" run --output "$CI_WORKSPACE/ci-valgrind/candidate" || result=$?
    printf '%s\n' "$(( $(date +%s) - started ))" > "$CI_WORKSPACE/candidate-valgrind-seconds.txt"
    [ "$result" -eq 0 ] || exit "$result"
fi

started="$(date +%s)"
GKYL_UNIT_RESULTS="$CI_WORKSPACE/candidate-unit-results.txt" make unit-run
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/unit-test-seconds.txt"
echo "Unit-test runtime: $elapsed seconds"
