#!/bin/bash -l
# Stellar CPU serial C-regression payload for the shared Slurm Pipeline.
#
# Jenkins supplies precompiled tests and the Slurm resource request.
# Write timing and comparison summaries to the shared CI_WORKSPACE.
set -euo pipefail

# Require the shared paths and the baseline-create/candidate-check mode.
: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"
: "${CI_BASELINE_DIR:?CI_BASELINE_DIR must name the trusted baseline checkout}"
: "${CI_BASELINE_PREFIX:?CI_BASELINE_PREFIX must name the baseline install prefix}"
: "${CI_REGRESSION_MODE:?CI_REGRESSION_MODE must be baseline-create or candidate-check}"
: "${CI_SERIAL_REGRESSION_JOBS:?CI_SERIAL_REGRESSION_JOBS must be set}"
: "${CI_SERIAL_REGRESSION_TEST_TIMEOUT:?CI_SERIAL_REGRESSION_TEST_TIMEOUT must be set}"

# Load the cluster environment before running either regression mode.
cd "$CI_WORKSPACE"
. machines/module_load.stellar-intel.sh
# Add site-specific module unloads or runtime exports here, before launching tests.

baseline_gkeyll="$CI_BASELINE_PREFIX/gkeyll/bin/gkeyll"

test -x "$baseline_gkeyll"

# On a cache miss, generate the baseline outputs used for later comparisons.
if [[ "$CI_REGRESSION_MODE" == baseline-create ]]; then
  cd "$CI_BASELINE_DIR"
  started="$(date +%s)"
  "$baseline_gkeyll" runregression run -c --execute-only create \
    --jobs "$CI_SERIAL_REGRESSION_JOBS" \
    --timeout "$CI_SERIAL_REGRESSION_TEST_TIMEOUT"
  elapsed="$(( $(date +%s) - started ))"
  printf '%s\n' "$elapsed" > "$CI_WORKSPACE/baseline-serial-c-regression-create-seconds.txt"
  echo "Baseline serial C-regression create runtime: $elapsed seconds"
  exit 0
fi

# Candidate comparison additionally needs its checkout, install and trusted checker.
[[ "$CI_REGRESSION_MODE" == candidate-check ]] || { echo "Unknown CI_REGRESSION_MODE: $CI_REGRESSION_MODE" >&2; exit 2; }
: "${CI_CANDIDATE_PREFIX:?CI_CANDIDATE_PREFIX must name the candidate install prefix}"
: "${CI_CANDIDATE_DIR:?CI_CANDIDATE_DIR must name the candidate checkout}"
: "${CI_TRUSTED_CHECKER:?CI_TRUSTED_CHECKER must name the reviewed regression checker}"
candidate_gkeyll="$CI_CANDIDATE_PREFIX/gkeyll/bin/gkeyll"
test -x "$candidate_gkeyll"
test -s "$CI_TRUSTED_CHECKER"

cd "$CI_CANDIDATE_DIR"
# Link each solver's accepted C baseline output into the candidate results.
for layer in moments vlasov gyrokinetic pkpm; do
  baseline_accepted="$CI_BASELINE_PREFIX/gkeyll-results/$layer/creg-accepted"
  candidate_accepted="$CI_CANDIDATE_PREFIX/gkeyll-results/$layer/creg-accepted"
  rm -rf "$candidate_accepted"
  if [[ -d "$baseline_accepted" ]]; then
    ln -s "$baseline_accepted" "$candidate_accepted"
  fi
done

# Run the candidate against the baseline outputs and record elapsed time.

started="$(date +%s)"
"$candidate_gkeyll" runregression run -c --execute-only check \
  --jobs "$CI_SERIAL_REGRESSION_JOBS" \
  --timeout "$CI_SERIAL_REGRESSION_TEST_TIMEOUT"
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/candidate-serial-c-regression-check-seconds.txt"
echo "Candidate serial C-regression check runtime: $elapsed seconds"

# Write the comparison summary; the trusted checker determines CI pass/fail.
cd "$CI_BASELINE_DIR"
"$baseline_gkeyll" -S "$CI_TRUSTED_CHECKER" \
  "$CI_CANDIDATE_PREFIX/gkeyll-results" \
  "$CI_CANDIDATE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_WORKSPACE/ci-serial-regression-summary.txt" \
  "$CI_BASELINE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_BASELINE_PREFIX/gkeyll-results"
