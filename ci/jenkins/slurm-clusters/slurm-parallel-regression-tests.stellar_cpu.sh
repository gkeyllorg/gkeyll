#!/bin/bash -l
# Stellar CPU parallel C-regression payload for the shared Slurm Pipeline.
# Jenkins supplies precompiled tests and the Slurm resource request.
# Write timing and comparison summaries to the shared CI_WORKSPACE.
set -euo pipefail

# Require the shared paths and the baseline-create/candidate-check mode.
: "${CI_WORKSPACE:?}"
: "${CI_BASELINE_DIR:?}"
: "${CI_BASELINE_PREFIX:?}"
: "${CI_REGRESSION_MODE:?}"

# Load the cluster environment before running either regression mode.
cd "$CI_WORKSPACE"
. machines/module_load.stellar-intel.sh
# Add site-specific module unloads or runtime exports here, before launching tests.

baseline_gkeyll="$CI_BASELINE_PREFIX/gkeyll/bin/gkeyll"

# The regression driver launches MPI tests; do not wrap it in another srun/mpirun.

# On a cache miss, generate the baseline outputs used for later comparisons.
if [[ "$CI_REGRESSION_MODE" == baseline-create ]]; then
  cd "$CI_BASELINE_DIR"
  started="$(date +%s)"
  "$baseline_gkeyll" runregression run -c --parallel --execute-only create
  elapsed="$(( $(date +%s) - started ))"
  printf '%s\n' "$elapsed" > "$CI_WORKSPACE/baseline-parallel-c-regression-create-seconds.txt"
  echo "Baseline parallel C-regression create runtime: $elapsed seconds"
  exit 0
fi

# Candidate comparison additionally needs its checkout, install and trusted checker.
[[ "$CI_REGRESSION_MODE" == candidate-check ]] || { echo "Unknown CI_REGRESSION_MODE: $CI_REGRESSION_MODE" >&2; exit 2; }
: "${CI_CANDIDATE_PREFIX:?}"
: "${CI_CANDIDATE_DIR:?}"
: "${CI_TRUSTED_CHECKER:?}"
candidate_gkeyll="$CI_CANDIDATE_PREFIX/gkeyll/bin/gkeyll"
test -s "$CI_TRUSTED_CHECKER"

cd "$CI_CANDIDATE_DIR"
# Link each solver's accepted baseline output into the candidate results.
for layer in moments vlasov gyrokinetic pkpm; do
  src="$CI_BASELINE_PREFIX/gkeyll-results/parallel-c-4/$layer/creg-accepted"
  dst="$CI_CANDIDATE_PREFIX/gkeyll-results/parallel-c-4/$layer/creg-accepted"
  rm -rf "$dst"
  if [[ -d "$src" ]]; then ln -s "$src" "$dst"; fi
done

# Run the candidate against the baseline outputs and record elapsed time.
started="$(date +%s)"
"$candidate_gkeyll" runregression run -c --parallel --execute-only check
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/candidate-parallel-c-regression-check-seconds.txt"
echo "Candidate parallel C-regression check runtime: $elapsed seconds"

# Write the comparison summary; the trusted checker determines CI pass/fail.
cd "$CI_BASELINE_DIR"
"$baseline_gkeyll" -S "$CI_TRUSTED_CHECKER" \
  "$CI_CANDIDATE_PREFIX/gkeyll-results/parallel-c-4" \
  "$CI_CANDIDATE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_WORKSPACE/ci-parallel-regression-summary.txt" \
  "$CI_BASELINE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_BASELINE_PREFIX/gkeyll-results/parallel-c-4"
