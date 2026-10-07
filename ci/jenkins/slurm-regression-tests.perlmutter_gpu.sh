#!/bin/bash -l
# All-C regression payload for ci/jenkins/jenkinsfile.perlmutter_gpu.
#
# Jenkins builds, installs, and compiles the candidate and trusted baseline on
# login node. This payload only executes their precompiled C tests, creates
# the baseline's outputs, checks the candidate against them, and turns SQLite-recorded
# failures into a Slurm/Jenkins failure.
set -euo pipefail

: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"
: "${CI_BASELINE_DIR:?CI_BASELINE_DIR must name the trusted baseline checkout}"
: "${CI_BASELINE_PREFIX:?CI_BASELINE_PREFIX must name the baseline install prefix}"
: "${CI_REGRESSION_MODE:?CI_REGRESSION_MODE must be baseline-create or candidate-check}"
: "${CI_REGRESSION_JOBS:?CI_REGRESSION_JOBS must be set}"
: "${CI_REGRESSION_TEST_TIMEOUT:?CI_REGRESSION_TEST_TIMEOUT must be set}"

cd "$CI_WORKSPACE"
. machines/module_load.perlmutter-gpu.sh
export SLURM_CPU_BIND=cores

baseline_gkeyll="$CI_BASELINE_PREFIX/gkeyll/bin/gkeyll"

test -x "$baseline_gkeyll"

if [[ "$CI_REGRESSION_MODE" == baseline-create ]]; then
  cd "$CI_BASELINE_DIR"
  started="$(date +%s)"
  srun --ntasks=1 --cpus-per-task=32 --gpus-per-task=1 --cpu-bind=cores \
    "$baseline_gkeyll" runregression run -c --execute-only create \
    --jobs "$CI_REGRESSION_JOBS" \
    --timeout "$CI_REGRESSION_TEST_TIMEOUT"
  elapsed="$(( $(date +%s) - started ))"
  printf '%s\n' "$elapsed" > "$CI_WORKSPACE/baseline-c-regression-create-seconds.txt"
  echo "Baseline C-regression create runtime: $elapsed seconds"
  exit 0
fi

[[ "$CI_REGRESSION_MODE" == candidate-check ]] || { echo "Unknown CI_REGRESSION_MODE: $CI_REGRESSION_MODE" >&2; exit 2; }
: "${CI_CANDIDATE_PREFIX:?CI_CANDIDATE_PREFIX must name the candidate install prefix}"
: "${CI_CANDIDATE_DIR:?CI_CANDIDATE_DIR must name the candidate checkout}"
candidate_gkeyll="$CI_CANDIDATE_PREFIX/gkeyll/bin/gkeyll"
test -x "$candidate_gkeyll"

cd "$CI_CANDIDATE_DIR"
# Use only C accepted output from the fixed baseline. Do not carry Lua
# baselines into this CPU C-regression job.
for layer in moments vlasov gyrokinetic pkpm; do
  baseline_accepted="$CI_BASELINE_PREFIX/gkeyll-results/$layer/creg-accepted"
  candidate_accepted="$CI_CANDIDATE_PREFIX/gkeyll-results/$layer/creg-accepted"
  rm -rf "$candidate_accepted"
  if [[ -d "$baseline_accepted" ]]; then
    ln -s "$baseline_accepted" "$candidate_accepted"
  fi
done

started="$(date +%s)"
srun --ntasks=1 --cpus-per-task=32 --gpus-per-task=1 --cpu-bind=cores \
  "$candidate_gkeyll" runregression run -c --execute-only check \
  --jobs "$CI_REGRESSION_JOBS" \
  --timeout "$CI_REGRESSION_TEST_TIMEOUT"
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/candidate-c-regression-check-seconds.txt"
echo "Candidate C-regression check runtime: $elapsed seconds"
"$baseline_gkeyll" "$CI_BASELINE_DIR/ci/jenkins/check_regression_results.lua" \
  "$CI_CANDIDATE_PREFIX/gkeyll-results" \
  "$CI_CANDIDATE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_WORKSPACE/ci-regression-summary.txt" \
  "$CI_BASELINE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_BASELINE_PREFIX/gkeyll-results"
