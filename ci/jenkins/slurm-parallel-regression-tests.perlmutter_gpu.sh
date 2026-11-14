#!/bin/bash -l
set -euo pipefail

: "${CI_WORKSPACE:?}"
: "${CI_BASELINE_DIR:?}"
: "${CI_BASELINE_PREFIX:?}"
: "${CI_REGRESSION_MODE:?}"

cd "$CI_WORKSPACE"
. machines/module_load.perlmutter-gpu.sh
export SLURM_CPU_BIND=cores

baseline_gkeyll="$CI_BASELINE_PREFIX/gkeyll/bin/gkeyll"

if [[ "$CI_REGRESSION_MODE" == baseline-create ]]; then
  cd "$CI_BASELINE_DIR"
  started="$(date +%s)"
  "$baseline_gkeyll" runregression run -c --parallel --execute-only create
  elapsed="$(( $(date +%s) - started ))"
  printf '%s\n' "$elapsed" > "$CI_WORKSPACE/baseline-parallel-c-regression-create-seconds.txt"
  echo "Baseline parallel C-regression create runtime: $elapsed seconds"
  exit 0
fi

[[ "$CI_REGRESSION_MODE" == candidate-check ]] || { echo "Unknown CI_REGRESSION_MODE: $CI_REGRESSION_MODE" >&2; exit 2; }
: "${CI_CANDIDATE_PREFIX:?}"
: "${CI_CANDIDATE_DIR:?}"
: "${CI_TRUSTED_CHECKER:?}"
candidate_gkeyll="$CI_CANDIDATE_PREFIX/gkeyll/bin/gkeyll"
test -s "$CI_TRUSTED_CHECKER"

cd "$CI_CANDIDATE_DIR"
for layer in moments vlasov gyrokinetic pkpm; do
  src="$CI_BASELINE_PREFIX/gkeyll-results/parallel-c-4/$layer/creg-accepted"
  dst="$CI_CANDIDATE_PREFIX/gkeyll-results/parallel-c-4/$layer/creg-accepted"
  rm -rf "$dst"
  if [[ -d "$src" ]]; then ln -s "$src" "$dst"; fi
done
started="$(date +%s)"
"$candidate_gkeyll" runregression run -c --parallel --execute-only check
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/candidate-parallel-c-regression-check-seconds.txt"
echo "Candidate parallel C-regression check runtime: $elapsed seconds"
cd "$CI_BASELINE_DIR"
"$baseline_gkeyll" -S "$CI_TRUSTED_CHECKER" \
  "$CI_CANDIDATE_PREFIX/gkeyll-results/parallel-c-4" \
  "$CI_CANDIDATE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_WORKSPACE/ci-parallel-regression-summary.txt" \
  "$CI_BASELINE_DIR/ci/jenkins/expected_regression_diffs.txt" \
  "$CI_BASELINE_PREFIX/gkeyll-results/parallel-c-4"
