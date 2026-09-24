#!/bin/bash -l
set -euo pipefail

: "${CI_WORKSPACE:?}"
: "${CI_BASELINE_DIR:?}"
: "${CI_BASELINE_PREFIX:?}"
: "${CI_CANDIDATE_PREFIX:?}"

cd "$CI_WORKSPACE"
. machines/module_load.perlmutter-gpu.sh
export SLURM_CPU_BIND=cores

baseline_gkeyll="$CI_BASELINE_PREFIX/gkeyll/bin/gkeyll"
candidate_gkeyll="$CI_CANDIDATE_PREFIX/gkeyll/bin/gkeyll"

cd "$CI_BASELINE_DIR"
started="$(date +%s)"
"$baseline_gkeyll" runregression run -c --parallel --execute-only create
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/baseline-parallel-c-regression-create-seconds.txt"
echo "Baseline parallel C-regression create runtime: $elapsed seconds"

cd "$CI_WORKSPACE"
for layer in moments vlasov gyrokinetic pkpm; do
  src="$CI_BASELINE_PREFIX/gkeyll-results/parallel-c-4/$layer/creg-accepted"
  dst="$CI_CANDIDATE_PREFIX/gkeyll-results/parallel-c-4/$layer/creg-accepted"
  rm -rf "$dst"
  if [[ -d "$src" ]]; then mv "$src" "$dst"; fi
done
started="$(date +%s)"
"$candidate_gkeyll" runregression run -c --parallel --execute-only check
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > candidate-parallel-c-regression-check-seconds.txt
echo "Candidate parallel C-regression check runtime: $elapsed seconds"
"$candidate_gkeyll" ci/jenkins/check_regression_results.lua \
  "$CI_CANDIDATE_PREFIX/gkeyll-results/parallel-c-4" \
  ci/jenkins/expected_regression_diffs.txt ci-parallel-regression-summary.txt
