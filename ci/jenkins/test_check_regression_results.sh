#!/usr/bin/env bash
# Focused acceptance test for candidate-versus-baseline regression acknowledgments.
# Usage: ci/jenkins/test_check_regression_results.sh /path/to/gkeyll
set -euo pipefail

gkeyll="${1:?Usage: $0 /path/to/gkeyll}"
test -x "$gkeyll"

readonly script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly work_dir="$(mktemp -d "${TMPDIR:-/tmp}/gkeyll-ci-ack-test.XXXXXX")"
trap 'rm -rf "$work_dir"' EXIT

results_dir="$work_dir/gkeyll-results"
baseline_results_dir="$work_dir/baseline-gkeyll-results"
mkdir -p "$results_dir/moments" "$baseline_results_dir/moments"
sqlite3 "$results_dir/moments/regressiondb" <<'SQL'
CREATE TABLE RegressionMeta (guid TEXT, run_mode TEXT);
INSERT INTO RegressionMeta VALUES ('fixture', 'check');
CREATE TABLE RegressionData (guid TEXT, name TEXT, status INTEGER, runlog TEXT);
INSERT INTO RegressionData VALUES ('fixture', 'moments/creg/rt_ack_delta', 0, '--- Comparison failures ---\nfixture difference');
SQL
sqlite3 "$baseline_results_dir/moments/regressiondb" <<'SQL'
CREATE TABLE RegressionMeta (guid TEXT, run_mode TEXT);
INSERT INTO RegressionMeta VALUES ('baseline', 'create');
CREATE TABLE RegressionData (guid TEXT, name TEXT, status INTEGER, runlog TEXT);
INSERT INTO RegressionData VALUES ('baseline', 'moments/creg/rt_ack_delta', -2, '');
SQL

run_case() {
  local name="$1" candidate_line="$2" baseline_line="$3" expected_status="$4"
  local expected_acked="$5" expected_unacked="$6" expected_candidate_only="$7"
  local candidate_file="$work_dir/${name}-candidate.txt"
  local baseline_file="$work_dir/${name}-baseline.txt"
  local summary_file="$work_dir/${name}-summary.txt"

  printf '%s\n' "$candidate_line" > "$candidate_file"
  printf '%s\n' "$baseline_line" > "$baseline_file"

  set +e
  "$gkeyll" "$script_dir/check_regression_results.lua" "$results_dir" \
    "$candidate_file" "$summary_file" "$baseline_file" "$baseline_results_dir" > "$work_dir/${name}.log" 2>&1
  local status=$?
  set -e
  test "$status" -eq "$expected_status"
  grep -Fx "c_regression_acknowledged=$expected_acked" "$summary_file"
  grep -Fx "c_regression_unacknowledged=$expected_unacked" "$summary_file"
  grep -Fx "c_regression_candidate_only=$expected_candidate_only" "$summary_file"
}

run_case added \
  'moments/creg/rt_ack_delta # intentional algorithm change' '' 0 1 0 0
run_case inherited \
  'moments/creg/rt_ack_delta # intentional algorithm change' \
  'moments/creg/rt_ack_delta # intentional algorithm change' 1 0 1 0
run_case removed '' \
  'moments/creg/rt_ack_delta # intentional algorithm change' 1 0 1 0
run_case refreshed \
  'moments/creg/rt_ack_delta # updated algorithm rationale' \
  'moments/creg/rt_ack_delta # intentional algorithm change' 0 1 0 0

sqlite3 "$results_dir/moments/regressiondb" <<'SQL'
DELETE FROM RegressionData;
INSERT INTO RegressionData VALUES ('fixture', 'moments/creg/rt_candidate_only', 0, '');
SQL
run_case candidate_only '' '' 0 0 0 1

echo 'check_regression_results acknowledgment tests passed'
