#!/usr/bin/env bash
# Focused acceptance test for candidate-versus-baseline regression acknowledgments.
# Usage: ci/jenkins/tests/test_check_regression_results.sh /path/to/gkeyll
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
  "$gkeyll" -S "$script_dir/../check_regression_results.lua" "$results_dir" \
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

sqlite3 "$results_dir/moments/regressiondb" <<'SQL'
UPDATE RegressionMeta SET run_mode = 'cpu_serial';
DELETE FROM RegressionData;
INSERT INTO RegressionData VALUES ('fixture', 'moments/creg/rt_ack_delta', 1, '');
INSERT INTO RegressionData VALUES ('fixture', 'moments/creg/rt_skip', -1, '');
INSERT INTO RegressionData VALUES ('older', 'moments/creg/rt_old', 1, '');
SQL
entries=$'rt_ack_delta # inherited\nmoments/rt_ack_delta cpu_serial # active\nmoments/creg/rt_ack_delta out.gkyl # file\nrt_ack_delta cpu_serial out.gkyl # scoped file\nrt_ack_delta gpu_serial # other mode\nrt_ack_delta cpu_parallel out.gkyl # other file mode\nrt_skip\nrt_unrun\nrt_old'
run_case passing "$entries" 'rt_ack_delta # inherited' 0 0 0 0
test "$(grep -c '^c_regression_listed_passing=' "$work_dir/passing-summary.txt")" -eq 4
grep -Fx 'c_regression_listed_passing=moments/creg/rt_ack_delta|inherited|rt_ack_delta # inherited' "$work_dir/passing-summary.txt"
grep -Fx 'c_regression_listed_passing=moments/creg/rt_ack_delta|new or updated|moments/rt_ack_delta cpu_serial # active' "$work_dir/passing-summary.txt"

# A failed comparison is not evidence that a listed test agrees.
sqlite3 "$results_dir/moments/regressiondb" "UPDATE RegressionData SET status=0 WHERE guid='fixture' AND name='moments/creg/rt_ack_delta';"
run_case failing "$entries" '' 0 1 0 0
test "$(grep -c '^c_regression_listed_passing=' "$work_dir/failing-summary.txt")" -eq 0

echo 'check_regression_results acknowledgment tests passed'
