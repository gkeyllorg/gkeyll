#!/usr/bin/env bash
# Focused acceptance test for persistent CI cache lifecycle helpers.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cache_tool="$repo_root/ci/jenkins/baseline_cache.sh"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/gkeyll-baseline-cache.XXXXXX")"
trap 'rm -rf "$work_dir"' EXIT

root="$work_dir/root"
platform=personal
export BUILD_TAG=fixture-run-1
unset GKEYLL_CI_PREBUILT_CONFIG CI_BASELINE_CACHE_CONTEXT GKEYLL_CI_CACHE_REVISION

fixture="$work_dir/fixture"
git init -q "$fixture"
git -C "$fixture" config user.email ci@example.invalid
git -C "$fixture" config user.name CI
touch "$fixture/placeholder"
git -C "$fixture" add placeholder
GIT_AUTHOR_DATE='2000-01-01T00:00:00Z' GIT_COMMITTER_DATE='2000-01-01T00:00:00Z' git -C "$fixture" commit -qm fixture
sha="$(git -C "$fixture" rev-parse HEAD)"
if "$cache_tool" prepare-candidate / "$platform" "$sha" "$fixture" "$sha" >/dev/null 2>&1; then
  echo 'filesystem root was accepted as a cache root' >&2
  exit 1
fi
if "$cache_tool" prepare-candidate "$root" ../other "$sha" "$fixture" "$sha" >/dev/null 2>&1; then
  echo 'invalid platform was accepted' >&2
  exit 1
fi

stage="$($cache_tool prepare-baseline "$root" "$platform" "$sha")"
test "$stage" = "$root/baseline-cache/$platform/$sha"
test -d "$stage/gkeyll"
if "$cache_tool" valid "$root" "$platform" "$sha"; then
  echo 'unpublished baseline was accepted as valid' >&2
  exit 1
fi
git clone -q "$fixture" "$stage/gkeyll"
mkdir -p "$stage/gkylsoft/gkeyll/bin"
for layer in moments vlasov gyrokinetic pkpm; do
  mkdir -p "$stage/gkylsoft/gkeyll-results/$layer/creg-accepted"
  mkdir -p "$stage/gkylsoft/gkeyll-results/parallel-c-4/$layer/creg-accepted"
done
touch "$stage/dependency"
printf '#!/usr/bin/env bash\ntest -f "%s/dependency"\n' "$stage" > "$stage/gkylsoft/gkeyll/bin/gkeyll"
chmod +x "$stage/gkylsoft/gkeyll/bin/gkeyll"
# Real SQLite fixtures also exercise create-run validation at publication.
make_results() {
  python3 - "$1" <<'PYTHON'
from pathlib import Path
import sqlite3
import sys
root = Path(sys.argv[1]) / 'gkylsoft/gkeyll-results'
for results in (root, root / 'parallel-c-4'):
    results.mkdir(parents=True, exist_ok=True)
    (root / 'runregression.config.lua').write_text('return {}\n')
    for layer in ('moments', 'vlasov', 'gyrokinetic', 'pkpm'):
        accepted = results / layer / 'creg-accepted'
        accepted.mkdir(parents=True, exist_ok=True)
        (accepted / 'field.gkyl').write_bytes(b'accepted fixture')
        with sqlite3.connect(results / layer / 'regressiondb') as db:
            db.executescript("""create table RegressionMeta (guid text, ntotal integer, nfail integer);
                create table RegressionData (guid text, status integer);
                insert into RegressionMeta values ('create', 1, 0);
                insert into RegressionData values ('create', -2);""")
PYTHON
}
make_results "$stage"

mkdir -p "$work_dir/baseline-workspace/ci-command-logs" "$stage/gkeyll/cuda-build"
printf 'baseline warning\n' > "$work_dir/baseline-workspace/ci-command-logs/baseline-unit-build.log"
printf '0\n' > "$work_dir/baseline-workspace/ci-command-logs/baseline-unit-build.log.exit"
printf 'baseline compile warning\n' > "$work_dir/baseline-workspace/baseline-c-compile.log"
printf '0\n' > "$work_dir/baseline-workspace/baseline-c-compile.log.exit"
printf '0\n' > "$work_dir/baseline-workspace/baseline-install.log.exit"
printf 'baseline install warning\n' > "$work_dir/baseline-workspace/baseline-install.log"
printf '12\n' > "$work_dir/baseline-workspace/baseline-install-seconds.txt"
printf 'baseline CUDA log\n' > "$stage/gkeyll/cuda-build/compile.log"
printf '\000baseline field\377' > "$stage/gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl"
printf 'return {}\n' > "$stage/gkylsoft/gkeyll-results/runregression.config.lua"
printf 'PREFIX=baseline\n' > "$stage/gkeyll/config.mak"
entry="$($cache_tool publish-baseline "$root" "$platform" "$sha" "$stage" "$work_dir/baseline-workspace")"
$cache_tool valid "$root" "$platform" "$sha"
test "$entry" = "$stage"
"$entry/gkylsoft/gkeyll/bin/gkeyll"

printf 'outside manifest\n' > "$work_dir/outside-manifest.txt"
ln -s "$work_dir/outside-manifest.txt" "$fixture/candidate-manifest.txt"
candidate="$($cache_tool prepare-candidate "$root" "$platform" "$sha" "$fixture" "$sha")"
test "$candidate" = "$root/runs/$platform/$BUILD_TAG/candidate/$sha"
test "$(git -C "$candidate/gkeyll" rev-parse HEAD)" = "$sha"
test -f "$candidate/gkeyll/placeholder"
test "$(cat "$work_dir/outside-manifest.txt")" = 'outside manifest'
mkdir -p "$candidate/gkylsoft/gkeyll/bin" "$candidate/gkylsoft/gkeyll-results/moments/creg-runs" "$candidate/gkeyll/build"
printf 'PREFIX=candidate\n' > "$candidate/gkeyll/config.mak"
dependencies="$root/runs/$platform/$BUILD_TAG/dependencies"
mkdir -p "$dependencies"
printf '{"source_sha256":"fixture"}\n' > "$dependencies/manifest.json"
printf 'PREFIX=original\n' > "$dependencies/source-config.mak"
printf 'export MPI_HOME=fixture\n' > "$dependencies/env.sh"
touch "$candidate/dependency" "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"
printf '#!/usr/bin/env bash\ntest -f "%s/dependency"\n' "$candidate" > "$candidate/gkylsoft/gkeyll/bin/gkeyll"
chmod +x "$candidate/gkylsoft/gkeyll/bin/gkeyll"
"$candidate/gkylsoft/gkeyll/bin/gkeyll"
printf 'failure details\n' > "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/_rr_failures.txt"
printf 'build log\n' > "$candidate/gkeyll/build/compile.log"
printf 'outside data\n' > "$work_dir/outside.log"
ln -s "$work_dir/outside.log" "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/outside.log"
mkdir -p "$work_dir/workspace/ci-command-logs" "$candidate/gkeyll/cuda-build"
printf 'candidate CUDA log\n' > "$candidate/gkeyll/cuda-build/compile.log"
printf 'candidate error\n' > "$work_dir/workspace/ci-command-logs/candidate-unit-build.log"
printf '1\n' > "$work_dir/workspace/ci-command-logs/candidate-unit-build.log.exit"
printf 'report part\n' > "$work_dir/workspace/ci-report-2.md"
printf '{}\n' > "$work_dir/workspace/ci-report.md.json"
printf 'unit results\n' > "$work_dir/workspace/candidate-unit-results.txt"
mkdir -p "$work_dir/workspace/ci-valgrind/candidate"
printf 'full memory error log\n' > "$work_dir/workspace/ci-valgrind/candidate/ctest_fixture.log"
rm "$candidate/candidate-manifest.txt"
ln -s "$work_dir/outside-manifest.txt" "$candidate/candidate-manifest.txt"
$cache_tool finalize-candidate "$root" "$platform" "$sha" "$sha" failure 'C regressions' "$work_dir/workspace"
"$candidate/gkylsoft/gkeyll/bin/gkeyll"
test "$(cat "$work_dir/outside-manifest.txt")" = 'outside manifest'
grep -q '^result=failure$' "$candidate/candidate-manifest.txt"
grep -q '^baseline_commit='"$sha"'$' "$candidate/candidate-manifest.txt"
test -f "$candidate/candidate-unit-results.txt"
cmp "$work_dir/workspace/ci-command-logs/candidate-unit-build.log" "$candidate/ci-command-logs/candidate-unit-build.log"
test -f "$candidate/ci-command-logs/candidate-unit-build.log.exit"
test -f "$candidate/ci-report-2.md"
test -f "$candidate/ci-report.md.json"
cmp "$work_dir/workspace/ci-valgrind/candidate/ctest_fixture.log" "$candidate/ci-valgrind/candidate/ctest_fixture.log"
ln -s "$work_dir" "$work_dir/workspace/gkylsoft"
ln -s "$work_dir" "$work_dir/workspace/build"
$cache_tool stage-candidate-artifacts "$root" "$platform" "$sha" "$work_dir/workspace" "$sha"
test -f "$work_dir/workspace/gkylsoft/gkeyll-results/moments/creg-runs/_rr_failures.txt"
test -f "$work_dir/workspace/build/compile.log"
test -f "$work_dir/workspace/cuda-build/compile.log"
test -f "$work_dir/workspace/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"
cmp "$stage/gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl" "$work_dir/workspace/_baseline/gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl"
cmp "$stage/gkylsoft/gkeyll-results/runregression.config.lua" "$work_dir/workspace/_baseline/gkylsoft/gkeyll-results/runregression.config.lua"
cmp "$stage/gkeyll/config.mak" "$work_dir/workspace/ci-baseline-config.mak"
cmp "$candidate/gkeyll/config.mak" "$work_dir/workspace/ci-candidate-config.mak"
cmp "$dependencies/manifest.json" "$work_dir/workspace/ci-dependencies/manifest.json"
cmp "$dependencies/env.sh" "$work_dir/workspace/ci-dependencies/env.sh"
test -f "$work_dir/workspace/_baseline/cuda-build/compile.log"
cmp "$work_dir/baseline-workspace/ci-command-logs/baseline-unit-build.log" "$work_dir/workspace/ci-command-logs/baseline-unit-build.log"
test "$(cat "$work_dir/workspace/ci-command-logs/baseline-unit-build.log.exit")" = 0
test "$(cat "$work_dir/workspace/baseline-install-seconds.txt")" = 12
printf 'current run log\n' > "$work_dir/workspace/ci-command-logs/baseline-unit-build.log"
$cache_tool stage-candidate-artifacts "$root" "$platform" "$sha" "$work_dir/workspace" "$sha"
test "$(cat "$work_dir/workspace/ci-command-logs/baseline-unit-build.log")" = 'current run log'
test ! -e "$work_dir/workspace/gkylsoft/gkeyll-results/moments/creg-runs/outside.log"
test -f "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"

if "$cache_tool" prepare-candidate "$root" "$platform" "$sha" "$fixture" "$sha" >/dev/null 2>&1; then
  echo 'an existing CI run was overwritten' >&2
  exit 1
fi
export BUILD_TAG=fixture-run-2
repeated="$($cache_tool prepare-candidate "$root" "$platform" "$sha" "$fixture" "$sha")"
test "$repeated" != "$candidate"
test -e "$candidate/dependency"
test -e "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"
grep -q '^result=building$' "$repeated/candidate-manifest.txt"

touch "$fixture/another-file"
git -C "$fixture" add another-file
GIT_AUTHOR_DATE='2000-01-02T00:00:00Z' GIT_COMMITTER_DATE='2000-01-02T00:00:00Z' git -C "$fixture" commit -qm another
next_sha="$(git -C "$fixture" rev-parse HEAD)"
export BUILD_TAG=fixture-run-3
next_candidate="$($cache_tool prepare-candidate "$root" "$platform" "$next_sha" "$fixture" "$sha")"
test "$next_candidate" = "$root/runs/$platform/$BUILD_TAG/candidate/$next_sha"
test -d "$candidate"
test -d "$repeated"

rm "$entry/cache-manifest.txt"
if "$cache_tool" valid "$root" "$platform" "$sha"; then
  echo 'incomplete baseline cache was accepted as valid' >&2
  exit 1
fi
test "$($cache_tool prepare-baseline "$root" "$platform" "$sha")" = "$entry"
test ! -e "$entry/dependency"
test -f "$candidate/_baseline/cuda-build/compile.log"
test -f "$candidate/_baseline/gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl"

# Prebuilt baselines use the same persistent path, independent of BUILD_TAG.
export GKEYLL_CI_PREBUILT_CONFIG="$work_dir/prebuilt-config.mak"
printf 'PREFIX=/fixture\n' > "$GKEYLL_CI_PREBUILT_CONFIG"
$cache_tool acquire "$root" "$platform"
if BUILD_TAG=other-run "$cache_tool" acquire "$root" "$platform" 2>/dev/null; then
  echo 'overlapping cache user was allowed' >&2; exit 1
fi
if BUILD_TAG=other-run "$cache_tool" release "$root" "$platform" 2>/dev/null; then
  echo 'another run released the cache lease' >&2; exit 1
fi

populate_baseline() {
  local commit="$1" target
  target="$($cache_tool prepare-baseline "$root" "$platform" "$commit")"
  git clone -q "$fixture" "$target/gkeyll"
  git -C "$target/gkeyll" checkout -q "$commit"
  mkdir -p "$target/gkylsoft/gkeyll/bin" "$target/dependencies/lib"
  printf 'PREFIX=baseline\n' > "$target/gkeyll/config.mak"
  printf 'dependency contents\n' > "$target/dependencies/lib/library.so"
  printf 'export BASELINE_LIBRARY=%q\n' "$target/dependencies/lib/library.so" > "$target/dependencies/env.sh"
  printf '{}\n' > "$target/dependencies/manifest.json"
  printf '#!/bin/bash\ntest -f "$BASELINE_LIBRARY"\n' > "$target/gkylsoft/gkeyll/bin/gkeyll"
  chmod +x "$target/gkylsoft/gkeyll/bin/gkeyll"
  make_results "$target"
}
populate_baseline "$sha"
entry="$($cache_tool publish-baseline "$root" "$platform" "$sha" "$entry" "$work_dir/baseline-workspace")"
$cache_tool release "$root" "$platform"

for run in fixture-run-4 fixture-run-5; do
  export BUILD_TAG="$run"
  $cache_tool acquire "$root" "$platform"
  $cache_tool valid "$root" "$platform" "$sha"
  test "$($cache_tool baseline-path "$root" "$platform" "$sha")" = "$entry"
  # No checkout, rebuild or baseline execution: restore directly into a new run.
  mkdir -p "$work_dir/$run"
  $cache_tool restore-diagnostics "$entry" "$work_dir/$run"
  cmp "$work_dir/baseline-workspace/baseline-c-compile.log.exit" "$work_dir/$run/baseline-c-compile.log.exit"
  cmp "$work_dir/baseline-workspace/baseline-install.log" "$work_dir/$run/baseline-install.log"
  ( . "$entry/dependencies/env.sh"; "$entry/gkylsoft/gkeyll/bin/gkeyll" )
  $cache_tool release "$root" "$platform"
done

# The trusted helper itself can be fetched to a different temporary directory.
mkdir "$work_dir/other-tools"
cp "$cache_tool" "$repo_root/ci/jenkins/prebuilt_config.py" "$repo_root/ci/jenkins/valgrind.py" "$work_dir/other-tools/"
test "$($cache_tool context)" = "$("$work_dir/other-tools/baseline_cache.sh" context)"

for file in baseline-install.log baseline-c-compile.log.exit gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl dependencies/lib/library.so; do
  mv "$entry/$file" "$work_dir/missing"
  if $cache_tool valid "$root" "$platform" "$sha" 2>/dev/null; then
    echo "missing $file was accepted" >&2; exit 1
  fi
  mv "$work_dir/missing" "$entry/$file"
done
printf 'corrupt output\n' >> "$entry/gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl"
if $cache_tool valid "$root" "$platform" "$sha"; then
  echo 'corrupt accepted output was accepted' >&2; exit 1
fi
populate_baseline "$sha"
$cache_tool publish-baseline "$root" "$platform" "$sha" "$entry" "$work_dir/baseline-workspace" >/dev/null
printf '# new settings\n' >> "$GKEYLL_CI_PREBUILT_CONFIG"
if $cache_tool valid "$root" "$platform" "$sha"; then
  echo 'changed prebuilt config was accepted' >&2; exit 1
fi
sed -i '$d' "$GKEYLL_CI_PREBUILT_CONFIG"
if GKEYLL_CI_CACHE_REVISION=changed $cache_tool valid "$root" "$platform" "$sha"; then
  echo 'changed toolchain revision was accepted' >&2; exit 1
fi
$cache_tool valid "$root" "$platform" "$sha"

# Failed creation never publishes or evicts the previous complete baseline.
populate_baseline "$next_sha"
next_entry="$($cache_tool baseline-path "$root" "$platform" "$next_sha")"
python3 - "$next_entry/gkylsoft/gkeyll-results/moments/regressiondb" <<'PYTHON'
import sqlite3, sys
with sqlite3.connect(sys.argv[1]) as db:
    db.execute('update RegressionData set status=-3')
PYTHON
if $cache_tool publish-baseline "$root" "$platform" "$next_sha" "$next_entry" "$work_dir/baseline-workspace" >/dev/null 2>&1; then
  echo 'failed create run was published' >&2; exit 1
fi
test ! -e "$next_entry/cache-manifest.txt"
$cache_tool valid "$root" "$platform" "$sha"
populate_baseline "$next_sha"
$cache_tool publish-baseline "$root" "$platform" "$next_sha" "$next_entry" "$work_dir/baseline-workspace" >/dev/null
$cache_tool valid "$root" "$platform" "$next_sha"
test ! -e "$entry"
# Results retained with an earlier candidate survive replacement of the cache.
test -f "$candidate/_baseline/gkylsoft/gkeyll-results/moments/creg-accepted/field.gkyl"

echo 'baseline cache helper test passed'
