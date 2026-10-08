#!/usr/bin/env bash
# Focused acceptance test for persistent CI cache lifecycle helpers.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cache_tool="$repo_root/ci/jenkins/baseline_cache.sh"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/gkeyll-baseline-cache.XXXXXX")"
trap 'rm -rf "$work_dir"' EXIT

root="$work_dir/root"
platform=personal

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
cat > "$stage/cache-manifest.txt" <<EOF
format_version=1
platform=$platform
baseline_commit=$sha
EOF

entry="$($cache_tool publish-baseline "$root" "$platform" "$sha" "$stage")"
$cache_tool valid "$root" "$platform" "$sha"
test "$entry" = "$stage"
"$entry/gkylsoft/gkeyll/bin/gkeyll"

printf 'outside manifest\n' > "$work_dir/outside-manifest.txt"
ln -s "$work_dir/outside-manifest.txt" "$fixture/candidate-manifest.txt"
candidate="$($cache_tool prepare-candidate "$root" "$platform" "$sha" "$fixture" "$sha")"
test "$candidate" = "$root/candidate-cache/$platform/$sha"
test "$(git -C "$candidate/gkeyll" rev-parse HEAD)" = "$sha"
test -f "$candidate/gkeyll/placeholder"
test "$(cat "$work_dir/outside-manifest.txt")" = 'outside manifest'
mkdir -p "$candidate/gkylsoft/gkeyll/bin" "$candidate/gkylsoft/gkeyll-results/moments/creg-runs" "$candidate/gkeyll/build"
touch "$candidate/dependency" "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"
printf '#!/usr/bin/env bash\ntest -f "%s/dependency"\n' "$candidate" > "$candidate/gkylsoft/gkeyll/bin/gkeyll"
chmod +x "$candidate/gkylsoft/gkeyll/bin/gkeyll"
"$candidate/gkylsoft/gkeyll/bin/gkeyll"
printf 'failure details\n' > "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/_rr_failures.txt"
printf 'build log\n' > "$candidate/gkeyll/build/compile.log"
printf 'outside data\n' > "$work_dir/outside.log"
ln -s "$work_dir/outside.log" "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/outside.log"
mkdir -p "$work_dir/workspace"
printf 'unit results\n' > "$work_dir/workspace/candidate-unit-results.txt"
rm "$candidate/candidate-manifest.txt"
ln -s "$work_dir/outside-manifest.txt" "$candidate/candidate-manifest.txt"
$cache_tool finalize-candidate "$root" "$platform" "$sha" "$sha" failure 'C regressions' "$work_dir/workspace"
"$candidate/gkylsoft/gkeyll/bin/gkeyll"
test "$(cat "$work_dir/outside-manifest.txt")" = 'outside manifest'
grep -q '^result=failure$' "$candidate/candidate-manifest.txt"
grep -q '^baseline_commit='"$sha"'$' "$candidate/candidate-manifest.txt"
test -f "$candidate/candidate-unit-results.txt"
ln -s "$work_dir" "$work_dir/workspace/gkylsoft"
ln -s "$work_dir" "$work_dir/workspace/build"
$cache_tool stage-candidate-artifacts "$root" "$platform" "$sha" "$work_dir/workspace"
test -f "$work_dir/workspace/gkylsoft/gkeyll-results/moments/creg-runs/_rr_failures.txt"
test -f "$work_dir/workspace/build/compile.log"
test ! -e "$work_dir/workspace/gkylsoft/gkeyll-results/moments/creg-runs/outside.log"
test -f "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"

test "$($cache_tool prepare-candidate "$root" "$platform" "$sha" "$fixture" "$sha")" = "$candidate"
test ! -e "$candidate/dependency"
test ! -e "$candidate/gkylsoft/gkeyll-results/moments/creg-runs/output.gkyl"
grep -q '^result=building$' "$candidate/candidate-manifest.txt"

touch "$fixture/another-file"
git -C "$fixture" add another-file
GIT_AUTHOR_DATE='2000-01-02T00:00:00Z' GIT_COMMITTER_DATE='2000-01-02T00:00:00Z' git -C "$fixture" commit -qm another
next_sha="$(git -C "$fixture" rev-parse HEAD)"
next_candidate="$($cache_tool prepare-candidate "$root" "$platform" "$next_sha" "$fixture" "$sha")"
test "$next_candidate" = "$root/candidate-cache/$platform/$next_sha"
test ! -e "$candidate"

rm "$entry/cache-manifest.txt"
if "$cache_tool" valid "$root" "$platform" "$sha"; then
  echo 'incomplete baseline cache was accepted as valid' >&2
  exit 1
fi
test "$($cache_tool prepare-baseline "$root" "$platform" "$sha")" = "$entry"
test ! -e "$entry/dependency"

echo 'baseline cache helper test passed'
