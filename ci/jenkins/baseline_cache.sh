#!/usr/bin/env bash
# Helpers for the persistent Jenkins baseline and candidate build trees.
set -euo pipefail

die() { echo "baseline-cache: $*" >&2; exit 2; }

cache_root() {
  local root="$1"
  while [[ "$root" == */ ]]; do root="${root%/}"; done
  [[ -n "$root" && "$root" = /* && ! "$root" =~ (^|/)\.\.(/|$) ]] || die 'GKEYLL_CI_ROOT must be an absolute directory other than /'
  printf '%s' "$root"
}

platform_dir() {
  local root
  case "$3" in
    personal|team-workstation|stellar-cpu|perlmutter-gpu) ;;
    *) die "invalid CI platform: $3" ;;
  esac
  root="$(cache_root "$1")" || return
  printf '%s/%s/%s' "$root" "$2" "$3"
}

run_path() {
  local tag="${BUILD_TAG:?BUILD_TAG must identify this CI run}" parent
  [[ "$tag" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && "$tag" != *..* ]] || die 'invalid BUILD_TAG'
  parent="$(platform_dir "$1" runs "$2")" || return
  printf '%s/%s' "$parent" "$tag"
}

baseline_path() {
  platform_dir "$1" baseline-cache "$2"
}

# Capture this before build helpers change MPI/runtime variables. An explicit
# revision also lets operators invalidate after an in-place toolchain upgrade.
cache_context() {
  local name
  {
    printf 'cache-format=2\n'
    for name in GKEYLL_CI_CACHE_REVISION GKEYLL_CI_PREBUILT_CONFIG PERSONAL_MKDEPS_SCRIPT PERSONAL_CONFIGURE_SCRIPT PERSONAL_MPIEXEC PERSONAL_MPI_HOME CC CXX CUDA_HOME; do
      printf '%s=%s\n' "$name" "${!name:-}"
    done
    if [[ -n "${GKEYLL_CI_PREBUILT_CONFIG:-}" ]]; then
      sha256sum < "$GKEYLL_CI_PREBUILT_CONFIG" || return
    fi
    sha256sum < "${BASH_SOURCE[0]}" || return
    local helper="$(dirname "${BASH_SOURCE[0]}")/prebuilt_config.py"
    if [[ -f "$helper" ]]; then sha256sum < "$helper" || return; fi
  } | sha256sum | cut -d' ' -f1
}

# Hold this lease until results/report artifacts have been copied. Never delete
# a tree that another job is building or reading (including multibranch jobs).
acquire_cache() {
  local parent
  parent="$(baseline_path "$1" "$2")"
  mkdir -p "$parent"
  : "${BUILD_TAG:?}"
  if ! mkdir "$parent/.lock" 2>/dev/null; then
    [[ -d "$parent/.lock" ]] || die "cannot acquire cache lease: $parent/.lock"
    echo "baseline-cache: waiting for $parent/.lock (owner: $(cat "$parent/.lock/owner" 2>/dev/null || true))" >&2
    return 75
  fi
  printf '%s\n' "$BUILD_TAG" > "$parent/.lock/owner"
}

release_cache() {
  local parent
  parent="$(baseline_path "$1" "$2")"
  [[ "$(cat "$parent/.lock/owner" 2>/dev/null)" == "${BUILD_TAG:?}" ]] || die 'cache lease owner does not match this run'
  rm "$parent/.lock/owner"
  rmdir "$parent/.lock"
}

candidate_path() {
  local run
  run="$(run_path "$1" "$2")" || return
  printf '%s/candidate' "$run"
}

cache_payload_valid() {
  local entry="$1" sha="$2" layer results
  [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || return 1
  [[ -d "$entry/gkeyll/.git" && -x "$entry/gkylsoft/gkeyll/bin/gkeyll" && -s "$entry/gkeyll/config.mak" ]] || return 1
  [[ "$(git -C "$entry/gkeyll" rev-parse HEAD 2>/dev/null)" == "$sha" ]] || return 1
  [[ -s "$entry/gkylsoft/gkeyll-results/runregression.config.lua" ]] || return 1
  for results in "$entry/gkylsoft/gkeyll-results" "$entry/gkylsoft/gkeyll-results/parallel-c-4"; do
    for layer in moments vlasov gyrokinetic pkpm; do
      [[ -d "$results/$layer/creg-accepted" && -s "$results/$layer/regressiondb" ]] || return 1
    done
  done
  if [[ -n "${GKEYLL_CI_PREBUILT_CONFIG:-}" ]]; then
    [[ -s "$entry/dependencies/env.sh" && -s "$entry/dependencies/manifest.json" ]] || return 1
  fi
}

cache_tree_valid() {
  local entry="$1" platform="$2" sha="$3" manifest context
  manifest="$entry/cache-manifest.txt"
  [[ -f "$manifest" && -s "$entry/cache-checksums.txt" ]] || return 1
  [[ "$(awk -F= '$1 == "format_version" {print $2}' "$manifest")" == 2 ]] || return 1
  [[ "$(awk -F= '$1 == "platform" {print $2}' "$manifest")" == "$platform" ]] || return 1
  [[ "$(awk -F= '$1 == "baseline_commit" {print $2}' "$manifest")" == "$sha" ]] || return 1
  context="${CI_BASELINE_CACHE_CONTEXT:-$(cache_context)}" || return
  [[ "$(awk -F= '$1 == "context" {print $2}' "$manifest")" == "$context" ]] || return 1
  cache_payload_valid "$entry" "$sha" || return 1
  (cd "$entry" && sha256sum --check --status cache-checksums.txt) || return 1
}

manifest_valid() {
  local root="$1" platform="$2" sha="$3"
  cache_tree_valid "$(baseline_path "$root" "$platform")/$sha" "$platform" "$sha"
}

prepare_baseline_stage() {
  local root platform sha parent target
  root="$1"; platform="$2"; sha="$3"
  [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || die 'baseline commit must be a full SHA'
  parent="$(baseline_path "$root" "$platform")"
  mkdir -p "$parent"
  target="$parent/$sha"
  # Installed binaries and libraries can embed this absolute path. Keep it fixed
  # from the build through every later cache hit.
  rm -rf "$target"
  mkdir -p "$target/gkeyll"
  printf '%s' "$target"
}

publish_baseline() {
  local root="$1" platform="$2" sha="$3" stage="$4" workspace="${5:-}" parent target file
  [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || die 'baseline commit must be a full SHA'
  parent="$(baseline_path "$root" "$platform")"
  target="$parent/$sha"
  [[ -d "$stage" ]] || die "baseline directory is missing: $stage"
  [[ "$stage" == "$target" ]] || die "baseline must be built at its final path: $target"
  rm -f "$target/cache-manifest.txt"
  cache_payload_valid "$target" "$sha" || die "baseline cache is incomplete or invalid: $target"
  [[ -d "$workspace" ]] || die 'baseline diagnostics workspace is required'
  # runregression records failures in SQLite and may still exit successfully.
  # Never turn a failed or interrupted create run into a reusable reference.
  python3 - "$target/gkylsoft/gkeyll-results" <<'PYTHON'
from pathlib import Path
import sqlite3
import sys

root = Path(sys.argv[1])
for results in (root, root / 'parallel-c-4'):
    created = 0
    for layer in ('moments', 'vlasov', 'gyrokinetic', 'pkpm'):
        path = results / layer / 'regressiondb'
        with sqlite3.connect(path.as_uri() + '?mode=ro', uri=True) as db:
            meta = db.execute('select guid, ntotal, nfail from RegressionMeta order by rowid desc limit 1').fetchone()
            if meta is None:  # Some layers have no selected parallel tests.
                if db.execute('select count(*) from RegressionData').fetchone()[0]:
                    sys.exit(f'baseline-cache: unfinished regression run in {path}')
                continue
            guid, total, failed = meta
            rows = db.execute('select status from RegressionData where guid=?', (guid,)).fetchall()
            if failed or len(rows) != total or any(status not in (-2, -1) for status, in rows):
                sys.exit(f'baseline-cache: incomplete or failed baseline creation in {path}')
            created += sum(status == -2 for status, in rows)
    if not created:
        sys.exit(f'baseline-cache: no accepted tests created in {results}')
PYTHON
  if [[ -n "$workspace" ]]; then
    mkdir -p "$target/ci-command-logs"
    for file in "$workspace"/baseline-*.log* "$workspace"/baseline-*-seconds.txt "$workspace"/baseline-unit-results.txt; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/"
    done
    for file in "$workspace"/ci-command-logs/baseline-*.log*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/ci-command-logs/"
    done
  fi
  if [[ "$platform" == stellar-cpu || "$platform" == perlmutter-gpu ]]; then
    for file in "$workspace"/slurm-regression-* "$workspace"/slurm-parallel-regression-*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/baseline-$(basename "$file")"
    done
    for file in "$workspace"/ci-command-logs/shared-slurm-*regression.log*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/ci-command-logs/baseline-${file##*/shared-}"
    done
  fi
  # Require the compiler diagnostics used to distinguish introduced warnings.
  for file in unit-build install c-compile; do
    local log="$target/baseline-$file.log"
    [[ -f "$log" ]] || log="$target/ci-command-logs/baseline-$file.log"
    [[ -f "$log" && "$(cat "$log.exit" 2>/dev/null)" == 0 ]] || die "missing or failed baseline $file diagnostics"
  done
  (
    cd "$target"
    # Include all saved diagnostics and numerical data; detect partial copies,
    # deletion, or corruption before trusting a hit. Build objects stay in place.
    find gkylsoft ci-command-logs -type f -print0
    find . -maxdepth 1 -type f -name 'baseline-*' -print0
    printf '%s\0' gkeyll/config.mak
    if [[ -f gkeyll/ci/jenkins/expected_regression_diffs.txt ]]; then
      printf '%s\0' gkeyll/ci/jenkins/expected_regression_diffs.txt
    fi
    if [[ -d dependencies ]]; then find dependencies -type f -print0; fi
  ) | (cd "$target" && sort -z | xargs -0 sha256sum) > "$target/cache-checksums.txt"
  # Only this final rename publishes a reusable cache. Interrupted builds,
  # failed diagnostic copies and old-format entries always remain misses.
  printf 'format_version=2\nplatform=%s\nbaseline_commit=%s\ncontext=%s\ncreated_utc=%s\n' \
    "$platform" "$sha" "${CI_BASELINE_CACHE_CONTEXT:-$(cache_context)}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$target/cache-manifest.txt.tmp"
  mv "$target/cache-manifest.txt.tmp" "$target/cache-manifest.txt"
  find "$parent" -mindepth 1 -maxdepth 1 -type d ! -name "$sha" ! -name '.lock' -exec rm -rf {} +
  printf '%s' "$target"
}

prepare_candidate() {
  local root="$1" platform="$2" sha="$3" source="$4" baseline_sha="$5" parent target head
  [[ "$sha" =~ ^[0-9a-f]{40}$ && "$baseline_sha" =~ ^[0-9a-f]{40}$ ]] || die 'candidate and baseline commits must be full SHAs'
  [[ -d "$source/.git" ]] || die "candidate checkout is missing: $source"
  head="$(git -C "$source" rev-parse HEAD 2>/dev/null || true)"
  [[ "$head" == "$sha" ]] || die "candidate checkout does not match $sha"
  parent="$(candidate_path "$root" "$platform")"
  target="$parent/$sha"
  [[ ! -e "$target" ]] || die "candidate already exists for this CI run: $target"
  mkdir -p "$target/gkeyll"
  # Copy source before building; installed libraries must stay at this path.
  cp -a "$source/." "$target/gkeyll/"
  [[ "$(git -C "$target/gkeyll" rev-parse HEAD 2>/dev/null || true)" == "$sha" ]] || die "candidate copy does not match $sha"
  rm -f "$target/candidate-manifest.txt"
  printf 'candidate_commit=%s\nbaseline_commit=%s\nresult=building\nbuild_tag=%s\n' "$sha" "$baseline_sha" "${BUILD_TAG:-manual}" > "$target/candidate-manifest.txt"
  printf '%s' "$target"
}

finalize_candidate() {
  local root="$1" platform="$2" sha="$3" baseline_sha="$4" result="$5" stage="$6" workspace="$7" target file destination
  [[ "$sha" =~ ^[0-9a-f]{40}$ && "$baseline_sha" =~ ^[0-9a-f]{40}$ ]] || die 'candidate and baseline commits must be full SHAs'
  target="$(candidate_path "$root" "$platform")/$sha"
  [[ -d "$target/gkeyll/.git" ]] || die "candidate checkout is missing: $target/gkeyll"
  [[ -d "$workspace" ]] || die "Jenkins workspace is missing: $workspace"
  [[ "$(git -C "$target/gkeyll" rev-parse HEAD 2>/dev/null || true)" == "$sha" ]] || die "candidate checkout does not match $sha"
  for file in "$workspace"/ci-*.txt "$workspace"/ci-*-config.mak "$workspace"/*-unit-results.txt "$workspace"/*-seconds.txt "$workspace"/candidate-*.log* "$workspace"/baseline-*.log* "$workspace"/baseline-slurm-* "$workspace"/slurm-*.out "$workspace"/ci-report*.md "$workspace"/ci-report*.md.json; do
    [[ -f "$file" && ! -L "$file" ]] || continue
    destination="$target/$(basename "$file")"
    rm -f "$destination"
    cp "$file" "$destination"
  done
  rm -rf "$target/ci-command-logs"
  if [[ -d "$workspace/ci-command-logs" && ! -L "$workspace/ci-command-logs" ]]; then
    mkdir -p "$target/ci-command-logs"
    for file in "$workspace"/ci-command-logs/*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/ci-command-logs/"
    done
  fi
  rm -f "$target/candidate-manifest.txt"
  printf 'candidate_commit=%s\nbaseline_commit=%s\nresult=%s\nterminal_stage=%s\nbuild_tag=%s\nfinished_utc=%s\n' \
    "$sha" "$baseline_sha" "$result" "$stage" "${BUILD_TAG:-manual}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$target/candidate-manifest.txt"
}

restore_diagnostics() {
  local baseline="$1" workspace="$2" file relative
    for file in "$baseline"/baseline-*.log* "$baseline"/baseline-*.txt "$baseline"/baseline-*.out "$baseline"/ci-command-logs/baseline-*.log*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      relative="${file#"$baseline"/}"
      [[ ! -e "$workspace/$relative" && ! -L "$workspace/$relative" ]] || continue
      mkdir -p "$workspace/$(dirname "$relative")"
      cp "$file" "$workspace/$relative"
    done
}

# Stream regular files in bulk, without following symlinks or starting a shell
# command for each file. NUL-delimited names preserve spaces and newlines.
copy_artifact_files() {
  local source="$1" destination="$2" selection="$3"
  [[ -d "$source" && ! -L "$source" ]] || return 0
  mkdir -p "$destination"
  (
    cd "$source"
    case "$selection" in
      all) find . -type f -print0 ;;
      diagnostics) find . -type f ! -name '*.gkyl' -print0 ;;
      logs) find . -type f -name '*.log' -print0 ;;
      *) die "invalid artifact selection: $selection" ;;
    esac | tar -c -f - --null -T -
  ) | tar -x -f - -C "$destination"
}

pack_numerical_outputs() {
  local source="$1" destination="$2" extension started=$SECONDS
  local -a compressor
  [[ -d "$source/gkylsoft/gkeyll-results" && ! -L "$source/gkylsoft" && ! -L "$source/gkylsoft/gkeyll-results" ]] || return 0
  # Bound compression workers on shared agents. Gzip keeps agents without zstd
  # usable, including personal macOS installations.
  if command -v zstd >/dev/null 2>&1; then
    compressor=(zstd -q -1 -T2 -c)
    extension=zst
  else
    compressor=(gzip -1 -c)
    extension=gz
  fi
  destination="$destination.tar.$extension"
  # Publish only a complete archive. Jenkins excludes the temporary filename.
  if (
    cd "$source" || exit 1
    find gkylsoft/gkeyll-results -type f -name '*.gkyl' -print0 |
      tar -c -f - --null -T - | "${compressor[@]}" > "$destination.tmp"
  ); then
    mv "$destination.tmp" "$destination"
  else
    rm -f "$destination.tmp"
    die "could not package numerical outputs from $source"
  fi
  echo "Packaged $destination: $(wc -c < "$destination") bytes in $((SECONDS - started)) s"
}

stage_candidate_artifacts() {
  local root="$1" platform="$2" sha="$3" workspace="$4" baseline_sha="${5:-}" target subtree file baseline started=$SECONDS
  [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || die 'candidate commit must be a full SHA'
  [[ -d "$workspace" ]] || die "Jenkins workspace is missing: $workspace"
  target="$(candidate_path "$root" "$platform")/$sha"
  [[ -d "$target" ]] || die "candidate directory is missing: $target"
  if [[ -f "$target/gkeyll/config.mak" && ! -L "$target/gkeyll/config.mak" ]]; then
    rm -f "$workspace/ci-candidate-config.mak"
    cp "$target/gkeyll/config.mak" "$workspace/ci-candidate-config.mak"
  fi
  local dependencies="$(run_path "$root" "$platform")/dependencies"
  rm -rf "$workspace/ci-dependencies"
  if [[ -d "$dependencies" ]]; then
    mkdir -p "$workspace/ci-dependencies"
    for file in manifest.json source-config.mak env.sh; do
      [[ -f "$dependencies/$file" && ! -L "$dependencies/$file" ]] || continue
      cp "$dependencies/$file" "$workspace/ci-dependencies/$file"
    done
  fi
  # The workspace came from the candidate checkout; recreate these archive
  # destinations so tracked symlinks cannot redirect diagnostic copies.
  rm -rf "$workspace/gkylsoft" "$workspace/build" "$workspace/cuda-build" "$workspace/_baseline" "$workspace/ci-numerical"
  mkdir -p "$workspace/ci-numerical"
  if [[ ! -L "$target/gkylsoft" ]]; then
    copy_artifact_files "$target/gkylsoft/gkeyll-results" "$workspace/gkylsoft/gkeyll-results" diagnostics
  fi
  for subtree in build cuda-build; do
    copy_artifact_files "$target/gkeyll/$subtree" "$workspace/$subtree" logs
  done
  if [[ -n "$baseline_sha" ]]; then
    [[ "$baseline_sha" =~ ^[0-9a-f]{40}$ ]] || die 'baseline commit must be a full SHA'
    baseline="$(baseline_path "$root" "$platform")/$baseline_sha"
    if [[ -f "$baseline/gkeyll/config.mak" && ! -L "$baseline/gkeyll/config.mak" ]]; then
      rm -f "$workspace/ci-baseline-config.mak"
      cp "$baseline/gkeyll/config.mak" "$workspace/ci-baseline-config.mak"
    fi
    # Preserve a complete independent snapshot before the shared cache changes.
    # Numerical data is copied once, directly into the retained run tree.
    rm -rf "$target/_baseline"
    if [[ ! -L "$baseline/gkylsoft" ]]; then
      copy_artifact_files "$baseline/gkylsoft/gkeyll-results" "$target/_baseline/gkylsoft/gkeyll-results" all
    fi
    copy_artifact_files "$target/_baseline/gkylsoft/gkeyll-results" "$workspace/_baseline/gkylsoft/gkeyll-results" diagnostics
    for subtree in build cuda-build; do
      copy_artifact_files "$baseline/gkeyll/$subtree" "$target/_baseline/$subtree" logs
      copy_artifact_files "$target/_baseline/$subtree" "$workspace/_baseline/$subtree" logs
    done
    # Restore completed baseline diagnostics on cache hits, without replacing
    # logs from a baseline built during this run.
    restore_diagnostics "$baseline" "$workspace"
  fi
  echo "Staged diagnostics and retained baseline snapshot in $((SECONDS - started)) s"
  # Package directly from retained trees: never stage thousands of loose .gkyl
  # files only to have Jenkins transfer and recreate them on the controller.
  pack_numerical_outputs "$target" "$workspace/ci-numerical/candidate"
  if [[ -n "$baseline_sha" ]]; then
    pack_numerical_outputs "$target/_baseline" "$workspace/ci-numerical/baseline"
  fi
}

case "${1:-}" in
  context) cache_context ;;
  acquire) acquire_cache "$2" "$3" ;;
  release) release_cache "$2" "$3" ;;
  restore-diagnostics) restore_diagnostics "$2" "$3" ;;
  run-path) run_path "$2" "$3"; echo ;;
  valid) manifest_valid "$2" "$3" "$4" ;;
  baseline-path) printf '%s\n' "$(baseline_path "$2" "$3")/$4" ;;
  prepare-baseline) prepare_baseline_stage "$2" "$3" "$4"; echo ;;
  publish-baseline) publish_baseline "$2" "$3" "$4" "$5" "${6:-}"; echo ;;
  prepare-candidate) prepare_candidate "$2" "$3" "$4" "$5" "$6"; echo ;;
  finalize-candidate) finalize_candidate "$2" "$3" "$4" "$5" "$6" "$7" "$8" ;;
  stage-candidate-artifacts) stage_candidate_artifacts "$2" "$3" "$4" "$5" "${6:-}" ;;
  *) die 'usage: baseline_cache.sh {context|acquire|release|restore-diagnostics|run-path|valid|baseline-path|prepare-baseline|publish-baseline|prepare-candidate|finalize-candidate|stage-candidate-artifacts} ...' ;;
esac
