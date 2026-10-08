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
  if [[ -n "${GKEYLL_CI_PREBUILT_CONFIG:-}" ]]; then
    # A cached executable embeds its dependency paths. Rebuild it against this
    # run's private dependencies rather than reuse another run's installation.
    local run
    run="$(run_path "$1" "$2")" || return
    printf '%s/baseline' "$run"
  else
    platform_dir "$1" baseline-cache "$2"
  fi
}
candidate_path() {
  local run
  run="$(run_path "$1" "$2")" || return
  printf '%s/candidate' "$run"
}

cache_tree_valid() {
  local entry="$1" platform="$2" sha="$3" manifest head layer
  [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || return 1
  manifest="$entry/cache-manifest.txt"
  [[ -d "$entry/gkeyll/.git" && -x "$entry/gkylsoft/gkeyll/bin/gkeyll" && -f "$manifest" && -d "$entry/gkylsoft/gkeyll-results" ]] || return 1
  [[ "$(awk -F= '$1 == "format_version" {print $2}' "$manifest")" == 1 ]] || return 1
  [[ "$(awk -F= '$1 == "platform" {print $2}' "$manifest")" == "$platform" ]] || return 1
  [[ "$(awk -F= '$1 == "baseline_commit" {print $2}' "$manifest")" == "$sha" ]] || return 1
  head="$(git -C "$entry/gkeyll" rev-parse HEAD 2>/dev/null || true)"
  [[ "$head" == "$sha" ]] || return 1
  for layer in moments vlasov gyrokinetic pkpm; do
    [[ -d "$entry/gkylsoft/gkeyll-results/$layer/creg-accepted" ]] || return 1
    [[ -d "$entry/gkylsoft/gkeyll-results/parallel-c-4/$layer/creg-accepted" ]] || return 1
  done
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
  if [[ -n "${GKEYLL_CI_PREBUILT_CONFIG:-}" ]]; then
    [[ ! -e "$target" ]] || die "baseline already exists for this CI run: $target"
  else
    rm -rf "$target"
  fi
  mkdir -p "$target/gkeyll"
  printf '%s' "$target"
}

publish_baseline() {
  local root="$1" platform="$2" sha="$3" stage="$4" workspace="${5:-}" parent target file
  parent="$(baseline_path "$root" "$platform")"
  target="$parent/$sha"
  [[ -d "$stage" ]] || die "baseline directory is missing: $stage"
  [[ "$stage" == "$target" ]] || die "baseline must be built at its final path: $target"
  cache_tree_valid "$target" "$platform" "$sha" || die "baseline cache is incomplete or invalid: $target"
  if [[ -n "$workspace" ]]; then
    mkdir -p "$target/ci-command-logs"
    for file in "$workspace"/baseline-*.log "$workspace"/baseline-*-seconds.txt "$workspace"/baseline-unit-results.txt; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/"
    done
    for file in "$workspace"/ci-command-logs/baseline-*.log*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      cp "$file" "$target/ci-command-logs/"
    done
  fi
  find "$parent" -mindepth 1 -maxdepth 1 -type d ! -name "$sha" -exec rm -rf {} +
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
  for file in "$workspace"/ci-*.txt "$workspace"/ci-*-config.mak "$workspace"/*-unit-results.txt "$workspace"/*-seconds.txt "$workspace"/candidate-*.log "$workspace"/baseline-*.log "$workspace"/slurm-*.out "$workspace"/ci-report*.md "$workspace"/ci-report*.md.json; do
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

stage_candidate_artifacts() {
  local root="$1" platform="$2" sha="$3" workspace="$4" baseline_sha="${5:-}" target subtree file relative baseline
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
  rm -rf "$workspace/gkylsoft" "$workspace/build" "$workspace/cuda-build"
  for subtree in gkylsoft/gkeyll-results gkeyll/build gkeyll/cuda-build; do
    [[ -d "$target/$subtree" ]] || continue
    while IFS= read -r -d '' file; do
      relative="${file#"$target"/}"
      relative="${relative#gkeyll/}"
      mkdir -p "$workspace/$(dirname "$relative")"
      cp "$file" "$workspace/$relative"
    done < <(find "$target/$subtree" -type f \( -path '*/gkeyll-results/*' -o -name '*.log' \) -print0)
  done
  if [[ -n "$baseline_sha" ]]; then
    [[ "$baseline_sha" =~ ^[0-9a-f]{40}$ ]] || die 'baseline commit must be a full SHA'
    baseline="$(baseline_path "$root" "$platform")/$baseline_sha"
    if [[ -f "$baseline/gkeyll/config.mak" && ! -L "$baseline/gkeyll/config.mak" ]]; then
      rm -f "$workspace/ci-baseline-config.mak"
      cp "$baseline/gkeyll/config.mak" "$workspace/ci-baseline-config.mak"
    fi
    rm -rf "$workspace/_baseline"
    for subtree in gkylsoft/gkeyll-results gkeyll/build gkeyll/cuda-build; do
      [[ -d "$baseline/$subtree" ]] || continue
      while IFS= read -r -d '' file; do
        relative="${file#"$baseline"/}"
        relative="${relative#gkeyll/}"
        mkdir -p "$workspace/_baseline/$(dirname "$relative")"
        cp "$file" "$workspace/_baseline/$relative"
      done < <(find "$baseline/$subtree" -type f \( -path '*/gkeyll-results/*' -o -name '*.log' \) -print0)
    done
    # Restore completed baseline diagnostics on cache hits, without replacing
    # logs from a baseline built during this run.
    for file in "$baseline"/baseline-*.log "$baseline"/baseline-*-seconds.txt "$baseline"/baseline-unit-results.txt "$baseline"/ci-command-logs/baseline-*.log*; do
      [[ -f "$file" && ! -L "$file" ]] || continue
      relative="${file#"$baseline"/}"
      [[ ! -e "$workspace/$relative" && ! -L "$workspace/$relative" ]] || continue
      mkdir -p "$workspace/$(dirname "$relative")"
      cp "$file" "$workspace/$relative"
    done
  fi
  # Preserve the baseline snapshot with this run even after the shared cache
  # changes. The staging loops above copy regular files only.
  if [[ -d "$workspace/_baseline" && ! -L "$workspace/_baseline" ]]; then
    rm -rf "$target/_baseline"
    cp -a "$workspace/_baseline" "$target/_baseline"
  fi
}

case "${1:-}" in
  run-path) run_path "$2" "$3"; echo ;;
  valid) manifest_valid "$2" "$3" "$4" ;;
  baseline-path) printf '%s\n' "$(baseline_path "$2" "$3")/$4" ;;
  prepare-baseline) prepare_baseline_stage "$2" "$3" "$4"; echo ;;
  publish-baseline) publish_baseline "$2" "$3" "$4" "$5" "${6:-}"; echo ;;
  prepare-candidate) prepare_candidate "$2" "$3" "$4" "$5" "$6"; echo ;;
  finalize-candidate) finalize_candidate "$2" "$3" "$4" "$5" "$6" "$7" "$8" ;;
  stage-candidate-artifacts) stage_candidate_artifacts "$2" "$3" "$4" "$5" "${6:-}" ;;
  *) die 'usage: baseline_cache.sh {run-path|valid|baseline-path|prepare-baseline|publish-baseline|prepare-candidate|finalize-candidate|stage-candidate-artifacts} ...' ;;
esac
