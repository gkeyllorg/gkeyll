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
  root="$(cache_root "$1")"
  printf '%s/%s/%s' "$root" "$2" "$3"
}

baseline_path() { platform_dir "$1" baseline-cache "$2"; }
candidate_path() { platform_dir "$1" candidate-cache "$2"; }

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
  rm -rf "$target"
  mkdir -p "$target/gkeyll"
  printf '%s' "$target"
}

publish_baseline() {
  local root="$1" platform="$2" sha="$3" stage="$4" parent target
  parent="$(baseline_path "$root" "$platform")"
  target="$parent/$sha"
  [[ -d "$stage" ]] || die "baseline directory is missing: $stage"
  [[ "$stage" == "$target" ]] || die "baseline must be built at its final path: $target"
  cache_tree_valid "$target" "$platform" "$sha" || die "baseline cache is incomplete or invalid: $target"
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
  rm -rf "$parent"
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
  for file in "$workspace"/ci-*.txt "$workspace"/*-unit-results.txt "$workspace"/*-seconds.txt "$workspace"/candidate-*.log "$workspace"/slurm-*.out "$workspace"/ci-report.md; do
    [[ -f "$file" && ! -L "$file" ]] || continue
    destination="$target/$(basename "$file")"
    rm -f "$destination"
    cp "$file" "$destination"
  done
  rm -f "$target/candidate-manifest.txt"
  printf 'candidate_commit=%s\nbaseline_commit=%s\nresult=%s\nterminal_stage=%s\nbuild_tag=%s\nfinished_utc=%s\n' \
    "$sha" "$baseline_sha" "$result" "$stage" "${BUILD_TAG:-manual}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$target/candidate-manifest.txt"
}

stage_candidate_artifacts() {
  local root="$1" platform="$2" sha="$3" workspace="$4" target subtree file relative
  [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || die 'candidate commit must be a full SHA'
  [[ -d "$workspace" ]] || die "Jenkins workspace is missing: $workspace"
  target="$(candidate_path "$root" "$platform")/$sha"
  [[ -d "$target" ]] || die "candidate directory is missing: $target"
  # The workspace came from the candidate checkout; recreate these archive
  # destinations so tracked symlinks cannot redirect diagnostic copies.
  rm -rf "$workspace/gkylsoft" "$workspace/build"
  for subtree in gkylsoft/gkeyll-results gkeyll/build; do
    [[ -d "$target/$subtree" ]] || continue
    while IFS= read -r -d '' file; do
      relative="${file#"$target"/}"
      relative="${relative#gkeyll/}"
      mkdir -p "$workspace/$(dirname "$relative")"
      cp "$file" "$workspace/$relative"
    done < <(find "$target/$subtree" -type f \( -name regressiondb -o -name _rr_failures.txt -o -name '*.log' \) -print0)
  done
}

case "${1:-}" in
  valid) manifest_valid "$2" "$3" "$4" ;;
  baseline-path) printf '%s\n' "$(baseline_path "$2" "$3")/$4" ;;
  prepare-baseline) prepare_baseline_stage "$2" "$3" "$4"; echo ;;
  publish-baseline) publish_baseline "$2" "$3" "$4" "$5"; echo ;;
  prepare-candidate) prepare_candidate "$2" "$3" "$4" "$5" "$6"; echo ;;
  finalize-candidate) finalize_candidate "$2" "$3" "$4" "$5" "$6" "$7" "$8" ;;
  stage-candidate-artifacts) stage_candidate_artifacts "$2" "$3" "$4" "$5" ;;
  *) die 'usage: baseline_cache.sh {valid|baseline-path|prepare-baseline|publish-baseline|prepare-candidate|finalize-candidate|stage-candidate-artifacts} ...' ;;
esac
