#!/usr/bin/env bash
# Run trailing-comma stripping and clang-format over every in-scope C/C++/CUDA
# file in the repo, in that order (a trailing comma before a closing '}'
# makes clang-format expand the whole list onto one line per element, so it
# has to be gone before clang-format sees the file).
#
# In scope: tracked *.c/*.h/*.cpp/*.hpp/*.cu/*.cuh files, excluding anything
# under a */ker/* directory (auto-generated DG kernel code) or under
# core/minus/ (vendored third-party libraries) -- see
# .agents/skills/code_style/SKILL.md.
#
# Usage:
#   ci/format-all.sh          # reformat all in-scope files in place
#   ci/format-all.sh --check  # fail (exit 1) if any in-scope file is not formatted
set -euo pipefail

PINNED_VERSION="18.1.8" # keep in sync with the `rev:` in .pre-commit-config.yaml

repo_root=$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)
cd "$repo_root"

mode=${1:-}

version=$(clang-format --version | grep -o '[0-9][0-9.]*' | head -1)
if [[ "$version" != "$PINNED_VERSION" ]]; then
  echo "warning: clang-format $version is on PATH, but pre-commit/CI enforce $PINNED_VERSION." >&2
  echo "         results may not exactly match the CI check. See README.md." >&2
fi

list_files() {
  git ls-files -z -- '*.c' '*.h' '*.cpp' '*.hpp' '*.cu' '*.cuh' \
    | grep -zv '/ker/' \
    | grep -zv '^core/minus/'
}

if [[ "$mode" == "--check" ]]; then
  status=0
  while IFS= read -r -d '' f; do
    if ! python3 ci/strip-trailing-commas.py --check "$f" > /dev/null 2>&1; then
      echo "trailing comma(s) found: $f"
      status=1
    fi
    if ! clang-format --dry-run --Werror "$f" > /dev/null 2>&1; then
      echo "not formatted: $f"
      status=1
    fi
  done < <(list_files)
  exit "$status"
else
  count=0
  while IFS= read -r -d '' f; do
    python3 ci/strip-trailing-commas.py "$f"
    clang-format -i "$f"
    count=$((count + 1))
  done < <(list_files)
  echo "Formatted $count files."
fi
