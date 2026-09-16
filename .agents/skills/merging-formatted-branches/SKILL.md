---
name: merging-formatted-branches
description: Resolve Gkeyll merges involving broad clang-format changes while preserving branch functionality, using normalized three-way comparisons.
---

# Merging formatted branches

Use when merging `agent_tools`, a branch containing it, or another broad formatting
change. Compare **base, ours, and theirs after applying the same formatting**.
Whitespace-ignore options alone miss inserted braces, removed trailing commas,
renamed tests, and moved code. Avoid blanket `ours`/`theirs` resolutions.

## Prepare comparisons

Record the original tips and protect existing local changes. Start the requested
merge with `git merge --no-commit <incoming>`; conflicts are expected.

Use the version pinned by `ci/format-all.sh` and `.pre-commit-config.yaml`, currently
**18.1.8**. A different installed clang-format can produce misleading differences.
If necessary, install the pinned version in a temporary Python venv.

Run from the repository root:

```sh
python3 .agents/skills/merging-formatted-branches/scripts/prepare_merge.py \
  --theirs <incoming> --clang-format /path/to/pinned/clang-format
```

The helper writes to a new temporary directory and leaves the index and working
tree untouched. It considers the C/C++/CUDA paths changed on **ours** since the
merge base, including cleanly merged paths and new files. It excludes `*/ker/*`
and `core/minus/*`. Symlinks and other file types still need ordinary merge review.

Read its `manifest.json`, then only the relevant `incoming-changes/*.diff` and
`merged/` files. It provides normalized `base/`, `ours/`, and `theirs/` copies;
strips trailing commas before formatting; formats to a stable result; and runs
`git merge-file --diff3`. `--ours` defaults to `HEAD`; `--theirs` defaults to
`MERGE_HEAD`. `--base` handles an explicitly chosen merge base. Formatting rules
come from `--style-ref` (default: theirs); choose a revision with the intended
formatter configuration if the branches use different rules.

## Resolve and stage

- If normalized base equals normalized theirs, the incoming edit is formatting
  only: the normalized ours version preserves the current branch's changes.
- Inspect remaining conflicts for actual edits. Keep compatible changes from
  both parents; use the user's requested branch for incompatible overlap.
  Preserve moved/refactored code and expanded tests when the other side merely
  formats their predecessors. Apply compatible API/test renames to the retained
  implementation and its callers.
- A delete/modify conflict can retain the deletion when the other side only
  formatted the obsolete file. Check its replacement and callers. Review
  add/add conflicts, real edits to deleted files, and renames individually.
- Copy reviewed candidates from `merged/` to their repository paths. Handle
  deletions explicitly; the helper does not apply them. Run trailing-comma
  stripping and the pinned formatter on the retained paths, then stage only
  reviewed paths. Include branch-added source files in formatting checks.

Compare the result with normalized ours to verify that additional differences
are intended incoming changes. Check for unmerged index entries and conflict
markers. `git diff --cached --check <incoming>` checks resolution whitespace
without reporting unrelated whitespace already present in the incoming branch.

## Verification shortcuts

Follow [compiling](../compiling/SKILL.md) and
[testing](../testing-and-verification/SKILL.md), with these practical details:

- Build the owning library first, e.g. `make -j"$NPROC" gyrokinetic`, to refresh
  its lower-layer dependencies. Then build individual executable targets
  sequentially. Requesting several executable targets together can launch
  recursive makes that compile the same shared objects concurrently; targeting
  only an executable can also reuse stale lower-layer libraries.
- For loss-cone/wall changes, use `ctest_loss_cone_mask_gyrokinetic`,
  `ctest_gk_species_wall_potential`, and appropriate mirror `-s1` smoke tests.
  The parallel loss-cone test needs **four ranks and `-M -E`**:
  `<matching-mpirun> -np 4 build/gyrokinetic/unit/mctest_loss_cone_mask_gyrokinetic -M -E`.
- A full format failure need not be a merge regression. At `agent_tools` tip
  `85781f6cc`, 52 files already fail `ci/format-all.sh --check`. Compare failing
  files with the actual incoming revision and report inherited failures;
  validate resolved files separately.

In the merge of `8a5d61c5d` with `85781f6cc`, normalization reduced 20 conflicted
files to five requiring review. These counts describe that merge, not future
branches. Complete the requested merge commit after verification.
