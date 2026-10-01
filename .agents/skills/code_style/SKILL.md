---
name: code_style
description: Gkeyll's C/C++/CUDA code style is enforced by clang-format. Use when writing or editing any .c/.h/.cpp/.hpp/.cu/.cuh file, to know how to format it and which files must never be reformatted.
user-invocable: true
---

# Instructions

* Operate relative to the repo root (detect via `git rev-parse --show-toplevel`).

# Code style

The style itself (indentation, brace placement, spacing, etc.) is defined by the
root `.clang-format` file plus one rule in `ci/format.py` (designated-initializer
layout, see below). Don't try to memorize or restate those rules -- just run the
formatter and let it decide:

```
python3 ci/format.py <file>...
```

Do not run bare `clang-format -i`: it packs designated initializers back onto one
line. `ci/format.py` needs clang-format 18.1.8 on PATH (`pip install clang-format==18.1.8`).

A `pre-commit` hook runs this automatically on commit, and CI (`.github/workflows/format-check.yml`)
double-checks on push/PR, so a file that hasn't been run through clang-format will
get reformatted or flagged regardless.

## The rules the config file can't express: never format `ker/` or `core/minus/`

Two kinds of files must never be run through clang-format, and never hand-formatted
to "match" the style either -- leave them exactly as they are:

- Any path matching `*/ker/*` (e.g. `core/ker/`, `gyrokinetic/ker/`) contains
  auto-generated DG kernel code (see `gkeyll_guide`'s note on `ker/`: generated with
  Maxima, never hand-edited).
- Anything under `core/minus/` is a vendored third-party library (e.g. sqlite3,
  kann, pcg_basic, STC), not Gkeyll's own code.

The pre-commit hook and CI already exclude both for this reason.

## Designated initializers: the script decides, not the trailing comma

A designated initializer that fits on one line stays on one line. One that does
not gets one member per line, first member below the `{`, closing `}` on its own
line. `ci/format.py` enforces this by stripping every trailing comma, running
clang-format, re-adding a trailing comma to each list that did not fit, and running
clang-format again. So don't add or remove trailing commas to steer the layout;
write the list any way you like and run the script.

## Two conventions the formatter enforces that you should know when reading code

- Every function *definition* has its return type on its own line, so
  `grep '^gkyl_foo('` finds the definition. Declarations keep the type on the same line.
- `case` labels are indented one level inside `switch`.
