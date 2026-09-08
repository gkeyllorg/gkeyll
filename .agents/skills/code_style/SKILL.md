---
name: code_style
description: Gkeyll's C/C++/CUDA code style is enforced by clang-format. Use when writing or editing any .c/.h/.cpp/.hpp/.cu/.cuh file, to know how to format it and which files must never be reformatted.
user-invocable: true
---

# Instructions

* Operate relative to the repo root (detect via `git rev-parse --show-toplevel`).

# Code style

The style itself (indentation, brace placement, spacing, etc.) is defined entirely
by the root `.clang-format` file. Don't try to memorize or restate those rules --
just run the formatter and let it decide:

```
clang-format -i <file>
```

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

## Trailing commas in initializer lists are stripped automatically

Don't add a trailing comma after the last element of a struct/array
initializer expecting it to force (or preserve) a multi-line layout --
clang-format expands every element onto its own line whenever one is
present, which is rarely what's wanted. `ci/strip-trailing-commas.py` removes
any comma directly before a closing `}` before clang-format runs, via the
same pre-commit hook and CI check as clang-format itself, so there's also no
need to manually remove existing ones -- the tooling already does it.
