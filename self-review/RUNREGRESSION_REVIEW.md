# Review: the `runregression` system

A design review of `gkeyll/lua/Tool/runregression.lua` (2817 lines) and its
companion `gkeyll/lua/Tool/test_costs.lua`, based on a full read of the source,
the committed cost table, and recent git history.

## What's well-designed

The bones are genuinely good. A few choices stand out as thoughtful:

- **Hierarchical layer model** (moments → vlasov → gyrokinetic → pkpm) with
  per-layer SQLite DBs and uniform Lua/C discovery. Clean and scales.
- **On-the-fly C compilation** from the installed `share/Makefile`, with a local
  `rt_arg_parse.h` overriding the installed copy, removes a whole
  `make regression` step.
- **Cost-ordered scheduling** (cheapest-first) to shorten makespan under
  `--jobs N`, with a robust cross-platform timeout layer: GNU `timeout` →
  `gtimeout` → a Perl fork + `setpgrp` + `alarm` fallback that correctly kills
  the entire process group (so a grandchild `gkeyll` launched by `sh` is reached).
- **Combined abs-AND-rel tolerance** in `compareFiles` — the correct way to
  avoid false failures both near zero (tiny abs, huge rel) and at large
  magnitudes (large abs, tiny rel).
- **Walking the *accepted* directory** (not the run dir) in `check_action`, so a
  test that crashes before writing its last frame is flagged as missing rather
  than silently skipped. A subtle, correct choice.

## Issues — most serious first

### 1. The step-cap mechanism has quietly turned the suite into a smoke test

The numbers tell the story. Of 554 tracked tests in `test_costs.lua`:

| num_steps cap | tests |
|--------------:|------:|
| 10            | 445   |
| 3             | 74    |
| 1             | 1     |
| 0             | 33    |

So by default `run check` exercises essentially *no* time evolution. A bug that
manifests at step 500 — instability, slow drift, boundary accumulation,
conservation violation — cannot be caught. The accepted baselines are *also*
created under the same cap (`create` flows through `stepArgFor` too), so a
10-step state is compared against a 10-step baseline.

This is a legitimate "make it run in 5 minutes" tradeoff (recent commits show
that was the explicit goal), but it is currently **invisible**: nothing in the
run summary says "this suite covers a median of 10 steps," and a reader would
reasonably assume `check` tests the physics. The coverage reduction lives in an
auto-generated file, not in the tests.

### 2. `num_steps = 0` converts crashing tests into passing tests

33 tests are pinned to `-s 0`. The code comment is explicit (runregression.lua
lines 469–479): a sim that aborts on its first step records 0, then runs zero
steps "so the test exits cleanly instead of failing." That means a test that
*crashes* is recorded as **green**: it compares initialization-only output
against an initialization-only baseline and reports pass.

Most of these are `rt_gr_*` / `rt_vacuum_einstein_*` / `rt_mhd_ot` — exactly the
family flagged as still broken (`gr-moments-uninit-bugs`). This is a
reward-hacking-shaped hole: the harness reports health it does not have. At
minimum these should carry a distinct status (`SKIP-CRASHES` / xfail) and must
never be counted as `pass`.

### 3. No process exit code reflects test failure

Every `os.exit(1)` in the file is a config/parse error path. `run_action`
returns normally whether 0 or 500 tests fail. There is no
`os.exit(nfail > 0 and 1 or 0)`. For a regression harness this is the single
biggest *automation* gap — CI cannot gate on it without scraping the log or
querying SQLite. One-line fix, outsized value.

### 4. Unit-test runners can false-pass and record nothing

`runLuaUnitTest` / `runCxxUnitTest` only `log` PASS/FAIL — no DB row, no counter,
no exit code. Worse, `runCxxUnitTest` greps stdout for `"FAILED"`: a binary that
**segfaults and prints nothing** is treated as passed. `runLuaUnitTest` keys on
the substring `"PASSED"` with the same fragility. `rununit` therefore cannot
gate anything and can hide crashes.

### 5. Unreadable output silently passes

In `compareFiles`, when `arrayNewFromFile` fails it returns `true` ("skipping …
unsupported file format"). A `.gkyl` file that cannot be read counts as a pass.
A corrupt or truncated output should fail, not be waved through.

### 6. Committed millisecond costs cause churn and aren't portable

`test_costs.lua` conflates two very different things:

- **`num_steps`** — genuine test configuration; belongs in the repo.
- **`cost`** — machine- and run-specific wall-clock floats.

Every `run create` rewrites the cost column, so the file diffs on every run, on
every machine, and the ordering a fast workstation produces won't match a CI
box. The recent history ("Update costs…", "Update test_costs and a t_end") is
this churn. Split them: keep `num_steps` (and any xfail/skip flags) committed as
test config; push the cost cache into the per-layer DB or a gitignored local
file.

### 7. Fixed batches leave cores idle; a worker pool would be strictly better

`executeBatch` launches a window of `jobCount`, then `wait`s for *all* of them
before starting the next window. One slow test in a batch idles every freed core
until the barrier clears. The elaborate cost-sorting is largely a workaround for
this — a dynamic pool (start the next test the moment any slot frees) would beat
fixed batches and make cost-ordering far less load-bearing.

### 8. A real shell-escaping bug in the Perl fallback

`scriptPath:gsub("'", "\\'")` (line 377) is incorrect for POSIX single-quoted
strings — inside `'…'` a backslash does not escape the quote; the correct idiom
is `'\''`. Harmless today because paths are developer-controlled and
quote-free, but it is a latent landmine, and the broader pattern of
`'%s'`-quoting everywhere will break on any path containing a quote.

## Features worth adding

- **Exit code + machine-readable summary** (JSON/TAP). Unlocks CI gating; highest
  leverage.
- **A `diff` subcommand.** The layout deliberately keeps `runs/` beside
  `accepted/` "so the developer can diff," but there is no tooling to show
  *which fields* diverged and by how much. The data is there; the command isn't.
- **In-test annotations via the existing `--!` magic comment.** That hook already
  exists but only does `numProc`-skip. Let tests declare their own `tol`,
  `timeout`, `run_to_completion` (never cap steps), and `xfail` — putting "what
  is tested and why" *in the test* instead of scattered across `test_costs.lua`,
  `ignore_*.lua`, and the 0-step hack (three separate, fragmented skip
  mechanisms today).
- **A small set of "full-physics" tests that run uncapped.** Even 5–10
  long-running tests run to completion would restore confidence that the suite
  catches late-emerging bugs, while keeping the 5-minute smoke pass for the bulk.
- **Baseline provenance.** Accepted `.gkyl` dirs carry no record of the
  commit/build/num_steps that produced them, so you can `check` against stale
  baselines and never know. A baseline manifest would catch this.
- **Flake handling.** History shows flaky tests are dealt with by *deletion*
  ("remove rt_dg_5m_mom_beach_p2, since it sporadically fails"). A
  retry-and-quarantine path would preserve coverage instead of dropping it.
- **MPI coverage.** `mpiExec` is configured but "not yet implemented," and
  MPI-marked tests are skipped outright — domain-decomposition correctness, a
  classic bug source in a physics code, is entirely uncovered.

## Bottom line

The architecture is solid and the engineering details (timeouts, process groups,
abs/rel tolerance, accepted-dir walking) show real care. The concerns cluster
around **integrity**: the step-cap and `num_steps = 0` mechanisms, added under
time pressure to hit a 5-minute target, have shifted the suite from "physics
regression" to "initialization smoke test" — and that shift is invisible in the
output, with 33 crashing tests reporting green.

Prioritize:

1. An honest status for the 0-step tests (never `pass`).
2. A real process exit code on failure.
3. Splitting machine-specific costs out of the committed file.

All three are small changes that restore trust in what a green run actually means.
