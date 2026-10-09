# CI regression tests

[CI overview](README.md)

## Numerical regression differences

Expected numerical changes require a reviewed, new or updated entry in
`expected_regression_diffs.txt`. For a candidate/baseline comparison, CI
honors only lines that are new or changed in the candidate file relative to
the baseline file. Unchanged inherited entries are inert, so a later PR can
safely remove stale lines. Updating an inherited line's reason explicitly
acknowledges a new intentional change for that test.

Entries accept four forms, with an optional `# reason` comment:

```text
test_name                            # entire test, all modes
test_name cpu_parallel               # entire test, CPU parallel only
test_name file_name                  # one file, all modes
test_name cpu_parallel file_name     # one file, CPU parallel only
```

Modes are `cpu_serial`, `cpu_parallel`, `gpu_serial`, and `gpu_parallel`.
Test names may use the full report name, `<layer>/<basename>`, or the bare
basename. Separate fields with spaces or tabs; filenames are matched exactly.
File entries acknowledge any comparison failure for the named file, including
numerical differences and missing or unreadable arrays. Every failing file
must match an entry. Execution failures (crashes, timeouts, compilation
failures, or no output) require an entire-test or entire-mode entry. Tests
still run, and all entries follow the same new-or-changed rule.

A C regression test introduced by the candidate is executed, but is not
numerically compared until it exists in a baseline. CI reports it as
candidate-only. It must still compile, finish without a timeout or crash, and
write output.

## Regression scheduling

Compilation and execution use separate platform settings: `*_BUILD_JOBS` and
`*_REGRESSION_JOBS`, with prefixes `PERSONAL`, `TEAM_WORKSTATION`, `STELLAR_CPU`,
or `PERLMUTTER_GPU`. For example, `PERSONAL_BUILD_JOBS=10` and
`PERSONAL_REGRESSION_JOBS=4` use ten compilation workers and four serial test
workers. Personal/team regression workers default to their build count (three).
Deploy updated trusted Jenkinsfiles for these defaults to take effect.

Pipelines run `run --c-only --jobs N compile`, then
`run --c-only --execute-only --jobs M create` or `check`. Outside CI,
`runregression run --jobs N` uses N workers for each phase; `--jobs 0` detects
CPU count. C and Lua tests share an asynchronous queue: workers take the next
test when free, with no runtime-cost ordering or cost table. Unequal durations
can leave workers idle as the queue drains.

MPI execution treats workers as a rank budget, reserving each test's manifest
rank count until completion. Eight workers allow two four-rank tests at once:

```sh
gkeyll runregression run --c-only --parallel --execute-only --jobs 8 check
```

A test requiring more ranks than the budget runs alone, so `--jobs 1` still
executes MPI tests sequentially. Compilation uses one worker per test.
The MPI launcher controls CPU/GPU placement; worker counts must fit allocated
devices and memory. Stellar allocates at least four MPI tasks, growing to the
regression budget, and uses `srun --exact` to reserve each test's resources.
Perlmutter caps the MPI budget at four allocated GPUs, so its four-rank tests
run one at a time. HPC compilation uses the configured build worker count.

See [scheduler and output tests](tests/README.md#regression-scheduling-and-output)
for verification commands, and [GitHub reporting](README.reporting.md#completed-reports)
for regression timing observations.
