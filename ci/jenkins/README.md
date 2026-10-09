# Gkeyll Jenkins CI

Gkeyll has independent Jenkins setups for machines that provide persistent or
specialized CI capacity. Each setup owns its Jenkins job, node configuration,
credentials, and command-line client.

- [Personal computer CI](README.personal.md): manually run a selected GitHub
  pull request or candidate/baseline branch or commit on a local computer.
- [Stellar CPU CI](README.stellar_cpu.md): manual CPU-only Slurm CI at
  Princeton Stellar.
- [Perlmutter GPU CI](README.perlmutter_gpu.md): manual GPU Slurm CI at NERSC.
- [Team workstation CI](README.team_workstation.md): periodically discover
  pull requests into `main` from every author, with optional selected runs.

Shared guides cover [installed dependencies](README.dependencies.md),
[GitHub reporting and controller setup](README.reporting.md),
[regression scheduling and expected differences](README.regressions.md),
[Valgrind memory checking](README.valgrind.md),
[storage and caching](README.storage.md), and [CI development tests](tests/README.md).

## Unified local command

Use `gkeyll-ci.sh` to run, query, and terminate CI on any platform:

```sh
./ci/jenkins/gkeyll-ci.sh --help
./ci/jenkins/gkeyll-ci.sh personal run --pr 1128 --follow
./ci/jenkins/gkeyll-ci.sh stellar_cpu recent
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu status --queue 42
./ci/jenkins/gkeyll-ci.sh stellar_cpu info --build 6
./ci/jenkins/gkeyll-ci.sh stellar_cpu artifact --build 6 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh team scan
```

The wrapper runs the selected client locally, passing through arguments and
environment variables. Use `-h` for short help or `PLATFORM --help` for that
client's commands and setup requirements.

`info --build NUMBER` shows build metadata, failure summaries, regression
failures, and retained artifacts. `artifact --build NUMBER` lists artifacts;
add `--fetch` to download them into a new `gkeyll-ci-build-NUMBER` directory.
Use `--only PATH[,PATH...]` and `--output-dir DIR` to select files and a new
destination. Downloads preserve the Jenkins-relative directory hierarchy.

## Results in GitHub

Each machine publishes a commit status while queued, running, and completed.
The controller listener reports queue position and running time/ETA, pins the
accepted candidate SHA, and cancels superseded PR runs before agent allocation.
Manual platforms need a new submission to test the replacement commit.

Completed runs post a report on the **exact tested commit**, including PR runs.
The status's **Details** link opens it. Reruns update one report per commit and
machine context. Reports include unit and regression results, compiler
diagnostics, warnings relative to the baseline, stage durations, the 20 slowest
regressions, and timing changes. Reporting failures do not change test results.
Optional Valgrind checks include complete memory-error logs in the report.

Source discovery and checkouts are anonymous. Each controller uses its existing
GitHub username/token credential for statuses and comments; see
[reporting credentials](README.reporting.md#reporting-credentials).
Do not put credentials in the repository or candidate branches. Reporting and
regression comparison use reviewed, trusted CI code, independently of the
candidate. See the [reporting guide](README.reporting.md) for installation,
trust configuration, and failure diagnostics.

## Candidate freshness

Before building, CI resolves both references to commit SHAs and requires the
candidate to contain the baseline commit. A behind candidate fails before
dependency builds or Slurm submission. Update the candidate with its baseline
before rerunning. For an intentional historical comparison, pass
`--allow-behind-candidate` (or enable `ALLOW_BEHIND_CANDIDATE` in Jenkins);
the artifact and report record the override.

## Numerical regression differences

Expected changes require a reviewed, new or updated entry in
`expected_regression_diffs.txt`. Only lines added or changed relative to the
baseline apply; unchanged inherited entries are inert. Entries can acknowledge
a whole test, one CPU/GPU serial/parallel mode, or individual output files.
See [entry syntax and scheduling](README.regressions.md).

New C regression tests run without numerical comparison until they exist in a
baseline. They must still compile, finish without a timeout or crash, and write
output; reports label them candidate-only.

## Dependencies and retained data

Set `GKEYLL_CI_PREBUILT_CONFIG` to reuse an installed dependency configuration,
or leave it unset to build dependencies with the platform's machine scripts.
The [dependency guide](README.dependencies.md) covers configuration, MPI,
and ADAS data.

Every controller needs a writable, agent-visible `GKEYLL_CI_ROOT`. Each build
retains its candidate source, executable, results, and baseline snapshot under
`runs/<platform>/<BUILD_TAG>/`; `ci-run-path.txt` records the path. A compatible
baseline is reused from `baseline-cache/<platform>/<baseline-sha>/`.

Jenkins archives full regression results, then removes its temporary workspace.
Runs and Jenkins artifacts have no automatic retention limit by default.
See [storage and caching](README.storage.md) for layout, cache invalidation,
concurrency, and cleanup.
