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

The Jenkinsfiles and clients in this directory publish machine-specific GitHub
commit statuses and a per-machine CI report comment (see below). Gkeyll is
public, so Pipeline source and candidate checkouts are anonymous. Each
controller stores its own classic GitHub PAT with only the `public_repo` scope
for authenticated GitHub API requests, status publication, and report
comments; an outside collaborator with push access can use this model without
Gkeyll organization membership. See the platform guide for the
credential owner. Do not place credentials in the repository or in candidate
branches.

## Results in GitHub

Every Jenkins controller is private, so the commit status alone would only show
red or green. Each run therefore also posts a Markdown report to GitHub, built
by `github_report.py` from the summaries the Pipeline archives: candidate and
baseline commits, the failing stage with the extracted compiler or test error,
the unacknowledged and acknowledged regression tests, and timings. For a PR
run the report is a pull-request comment; for a candidate/baseline run it is a
comment on the candidate commit. The comment carries a hidden marker keyed on
the machine's status context, so each machine owns one comment per PR and
later runs update it in place instead of adding new ones.

Publishing is best effort: if GitHub is unreachable the build result is
unchanged and the status description ends with "(report not posted)". The
Pipeline runs `github_report.py` only from reviewed code (the trusted
team-workstation checkout, otherwise `main`, or `GKEYLL_CI_TRUSTED_REF` when
staging a CI change), never from the candidate checkout, because the GitHub
token is bound while it runs.
The regression-result checker also comes from the reviewed trusted CI commit,
recorded in `ci-trusted-checker-commit.txt`. CI runs it with the baseline
executable from the baseline source directory, using `-S` because this Lua
check does not require MPI.

## Numerical regression differences

Expected numerical changes require a reviewed, new or updated entry in
`expected_regression_diffs.txt`. For a candidate/baseline comparison, CI
honors only lines that are new or changed in the candidate file relative to
the baseline file. Unchanged inherited entries are inert, so a later PR can
safely remove stale lines. Updating an inherited line's reason explicitly
acknowledges a new intentional change for that test.

A C regression test introduced by the candidate is executed, but is not
numerically compared until it exists in a baseline. CI reports it as
candidate-only. It must still compile, finish without a timeout or crash, and
write output.

## Unified local command

### Regression worker budget

`gkeyll runregression run --jobs N` uses a worker pool for C compilation and
test execution. Each finished job releases its slots immediately so the next
test can start without waiting for a whole batch. Compilation finishes before
execution starts, with each compiler invocation limited to `make -j1`.

Serial tests use one slot; `--parallel` tests reserve one slot per MPI rank.
For example, `--parallel --jobs 12` can launch a four-rank test and an eight-rank
test concurrently, or three four-rank tests.
A test requiring more ranks than the budget runs alone. The default remains
`--jobs 1`; `--jobs 0` detects the available CPU count.

Worker slots only account for concurrency within one runregression invocation.
They do not assign CPU IDs or reserve physical cores, and separate invocations
have independent budgets. The MPI launcher controls CPU and GPU placement;
GPU jobs share the available devices. Set the budget to fit your allocation,
and configure the launcher for concurrent tests as described below.

### Concurrent MPI tests and placement

Configure launcher options with `runregression configure --mpi-args`; these
options apply to each MPI launch. Keep your existing source, prefix, results
subdirectory, and other configuration options when configuring your suite.

For Open MPI, independent `mpirun` launches can bind ranks to the same cores.
For concurrent tests on a shared CPU allocation, disabling binding lets the
operating system schedule ranks across the available CPUs:

```sh
gkeyll runregression configure --source-dir /path/to/gkeyll \
  --mpiexec /path/to/openmpi/bin/mpirun --mpi-args="--bind-to none"
gkeyll runregression run --parallel --jobs 12 check
```

`--bind-to none` does not give tests exclusive cores. Guaranteed disjoint core
placement requires a resource manager or launcher wrapper that assigns separate
CPU sets to concurrent tests; runregression does not allocate those sets.
See the [Open MPI launcher documentation](https://docs.open-mpi.org/en/v5.0.10/man-openmpi/man1/mpirun.1.html).

For Slurm, use `--exact` when launching concurrent steps inside an existing job
allocation. By default, a step can consume all non-GRES resources on its nodes,
causing another step to wait even when its rank count fits the worker budget.
`--exact` limits each step to its requested resources:

```sh
gkeyll runregression configure --source-dir /path/to/gkeyll \
  --mpiexec /path/to/srun --mpi-args="--exact"
# Run inside a Slurm allocation with sufficient resources for these ranks.
gkeyll runregression run --parallel --jobs 12 check
```

Slurm placement still depends on the allocation, site configuration, and other
resource requests, including CPUs per task, node count, and GPUs. `--exact`
does not guarantee that every step can start immediately. See the
[Slurm srun documentation](https://slurm.schedmd.com/srun.html).

`--timeout` measures the whole launcher invocation, including startup and time
spent waiting for a Slurm step allocation. Waiting for runregression's own worker
slots occurs before this timer starts.

Use `--parallel --jobs 1` to serialize MPI tests when concurrent placement is
unsuitable. Each test still launches its full manifest-defined MPI rank count;
`--jobs 1` does not reduce it to one rank. This prevents overlap between tests
within that invocation only.

### CI client

The script `gkeyll-ci.sh` provides a CLI to run, query and terminate CI. See

```sh
./ci/jenkins/gkeyll-ci.sh -h
./ci/jenkins/gkeyll-ci.sh --help
```

## Candidate freshness

Before building, CI resolves both references to commit SHAs and requires the
candidate to contain the baseline commit. A behind candidate fails immediately,
before dependency builds or Slurm submission. Update the candidate with its
baseline before rerunning. For an intentional historical comparison, pass
`--allow-behind-candidate` (or enable `ALLOW_BEHIND_CANDIDATE` in Jenkins); the
CI artifact and GitHub report record that override.

## Persistent regression data

Every Jenkins controller requires a writable, agent-visible `GKEYLL_CI_ROOT`.
A complete baseline lives at
`$GKEYLL_CI_ROOT/baseline-cache/<platform>/<baseline-sha>/` and is reused
until the resolved baseline SHA changes. A failed or incomplete baseline build
is removed after diagnostics are collected. A later candidate failure leaves
a valid baseline intact. The latest candidate attempt, including its source,
build, and results, stays at
`$GKEYLL_CI_ROOT/candidate-cache/<platform>/<build-tag>/<candidate-sha>/`
until the next run starts. Every run rebuilds the candidate, including when
its SHA is unchanged. Within each build tree, `gkeyll/` is the source checkout
and `gkylsoft/` contains the install and regression results.

At startup, `cache_manager.py` removes inactive prior candidates and, on
Stellar and Perlmutter, inactive build-tag workspaces. It preserves trees used
by a running Jenkins build or surviving Slurm job. Concurrent runs can
temporarily leave multiple candidate and baseline trees. When their cleanup
stages finish, the candidate from the most recently started run remains,
regardless of which run finished last. A newly started run removes that
completed candidate; at most one completed candidate and one valid baseline
remain after all runs finish. The cleanup uses a
filesystem lock and queries Jenkins and Slurm before deleting an old tree.
If either required query fails, CI stops cleanup and preserves the files.
`cache-runs/` stores small ownership records and `cache-control/` contains
reviewed cache scripts for each run. Jenkins still archives summaries and
small diagnostics and retains its configured build history separately.

The Jenkins agent must be able to read the Jenkins build-state API using
`CI_JENKINS_READ_AUTH_FILE`, a mode-600 `user:api-token` file, or anonymous
read access when that variable is empty. Its default is the platform CLI
credential file. Set it explicitly when the agent runs under a different Unix
account. On cluster systems, keep the root on shared storage that supports
cross-process `flock`, and provide enough Jenkins executors for the desired
parallel builds. Drain old CI runs before first using this layout so the
legacy SHA-only candidate cache can be removed safely.
The `<platform>` directory is `personal`, `team-workstation`, `stellar-cpu`,
or `perlmutter-gpu`, depending on the Jenkins job.

It can be used with any of the machines listed above, for example:

```sh
./ci/jenkins/gkeyll-ci.sh personal run --pr 1128 --follow
./ci/jenkins/gkeyll-ci.sh stellar_cpu recent
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu status --queue 42
./ci/jenkins/gkeyll-ci.sh stellar_cpu info --build 6
./ci/jenkins/gkeyll-ci.sh stellar_cpu artifact --build 6 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh team scan
```

The wrapper runs the chosen platform client locally and passes through all
arguments and environment variables. Use `./ci/jenkins/gkeyll-ci.sh --help`
to list platforms, or `./ci/jenkins/gkeyll-ci.sh PLATFORM --help` for the
selected client's commands and setup requirements.

### Build details and artifacts

`info --build NUMBER` shows retained build metadata, its archived failure
summary when present, queryable regression failures, and every retained
artifact. `artifact --build NUMBER` lists those artifacts without opening a
browser. Add `--fetch` to download all artifacts into a new
`gkeyll-ci-build-NUMBER` directory, or use `--only PATH[,PATH...]` and
`--output-dir DIR` to select artifacts and a new destination. Downloaded
artifacts retain their Jenkins-relative directory hierarchy.
