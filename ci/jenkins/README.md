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

The Jenkinsfiles publish machine-specific GitHub commit statuses for queue,
stage, progress, and final results. Completed reports are readable PR or commit
comments, updated in place, and linked from the status's **Details** button.
Reporting uses the existing username/token credential. Public source discovery
and checkouts remain anonymous. See [Reporting credentials](#reporting-credentials)
and the platform guide for setup. Do not place credentials in the repository
or in candidate branches.

## Reusing installed dependencies

All four workflows accept the optional Jenkins environment variable
`GKEYLL_CI_PREBUILT_CONFIG`. Set it in **Manage Jenkins → System → Global
properties → Environment variables** to the absolute path of an existing
Gkeyll `config.mak` on the build agent, for example:

```text
GKEYLL_CI_PREBUILT_CONFIG=/opt/gkylsoft/gkeyll/share/config.mak
```

You can use either the `config.mak` from a configured source checkout or the
copy saved by `make install` in `<prefix>/gkeyll/share/config.mak`. Keep the
config and dependencies outside Jenkins workspaces, which are deleted during
checkout. The Jenkins agent must be able to read them; on Slurm machines the
libraries must also be visible at the same paths on compute nodes.

Specify each dependency's include and library paths in this config. A generated
config already contains these settings; for dependencies installed elsewhere,
make a dedicated CI copy and edit the appropriate entries:

| Dependency | Path variables in `config.mak` |
| --- | --- |
| BLAS/LAPACK | `LAPACK_INC_DIR`, `LAPACK_LIB_DIR` (and `LAPACK_LIB_NAME`) |
| SuperLU | `SUPERLU_INC_DIR`, `SUPERLU_LIB_DIR` (and `SUPERLU_LIB_NAME`) |
| MPI | `CONF_MPI_INC_DIR`, `CONF_MPI_LIB_DIR` |
| LuaJIT | `CONF_LUA_INC_DIR`, `CONF_LUA_LIB_DIR` (and `CONF_LUA_LIB`) |
| NCCL | `CONF_NCCL_INC_DIR`, `CONF_NCCL_LIB_DIR` |
| cuDSS | `CONF_CUDSS_INC_DIR`, `CONF_CUDSS_LIB_DIR` |
| CUDA math libraries | `CUDAMATH_LIB_DIR` |

Use absolute paths for enabled dependencies, or expressions relative to the
config's original `PREFIX`, such as `$(PREFIX)/superlu/lib`. Supply a complete,
self-contained generated config with `PREFIX` defined; relative includes and
paths relative to the original source checkout are not supported. Its compiler,
architecture, solver and feature settings are reused too. Choose settings
compatible with the machine's existing CI lanes: Lua and MPI are needed by the
regression runner, and Perlmutter requires the CUDA/NCCL configuration. Existing
HPC module-loading steps still run; dependencies must match those modules.

For personal CI also specify `PERSONAL_MPIEXEC=/opt/mpi/bin/mpiexec` or
`PERSONAL_MPI_HOME=/opt/mpi`, matching the MPI in the config. For team-workstation
CI use `TEAM_WORKSTATION_MPIEXEC` or `WORKSTATION_MPI_HOME`. HPC workflows retain
their existing Slurm launchers. `PERSONAL_MKDEPS_SCRIPT` and
`PERSONAL_CONFIGURE_SCRIPT` (or their `TEAM_WORKSTATION_*` equivalents) are not
required in prebuilt mode.

With this variable set, candidate and baseline both skip the machine dependency
and configure scripts. Each gets its own copy of the config with Gkeyll's
`PREFIX` and `INSTALL_PREFIX` redirected to its workspace `gkylsoft` directory.
Dependency references to the old prefix are preserved. The supplied config and
installed dependencies are not modified. Gkeyll itself is still rebuilt and
all existing test lanes run. An invalid config path fails the build rather than
silently rebuilding dependencies. Leave the variable unset or empty to retain
the existing dependency-build behavior. Deploy the updated trusted Jenkinsfiles
to enable this option; candidate and baseline refs can predate it.

The offline config test uses Groovy 2.4 and GNU Make:

```sh
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/test_prebuilt_config.groovy
```

## Results in GitHub

### Reporting credentials

Use the platform's existing `*_GITHUB_CREDENTIAL_ID` Jenkins **Username with
password** credential, with the GitHub username and token as its password.
An existing classic PAT with `public_repo` scope covers status and comment
publication on this public repository. A fine-grained token needs **Commit
statuses: write**, **Pull requests: write** (PR comments), and **Contents: write**
(commit comments) on `gkeyllorg/gkeyll`; its owner needs push access for statuses.
There is no additional credential or reporting plugin to configure. See GitHub's
[commit status API](https://docs.github.com/en/rest/commits/statuses),
[PR comment API](https://docs.github.com/en/rest/issues/comments), and
[commit comment API](https://docs.github.com/en/rest/commits/comments).

Each controller must use its own stable status context. Keep existing contexts
so branch protection continues to use the same names. Retries update the same
report comment for that context. Delivery diagnostics are archived in
`ci-report-delivery.json`; reporting outages do not change test results.

### Queued builds and superseded PR commits

Install [queue-status.groovy](queue-status.groovy) on **each Jenkins controller**
to publish a yellow `pending` commit status as soon as a job is accepted into
its queue. This covers CLI and browser submissions on personal, Stellar CPU, and
Perlmutter GPU controllers, plus automatic and manual team-workstation builds.
The listener runs without a build executor, including when
`disableConcurrentBuilds()` blocks another run of the same job. It also covers
Pipelines waiting for their first `node()` allocation.

The pending description includes the position, for example
`Gkeyll CI queued: position 3 of 10 (Jenkins queue #42).` Positions count waiting
Gkeyll builds for the same platform on that controller in original submission
order, including the first agent wait, and exclude builds already allocated an
agent. They refresh after queue changes and every minute; unchanged positions
are not reposted. This is a submission-order position, not a guaranteed execution
order: Jenkins can skip blocked jobs or allocate multiple available executors.

Once an agent is allocated, the same status refreshes every minute with elapsed
build time, an approximate percentage, and estimated remaining time, for example
`Gkeyll CI running: started 29 m ago; ~20%; est. remaining 1 hr 50 m (Jenkins queue #42).`
The estimate comes from Jenkins' historical build duration, not completed test
counts; elapsed time uses Jenkins' build start time (including the Pipeline's
initial agent wait). Runs with no estimate show `ETA unavailable`; runs exceeding
the estimate say so instead of claiming completion. Unchanged descriptions are
not reposted, and final results replace progress. Stage names are shortened
before timing information so the ETA remains readable. Queue position does not
predict a start time: blocked agents and Slurm scheduling make that uncertain.
Install the updated controller hook and trusted Jenkinsfiles together to enable
progress for all four platforms.

The listener resolves and records the candidate SHA before releasing the job
to run. The Pipeline checks out that exact SHA and updates the same status to
running and then its final result. A queued build therefore cannot silently
switch to a different commit when a branch moves.

For PR selections, the controller checks for newer head commits after
submissions and every minute while work is queued. It cancels superseded runs
before they acquire their first agent, and updates the old commit's status to
`error` with `Gkeyll CI cancelled while queued; superseded by <commit>.` It does
not abort work that already acquired an agent, cancel a different PR, or cancel
a second request for the same commit. Branch/SHA selections remain pinned and
are not automatically cancelled. Cancellation does not submit a replacement:
team-workstation discovery supplies new runs; manual platforms need a new
submission. GitHub commit statuses have no separate `cancelled` state, so
cancellation clears yellow using `error` and an explanatory description.

Manual queue cancellation and failures before the Pipeline can report also
close the pending status. Cleanup checks the current GitHub status before
writing, to preserve a newer run's status and the Pipeline's detailed results.
GitHub outages are logged and do not prevent builds from starting; notification
and supersession checks require a working GitHub API and valid credentials.

#### Controller installation

First deploy the updated trusted Jenkinsfiles along with the listener. It uses
the existing **Credentials** plugin and the platform's existing username/PAT
credential. Set the credential ID in Jenkins **global environment variables**
as described in the platform guide; node-only environment variables are not
available while waiting for an executor. Keep the team GitHub source anonymous:
this listener owns queue notifications, independently of GitHub Branch Source.

On personal controllers, explicitly set `PERSONAL_STATUS_CONTEXT` to the
existing context, e.g. `continuous-integration/jenkins/personal-mycomputer`.
The controller cannot infer the build agent's hostname before allocation.

From a **reviewed, trusted checkout**, as the Jenkins service owner, install:

```sh
install -d -m 700 "$JENKINS_HOME/init.groovy.d"
install -m 600 ci/jenkins/queue-status.groovy \
  "$JENKINS_HOME/init.groovy.d/gkeyll-queue-status.groovy"
```

`JENKINS_HOME` is the actual controller home, not a workspace: typically
`/var/lib/jenkins` for the workstation service, or
`$GKEYLL_CI_ROOT/jenkins_home` for Stellar/Perlmutter. For system services, use
an administrator to install the file with the Jenkins service user's ownership.
Restart the controller during an idle maintenance window. Check its log for
`Gkeyll GitHub queue status listener installed`. Install the hook once; do not
also execute it in the Script Console. Updating the repository alone does not
update the installed controller hook. To remove it, delete the installed hook
and restart; the Pipelines retain their existing running/final reporting.

The default job names come from the platform guides. For renamed or staging
jobs, override the appropriate global variable with the **full Jenkins job
name**, including folders:

| Platform | Optional job-name override |
| --- | --- |
| Personal | `GKEYLL_CI_QUEUE_PERSONAL_JOB` |
| Stellar CPU | `GKEYLL_CI_QUEUE_STELLAR_CPU_JOB` |
| Perlmutter GPU | `GKEYLL_CI_QUEUE_PERLMUTTER_GPU_JOB` |
| Team workstation | `GKEYLL_CI_QUEUE_TEAM_WORKSTATION_JOB` (multibranch parent) |

The hook only handles these configured jobs (direct children for the team
multibranch job). `PERSONAL_STATUS_CONTEXT` and
`TEAM_WORKSTATION_STATUS_CONTEXT` must match the Pipeline configuration.
Queue metadata is saved under `$JENKINS_HOME/gkeyll-queue-status/` and transferred
to internal `CI_QUEUE_*` build parameters when the Pipeline starts. This keeps
the accepted SHA across restarts without disrupting Jenkins' queue deduplication.
Do not add these internal parameters to the browser's input form.

Validate installation by keeping the selected agent unavailable, submitting a
PR, and confirming yellow appears before any candidate checkout. Push another
commit to that PR and confirm the old queue entry disappears within the next
successful sweep and its GitHub description says it was superseded. Then test
manual cancellation and a normal running/completed build.

### Completed build reports

Every completed or failed run attempts to publish a Markdown comment on the PR,
or on the exact candidate commit for branch/SHA selections. Stage and command
transitions update the commit status immediately; the controller refreshes
queue positions and elapsed-time/ETA indicators every minute. Progress does not
create PR comment notifications. Cancellation and timeout descriptions are
explicit; both use GitHub's `error` status. Reports identify the candidate,
baseline, machine, build and queue IDs, and status context, with unit/regression
results, timings, and collapsible sections for:

- **New warnings vs main** (or the explicitly selected baseline).
- **All warnings** and **all errors**, with source log names, line numbers, and
  nearby diagnostic context.
- **Failed build log — last 100 lines**, immediately below the failed stage,
  in a collapsed dropdown preserving compiler commands, source lines, and carets.
  This excerpt preserves the original output order (including any warnings);
  extracted warnings and errors remain in their own separate dropdowns.
- An inventory of captured logs.

Reports also show the full Jenkinsfile, reporting-tool, and regression-checker
commits when recorded. Personal CI's `CI_REF` selects an implementation
independently of candidate and baseline; see the
[personal guide](README.personal.md#cli-launch). The loader logs its resolved
commit before running that implementation and archives it on completion. The controller queue
listener is installed separately and is not replaced by a per-run selection.

The reporter reads build/configuration logs, Slurm output, and the individual
compiler/runtime logs stored in regression databases. Warning comparison uses
matching completed baseline steps from the same run. It ignores checkout paths,
ANSI colors, and source line/column changes. Missing or failed baseline steps are
reported as unavailable, never treated as a clean baseline. Warnings from shared
or unmatched logs remain visible in **All warnings**. A selected baseline other
than `main` is labeled explicitly. HPC runs also compile baseline unit tests so
candidate unit-build warnings have a corresponding baseline.

Large diagnostic sections continue in linked report comments, published before
the main report is updated. Pages leave room within GitHub's size limit and keep
Markdown fences and collapsed sections balanced. The reporter only edits
comments written by the authenticated account with the exact context marker;
it clears obsolete continuation pages when a later report is shorter. Accepted
queue IDs (or build numbers within a job) prevent an older run from replacing
a newer report already on the PR. The archived `ci-report.md` contains the
entire report; `ci-report.md.json` contains the report pages. Raw logs remain
available through the artifact command below. `ci-stage-history.txt` records
stage changes in UTC, including stages before candidate checkout.

Reporting runs during failure cleanup. If no exact candidate SHA was resolved,
the generated report stays in Jenkins artifacts. Timing-summary errors do not
skip publication. Network/timeouts and transient server errors receive up to
three attempts; each comment retry discovers existing pages first to recover
from a lost POST response. Status publication is attempted even if comments
fail, and vice versa. The delivery artifact records which operation failed.
The controller closes pending statuses for early failures and preserves detailed
terminal results against late progress updates. A persistent GitHub outage or
an unavailable agent can still prevent delivery; consult the build log and
archived report in that case.

Report generation and GitHub publication live in `github_report.py`. Shared
Pipeline stage/command reporting lives in `jenkins_reporting.groovy`; the
controller only needs `queue-status.groovy`. The older standalone `build` and
`publish` commands remain supported. Pipelines load both shared report files
before candidate checkout from reviewed code (`main`, the pinned trusted team
checkout, or `GKEYLL_CI_TRUSTED_REF` while staging changes). Deploy those files
with the updated Jenkinsfiles. Diagnostic collection runs before binding the
GitHub token; credential-bearing shell steps are not captured in published logs.

The regression-result checker also comes from the reviewed trusted CI commit,
recorded in `ci-trusted-checker-commit.txt`. CI runs it with the baseline
executable from the baseline source directory, using `-S` because this Lua
check does not require MPI.

Run offline reporter tests with:

```sh
python3 -m unittest discover -s ci/jenkins -p 'test_*.py'
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/test_jenkins_reporting.groovy
```

The queue listener's integration test runs on a disposable Jenkins controller
with mocked GitHub responses. It checks all four configurations, supersession
before Pipeline start and during agent wait, manual cancellation, active-build
preservation, queue persistence data, and reporting failures. Install Credentials,
Folders, Pipeline: Job, Pipeline: Groovy, Pipeline: Basic Steps, and Pipeline:
Nodes and Processes, and Credentials Binding (including dependencies) in that
test controller. With a
Jenkins WAR and those plugin archives available locally:

```sh
queue_test_home=$(mktemp -d)
mkdir -p "$queue_test_home/plugins" "$queue_test_home/init.groovy.d"
cp /path/to/test-plugin-archives/*.jpi "$queue_test_home/plugins/"
cp ci/jenkins/test_queue_status.groovy "$queue_test_home/init.groovy.d/90-queue-test.groovy"
JENKINS_HOME="$queue_test_home" java \
  -Djenkins.install.runSetupWizard=false -Dgkeyll.queue.test=true \
  -Dgkeyll.queue.source="$PWD/ci/jenkins" \
  -jar /path/to/jenkins.war --httpListenAddress=127.0.0.1 --httpPort=18089
cat "$queue_test_home/queue-test-result.txt"
```

This test requires an empty job directory and shuts down its test JVM on
completion. Never install `test_queue_status.groovy` on a production controller.

## Regression scheduling

Jenkins configures compilation and execution separately using each platform's
`*_BUILD_JOBS` and `*_REGRESSION_JOBS` environment variables. For example, set
`PERSONAL_BUILD_JOBS=10` and `PERSONAL_REGRESSION_JOBS=4` for 10 compilation
workers and 4 concurrent serial regression runs. Use the `TEAM_WORKSTATION`,
`STELLAR_CPU`, or `PERLMUTTER_GPU` prefix for the other workflows. The pipelines
invoke `run --c-only --jobs N compile` first, then `run --c-only --execute-only
--jobs M create` or `check`. MPI regressions continue to run one test at a time
with four MPI ranks per test.

`runregression run --jobs N` uses up to N compilation workers, then up to N
execution workers. `--jobs 0` detects the available CPU count. C and Lua tests
share an execution queue without runtime-cost ordering. Workers take the next
test as soon as a slot is free, without waiting for a batch to finish. No cost
table or scheduling artifacts need maintaining. Unequal test durations can
leave workers idle as the queue drains at the end.

MPI regression collectives still execute one test at a time; `--jobs`
parallelizes their compilation. GPU worker counts must fit the devices and
memory allocated to the job.

HPC compile stages pass their configured build worker count to runregression;
execution uses the existing regression worker allocation. Personal and team
workstation regression workers default to their build worker count (three),
with the existing regression-jobs setting available as an override. Updated
trusted Jenkinsfiles must be deployed for these CI defaults to take effect.

The scheduler fixtures use tiny shell/make jobs, not plasma simulations:

```sh
LUAJIT=/path/to/luajit python3 -m unittest discover -s ci/jenkins -p 'test_*.py'
```

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
CI retains one complete baseline per platform in
`$GKEYLL_CI_ROOT/baseline-cache/<platform>/<baseline-sha>/`. A repeated run
against the same resolved baseline SHA reuses that tree; a changed SHA rebuilds
and replaces it. The candidate source, build, installed executable, and
regression results are kept together in
`$GKEYLL_CI_ROOT/candidate-cache/<platform>/<candidate-sha>/`. CI rebuilds the
candidate from scratch on every run, even when its SHA is unchanged, and keeps
only the latest candidate tree per platform.
Within either SHA directory, `gkeyll/` is the source checkout and `gkylsoft/`
is its sibling install and results directory. For example, the executable is
`<sha>/gkylsoft/gkeyll/bin/gkeyll` and regression output is under
`<sha>/gkylsoft/gkeyll-results/`.

Both trees are built at their final paths so installed libraries retain valid
absolute paths. An incomplete baseline build has no valid cache manifest and
is rebuilt on the next run. Jenkins keeps summaries and small diagnostic
artifacts in its workspace for reporting, then removes that workspace; the
complete candidate build and results remain under `candidate-cache`.
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
