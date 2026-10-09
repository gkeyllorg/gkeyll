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
stage, progress, and final results. Completed reports are comments on the exact
tested commit, including PR runs, updated in place per commit and status context,
and linked from the status's **Details** button.
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
checkout. The Jenkins agent must be able to read them. On Slurm machines the
CI run directory must be visible at the same path on compute nodes.

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
and configure scripts. CI copies the configured dependency installations into
`$GKEYLL_CI_ROOT/runs/<platform>/<BUILD_TAG>/dependencies/`, following symlinks
to keep those copies independent of the originals. The candidate gets a new
copy each run. The baseline owns a separate persistent copy under
`baseline-cache/<platform>/<baseline-sha>/dependencies/`, reused with its
executable and results. The supplied config and original dependencies remain untouched. OS libraries under
`/usr`, `/lib`, `/lib64`, and `/System`, compilers, and platform SDKs remain
machine-provided; they are not copied into the run.

ADAS runtime data is copied from the original `PREFIX/gkeyll/share/adas` into
`dependencies/data/adas`, then copied into each build's own
`gkylsoft/gkeyll/share/adas` before compilation and tests. These are independent
directories so installing a revision's radiation fits cannot modify the other
build, the dependency snapshot, or the original installation. Prebuilt
gyrokinetic and PKPM configs require ADAS `.npy` tables in the original
installation; setup fails early if none are present. Prepare that installation
with `--build-adas=yes` before using it for CI. Run `mkdeps.sh` from
`install-deps/` and pass the same absolute prefix used by the supplied config;
otherwise it defaults to the invoking user's `$HOME/gkylsoft`. For an existing
dependency installation, only the data needs preparing:

```sh
# Run from the repository root, with Python and NumPy available.
set -e
prebuilt_prefix=/opt/gkylsoft  # Replace with PREFIX from the prebuilt config.
(
    cd install-deps
    ./mkdeps.sh --prefix="$prebuilt_prefix" --build-adas=yes
)
# Check every generated table; a lone .npy file is not a complete installation.
for element in h he li be b c n o ar; do
    for table in ioniz recomb logT logN; do
        test -s "$prebuilt_prefix/gkeyll/share/adas/${table}_${element}.npy"
    done
done
```

When starting from scratch, build the libraries with the appropriate machine
dependency script (including ADAS), then configure Gkeyll with that same prefix
and the intended compiler, MPI, Lua, and CPU/GPU options. Retain the generated
`config.mak` outside disposable workspaces. Prebuilding a Gkeyll executable is
not required when supplying the source checkout's config. The agent must see
the config's absolute dependency paths; container mounts must preserve them.

The runtime path is compiled into Gkeyll through `GKYL_SHARE_DIR`. CI explicitly
overrides this to each build's `gkylsoft/gkeyll/share`, including when an older
prebuilt config pins it to the reference installation. Generic dependency-path
rewriting alone would point such an override into the excluded
`dependencies/gkylsoft/gkeyll` tree, despite successful ADAS copying. The unit
build installs radiation fits from the candidate or baseline source revision
into its own share directory before running tests. These fits are versioned
source data, separate from the downloaded reaction tables.

After updating the trusted `prebuilt_config.py`, start a new CI run so both
builds compile with the corrected path. Copying data or editing `config.mak`
after compilation does not change the path embedded in an existing binary.

Each build gets a config pointing to the copied libraries, with `PREFIX` and
`INSTALL_PREFIX` set to its own `gkylsoft` directory. CI uses the copied MPI
launcher and library paths; `dependencies/env.sh` saves that runtime environment
for later inspection or reruns. Source it before manually running a retained
executable. Both installations stay at their final paths because binaries embed
those paths. A cache hit skips baseline checkout, dependency copying, compilation,
unit tests, and regression creation. Candidate test lanes still run.
Invalid paths or a source config changed midway through a run fail the build.
Leave the variable unset or empty to build dependencies from machine scripts;
the baseline cache is shared across runs in either mode. Deploy the updated trusted Jenkinsfiles and helpers
to enable this option; candidate and baseline refs can predate it.

The trusted Python helper returns literal `KEY=value` lines to the Pipeline,
which reads them with sandbox-approved string operations. Prebuilt setup needs
no additional script approvals or Pipeline Utility Steps plugin. The helper's
default CLI output remains JSON; Jenkins selects `--format env`.

Regression configuration is stored at
`<installation-prefix>/gkeyll-results/runregression.config.lua`. The regression
tools derive this location from the installed executable. The trusted prebuilt
config helper also retains a literal `PREFIX=` entry for historical runners
that cannot read Make's `override PREFIX :=` syntax. Deploy the updated trusted
Jenkinsfiles to fix those older baseline refs; updating only the candidate's
runner does not change the baseline's runner.

The offline config tests use Python 3, Groovy 2.4, GNU Make, and a C compiler:

```sh
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_prebuilt_config.groovy
python3 -m unittest -v ci.jenkins.tests.test_prebuilt_config
```

The sandbox integration test uses a fresh Jenkins home with Pipeline: Job,
Pipeline: Groovy, Pipeline: Basic Steps, and Pipeline: Nodes and Processes
(including dependencies). With the WAR and plugin archives available locally:

```sh
prebuilt_test_home=$(mktemp -d)
mkdir -p "$prebuilt_test_home/plugins" "$prebuilt_test_home/init.groovy.d"
cp /path/to/test-plugin-archives/*.jpi "$prebuilt_test_home/plugins/"
cp ci/jenkins/tests/test_prebuilt_config_sandbox.groovy "$prebuilt_test_home/init.groovy.d/90-prebuilt-test.groovy"
JENKINS_HOME="$prebuilt_test_home" java \
  -Djenkins.install.runSetupWizard=false -Dgkeyll.prebuilt.test=true \
  -Dgkeyll.prebuilt.source="$PWD/ci/jenkins" \
  -jar /path/to/jenkins.war --httpPort=-1
cat "$prebuilt_test_home/prebuilt-test-result.txt"
```

This runs every helper in the default Groovy sandbox for both candidate and
baseline, with and without MPI. It uses small dependency fixtures and shuts
down its test JVM on completion. Use only a disposable controller; the offline
Groovy tests alone do not check Jenkins sandbox permissions.

The installed-tool fixture checks configure/load, installation isolation, and
Lua failure exit codes without running simulations:

```sh
GKEYLL=/path/to/gkeyll/bin/gkeyll python3 -m unittest -v ci.jenkins.tests.test_regression_config
```

## Results in GitHub

### Reporting credentials

Use the platform's existing `*_GITHUB_CREDENTIAL_ID` Jenkins **Username with
password** credential, with the GitHub username and token as its password.
An existing classic PAT with `public_repo` scope covers status and comment
publication on this public repository. A fine-grained token needs **Commit
statuses: write** and **Contents: write** (to create and update commit comments)
on `gkeyllorg/gkeyll`; its owner needs push access for statuses.
There is no additional credential or reporting plugin to configure. See GitHub's
[commit status API](https://docs.github.com/en/rest/commits/statuses) and
[commit comment API](https://docs.github.com/en/rest/commits/comments).

Each controller must use its own stable status context. Keep existing contexts
so branch protection continues to use the same names. Retries update the same
report comment for that commit and context. Delivery diagnostics are archived in
`ci-report-delivery.json`; reporting outages do not change test results.

### Queued builds and superseded PR commits

Install [queue-status.groovy](queue-status.groovy) on **each Jenkins controller**
to publish a yellow `pending` commit status as soon as a job is accepted into
its queue. This covers CLI and browser submissions on personal, Stellar CPU, and
Perlmutter GPU controllers, plus automatic and manual team-workstation builds.
The listener runs without a build executor, including when
`disableConcurrentBuilds()` blocks another run of the same job. It also covers
Pipelines waiting for their first `node()` allocation.

The pending description counts other jobs ahead, for example
`Gkeyll CI queued: waiting behind 3 other jobs (Jenkins ID 304).` With none ahead,
it says `no jobs ahead; waiting for an executor`. The count includes waiting
Gkeyll builds for the same platform on that controller in original submission
order, including the first agent wait, and excludes builds already allocated an
agent. It refreshes after queue changes and every minute; unchanged counts
are not reposted. Jenkins can skip blocked jobs or allocate multiple available
executors, so submission order does not guarantee execution order. The trailing
ID identifies the original submission for rerun tracking; it is never the count
of jobs ahead.

Once an agent is allocated, the same status refreshes every minute with the
current stage or command, elapsed build time, Jenkins' estimated remaining time,
and the build number, for example
`Gkeyll CI running: candidate unit build; elapsed 29 m; ETA ~1 hr 50 m (Jenkins build #111; ID 304).`
The controller reads Jenkins' `estimatedDuration` and build start time each
minute; elapsed time includes the Pipeline's initial agent wait. This is a
historical duration estimate, not a count of completed tests. Runs with no
estimate show `ETA unavailable`; overruns say `ETA unavailable (estimate exceeded)`.
Stage transitions retain the latest timing snapshot until the next controller
update. Unchanged descriptions are not reposted, and final results replace
progress. Stage names are shortened before timing information so the ETA remains
readable. GitHub's commit-status API uses yellow `pending` for both queued and
running work; the description distinguishes them. Queue position does not
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

#### Setup resource failures and early reports

Containerized build agents must run an init process to reap orphaned children.
See [containerized build agents](README.personal.md#containerized-build-agents)
for a reusable Podman Quadlet drop-in, installation commands, and verification.
This host configuration applies to personal and team workers; Jenkinsfile
retries cannot clear zombies owned by container PID 1.

Personal CI records pipeline selection and reporter loading before the shared
reporter is available. A Git setup command that reports process/thread exhaustion
or an allocation failure gets up to **four total attempts**, with Jenkins-managed
waits of **15, 30, and 60 seconds**, within a **10-minute timeout per setup command**.
Other failures (including an invalid ref with exit 128) fail immediately. User
cancellation and timeouts are not retried. This does not retry compiler errors or
failed tests, or change the configured compilation worker count.

Each attempt archives its output, exit code, command, elapsed milliseconds, and
best-effort resource snapshots under `ci-bootstrap/`. Records live outside the
checkout and are isolated by build number. Snapshots include process limits,
memory, visible cgroup v2 PID/memory limits and ancestors, and process thread
counts. Containers can hide ancestor limits; inspect those on the host when the
visible limits do not explain exhaustion. CPU count alone does not establish
available process/thread capacity. Repeated failures need the actual resource
limit or competing workload corrected; retries only address temporary contention.

If the Pipeline cannot publish a detailed result, the controller posts a bounded
failure comment on the candidate commit accepted at queue time and links the
terminal status to it. This works before checkout and without a working agent
shell. It preserves detailed reports and newer runs' statuses. A per-run comment
marker prevents duplicate fallback comments. Full command logs remain Jenkins
artifacts; the fallback Markdown and delivery outcome are also retained in the
controller's build directory as `gkeyll-fallback-report.md` and
`gkeyll-fallback-delivery.json`. A GitHub delivery failure is logged and does not
replace the original build failure.

Terminal statuses include seconds, for example `Failed after 7 s: Load CI
pipeline; Process/thread creation failed ...`. Reports distinguish queue wait
from execution elapsed, and list setup retry waiting as part of execution.
Failed commands retain timing just like successful commands. Shared-reporter
stage durations and outcomes are archived in `ci-stage-timings.json` and included
in the report. Older controllers without queue timestamps report the timing
that is available rather than inventing a queue duration.

Both the trusted entry Jenkinsfile and the **installed controller hook** must be
updated. Changing only the repository does not update `init.groovy.d`; follow the
installation instructions below and restart the controller while idle.

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

Every completed or failed run attempts to publish a Markdown comment on the
exact candidate commit, including PR runs. A new commit gets its own report;
reruns update the report for the same commit and status context. Stage and command
transitions update the commit status immediately; the controller refreshes
queue positions and elapsed-time/ETA indicators every minute. Progress does not
create comment notifications. Cancellation and timeout descriptions are
explicit; both use GitHub's `error` status. Reports identify the candidate,
baseline, machine, build and queue IDs, and status context, with unit/regression
results and diagnostics organized into dropdowns that are collapsed by default.
Only the overall result heading is expanded; dropdown labels show test totals
and the current or failed stage. Sections include:

- **Run details**, including commit selection, provenance, machine, and timing metadata.
- **Candidate/Baseline unit tests** and **Serial/Parallel C regressions**, with
  per-layer tables and nested failure details, including when tests fail.
- **Stage durations** and **Timings**.
- **New warnings vs main** (or the explicitly selected baseline).
- **All warnings** and **all errors**, with source log names, line numbers, and
  nearby diagnostic context.
- **Failed build log — last 100 lines**, immediately below the failed stage,
  in a collapsed dropdown preserving compiler commands, source lines, and carets.
  This excerpt preserves the original output order (including any warnings);
  extracted warnings and errors remain in their own separate dropdowns.
  Absolute and relative paths to the same log share one excerpt, with the
  command and exit code combined when available.
- An inventory of captured logs.
- Stage history, Slurm allocation details, and instructions for fetching artifacts.

Run details include UTC start/end timestamps and elapsed time. The Timings table
ends with total wall-clock time from the Jenkins build start through report
generation, including the initial agent wait. Individual timings can overlap,
so this total is not the sum of the step durations.
Displayed log excerpts wrap automatically to the available width using inline
code formatting, preserving original line breaks and indentation. Archived raw
logs retain their original lines.

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
a newer report for the same commit and context. Reports on other commits remain
intact even when an older run finishes later. The archived `ci-report.md` contains
the entire report; `ci-report.md.json` contains the report pages. Raw logs remain
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
`publish` commands remain supported. Publication requires the full tested SHA
in `--commit`; `--pr` and `--ref` are accepted for compatibility but never select
the comment target. If checkout failed before recording a SHA, keep the report
as an artifact instead of resolving a branch or PR that may have moved.
Pipelines load both shared report files before candidate checkout from reviewed
code (`main`, the pinned trusted team
checkout, or `GKEYLL_CI_TRUSTED_REF` while staging changes). Deploy those files
with the updated Jenkinsfiles. Diagnostic collection runs before binding the
GitHub token; credential-bearing shell steps are not captured in published logs.

The regression-result checker also comes from the reviewed trusted CI commit,
recorded in `ci-trusted-checker-commit.txt`. CI runs it with the baseline
executable from the baseline source directory, using `-S` because this Lua
check does not require MPI.

Run offline reporter tests with:

```sh
python3 -m unittest discover -s ci/jenkins/tests -p 'test_*.py'
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_jenkins_reporting.groovy
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_bootstrap_reporting.groovy
```

The queue listener's integration test runs on a disposable Jenkins controller
with mocked GitHub responses. It checks all four configurations, supersession
before Pipeline start and during agent wait, manual cancellation, active-build
preservation, queue persistence data, and reporting failures. It also verifies
job counts ahead and progress polling with anonymous job access disabled.
Install Credentials,
Folders, Pipeline: Job, Pipeline: Groovy, Pipeline: Basic Steps, and Pipeline:
Nodes and Processes, and Credentials Binding (including dependencies) in that
test controller. With a
Jenkins WAR and those plugin archives available locally:

```sh
queue_test_home=$(mktemp -d)
mkdir -p "$queue_test_home/plugins" "$queue_test_home/init.groovy.d"
cp /path/to/test-plugin-archives/*.jpi "$queue_test_home/plugins/"
cp ci/jenkins/tests/test_queue_status.groovy "$queue_test_home/init.groovy.d/90-queue-test.groovy"
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
LUAJIT=/path/to/luajit python3 -m unittest discover -s ci/jenkins/tests -p 'test_*.py'
```

The BGK output fixture also exercises the real CBC input for one step in serial
and on four MPI ranks, and compares the source rate and equilibrium arrays.
Set `GKEYLL` to the installed executable, `GKEYLL_CBC` to the built
`rt_gk_cbc_3x2v_p1`, and `MPIEXEC` to the matching MPI launcher, then run
`python3 -m unittest ci.jenkins.tests.test_bgk_source_io`. When using build-tree
executables, expose their shared libraries through `LD_LIBRARY_PATH`.

## Numerical regression differences

Expected numerical changes require a reviewed, new or updated entry in
`expected_regression_diffs.txt`. For a candidate/baseline comparison, CI
honors only lines that are new or changed in the candidate file relative to
the baseline file. Unchanged inherited entries are inert, so a later PR can
safely remove stale lines. Updating an inherited line's reason explicitly
acknowledges a new intentional change for that test.

For a writer fix whose old baseline array is corrupt, use the narrower form
`<full-test-name> <run-mode> <filename> baseline-array-read-failed # reason`.
This acknowledges only that file's baseline read failure, and only when the
candidate array is readable. Every failing file must have a matching entry;
numerical differences, missing or unreadable candidate files, and execution
failures still fail CI. These entries follow the same new-or-changed rule.

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
Each build keeps a separate directory, including on failure:

```text
$GKEYLL_CI_ROOT/runs/<platform>/<BUILD_TAG>/
  dependencies/                  private prebuilt libraries and saved environment
  candidate/<candidate-sha>/
    gkeyll/                      source, build, and config.mak
    gkylsoft/gkeyll/              installed Gkeyll
    gkylsoft/gkeyll-results/      databases, configs, logs, and .gkyl outputs
    _baseline/                   snapshot of baseline results and diagnostics
```

The `ci-run-path.txt` artifact records the exact directory. CI builds at these
final paths and does not delete previous run directories, including repeated
runs of the same commit. The original prebuilt installation receives no results.
CI reuses `$GKEYLL_CI_ROOT/baseline-cache/<platform>/<baseline-sha>/` in both
prebuilt and machine-script modes. It retains source/build trees, installed
Gkeyll, its dependencies and runtime environment, serial and parallel accepted
outputs, regression databases/configuration, compiler logs with exit codes,
unit results, Slurm diagnostics, and timings. Each candidate retains a snapshot
of the baseline results and diagnostics.

A matching full baseline SHA and build configuration loads the cache. Changing
the baseline builds and publishes a replacement; only successful publication
removes the previous entry. Missing or corrupt saved artifacts, failed regression
creation, or an interrupted publication cannot produce a cache hit. A manifest
is published atomically after validation and diagnostic copies finish. Old
manifest formats trigger a one-time rebuild. The GitHub report labels the
baseline `(loaded from cache)` or `(saved to cache)`; a failed attempt remains
`(cache miss; not saved)`. Restored baseline timings describe the original build,
while the report's total elapsed time describes the current run.

Compatibility includes the supplied prebuilt config content, trusted cache and
dependency helpers, configured local machine scripts/MPI paths, and
`GKEYLL_CI_CACHE_REVISION`. Change that revision on the agent after an in-place
compiler, module, SDK, or dependency upgrade (for example, `2026-10-toolchain-2`).
System toolchains and dependency contents at the original installation are not
rehashed on every hit. A platform's jobs share a filesystem lease until artifact
staging finishes, so concurrent jobs wait instead of replacing a live baseline.
The lease is released on success, failure, and normal cancellation. After an
unrecoverable controller/agent loss, remove `baseline-cache/<platform>/.lock`
only after confirming the recorded owner is no longer running. Waiting for a
lease times out after 12 hours.

Deploy the updated trusted Jenkinsfiles together with `jenkins_reporting.groovy`,
`baseline_cache.sh`, `prebuilt_config.py`, and `github_report.py`. HPC payloads
also source the baseline's saved runtime when invoking the trusted comparator.

Jenkins archives the complete regression result trees for candidate and
baseline, including `.gkyl` outputs, databases, configuration, and logs. The
temporary Jenkins workspace is removed after archival; the run directory
remains. Accepted baseline outputs are copied into the candidate's result tree
so inspecting an old run does not depend on a surviving shared-cache symlink.
Jenkins builds and artifacts have no automatic count limit by default. Set
`GKEYLL_CI_BUILDS_TO_KEEP` to a positive build count to opt into Jenkins retention
limits; this does not delete directories under `runs/`. Remove retained run
directories explicitly when their data is no longer needed.

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
