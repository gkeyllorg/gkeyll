# GitHub reporting

[CI overview](README.md)

Each machine publishes commit statuses for queued, running, and completed CI.
Completed reports are comments on the exact tested commit, including PR runs,
linked from the status's **Details** button. Reruns update the report for that
commit and status context. Reporting outages do not change test results;
`ci-report-delivery.json` records delivery failures.

## Reporting credentials

Use the platform's existing `*_GITHUB_CREDENTIAL_ID` **Username with password**
credential: GitHub username and token as password. A classic PAT needs
`public_repo`; a fine-grained token needs **Commit statuses: write** and
**Contents: write** on `gkeyllorg/gkeyll`. The owner needs push access for
statuses. See GitHub's [status API](https://docs.github.com/en/rest/commits/statuses)
and [comment API](https://docs.github.com/en/rest/commits/comments).

No extra credential or reporting plugin is required. Keep each controller's
status context stable so branch protection and report updates retain their
existing names. Public source discovery and checkouts remain anonymous.

## Controller installation

Install [queue-status.groovy](queue-status.groovy) on each controller and deploy
the trusted Jenkinsfiles with `jenkins_reporting.groovy` and `github_report.py`.
The hook uses the existing **Credentials** plugin. Set credential IDs in Jenkins
**global environment variables**; node-only variables are unavailable before
agent allocation. Keep the team GitHub source anonymous.

Set `PERSONAL_STATUS_CONTEXT` explicitly, for example to
`continuous-integration/jenkins/personal-mycomputer`; the controller cannot infer
an agent hostname while queued. Personal and team status contexts must match
their Pipeline configurations.

From a reviewed checkout, as the Jenkins service owner:

```sh
install -d -m 700 "$JENKINS_HOME/init.groovy.d"
install -m 600 ci/jenkins/queue-status.groovy \
  "$JENKINS_HOME/init.groovy.d/gkeyll-queue-status.groovy"
```

Use the actual controller home: typically `/var/lib/jenkins` for a workstation
service or `$GKEYLL_CI_ROOT/jenkins_home` on HPC. System services may need an
administrator to install the file with Jenkins service ownership. Restart while
idle and check for `Gkeyll GitHub queue status listener installed` in the log.
Install once; do not also run the hook in Script Console. Repository updates
do not update this installed copy. To remove it, delete the hook and restart;
Pipeline running/final reporting remains available.

For renamed or staging jobs, set the full job name (including folders) globally:

| Platform | Job-name override |
| --- | --- |
| Personal | `GKEYLL_CI_QUEUE_PERSONAL_JOB` |
| Stellar CPU | `GKEYLL_CI_QUEUE_STELLAR_CPU_JOB` |
| Perlmutter GPU | `GKEYLL_CI_QUEUE_PERLMUTTER_GPU_JOB` |
| Team workstation | `GKEYLL_CI_QUEUE_TEAM_WORKSTATION_JOB` (multibranch parent) |

Defaults match the platform guides. The hook handles only configured jobs and
direct children of the team parent. Queue metadata persists in
`$JENKINS_HOME/gkeyll-queue-status/` across restarts and becomes internal
`CI_QUEUE_*` build parameters without disrupting queue deduplication. Do not
add these parameters to the browser form.

To verify, submit a PR with its agent unavailable and check that a yellow status
appears before checkout. Push a newer commit and check that the old queued run
is cancelled on the next successful sweep. Also check manual cancellation and
a normal running/completed build. [Integration tests](tests/README.md) use a
disposable controller.

## Queue and progress

The listener publishes yellow `pending` when Jenkins accepts a job, without an
executor. This covers CLI, browser, and automatic submissions, jobs blocked by
`disableConcurrentBuilds()`, and the first `node()` wait. It resolves and pins
the candidate SHA; the Pipeline checks out that exact commit.

Queue descriptions count earlier waiting Gkeyll jobs for the same platform on
that controller, including initial agent waits and excluding allocated builds.
Counts refresh on queue changes and every minute. The trailing Jenkins ID is
the original submission ID, not the count. Blocked jobs and multiple executors
can change execution order; queue position does not predict start time.

After allocation, statuses show the stage/command, elapsed build time, Jenkins'
estimated remaining time, build number, and submission ID. Stage changes update
immediately; timing refreshes every minute, retaining the previous snapshot
between updates. Elapsed time includes the initial agent wait. ETA uses historical
`estimatedDuration`, with `unavailable` or `estimate exceeded` when appropriate.
Long stage names are shortened to leave timing readable. Unchanged descriptions
are not reposted; progress creates no report-comment notifications.

For PR selections, the listener checks for newer heads after submissions and
every minute while queued. It cancels superseded runs before their first agent
allocation and marks the old commit `error` with the replacement SHA. It leaves
active runs, other PRs, repeat requests for the same commit, and branch/SHA
selections alone. It does not submit replacements: team discovery supplies new
runs; manual platforms require resubmission.

Cancellation, timeout, and early failures close pending statuses using `error`
and an explanation. Cleanup preserves newer statuses and detailed Pipeline
results. GitHub outages do not block execution, but notifications and
supersession checks need a working API and valid credentials.

## Completed reports

Reports use collapsed sections for:

- Run details: candidate/baseline, machine, build/queue IDs, status context,
  trusted code commits, UTC timestamps, and queue/execution times.
- Candidate/baseline unit tests and serial/parallel C regressions, with layer
  totals and failure details.
- Stage durations, timings, the 20 slowest candidate regressions, and regression
  timing changes.
- New warnings relative to the selected baseline, all warnings/errors with log
  locations and context, and the failed build log's last 100 lines.
- Stage history, Slurm details, and artifact retrieval instructions.

The failed-log excerpt preserves output order, compiler commands, source lines,
and carets; duplicate paths to one log share an excerpt. Displayed logs wrap
while preserving line breaks and indentation; archived logs remain raw.

Warning comparison uses matching completed baseline steps, ignoring checkout
paths, ANSI colors, and source line/column changes. Missing or failed steps show
comparison as unavailable; unmatched/shared warnings remain in **All warnings**.
HPC builds baseline unit tests too, providing corresponding compiler diagnostics.

Per-test timing uses SQLite `runtime` from the latest finalized invocation per
suite/layer, excluding compilation and comparison. Incomplete, invalid, or
ambiguous observations are omitted. The slowest list includes numerical
failures; timing comparisons require passed/baseline-created results matched
by suite/layer, test name, type, and CPU/GPU serial/parallel mode. They show
baseline/candidate seconds, signed differences, and ratios when both runs took
at least 5 seconds and differ by at least 2x. Missing data is unavailable, not
zero. Timing observations do not affect CI status.

These are single-run measurements affected by worker counts, GPU sharing,
machine load, I/O, and cache age; cached baseline timings are labeled. Confirm
performance changes with repeated, interleaved runs on the same hardware and
fixed settings/load, comparing medians and variability. Changed step counts or
output frequency also affect total cost. Report total time runs from Jenkins
build start through report generation; overlapping steps need not sum to it.

Large reports continue in linked comments, published before updating the main
comment. Only comments owned by the authenticated account with the exact
context marker are edited; obsolete continuation pages are cleared. Queue IDs
(or build numbers within a job) prevent older runs from replacing newer reports
on the same commit/context. Other commits' reports remain intact.

Artifacts include `ci-report.md` (full report), `ci-report.md.json` (pages),
`ci-stage-history.txt` (UTC transitions, including before checkout), and
`ci-stage-timings.json`. The captured-log inventory is artifact-only, in
`ci-diagnostic-summary.json` as `captured_log_paths` and count `captured_logs`.

## Setup resource failures and early reports

Containerized personal/team agents need an init process to reap orphaned
children; see [container setup](README.personal.md#containerized-build-agents).
Pipeline retries cannot clear zombies owned by container PID 1.

Personal CI records pipeline selection and reporter loading before the shared
reporter is available. Git setup commands reporting process/thread exhaustion
or allocation failure get up to four attempts, with 15, 30, and 60 second waits
and a 10-minute timeout per setup command. Other errors, cancellation, and
timeouts are not retried. This does not retry builds/tests or change worker counts.

Per-attempt logs, commands, exit codes, elapsed milliseconds, and resource
snapshots live under `ci-bootstrap/`, outside checkout and isolated by build
number. Snapshots include process/thread counts, process limits, memory, and
visible cgroup v2 limits/ancestors. Inspect host limits when containers hide
them; repeated failures require correcting limits or competing load.

If the Pipeline cannot report, the controller posts a bounded failure comment
on the queue-pinned commit and links the terminal status to it, even before
checkout or without an agent shell. A per-run marker prevents duplicates;
detailed reports and newer statuses are preserved. The controller's build
directory retains `gkeyll-fallback-report.md` and `gkeyll-fallback-delivery.json`;
full logs remain artifacts. Delivery failure does not replace the build failure.
Terminal statuses include seconds; reports separate queue wait from execution
and include setup retry waits and failed-command timings. Older controllers
report only available timing data. Update both the trusted entry Jenkinsfile
and installed hook to enable this behavior.

## Trusted code and delivery

Pipelines load `github_report.py` and `jenkins_reporting.groovy` before candidate
checkout from reviewed code: `main`, the pinned trusted team checkout, or
`GKEYLL_CI_TRUSTED_REF` for staging. Personal `CI_REF` independently selects its
implementation; see [personal CLI selection](README.personal.md#cli-launch).
The loader logs and archives its resolved commit. Per-run selections do not
replace the installed controller hook. Diagnostics are collected before binding
the GitHub token; credential-bearing shell steps are not published.

The regression checker comes from the trusted CI commit recorded in
`ci-trusted-checker-commit.txt`. It runs with the baseline executable from the
baseline source directory, using `-S` because this Lua check does not need MPI.

Reporting runs during failure cleanup, even if timing summaries fail. Without
an exact candidate SHA, reports remain artifacts. The reporter's standalone
`build` and `publish` commands remain supported; publication requires the full
tested SHA in `--commit`. Compatibility options `--pr` and `--ref` never select
the comment target. Do not resolve moving refs after a checkout failure.

Network/timeouts and transient server errors get up to three attempts. Comment
retries discover existing pages to recover lost POST responses; status and
comment publication are attempted independently. Late progress cannot overwrite
terminal results. Persistent outages or unavailable agents can still prevent
delivery; inspect the build log, archived report, and delivery diagnostics.
