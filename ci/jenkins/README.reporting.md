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

All main report content must be inside collapsed dropdowns, including content
on continuation pages. Only the report heading and page navigation remain
outside. Prefer sibling sections to keep the report shallow; avoid adding
wrapper dropdowns around existing groups or dropdowns for individual regression
test previews.

Reports use collapsed sections for:

- Run details: candidate/baseline, machine, build/queue IDs, status context,
  trusted code commits, UTC timestamps, and queue/execution times.
- Candidate/baseline unit tests, with layer totals and failure details.
- Baseline/candidate regression execution failures, near the top of failed
  reports: test names, CPU/GPU and serial/parallel mode, failure status, fatal
  messages with line numbers, and the database containing the full output.
  These come directly from the latest finalized regression invocation, even
  when the runner exits zero or no failed command log was captured. Baseline
  failures explain why the results cannot be cached. This reports observed
  failures; identifying the responsible source change still requires diagnosis.
- Serial/parallel regression totals, layer counts, candidate-only tests, and
  expected-diff entries whose tests pass, in a summary dropdown for each suite.
- Acknowledged and unacknowledged regression tests, grouped separately for each
  serial/parallel suite, with every failing test's status and differing or
  missing files. These groups are siblings of the regression summary dropdown;
  individual test previews do not have dropdowns.
- Stage durations, timings, the 20 slowest candidate regressions, and regression
  timing changes (only when differences meet the reporting thresholds).
- Candidate regressions outputting more than one frame after initial-condition
  frame 0, only when any are found. Counts use distinct numbered `.gkyl` outputs
  across fields and species, so frames 0 and 1 alone do not trigger this section
  and tests without field output are included. Staging preserves a filename
  inventory for this check alongside the diagnostic artifacts.
- New/all compiler and other warnings, separate new/all time-step warnings,
  all errors with log locations and context, and the failed build log's last
  100 lines.
- Stage history and Slurm details.

The unacknowledged group includes copyable lines for
`ci/jenkins/expected_regression_diffs.txt` and explains how to add a reason,
refresh inherited acknowledgments, or narrow them to a run mode or file.
Each test summary shows its status, differing/missing file counts, the largest
absolute and relative differences with their filenames, and counts by comparison
reason (for example, numerical differences or incompatible dynvec outputs).
Maxima cover every file, including files not named in the summary. At most four
reason types are shown, with an explicit count of additional types. These are
observed comparison failures, not inferred physical root causes. Accepted
differences are explicitly identified as not failing CI.

Complete file lists and original comparison details are retained in
`ci-report-regressions.md` in Jenkins artifacts, rather than repeated in GitHub
continuation comments. Regression-check command tails are replaced by log paths
and exit codes when a checker summary is available; incomplete checks and other
failed commands retain their tails. Raw logs remain archived. Timing-data
omissions are summarized by candidate/baseline diagnostic counts on GitHub;
every diagnostic is preserved in `ci-report-timing-details.md`. Existing
`ci-report*.md` archive patterns include both new artifacts on all platforms.

Deploy both `github_report.py` and `check_regression_results.lua` in the trusted
CI revision for complete acknowledged-file details. Older summaries remain
readable, but details that the old checker omitted cannot be recovered from
those summaries.

The failed-log excerpt preserves output order, compiler commands, source lines,
and carets; duplicate paths to one log share an excerpt. Displayed logs wrap
while preserving line breaks and indentation; archived logs remain raw.

Warning comparison uses matching completed baseline steps, ignoring checkout
paths, ANSI colors, and source line/column changes. Missing or failed steps show
comparison as unavailable; unmatched/shared warnings remain in their **All** section.
HPC builds baseline unit tests too, providing corresponding compiler diagnostics.

Compiler warnings whose source path contains a `ker/` directory component are
excluded from warning sections, comparisons, counts, and `ci-report-warnings.md`.
This applies to all layers and both candidate and baseline builds. Kernel errors
are still reported, and raw logs (including failure excerpts) remain unfiltered.
Only identifiable source locations are filtered; diagnostics whose paths were
destroyed by interleaved compiler output remain visible rather than guessing
their origin. Paths mentioned in messages or include context do not suppress
warnings originating outside `ker/`.

Warnings tagged `[-Wnan-infinity-disabled]` are omitted from GitHub warning
sections and the counts in `ci-diagnostic-summary.json`, including new and
unclassified warning counts. They remain in `ci-report-warnings.md` (unless
excluded by the kernel rule above) and raw logs. This reporting filter does not
change compiler flags or suppress runtime NaN/infinity diagnostics or errors.

Warnings are split into **Compiler and other warnings** and **Time-step warnings**,
each with **New ... vs baseline** and **All ...** collapsed sections. Small-time-step
warnings (`WARNING: Time-step ...`) appear only in the time-step sections, including
messages with nonfinite time steps. Actual abort diagnostics (`ERROR: Time-step ...`)
remain in **All errors**.

**All** sections include occurrences across candidate, baseline, and shared logs,
without a preview limit. **New** sections retain distinct new messages per matching
build/test step. Repeated time-step warnings are condensed to one entry per log/test
and threshold, showing the line range, time-step range, and maximum `num_failures`.
Counts describe occurrences in **All** sections and distinct new messages in **New**
sections. Grouping changes display only; baseline comparison and aggregate warning
counts in `ci-diagnostic-summary.json` stay the same. The complete new/all warning
sections, split into the same categories and including every nonexcluded occurrence
and its context, are retained in `ci-report-warnings.md` in Jenkins artifacts.
Many compiler warnings or errors can still produce continuation comments; none
are omitted or truncated to shorten the report.

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
context marker are edited; obsolete continuation pages are deleted after the
replacement report is published. Queue IDs
(or build numbers within a job) prevent older runs from replacing newer reports
on the same commit/context. Other commits' reports remain intact.

Artifacts include `ci-report.md` (posted report), `ci-report.md.json` (pages),
`ci-report-warnings.md` (complete new/all warning details),
`ci-report-regressions.md` (complete regression comparisons),
`ci-report-timing-details.md` (timing-data diagnostics),
`ci-stage-history.txt` (UTC transitions, including before checkout), and
`ci-stage-timings.json`. The captured-log inventory is artifact-only, in
`ci-diagnostic-summary.json` as `captured_log_paths` and count `captured_logs`.

Small timing, command-exit, Slurm job-ID/status, regression/failure-summary,
provenance, cache, and report-delivery files are consolidated into
`ci-statistics.md` before archiving. The original files remain in the workspace
for CI and cache readers, but are excluded from Jenkins artifacts only after
consolidation succeeds. Full logs, result databases, and detailed report JSON
remain separate artifacts.

Deploy both `jenkins_reporting.groovy` and `github_report.py` in the trusted
reporting revision to enable consolidation; updating only the candidate checkout
is insufficient. If the trusted reporter is older or consolidation fails, Jenkins
logs a warning and archives the original files to preserve diagnostics. Existing
build artifacts are not rewritten by this change.

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
