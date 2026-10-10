# CI storage and baseline caching

[CI overview](README.md)

Every controller requires a writable, agent-visible `GKEYLL_CI_ROOT`. Each
build keeps a separate directory, including failures and reruns of the same SHA:

```text
$GKEYLL_CI_ROOT/runs/<platform>/<BUILD_TAG>/
  dependencies/                  prebuilt libraries and saved environment
  candidate/<candidate-sha>/
    gkeyll/                      source, build, and config.mak
    gkylsoft/gkeyll/              installed executable and libraries
    gkylsoft/gkeyll-results/      databases, configs, logs, and .gkyl outputs
    _baseline/                   baseline results and diagnostics snapshot
```

`ci-run-path.txt` records the exact path. Platform names are `personal`,
`team-workstation`, `stellar-cpu`, and `perlmutter-gpu`. Builds use their final
paths so embedded library and data paths remain valid. The original prebuilt
installation receives no results.

## Baseline cache

Both dependency modes reuse
`$GKEYLL_CI_ROOT/baseline-cache/<platform>/<baseline-sha>/`. The cache retains
source/build trees, installed Gkeyll, dependencies and runtime environment,
serial/parallel accepted outputs, regression databases/configs, compiler logs
and exit codes, unit results, Slurm diagnostics, and timings.

A hit requires the same full baseline SHA and compatible build configuration.
It skips baseline checkout, dependency copying, compilation, unit tests, and
regression creation. Candidate builds and tests still run from scratch.
Changed baselines replace the previous entry only after successful publication.
Missing or corrupt artifacts, failed regression creation, and interrupted
publication cannot produce a hit; the manifest is published atomically after
validation and diagnostic copying. Old manifest formats rebuild once.

Compatibility includes prebuilt config content, trusted cache/dependency
helpers, local machine scripts/MPI paths, and `GKEYLL_CI_CACHE_REVISION`.
Change that revision on the agent after an in-place compiler, module, SDK, or
dependency upgrade, for example to `2026-10-toolchain-2`. System toolchains and
original dependency contents are not rehashed on every hit.

Reports label baselines `loaded from cache`, `saved to cache`, or
`cache miss; not saved`. On a cache hit, `baseline_*_seconds` in
`ci-timing-summary.txt` are zero for skipped work. Historical measurements use
`cached_baseline_*_seconds` and appear separately in the report; they do not
contribute to current-run compile totals. Restored raw `baseline-*-seconds.txt`
files describe the original build. Total elapsed time describes the current run.
Each candidate retains baseline
diagnostics and copies of accepted outputs, so old runs do not depend on a
surviving cache symlink.

Jobs for one platform share a filesystem lease until artifact staging finishes.
Concurrent jobs wait up to 12 hours; success, failure, and normal cancellation
release the lease. After controller/agent loss, remove
`baseline-cache/<platform>/.lock` only after confirming its recorded owner is
no longer running.

Deploy trusted Jenkinsfiles with `jenkins_reporting.groovy`, `baseline_cache.sh`,
`prebuilt_config.py`, and `github_report.py`. HPC payloads also source the
baseline's saved runtime when invoking the trusted comparator.

## Retention

Jenkins archives complete candidate and baseline regression trees, including
`.gkyl` outputs, databases, configs, and logs, then removes its temporary
workspace. Run directories remain.

Jenkins builds and artifacts have no automatic count limit by default. Set
`GKEYLL_CI_BUILDS_TO_KEEP` to a positive count to limit Jenkins retention.
This does not remove `runs/` directories; delete those explicitly when no longer
needed. Use the [artifact command](README.md#unified-local-command) to retrieve
archived files.
