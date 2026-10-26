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
`team-workstation`, `stellar-cpu`, `stellar-ai-cpu`, and `perlmutter-gpu`. Builds use their final
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

Personal, team workstation, Stellar CPU, and Perlmutter use `cache_manager.py`
to coordinate concurrent runs under a filesystem lock. A run claims its baseline
until diagnostics and numerical snapshots have been staged. Builds waiting for
the same baseline retry for up to 12 hours; incompatible configurations wait
until existing readers finish before rebuilding at the fixed install path.
Cleanup checks Jenkins build state and surviving Slurm jobs before deleting
inactive workspaces or baselines. Failed API/Slurm queries preserve files and
fail cleanup. Failed baseline builds are removed after diagnostics are collected.
Ownership records live under `cache-runs/`; trusted helpers use `cache-control/`.

Agents need Jenkins API read access through `CI_JENKINS_READ_AUTH_FILE`, a
mode-600 `user:api-token` file. An explicitly empty value enables anonymous
reads. Cluster agents can set `CI_JENKINS_LOCAL_URL` (default
`http://127.0.0.1:8080`) to the controller URL reachable from the login node.
Shared storage must support `flock`; concurrent builds need enough executors.
Drain old CI runs before deployment, and remove a legacy
`baseline-cache/<platform>/.lock` only after its owner has stopped.
Stellar AI retains the existing platform-wide lease workflow.

Deploy trusted Jenkinsfiles with `jenkins_reporting.groovy`, `baseline_cache.sh`,
`prebuilt_config.py`, `cache_manager.py`, and `github_report.py`. HPC payloads also source the
baseline's saved runtime when invoking the trusted comparator.

## Retention

Jenkins retains all candidate and baseline regression results, then removes its
temporary workspace. By default the two newest builds per job retain their
entire run directories, including source, build objects, installed binaries,
run-local dependencies, numerical outputs, and independent baseline snapshots.
Active runs are preserved until they finish.

Databases, configs, logs, and other non-`.gkyl` files remain individually
downloadable at their existing artifact paths. The staging helper copies these
files in bulk. It packages all `.gkyl` files directly from the retained trees
into two artifacts under `ci-numerical/`: `candidate.tar.zst` and
`baseline.tar.zst`. Each contains paths rooted at `gkylsoft/gkeyll-results/`,
including serial and parallel results. Symlinks and special files are excluded
from both staging and packaging. There is no filtering of successful tests.

Packaging uses zstd level 1 with the serial regression worker budget when
`zstd` is on the agent's PATH; otherwise it uses parallel pigz or serial gzip
level 1 and the names end in `.tar.gz`.
Install `zstd` on agents for faster compression. A missing results tree (for
example, after an early build failure) produces no numerical archive. The helper
logs staging time, packaging time, and compressed sizes separately.

For example, download and unpack a personal build's numerical outputs:

```sh
./ci/jenkins/gkeyll-ci.sh personal artifact --build 200 --fetch \
  --only ci-numerical/candidate.tar.zst,ci-numerical/baseline.tar.zst \
  --output-dir results-200
cd results-200
mkdir -p _baseline
zstd -dc ci-numerical/candidate.tar.zst | tar -xf -
zstd -dc ci-numerical/baseline.tar.zst | tar -xf - -C _baseline
```

Use the build's actual archive names (`artifact --build NUMBER --list` lists
them). For `.tar.gz` archives, substitute `gzip -dc` and `.tar.gz` above. Builds
created before this packaging change retain their original loose `.gkyl` files.

Personal, team workstation, Stellar CPU, Stellar AI CPU, and Perlmutter default
to retaining two Jenkins builds, their artifacts, and their **entire run
directories**. `GKEYLL_CI_BUILDS_TO_KEEP` sets this positive count for both
Jenkins and disk retention (`--runs-to-keep` in the cache manager).
Jenkins may preserve explicitly protected builds in addition to this count;
that protection does not pin their external run directories.

With the default `GKEYLL_CI_RETAIN_RUNS=0`, startup reserves a slot for the
incoming build and deletes older inactive run directories. Completion retries
cleanup. Runs are ordered by build number within each job, not modification
time. Deletion removes the whole directory, including `dependencies/`,
`candidate/`, and any `run.tar.*` archive. Jenkins and Slurm activity checks
protect active runs, which can temporarily exceed the limit. Query failures
preserve files and fail cleanup. Keep the Jenkins read credentials configured.

The next updated build also prunes existing directories whose tags identify
older builds of the same job, even when Jenkins has discarded their records.
Unknown owners are never deleted automatically. Managed platforms stop on
unknown owners; Stellar AI only considers tags belonging to the current job.
Stellar AI prunes while holding its platform lease; completion pruning requires
successful artifact archival and Slurm finalization. Deploy the Jenkinsfile and
`cache_manager.py` together through the trusted CI revision.

`GKEYLL_CI_RETAIN_RUNS=1` explicitly disables run-directory pruning and retains
unlimited history, overriding the disk count. Leave it unset or `0` to bound
disk usage. Jenkins history remains bounded by `GKEYLL_CI_BUILDS_TO_KEEP`.
The reusable baseline cache independently retains one valid entry, deleting the
outdated entry after its replacement is successfully published; entries still
used by active builds are removed once those builds finish.
Use the [artifact command](README.md#unified-local-command) to retrieve
archived files.

## Optional run compaction

Run compaction is disabled by default. Set `GKEYLL_CI_COMPACT_RUNS=1` to
compact finalized runs after artifact archival. Personal, team workstation,
Stellar, and Perlmutter produce `runs/<platform>/<BUILD_TAG>/run.tar.zst`
(or `run.tar.gz`), containing the candidate tree and any run-local dependencies.
This includes source, build objects, installed binaries, diagnostics, numerical
outputs, and the independent baseline snapshot. Cluster pipelines require
successful Slurm cleanup before compaction.

The helper compares the archive against the original files before publishing
it and removing the loose trees. Creation or verification failures preserve
the originals. Compaction requires GNU tar (`tar`, or `gtar` on macOS), and
uses zstd level 1, falling back to pigz level 1 and then serial gzip level 1.
Both compaction and numerical artifact packaging use the serial regression
worker count: `STELLAR_CPU_REGRESSION_JOBS` (default 16),
`PERLMUTTER_GPU_REGRESSION_JOBS` (default 4), or `PERSONAL_REGRESSION_JOBS`
(defaulting to `PERSONAL_BUILD_JOBS`, then 3). Team workstation settings are
forwarded to the personal workflow. Manual helper invocations accept
`CI_COMPRESSION_JOBS` (default 2).

Compaction needs room for the archive alongside the original files. It reduces
retained file counts, not peak compilation/test storage. The baseline cache
remains unpacked. With `GKEYLL_CI_RETAIN_RUNS=1`, run archives have no automatic count limit;
Jenkins build retention does not delete them. Bounded run retention (`0`)
removes inactive archived runs under the same policy as unpacked runs.

Restore in the original run directory to preserve embedded absolute paths:

```sh
cd "$GKEYLL_CI_ROOT/runs/stellar-cpu/<BUILD_TAG>"
zstd -dc run.tar.zst | tar -xf -
# For the gzip fallback: tar -xzf run.tar.gz
```

An existing finalized run can also be compacted explicitly after confirming
Jenkins and Slurm have finished using it:

```sh
BUILD_TAG=jenkins-gkeyll-ci-stellar_cpu-9 \
  ci/jenkins/baseline_cache.sh compact-run "$GKEYLL_CI_ROOT" stellar-cpu
```

The helper refuses unfinished manifests and another build's active cache lease.
An interrupted compaction can leave `.compact-lock` or `run.tar.*.tmp`; remove
these only when no compactor is running. If both the published archive and
payload remain, verify the archive before removing either.
