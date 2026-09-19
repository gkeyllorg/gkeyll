# Jenkins workstation handoff — 2026-09-17, 14:08 UTC

## User instructions and current result

The user wants this RTX 3090 workstation shared by all Gkeyll contributors,
using the existing team-workstation PR discovery workflow. The immediate
validation target is PR 1132 (`agent_tools`). The user explicitly instructed:
**do not commit anything to this branch**, then asked to conclude and leave
a detailed handoff because more contributions are coming soon.

No commits, pushes, or PRs were created. HEAD remains
`91ef9ed70d4636ce8bb4da32ed5fa7e4a1eae0b5`. Preserve all working-tree edits:
several provisioning files and personal-CI fixes already existed at session
start, and this session extended them. Solver source files were not edited.

The shared node is operational, but PR 1132 is **not passing the full CI**.
Its shared Jenkins build #2 finished with FAILURE during candidate unit-test
linking. Jenkins and the GPU agent remain active; the queue was empty at
handoff. Automatic two-minute discovery remains enabled. No build was aborted
in response to the user's conclusion request; the build completed naturally.

## Immediate source failure to investigate

Build: <https://bazzite.taile273b2.ts.net/job/gkeyll-ci-team-workstation/job/PR-1132/2/>

Candidate: `91ef9ed70d4636ce8bb4da32ed5fa7e4a1eae0b5`.
Selected baseline: `main`, which resolved earlier to
`0e4377b7687721d1a9cc8c2fc4425add98a43514`; the baseline build was never reached.

The dependency scripts completed OpenMPI, LuaJIT, OpenBLAS, and SuperLU.
The subsequent CUDA unit build failed linking
`cuda-build/core/unit/ctest_linsolvers` with undefined references to:

```text
test_cusolver_qr_dev
test_cusolver_rf_dev
test_cusolver_ops_dev
test_cusolver_ops_multiple_rhs_dev
test_cusolver_ops_multiple_prob_dev
```

The naming mismatch is visible directly in the checked-out source:

- `core/unit/ctest_linsolvers.c`: declarations near line 10 and registrations
  near line 632 use `_dev` suffixes.
- `core/unit/ctest_cusolver.cu`: the `extern "C"` declarations near line 18
  and function definitions use the corresponding names **without** `_dev`.

No fix was attempted after the user requested conclusion. Start with a
focused build/test of this target before spending time on another full run.
The full unit suite, baseline build, and regression comparison remain
unverified. The earlier CUDA launch-delimiter issue is fixed sufficiently for
the focused CUDA tests and the core library to compile.

## Services, access, and credentials

- Controller service/container: `gkeyll-jenkins`.
- GPU agent service/container: `gkeyll-gpu-agent`.
- Node: `bazzite-rtx3090`; one executor, exclusive labels
  `linux bazzite cuda rtx3090`; controller has zero executors.
- Shared job: `gkeyll-ci-team-workstation`, with `main` and `PR-*` children.
- Old `gkeyll-ci-personal` job is retained.
- Jenkins: <https://bazzite.taile273b2.ts.net/> or local port 8080.
- Cockpit: <https://bazzite.taile273b2.ts.net:8443/>.
- These HTTPS endpoints are tailnet-only. Contributor PR discovery is
  independent of dashboard access; colleagues still need appropriate
  Tailscale access and individual Jenkins accounts to use the dashboard/CLI.
- User lingering is enabled and suspend targets are masked. Lock/display-off
  allows CI to continue. Actual suspend does not.
- Controller and agent survived a controller restart during verification.
  Remote access from a second physical device was not tested in this session.

Jenkins login is `mrosen`. Password file:
`~/.local/share/gkeyll-jenkins/admin-password`.
CLI authentication file: `~/.config/gkeyll/jenkins-cli.auth`, mode 600.
Do not print either credential into logs or conversation.

GitHub credential ID: `gkeyll-github`. The original fine-grained token failed
because its lifetime exceeded the organization's 366-day maximum. After
lifetime/permission adjustments it could read PRs, but status writes still
returned HTTP 403. The user replaced it with a classic token; PR reads then
returned HTTP 200 and commit-status creation HTTP 201. Build #2 successfully
published its pending and failure statuses.

**Outstanding credential issue:** the successful classic token's actual
`X-OAuth-Scopes` response included broad repository and administration scopes,
not just `repo:status`. The user was told to narrow it. Narrowing was not
confirmed before handoff. Inspect scope metadata without exposing the token.

Matrix authorization now reserves administrator access for `mrosen` and
grants authenticated colleague accounts read/build/cancel permissions.
Self-registration and anonymous Jenkins access remain disabled.

## What is deployed

The controller uses the standard Pipeline: Multibranch with defaults and
Config File Provider plugins to run a sandboxed, controller-managed snapshot
of this repository's team and personal Pipelines. This is needed because
`main` does not yet contain `ci/jenkins`: existing PRs should not need special
CI files, and their source should not replace the controlling Pipeline.

The shared scan discovered `main` and all 24 open PRs. Twenty-one target
`main`; the wrapper skips the other targets. Discovery uses the GitHub
credential for API access because anonymous discovery was throttled by the
60-request/hour quota. Candidate/baseline Git checkouts remain anonymous.

The actual build uses existing machine scripts:

```text
TEAM_WORKSTATION_MKDEPS_SCRIPT=mkdeps.linux.sh
TEAM_WORKSTATION_CONFIGURE_SCRIPT=configure.linux.gpu.sh
TEAM_WORKSTATION_BUILD_JOBS=8
TEAM_WORKSTATION_REGRESSION_JOBS=1
MAKEFLAGS=CUDA_ARCH=86
```

MPI/NCCL/cuDSS are not enabled by the configure script, though the dependency
script builds OpenMPI. Dependency scripts have their own hard-coded build
parallelism, independent of the pipeline's eight jobs.

CUDA 12.6, Compute Sanitizer 2024.3.0, and Valgrind 3.22.0 are installed in the
persistent GPU-agent image. The tools are available via `podman exec`; they
were not installed into Bazzite's immutable host filesystem. The default
Jenkins pipeline does not automatically run memory checkers.

## Uncommitted changes and deployment distinction

Repository CLI/Pipeline fixes made or retained locally include:

- Linux/macOS credential-file mode check order in personal/team clients.
- Queue parsing that preserves an empty build-number field.
- Team scan support for Jenkins' HTTP 302 folder redirect and indexing console
  endpoint; the assumed computation JSON API does not exist here.
- Team wrapper uses `evaluate(readTrusted(...))` without requiring an outer
  executor/workspace for `load`.
- Sequential unit compilation/execution, also needed for the older baseline's
  parallel `make check` dependency ordering.
- Trusted regression checker available to candidate branches predating CI.
- Selection metadata written after checkout cleanup, CUDA log archiving,
  useful GitHub error bodies, and workspace cleanup after status failures.
- Provisioning files now configure shared discovery, restricted colleague
  permissions, and persistent memory-checking tools.
- Updated personal/team/Bazzite setup documentation.

**One important change is local and tested but NOT deployed:**
`ci/jenkins/check_regression_results.lua` now also checks `gpu_status` and
`cpu_gpu_diff`, and supports legacy databases without these columns. The
old checker ignored GPU failures. Ten fixtures passed for CPU/GPU success,
failure, GPU crash/timeout, CPU-GPU mismatch, acknowledged differences, and
legacy schemas. The current controller image/snapshot still contains the
older checker. Review this change, then rebuild/redeploy the controller when
idle. The existing acknowledgment-file semantics are preserved; this was not
a complete audit of missing/empty databases or acknowledged failure policy.

Do not assume local file edits are automatically active in Jenkins: the
controller image and `bootstrap.groovy` create the managed pipeline snapshot.
`install-user.sh` rebuilds both images and restarts both services, so check
for running work first. No additional plugin is needed to deploy the checker.

## Verification and retained evidence

At the exact PR 1132 commit above:

- Array/matrix CUDA tests: **35 cases passed, including 9 device cases**.
- Compute Sanitizer: both executables reported **0 errors, 0 leaked bytes**.
- Separate CPU array/matrix tests under Valgrind: **26 cases passed**, zero
  errors, and all heap allocations freed.
- Seven CLI fixtures passed, including waiting/cancelled queues, mode 644
  rejection, and the real multibranch scan response shape.
- Ten regression-checker database fixtures passed.
- `bash -n ci/jenkins/*.sh` and `git diff --check` passed.

The focused diagnostic build used the agent image's OpenBLAS/LAPACKE and
SuperLU. Ubuntu splits LAPACKE into a separate library, so it required
`LAPACK_LIBS='-llapacke -lopenblas'`; this does not alter the Jenkins builds,
which build repository dependencies themselves.

For Valgrind on this Ryzen 7700X, use a separate CPU build configured with
`ARCH_FLAGS=-march=x86-64-v3` and `OPENBLAS_CORETYPE=Haswell`; Valgrind 3.22
could not decode AVX-512 from the initial native build. Use `--no-exec` for
the unit executable and `--error-exitcode=1` for the checker.

Diagnostic checkout/binaries remain inside the agent at
`/home/jenkins/agent/pr1132-verification` (host path
`~/.local/share/gkeyll-jenkins/agent/pr1132-verification`). This is a source
archive of the exact commit, not a Git clone. `cuda-build`, `build`, and
`valgrind-build` are separate directories. The normal Jenkins workspace was
deleted by its pipeline cleanup after failure.

Evidence is saved outside the repository under
`~/.local/share/gkeyll-jenkins/verification/`:

```text
pr1132-gpu-checks.log
pr1132-valgrind.log
pr1132-build-2-console-snapshot.log
test-jenkins-clients.py
test-regression-checker.py
```

The GPU log includes initial diagnostic link errors and a native-AVX-512
Valgrind attempt; later sections contain the successful CUDA checks. The
separate Valgrind log contains the successful compatible CPU run.
The regression fixture script requires an available LuaJIT executable
(`GKEYLL_LUAJIT`), source tree (`GKEYLL_SOURCE`), and checker path
(`GKEYLL_CHECKER`). Its successful run used LuaJIT from the now-cleaned Jenkins
workspace. Those diagnostic scripts were saved for reference, not committed
as a portable repository test suite.

## Resume commands

```sh
source "$HOME/.config/gkeyll/jenkins.env"
./ci/jenkins/gkeyll-ci.sh team scan
JENKINS_JOB=gkeyll-ci-team-workstation/PR-1132 ./ci/jenkins/gkeyll-ci.sh team recent
JENKINS_JOB=gkeyll-ci-team-workstation/PR-1132 ./ci/jenkins/gkeyll-ci.sh team follow --build 2
./ci/jenkins/gkeyll-ci.sh team run --pr 1132 --follow
```

The last command submits through the trusted `main` child. `active`, `recent`,
and build numbers refer to the selected job, not every child in the folder;
use the dashboard or set `JENKINS_JOB` to inspect automatic PR builds.
Unchanged PRs which failed during the initial credential setup are not
automatically retried just by scanning. New revisions are discovered normally;
explicitly resubmit an unchanged PR when appropriate.

## Optional plugin action rejected by automatic review

An attempted installation of **Skip Notifications Trait** was rejected by
automatic approval review because installing another plugin expands the
shared controller's trusted code surface and was not explicitly requested.
It was not installed, and no workaround was attempted. Default GitHub Branch
Source notifications therefore remain in addition to the pipeline's own
team-workstation status. Suppressing those duplicate notifications is optional
and would require approval for the plugin installation if pursued later.
