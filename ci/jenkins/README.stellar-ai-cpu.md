# Stellar AI CPU Jenkins CI

The `stellar-ai-cpu` job runs CPU unit tests, serial C regressions, and
manifest-selected MPI C regressions on Stellar AI's `cpu` partition. It uses
`machines/*.stellarai-cpu.sh`, with Open MPI 5.0.10 and no GPU requests.
Anaconda is loaded only for dependency downloads; its MPI must not shadow
Open MPI in build or simulation environments. MPI regression payloads select
`SLURM_MPI_TYPE=pmix` explicitly: Stellar AI provides the PMIx plugin but does
not configure a default MPI plugin for `srun`.

Slurm submissions start with a minimal environment and explicit CI paths,
including the dependency environment file when configured. Batch scripts skip
user login profiles, initialize Modules on the compute node, and load the CPU
toolchain there. This avoids exporting visualization-host GPU module state
whose modulefiles may not exist on CPU nodes. `SLURM_EXPORT_ENV=ALL` lets MPI
steps inherit the resulting compute-node environment.

The CI pipeline sets `ARCH_FLAGS="-march=x86-64-v3 -mtune=znver5"` in each
run's generated `config.mak`, for both candidate and baseline builds (including
copies of prebuilt configurations). Machine files and the original prebuilt
configuration are unchanged. The visualization nodes are Intel Granite Rapids,
while the CPU partition is AMD Turin (Zen 5). The shared instruction set lets
Jenkins run the installed regression driver on the visualization host and tunes
code for the AMD compute nodes. AVX2 avoids a reproduced GCC 14.3.1 AVX-512
SLP miscompilation of `binop_cross_mul_1d_3d_ser_p1`: coefficients 2 and 3
were incorrect in build #7. Disabling fast math did not fix it; compiling for
AVX2 did. Keep this setting until a compiler upgrade passes the kernel codegen
fixture and unit tests. It also avoids GCC 14.3's `-march=native` AVX10 alias
warning. The effective CI architecture flags are included in the baseline
cache fingerprint, so changing them forces a fresh baseline.

## Install and access

From a reviewed checkout on a Stellar AI login or visualization host:

```sh
./ci/jenkins/install-stellar-ai-cpu.sh
```

The installer uses Jenkins 2.580.1 LTS, verifies the WAR SHA-256, installs
Pipeline/Git/Credentials Binding plugins, and creates a private authenticated
controller. It refuses to overwrite an existing configured controller. Java 21
and tmux must be installed. The defaults are:

- Shared root: `/scratch/gpfs/CSD/$USER/gkeyll_ci_stellar_ai_cpu`
- Java: `/usr/lib/jvm/jre-21-openjdk`
- tmux session: `gkeyll_ci_stellar_ai_cpu`
- Job and node label: `gkeyll-ci-stellar-ai-cpu` and `stellar-ai-cpu`
- Slurm: partition `cpu`, account `csd`, 3 GiB of memory per allocated CPU core
- Unit tests: QoS `test`, 30-minute allocations (`00:30:00`)
- Serial and MPI regressions: QoS `short`, four-hour allocations (`04:00:00`)
- Eight compilation workers and sixteen regression worker/rank slots

Override `GKEYLL_CI_ROOT` or `JAVA_HOME` before installation if necessary.
Change the account, QoS, or worker budgets in Jenkins global environment
variables (`STELLAR_AI_CPU_SLURM_ACCOUNT`, `STELLAR_AI_CPU_SLURM_QOS`
for regressions, `STELLAR_AI_CPU_UNIT_QOS` for unit tests,
`STELLAR_AI_CPU_BUILD_JOBS`, `STELLAR_AI_CPU_REGRESSION_JOBS`).

The [StellarAI scheduling documentation](https://researchcomputing.princeton.edu/systems/stellarai)
lists CPU QoS limits of 30 minutes for `test` and four hours for `short`.
`STELLAR_AI_CPU_UNIT_TIME` and `STELLAR_AI_CPU_REGRESSION_TIME` override the
unit and regression allocation limits; keep them within the selected QoS
limit. Unit QoS is independent of `STELLAR_AI_CPU_SLURM_QOS`. CPU jobs must
request memory explicitly; the pipeline requests 3 GiB per core, below the
documented 370 GB per 96-core node.

For an existing controller, set `STELLAR_AI_CPU_SLURM_QOS=short` and
`STELLAR_AI_CPU_UNIT_QOS=test` in
**Manage Jenkins → System → Global properties → Environment variables**.
Remove existing `STELLAR_AI_CPU_UNIT_TIME` and
`STELLAR_AI_CPU_REGRESSION_TIME` overrides to use the defaults, or set them
to `00:30:00` and `04:00:00`, respectively.
The bootstrap hook preserves existing settings, so updating the checkout or
restarting Jenkins alone does not replace a configured QoS or time override.

On a virtual desktop on the **same host as Jenkins**, open
`http://127.0.0.1:8080`. Sign in as `admin`; read the password locally:

```sh
cat "$GKEYLL_CI_ROOT/jenkins_home/admin-password"
```

If `GKEYLL_CI_ROOT` is unset, substitute the expanded default path above.
From another machine, tunnel to the exact host running the controller, e.g.:

```sh
ssh -N -L 127.0.0.1:8083:127.0.0.1:8080 USER@stellarai-vis2.princeton.edu
```

Then open `http://127.0.0.1:8083`. Jenkins listens only on loopback and denies
anonymous job access. The generated API token is stored with mode 600 in
`jenkins_home/jenkins-cli.auth`; the CLI reads it without printing it.

In **Manage Jenkins → Credentials → System → Global credentials**, add a
**Username with password** credential. Use your GitHub username and token,
and ID **`gkeyll-github-stellar-ai-cpu`**. See
[GitHub reporting](README.reporting.md#reporting-credentials) for permissions.
Do not paste tokens into source files or build parameters.

## Run and restart

```sh
./ci/jenkins/gkeyll-ci.sh stellar-ai-cpu run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh stellar-ai-cpu run --candidate-ref main --baseline-ref main
./ci/jenkins/gkeyll-ci.sh stellar-ai-cpu recent
./ci/jenkins/gkeyll-ci.sh stellar-ai-cpu info --build 1
./ci/jenkins/gkeyll-ci.sh stellar-ai-cpu start
```

The controller survives terminal disconnection in tmux. It does not restart
automatically after a host reboot; run `start` on the same host to restart it.
The persistent launcher is also under `$GKEYLL_CI_ROOT/trusted/ci/jenkins/`.
Jenkins data is in scratch and remains subject to cluster storage policies.

## Workflow branch and caching

This installation always loads `ci/jenkins/jenkinsfile.stellar-ai-cpu` from
**`jenkins-github-more-reporting`** using Jenkins **Pipeline script from SCM**.
`GKEYLL_CI_TRUSTED_REF=jenkins-github-more-reporting` selects the shared reporting,
cache, comparison, Slurm, and machine scripts. Each run fetches that branch once
into a separate temporary directory and freezes its SHA for all those helpers.
The tested PR, candidate commit, and baseline commit never select workflow code.
Reports record the trusted CI SHA independently of the two tested SHAs.

The branch must contain the Stellar AI files before the first build. Push
reviewed workflow updates to that branch; subsequent builds pick them up without
editing the job. The locally installed controller queue hook and launcher are
separate: update those from a reviewed checkout and restart Jenkins while idle.
The bootstrap hook does not overwrite administrator changes on restart.

The job retains **two Jenkins builds, their artifacts, and their entire run
directories**. `GKEYLL_CI_BUILDS_TO_KEEP` overrides the count for both Jenkins
and disk retention. Startup reserves a slot for the incoming build; completion
prunes again after artifact archival and Slurm cleanup. Older inactive run
trees are deleted in full, including dependencies and compacted archives.
Jenkins/Slurm checks preserve active runs; failed checks preserve files and
fail cleanup. The shared baseline cache is separate. Leave
`GKEYLL_CI_RETAIN_RUNS` unset or `0`; `1` opts into unlimited disk history.

Pruning reads Jenkins using `$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth` and
`http://127.0.0.1:8080`. Override `CI_JENKINS_READ_AUTH_FILE` or
`CI_JENKINS_LOCAL_URL` if needed; an explicitly empty auth-file value permits
anonymous reads. Deploy both this Jenkinsfile and `cache_manager.py` on the
trusted CI branch; the next build also cleans up older runs. See
[storage](README.storage.md#retention).

The baseline lives outside Jenkins workspaces and build records at
`$GKEYLL_CI_ROOT/baseline-cache/stellar-ai-cpu/<baseline-sha>/`. The first run
creates it; later compatible runs reuse the executable, serial/MPI results,
and diagnostics. Publication requires successful, complete regression data.
An atomic manifest, checksums, and a platform lease protect against partial
publication, corruption, and concurrent replacement. The cache key includes
dependency configuration, all three machine scripts, Slurm payloads, and the
MPI reduction policy. Set `GKEYLL_CI_CACHE_REVISION` after an in-place module
or dependency upgrade. The new platform cannot reuse old Stellar CPU caches.

Cache status is recorded in `ci-baseline-cache.txt` and consolidated into
`ci-statistics.md`. A cache fixture pass verifies the lifecycle; an actual
numerical baseline is populated only after a successful regression create run.

## Verification

```sh
TEST_CACHE_PLATFORM=stellar-ai-cpu bash ci/jenkins/tests/test_baseline_cache.sh
python3 -m unittest -v ci.jenkins.tests.test_stellar_ai_mpi_policy ci.jenkins.tests.test_stellar_ai_batch_environment
python3 -m unittest -v ci.jenkins.tests.test_stellar_ai_codegen
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_stellar_resources.groovy
```

Before the first full run, verify regression timing and four-rank `srun --exact`
on a CPU allocation with the same modules, account, and QoS. The regression
runner probes for numeric timestamps, preferring GNU `date +%s.%N`, then Perl's
`Time::HiRes`, then whole-second `date +%s`. Perl being installed does not imply
that `Time::HiRes` is available on compute nodes. Build #9 ran all 320 serial
simulations successfully but misreported them as crashes because missing
`Time::HiRes` left empty timing markers and error text after their exit markers.
The clock fallback avoids requiring additional Perl modules on these nodes.
