# Jenkins CI on Slurm clusters

Use `gkeyll-ci.sh slurm COMMAND` for any CPU or GPU cluster. Select the cluster
with `GKEYLL_CI_CLUSTER`; for example, `perlmutter_gpu` selects
`slurm-clusters/slurm_ci_profile_perlmutter_gpu.sh`.

A profile contains Bash settings for resources, paths and build commands. Each
cluster keeps its own submitted Slurm scripts for module loading, environment
preparation and test launching. The shared Jenkins Pipeline handles builds,
submission, comparison with a baseline, cancellation and reports.

- First setup: [choose a profile](#choose-a-profile), [install the controller](#install-the-controller), then [configure the trusted job](#configure-the-trusted-job).
- Daily use: [run and inspect builds](#operating-the-cluster-ci).
- Configuration: [change settings](#changing-settings), including walltime, or consult the [complete profile reference](slurm-clusters/README.md#shell-profile-reference).
- New sites: [add a cluster](#add-a-cluster) and follow the [site script contract](#site-script-contract).
- Troubleshooting: [retained data and cleanup](#retained-data-and-cleanup).

## Choose a profile

Run on a permitted login-side host with Bash, Python 3, Git, curl, tmux, Linux
`ss`, Slurm client tools and a Java installation compatible with your Jenkins
release (the shipped setups use Java 21). Python must be available before site
modules load. The Jenkins agent must be able to submit jobs under your account.
Use shared storage with filesystem locks: workspaces, installations and trusted
CI files must have the same absolute paths on the agent and compute nodes.

These instructions assume a fresh installation with a new empty CI root. Stop
old builds and controllers before reusing their ports. Old cluster clients,
environment aliases and caches are not migrated. Back up Jenkins credentials
and configuration outside purgeable scratch storage.

| Profile | Unit / serial resources | Parallel resources | Compilation |
| --- | --- | --- | --- |
| [Stellar CPU](slurm-clusters/slurm_ci_profile_stellar_cpu.sh) | 1 task; 4 / 8 CPUs; no GPUs | 4 tasks, 1 CPU each | Agent; Stellar Intel toolchain |
| [Perlmutter GPU](slurm-clusters/slurm_ci_profile_perlmutter_gpu.sh) | 1 task, 32 CPUs, 1 GPU | 4 tasks, 8 CPUs and 1 GPU per task | Agent; Perlmutter CUDA/NCCL toolchain |
| [StellarAI GPU](slurm-clusters/slurm_ci_profile_stellarai_gpu.sh) | 1 task, 4 CPUs, 1 GPU | 4 tasks, 1 CPU each; 4 GPUs on one node | Agent; StellarAI CUDA/NCCL toolchain |

All profiles use one node. Their default walltimes are 30 minutes for unit
jobs and 4 hours for each serial or parallel regression job. For another site,
[create a profile](#add-a-cluster) before continuing.

### Stellar CPU example

In your authenticated Stellar login session, set:

```sh
export GKEYLL_CI_CLUSTER=stellar_cpu
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.13.0.9-1.1.el8.x86_64
export SLURM_CI_QOS='your-qos'
# Set this only if your allocation requires an account:
export SLURM_CI_ACCOUNT='your-account'
```

The root and Java path above are also the profile defaults. Override Java if
the installed path changes. Choose a QoS/account permitted for your allocation.

### Perlmutter GPU example

In your authenticated Perlmutter login session, set:

```sh
export GKEYLL_CI_CLUSTER=perlmutter_gpu
export GKEYLL_CI_ROOT=/pscratch/sd/${USER:0:1}/$USER/gkeyll_ci
export JAVA_HOME=/absolute/path/to/your/java21
export SLURM_CI_ACCOUNT='your-project'
```

Perlmutter has no profile default for the root or Java path; supply both.
If you use an OpenJDK installation managed by the site's Spack module, load it
before setting `JAVA_HOME`, for example:

```sh
module load spack
spack load openjdk@21
export JAVA_HOME="$(spack location --install-dir openjdk@21)"
```

The profile requests `shared` QoS and the `gpu` constraint. Account, resource
and walltime requests must fit your site's scheduler policy.

### StellarAI GPU example

On a permitted StellarAI host with the GPU build toolchain, set:

```sh
export GKEYLL_CI_CLUSTER=stellarai_gpu
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci_stellarai_gpu
export JAVA_HOME=/absolute/path/to/your/java21
export SLURM_CI_QOS='your-gpu-qos'
# Set this only if your allocation requires an account:
export SLURM_CI_ACCOUNT='your-account'
./ci/jenkins/gkeyll-ci.sh slurm validate
```

Supply the shared root and Java path; this profile has no defaults for them.
Choose a GPU QoS that permits the requested walltimes. The profile selects the
`gpu` partition and uses the three `machines/*stellarai-gpu.sh` build/environment
files. Compilation runs on the Jenkins agent; use a host permitted for builds
that can also execute the resulting Emerald Rapids CPU code. The configure
script targets CUDA architecture 100. Check current host access and GPU QoS
names in the [StellarAI site guide](https://researchcomputing.princeton.edu/systems/stellarai).

The submitted scripts unload `anaconda3/2026.7` after loading the environment
and set `UCX_TLS=sm,cuda_copy,self`, matching `machines/jobScript.stellarai-gpu.sh`.
Anaconda remains available during dependency installation for ADAS downloads.
Parallel regressions use `mpirun`, four ranks and four GPUs on one node; keep
this allocation on one node because these UCX transports do not span nodes.
Serial regressions run one at a time on one GPU (`SLURM_CI_SERIAL_JOBS=1`).

Follow the common controller and trusted-job setup below, using job name
`gkeyll-ci-stellarai_gpu` and node label `stellarai_gpu`. The trusted CI revision
must contain this profile and its payloads, and both candidate and baseline
must contain the StellarAI machine files. A local profile alone does not make
it available to an existing Jenkins job pinned to another revision.

## Install the controller

From the repository root of a reviewed checkout, after setting the environment
for your profile:

```sh
"$JAVA_HOME/bin/java" -version
./ci/jenkins/gkeyll-ci.sh slurm validate
mkdir -p "$GKEYLL_CI_ROOT"
curl -fL -o "$GKEYLL_CI_ROOT/jenkins.war" https://get.jenkins.io/war-stable/latest/jenkins.war
./ci/jenkins/gkeyll-ci.sh slurm start
hostname -f
```

`validate` checks local profile settings and script syntax and prints resolved
settings plus `sbatch` arguments. It does not contact Jenkins or Slurm or run
site scripts. `start` runs Jenkins in a detached tmux session, listening on
`127.0.0.1:8080` by default. The controller survives logout.

From your laptop, tunnel to the **exact host** printed by `hostname -f`:

```sh
ssh -N -o ExitOnForwardFailure=yes \
  -L 127.0.0.1:8083:127.0.0.1:8080 your-user@controller.example.edu
```

On Perlmutter, use its public address as a jump host if needed; replace the
final hostname with the recorded login host:

```sh
ssh -N -o ExitOnForwardFailure=yes \
  -J your-user@perlmutter.nersc.gov \
  -L 127.0.0.1:8083:127.0.0.1:8080 your-user@login20.chn.perlmutter.nersc.gov
```

Open `http://127.0.0.1:8083`, unlock Jenkins with the password in
`$GKEYLL_CI_ROOT/jenkins_home/secrets/initialAdminPassword`, and install Pipeline,
Git and Credentials Binding with their dependencies.

If port 8080 is occupied on the cluster host, set `JENKINS_PORT` before
starting Jenkins and change the tunnel's destination port to match. Set
`CI_JENKINS_LOCAL_URL` in Jenkins to that host-side URL. The laptop's forwarded
port is independent. Restore your cluster/root/Java/port settings in new shells.
Separate controllers on one host need distinct roots, homes, ports and tmux
sessions (`JENKINS_SESSION`).

Create a Jenkins API token for your Jenkins user, then store it on the cluster:

```sh
umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' > "$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth"
chmod 600 "$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth"
```

`JENKINS_CLI_AUTH_FILE` can override this client credential path. Separately,
create a Jenkins **Username with password** credential containing your GitHub
username and PAT; record its credential ID for the next step. Current CI uses
commit statuses and report comments, so the PAT needs the permissions described
in [Results in GitHub](README.md#results-in-github). Tokens never go in profiles.

## Configure the trusted job

For a controller dedicated to one cluster profile, open **Manage Jenkins →
System → Global properties → Environment variables** and add these settings:

| Variable | Value / purpose |
| --- | --- |
| `GKEYLL_CI_CLUSTER` | Profile ID, such as `perlmutter_gpu`. |
| `GKEYLL_CI_ROOT` | The expanded absolute shared root used by the CLI. |
| `SLURM_CI_GITHUB_CREDENTIAL_ID` | The GitHub credential ID created above. |
| `SLURM_CI_ACCOUNT`, `SLURM_CI_QOS` | Account/QoS overrides required by your profile. |
| `SLURM_CI_NODE_LABEL` | Optional agent label; defaults to the profile ID. |
| `GKEYLL_CI_TRUSTED_REF` | Reviewed branch or full commit SHA for CI code; defaults to `main`. |
| `CI_JENKINS_LOCAL_URL` | Agent-reachable Jenkins URL; defaults to `http://127.0.0.1:8080`, not the laptop tunnel URL. |
| `CI_JENKINS_READ_AUTH_FILE` | Agent-readable API credential file for cache activity queries; defaults to `$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth`. Empty only if Jenkins permits anonymous build-state reads. |

Save the settings. The profile ID and node label must be available before agent
allocation; other overrides may also come from the agent environment. These
installation variables are **not profile-file options**, except the account
and QoS settings also listed in the reference below. Shell exports used for the
CLI are not automatically forwarded to Jenkins builds.

Give the login-side agent the selected label. Create a Pipeline job named
`gkeyll-ci-PROFILE` (or your profile's `SLURM_CI_JOB_NAME`) and configure:

- Definition: **Pipeline script from SCM**; SCM: **Git**.
- Repository: `https://github.com/gkeyllorg/gkeyll.git`; no SCM credential.
- Branch: reviewed `main`, or the reviewed CI branch you are staging.
- Script path: `ci/jenkins/jenkinsfile.slurm_cluster`.

When staging CI changes, set `GKEYLL_CI_TRUSTED_REF` to the reviewed revision
containing them too. Never select the candidate PR as the trusted CI source.
The Jenkinsfile pins the profile, helpers and site scripts, then loads its own
copy from that snapshot to run the workflow. Do not set the internal
`CI_SLURM_BOOTSTRAPPED` flag yourself.

Run **Build Now** once to register parameters; its empty-selector failure is
expected and submits no Slurm jobs. Then use **Build with Parameters** or the CLI.

For several profiles on one controller, use distinct job/cache names and set
the profile and its overrides in each job's trusted inline bootstrap instead
of global settings. For example:

```groovy
node('example_cpu') {
    dir('_slurm_entry') {
        git branch: 'main', url: 'https://github.com/gkeyllorg/gkeyll.git'
        withEnv(['GKEYLL_CI_CLUSTER=example_cpu', 'CI_SLURM_ENTRY_NODE=true']) {
            load 'ci/jenkins/jenkinsfile.slurm_cluster'
        }
    }
}
```

## Operating the cluster CI

Run from the repository root on the cluster with your installation environment
restored. `GKEYLL_CI_CLUSTER` is required; no hostname detection is performed.

```sh
./ci/jenkins/gkeyll-ci.sh slurm run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh slurm run --candidate-ref feature --baseline-ref main
./ci/jenkins/gkeyll-ci.sh slurm active
./ci/jenkins/gkeyll-ci.sh slurm recent --limit 5
./ci/jenkins/gkeyll-ci.sh slurm status --queue 42
./ci/jenkins/gkeyll-ci.sh slurm follow --queue 42
./ci/jenkins/gkeyll-ci.sh slurm info --build 42
./ci/jenkins/gkeyll-ci.sh slurm artifact --build 42 --fetch
./ci/jenkins/gkeyll-ci.sh slurm abort --build 42
```

`run` returns a queue ID; `--follow` streams the console. Ctrl-C stops following
without aborting the build. Cancel with `abort --queue ID` or `abort --build
NUMBER`. Use `gkeyll-ci.sh slurm COMMAND --help` for command-specific flags.

In **Build with Parameters**, set `CANDIDATE_PR`, or both `CANDIDATE_REF` and
`BASELINE_REF`. Candidates must contain their baseline; use
`--allow-behind-candidate` / `ALLOW_BEHIND_CANDIDATE` only for intentional
historical comparisons. Accepted regression differences still use the tracked
file described in [Numerical regression differences](README.md#numerical-regression-differences).

## Changing settings

Edit `ci/jenkins/slurm-clusters/slurm_ci_profile_PROFILE.sh` to change defaults
for that profile; see the [available options and defaults](slurm-clusters/README.md#shell-profile-reference).
Commit the change to the trusted CI revision Jenkins loads
(`main` unless `GKEYLL_CI_TRUSTED_REF` selects another revision). An uncommitted
local edit, or an edit only in the candidate PR, does not change Jenkins's profile.

For an installation-specific override, add the scalar variable with the same
name in Jenkins's environment, using the location described in
[Configure the trusted job](#configure-the-trusted-job). Global overrides apply
to all profiles on that controller. Start a new build to use the new settings;
running builds retain their resolved configuration.

### Example: increase Perlmutter walltime

In `slurm-clusters/slurm_ci_profile_perlmutter_gpu.sh`, add:

```sh
SLURM_CI_UNIT_TIME=01:00:00
SLURM_CI_SERIAL_TIME=08:00:00
SLURM_CI_PARALLEL_TIME=08:00:00
```

This requests one hour for each unit job and eight hours for each regression
job, including separate baseline-creation and candidate-comparison submissions.
It is not an overall Pipeline timeout. Use `SLURM_CI_BUILD_TIME` for compilation
jobs only when `SLURM_CI_BUILD_LOCATION=slurm`. `SLURM_CI_SERIAL_TEST_TIMEOUT`
is a separate timeout for an individual serial regression test, in seconds.

Validate the local edit with your account/QoS environment set:

```sh
GKEYLL_CI_CLUSTER=perlmutter_gpu ./ci/jenkins/gkeyll-ci.sh slurm validate
```

Check `UNIT_TIME`, `SERIAL_TIME`, `PARALLEL_TIME` and the `--time=...` allocation
arguments. After deploying the change, inspect the new Jenkins build's
`ci-slurm-profile.json` artifact to confirm its actual settings. Local validation
does not read Jenkins's environment overrides.

## Add a cluster

1. Choose a profile ID using letters, digits, underscores or hyphens, starting
   with a letter or digit. `personal`, `team` and `slurm` are reserved.
2. Copy a shipped profile to `ci/jenkins/slurm-clusters/slurm_ci_profile_PROFILE.sh`.
   Update its display name, root, Java and machine scripts. Job/cache/status
   names default to the new ID; remove or replace copied identity overrides.
3. Create the site's unit, serial and parallel scripts in the same folder,
   starting from the shipped examples. Set `SLURM_CI_UNIT_SCRIPT`,
   `SLURM_CI_SERIAL_SCRIPT` and `SLURM_CI_PARALLEL_SCRIPT` to their paths.
   Customize preparation and launching, and set allocations in the profile.
   Remove inherited GPU requests/constraints for a CPU profile.
4. If compilation is not permitted on the agent, set
   `SLURM_CI_BUILD_LOCATION=slurm`, supply `SLURM_CI_BUILD_SCRIPT` and configure
   its resources using the contract below.
5. Export `GKEYLL_CI_CLUSTER=PROFILE` and required account/QoS values. Run
   `./ci/jenkins/gkeyll-ci.sh slurm validate`, then follow the installation
   steps. The trusted CI revision must contain the new profile and site scripts;
   machine scripts must exist in both tested source revisions.
6. Complete [verification](#verification-before-enabling-a-profile).

No shared-code edits or cluster-name switch are needed for a new profile.
Use separate profiles and cache namespaces for different CPU/GPU toolchains.

## Site script contract

### Scripts submitted to Slurm

`slurm-submit.sh` runs on the Jenkins agent and passes the profile's resource
arguments and selected site script to `sbatch`. The scripts below execute on
compute nodes; each regression script runs for baseline creation on a cache
miss and again for candidate comparison.

| Profile option | Stellar CPU script | Perlmutter GPU script | StellarAI GPU script |
| --- | --- | --- | --- |
| `SLURM_CI_UNIT_SCRIPT` | [slurm-unit-tests.stellar_cpu.sh](slurm-clusters/slurm-unit-tests.stellar_cpu.sh) | [slurm-unit-tests.perlmutter_gpu.sh](slurm-clusters/slurm-unit-tests.perlmutter_gpu.sh) | [slurm-unit-tests.stellarai_gpu.sh](slurm-clusters/slurm-unit-tests.stellarai_gpu.sh) |
| `SLURM_CI_SERIAL_SCRIPT` | [slurm-serial-regression-tests.stellar_cpu.sh](slurm-clusters/slurm-serial-regression-tests.stellar_cpu.sh) | [slurm-serial-regression-tests.perlmutter_gpu.sh](slurm-clusters/slurm-serial-regression-tests.perlmutter_gpu.sh) | [slurm-serial-regression-tests.stellarai_gpu.sh](slurm-clusters/slurm-serial-regression-tests.stellarai_gpu.sh) |
| `SLURM_CI_PARALLEL_SCRIPT` | [slurm-parallel-regression-tests.stellar_cpu.sh](slurm-clusters/slurm-parallel-regression-tests.stellar_cpu.sh) | [slurm-parallel-regression-tests.perlmutter_gpu.sh](slurm-clusters/slurm-parallel-regression-tests.perlmutter_gpu.sh) | [slurm-parallel-regression-tests.stellarai_gpu.sh](slurm-clusters/slurm-parallel-regression-tests.stellarai_gpu.sh) |

Keep site module loads, unloads, runtime exports and launch commands in these
scripts, before test execution. Do not duplicate profile resources in `#SBATCH`
directives. Preserve `set -euo pipefail` and propagate preparation/test failures.
Parallel regression drivers launch MPI tests using manifest cuts: provide enough
resources, retain the `parallel-c-4` results directory, and do not wrap the
driver in an outer MPI launch.

### Inputs and outputs

The Pipeline exports `CI_WORKSPACE`, candidate/baseline source directories
(`CI_CANDIDATE_DIR`, `CI_BASELINE_DIR`) and installation prefixes
(`CI_CANDIDATE_PREFIX`, `CI_BASELINE_PREFIX`) as appropriate for the job type.
Regression scripts receive `CI_REGRESSION_MODE` (`baseline-create` or
`candidate-check`) and `CI_TRUSTED_CHECKER`; serial scripts also receive
`CI_SERIAL_REGRESSION_JOBS` and `CI_SERIAL_REGRESSION_TEST_TIMEOUT`.
`CI_SLURM_CPUS_PER_TASK` and `CI_SLURM_GPUS_PER_TASK` describe the allocation.
All resolved scalar `SLURM_CI_*` settings are inherited.

Submitted scripts and the checker come from the pinned trusted snapshot.
Preserve the shipped scripts' summary/timing filenames, including
`ci-serial-regression-summary.txt`, `ci-parallel-regression-summary.txt`, and
candidate/baseline `*-serial-c-*` / `*-parallel-c-*` timing files. The checker
receives the candidate/baseline acknowledgment files. Use the shipped scripts as
executable examples of the full contract.

### Compilation inside Slurm

With `SLURM_CI_BUILD_LOCATION=slurm`, each build operation submits your
`SLURM_CI_BUILD_SCRIPT`. That script prepares its environment and calls the
shared helper:

```sh
#!/bin/bash -l
set -euo pipefail
# Site-specific preparation here.
bash "$CI_BUILD_HELPER"
```

The helper receives `CI_BUILD_ACTION` (`deps`, `configure`, `unit`, `install`,
`serial-compile`, `parallel-compile`), `CI_BUILD_SOURCE`, `CI_BUILD_PREFIX`,
`CI_BUILD_LABEL` and `CI_WORKSPACE`. It sources `SLURM_CI_ENV_SCRIPT` and uses
the machine scripts in the corresponding candidate/baseline tree. The environment
script must initialize any required module functions. Each operation completes
before its dependent stage.

## Retained data and cleanup

Baseline caches are checked against a fingerprint of the build/test settings
and trusted CI revision. A configuration change can cause a rebuild; active
baseline users finish before their data is replaced. Controller paths/ports and
display identity are excluded from the fingerprint. Each build archives its
resolved `ci-slurm-profile.json` and keeps trusted files under
`$GKEYLL_CI_ROOT/cache-control/BUILD_TAG/trusted` until activity-checked cleanup.

Cache cleanup queries Jenkins and Slurm before deleting old data. Unknown
activity preserves files. Workspaces remain available for later cleanup; do not
manually remove a tree used by a running job. Jenkins artifacts and the GitHub
report contain the final outcome. Retained candidate diagnostics are captured
before cache finalization; the final report follows cleanup so it includes failures.

### Cancellation and process cleanup

The helper records a unique ownership token before submission. It and the
Pipeline finalizer use one cleanup routine to verify jobs by token, user and
workspace, even if submission was interrupted before its numeric ID was saved.
Scheduler checks do not rely on old or incomplete status files. Failed or
interrupted submissions attempt cancellation; successful runs verify that
allocations have ended. Ownership/status records are archived.

[Jenkins handles agent-side build abortion](https://www.jenkins.io/doc/book/using/aborting-a-build/);
[`scancel` terminates compute allocations](https://slurm.schedmd.com/scancel.html).
Site scripts must keep work inside their allocation and wait for background
tasks. The persistent controller and tmux session remain running between builds.

Forced Jenkins termination, an unavailable scheduler, or acceptance after
recovery polling ends can prevent immediate cleanup. Ownership records and
workspaces are retained for retry. After restoring access, run:

```sh
CI_WORKSPACE=/absolute/path/to/retained/workspace \
  bash /absolute/path/to/reviewed/ci/jenkins/slurm-submit.sh cleanup
```

For diagnosis, inspect Jenkins artifacts, `squeue --me`, and the controller log
at `$GKEYLL_CI_ROOT/logs/jenkins.log`; attach with
`tmux attach -t "${JENKINS_SESSION:-gkeyll_ci}"`.

## Verification before enabling a profile

Run focused local checks from the repository root:

```sh
python3 -m unittest discover -s ci/jenkins -p 'test_*.py'
bash ci/jenkins/test_baseline_cache.sh
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/test_slurm_pipeline.groovy
```

The Groovy fixtures exercise bootstrap, cache reuse/builds, failures and cleanup
with mocked Pipeline steps, and syntax-check reached shell steps. Use Jenkins's
Groovy JAR. These checks do not validate site modules, storage visibility,
scheduler policy or actual CPU/GPU execution.

Before enabling a new site, try its machine/build scripts in a disposable shared
checkout. In live Jenkins, run a same-SHA candidate/baseline comparison through
unit, serial and parallel jobs, then repeat to confirm a baseline-cache hit.
Abort short queued and running jobs and an agent-side compilation. Verify that
allocations and test/compiler processes disappear after completion or abort.
