# Slurm cluster profiles and job scripts

This directory contains the settings and submitted job scripts for each cluster.
For installation, CLI use and adding a new cluster, see the
[Slurm CI setup guide](../README.slurm_cluster.md). The
[profile reference](#shell-profile-reference) below lists every editable option.

## Four file types

`PROFILE` identifies a CI environment, such as `stellar_cpu` or `perlmutter_gpu`.
A cluster can have separate CPU and GPU profiles.

| File pattern | Purpose |
| --- | --- |
| `slurm_ci_profile_PROFILE.sh` | Bash assignments selecting machine scripts, submitted scripts, resource requests and controller defaults. Read by the shared configuration loader; not submitted to Slurm. |
| `slurm-unit-tests.PROFILE.sh` | Submitted job script that prepares the site environment and runs candidate unit tests. |
| `slurm-serial-regression-tests.PROFILE.sh` | Submitted job script that creates serial regression baselines or checks candidate results against them. Serial tests may run concurrently within the allocation. |
| `slurm-parallel-regression-tests.PROFILE.sh` | Submitted job script that creates or checks parallel regression results; the regression driver launches each MPI test. |

The shipped profiles are [Stellar CPU](slurm_ci_profile_stellar_cpu.sh),
[Perlmutter GPU](slurm_ci_profile_perlmutter_gpu.sh), and
[StellarAI GPU](slurm_ci_profile_stellarai_gpu.sh). The profile supplies
allocations to the shared submission helper. The three test scripts own site
module loads, unloads, runtime exports and launch commands; edit those scripts
to add preparation before test execution. Keep resource requests in the profile
rather than duplicating them in `#SBATCH` directives.

All shipped profiles compile on the Jenkins agent. A cluster that compiles
inside Slurm also supplies its own `SLURM_CI_BUILD_SCRIPT`; see the
[build script contract](../README.slurm_cluster.md#compilation-inside-slurm).
Shared Jenkins code and submission/build helpers remain in the parent directory.

## Shell profile reference

This is the complete editable vocabulary accepted by [`slurm_config.py`](../slurm_config.py).
Use Bash assignments and arrays; quoting and `${USER}` expansion work. Keep
module loading and other commands in site scripts. All paths ending in
`_SCRIPT` are repository-relative, must exist, and must pass `bash -n`.
Root/Java paths must be absolute; path traversal, control characters and shell
metacharacters rejected by the loader are not allowed. Unknown `SLURM_CI_*`
variables assigned in a profile are errors.

Values below are **shared defaults**; the shipped profiles override resource
counts and machine paths. Only set the values your profile needs to change.

### Identity, controller and build options

| Profile option | Default | Meaning |
| --- | --- | --- |
| `SLURM_CI_CONFIG_VERSION` | `1` | Supported profile format version. |
| `SLURM_CI_DISPLAY_NAME` | Profile ID | Descriptive name recorded in configuration. |
| `SLURM_CI_JOB_NAME` | `gkeyll-ci-PROFILE` | Jenkins job the CLI addresses; create the job with this name. |
| `SLURM_CI_CACHE_NAMESPACE` | Profile ID | Cache directory name; distinct toolchains need distinct namespaces. |
| `SLURM_CI_STATUS_CONTEXT` | `continuous-integration/jenkins/PROFILE` | GitHub commit-status context. |
| `SLURM_CI_ROOT` | Empty | Shared CI root; required to run CI. |
| `SLURM_CI_JAVA_HOME` | Empty | Java installation directory; required to start the controller. |
| `SLURM_CI_PORT` | `8080` | Controller HTTP port, 1–65535. |
| `SLURM_CI_SESSION` | `gkeyll_ci` | Controller tmux session name. |
| `SLURM_CI_ENV_SCRIPT` | Empty; required | Machine environment script sourced for compilation. |
| `SLURM_CI_MKDEPS_SCRIPT` | Empty; required | Machine dependency-build script. |
| `SLURM_CI_CONFIGURE_SCRIPT` | Empty; required | Machine configure script; selects CPU/CUDA toolchain. |
| `SLURM_CI_BUILD_LOCATION` | `agent` | Compile on the Jenkins agent, or use `slurm` allocations. |
| `SLURM_CI_BUILD_JOBS` | `3` | Positive concurrency for unit builds, installation, and C regression compilation. |
| `SLURM_CI_SERIAL_JOBS` | `4` | Positive number of serial regression tests run concurrently. |
| `SLURM_CI_SERIAL_TEST_TIMEOUT` | `900` | Positive timeout in seconds per serial regression test. |
| `SLURM_CI_MPIEXEC` | `srun` | Parallel regression launcher. |
| `SLURM_CI_MPI_ARGS` | Empty | Additional launcher arguments used by parallel regression configuration. |

Job, cache and session names use letters, digits, underscores or hyphens,
starting with a letter or digit. Status contexts also allow dots and slashes.

### Common scheduler options

| Profile option | Default | Meaning |
| --- | --- | --- |
| `SLURM_CI_ACCOUNT` | Empty | Account for all submitted jobs. |
| `SLURM_CI_QOS` | Empty | QoS for all submitted jobs. |
| `SLURM_CI_PARTITION` | Empty | Partition for all submitted jobs. |
| `SLURM_CI_CONSTRAINT` | Empty | Node constraint for all submitted jobs. |
| `SLURM_CI_ACCOUNT_REQUIRED` | `false` | Set `true` to require an account for every submitted job type. |
| `SLURM_CI_QOS_REQUIRED` | `false` | Set `true` to require a QoS for every submitted job type. |
| `SLURM_CI_SBATCH_EXTRA_ARGS` | `()` | Bash array of extra resource options applied to all jobs. |

Empty optional scheduler settings omit the corresponding `sbatch` flag.
Stellar requires a QoS; Perlmutter requires an account and defaults to `shared`
QoS. There are no common memory, GPU-count or walltime variables; set those per
job type below.

### Options for each job type

Replace `<LANE>` with **`UNIT`, `SERIAL`, `PARALLEL` or `BUILD`**. Each row
therefore defines four options: for example, `SLURM_CI_UNIT_NODES`,
`SLURM_CI_SERIAL_NODES`, `SLURM_CI_PARALLEL_NODES`, `SLURM_CI_BUILD_NODES`.
`BUILD` resources apply only when `SLURM_CI_BUILD_LOCATION=slurm`.

| Profile option | Default | Meaning |
| --- | --- | --- |
| `SLURM_CI_<LANE>_SCRIPT` | Empty | Submitted site script; required for unit/serial/parallel, and for Slurm compilation. |
| `SLURM_CI_<LANE>_NODES` | `1` | Positive number of nodes (`--nodes`). |
| `SLURM_CI_<LANE>_TASKS` | `4` for parallel; `1` otherwise | Positive number of tasks (`--ntasks`). |
| `SLURM_CI_<LANE>_CPUS_PER_TASK` | `1` | Positive CPUs per task (`--cpus-per-task`). |
| `SLURM_CI_<LANE>_GPUS_PER_TASK` | `0` | GPUs per task; zero omits `--gpus-per-task`. |
| `SLURM_CI_UNIT_TIME` | `00:30:00` | Unit-job walltime (`--time`). |
| `SLURM_CI_SERIAL_TIME` | `04:00:00` | Serial-regression job walltime. |
| `SLURM_CI_PARALLEL_TIME` | `04:00:00` | Parallel-regression job walltime. |
| `SLURM_CI_BUILD_TIME` | `04:00:00` | Walltime for each Slurm compilation operation. |
| `SLURM_CI_<LANE>_MEM` | Empty | Memory per node (`--mem`), e.g. `32G`. |
| `SLURM_CI_<LANE>_GRES` | Empty | Generic resources (`--gres`), e.g. `gpu:2`; cannot combine with nonzero GPUs per task. |
| `SLURM_CI_<LANE>_ACCOUNT` | Empty | Nonempty value overrides the common account. |
| `SLURM_CI_<LANE>_QOS` | Empty | Nonempty value overrides the common QoS. |
| `SLURM_CI_<LANE>_PARTITION` | Empty | Nonempty value overrides the common partition. |
| `SLURM_CI_<LANE>_CONSTRAINT` | Empty | Nonempty value overrides the common node constraint. |
| `SLURM_CI_<LANE>_EXTRA_ARGS` | `()` | Bash array of extra options appended to common extras. |

Walltime syntax is `HH:MM:SS` or `D-HH:MM:SS`. GPU counts must be nonnegative
integers. Requesting GPUs does not enable CUDA: machine scripts and site
launch commands must agree with the resources.

Extra arrays support only `--reservation`, `--licenses`, `--hint`, `--cpu-bind`,
`--gpu-bind`, `--gpus-per-node`, `--ntasks-per-node`, `--exclusive`,
`--mem-per-cpu` and `--gres-flags`. Use `--option=value`, except `--exclusive`:

```sh
SLURM_CI_SBATCH_EXTRA_ARGS=(--reservation=my-reservation)
SLURM_CI_PARALLEL_EXTRA_ARGS=(--gpu-bind=closest)
```

Duplicate options across common/per-job extras are rejected. Do not combine
`--gpus-per-node` with GRES or GPUs per task, or `--mem-per-cpu` with per-node
memory. Working directories, logging, exports and payload selection are managed
by the shared helper and cannot be overridden through these arrays.

### Override precedence

The order is shared defaults → profile assignments → same-name environment
overrides. Arrays are edited in the profile; they are not environment overrides.
These installation variables take precedence over both the profile and its
same-name environment variable:

| Installation variable | Profile setting overridden |
| --- | --- |
| `GKEYLL_CI_ROOT` | `SLURM_CI_ROOT` |
| `JAVA_HOME` | `SLURM_CI_JAVA_HOME` |
| `JENKINS_PORT` | `SLURM_CI_PORT` |
| `JENKINS_SESSION` | `SLURM_CI_SESSION` |
| `JENKINS_JOB` | `SLURM_CI_JOB_NAME` |

Keep credentials and Jenkins-only variables in the installation environment,
not the profile. The JSON snapshot is generated automatically; you do not edit it.

