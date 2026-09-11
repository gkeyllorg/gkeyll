# Manual Jenkins CI on Princeton Stellar Intel

This is the Level-1, CPU-only Gkeyll CI setup for Princeton's Stellar Intel
cluster. A person authenticates to Stellar with SSH/Duo, then manually starts
a Jenkins build for a GitHub pull-request number. Jenkins itself neither logs
in through Duo nor accepts GitHub webhooks.

The job runs `ci/jenkins/Jenkinsfile.stellar-intel` from a reviewed, trusted
CI branch during bring-up. It fetches the requested `refs/pull/<number>/head`
commit, builds Gkeyll and its unit-test executables on the login node, then
submits `make unit-run` as a one-core CPU Slurm job. Regression, MPI, and GPU
testing are deliberately not part of this first level.

## 0. Publish these CI files first

Create the dedicated branch `agent_tools-jenkins-stellar_intel`, review, and
push the four
implementation files to it:

```text
ci/jenkins/Jenkinsfile.stellar-intel
ci/jenkins/slurm-unit-tests.stellar-intel.sh
ci/jenkins/README.stellar-intel.md
machines/configure.stellar-intel.sh
```

Configure the initial Jenkins job to load its Pipeline script from this branch.
It is the trusted orchestration branch; do not point the job at an unreviewed
PR branch. After the pipeline has been validated and the change is merged,
switch its SCM branch specifier to `*/main`. No pipeline-code change is needed
at that point.

## 1. Obtain approval and choose locations

Before installing anything, ask Princeton Research Computing whether a
personal persistent Jenkins service is permitted on the Stellar Intel login
node and whether binding it to localhost is acceptable. Do not proceed on an
assumption: the login nodes are shared and must only be used for interactive
work such as compilation and job submission.

Choose a shared workspace root that compute nodes can access. Set it under the
project's `/scratch/gpfs` allocation, for example:

```sh
/scratch/gpfs/<project>/jenkins-workspaces
```

Do not use `/tmp`, which is node-local, and do not use an unbounded directory
in `/home`. The pipeline creates an isolated workspace below this root for
each Jenkins build. Arrange a site-appropriate scratch cleanup policy.

Record the Slurm values for your group:

```sh
qos
sshare
```

`STELLAR_SLURM_QOS` is required. `STELLAR_SLURM_ACCOUNT` is optional for PU
users but required for PPPL/CIMES users when their project policy requires
`--account`.

The initial job intentionally requests one CPU core because `make unit-run`
does not distribute the full unit suite as a Slurm MPI job. Stellar may place
requests of 47 cores or fewer in its low-priority serial queue. Do not request
an entire 96-core node merely to bypass that queue; add a genuinely parallel
test profile first.

## 2. Validate the cluster setup by hand

SSH to the Intel side and approve Duo:

```sh
ssh <NetID>@stellar.princeton.edu
```

From a disposable checkout on the selected shared filesystem, validate the
same environment the job will use:

```sh
module purge
PREFIX="$PWD/gkylsoft" ./machines/mkdeps.stellar-intel.sh
module purge
PREFIX="$PWD/gkylsoft" ./machines/configure.stellar-intel.sh
make -j3 unit
sbatch --wait --qos <your-qos> --nodes 1 --ntasks 1 --cpus-per-task 1 \
  --time 00:30:00 --chdir "$PWD" \
  ci/jenkins/slurm-unit-tests.stellar-intel.sh
```

For PPPL/CIMES, add `--account <your-account>` to the `sbatch` command. This
must complete successfully before introducing Jenkins. The configuration
targets `-march=cascadelake`, the documented Intel compute-node architecture;
it intentionally does not use login-node `-march=native` detection.

## 3. Install and expose Jenkins privately

Follow the Research Computing-approved installation method and run Jenkins as
the account that owns the Slurm allocation. Bind it to `127.0.0.1` only.
Do not expose a public port or configure a GitHub webhook for this Level-1
setup.

After SSH/Duo authentication from the laptop, open a tunnel. The first port is
on the laptop and can be any unused port; the second is Jenkins' listening
port on Stellar. If Jenkins uses its default remote port 8080 and the Mac's
local Jenkins already owns local port 8080, use 8081 locally:

```sh
ssh -N -L 8081:127.0.0.1:8080 <NetID>@stellar.princeton.edu
```

Then open `http://localhost:8081` in the laptop browser. The Mac Jenkins
remains at `http://localhost:8080`; both interfaces can be open at once. If
you configure the Stellar Jenkins service itself to listen on a non-default
remote port, replace only the final `8080` in the tunnel command.

Complete Jenkins initial setup and create an administrator account.

Install these Jenkins plugins if they are not already present: Pipeline, Git,
Credentials Binding, and the Git client plugin. The Pipeline and Git plugins
are normally part of Jenkins' suggested-plugin installation.

## 4. Create the GitHub credential

Create a GitHub fine-grained token with read-only access to the Gkeyll
repository (and pull requests/contents as required by Git). In Jenkins add it
as a Username-with-password credential, using the GitHub username and token,
and give it an ID such as `gkeyll-github-read`.

This credential only fetches source. Do not put Duo secrets, a personal SSH
key, or unrelated credentials in the Jenkins account.

## 5. Configure the Jenkins node and global environment

Use the controller's built-in node or a local agent. Give it the label
`stellar-intel` (or select another label and set
`STELLAR_INTEL_NODE_LABEL` below). The agent must run on the Stellar Intel
login side and as the same Unix account that can submit Slurm jobs.

In **Manage Jenkins → System → Global properties → Environment variables**,
set:

| Name | Required value |
| --- | --- |
| `CI_STELLAR_WORKSPACE_ROOT` | Shared project path, e.g. `/scratch/gpfs/<project>/jenkins-workspaces` |
| `STELLAR_GITHUB_CREDENTIAL_ID` | Jenkins credential ID, e.g. `gkeyll-github-read` |
| `STELLAR_SLURM_QOS` | Your valid CPU QoS, e.g. `pppl-short` |
| `STELLAR_SLURM_ACCOUNT` | Project account, if required; otherwise omit it |
| `STELLAR_SLURM_TIME` | Optional time limit; defaults to `00:30:00` |
| `STELLAR_BUILD_JOBS` | Optional login-node compile parallelism; defaults to `3` |
| `STELLAR_INTEL_NODE_LABEL` | Optional agent label; defaults to `stellar-intel` |

Do not set a broad global `PATH` to an interactive shell configuration. The
pipeline explicitly initializes the Stellar modules for every build and Slurm
job.

## 6. Create the one parameterized Pipeline job

Create **New Item → Pipeline** named `gkeyll-stellar-intel-ci`.

Configure its pipeline definition as **Pipeline script from SCM**:

- SCM: Git
- Repository: `https://github.com/gkeyllorg/gkeyll.git`
- Credentials: the read-only GitHub credential
- Branch specifier: `*/agent_tools-jenkins-stellar_intel`
- Script path: `ci/jenkins/Jenkinsfile.stellar-intel`

The CI-branch selection is intentional during bring-up. Do not point this job
at a PR branch or let the requested PR select its own Jenkinsfile. Once this
pipeline is validated and merged, change this branch specifier to `*/main`.

Save the job. It has one parameter, `PR_NUMBER`.

## 7. Run and inspect a build

After reaching Jenkins through the SSH tunnel, select **Build with
Parameters**, enter a pull-request number (for example `1104`), and start the
build. Jenkins records the exact commit in `ci-pr-commit.txt`, builds unit
tests on the login node, waits for Slurm, and archives the Slurm output as
`slurm-<jobid>.out`.

If a build queues too long or must be stopped, cancel it in Jenkins and run
`squeue --me`/`scancel <jobid>` if necessary. Automatic Slurm cancellation on
Jenkins abort is a later hardening stage, not part of Level 1.

## What this does not yet do

- automatic GitHub polling or webhooks;
- MPI regression, MOAT baseline, or GPU tests;
- Slurm job-ID polling, state classification, and automatic cancellation.

Those are future stages, after this manually-triggered CPU unit-test path is
reliable.
