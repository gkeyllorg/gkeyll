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

## 1. Choose the CI root and Slurm settings

Keep the Jenkins controller, its installation, its workspaces, and all CI
artifacts under one shared, compute-node-visible directory:

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
```

Stellar provides the username through `$USER`. The setup below creates:

```text
$GKEYLL_CI_ROOT/jenkins.war          Jenkins installation archive
$GKEYLL_CI_ROOT/jenkins_home/        JENKINS_HOME: controller configuration, plugins, jobs
$GKEYLL_CI_ROOT/jenkins_webroot/     unpacked Jenkins web application
$GKEYLL_CI_ROOT/tmp/                 controller temporary files
$GKEYLL_CI_ROOT/logs/                controller logs
$GKEYLL_CI_ROOT/workspaces/          Pipeline build workspaces and Slurm output
```

Do not use `/tmp`, `/home`, or a project directory for these CI files. The
pipeline creates an isolated workspace below `$GKEYLL_CI_ROOT/workspaces` for every
build. Since scratch storage can be purged, do not treat Jenkins build history
or credentials stored there as durable backups.

Record the Slurm values for your group:

```sh
qos
sshare
```

`STELLAR_SLURM_QOS` is required. `STELLAR_SLURM_ACCOUNT` is optional for PU
users but required for PPPL/CIMES users when their project policy requires
`--account`.

The initial job requests one task with four CPUs. `make unit-run` is launched
only once and does not distribute the full suite as a Slurm MPI job; the
four-CPU allocation is needed for its roughly 30 GB default memory allocation
(Stellar allocates 7.5 GB per core by default). In particular,
`test_dg_interpolate_3x2v_gk_ho` did not fit in a one-core allocation. Stellar
may place requests of 47 cores or fewer in its low-priority serial queue. Do
not request an entire 96-core node merely to bypass that queue; add a genuinely
parallel test profile first.

## 2. Validate the cluster setup by hand

SSH to the Intel side and approve Duo:

```sh
ssh <NetID>@stellar.princeton.edu
```

Create a disposable checkout below the CI root. This branch contains the
Stellar pipeline files; do not use an unrelated checkout that might lack the
`machines/` configuration scripts or Slurm test payload.

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
mkdir "$GKEYLL_CI_ROOT"
cd "$GKEYLL_CI_ROOT"
git clone --branch agent_tools-jenkins-stellar_intel --single-branch \
  https://github.com/gkeyllorg/gkeyll.git gkeyll
cd gkeyll
```

From this `gkeyll/` checkout, validate the same environment the job will use:
`PREFIX="$PWD/../gkylsoft"` places the dependencies at
`$GKEYLL_CI_ROOT/gkylsoft`, beside the disposable source checkout.

```sh
module purge
PREFIX="$PWD/../gkylsoft" ./machines/mkdeps.stellar-intel.sh
module purge
PREFIX="$PWD/../gkylsoft" ./machines/configure.stellar-intel.sh
make -j32 unit
sbatch --wait --qos pppl-short --nodes 1 --ntasks 1 --cpus-per-task 4 \
  --time 00:30:00 --chdir "$PWD" \
  --export=ALL,CI_WORKSPACE="$PWD" \
  ci/jenkins/slurm-unit-tests.stellar-intel.sh
```

For PPPL/CIMES, add `--account <your-account>` to the `sbatch` command. This
must complete successfully before introducing Jenkins. `CI_WORKSPACE` is the
shared checkout path that the batch payload uses after it starts on a compute
node; Jenkins supplies the same variable when it submits this job.

## 3. Install and run Jenkins privately

Stellar Intel is Red Hat Enterprise Linux 8.10. Because this is an
unprivileged personal setup, use the Jenkins WAR directly rather than a
system RPM/service. This procedure requires no root access: do not use
`sudo`, `dnf`, `rpm`, or `systemctl`. The login node's default Java is Java
17; current Jenkins requires Java 21 or later. Stellar's Java 21 location is
`/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64`. Do not use the
default `java` command unless `java -version` confirms it is Java 21 or newer.

After logging in through Duo, set up the scratch-only controller directory:

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64
export PATH="$JAVA_HOME/bin:$PATH"
export JENKINS_HOME="$GKEYLL_CI_ROOT/jenkins_home"
export JENKINS_WEBROOT="$GKEYLL_CI_ROOT/jenkins_webroot"
export TMPDIR="$GKEYLL_CI_ROOT/tmp"

mkdir -p "$JENKINS_HOME" "$JENKINS_WEBROOT" "$TMPDIR" \
  "$GKEYLL_CI_ROOT/logs" "$GKEYLL_CI_ROOT/workspaces"
java -version
```

The final command must report Java 21 or newer. Download the Jenkins LTS WAR
into the same directory; if Stellar cannot reach the download URL, download it
on the laptop and copy it into `$GKEYLL_CI_ROOT` through the authenticated SSH
connection instead:

```sh
cd "$GKEYLL_CI_ROOT"
curl -fL -o jenkins.war https://get.jenkins.io/war-stable/latest/jenkins.war
```

Start Jenkins in a named tmux session. This is the initial operational mode:
it keeps the controller available after SSH disconnects without installing a
system service. Stop it after a test session if it is not needed.

```sh
tmux new -s gkeyll-jenkins

export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64
export PATH="$JAVA_HOME/bin:$PATH"
export JENKINS_HOME="$GKEYLL_CI_ROOT/jenkins_home"
export JENKINS_WEBROOT="$GKEYLL_CI_ROOT/jenkins_webroot"
export TMPDIR="$GKEYLL_CI_ROOT/tmp"

java -Djava.io.tmpdir="$TMPDIR" -jar "$GKEYLL_CI_ROOT/jenkins.war" \
  --webroot="$JENKINS_WEBROOT" \
  --httpListenAddress=127.0.0.1 \
  --httpPort=8080 \
  2>&1 | tee -a "$GKEYLL_CI_ROOT/logs/jenkins.log"
```

In another Stellar shell, use these commands to manage that controller:

```sh
tmux attach -t gkeyll-jenkins   # return to the Jenkins console
tmux ls                         # verify the session exists
tmux kill-session -t gkeyll-jenkins  # stop Jenkins after testing
```

The controller must run as the same Unix account that owns the Slurm
allocation. It binds only to `127.0.0.1`; do not expose a public port or
configure a GitHub webhook for this Level-1 setup.

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

On its first start, get the unlock password with:

```sh
cat "$JENKINS_HOME/secrets/initialAdminPassword"
```

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
| `CI_STELLAR_WORKSPACE_ROOT` | `/scratch/gpfs/$USER/gkeyll_ci/workspaces` |
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

Create **New Item → Pipeline** named `gkeyll-ci-stellar-intel`.

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
