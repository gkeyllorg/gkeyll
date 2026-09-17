# Gkeyll Jenkins CI on Princeton Stellar CPU

This guide sets up Gkeyll's private, CPU-only Jenkins CI on Princeton Stellar.
After SSH/Duo authentication, a developer manually starts a build for a GitHub
pull request or an explicit candidate/baseline comparison. Jenkins runs on the
login side and submits unit tests and C regressions to Slurm CPU nodes.

The job uses a reviewed Pipeline and fixed reviewed baseline. It does not use
GitHub webhooks or run Lua, MOAT-only, or GPU regressions. Alongside the
serial C suite it runs manifest-selected C MPI regressions on four ranks. Jenkins is
available only through the authenticated Stellar login session or an SSH tunnel.

# Installation

## Install Jenkins

SSH to Stellar's CPU side and approve Duo:

```sh
ssh <NetID>@stellar.princeton.edu
```

This is a private controller owned by one Unix account. Do not use `sudo`,
`dnf`, `rpm`, or `systemctl`. Keep the controller, workspaces, Slurm output,
and Jenkins configuration in shared scratch storage visible to compute nodes;
do not use `/tmp`, `/home`, or a project checkout.

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64
export PATH="$JAVA_HOME/bin:$PATH"

mkdir -p "$GKEYLL_CI_ROOT"
java -version
```

The final command must report Java 21 or newer. Stellar's default Java may be
too old. The controller stores its WAR, `jenkins_home`, webroot, temporary
files, logs, and build workspaces below `GKEYLL_CI_ROOT`; scratch purging can
remove this data, so do not treat its credentials or build history as durable.

Download the Jenkins LTS WAR. If Stellar cannot reach the download site,
download it on your laptop and copy it to the same path through SSH.

```sh
cd "$GKEYLL_CI_ROOT"
curl -fL -o jenkins.war https://get.jenkins.io/war-stable/latest/jenkins.war
```

From the reviewed Gkeyll checkout, start the controller (see
`./ci/jenkins/gkeyll-ci.sh -h`):

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64
./ci/jenkins/gkeyll-ci.sh stellar_cpu start
```

The launcher creates detached tmux session `gkeyll_ci` and binds Jenkins only
to `127.0.0.1:8080`. It is safe to disconnect after startup. To inspect or
stop the controller later:

```sh
tmux attach -t gkeyll_ci
tmux ls
tmux kill-session -t gkeyll_ci
```

The controller and its Jenkins agent must run as the same Unix account that
submits Slurm allocations.

## Open Jenkins browser

The CLI is the normal way to run CI, but the browser is needed for first-time
setup and remains useful for inspection. From your laptop, after SSH/Duo
authentication, create a tunnel:

```sh
ssh -N -L 8081:127.0.0.1:8080 <NetID>@stellar.princeton.edu
```

Open `http://localhost:8081`. Choose another first port if 8081 is occupied;
the final `8080` is the remote Jenkins port and normally remains unchanged.

On Jenkins' first start, obtain the unlock password on Stellar:

```sh
cat "$GKEYLL_CI_ROOT/jenkins_home/secrets/initialAdminPassword"
```

Complete initial setup and create an administrator account. Install Pipeline,
Git, Credentials Binding, and Git client if Jenkins did not install them as
suggested plugins.

## Set up Jenkins

### Create a Jenkins API token for the launcher

In Jenkins, sign in as the user who will launch CI builds and create an API
token in that user's security settings. On Stellar, store the Jenkins username
and token in a file owned and readable only by that Unix user:

```sh
umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' \
  > "$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth"
chmod 600 "$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth"
```

This is a Jenkins credential, not the GitHub credential below. The launcher
uses it only against loopback Jenkins and does not disable CSRF protection.

### Create the GitHub credential

Create a classic GitHub PAT with only the `repo:status` scope and a short
expiration. Its owner must have push access to `gkeyllorg/gkeyll`, which GitHub
requires to publish commit statuses. In Jenkins, add it as a **Username with
password** credential, using the owner's GitHub username and the PAT. Give it
an ID such as `gkeyll-github-stellar-cpu`. Organization membership is not
required.

The Pipeline fetches public candidates anonymously; this credential is used
only for authenticated GitHub API requests and status publication. Do not
share it or store Duo secrets, personal SSH keys, or unrelated credentials in
Jenkins. Revoke or replace it when this controller or its owner changes.

### Configure the Jenkins node and global environment

Use the controller's built-in node or a local agent on the Stellar CPU login
side. Give it label `stellar_cpu`. It must use the Unix account that can submit
Slurm jobs.

In **Manage Jenkins → System → Global properties → Environment variables**,
set the following values. Paste an expanded scratch path, not a literal `$USER`.

| Name | Required value |
| --- | --- |
| `GKEYLL_CI_ROOT` | `/scratch/gpfs/<your-user>/gkeyll_ci` |
| `STELLAR_CPU_GITHUB_CREDENTIAL_ID` | GitHub status/API credential ID, e.g. `gkeyll-github-stellar-cpu` |
| `STELLAR_CPU_SLURM_QOS` | A valid CPU QoS for your group |
| `STELLAR_CPU_SLURM_ACCOUNT` | Project account when required; otherwise omit it |
| `STELLAR_CPU_BUILD_JOBS` | Optional login-node compile parallelism; default `3` |
| `STELLAR_CPU_UNIT_TIME` | Optional unit-test allocation limit; default `00:30:00` |
| `STELLAR_CPU_REGRESSION_TIME` | Optional C-regression allocation limit; default `04:00:00` |
| `STELLAR_CPU_REGRESSION_JOBS` | Optional concurrent C test runs; default `4` |
| `STELLAR_CPU_REGRESSION_TEST_TIMEOUT` | Optional per-test limit in seconds; default `900` |
| `STELLAR_CPU_NODE_LABEL` | Optional node label; default `stellar_cpu` |

Do not set a broad global `PATH` to an interactive shell configuration. The
Pipeline initializes Stellar modules for each build and Slurm job. Its
`GKEYLL_CI_ROOT` must exactly match the controller's root.

### Create the one parameterized Pipeline job

Create **New Item → Pipeline** named `gkeyll-ci-stellar_cpu`. Configure
**Pipeline script from SCM** with:

- SCM: Git
- Repository: `https://github.com/gkeyllorg/gkeyll.git`
- Credentials: leave empty; the public repository is fetched anonymously
- Branch specifier: `*/main`
- Script path: `ci/jenkins/Jenkinsfile.stellar_cpu`

The job must load its Pipeline from reviewed `main`, not from a pull request
being tested.

Leave **This project is parameterized** unchecked: the trusted Pipeline
declares the parameters itself. On a newly created job, click **Build Now**
once. The expected initial build stops immediately because no selector is set,
but it registers `CANDIDATE_PR`, `CANDIDATE_REF`, and `BASELINE_REF` with
Jenkins. The page then offers **Build with Parameters**. This initial build
does not submit a Slurm job.

# Launching CI jobs

## CLI launch

After SSH/Duo authentication, run these commands from the reviewed Gkeyll
checkout. The launcher starts the detached controller when necessary and uses
the loopback Jenkins API.

```sh
# Queue a pull-request build and return after Jenkins accepts it.
./ci/jenkins/gkeyll-ci.sh stellar_cpu run --pr 1104

# Compare a branch or commit with a selected baseline and follow its result.
./ci/jenkins/gkeyll-ci.sh stellar_cpu run \
  --candidate-ref feature/new-solver --baseline-ref main --follow
```

The default `run` output includes a Jenkins queue ID. Use it to monitor or
stop work after reconnecting:

```sh
./ci/jenkins/gkeyll-ci.sh stellar_cpu follow --queue 42
./ci/jenkins/gkeyll-ci.sh stellar_cpu status --build 187
./ci/jenkins/gkeyll-ci.sh stellar_cpu active
./ci/jenkins/gkeyll-ci.sh stellar_cpu recent --limit 5
./ci/jenkins/gkeyll-ci.sh stellar_cpu info --build 187
./ci/jenkins/gkeyll-ci.sh stellar_cpu artifact --build 187 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh stellar_cpu abort --queue 42
```

`--follow` streams Jenkins' console through its terminal result and returns
zero only for `SUCCESS`. Ctrl-C stops local monitoring; it does not abort
Jenkins or Slurm work. Jenkins serializes builds for this job, so a later
request waits in its queue.

## Browser launch

Through the SSH tunnel, open `gkeyll-ci-stellar_cpu` and select **Build with
Parameters**. For a pull request, set only `CANDIDATE_PR`. For a direct
comparison, set both `CANDIDATE_REF` and `BASELINE_REF`. Do not combine a PR
number with either reference field.

Jenkins records the selectors and exact commits, posts the
`continuous-integration/jenkins/stellar_cpu` GitHub status, compiles on the
login node, and submits separate unit and C-regression Slurm jobs. Use the
build page to inspect its console and artifacts, or **Abort** to stop it.

# Troubleshooting

## Controller and Jenkins API

If the controller does not become ready, inspect its console and log:

```sh
tmux attach -t gkeyll_ci
tail -n 100 "$GKEYLL_CI_ROOT/logs/jenkins.log"
```

If Jenkins reports a missing job, confirm the job name is
`gkeyll-ci-stellar_cpu`, the API-token file has mode 600, and its Jenkins user
can read and build that job. `recent` and `active` work only after the job
exists and is visible to that user.

## Slurm jobs and controller crashes

To diagnose Slurm work or clean up after a controller crash:

```sh
squeue --me
scancel <jobid>
```

Do not remove a workspace while its Slurm job appears in `squeue`. After it
stops, it is safe to remove the abandoned workspace below
`$GKEYLL_CI_ROOT/workspaces` and restart Jenkins. Jenkins retains 20 recent
build records and archives artifacts before removing normal completed workspaces.

## Numerical regression differences

An unlisted numerical regression difference fails CI. Review a genuine change
and add its accepted difference to `ci/jenkins/expected_regression_diffs.txt`;
do not acknowledge output differences merely to make a build pass.

## Validate manually

For a first-installation failure, validate the same cluster/toolchain workflow
before changing Jenkins. Create disposable candidate and baseline checkouts
below the CI root. Replace `<qos>` and, where required, `<account>` with the
values used for your Jenkins job.

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
mkdir -p "$GKEYLL_CI_ROOT"
cd "$GKEYLL_CI_ROOT"
git clone --branch main --single-branch https://github.com/gkeyllorg/gkeyll.git gkeyll
git clone --branch main --single-branch https://github.com/gkeyllorg/gkeyll.git gkeyll-baseline
```

Build and install the candidate, then submit its unit-test Slurm script:

```sh
cd "$GKEYLL_CI_ROOT/gkeyll"
PREFIX="$PWD/../gkylsoft" ./machines/mkdeps.stellar-intel.sh
PREFIX="$PWD/../gkylsoft" ./machines/configure.stellar-intel.sh
. ./machines/module_load.stellar-intel.sh
make -j32 unit
make -j32 install
sbatch --wait --qos <qos> --nodes 1 --ntasks 1 --cpus-per-task 4 \
  --time 00:30:00 --chdir "$PWD" --export=ALL,CI_WORKSPACE="$PWD" \
  ci/jenkins/slurm-unit-tests.stellar_cpu.sh
```

For PPPL/CIMES, add `--account <account>` to each `sbatch` command. Build and
install the baseline with its own prefix, then configure and compile C
regressions in both checkouts:

```sh
cd "$GKEYLL_CI_ROOT/gkeyll-baseline"
PREFIX="$PWD/gkylsoft" ./machines/mkdeps.stellar-intel.sh
PREFIX="$PWD/gkylsoft" ./machines/configure.stellar-intel.sh
. ../gkeyll/machines/module_load.stellar-intel.sh
make -j32 install
"$PWD/gkylsoft/gkeyll/bin/gkeyll" runregression configure --source-dir "$PWD" --prefix "$PWD/gkylsoft"
"$PWD/gkylsoft/gkeyll/bin/gkeyll" runregression run -c compile

cd "$GKEYLL_CI_ROOT/gkeyll"
. machines/module_load.stellar-intel.sh
make -j32 install
"$PWD/gkylsoft/gkeyll/bin/gkeyll" runregression configure --source-dir "$PWD" --prefix "$PWD/gkylsoft"
"$PWD/gkylsoft/gkeyll/bin/gkeyll" runregression run -c compile
```

Finally submit the execution-only comparison:

```sh
cd "$GKEYLL_CI_ROOT/gkeyll"
sbatch --wait --qos <qos> --nodes 1 --ntasks 1 --cpus-per-task 8 \
  --time 04:00:00 --chdir "$PWD" \
  --export=ALL,CI_WORKSPACE="$PWD",CI_BASELINE_DIR="$GKEYLL_CI_ROOT/gkeyll-baseline",CI_BASELINE_PREFIX="$GKEYLL_CI_ROOT/gkeyll-baseline/gkylsoft",CI_CANDIDATE_PREFIX="$GKEYLL_CI_ROOT/gkeyll/gkylsoft",CI_REGRESSION_JOBS=4,CI_REGRESSION_TEST_TIMEOUT=900 \
  ci/jenkins/slurm-regression-tests.stellar_cpu.sh
```

The baseline's accepted C-regression output is copied into the candidate
results tree and is disposable. Do not use a persistent accepted-output cache.
