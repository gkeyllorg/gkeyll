# Gkeyll Jenkins CI on Princeton Stellar CPU

To reuse an existing installation’s dependency paths and `config.mak`, see
[Reusing installed dependencies](README.dependencies.md).

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

Verify the private listener after startup:

```sh
curl --fail --output /dev/null http://127.0.0.1:8080/login
ss -ltn | grep '127.0.0.1:8080'
```

The controller and its Jenkins agent must run as the same Unix account that
submits Slurm allocations.

## Open Jenkins browser

The CLI is the normal way to run CI, but the browser is needed for first-time
setup and remains useful for inspection. From your laptop, after SSH/Duo
authentication, create a tunnel:

```sh
ssh -N -o ExitOnForwardFailure=yes \
  -L 127.0.0.1:8083:127.0.0.1:8080 <NetID>@stellar.princeton.edu
```

Open `http://127.0.0.1:8083`. The first `8083` is the local browser port; the
final `8080` is the remote Jenkins port and normally remains unchanged. Choose
another unused local port if necessary.

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

The [shared reporter](README.reporting.md#reporting-credentials) uses this existing
username/token credential for statuses and report comments on the tested commit.

Create a classic GitHub PAT with only the `public_repo` scope and a short
expiration. That scope covers commit statuses and commit comments on the public
repository; a fine-grained token needs **Commit
statuses: write** and **Contents: write** (to create and update commit
comments) on `gkeyllorg/gkeyll` instead. Its owner must have push access to
`gkeyllorg/gkeyll`, which GitHub requires to publish commit statuses. In Jenkins, add it as a **Username with
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
| `STELLAR_CPU_BUILD_JOBS` | Optional login-node compilation workers, including C regressions; default `8` |
| `STELLAR_CPU_UNIT_TIME` | Optional unit-test allocation limit; default `00:30:00` |
| `STELLAR_CPU_REGRESSION_TIME` | Optional C-regression allocation limit; default `04:00:00` |
| `STELLAR_CPU_REGRESSION_JOBS` | Optional regression worker budget (one per serial test or MPI rank); default `16` |
| `STELLAR_CPU_REGRESSION_TEST_TIMEOUT` | Optional per-test limit in seconds; default `900` |
| `STELLAR_CPU_NODE_LABEL` | Optional node label; default `stellar_cpu` |
| `PERL5LIB` | Shared Perl library path when compute nodes lack `Time::HiRes`; see [runtime verification](#verify-the-compute-node-perl-runtime) |
| `GKEYLL_CI_TRUSTED_REF` | Optional reviewed branch or full SHA for reporting, regression checking, and cache tools; default `main`. Set it only while staging a CI change |
| `CI_JENKINS_READ_AUTH_FILE` | Optional agent-readable Jenkins `user:api-token` file; default `$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth`; empty for anonymous reads |
| `CI_JENKINS_LOCAL_URL` | Optional controller URL reachable from the login agent; default `http://127.0.0.1:8080` |
| `GKEYLL_CI_BUILDS_TO_KEEP` | Positive count of Jenkins builds, artifacts, and entire run directories to retain per job; default `2` |
| `GKEYLL_CI_RETAIN_RUNS` | Default `0` enables bounded retention after Jenkins/Slurm checks; `1` retains unlimited disk history |

Do not set a broad global `PATH` to an interactive shell configuration. The
Pipeline initializes Stellar modules for each build and Slurm job. Its
`GKEYLL_CI_ROOT` must exactly match the controller's root.
It retains a separate directory for each build under `runs/stellar-cpu/`;
see [persistent regression data](README.storage.md). After an intentional in-place toolchain or dependency upgrade, change
`GKEYLL_CI_CACHE_REVISION` to force baseline regeneration in either dependency mode.
Jenkins defaults to retaining two builds, their artifacts, and their entire run
directories; `GKEYLL_CI_BUILDS_TO_KEEP` overrides all three counts. Older inactive
run trees, including dependencies and archives, are deleted at startup and
completion. Active runs remain protected.

Serial regressions request one CPU per worker: the default is 16 concurrent
tests on 16 CPUs. MPI regressions request one CPU per rank slot (at least four
slots) and reserve each test's manifest rank count until it finishes. With the
default 16-slot budget, four four-rank tests can run concurrently; available
slots are reused as tests finish. These settings apply to both baseline and
candidate runs.

Parallel regressions disable Open MPI's HCOLL component and select tuned
recursive-doubling allreduce (`coll_tuned_allreduce_algorithm=3`) for both baseline
creation and candidate checks. This fixes the rank-based summation order for a
given rank count; reduction roundoff can otherwise affect integrated diagnostics
and feed back into adaptive sources. The parallel Slurm payload is included in
the baseline-cache fingerprint, so changing this policy regenerates the baseline
once. Existing baselines made with another reduction policy must not be reused.

### Verify the compute-node Perl runtime

The regression runner uses Perl's `Time::HiRes` module for execution timings.
Stellar login and compute nodes can have different Perl packages: a successful
login-node check does not establish that this module is available in Slurm.
Verify it before the first CI build and after runtime upgrades, using the same
account, QoS, modules, and `PERL5LIB` as Jenkins.

From the reviewed checkout, create and submit this preflight script. Replace
`<qos>` and `<account>` with the Jenkins settings (omit `--account` if unused):

```sh
mkdir -p "$GKEYLL_CI_ROOT/tmp" "$GKEYLL_CI_ROOT/logs"
cat > "$GKEYLL_CI_ROOT/tmp/perl-preflight.sh" <<'EOF'
#!/bin/bash -l
set -euo pipefail
. machines/module_load.stellar-intel.sh
hostname
perl -V:version -V:archname
perl -MTime::HiRes=time -e 'printf "%.6f\n", time'
EOF
sbatch --wait --account <account> --qos <qos> \
  --nodes 1 --ntasks 1 --cpus-per-task 1 --time 00:05:00 \
  --chdir "$PWD" --export=ALL \
  --output "$GKEYLL_CI_ROOT/logs/perl-preflight-%j.log" \
  "$GKEYLL_CI_ROOT/tmp/perl-preflight.sh"
```

Require a successful job and a numeric timestamp in its log. If it reports
`Can't locate Time/HiRes.pm`, install `Time::HiRes` for the compute node's Perl
in shared storage, such as `$GKEYLL_CI_ROOT/runtime/perl5`. Do not install it
only on the login node or in a node-local temporary directory.

For matching system Perl versions and architectures, the login node's existing
module can be copied into this account-local directory. Compare the preflight's
`version` and `archname` with `perl -V:version -V:archname` on the login node
first. The following paths apply to Stellar's system Perl 5.26.3 installation;
use a compatible module installation if the paths or runtime differ:

```sh
mkdir -p "$GKEYLL_CI_ROOT/runtime/perl5/Time" \
  "$GKEYLL_CI_ROOT/runtime/perl5/auto/Time/HiRes"
cp /usr/lib64/perl5/vendor_perl/Time/HiRes.pm \
  "$GKEYLL_CI_ROOT/runtime/perl5/Time/"
cp /usr/lib64/perl5/vendor_perl/auto/Time/HiRes/HiRes.so \
  "$GKEYLL_CI_ROOT/runtime/perl5/auto/Time/HiRes/"
export PERL5LIB="$GKEYLL_CI_ROOT/runtime/perl5${PERL5LIB:+:$PERL5LIB}"
```

Repeat the Slurm preflight with this environment. Once it passes, set Jenkins's
**global environment variable** `PERL5LIB` to the same expanded absolute path
(or colon-separated paths if preserving an existing value). A terminal-only
`export` does not update an already running Jenkins controller. The Pipeline
passes its environment to Slurm with `--export=ALL`. Revalidate the module,
including its compiled `HiRes.so`, after a system Perl upgrade.

### Create the one parameterized Pipeline job

Install the [controller queue listener](README.reporting.md#controller-installation) in
`$GKEYLL_CI_ROOT/jenkins_home/init.groovy.d/` to report pending while jobs wait
and cancel superseded queued PR commits. It uses the global
`STELLAR_CPU_GITHUB_CREDENTIAL_ID` and requires no Slurm allocation.

Create **New Item → Pipeline** named `gkeyll-ci-stellar_cpu`. Configure
**Pipeline script from SCM** with:

- SCM: Git
- Repository: `https://github.com/gkeyllorg/gkeyll.git`
- Credentials: leave empty; the public repository is fetched anonymously
- Branch specifier: `*/main`
- Script path: `ci/jenkins/jenkinsfile.stellar_cpu`

The job must load its Pipeline from reviewed `main`, not from a pull request
being tested.

Leave **This project is parameterized** unchecked: the trusted Pipeline
declares the parameters itself. On a newly created job, click **Build Now**
once. The expected initial build stops immediately because no selector is set,
but it registers `CANDIDATE_PR`, `CANDIDATE_REF`, `BASELINE_REF`, and
`ALLOW_BEHIND_CANDIDATE` with
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

## Candidate freshness

Before dependency builds or Slurm submission, CI resolves both references to
commits and rejects a candidate that does not contain its baseline. Merge or
rebase the baseline into the candidate before rerunning. For an intentional
historical comparison, use `--allow-behind-candidate` or enable `ALLOW_BEHIND_CANDIDATE` in
**Build with Parameters**; the resulting artifact
and GitHub report record the override.

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
`$GKEYLL_CI_ROOT/workspaces` and restart Jenkins. Jenkins archives artifacts
before removing completed workspaces; see [retention settings](README.storage.md#retention).

## Every regression reported as a crash

Inspect the per-test `runlog` in the retained regression databases. Missing
`Time::HiRes` can produce empty `__START__:`/`__END__:` markers and Perl errors
after a successful `__EXIT__:0`. The runner's exit-marker parsing can then
misclassify successful executions as crashes, with no completed timings.
Passing unit tests do not validate this regression-runner dependency.

Follow [runtime verification](#verify-the-compute-node-perl-runtime), persist
the repaired environment in Jenkins, and rerun the PR. Do not acknowledge
these failures in `expected_regression_diffs.txt`; numerical comparisons still
need to run successfully after the runtime is repaired.

## Numerical regression differences

An unlisted numerical regression difference fails CI. Review a genuine change
and add a new or updated accepted difference to
`ci/jenkins/expected_regression_diffs.txt`; do not acknowledge output
differences merely to make a build pass. Entries unchanged from the baseline
are inert, so stale entries may safely be removed in any later PR. Updating a
baseline entry's reason acknowledges a new intentional change for that test.

## Validate manually

For a first-installation failure, validate the same cluster/toolchain workflow
before changing Jenkins. Create disposable candidate and baseline checkouts
below the CI root. Replace `<qos>` and, where required, `<account>` with the
values used for your Jenkins job.

```sh
export GKEYLL_CI_ROOT=/scratch/gpfs/$USER/gkeyll_ci
mkdir -p "$GKEYLL_CI_ROOT"
cd "$GKEYLL_CI_ROOT"
sha="$(git ls-remote https://github.com/gkeyllorg/gkeyll.git refs/heads/main | cut -f1)"
candidate_root="$GKEYLL_CI_ROOT/manual-candidate/$sha"
baseline_root="$GKEYLL_CI_ROOT/manual-baseline/$sha"
mkdir -p "$candidate_root" "$baseline_root"
git clone --branch main --single-branch https://github.com/gkeyllorg/gkeyll.git "$candidate_root/gkeyll"
git clone --branch main --single-branch https://github.com/gkeyllorg/gkeyll.git "$baseline_root/gkeyll"
test "$(git -C "$candidate_root/gkeyll" rev-parse HEAD)" = "$sha"
test "$(git -C "$baseline_root/gkeyll" rev-parse HEAD)" = "$sha"
```

Build and install the candidate, then submit its unit-test Slurm script:

```sh
cd "$candidate_root/gkeyll"
PREFIX="$PWD/../gkylsoft" ./machines/mkdeps.stellar-intel.sh
PREFIX="$PWD/../gkylsoft" ./machines/configure.stellar-intel.sh
. ./machines/module_load.stellar-intel.sh
make -j32 unit
make -j32 install
sbatch --wait --qos <qos> --nodes 1 --ntasks 1 --cpus-per-task 4 \
  --time 00:30:00 --chdir "$PWD" --export=ALL,CI_WORKSPACE="$PWD",CI_CANDIDATE_DIR="$PWD" \
  ci/jenkins/slurm-unit-tests.stellar_cpu.sh
```

For PPPL/CIMES, add `--account <account>` to each `sbatch` command. Build and
install the baseline with its own prefix, then configure and compile C
regressions in both checkouts:

```sh
cd "$baseline_root/gkeyll"
PREFIX="$PWD/../gkylsoft" ./machines/mkdeps.stellar-intel.sh
PREFIX="$PWD/../gkylsoft" ./machines/configure.stellar-intel.sh
. machines/module_load.stellar-intel.sh
make -j32 install
"$PWD/../gkylsoft/gkeyll/bin/gkeyll" runregression configure --source-dir "$PWD" --prefix "$PWD/../gkylsoft"
"$PWD/../gkylsoft/gkeyll/bin/gkeyll" runregression run -c compile

cd "$candidate_root/gkeyll"
. machines/module_load.stellar-intel.sh
make -j32 install
"$PWD/../gkylsoft/gkeyll/bin/gkeyll" runregression configure --source-dir "$PWD" --prefix "$PWD/../gkylsoft"
"$PWD/../gkylsoft/gkeyll/bin/gkeyll" runregression run -c compile
```

Finally submit the execution-only comparison:

```sh
cd "$candidate_root/gkeyll"
sbatch --wait --qos <qos> --nodes 1 --ntasks 1 --cpus-per-task 8 \
  --time 04:00:00 --chdir "$PWD" \
  --export=ALL,CI_WORKSPACE="$PWD",CI_CANDIDATE_DIR="$PWD",CI_BASELINE_DIR="$baseline_root/gkeyll",CI_BASELINE_PREFIX="$baseline_root/gkylsoft",CI_CANDIDATE_PREFIX="$candidate_root/gkylsoft",CI_TRUSTED_CHECKER="$baseline_root/gkeyll/ci/jenkins/check_regression_results.lua",CI_REGRESSION_MODE=candidate-check,CI_REGRESSION_JOBS=4,CI_REGRESSION_TEST_TIMEOUT=900 \
  ci/jenkins/slurm-regression-tests.stellar_cpu.sh
```

The baseline's accepted C-regression output is linked into the candidate
results tree for checking. These manual checkouts are disposable.

Run compaction is disabled by default. Set `GKEYLL_CI_COMPACT_RUNS=1` to
archive completed run trees; see [compaction](README.storage.md#optional-run-compaction).
