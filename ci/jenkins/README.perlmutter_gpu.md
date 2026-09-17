# Gkeyll Jenkins CI on NERSC Perlmutter GPU

This private, manually triggered CUDA CI builds candidate and baseline CUDA/NCCL
installations on the login node, then submits unit and C-regression work to GPU
nodes. The trusted Pipeline never comes from the candidate PR. CUDA unit tests
run automatically; GPU-capable C regressions are compared with the CUDA/NCCL
baseline built by the Pipeline.
A separate manifest-selected C MPI lane uses four GPU ranks for Vlasov,
gyrokinetic, and PKPM, plus a moments CPU MPI test.

# Installation

## Install Jenkins

Use a controller only in an authenticated Perlmutter session and as permitted
by NERSC policy. Set a project scratch root visible to compute nodes; do not
use `/tmp` or retain durable credentials only there.

```sh
export GKEYLL_CI_ROOT=/pscratch/sd/<first-letter>/<username>/gkeyll_ci
export JAVA_HOME=<java-21-or-newer-installation>
mkdir -p "$GKEYLL_CI_ROOT"
java -version
cd "$GKEYLL_CI_ROOT"
curl -fL -o jenkins.war https://get.jenkins.io/war-stable/latest/jenkins.war
```

From the reviewed Gkeyll checkout, start the controller (see
`./ci/jenkins/gkeyll-ci.sh -h`):

```sh
export GKEYLL_CI_ROOT=/pscratch/sd/<first-letter>/<username>/gkeyll_ci
export JAVA_HOME=<java-21-or-newer-installation>
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu start
```

It creates detached tmux session `gkeyll_ci` and binds only to loopback. Use
`tmux attach -t gkeyll_ci` to inspect it or `tmux kill-session -t gkeyll_ci`
to stop it.

Verify the private listener after startup:

```sh
curl --fail --output /dev/null http://127.0.0.1:8080/login
ss -ltn | grep '127.0.0.1:8080'
```

## Open Jenkins browser

From your laptop, tunnel a local port to the controller and open the resulting
local URL:

```sh
ssh -N -o ExitOnForwardFailure=yes \
  -L 127.0.0.1:8084:127.0.0.1:8080 <username>@perlmutter.nersc.gov
```

Open `http://127.0.0.1:8084`. The first `8084` is the local browser port; the
final `8080` is the remote Jenkins port and normally remains unchanged. Choose
another unused local port if necessary. On first start, read
`$GKEYLL_CI_ROOT/jenkins_home/secrets/initialAdminPassword`, create an admin
account, and install Pipeline, Git, Credentials Binding, Git client, and
GitHub plugins.

## Set up Jenkins

### Create a Jenkins API token for the launcher

Create an API token for the Jenkins user that will launch builds, then save it
on Perlmutter with mode 600:

```sh
umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' > "$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth"
chmod 600 "$GKEYLL_CI_ROOT/jenkins_home/jenkins-cli.auth"
```

### Create the GitHub credential

Create a classic GitHub PAT with only the `repo:status` scope and a short
expiration. Its owner must have push access to `gkeyllorg/gkeyll`, which GitHub
requires to publish commit statuses. In **Manage Jenkins → Credentials**, add
it to this controller as a **Username with password** credential: use the
owner's GitHub username and the PAT as the password, then record its ID.
Organization membership is not required.

The Pipeline fetches public candidates anonymously; this PAT is used only for
authenticated GitHub API requests and status publication. Do not share it or
store it outside this controller. Revoke or replace it when the controller or
its owner changes.

### Configure the Jenkins node and global environment

Use a login-side node/agent labelled `perlmutter_gpu`. In **Manage Jenkins →
System → Global properties → Environment variables**, set:

| Name | Value |
| --- | --- |
| `GKEYLL_CI_ROOT` | Expanded shared project-scratch root |
| `PERLMUTTER_GPU_GITHUB_CREDENTIAL_ID` | GitHub status/API credential ID |
| `PERLMUTTER_GPU_SLURM_ACCOUNT` | Required NERSC project/account |
| `PERLMUTTER_GPU_NODE_LABEL` | Optional; default `perlmutter_gpu` |
| `PERLMUTTER_GPU_SLURM_QOS` | Optional; default `shared` |
| `PERLMUTTER_GPU_BUILD_JOBS` | Optional; default `3` |
| `PERLMUTTER_GPU_UNIT_TIME` | Optional; default `00:30:00` |
| `PERLMUTTER_GPU_REGRESSION_TIME` | Optional; default `04:00:00` |
| `PERLMUTTER_GPU_REGRESSION_JOBS` | Optional; default `4` |
| `PERLMUTTER_GPU_REGRESSION_TEST_TIMEOUT` | Optional; default `900` |

### Create the one parameterized Pipeline job

Create **New Item → Pipeline** named `gkeyll-ci-perlmutter_gpu`. Use
**Pipeline script from SCM** with repository `https://github.com/gkeyllorg/gkeyll.git`,
no SCM credential, branch `*/main`, and script path
`ci/jenkins/Jenkinsfile.perlmutter_gpu`. The public repository is fetched
anonymously; use node label `perlmutter_gpu`.

Leave **This project is parameterized** unchecked. Click **Build Now** once to
register the Pipeline parameters; its expected empty-selector failure submits
no Slurm job. Do not let a PR supply this Pipeline.

# Launching CI jobs

## CLI launch

```sh
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu run --pr 1234
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu run --candidate-ref feature --baseline-ref main --follow
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu follow --queue 42
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu active
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu info --build 42
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu artifact --build 42 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu abort --build 42
```

## Browser launch

Open the job through the tunnel, select **Build with Parameters**, and set
either `CANDIDATE_PR` or both `CANDIDATE_REF` and `BASELINE_REF`.

# Troubleshooting

## Controller and Slurm jobs

Inspect `tmux attach -t gkeyll_ci`, `$GKEYLL_CI_ROOT/logs/jenkins.log`, and
`squeue --me`. After a controller crash, use `scancel <jobid>` before removing
an abandoned workspace.

## Validate manually

Before changing Jenkins, build a disposable main checkout with
`machines/mkdeps.perlmutter.gpu.sh` and `machines/configure.perlmutter.gpu.sh`,
then submit `slurm-unit-tests.perlmutter_gpu.sh` using your account, `shared`
QoS, `--constraint gpu`, 32 CPUs, and one GPU. Build a separate baseline prefix
and compile C regressions in both trees before submitting the regression
Slurm script. This separates NERSC toolchain/allocation failures from Jenkins setup.
