# Gkeyll Jenkins CI on a Personal Computer

To reuse an existing installation’s dependency paths and `config.mak`, see
[Reusing installed dependencies](README.md#reusing-installed-dependencies).

This private CI lets the computer owner explicitly test a Gkeyll PR or
candidate/baseline comparison. It does not poll GitHub or execute unselected
contributor code.

# Installation

## Install Jenkins

Install Jenkins LTS and Java 21 or newer with the normal local service
mechanism (for example `brew install jenkins-lts` and `brew services start
jenkins-lts` on macOS).

## Open Jenkins browser

Open the local Jenkins URL, normally `http://127.0.0.1:8080`, complete first
start setup, and create an administrator account. During initial setup, install
the Pipeline, Git, Credentials Binding, Git client, and GitHub Branch Source
plugins.

## Set up Jenkins

### Create a Jenkins API token for the launcher

Create an API token for the user who will run the client and save it locally:

```sh
mkdir -p "$HOME/.config/gkeyll/jenkins"; umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' > "$HOME/.config/gkeyll/jenkins/personal.auth"
chmod 600 "$HOME/.config/gkeyll/jenkins/personal.auth"
```

### Create the GitHub credential

The [shared reporter](README.md#reporting-credentials) uses this existing
username/token credential for commit statuses and PR/commit comments.

Create a classic GitHub PAT with only the `public_repo` scope and a short
expiration. That scope covers commit statuses and pull-request or commit
comments on the public repository; a fine-grained token needs **Commit
statuses: write**, **Pull requests: write**, and **Contents: write** (commit
comments) on `gkeyllorg/gkeyll` instead. Its owner must have push access to
`gkeyllorg/gkeyll`, which GitHub requires to publish commit statuses. In **Manage Jenkins → Credentials**, add
it to this controller as a **Username with password** credential: use the
owner's GitHub username and the PAT as the password, then record its ID.
Organization membership is not required.

The Pipeline reads public source anonymously; this PAT is used only for
authenticated GitHub API requests and status publication. Do not share it or
store it outside this controller. Revoke or replace it when the controller or
its owner changes.

### Configure the Jenkins node and global environment

Label the build node and set its non-interactive PATH to compilers, `cmake`,
and Python with NumPy. Set these global environment variables:

| Name | Meaning |
| --- | --- |
| `PERSONAL_NODE_LABEL` | Local Jenkins build-node label |
| `PERSONAL_MKDEPS_SCRIPT` | `machines/` dependency script |
| `PERSONAL_CONFIGURE_SCRIPT` | `machines/` configure script |
| `PERSONAL_GITHUB_CREDENTIAL_ID` | GitHub status/API credential ID |
| `PERSONAL_BUILD_JOBS` | Optional compilation workers, including C regressions; default `3` |
| `PERSONAL_REGRESSION_JOBS` | Optional concurrent serial C test runs; defaults to `PERSONAL_BUILD_JOBS` (default `3`) |
| `PERSONAL_MPI_HOME` | Optional MPI installation path for both trees; default is each tree's `gkylsoft/openmpi` |
| `PERSONAL_MPIEXEC` | Optional launcher override for both trees; default is `bin/mpiexec` under the selected MPI installation |
| `PERSONAL_STATUS_CONTEXT` | Set explicitly for queue reporting, e.g. `continuous-integration/jenkins/personal-<hostname>`; use a distinct context for each computer |
| `GKEYLL_CI_ROOT` | Required persistent writable root for baseline and candidate source, build, and regression data |
| `GKEYLL_CI_TRUSTED_REF` | Default branch or full SHA for the personal Jenkinsfile and CI helpers; default `main`. A run's `CI_REF` overrides it |

The selected dependency script must pass `--build-adas=yes` to
`install-deps/mkdeps.sh` so ADAS data is available before unit tests run.

The workflow sets `MPI_HOME` to `PERSONAL_MPI_HOME` while building each tree.
If `PERSONAL_MPI_HOME` is unset or blank, it uses that tree's SHA-local
`gkylsoft/openmpi`. Set `PERSONAL_MPI_HOME` (for example, `/opt/openmpi`) to use
an existing MPI installation for both builds and parallel regressions.
An ambient `MPI_HOME` does not select the CI installation.

### Create the one parameterized Pipeline job

Create Pipeline `gkeyll-ci-personal` from SCM repository
`https://github.com/gkeyllorg/gkeyll.git`, branch `*/main`, and script path
`ci/jenkins/jenkinsfile.personal`. Leave the SCM **Credentials** field empty:
Gkeyll is public and the status credential is not a Git checkout credential.
Do not let a selected PR provide its Pipeline. Run it once without selectors
to register parameters.

The SCM Jenkinsfile supplies the initial loader. It resolves `CI_REF` (or the
default above) to one commit, then loads `jenkinsfile.personal` and its CI
helpers from that commit. Candidate selection does not select the CI
implementation. When staging a change to the loader itself, point the job's
SCM branch and refspec to your CI feature branch too.

Install the [controller queue listener](README.md#controller-installation) to
report pending before an executor is available and cancel superseded queued PR
commits. Set `PERSONAL_STATUS_CONTEXT` globally to the context already used by
this machine, so queued and final statuses update the same GitHub check.

# Launching CI jobs

## CLI launch

```sh
./ci/jenkins/gkeyll-ci.sh personal run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh personal run --candidate-ref feature/new-solver --baseline-ref main
./ci/jenkins/gkeyll-ci.sh personal active
./ci/jenkins/gkeyll-ci.sh personal info --build 42
./ci/jenkins/gkeyll-ci.sh personal artifact --build 42 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh personal abort --build 42
```

To compare `branch1` with `main` using the CI implementation from `branch2`:

```sh
./ci/jenkins/gkeyll-ci.sh personal run \
  --candidate-ref branch1 --baseline-ref main --ci-ref branch2 --follow
```

`--ci-ref` also accepts a full commit SHA and can be combined with `--pr`.
The branch must exist in `gkeyllorg/gkeyll`. This selector executes CI code with
the job's normal access, so select a CI revision you trust. Team CI continues to
use its administrator-controlled `TEAM_WORKSTATION_TRUSTED_CI_REF`.

The console prints the selected ref and full Jenkinsfile commit before the
build starts. `ci-pipeline-source.txt` and `ci-trusted-ci-commit.txt` retain that
selection, including when an older selected Jenkinsfile deletes its workspace
or fails. Current reporting code also puts Jenkinsfile, reporting-tool, and
regression-checker commits in the GitHub report. Branch movement during the run
does not change the pinned CI implementation. Older Jenkinsfiles must support
`CI_NODE_ALREADY_ALLOCATED` and `CI_TRUSTED_CI_COMMIT` to run under this loader.

Set `JENKINS_CLI_AUTH_FILE` only to use a credential file at a different path.

## Browser launch

Open `gkeyll-ci-personal`, select **Build with Parameters**, and set either a
PR number or both candidate and baseline references. Optionally set **CI_REF**
to choose the Jenkinsfile and CI helpers independently.

# Troubleshooting

## Jenkins API and builds

Confirm the token file is user-owned and mode 600, Jenkins is running, and the
user can read/build `gkeyll-ci-personal`. `follow` and `status` accept queue or
build IDs; `abort --queue` cancels waiting work and `abort --build` stops it.
CI rejects a candidate that does not contain its baseline before building. Use
`--allow-behind-candidate` only for an intentional historical comparison.

## Numerical regression differences

Expected numerical changes require a reviewed, new or updated entry in
`ci/jenkins/expected_regression_diffs.txt`; unlisted differences fail CI.
Entries unchanged from the baseline are inert, so stale entries may safely be
removed in any later PR. Update the reason on a baseline entry to acknowledge
a new intentional change for that same test.
