# Gkeyll Jenkins CI on a Personal Computer

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

Create a classic GitHub PAT with only the `repo:status` scope and a short
expiration. Its owner must have push access to `gkeyllorg/gkeyll`, which GitHub
requires to publish commit statuses. In **Manage Jenkins → Credentials**, add
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
| `PERSONAL_BUILD_JOBS` | Optional; default `3` |
| `PERSONAL_REGRESSION_JOBS` | Optional; default `1` |
| `PERSONAL_MPIEXEC` | Optional launcher override for both trees; otherwise each tree uses its own `gkylsoft/openmpi/bin/mpiexec` |

### Create the one parameterized Pipeline job

Create Pipeline `gkeyll-ci-personal` from SCM repository
`https://github.com/gkeyllorg/gkeyll.git`, branch `*/main`, and script path
`ci/jenkins/jenkinsfile.personal`. Leave the SCM **Credentials** field empty:
Gkeyll is public and the status credential is not a Git checkout credential.
Do not let a selected PR provide its Pipeline. Run it once without selectors
to register parameters.

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

Set `JENKINS_CLI_AUTH_FILE` only to use a credential file at a different path.

## Browser launch

Open `gkeyll-ci-personal`, select **Build with Parameters**, and set either a
PR number or both candidate and baseline references.

# Troubleshooting

## Jenkins API and builds

Confirm the token file is user-owned and mode 600, Jenkins is running, and the
user can read/build `gkeyll-ci-personal`. `follow` and `status` accept queue or
build IDs; `abort --queue` cancels waiting work and `abort --build` stops it.

## Numerical regression differences

Expected numerical changes require a reviewed entry in
`ci/jenkins/expected_regression_diffs.txt`; unlisted differences fail CI.
