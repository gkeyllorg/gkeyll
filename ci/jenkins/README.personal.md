# Gkeyll Jenkins CI on a Personal Computer

This private CI lets the computer owner explicitly test a Gkeyll PR or
candidate/baseline comparison. It does not poll GitHub or execute unselected
contributor code.

# Installation

## Install Jenkins

Install Jenkins LTS and Java 21 or newer with the normal local service
mechanism (for example `brew install jenkins-lts` and `brew services start
jenkins-lts` on macOS). Install Pipeline, Git, Credentials Binding, Git client,
and GitHub Branch Source plugins.

## Open Jenkins browser

Open the local Jenkins URL, normally `http://127.0.0.1:8080`, complete first
start setup, and create an administrator account.

## Set up Jenkins

### Create a Jenkins API token for the launcher

Create an API token for the user who will run the client and save it locally:

```sh
mkdir -p "$HOME/.config/gkeyll"; umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' > "$HOME/.config/gkeyll/jenkins-cli.auth"
chmod 600 "$HOME/.config/gkeyll/jenkins-cli.auth"
```

### Create the GitHub credential

Create a fine-grained `gkeyllorg/gkeyll` token with Contents and Pull requests
read access plus Commit statuses read/write. Add it as a Jenkins **Username
with password** credential and record its ID.

### Configure the Jenkins node and global environment

Label the build node and set its non-interactive PATH to compilers, `cmake`,
and Python with NumPy. Set these global environment variables:

| Name | Meaning |
| --- | --- |
| `PERSONAL_NODE_LABEL` | Local Jenkins build-node label |
| `PERSONAL_MKDEPS_SCRIPT` | `machines/` dependency script |
| `PERSONAL_CONFIGURE_SCRIPT` | `machines/` configure script |
| `PERSONAL_GITHUB_CREDENTIAL_ID` | GitHub credential ID |
| `PERSONAL_BUILD_JOBS` | Optional; default `3` |
| `PERSONAL_REGRESSION_JOBS` | Optional; default `1` |

### Create the one parameterized Pipeline job

Create Pipeline `gkeyll-ci-personal` from SCM repository
`https://github.com/gkeyllorg/gkeyll.git`, branch `*/main`, and script path
`ci/jenkins/jenkinsfile.personal`. Do not let a selected PR provide its
Pipeline. Run it once without selectors to register parameters.

# Launching CI jobs

## CLI launch

```sh
export JENKINS_CLI_AUTH_FILE="$HOME/.config/gkeyll/jenkins-cli.auth"
./ci/jenkins/gkeyll-ci.sh personal run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh personal run --candidate-ref feature/new-solver --baseline-ref main
./ci/jenkins/gkeyll-ci.sh personal active
./ci/jenkins/gkeyll-ci.sh personal abort --build 42
```

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
