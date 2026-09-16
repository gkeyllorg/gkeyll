# Gkeyll Jenkins CI on a Team Workstation

This shared-workstation CI discovers `main` and PRs targeting `main`, and lets
maintainers launch explicit PR or candidate/baseline comparisons.

# Installation

## Install Jenkins

Install Jenkins LTS and Java 21 or newer as a local service. Install Pipeline,
Git, Credentials Binding, Git client, and GitHub Branch Source plugins.

## Open Jenkins browser

Open the workstation's Jenkins URL, complete first-start setup, and create an
administrator account.

## Set up Jenkins

### Create a Jenkins API token for the launcher

Create an API token for the CLI user and save a user-owned mode-600
`jenkins-user:api-token` file, then set `JENKINS_CLI_AUTH_FILE` to it.

### Create the GitHub credential

Create a `gkeyllorg/gkeyll` token with repository/PR read and commit-status
write access. Add it as a Jenkins credential and record its ID.

### Configure the Jenkins node and global environment

Label the node, configure its non-interactive PATH with the toolchain, `cmake`,
and Python/NumPy, then set:

| Name | Meaning |
| --- | --- |
| `TEAM_WORKSTATION_NODE_LABEL` | Workstation node label |
| `TEAM_WORKSTATION_MKDEPS_SCRIPT` | `machines/` dependency script |
| `TEAM_WORKSTATION_CONFIGURE_SCRIPT` | `machines/` configure script |
| `TEAM_WORKSTATION_GITHUB_CREDENTIAL_ID` | GitHub credential ID |
| `TEAM_WORKSTATION_BUILD_JOBS` | Optional; default `3` |
| `TEAM_WORKSTATION_REGRESSION_JOBS` | Optional; default `1` |

### Create the multibranch Pipeline job

Create Multibranch Pipeline `gkeyll-ci-team-workstation` from GitHub source
`gkeyllorg/gkeyll`. Discover `main` and pull requests, exclude ordinary
branches that are also PRs, use script path
`ci/jenkins/jenkinsfile.team_workstation`, and enable two-minute scans. Run an
initial scan to create the trusted `main` child used by the CLI.

# Launching CI jobs

## CLI launch

```sh
export JENKINS_CLI_AUTH_FILE="$HOME/.config/gkeyll/jenkins-cli.auth"
./ci/jenkins/gkeyll-ci.sh team scan
./ci/jenkins/gkeyll-ci.sh team run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh team active
./ci/jenkins/gkeyll-ci.sh team abort --build 42
```

`scan` discovers `main` and eligible PR jobs, then reports builds Jenkins
scheduled. The client never starts Jenkins or tmux.

## Browser launch

Use the multibranch job and its `main` child for explicit parameterized runs;
PR child jobs are used for automatic PR builds.

# Troubleshooting

## Jenkins API and discovery

Confirm the token file mode and job permissions. If `main` is missing, run
`team scan` and inspect the multibranch job's indexing log.

## Numerical regression differences

Expected numerical changes require a reviewed entry in
`ci/jenkins/expected_regression_diffs.txt`; unlisted differences fail CI.
