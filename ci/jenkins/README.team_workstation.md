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

Create an API token for the CLI user and save it in the user-owned mode-600
default credential file:

```sh
mkdir -p "$HOME/.config/gkeyll/jenkins"; umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' > "$HOME/.config/gkeyll/jenkins/team-workstation.auth"
chmod 600 "$HOME/.config/gkeyll/jenkins/team-workstation.auth"
```

### Create the GitHub credential

Designate a workstation CI maintainer who has push access to
`gkeyllorg/gkeyll`. That maintainer creates a classic GitHub PAT with only the
`repo:status` scope and a short expiration. In **Manage Jenkins →
Credentials**, add it as a **Username with password** credential: use the
maintainer's GitHub username and the PAT as the password, then record its ID.
Organization membership is not required.

The Pipeline fetches public source anonymously; this PAT is used only for
authenticated GitHub API requests and status publication. Do not share it.
Revoke or replace it when the maintainer or workstation ownership changes.

### Configure the Jenkins node and global environment

Label the node, configure its non-interactive PATH with the toolchain, `cmake`,
and Python/NumPy, then set:

| Name | Meaning |
| --- | --- |
| `TEAM_WORKSTATION_NODE_LABEL` | Workstation node label |
| `TEAM_WORKSTATION_MKDEPS_SCRIPT` | `machines/` dependency script |
| `TEAM_WORKSTATION_CONFIGURE_SCRIPT` | `machines/` configure script |
| `TEAM_WORKSTATION_GITHUB_CREDENTIAL_ID` | GitHub status/API credential ID |
| `TEAM_WORKSTATION_BUILD_JOBS` | Optional; default `3` |
| `TEAM_WORKSTATION_REGRESSION_JOBS` | Optional; default `1` |

### Create the multibranch Pipeline job

Create Multibranch Pipeline `gkeyll-ci-team-workstation` from GitHub source
`gkeyllorg/gkeyll`. Discover `main` and pull requests, exclude ordinary
branches that are also PRs, use script path
`ci/jenkins/jenkinsfile.team_workstation`, and enable two-minute scans. Run an
initial scan to create the trusted `main` child used by the CLI. Configure the
GitHub source anonymously; do not use the status PAT for source discovery.

# Launching CI jobs

## CLI launch

```sh
./ci/jenkins/gkeyll-ci.sh team scan
./ci/jenkins/gkeyll-ci.sh team run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh team active
./ci/jenkins/gkeyll-ci.sh team abort --build 42
```

Set `JENKINS_CLI_AUTH_FILE` only to use a credential file at a different path.

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
