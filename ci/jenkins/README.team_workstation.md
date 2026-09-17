# Gkeyll Jenkins CI on a Team Workstation

This CI automatically discovers `main` and pull requests targeting `main`.
Each PR is built in its own multibranch child job, but no PR supplies the
Pipeline that controls the build: Jenkins runs a small controller-owned
bootstrap, which checks out reviewed `main` and loads the workflow from that
separate checkout. The workflow then separately fetches the PR candidate and
its baseline.

Candidate code is still arbitrary code under test. Run the Jenkins agent as a
dedicated, low-privilege account or in an isolated environment, with no
personal SSH keys, Jenkins home access, or unrelated credentials.

# Installation

## Install Jenkins

Install Jenkins LTS and a supported Java runtime as a local `systemd` service.
Follow Jenkins's [official Linux installation guide](https://www.jenkins.io/doc/book/installing/linux/)
for the workstation distribution.

## Open Jenkins browser

Open `http://localhost:8080` on the workstation, unless its Jenkins service was
configured with a different host or port. Complete first-start setup and create
an administrator account. During initial setup, install Pipeline, Git,
Credentials Binding, Git client, GitHub Branch Source, and **Multibranch
Pipeline Inline Definition**. The inline-definition plugin must support the
installed Jenkins version. Restrict job configuration and credential-management
permissions to CI administrators. The initial unlock password is available on
the workstation with:

```sh
sudo cat /var/lib/jenkins/secrets/initialAdminPassword
```

## Set up Jenkins

### Create a Jenkins API token for the launcher

Create a Jenkins API token for the CLI user and store it in a mode-600 file:

```sh
mkdir -p "$HOME/.config/gkeyll/jenkins"; umask 077
printf '%s:%s\n' '<jenkins-user>' '<jenkins-api-token>' > "$HOME/.config/gkeyll/jenkins/team-workstation.auth"
chmod 600 "$HOME/.config/gkeyll/jenkins/team-workstation.auth"
```

### Create the GitHub credential

Create a short-lived GitHub token limited to commit-status access and add it as
a Jenkins **Username with password** credential. The workflow uses it only to
read PR metadata and publish status; source discovery and Git checkouts are
anonymous. Revoke or replace it when workstation ownership changes.

### Configure the Jenkins node and global environment

Label the isolated agent and configure its non-interactive PATH with the
toolchain, `cmake`, and Python/NumPy. Set these global environment variables:

| Name | Meaning |
| --- | --- |
| `TEAM_WORKSTATION_NODE_LABEL` | Required Jenkins agent label |
| `TEAM_WORKSTATION_MKDEPS_SCRIPT` | Required `machines/` dependency script |
| `TEAM_WORKSTATION_CONFIGURE_SCRIPT` | Required `machines/` configure script |
| `TEAM_WORKSTATION_GITHUB_CREDENTIAL_ID` | GitHub status/API credential ID |
| `TEAM_WORKSTATION_BUILD_JOBS` | Optional; default `3` |
| `TEAM_WORKSTATION_REGRESSION_JOBS` | Optional; default `1` |
| `TEAM_WORKSTATION_MPIEXEC` | Optional parallel-test launcher |
| `TEAM_WORKSTATION_STATUS_CONTEXT` | Optional status context; default team-workstation |
| `TEAM_WORKSTATION_TRUSTED_CI_REF` | Trusted workflow branch/SHA; production value `main` |

### Create the one centralized multibranch Pipeline job

Create Multibranch Pipeline `gkeyll-ci-team-workstation` from GitHub source
`gkeyllorg/gkeyll`. Discover `main` and pull requests, exclude ordinary
branches that are also PRs, and enable periodic scans or GitHub webhooks.
Configure the GitHub source anonymously. Do not configure an SCM Jenkinsfile
path: choose **Common pipeline definition for markerfile**, set the marker to
`ci/jenkins/.team-workstation-marker`, and paste this controller-owned
bootstrap:

```groovy
def nodeLabel = env.TEAM_WORKSTATION_NODE_LABEL?.trim()
if (!nodeLabel) error('TEAM_WORKSTATION_NODE_LABEL is not configured.')

def trustedRef = (env.TEAM_WORKSTATION_TRUSTED_CI_REF ?: 'main').trim()
def isSha = trustedRef ==~ /[0-9a-fA-F]{40}/
if (!isSha && (!(trustedRef ==~ /[A-Za-z0-9][A-Za-z0-9._\/-]*/)
    || trustedRef.contains('..') || trustedRef.contains('//')
    || trustedRef.endsWith('/') || trustedRef.endsWith('.') || trustedRef.endsWith('.lock'))) {
    error("TEAM_WORKSTATION_TRUSTED_CI_REF must be a safe branch or full SHA, got '${trustedRef}'.")
}

def trustedDir = '_trusted_ci'
def trustedRefspec = isSha
    ? "+${trustedRef}:refs/remotes/origin/trusted"
    : "+refs/heads/${trustedRef}:refs/remotes/origin/trusted"
node(nodeLabel) {
    def trustedCommit = ''
    dir(trustedDir) {
        deleteDir()
        checkout([
            $class: 'GitSCM',
            branches: [[name: 'refs/remotes/origin/trusted']],
            doGenerateSubmoduleConfigurations: false,
            extensions: [
                [$class: 'CleanBeforeCheckout'],
                [$class: 'CloneOption', depth: 1, honorRefspec: true, noTags: true, shallow: true]
            ],
            userRemoteConfigs: [[
                refspec: trustedRefspec,
                url: 'https://github.com/gkeyllorg/gkeyll.git'
            ]]
        ])
        trustedCommit = sh(returnStdout: true, script: 'git rev-parse HEAD').trim()
    }
    withEnv([
        'CI_TEAM_WORKSTATION_BOOTSTRAPPED=true',
        "CI_TRUSTED_CI_COMMIT=${trustedCommit}",
        "TRUSTED_CI_DIR=${env.WORKSPACE}/${trustedDir}"
    ]) {
        dir(env.WORKSPACE) {
            load "${env.WORKSPACE}/${trustedDir}/ci/jenkins/jenkinsfile.team_workstation"
        }
    }
}
```

The bootstrap must remain Jenkins configuration. It is the trust boundary: do
not replace it with `ci/jenkins/jenkinsfile.team_workstation` as an SCM script,
and do not add a candidate `load` to it. The workflow archives the trusted CI
commit alongside candidate and baseline commits.

# Launching CI jobs

## CLI launch

PR scans create `PR-<number>` jobs automatically for PRs targeting `main`.
Use the CLI to inspect that work or to start an explicit run through the
trusted `main` child:

```sh
./ci/jenkins/gkeyll-ci.sh team scan
./ci/jenkins/gkeyll-ci.sh team run --pr 1234 --follow
./ci/jenkins/gkeyll-ci.sh team run --candidate-ref feature/new-solver --baseline-ref main --follow
./ci/jenkins/gkeyll-ci.sh team active
./ci/jenkins/gkeyll-ci.sh team info --build 42
./ci/jenkins/gkeyll-ci.sh team artifact --build 42 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh team abort --build 42
```

`--candidate-ref` and `--baseline-ref` each accept a safe branch name or a
full 40-character commit SHA. Supply both together; do not combine either with
`--pr`. The browser offers the same parameters on the `main` child. Set
`JENKINS_CLI_AUTH_FILE` only when using a non-default CLI credential path.

## Browser launch

Open the multibranch job. PR child jobs show automatic builds; use the `main`
child's **Build with Parameters** page for explicit PR or candidate/baseline
comparisons.

# Troubleshooting

## Jenkins API and builds

If a PR child does not appear, run `team scan` and inspect the multibranch
indexing log. Confirm that `ci/jenkins/.team-workstation-marker` exists in the
revision, the GitHub source discovers pull requests, and the PR targets `main`.

For a build failure, inspect `ci-trusted-ci-commit.txt`,
`ci-candidate-commit.txt`, and `ci-baseline-commit.txt` artifacts first. A
missing trusted-CI artifact indicates a bootstrap/checkout problem; an
invalid selector or missing status credential fails before candidate testing.
Review agent capacity and Jenkins queue state when several PRs are active.

GitHub branch protection requiring
`continuous-integration/jenkins/team-workstation` is recommended when this
workstation becomes an active merge gate. It is intentionally not required for
the present infrastructure setup.

Production PR jobs execute only the merged `main` versions of
`jenkinsfile.team_workstation` and `jenkinsfile.personal`. A PR changing either
file is tested by the current production workflow, not its own modified one.

To validate such a change before merge, an administrator creates or temporarily
reconfigures a separate staging multibranch job with the same bootstrap, sets
`TEAM_WORKSTATION_TRUSTED_CI_REF` to the reviewed CI-change branch or full
commit SHA, and sets `TEAM_WORKSTATION_STATUS_CONTEXT` to a distinct staging
context. Run the intended PR/reference builds, then restore the staging job's
trusted ref to `main`. Merge only after that staging run succeeds.

## Numerical regression differences

Expected numerical changes require a reviewed entry in
`ci/jenkins/expected_regression_diffs.txt`; unlisted differences fail CI.
