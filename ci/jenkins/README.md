# Jenkins CI for Gkeyll

This sets up Jenkins (following the pattern used by
[SUNDIALS](https://github.com/llnl/sundials/tree/main/test/jenkins), adapted
for GitHub instead of Bitbucket) to build Gkeyll and run its unit tests on our
own persistent machines any time a pull request into `main` is opened or
updated.

Like `.github/workflows/mac_build.yml`, every build starts from a clean git
checkout and builds its dependencies (`gkylsoft/`) from scratch into its own
Jenkins workspace via `machines/mkdeps.macos.sh` — but on our own persistent
hardware instead of an ephemeral GitHub-hosted runner, and it actually
executes the unit tests rather than just building them.

Each Jenkins controller runs on whichever machine you install it on —
laptop, desktop, or workstation. `manauref_lt1` (manauref's laptop) was the
first one set up this way. Build nodes — a controller's own machine, plus any
others registered later, personal or shared — are added and configured the
same way (see step 4).

This doc covers both adding a node under an *existing* controller (e.g.
`manauref_lt1`'s), and setting up your *own*, fully independent controller on
your own machine (step 3a) if you'd rather not grant that existing
controller's admin SSH access to your machine.

## 1. Install Jenkins

**macOS:**

```
brew install jenkins-lts
brew services start jenkins-lts
```

This pulls in `openjdk@21` as a dependency. `JENKINS_HOME` is
`~/.jenkins` (Homebrew's default on macOS — not
`/opt/homebrew/var/lib/jenkins`, despite what some older docs say). Get the
initial admin password with:

```
cat ~/.jenkins/secrets/initialAdminPassword
```

**Linux:**

Follow the [official Jenkins Linux install
instructions](https://www.jenkins.io/doc/book/installing/linux/) for your
distro (e.g. the `apt`/`yum` package repo), then start it via `systemctl
start jenkins` (most distros enable/start it automatically on install).
`JENKINS_HOME` is `/var/lib/jenkins` by default. Get the initial admin
password with:

```
sudo cat /var/lib/jenkins/secrets/initialAdminPassword
```

**Both:**

Jenkins listens on `http://localhost:8080`. Open it, paste the
password, and install the "suggested plugins" set when prompted.

## 2. Install additional plugins

Manage Jenkins → Plugins → Available plugins, install:

- **GitHub Branch Source** — lets Jenkins discover branches/PRs on a GitHub
  repo and is the GitHub equivalent of the Bitbucket Branch Source plugin
  SUNDIALS uses.
- **Collapsing Console Sections** (optional) — collapses long build log
  sections for readability, same as SUNDIALS' setup.

## 3. Add a GitHub credential

Jenkins needs read access to `gkeyllorg/gkeyll` to poll for branches and PRs:

1. Create a GitHub Personal Access Token with `repo` scope (Settings →
   Developer settings → Personal access tokens on GitHub).
2. In Jenkins: Manage Jenkins → Credentials → System → Global credentials →
   Add Credentials → kind "GitHub personal access token" (or "Username with
   password", username = your GitHub username, password = the token).

## 3a. Scope your controller to its own nodes

`ci/jenkins/Jenkinsfile` has a single `nodes` map listing every known build
node across every independent controller that builds this repo (yours and
anyone else's). This keeps that map a shared, PR-reviewable source of truth
— but it also means your controller must be told which of those labels it
actually owns agents for, or it will try to start a `node(label)` block for
someone else's node and hang forever waiting for an executor that will never
exist under your controller.

Manage Jenkins → System → Global properties → check "Environment variables"
→ Add:

- Name: `CI_OWNED_NODE_LABELS`
- Value: a comma-separated list of the node label(s) *this controller* owns,
  e.g. `manauref_lt1`, or `alice_laptop,alice_workstation` if you register
  more than one node under your own controller.

This is a controller-wide setting — set it once, regardless of how many
nodes your controller owns. If it's unset (or blank), the pipeline fails
fast with an error telling you to set it, rather than hanging; that's
expected the first time you set up a new controller.

## 4. Add a build node

Every machine that runs builds — including the controller's own machine — is
a Jenkins "node" with a label matching one used in `ci/jenkins/Jenkinsfile`'s
`nodes` map. This procedure is the same whether you're setting up
`manauref_lt1` (the first node) or adding another machine later.

### 4.1. Register the node and label it

- **This machine** (the controller itself, e.g. `manauref_lt1`): Manage
  Jenkins → Nodes → "Built-In Node" → Configure → Labels: add the label.
- **A separate machine** (e.g. a workstation): Manage Jenkins → Nodes → New
  Node → "Permanent Agent" → Launch method "Launch agents via SSH",
  host/credentials for that machine, and give it a label. Nothing to install
  there yourself — Jenkins pushes its agent jar over SSH. Make sure Java, a
  C/C++/Fortran toolchain, and `cmake` are present; `mkdeps` builds everything
  else from scratch per build.

This label is what the Jenkinsfile's `node('<label>')` matches.

### 4.2. Fix the node's `PATH`

On **macOS**, Jenkins runs as a `launchd` service (via `brew services`),
which does **not** source your shell profile (`~/.zshrc`, `~/.zprofile`) —
so it starts with launchd's bare-bones default `PATH`
(`/usr/bin:/bin:/usr/sbin:/sbin`) and can't find Homebrew tools like `cmake`,
even though `which cmake` works fine in your terminal.

On **Linux**, the built-in node's agent process is started by `systemd`
(via the `jenkins` service), which likewise does not source your shell
profile — it uses whatever minimal `PATH` is baked into the unit's
environment, which typically misses anything installed under `/opt`, a
package manager outside the base distro repos, or `$HOME`. An SSH-launched
agent on another machine (either OS) has its own, similar set of environment
quirks to sort out, since `sshd` non-interactive sessions don't source
`~/.bashrc`/`~/.zshrc` either.

Fix this per-node in Jenkins itself (not by editing the launchd plist,
systemd unit, or a shell profile) so it works the same way regardless of how
each node's agent process is started, and each node keeps its own value —
this generalizes cleanly as you add more, across OSes.

Manage Jenkins → Nodes → `<label>` → Configure → Node Properties → check
"Environment variables" → Add:

- Name: `PATH`
- Value: wherever this node's build tools (`cmake`, compilers, MPI, a
  `python` with `numpy`, etc.) actually live, e.g.:
  - macOS (`manauref_lt1`):
    `/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin`
  - Linux: typically
    `/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin`, plus anywhere your
    toolchain/MPI/Python were installed outside the base distro repos (e.g.
    `/opt/<something>/bin`, a Conda/venv `bin`, or `/usr/local/cuda/bin`).

Don't assume this matches another node's value — especially across OSes.

### 4.3. Personal machine or team machine?

- **Personal computer** (like `manauref_lt1`): the multibranch job builds —
  i.e. executes — the code of any PR opened against the repo, on whichever
  node picks it up. On a personal machine you probably don't want to run
  other contributors' PR code, so restrict this node to PRs you authored: in
  `ci/jenkins/Jenkinsfile`, add `allowedPrAuthors: ['<your-github-username>']`
  to this node's entry in the `nodes` map. PR builds from any other author
  are skipped on this node (logged in the console, not failed); builds of
  `main` itself are unaffected, since there's no PR author to check there.
  Note this is a Jenkinsfile-level control, not a GitHub-side protection — it
  only takes effect once Jenkins has already scanned and started that node's
  branch of the pipeline.
- **Shared/team machine**: omit `allowedPrAuthors` from its entry — it will
  build PRs from anyone.

### 4.4. Register the node in the Jenkinsfile

Add (or confirm) an entry in the `nodes` map in `ci/jenkins/Jenkinsfile`,
using whichever `mkdeps`/`configure` scripts under `machines/` match that
node's OS:

```groovy
def nodes = [
    'manauref_lt1': [mkdeps: 'mkdeps.macos.sh', configure: 'configure.macos.sh', allowedPrAuthors: ['manauref']],
    'workstation-node': [mkdeps: 'mkdeps.linux.sh', configure: 'configure.linux.cpu.sh'],
]
```

Commit and push (a normal PR against `main`, same as any other change).
Adding a node's entry to this shared map doesn't by itself make any
controller build on it — only a controller whose `CI_OWNED_NODE_LABELS`
(step 3a) includes this label will. The next PR build will then run on every
node your controller owns, in parallel (subject to each node's
`allowedPrAuthors`, if set).

## 5. Nothing to pre-build

Unlike an earlier version of this setup, the pipeline does **not** assume a
pre-built `gkylsoft/` on any node. Each build passes
`PREFIX=$WORKSPACE/gkylsoft` to that node's `mkdeps`/`configure` scripts, so
dependencies are built from scratch inside that build's own Jenkins
workspace — isolated from any manual checkout you have elsewhere, and safe
for concurrent builds of different branches/PRs on the same node (each gets
its own workspace, so its own `gkylsoft/`).

The tradeoff is build time: expect each build to spend several extra minutes
building SuperLU/LuaJIT from source before it even gets to compiling Gkeyll,
same as `.github/workflows/mac_build.yml` does on GitHub's runners.

You still need a working toolchain, `cmake`, and a `python` on `PATH` with
`numpy` installed (used to process the ADAS atomic data some gyrokinetic
unit tests, e.g. `ctest_dg_iz`/`ctest_dg_recomb`, need at runtime) on each
node for its `mkdeps` script to succeed — if you've built Gkeyll manually on
that machine before, you already have most of these:

- **macOS**: Xcode command line tools (`xcode-select --install`) for the
  C/C++/Fortran toolchain, plus `cmake`/`python` (e.g. via Homebrew).
- **Linux**: a C/C++/Fortran toolchain (e.g. `build-essential` + `gfortran`
  on Debian/Ubuntu, or the equivalent group on your distro), plus `cmake` and
  a `python` with `numpy` from your package manager (or a Conda/venv install,
  as long as it's on the `PATH` set in step 4.2).

## 6. Create the multibranch pipeline job

New Item → name it (e.g. `gkeyll`) → **Multibranch Pipeline**.

- **Branch Sources** → Add source → GitHub:
  - Credentials: the one added in step 3.
  - Repository HTTPS URL: `https://github.com/gkeyllorg/gkeyll`
  - Behaviors: "Discover branches" — set its strategy to **"Exclude
    branches that are also filed as PRs"**, not "All branches". Also add
    "Discover pull requests from origin". Add "Discover pull requests from
    forks" only if external contributors submit PRs from forks.

    Our PRs are opened from branches within the same repo (not forks), so
    with "All branches" the scanner indexes each such branch twice — once as
    a plain branch, once as a PR head — building it twice and posting two
    separate GitHub status checks for the same commit
    (`continuous-integration/jenkins/branch` and `.../pr-head`). Excluding
    branches that are also filed as PRs leaves `main` (and any other branch
    with no open PR) reporting as `.../branch`, while PR branches report
    only as `.../pr-head`.
- **Build Configuration** → Mode: "by Jenkinsfile", Script Path:
  `ci/jenkins/Jenkinsfile`.
- **Scan Multibranch Pipeline Triggers** → check "Periodically if not
  otherwise run" → interval **2 minutes**.

If your controller is behind NAT (as `manauref_lt1`'s is), GitHub can't push
a webhook to it. The periodic scan above is how Jenkins learns about
new/updated PRs instead — every 2 minutes it asks GitHub for changes and
kicks off a build if there are any. This means there's up to a ~2 minute
delay before a build starts, which is an acceptable tradeoff for not having
to run a public tunnel (e.g. ngrok) just to receive webhooks. If your
controller has a reachable public address, you can set up a GitHub webhook
instead for near-instant triggering — not covered here since none of our
controllers currently need it.

Save. Jenkins will scan the repo, find `main` and any open PRs with a
Jenkinsfile, and start building them.

## What the pipeline does

First, the pipeline filters the shared `nodes` map down to just the label(s)
in this controller's own `CI_OWNED_NODE_LABELS` (step 3a) — failing fast if
that's not set. Then, per node this controller owns, in parallel (skipping a
node entirely if the PR's author isn't in that node's `allowedPrAuthors`,
when set):

1. `printenv` — for debugging the build environment.
2. `make clean` and `rm -rf $WORKSPACE/gkylsoft $WORKSPACE/_main_baseline` —
   the workspace persists across builds for a given branch/PR, so without
   this, a stale or (e.g. from a prior interrupted or racing build)
   corrupted file can look up-to-date and never get rebuilt/reinstalled,
   silently breaking the link or mixing old and new dependency files. Wiping
   everything every run makes each build genuinely from scratch.
3. That node's `mkdeps` script (e.g. `machines/mkdeps.macos.sh`) — builds
   `gkylsoft/` from scratch into `$WORKSPACE/gkylsoft`.
4. That node's `configure` script (e.g. `machines/configure.macos.sh`) —
   generates `config.mak` pointing at that freshly-built `gkylsoft/`.
5. `make -j3 check` — builds **and runs** all unit tests (`core`, `moments`,
   `vlasov`, `gyrokinetic`, `pkpm`); a failing unit test fails the build.
6. `make -j3 regression` — builds (does not execute) all regression tests,
   matching today's `.github/workflows/mac_build.yml` scope. This is a cheap
   compile-only smoke check across the full regression-test corpus; it
   doesn't overlap with step 8 below, which actually runs a curated subset.
7. `make -j3 gkeyll-install` — installs the built tree into
   `$WORKSPACE/gkylsoft/gkeyll`, which the regression-test runner (step 8)
   needs: `gkeyll` is both the interpreter it runs as, and its C tests are
   compiled on the fly using the installed `share/Makefile`.
8. **PR builds only** — actually runs the MOAT ("Mother Of All Tests")
   curated regression-test subset (34 tests across `moments`/`vlasov`/
   `gyrokinetic`/`pkpm`) via `gkeyll runregression`, and fails the build on
   any unacknowledged failure or diff:
   1. Builds the PR's target branch (`main`, normally) from scratch in a
      side workspace (`$WORKSPACE/_main_baseline`) and runs
      `gkeyll runregression run --moat create` there, producing a baseline
      generated on this same node/build — never a stale or
      cross-machine-generated one (the default comparison tolerance is
      `1e-12`, tight enough that results aren't guaranteed to reproduce
      bit-for-bit across different machines/toolchains/BLAS builds).
   2. Configures `runregression` for the PR's own install, moves that
      baseline's accepted outputs into the PR's own `gkeyll-results/` tree,
      and runs `gkeyll runregression run --moat check` there, executing the
      34 tests and diffing their output against it.
   3. Runs `ci/jenkins/check_regression_results.lua` against the resulting
      per-layer `regressiondb` files, which fails the build if any test
      didn't pass **unless** it's listed in
      `ci/jenkins/expected_regression_diffs.txt` — see below. `main`-branch
      builds skip this whole stage (there's no target branch to diff
      against) but do a cheap check that this file has been cleared out.

   Building the PR's target branch from scratch roughly doubles a PR
   build's total time; that's the tradeoff for a same-session baseline
   instead of a committed or stale one.
9. Archives `build/**/*.log` and, for the regression-test stage,
   `gkylsoft/gkeyll-results/**/*.txt`/`regressiondb`, so logs and results are
   downloadable from the Jenkins build page even on failure.

### Acknowledging an expected regression diff

If a PR intentionally changes a MOAT test's output (e.g. a physics or
algorithm change), the build's `check_regression_results.lua` step still
fails until you tell it that's expected: add a
`<layer>/<test-name>` line to `ci/jenkins/expected_regression_diffs.txt`
(the failing build's console output prints the exact line to add) with a
short reason, then push again. **Clear those lines back out once the PR is
merged** — an entry left behind could silently hide a real future regression
on that same test the next time it fails for an unrelated reason. Leaving
the file non-empty on `main` prints a (non-blocking) warning as a reminder.
