# Gkeyll Jenkins CI

Gkeyll has independent Jenkins setups for machines that provide persistent or
specialized CI capacity. Each setup owns its Jenkins job, node configuration,
credentials, and command-line client.

- [Personal computer CI](README.personal.md): manually run a selected GitHub
  pull request or candidate/baseline branch or commit on a local computer.
- [Stellar CPU CI](README.stellar_cpu.md): manual CPU-only Slurm CI at
  Princeton Stellar.
- [Perlmutter GPU CI](README.perlmutter_gpu.md): manual GPU Slurm CI at NERSC.
- [Team workstation CI](README.team_workstation.md): periodically discover
  pull requests into `main` from every author, with optional selected runs.

The Jenkinsfiles and clients in this directory publish machine-specific GitHub
commit statuses and a per-machine CI report comment (see below). Gkeyll is
public, so Pipeline source and candidate checkouts are anonymous. Each
controller stores its own classic GitHub PAT with only the `public_repo` scope
for authenticated GitHub API requests, status publication, and report
comments; an outside collaborator with push access can use this model without
Gkeyll organization membership. See the platform guide for the
credential owner. Do not place credentials in the repository or in candidate
branches.

## Results in GitHub

Every completed or failed run attempts to post a Markdown report to GitHub.
PR runs update a pull-request comment; selected branch/commit runs update a
comment on the candidate commit. Reports include the failed stage and command
output, unit and regression results, timings, and collapsible sections for:

- **New warnings vs main** (or the explicitly selected baseline).
- **All warnings** and **all errors**, with source log names, line numbers, and
  nearby diagnostic context.
- Failed command output and an inventory of captured logs.

The reporter reads build/configuration logs, Slurm output, and the individual
compiler/runtime logs stored in regression databases. Warning comparison uses
matching completed baseline steps from the same run. It ignores checkout paths,
ANSI colors, and source line/column changes. Missing or failed baseline steps are
reported as unavailable, never treated as a clean baseline. Warnings from shared
or unmatched logs remain visible in **All warnings**. A selected baseline other
than `main` is labeled explicitly. HPC runs also compile baseline unit tests so
candidate unit-build warnings have a corresponding baseline.

Large diagnostic sections continue in additional comments rather than being
truncated. Hidden markers keyed on the machine's status context allow later runs
to update those comments and clear obsolete continuation pages. The archived
`ci-report.md` contains the entire report; `ci-report.md.json` contains the exact
comment pages. Raw logs remain available through the artifact command below.

Reporting runs during failure cleanup, including checkout failures when a PR or
candidate ref is known. Timing-summary errors do not skip publication. Publishing
retries three times; GitHub/network outages, missing credentials, or an unavailable
Jenkins agent can still prevent delivery. Publication failure leaves the original
build result unchanged and adds "(report not posted)" to the status description
when status publication is possible.

The Pipeline runs `github_report.py` only from reviewed code (the trusted
team-workstation checkout, otherwise `main`, or `GKEYLL_CI_TRUSTED_REF` when
staging a CI change), never from the candidate checkout. Diagnostic collection
runs before binding the GitHub token. Credential-bearing shell steps are not
captured in published logs.

Run offline reporter tests with:

```sh
python3 -m unittest discover -s ci/jenkins -p 'test_*.py'
```

## Unified local command

The script `gkeyll-ci.sh` provides a CLI to run, query and terminate CI. See

```sh
./ci/jenkins/gkeyll-ci.sh -h
./ci/jenkins/gkeyll-ci.sh --help
```

It can be used with any of the machines listed above, for example:

```sh
./ci/jenkins/gkeyll-ci.sh personal run --pr 1128 --follow
./ci/jenkins/gkeyll-ci.sh stellar_cpu recent
./ci/jenkins/gkeyll-ci.sh perlmutter_gpu status --queue 42
./ci/jenkins/gkeyll-ci.sh stellar_cpu info --build 6
./ci/jenkins/gkeyll-ci.sh stellar_cpu artifact --build 6 --fetch --only ci-regression-summary.txt
./ci/jenkins/gkeyll-ci.sh team scan
```

The wrapper runs the chosen platform client locally and passes through all
arguments and environment variables. Use `./ci/jenkins/gkeyll-ci.sh --help`
to list platforms, or `./ci/jenkins/gkeyll-ci.sh PLATFORM --help` for the
selected client's commands and setup requirements.

### Build details and artifacts

`info --build NUMBER` shows retained build metadata, its archived failure
summary when present, queryable regression failures, and every retained
artifact. `artifact --build NUMBER` lists those artifacts without opening a
browser. Add `--fetch` to download all artifacts into a new
`gkeyll-ci-build-NUMBER` directory, or use `--only PATH[,PATH...]` and
`--output-dir DIR` to select artifacts and a new destination. Downloaded
artifacts retain their Jenkins-relative directory hierarchy.
