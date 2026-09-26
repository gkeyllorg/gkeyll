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

Every Jenkins controller is private, so the commit status alone would only show
red or green. Each run therefore also posts a Markdown report to GitHub, built
by `github_report.py` from the summaries the Pipeline archives: candidate and
baseline commits, the failing stage with the extracted compiler or test error,
the unacknowledged and acknowledged regression tests, and timings. For a PR
run the report is a pull-request comment; for a candidate/baseline run it is a
comment on the candidate commit. The comment carries a hidden marker keyed on
the machine's status context, so each machine owns one comment per PR and
later runs update it in place instead of adding new ones.

Publishing is best effort: if GitHub is unreachable the build result is
unchanged and the status description ends with "(report not posted)". The
Pipeline runs `github_report.py` only from reviewed code (the trusted
team-workstation checkout, otherwise `main`, or `GKEYLL_CI_TRUSTED_REF` when
staging a CI change), never from the candidate checkout, because the GitHub
token is bound while it runs.

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
