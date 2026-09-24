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
commit statuses. Gkeyll is public, so Pipeline source and candidate checkouts
are anonymous. Each controller stores its own classic GitHub PAT with only the
`repo:status` scope for authenticated GitHub API requests and status
publication; an outside collaborator with push access can use this model
without Gkeyll organization membership. See the platform guide for the
credential owner. Do not place credentials in the repository or in candidate
branches.

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
