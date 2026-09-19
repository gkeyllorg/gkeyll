# Bazzite shared GPU workstation

This provisions the existing [team workstation CI](README.team_workstation.md)
on the Bazzite workstation with an NVIDIA RTX 3090. Jenkins discovers PRs from
all authors every two minutes and tests those targeting `main`. The controller
and CUDA agent are rootless Podman services. Tailscale Serve supplies HTTPS
access to Jenkins and the host's Cockpit terminal.

**Lock the screen or turn off the display when leaving the computer.** Suspend
and hibernate stop Jenkins and GPU work. The host setup disables those power
states; services start at boot and survive logout through user lingering.

## Install

From the repository root, as the desktop user:

```sh
bash ci/jenkins/install-user.sh
sudo bash ci/jenkins/setup-host.sh mrosen
tailscale up --hostname=bazzite-gkeyll
sudo bash ci/jenkins/setup-tailnet.sh
```

The first command builds the container images and restarts the services; run
it when no CI jobs are active. The host commands require local administrator
access and Tailscale sign-in. Existing Jenkins data and credentials survive
image rebuilds. This provisioner uses the desktop/Jenkins account `mrosen`.

Jenkins binds to `127.0.0.1:8080`. Cockpit and its SSH authentication backend
bind to loopback ports 9090 and 22. Remote access uses Tailscale HTTPS ports
443 for Jenkins and 8443 for Cockpit; no router forwarding is needed.
`tailscale serve status` prints the actual URLs.

## Credentials and shared access

Open <http://localhost:8080/> as `mrosen`. The generated administrator password
is in `~/.local/share/gkeyll-jenkins/admin-password`.

Add a **Username with password** credential with ID `gkeyll-github`, the
maintainer's GitHub username, and a classic PAT with only `repo:status` scope
and a short expiration. Its owner needs push access to `gkeyllorg/gkeyll`.
The controller uses it for GitHub API discovery and the Pipeline uses it to
publish statuses. Candidate and baseline Git checkouts are anonymous.

A fine-grained token instead needs access to `gkeyllorg/gkeyll`, Pull requests
read permission, Commit statuses read/write, and any required organization
approval. Outside collaborators should use the classic token, as described in
the [GitHub token limitations](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#fine-grained-personal-access-tokens-limitations).
Reading PRs successfully does not prove that status writes work.
GitHub rejected the original token because its lifetime exceeded the
organization's 366-day limit. Store tokens only in Jenkins.

Jenkins self-registration and anonymous access are disabled. The owner can
create separate Jenkins accounts for colleagues under **Manage Jenkins →
Users**. Authenticated colleagues may read, build, and cancel jobs; only
`mrosen` has administrator access by default. Each CLI user creates their own
Jenkins API token and mode-600 authentication file.

PR discovery covers all contributors independently of dashboard access.
Remote dashboard users also need Tailscale access to this machine. Configure
tailnet sharing/grants for the intended colleagues; keep Cockpit access to
host administrators. An external-device login must be tested from that device.

## Pipeline and node configuration

The provisioned multibranch job is `gkeyll-ci-team-workstation`. Its `main`
child accepts explicit PR or candidate/baseline requests; PR children run
automatically. The previous `gkeyll-ci-personal` job remains available.

The standard **Pipeline: Multibranch with defaults** and **Config File
Provider** plugins supply a controller-managed, sandboxed snapshot of
`jenkinsfile.team_workstation`, `jenkinsfile.personal`, and the regression
checker. The snapshot is assembled from the files copied into the controller
image. Candidate branches cannot replace it, and branches predating Jenkins
need no special files. Rebuild the controller to deploy reviewed CI changes.

The node has one executor. Global `TEAM_WORKSTATION_*` variables select the
existing `mkdeps.linux.sh` and `configure.linux.gpu.sh`, eight build processes,
and one regression process. `MAKEFLAGS=CUDA_ARCH=86` targets the RTX 3090.
Candidate and baseline dependencies are built separately. MPI, NCCL and cuDSS
are not enabled by that configure script. Dependency scripts contain their
own parallelism settings.

The agent includes Java 21, CUDA 12.6, compilers, CMake, rsync, Python/NumPy,
Valgrind, and Compute Sanitizer. Memory checkers run inside this same agent:

```sh
podman exec gkeyll-gpu-agent valgrind --version
podman exec gkeyll-gpu-agent compute-sanitizer --version
podman exec gkeyll-gpu-agent nvidia-smi
```

For CPU memory checks on this Ryzen 7700X, configure a separate CPU build with
`ARCH_FLAGS=-march=x86-64-v3`: Valgrind 3.22 cannot decode the AVX-512
instructions emitted by `-march=native`. When using the image's dynamic
OpenBLAS, set `OPENBLAS_CORETYPE=Haswell` for Valgrind runs. Use the unit-test
option `--no-exec` so memory errors are reported directly by Valgrind. For
example, from the test checkout inside the agent:

```sh
OPENBLAS_CORETYPE=Haswell valgrind --leak-check=full \
  --errors-for-leak-kinds=definite,indirect --error-exitcode=1 \
  ./build/core/unit/ctest_array --no-exec
compute-sanitizer --tool memcheck --leak-check full --error-exitcode 1 \
  ./cuda-build/core/unit/ctest_array --no-exec
```

The existing CI Pipeline does not automatically run these memory checkers.

## Run and inspect CI

```sh
source "$HOME/.config/gkeyll/jenkins.env"
./ci/jenkins/gkeyll-ci.sh team scan
./ci/jenkins/gkeyll-ci.sh team run --pr 1132 --follow
./ci/jenkins/gkeyll-ci.sh team recent
./ci/jenkins/gkeyll-ci.sh team active
./ci/jenkins/gkeyll-ci.sh team abort --build 42
```

These build commands select the `main` child. For an automatic PR build:

```sh
JENKINS_JOB=gkeyll-ci-team-workstation/PR-1132 ./ci/jenkins/gkeyll-ci.sh team recent
JENKINS_JOB=gkeyll-ci-team-workstation/PR-1132 ./ci/jenkins/gkeyll-ci.sh team follow --build 1
```

The environment file references the owner's mode-600 API-token file at
`~/.config/gkeyll/jenkins-cli.auth` and a Java runtime extracted from the
Jenkins image. A remote CLI uses its own Java installation, Jenkins API token,
and HTTPS `JENKINS_URL`.

Full comparisons build both trees and run the C regression suite; they can
take a long time. The one executor serializes work across all PR jobs.

## Verification on 2026-09-17

Both services survived a controller restart and the GPU agent reconnected.
User lingering is enabled, suspend targets are masked, and Tailscale Serve
is configured. The team scan discovered `main` and all 24 open PRs; 21 PRs
targeted `main`. The existing CLI successfully monitored a scan and submitted
PR 1132. Seven isolated CLI checks covered Linux credential permissions,
waiting/cancelled queue entries, and the multibranch scan's HTTP 302 response.

At PR 1132 commit `91ef9ed70d4636ce8bb4da32ed5fa7e4a1eae0b5`, the focused
array/matrix CUDA build passed all 35 cases, including nine device cases.
Compute Sanitizer reported zero errors and zero leaked bytes. The separate
CPU build passed all 26 cases under Valgrind with zero errors and all heap
allocations freed. These diagnostic builds used the agent image's
OpenBLAS/LAPACKE and SuperLU; full Jenkins jobs build their own dependencies.

Shared PR job `PR-1132` build #2 reached the candidate CUDA unit build, then
failed linking `ctest_linsolvers`: its C registrations reference five
`test_cusolver_*_dev` functions, while `ctest_cusolver.cu` defines the names
without `_dev`. The baseline and regression stages were not reached. See
[the handoff](HANDOFF.bazzite.md) for the remaining work and deployment state.

## Services and recovery

```sh
systemctl --user status gkeyll-jenkins gkeyll-gpu-agent
journalctl --user -u gkeyll-jenkins -u gkeyll-gpu-agent -n 60
loginctl show-user mrosen -p Linger
systemctl is-enabled sleep.target suspend.target
tailscale serve status
```

Jenkins state, credentials, and workspaces live in
`~/.local/share/gkeyll-jenkins/`. Quadlets live in
`~/.config/containers/systemd/`. Back up controller state and protected
credential files while Jenkins is stopped. Base images are digest-pinned;
update them deliberately and confirm the agent reconnects after rebuilding.

To restore normal suspend behavior:

```sh
sudo systemctl unmask sleep.target suspend.target hibernate.target hybrid-sleep.target suspend-then-hibernate.target
```

## References

- [Bazzite containers and Quadlet](https://docs.bazzite.gg/Installing_and_Managing_Software/Containers/)
- [NVIDIA CDI with Podman](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/cdi-support.html)
- [Pipeline: Multibranch with defaults](https://plugins.jenkins.io/pipeline-multibranch-defaults/)
- [Tailscale Serve](https://tailscale.com/docs/features/tailscale-serve)
- [Compute Sanitizer](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html)
