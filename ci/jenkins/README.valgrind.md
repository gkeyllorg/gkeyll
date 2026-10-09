# Valgrind memory checking

[CI overview](README.md)

Use `--use-valgrind=yes` with `install-deps/mkdeps.sh` and `configure` to enable
Valgrind-compatible CPU builds. Both accept `yes` or `no` (default: `no`), matching
the other build options. CUDA, NCCL, and cuDSS builds cannot use this profile.
This option prepares code for memory checking; install Valgrind separately and
use `make valcheck` to run it.

## Compiler and dependency profiles

Normal builds retain their configured compiler flags and OpenBLAS CPU detection,
including AVX512 on supported hardware. Valgrind builds use the compiler's target
triple to select a conservative instruction set:

| Compiler target | Gkeyll architecture flags | OpenBLAS target |
| --- | --- | --- |
| x86-64 | `-march=x86-64 -mtune=generic` | `GENERIC`, AVX/AVX2/AVX512 disabled |
| ARM64 | `-march=armv8-a -mtune=generic` | `ARMV8` |

Other compiler targets fail with an explicit diagnostic. This is a baseline
profile, not Skylake: Skylake retains extensions beyond baseline x86-64 even
when AVX is disabled. The profile replaces existing `-m*` compiler options,
preserves optimization and project defines, and adds `-g -gdwarf-4`. It covers
library objects, executables, and SQLite. OpenBLAS runtime CPU dispatch is disabled
in the Valgrind profile so it cannot select newer instructions at runtime.
See the [Valgrind instruction limitations](https://valgrind.org/docs/manual/manual-core.html)
and [OpenBLAS build options](https://www.openmathlib.org/OpenBLAS/docs/build_system/).

Rebuild OpenBLAS in a separate prefix before using the modifier; compiler flags
cannot change an already installed library. From `install-deps`:

```sh
# Default, native dependency build (AVX512 remains available).
./mkdeps.sh --prefix="$HOME/gkylsoft" --build-openblas=yes

# Conservative dependency build, kept separate from normal builds.
./mkdeps.sh --use-valgrind=yes \
  --prefix="$HOME/gkylsoft-valgrind" --build-openblas=yes
```

The selection is saved in `build-opts.sh` and `config.mak`, respectively, so later
builds use the same profile without repeating the option. To return to native
builds, rerun each script with `--use-valgrind=no`. Use fresh Gkeyll build
directories when changing profiles; make does not track compiler flags.
The OpenBLAS script extracts fresh sources to avoid reusing objects from the
other profile. Job counts default to half the available processors; `NPROC`
can select a smaller count.

Copy your CPU prebuilt `config.mak` and point `LAPACK_INC_DIR` and
`LAPACK_LIB_DIR` at the conservative OpenBLAS installation. Retain the remaining
dependency paths and `PREFIX` so MPI, Lua, SuperLU, and ADAS still resolve.
Other dependencies must also use compatible instructions; rebuild any that were
compiled with native CPU extensions. Select the copied config with
`GKEYLL_CI_PREBUILT_CONFIG`; see [prebuilt dependencies](README.dependencies.md).
Do not replace dependencies used by an active job.

## Coverage and local use

CI checks every serial C unit test (`ctest_*.c`) in `BUILD_APP` and its dependency
layers: core → moments → vlasov → gyrokinetic → pkpm. An unset or `all` selection
includes all five layers. MPI-launched tests, Lua tests, and regressions are not
run under Valgrind. Ordinary CI checks still run.

After configuring compatible CPU dependencies:

```sh
# From the repository root; add your usual solver and dependency options.
./configure --use-valgrind=yes --prefix="$HOME/gkylsoft-valgrind"
make BUILD_DIR=build-valgrind valcheck
# A single layer can also be checked:
make BUILD_DIR=build-valgrind core-valcheck
```

`valcheck` builds and checks the configured stack, continues after failing tests
or layers, and returns failure for memory errors, definite/indirect/possible
leaks, assertions, and crashes. Logs are stored beside each executable as
`_val_err`, with process status in `_val_err.exit`. Still-reachable allocations
are reported but do not fail the check by themselves.

## Jenkins

Existing scripts can still supply `GKEYLL_USE_VALGRIND=0|1` in the environment
or as an argument to `mkdeps.sh` and `configure`. Explicit script arguments take
precedence over the environment; the last argument wins. Generated settings
retain that selection even if the environment later changes. A make command-line
assignment, such as `make GKEYLL_USE_VALGRIND=1`, can override `config.mak`.

Set `GKEYLL_USE_VALGRIND=1` in the Jenkins node or global environment and select
the compatible prebuilt config. Personal and team workflows run on the worker;
Stellar runs inside its unit-test Slurm allocation. Increase
`STELLAR_CPU_UNIT_TIME` for the larger, slower suite. Perlmutter GPU is unsupported.
The baseline uses matching compiler flags for comparison but is not memory-checked.
The flag participates in the baseline cache key.

CI writes complete logs under `ci-valgrind/candidate/<layer>/`, with per-test
status in `summary.json`. Layer directories prevent duplicate test names from
colliding. Missing tools, executables, logs, or Valgrind summaries fail the check.
Interrupted logs are included in reports even without a completed summary.
Jenkins archives these artifacts and retains them in the persistent candidate run.

The GitHub report includes every failing test's complete stdout/stderr, paginated
across comments when necessary. Passing logs remain available as artifacts.
Deploy `valgrind.py`, `core/minus/valcheck.sh`, the shared reporter, and Jenkinsfiles
together in the trusted CI revision. The trusted helper supplies compiler flags
and the runner even for candidate/baseline revisions predating this feature.
Historical dependency scripts require prepared prebuilt dependencies.

## Optional container workers

Containers are optional. On clusters and native Ubuntu workers, install or load
Valgrind in the execution environment. On Bazzite or other container hosts, install
it in the worker image. The host operating system need not match the image.

Each machine's deployer manages its worker image, Valgrind installation,
mounts, credentials, and scheduling. This repository does not provide a worker
container build. Switch an existing worker to a new image when it is idle.

On macOS, use a Linux VM/container worker for this workflow, including ARM64
Linux on Apple Silicon. These checks exercise the Linux build, not native macOS
or Apple's Accelerate framework. See [Valgrind supported platforms](https://valgrind.org/info/platforms.html).
No personal paths, local image names, worker names, or fixed UIDs are required.

## Verification

```sh
python3 -m unittest -v ci.jenkins.tests.test_valgrind
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain \
  ci/jenkins/tests/test_valgrind.groovy
```

Fixtures cover default and conservative dependency builds, compiler flags,
configured-layer coverage, duplicate names, GPU rejection, cache separation,
failed builds, crashes, assertions, and complete reports. With Valgrind and a C
compiler installed, they also run programs containing intentional memory errors.
