# Reusing installed dependencies

[CI overview](README.md)

All four workflows accept `GKEYLL_CI_PREBUILT_CONFIG`. In **Manage Jenkins →
System → Global properties → Environment variables**, set it to the absolute
path of an existing `config.mak` on the build agent:

```text
GKEYLL_CI_PREBUILT_CONFIG=/opt/gkylsoft/gkeyll/share/config.mak
```

Use the config from a configured source checkout or the copy saved by
`make install` in `<prefix>/gkeyll/share/config.mak`. Keep it and its dependencies
outside disposable Jenkins workspaces, readable by the agent. Container mounts
must preserve absolute dependency paths; Slurm compute nodes must see the run
directory at the same path. Leave the variable empty to use machine scripts.

## Configuration

Supply a complete, self-contained generated config with `PREFIX` defined.
Dependency paths must be absolute or relative to that original prefix, such as
`$(PREFIX)/superlu/lib`; relative includes and source-relative paths are unsupported.
For dependencies installed elsewhere, edit a dedicated CI copy:

| Dependency | Config variables |
| --- | --- |
| BLAS/LAPACK | `LAPACK_INC_DIR`, `LAPACK_LIB_DIR`, `LAPACK_LIB_NAME` |
| SuperLU | `SUPERLU_INC_DIR`, `SUPERLU_LIB_DIR`, `SUPERLU_LIB_NAME` |
| MPI | `CONF_MPI_INC_DIR`, `CONF_MPI_LIB_DIR` |
| LuaJIT | `CONF_LUA_INC_DIR`, `CONF_LUA_LIB_DIR`, `CONF_LUA_LIB` |
| NCCL | `CONF_NCCL_INC_DIR`, `CONF_NCCL_LIB_DIR` |
| cuDSS | `CONF_CUDSS_INC_DIR`, `CONF_CUDSS_LIB_DIR` |
| CUDA math libraries | `CUDAMATH_LIB_DIR` |

The config also supplies compiler, architecture, solver, and feature settings.
Use settings compatible with the CI lanes: regression runs need Lua and MPI;
Perlmutter needs CUDA/NCCL. HPC module-loading steps still run, so dependencies
must match those modules.

For personal CI, set `PERSONAL_MPIEXEC=/opt/mpi/bin/mpiexec` or
`PERSONAL_MPI_HOME=/opt/mpi` to match the config's MPI. For team CI, use
`TEAM_WORKSTATION_MPIEXEC` or `WORKSTATION_MPI_HOME`. HPC retains its Slurm
launchers. Prebuilt mode skips both dependency and configure scripts for
candidate and baseline; the corresponding `PERSONAL_*` or `TEAM_WORKSTATION_*`
script settings are not required.

## Copies and runtime paths

CI copies dependencies, following symlinks, into each run's `dependencies/`
directory. The baseline keeps a persistent copy under
`baseline-cache/<platform>/<baseline-sha>/dependencies/`, reused with its
executable and results; see [cache behavior](README.storage.md). The config and installation are
unchanged. OS libraries under `/usr`, `/lib`, `/lib64`, and `/System`, compilers,
and platform SDKs remain machine-provided.

Each build's config points to its copied libraries and sets `PREFIX` and
`INSTALL_PREFIX` to its own `gkylsoft`. CI uses the copied MPI launcher and
libraries; source `dependencies/env.sh` before manually running a retained
executable. Builds stay at their final paths because binaries embed those paths.
Invalid paths or a source config changed during setup fail the build.

## ADAS data

Gyrokinetic and PKPM configs require ADAS `.npy` tables in the original
`PREFIX/gkeyll/share/adas`; setup fails early if none exist. To prepare the data
for an existing installation, run from the repository root with Python and NumPy:

```sh
set -e
prebuilt_prefix=/opt/gkylsoft  # Use PREFIX from the supplied config.
(
    cd install-deps
    ./mkdeps.sh --prefix="$prebuilt_prefix" --build-adas=yes
)
for element in h he li be b c n o ar; do
    for table in ioniz recomb logT logN; do
        test -s "$prebuilt_prefix/gkeyll/share/adas/${table}_${element}.npy"
    done
done
```

Pass the absolute prefix explicitly; `mkdeps.sh` otherwise defaults to
`$HOME/gkylsoft`. A lone `.npy` file is not a complete installation. Starting
from scratch, build dependencies with the machine script (including ADAS),
then configure Gkeyll with the same prefix and intended options. A prebuilt
Gkeyll executable is not required when using a source checkout's config.

CI snapshots ADAS into `dependencies/data/adas`, then copies it to each build's
own `gkylsoft/gkeyll/share/adas`. It overrides `GKYL_SHARE_DIR` to that build's
share directory, even if the supplied config pins another path. Unit builds
install radiation fits from their own source revision there before tests.
The fits and downloaded tables remain separate from the original installation
and the other build. Fixing data or config after compilation requires a new run
to change paths embedded in the binary.

## Trusted helpers and older revisions

Deploy the trusted Jenkinsfiles and `prebuilt_config.py` together; candidate
and baseline refs may predate prebuilt support. No extra script approvals or
Pipeline Utility Steps plugin are needed. The helper defaults to JSON output;
Jenkins uses `--format env` and reads literal `KEY=value` lines.

Regression tools derive
`<installation-prefix>/gkeyll-results/runregression.config.lua` from the installed
executable. The helper preserves a
literal `PREFIX=` for older runners that cannot parse `override PREFIX :=`.
Updating only the candidate runner cannot fix an older baseline's configuration.
See [CI development tests](tests/README.md#prebuilt-dependencies) for verification.
