---
name: gkeyll-compiling
description: Configure and compile Gkeyll libraries, executables, and specific unit or regression test targets on CPU or CUDA. Use for build tasks and build failures.
---

# Compiling Gkeyll

Operate relative to the repository root (`git rev-parse --show-toplevel`).
Dependencies are typically already installed in the sibling `gkylsoft/`
directory. Check the existing configuration before building. If configuration
is missing or needs changing, consult `machines/configure.<machine>.sh` for
required modules and library paths, or use `./configure` with the intended
installation prefix and solver.

Use half the available processors for parallel builds, rounded down with a
minimum of one job. Set `NPROC` using `nproc` before building:

```sh
NPROC=$(nproc)
NPROC=$((NPROC > 1 ? NPROC / 2 : 1))
```

Build and install with the existing configuration:

```sh
make -j"$NPROC" install
```

Compile a specific test by giving make its executable target:

```sh
make -j"$NPROC" build/core/unit/ctest_array
make -j"$NPROC" build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1
```

CPU builds use `build/`; CUDA builds with `CC=nvcc` use `cuda-build/`.
Match test targets to the configured build directory. Single-layer targets
include `core`, `moments`, `vlasov`, `gyrokinetic`, and `pkpm`; `gkeyll` builds
the Lua executable. Preserve the selected compiler, CUDA architecture, and MPI
configuration when rebuilding.

For executing the resulting tests, see
[gkeyll-running-tests](../gkeyll-running-tests/SKILL.md).
