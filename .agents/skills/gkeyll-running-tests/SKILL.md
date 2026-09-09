---
name: gkeyll-running-tests
description: Run Gkeyll unit tests and C or Lua regression simulations, including MPI and memory checks. Use when validating changes or reproducing simulation failures.
---

# Running Gkeyll tests

Operate relative to the repository root (`git rev-parse --show-toplevel`).
Choose `unit/` tests for individual modules, `creg/` for C regression
simulations, and `luareg/` for Lua regression simulations.

Build the relevant executable first; see
[gkeyll-compiling](../gkeyll-compiling/SKILL.md) for build commands.
Run compiled tests from the repository root:

```sh
./build/core/unit/ctest_array
./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -s1
```

`-s1` limits the regression to a one-step smoke test. Omit it for the full
simulation; a smoke test alone does not validate long-time behavior. Use the
executable's `-h` output to confirm supported options. CUDA builds use
`cuda-build/` paths and regression tests use `-g` to select GPU execution.

Run Lua inputs using the installed executable (adjust the prefix if needed):

```sh
../gkylsoft/gkeyll/bin/gkeyll <Lua_input_file>
```

## MPI regression runs

Use an MPI-enabled build and preferably the OpenMPI installation used to build
Gkeyll, typically `../gkylsoft/openmpi/bin/mpirun`. Regression executables need
`-M` for MPI. The flags `-c X -d Y -e Z` specify partition counts in the first,
second, and third configuration-space directions; their product must match
the number of MPI ranks.

For gyrokinetic partition restrictions and an example, read
[gkeyll-gyrokinetics](../gkeyll-gyrokinetics/SKILL.md).

## Memory checks

All code must be free of memory errors and leaks. Run the relevant CPU test
under Valgrind before committing:

```sh
valgrind --leak-check=full ./build/core/unit/ctest_array
```

For GPU changes, run the relevant executable under Compute Sanitizer:

```sh
compute-sanitizer --tool memcheck --leak-check full <executable>
```

Supply the executable's GPU options where required. Report unavailable tools
or hardware and distinguish completed checks from checks that could not run.
