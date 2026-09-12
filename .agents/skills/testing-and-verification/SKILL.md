---
name: testing-and-verification
description: How to run Gkeyll unit and regression tests, as well as performing memory checks. Use when testing changes and verifying that changes didn't break the code or change simulation results.
---

# Instructions

* For code changes, run relevant unit tests and focused regression tests when
  simulation behavior is affected. Choose tests based on the changed component.
* Do not run the entire regression suite unless you are asked to; it can take a long time.
  Use the task's scope and available resources to choose the suitable testing scope.
* For documentation or agent configuration changes, validate the changed files
  and references; solver tests and memory checks are unnecessary.
* Report what was checked and any checks that could not run. If required context
  is missing or unreadable, report the exact path rather than silently skipping it.

# Running unit tests

Operate relative to the repository root (`git rev-parse --show-toplevel`).

In order to run all Gkeyll unit tests:
1. Compile all unit tests; see [compilation](../compiling/SKILL.md).
2. Run all unit tests with
```sh
make -j"$NPROC" unit-run
```

In order to run all unit tests for a specific solver, e.g. gyrokinetic:
1. Compile all unit tests for that solver; see [compilation](../compiling/SKILL.md).
2. Run all unit tests for that solver, e.g. for gyrokinetic use
```sh
make -j"$NPROC" gyrokinetic-unit-run
```

One can also compile and run a single unit test of interest. Simply:
1. Compile the unit; see [compilation](../compiling/SKILL.md).
2. Run the executable from the repository root:
For example, in order to run `core/unit/ctest_array.c` do

```sh
make -j"$NPROC" ./build/core/unit/ctest_array
./build/core/unit/ctest_array
```

If the Gkeyll was build on a CPU-only machine or without GPU support, unit tests will simply ignore
GPU tests (via pre-processor if-statements). When running on a machine with a supported GPU, and if
Gkeyll is configured to build GPU code, running these tests will automatically run the corresponding
GPU tests. If a GPU is not available but Gkeyll was built to use GPUs, the GPU tests may simply fail
or exit (e.g. due to GPU memory allocation errors).

### Parallel unit tests

There are a few parallel unit tests that use MPI, whose names begin with `mctest`. These tests are
compiled the same way as serial tests, but must be run with the appropriate MPI execution command.

For example if Gkeyll was built and configured with the MPI in `gkylsoft/openmpi/`, then a parallel
unit test (for example `mctest_mpi_comm.c`) may be run with

```sh
<mpi_install_prefix>/bin/mpirun -np 4 build/core/unit/mctest_mpi_comm
```

# Running regression tests

## C regression tests

Regression tests written in C and located in `creg/` directories need to be compiled first,
see [compilation](../compiling/SKILL.md).

Once compiled the regression test may be run serially with, for
`gyrokinetic/creg/rt_gk_sheath_2x2v_p1.c` for example:
```sh
./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1
```

See 
```sh
./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -h
```
for command line arguments that may be available (not all are actually supported). For
example, the following command
```sh
./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -s1
```
limits the regression to a one-step smoke test.

CUDA builds place the executable in `cuda-build/` instead of `build/`, and regression tests use `-g`
to indicate GPU execution.

## Parallel regression tests

Regression executables take the `-M` command line argument to indicate a parallel run (requiring
MPI). The flags `-c X -d Y -e Z` specify that the first, second and third configuration-space
dimensions are to be subdivided amongs X, Y and Z MPI processes, respecitively. The product `X*Y*Z`
must match the number of cores available for this run. If one of `-c`, `-d` or `-e` is not given,
it is assumed to be 1.

For example, we may run `gyrokinetic/creg/rt_gk_sheath_2x2v_p1` using 2 cores in the second
dimension with
```sh
<mpi_install_prefix>/bin/mpirun -np 2 ./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -M -d 2
```

For gyrokinetic partition restrictions and an example, read
[gyrokinetic-details](../gyrokinetic-details/SKILL.md).

It is also possible to run with multiple GPUs. Gkeyll's model is to match each MPI process to a
single GPU. The procedure is similar as for multiple GPUs, but the additional `-g` flag is needed.
For example, to run `gyrokinetic/creg/rt_gk_sheath_2x2v_p1` using 2 GPUs in the second dimension
with
```sh
<mpi_install_prefix>/bin/mpirun -np 2 ./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -g -M -d 2
```

## Lua regression tests

Run Lua inputs using the installed executable:

```sh
<gkeyll_install_prefix>/bin/gkeyll <Lua_input_file>
```

## The runregression tool

After installing the Gkeyll library and executable one may use the `runregression` tool to run
regression tests, see

```sh
<gkeyll_install_prefix>/bin/gkeyll runregression -h
```

After configuring runregression, and having generated accepted results with
`runregression run create`, one may check the test using `runregression run check`. Note that the
flag `-r` allows check a specific test or set of tests.

# Memory checks

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
