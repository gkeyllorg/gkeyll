---
name: compiling
description: How to compile Gkeyll libraries, executables, and specific unit or regression test targets on CPU or CUDA.
---

# Instructions

Use this skill whenever you need to compile
* Any of the Gkeyll libraries.
* A C unit or regression test.
* A C input file.

# Compiling in Gkeyll

We assume the dependencies have already been installed (e.g. via a machines/mkdeps
file) and configured (e.g. via a machines/configure file).

Operate relative to the repository root (`git rev-parse --show-toplevel`).
Dependencies are typically already installed in the sibling `gkylsoft/`
directory. Check the existing configuration before building. If configuration
is missing or needs changing, consult `machines/configure.<machine>.sh` for
required modules and library paths, or use `./configure` with the intended
installation prefix and solver.

## Compiling the Gkeyll library

We compile using a Makefile. The `make` command allows you to specify the number
of cores to use via `-j`. We will use no more than half the available processors for parallel builds.
So before invoking `make`, set the `NPROC` environment variable on a Mac using
```
export NPROC=$(sysctl -n hw.physicalcpu)
NPROC=$((NPROC > 1 ? NPROC / 2 : 1))
```
or on Linux using

```sh
NPROC=$(nproc)
NPROC=$((NPROC > 1 ? NPROC / 2 : 1))
```

Then, you can build the Gkeyll library with

```sh
make -j"$NPROC" install
```

## Compiling a C unit or regression test

Unit tests are located in the `unit/` directory of each solver
(core/moments/vlasov/gyrokinetic/pkpm), while C regression tests are located in `creg/` folders and
have names beginning with `rt_`.

Compile a specific test by giving make its executable target:

```sh
make -j"$NPROC" build/core/unit/ctest_array
make -j"$NPROC" build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1
```

CPU targets are located in `build/`, and  GPU targets are in `cuda-build/`.
Match test targets to the configured build directory.

For executing the resulting tests, see
[testing-and-verification](../testing-and-verification/SKILL.md).

### Compiling all or groups of C unit tests

The Makefile supports compiling groups of unit tests all at once. For example, to compile all the
unit tests at once run

```sh
make -j"$NPROC" unit
```

Alternatively, you can compile all the unit tests for a specific solver. For example, to compile all
the unit tests for the gyrokinetic solver use

```sh
make -j"$NPROC" gyrokinetic-unit
```

## Compiling C input files

Input files written in C work just like C regression tests. They can be compiled similarly if they
are located in the `creg/` folder. For example, if `my_sim.c` is a gyrokinetic input file, and located
in `gyrokinetic/creg/`, it can be compiled with

```sh
make -j"$NPROC" build/gyrokinetic/creg/my_sim
```

If the input file is outside of the repository, use the shared makefile in
`<gkeyll_install_prefix>/share/`. Run it from the directory containing the input
file:

```sh
cd <input_file_directory>
cp <gkeyll_install_prefix>/share/Makefile ./
make
```
