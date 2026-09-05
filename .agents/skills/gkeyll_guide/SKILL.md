---
name: gkeyll_guide
description: Gkeyll solves partial differential equations (PDEs) to model plasmas. Gkeyll has components common to all solvers, and solver-specific components. This skill helps guide the agent more directly to where it needs to go. Use when the user is doing anything with the gkeyll codebase — creating, modifying, running, testing, debugging, or exploring code.
user-invocable: true
---

# Instructions

* Operate relative to the repo root (detect via git rev-parse --show-toplevel).

# Gkeyll guide

## Directory structure

Gkeyll has four PDE solvers: 
* Moments or fluid solver.
* Vlasov solver.
* Gyrokinetic solver.
* PKPM solver.

Correspondingly, these solvers are organized into four separate folders, and they
share some common functionality in a fifth folder. Gkeyll is mostly organized in:
* core/: functionality common to all solvers.
* moments/: files for the moments solver.
* vlasov/: files for the Vlasov solver.
* gyrokinetic/: files for the gyrokinetic solver.
* pkpm/: files for the PKPM solver.

Note that the four solvers are not independent. They have the following
dependencies:
* moments depends on core.
* vlasov depends on moments.
* gyrokinetic depends on vlasov.
* pkpm depends on gyrokinetic.

Each of the solvers' folders have the sub-folders:
* ker/: C kernels generated with Maxima stored in the gkylcas repository. Do
  not modify, and if you think you should, please ask for permission.
* zero/: C and CUDA functions or modules (sometimes we call them updaters),
  some of which call kernels in ker/.
* data/: data needed for some simulations.
* unit/: unit tests of specific components in zero/.
* apps/: apps are called by input files or regression tests, and they organize
  solver workflow or simulations and call modules in zero/.
* creg/: C regression tests or short simulations to ensure solvers work (these
  are also examples of C input files).
* luareg/: Lua regression tests or short simulations to ensure solvers work
  (these are also examples of Lua input files).

## Library dependencies

The Gkeyll source code in gkeyll/ depends on other libraries installed in gkylsoft/.
Most of the time you shouldn't need to read files in gkylsoft/ unless prompted.

## CPU - GPU organization

Gkeyll runs on both CPUs and GPUs, the latter using CUDA.
- CPU-only implementations live in .c files.
- C files call host-side wrappers of CUDA kernels, both of which live in _cu.cu
  files.

## Search hierarchy

Rather than looking into the entire codebase without prior knowledge, look into
folders in the following order:
- Look into a folder (ending in /) the user referenced.
- If the user referenced a specific solver (e.g. gyrokinetic), look into its
  folder.
- Look into folders for solvers that the specified solver depends on, following
  the dependency chain up through core/.

## Naming conventions.

### Files.

Gkeyll is written in C and CUDA, in addition to some Lua input files.
- Public header files have names starting with gkyl_ and ending in .h.
- Private header files have names starting with gkyl_ and ending in _priv.h.
- CUDA files have names ending in _cu.cu. 
- luareg/ folders have Lua input files whose names end with .lua. The Lua
  wrappers are in the apps/ folders and have names ending in _lw.c

### Functions

- Public functions (defined in public header files) should have a name that starts
  with gkyl_.
- Public functions in files in zero/ folders should have a name that starts
  with the name of the file that contains it.
- Private functions (defined in private header files) or static (and not defined
  in private headers) in files in zero/ folders should have a name that starts
  with an abbreviated version of the name of the file that contains it.

### Variables

- Do not use single letter names for variables whose scope spans more than ~15
  lines.

## Module structure and pattern

Most modules in zero/ or apps/ consist of 3 public functions:
1. A creation function, typically called gkyl_<module>_new or gkyl_<module>_init.
2. An execution function, typically called gkyl_<module>_advance or gkyl_<module>_apply.
3. A deletion function, typically called gkyl_<module>_release.

There may also be some additional auxiliary private or public functions.

## Building gkeyll

Typicaly the library dependencies will already be installed, and to build Gkeyll
all you have to do is

make -j3 install

from the gkeyll/ folder.

In order to compile unit regression tests you must give the makefile the target
executable, for example:

make -j3 build/core/unit/ctest_gkyl_array

Similarly, to build a regression test you would do, for example

make -j3 build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1

## Running unit and regression tests

Once the corresponding unit or regression test has been compiled, you can run
them from the gkeyll/ folder like:

./build/core/unit/ctest_gkyl_array
./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1

Regression tests can also run in parallel, preferably using the OpenMPI that
Gkeyll installed in gkylsoft/. For example, to run the rt_gk_sheath_2x2v_p1 on
4 cores you'd use (assuming you are in gkeyll/)

../gkylsoft/openmpi/bin/mpirun -np 4 ./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -M -d 4

Here -M is needed to parallelize, and -c X -d Y -e Z are used to tell the
program how many cores to use in each direction. But note that gyrokinetic tests
only parallelize along the last dimension, so use the flags in the following table:

| Dimensionality | Flag to use |
|---|---|
| 1x2v | `-c X` |
| 2x2v | `-d X` |
| 3x2v | `-e X` |

There are also Lua regression tests you can run using

../gkylsoft/gkeyll/bin/gkeyll <Lua_input_file>
