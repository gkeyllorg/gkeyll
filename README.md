# About

This is the main repository for `Gkeyll`: a multi-scale, multi-physics simulation
framework developed for a variety of applications in plasma physics, space physics,
general relativity, and high-energy astrophysics. The name is pronounced as in the book 
"The Strange Case of Dr. Jekyll and Mr. Hyde" which is required reading for all
members of the `Gkeyll` Team. The core of `Gkeyll` is written in C, featuring a
lightweight, modular design and minimal external dependencies, plus an additional Lua
scripting layer (for specifying simulation input parameters), with trace quantities of C++
for GPU support. `Gkeyll` is developed and maintained at the Princeton Plasma Physics
Laboratory (PPPL) and is copyrighted 2020-2025 by Ammar Hakim and the `Gkeyll` Team.

**NOTE**: Though the framework is called `Gkeyll`, the abbreviation `gkyl` is used
commonly within the code itself.

Documentation is available at http://gkeyll.rtfd.io.

# Code structure

`Gkeyll` is designed around a fully modular architecture, with the entire code being
built hierarchically, with the following layers:

* `core`: essential infrastructure (grid-generation, parallelism, file input/output, Lua
  tools, array operations, fundamental discontinuous Galerkin operations, etc.)
* `moments`: finite-volume solvers for hyperbolic PDE systems, including (relativistic)
  multi-fluid equations, Maxwell's equations, Einstein's equations, etc., plus
  infrastructure for specifying curved spacetime geometries (e.g. for black holes and 
  neutron stars)
* `vlasov`: modal discontinous Galerkin solvers for Vlasov-Maxwell and Vlasov-Poisson
  systems, as well as general Hamiltonian systems via their canonical Poisson bracket
  structure
* `gyrokinetic`: modal discontinuous Galerkin solvers for the full-$f$, long-wavelength
  gyrokinetic equation coupled to the gyrokinetic Poisson equation, plus infrastructure 
  for specifying tokamak and mirror geometries
* `pkpm`: modal discontinuous Galerkin solvers for the
  Parallel-Kinetic-Perpendicular-Moment system

Although `Gkeyll` as a whole comprises all of these layers, each layer can be built
individually, with dependencies only upon the layer below it (e.g. `moments` depends
only on `core`, `vlasov` depends only upon `moments`, etc.).

# Installation

Installing `Gkeyll` consists of three steps: building dependencies, configuring, and
compiling. How you do the first two of these steps varies depending upon whether you are
installing on a supercomputer or operating system for which we already have installation
scripts (i.e. "machine files"), or you are installing on an entirely new supercomputer or
operating system. Since the former is the most common, we describe that first.

## On a known operating system or supercomputer, using machine files

We provide a set of "machine files" to ease the build process. These are stored in the
`/machines` directory. For example, to build on macOS, you simply run:

```
./machines/mkdeps.macos.sh
./machines/configure.macos.sh
```

whereas to build on the Traverse supercomputer, you simply run:

```
./machines/mkdeps.traverse.sh
./machines/configure.traverse.sh
```

At this point, if you are building on a supercomputer/cluster (i.e. not on a personal
computer), you may need to load certain environment modules manually. The requisite
environment modules can be found at the top of the corresponding
`/machines/configure.<machine name>.sh` file for that computer. For example, for
Stellar-AMD, we see from the top of the `/machines/configure.stellar-amd.sh` configure
script that we must run:

```
module load cudatoolkit/12.4
module load openmpi/cuda-11.1/gcc/4.1.1
```

before proceeding. Next, simply compile the code with:

```
make install -j #
```

where # is the number of cores you wish to use for compilation (if working on the login
node of a cluster, it may be diplomatic to use fewer than the total number of cores
available, in the interests of preserving friendly relations with the cluster admins and
other users).

All of the `/machines` files default to installing dependencies in the `$HOME/gkylsoft/`
directory, and to searching for dependencies there too. If you wish to install somewhere
else, you'll need to indicate the desired installation directory using the `--prefix=`
flag in both the `/machines/mkdeps.<machine name>.sh` and
`/machines/configure.<machine name>.sh` machine files, and also ensure that the paths to
each dependency are correct in the latter file.

### Installing a single solver

The machine files in the repository build all solvers (e.g. moments, vlasov, gyrokinetic),
but if you only want to build one of them, say vlasov, use the additional `--app` flag
in the configure command:

```
./configure --app=vlasov
```

## On a new computer (no machine files available)

When installing on a new operating system or cluster that we don't yet have machine files
for, you could choose to build machine files using existing ones as guides, or you could do
each step manually. We expand on each of these steps below.

### Specifying/building dependencies

`Gkeyll` has a few dependencies (e.g. OpenBLAS and SuperLU, optionally MPI, CUDA, and
NCCL, etc.). If you are working on a cluster, some or all of these libraries are likely
already installed. To install dependencies (either all of them, or only the ones your
cluster doesn't have), use the `mkdeps.sh` script in the `/install-deps` directory, as
follows:

```
cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes
```

In this example we opted to build OpenBLAS and SuperLU, assuming that neither is available
on this new computer. This will install OpenBLAS and SuperLU in the ```$HOME/gkylsoft```
directory; in order to install in a different directory one must specify it with the
`--prefix=` flag.

Note that OpenBLAS **requires you to have gfortran** installed. Likewise, SuperLU
**requires you to have cmake** installed. **You** are responsible for installing this on
your machine. Also keep in mind that, on macOS, you should not need to do anything special
to use LAPACK/BLAS, as these come pre-packaged as part of the Apple developer tools.

At this stage, you can also download ADAS data for neutral interactions (useful for some
gyrokinetic simulations) by adding the flag ```--build-adas=yes``` after the other
dependencies to be installed.

### Configuring

The (optional) configuration step is used when we need to specify the use of specific
compilers, dependency paths, or other options. For example, if a dependency is already
available in your computer/cluster, you should use the configure script to specify the
location of the include, and the full path to the static libraries, for these
dependencies. Run:

```
./configure --help
```

for more details. You can also see machine files in ```machines/``` for examples of how
this script is used.

### Compiling

Once you are done with installing/specifying the compiler and dependencies, simply type:

```
 make -j #
```

in the top-level directory (# is the number of cores; see previous comment on this).

If you do not run configure (see previous section/step) you can also specify the compiler
when running make. For example:

```
 make CC=icc -j #
```

to build with the Intel C compiler (# is the number of cores; see previous comment on
this). If you use the NVIDIA `nvcc` compiler, then the CUDA specific parts of the code
will be built:

```
make CC=nvcc -j #
```

(# is the number of cores; see previous comment on this).

On the other hand, if you do not wish to compile the entire `Gkeyll` system, but only
one specific layer (e.g. `moments`), use the `--app` flag in the configure step and
in the compilation step just type:

```
make moments -j #
```

or its equivalent. To compile the `Gkeyll` executable itself (e.g. to run Lua simulation
files), simply type:

```
make gkeyll -j #
```

If in doubt, just run:

```
make everything -j #
```

to build every part of the `Gkeyll` system (including `Gkeyll` executable, and all unit
and regression tests).

Whatever you choose, in order to test that your installation worked, you can compile one
of the unit or regression tests and run it. See instructions for that below.

# Running a `Gkeyll` input file

There are Lua or C input files for `Gkeyll`. Whenever possible, if the feature you need is
supported, we recommend that you run `Gkeyll` simulations using a Lua input file. To do so
simply pass the (Lua) input file name to the `gkeyll` command installed in
`gkylsoft/gkeyll/bin/`:

```
gkylsoft/gkeyll/bin/gkeyll vlasov_two_stream_p1.lua
```

Some features may only be accessible via a C input file, which needs to be compiled. For
example, if you have the input file `gk_d3d_3x2v_p1.c`:

- Copy `gkylsoft/gkeyll/Makefile` to the directory where `gk_d3d_3x2v_p1.c` is.
- Run `make`, this will compile the input file and create the executable
  `gk_d3d_3x2v_p1`.
- Run the executable. How you run it depends on the cluster and desired runtime
  options, see `gk_d3d_3x2v_p1 -h` and sample SLURM scripts in `machines` starting with
  `jobScript`.

If you are without GPUs (on CPU only), you can alternatively combine the compile and run
the C input file with a single command as we do with Lua:

```
gkylsoft/gkeyll/bin/cgkeyll gk_d3d_3x2v_p1.c
```

# Developing for `Gkeyll`

Development philosophy
---------------------

Out goal is to keep `Gkeyll` as simple and dependency-free as possible. Some dependencies
are unavoidable, like MPI and linear algebra libraries. However, we must avoid an
exponentially increasing dependency chain. Another goal is that `Gkeyll` itself should be
written in pure, modern (C99/C11) C. Some parts of the code need C++ (since `nvcc` is a C++
compiler for CUDA), but the core code itself should be in C.

Developing in C (and C++) requires very strong focus and discipline. **Please consult**
https://en.cppreference.com/w/ for standards documentation for these languages and their
libraries. **Please use valgrind** to make sure all code is "valgrind clean" and
"compute-sanitizer clean" (see section on memory errors below). Pay attention to all
compiler warnings.

Most importantly, **please internalize and follow** the programming philosophy outlined in
the document
["A Minimalist Approach to Software"](https://www.ammar-hakim.org/sj/pn/pn0/pn0-minimalism.html).

When opening issues and pull-requests, you can consider using a template as those in:
[Suggested templates](https://gkeyll.readthedocs.io/en/latest/dev/suggested-templates.html).

Built-in tests
--------------

`Gkeyll` has built-in unit tests (for testing individual units of functionality) and
regression tests (for testing full-scale simulations, either from C or from Lua) that are
run automatically by our CI system, but can also be run manually by developers. These
tests are not compiled by the default `make` command used to compile the code (see
previous sections), but can be compiled by building specific targets. For example, to
build the `ctest_array` unit test in the `/core` directory, run:

```
make build/core/unit/ctest_array
```

or, to build the same unit test on an NVIDIA GPU-accelerated node:

```
make cuda-build/core/unit/ctest_array
```

These produce an executable in the `/build` or `/cuda-build` directories,
respectively, which can then be run directly. For example, to run the `ctest_array`
executable we just compiled, use:

```
./build/core/unit/ctest_array
```

or

```
./cuda-build/core/unit/ctest_array
```

Alternatively you could choose to build all unit tests with:

```
make unit -j #
```

or all unit tests for a specific layer (e.g. `moments`) with:

```
make moments-unit -j #
```

(# is the number of cores; see previous comment on this). A similar procedure
should be followed to compile and run C regression tests in `/<layer>/creg`, e.g. to
build the `rt_10m_gem` C regression test in the `/moments` directory, run:

```
make build/moments/creg/rt_10m_gem
```

or

```
make cuda-build/moments/creg/rt_10m_gem
```

followed by:

```
./build/moments/creg/rt_10m_gem
```

or

```
./cuda-build/moments/creg/rt_10m_gem
```

Regression test executables have a number of run-time options which you can list with the
`-h` flag, i.e. ```<executable> -h```. Likewise, you can also compile all C regression
tests with:

```
make regression -j #
```

or all C regression tests for a specific layer (e.g. `moments`) with:

```
make moments-regression -j #
```

You may also compile *and run* all the unit tests with a single command:

```
make check -j #
```

or, for a specific layer (e.g. `moments`):

```
make moments-check -j #
```

Checking for memory errors
--------------------------

We strive for `Gkeyll` to have no memory errors, even if they are silent or benign. To
ensure this is the case, please check that new solvers or simulations are valgrind clean
on CPUs with e.g.

```
valgrind --leak-check=full <executable> <command_line_arguments>
```

and compute-sanitizer clean on one GPU with e.g.

```
compute-sanitizer --tool memcheck --leak-check full <executable> <command_line_arguments>
```

To check using multiple GPUs on a Perlmuttter node with SLURM e.g.

```
srun -N 1 --ntasks=2 --gpus-per-task=1 --gpu-bind=closest compute-sanitizer --tool memcheck --leak-check full <executable> <command_line_arguments>
```

where the command line arguments must at least contain `-g -M -direction 2` with `direction` being the direction along which the domain is decomposed, e.g. `-e 2` for 3x2v gk simulations.

## Code style

C/C++/CUDA source is formatted by [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
using the style defined in the root `.clang-format` (a Linux-kernel-flavored style with
2-space indentation). Two things are excluded and must never be hand-formatted: files
under any `*/ker/*` directory (auto-generated DG kernel code) and anything under
`core/minus/` (third-party libraries).

One-time setup:

```
pip install pre-commit
pre-commit install
```

After that, `git commit` automatically reformats any staged C/C++/CUDA file that
isn't already formatted (excluding `ker/` and `core/minus/`). CI
(`.github/workflows/format-check.yml`) re-checks the same thing on push/PR as a backstop.

clang-format's output differs slightly across versions, so `pre-commit` and CI both
pin `clang-format` `18.1.8` (see the `rev:` in `.pre-commit-config.yaml`). To run
`ci/format-all.sh` or configure your editor's format-on-save, install the same
version so you don't fight the pinned one:

```
pip install clang-format==18.1.8
```

To format everything by hand (e.g. after pulling changes), run `ci/format-all.sh`,
or `ci/format-all.sh --check` to only check without modifying files.

clang-format expands every element of a braced initializer list onto its own
line whenever the list already ends with a trailing comma, regardless of
column limit -- so the same struct/array literal formats differently
depending on whether its last author happened to add one. `ci/strip-trailing-commas.py`
removes any comma directly before a closing `}` before clang-format ever runs,
so don't rely on a trailing comma to keep an initializer list multi-line; that's
clang-format's call once the comma's gone. This runs automatically as part of
the same pre-commit hook and CI check as clang-format (and as part of
`ci/format-all.sh`), so you shouldn't need to run it by hand.

Since the repo's history includes a single large reformatting commit, run:

```
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

so `git blame` attributes lines to their original author instead of that commit.
