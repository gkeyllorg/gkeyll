---
name: gkeyll-file-structure
description: Locate Gkeyll solver code, updaters, apps, kernels, and tests. Use when exploring the repository or deciding where implementation belongs.
---

# Gkeyll file structure

Operate relative to the repository root (`git rev-parse --show-toplevel`).

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
* The current PKPM build also links gyrokinetic (see `pkpm/Makefile-pkpm`).
  Keep lower layers independent of higher layers.

Each of the solvers' folders have the sub-folders:
* ker/: C kernels generated with Maxima stored in the gkylcas repository. Do
  not hand-edit generated kernels; modify their source templates instead.
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
