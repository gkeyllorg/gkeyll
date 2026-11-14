---
name: gyrokinetic-details
description: Apply gyrokinetic solver conventions when working on gyrokinetic-related source code or input files.
---

# Instructions

* Use when working on the gyrokinetic/ code, or with an input file that uses
  the gyrokinetic solver.

# Gyrokinetics

## MPI decomposition

Gyrokinetic regression tests parallelize along the last configuration-space
dimension. Choose the partition flag by dimensionality:

| Dimensionality | Flag for N partitions |
|---|---|
| 1x2v | `-c N` |
| 2x2v | `-d N` |
| 3x2v | `-e N` |

For example, run the compiled 2x2v sheath regression on four MPI ranks:

```sh
<mpi_install_prefix>/bin/mpirun -np 4 ./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -M -d 4
```

Use the MPI installation matching the build. `-M` enables MPI in the
regression executable; `-d 4` partitions its second configuration dimension.
Add `-s1` for a one-step smoke test.

See [compiling](../compiling/SKILL.md) when a build is needed,
and [testing-and-verification](../testing-and-verification/SKILL.md) for general test
execution and memory checks.
