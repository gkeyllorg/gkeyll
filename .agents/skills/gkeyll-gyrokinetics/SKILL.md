---
name: gkeyll-gyrokinetics
description: Navigate Gkeyll gyrokinetic solver code and select configuration-space MPI partitions for gyrokinetic simulations. Use when working on gyrokinetic code or inputs.
---

# Gkeyll gyrokinetics

Operate relative to the repository root (`git rev-parse --show-toplevel`).
Start in `gyrokinetic/`: `apps/` orchestrates simulations, `zero/` contains
updaters, `creg/` contains C regression inputs, and `luareg/` contains Lua
regression inputs. Follow shared implementation dependencies through
`vlasov/`, `moments/`, and `core/` as needed.

Keep computation and internal state in updaters; apps should delegate through
their API. Do not hand-edit generated `ker/` kernels; their Maxima templates
live in the gkylcas repository.

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
../gkylsoft/openmpi/bin/mpirun -np 4 ./build/gyrokinetic/creg/rt_gk_sheath_2x2v_p1 -M -d 4
```

Use the MPI installation matching the build. `-M` enables MPI in the
regression executable; `-d 4` partitions its second configuration dimension.
Add `-s1` for a one-step smoke test.

See [gkeyll-compiling](../gkeyll-compiling/SKILL.md) when a build is needed,
and [gkeyll-running-tests](../gkeyll-running-tests/SKILL.md) for general test
execution and memory checks.
