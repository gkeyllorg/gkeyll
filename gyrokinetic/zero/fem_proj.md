# FEM mass projection

`gkyl_fem_proj` solves the weak equation

```
integral(v * weight_left * u) = integral(v * weight_right * f)
```

in a continuous nodal FEM space and returns `u` as modal DG coefficients.
Equal weights give a weighted L2 projection. Either weight may independently
be NULL (unity). The consistent mass matrix is used without mass lumping;
weighted products are integrated without truncation into the input DG space.

## Interface

Pass the grid, solve range, serendipity basis, boundary specification, optional
weights, and `use_gpu` through `gkyl_fem_proj_inp` to `gkyl_fem_proj_inew`.
Call `gkyl_fem_proj_set_rhs` before each `gkyl_fem_proj_solve`, then release the
updater with `gkyl_fem_proj_release`. Input/output arrays contain one scalar DG
field, with `basis.num_basis` components and the parent range's array layout.
The output can alias the input; entries outside the solve range are untouched.

The object copies grid/range/basis/boundary metadata and acquires the weight
arrays. Both weights must remain unchanged until release. Recreate the updater
to change either weight. The caller must supply a positive left weight.
The matrix is constructed/factored once. All per-RHS GPU work remains on device.

Boundary types are independent of the Poisson solver:

- `GKYL_FEM_PROJ_NONE`: leave boundary nodes unconstrained (the zero default).
- `GKYL_FEM_PROJ_PERIODIC`: identify nodes; specify both ends of a direction.
- `GKYL_FEM_PROJ_DIRICHLET`: prescribe the constructor's constant value.
- `GKYL_FEM_PROJ_DIRICHLET_VARYING`: evaluate the DG boundary-value array passed
  to `set_rhs` at the physical boundary nodes of skin cells.

Prescribed values must agree at shared nodes, including periodic copies and
intersections of Dirichlet faces. This is a caller requirement. Specialized
source kernels assign prescribed values during cell assembly; all cells sharing
a constrained node assign rather than add to its RHS entry. CUDA uses atomic
assignments for these entries and atomic additions for unconstrained entries.
Dirichlet matrix rows are identity rows; free rows retain their boundary-column
entries. There are no derivative BCs or periodic nullspace/mean corrections:
the positive mass matrix has no constant nullspace.

Optional `bias_plane_list` data prescribes constant values on grid-aligned
planes, including planes inside the solve range. Each plane specifies `dir`,
`loc`, and `val`; its list is copied at construction. Physical coordinates use
the supplied full grid, with solve indices in `[1, grid.cells[d]]`. Planes
outside the solve range are ignored, and alignment is checked within `1e-9`
cell widths. Existing unbiased arbitrary-index ranges remain supported.

Bias constraints must agree at intersections, periodic copies, and shared
Dirichlet nodes. Contradictory constant data is rejected at construction;
compatibility with varying boundary data remains the caller's responsibility.
Bias rows are replaced by unit equations after mass assembly, visiting all
incident cells. Bias values are assigned after RHS assembly, using a separate
CUDA launch on GPU. Only compact incident-face ranges and values are stored;
there is no per-node constraint table or per-advance host transfer. Boundary
biases may instead be supplied through the ordinary Dirichlet BC interface,
as the IWL caller does for an LCFS bias.

## Supported configurations and generation

The initial implementation supports 1D/2D serendipity p1/p2. It has no MPI
collectives: each instance owns its entire solve range. The twist-shift caller
already provides the full perpendicular domain on each parallel partition.

The independent indexing generator is tested through 3D p1/p2, including the
20-node 3D p2 basis. Normal generation and the installed kernels remain 1D/2D;
3D indexing kernels are generated only in temporary test directories. A 3D
constructor call still fails explicitly. Tensor bases are unsupported.

The generator lives in the sibling `gkylcas/gyrokinetic/fem_proj/` directory.
Run from the CAS repository root, with an existing output directory (default
`~/max-out/`):

```sh
maxima --very-quiet -b gyrokinetic/fem_proj/ms-fem-proj.mac
maxima --very-quiet -b gyrokinetic/fem_proj/ms-fem-proj-header.mac
```

Set `fem_proj_out` before loading the drivers to override that directory. Copy generated files
into `gyrokinetic/ker/fem_proj/`; never edit or format them manually. The new
operator uses general nodal utilities and sparse solvers, but does not call
Poisson kernels or the Poisson updater.

## Direct global node numbering

The local-to-global kernels compute compact FEM indices directly, following
Poisson's cell-ownership ordering: the last coordinate varies fastest, and
owned nodes retain their local basis order. Upper-face nodes belong to the next
cell, except at a nonperiodic domain boundary; periodic boundaries wrap the
owner to the first cell. The Maxima implementation is independent of the
Poisson generator and derives row/plane offsets using the same dimensional
algorithm in 1D, 2D and 3D.

The updater computes these indices when assembling the matrix/RHS and when
reconstructing the solution. It retains no per-cell `l2g_map` array and makes
no host transfers for index lookup during GPU solves. Kernel selection is
centralized in `fem_proj_choose_kernels`; cell boundary flags select an entry
in the local-to-global table.

Mass and source assembly follow Poisson's location-based kernel dispatch.
Each direction has interior, lower, upper, and both-boundaries locations;
the fourth location handles a direction with just one cell. The updater selects
location-indexed kernel arrays once, using its fixed boundary types and weights.
There are no per-node boundary records. Ordinary boundary conditions are
handled within assembly; optional internal biases use the separate passes above.
Constant boundary values occupy a fixed-size array copied to CUDA at initialization.

Mass kernels accumulate free rows directly into matrix triples and insert a
unit diagonal at constrained rows. They do not insert off-diagonal zeros,
because periodic local nodes can alias the same global diagonal. Source kernels
integrate the volume term only at free nodes and evaluate prescribed values only
at constrained nodes. At intersecting faces the lowest coordinate direction is
used; the caller's consistency requirement makes this choice immaterial.
Matrix construction remains on the host. A GPU RHS update clears the RHS and
launches one cell-assembly kernel, followed by a bias kernel only when bias
planes intersect the solve range.

## Verification

Build with the configured CPU/CUDA build directory:

```sh
make -j2 LIBRARY_TARGET=library gyrokinetic
make -j2 build/gyrokinetic/unit/ctest_fem_proj \
  build/gyrokinetic/unit/ctest_bc_twistshift \
  build/gyrokinetic/creg/rt_gk_ts_fem_perpproj
./build/gyrokinetic/unit/ctest_fem_proj
./build/gyrokinetic/unit/ctest_bc_twistshift
```

Run the independent generator regression from the CAS repository root:

```sh
python3 gyrokinetic/fem_proj/check-fem-proj-indexing.py
```

It generates and compiles indexing kernels into a temporary directory, then
checks both the new generator and the Poisson generator against independently
enumerated cell owners. The 532 cases cover p1/p2, 1D/2D/3D, every periodicity
combination, all grid shapes with 1–3 cells per direction, and an additional
asymmetric 2x3x4 grid. It checks exact ordering, shared and periodic nodes,
contiguous numbering and bounds. `--work-dir PATH` retains sources and logs.

This test exposed a Poisson 3D preceding-plane offset error. On a nonperiodic
3x2x2 p1 grid, cell (1,0,0) previously produced index 40 despite there being
only 36 global nodes; its first node was 15 instead of 9. The corrected
Poisson generator counts the upper-y/interior-z and upper-y/upper-z portions
without changing the x boundary classification. Its 1D/2D output is unchanged.

`ctest_fem_proj` also checks installed indexing kernels and CPU dispatch on
shifted subranges for all 1D/2D periodicity combinations and cell counts 1–3.
It exercises 5,040 combinations of dimension/order, all supported lower/upper
boundary pairs, independent weights, and grids with one or three cells per
direction. Each uses random DG, continuous, and constant input fields. Another
120 cases check distinct lower/upper values in each direction, including constant
and varying data on opposite faces. Verification uses independent physical-coordinate
node matching and four-point Gauss quadrature to check the weak residual at free
nodes and the prescribed values at constrained nodes; there is no reference solver. Checks include idempotence where
appropriate, weighted integral preservation, face continuity, varying boundary
data, repeated RHS solves, aliased input/output, shifted subrange indices,
untouched ghost entries, unsupported-input rejection, and comparison with the
previous Poisson-based projection. CUDA variants use the same reference.

Run the application regression in isolated output directories, serially and on
four MPI ranks (`-M -e 4`; add `-g` for GPU execution). Compare field coefficients
for default FEM, explicit FEM, and explicit DG settings. Use Valgrind or macOS
leaks on CPU and Compute Sanitizer for CUDA validation when available.

### Local validation of the initial implementation

- CPU library build and all three `ctest_fem_proj` entries passed, including
  all 288 reference configurations and comparison with the Poisson projector.
- All 12 twist-shift unit tests passed, including the projection convergence test.
- Application regression passed serially and on four MPI ranks. Serial/MPI field
  coefficients were identical; replacing Poisson changed the manufactured FEM
  result by at most 1.11e-16. The explicit DG path was unchanged.
- A reduced 8x8x8, 8x4-velocity CBC case was compared against an isolated library
  using the original Poisson projection. Initial potential coefficients differed
  by at most 4.44e-16. After one step, the maximum potential difference was
  5.81e-12 (5.55e-14 relative to the maximum coefficient). Both took the same
  1.1880724244388795e-7 timestep. Ion distribution and density differences were
  also within the scaled 1e-11 tolerance.
- Separate Maxima regeneration produced identical kernel files.
- macOS leaks reported zero leaked bytes for the host reference and Poisson
  comparison tests. Valgrind, CUDA, and Compute Sanitizer were unavailable;
  device builds and tests remain to be verified on a CUDA machine.

### Validation after removing the stored node map

- Both generators passed all 532 independent indexing cases through 3D p2.
- Existing Poisson 1D/2D indexing output remained byte-for-byte unchanged.
- The four CPU `ctest_fem_proj` entries and all 12 twist-shift unit tests passed.
- No 3D FEM projection kernels were added to Gkeyll; mass/source/transform
  kernels are unchanged by the numbering update.
- macOS leaks reported zero leaked bytes for the indexing, reference projection,
  and Poisson-comparison tests.
- CUDA compilation and device validation remain unavailable locally.

### Validation after specializing boundary assembly

- The four CPU `ctest_fem_proj` entries passed, covering 5,160 projection
  configurations, direct-index checks, Poisson comparison, and invalid inputs.
- All 12 twist-shift unit tests passed, including field projection.
- Independent Maxima regeneration reproduced all five installed files exactly.
- macOS leaks reported zero leaked bytes for indexing, projection, and Poisson
  comparison tests. The invalid-input test, which intentionally forks, was excluded.
- CUDA compilation, runtime tests, and Compute Sanitizer remain unverified locally;
  no GPU performance measurement is claimed.

## IWL verification

The partial-range tests add 720 configurations per backend across 1D/2D,
p1/p2, both coordinate directions, independent left/right weights, and one-cell
or three-cell solve extents. Checks cover unconstrained quadrature residuals,
continuity, prescribed values, idempotence for equal weights, and unchanged
cells outside the solve range. Scenarios include internal, outer, crossing,
distinct parallel, outside-range, and periodically identified bias planes.
Invalid-input tests cover misalignment, conflicting constraints, and directions.
