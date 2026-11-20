# Twist-shift low-pass filter breaks the zonal (ky = 0) consistency on closed field lines

Status: findings to start an investigation on `main`. Nothing here is fixed yet.
Scope: `main` only (gyrokinetic, electrostatic). Data and scripts are listed at the end.

## Summary

Since PR #1088 (`gk_tsbc_filter_fix`, merged into `main` on 2026-10-01 as `e5134d8d1`), every
twist-shift (TS) boundary condition goes through a supersample, shift, low-pass filter and restrict
pipeline by default. The filter acts along the shear direction x. It is applied to the species
distribution functions and to the parallel smoothing of phi. It is on unless all three
`core_parallel_bcs.ts_*` parameters are set to -1.

An exact twist-shift is a shift in y only, so it must leave the y-average (ky = 0 component) of any
field unchanged. Filtering in x also changes the ky = 0 component: the y-averaged profile that comes
back through the TS boundary is a radially smoothed copy of the one that left. Observed consequences:

1. **Fields.** Inside the LCFS, the y-averaged phi is no longer continuous across z = +-pi. The jump
   is 1-16 % instead of ~1e-15, and it is largest at the radial edges of the closed-field-line region.
2. **Distribution functions.** f is smoothed in x every time it crosses z = +-pi, which acts as a
   numerical radial diffusion of the core profiles on closed field lines. A rough estimate for
   electrons is D ~ dx^2 * f_transit ~ (5 mm)^2 * 3.5e5 1/s ~ 9 m^2/s; ions are slower by ~sqrt(m_i/m_e).
3. **Physics.** In a TCV inner-wall-limited (IWL) electrostatic run:
   - the core Er well is about half as deep;
   - phi inside the LCFS turns positive;
   - n_e drops and T_e is smeared outward across the LCFS;
   - fluctuations inside the LCFS start earlier and are 2-3x larger.

Workaround: disable the filter (all three `ts_*` = -1). Proposed fix: apply the filter only to the
y-fluctuating part of the twist-shifted plane, which preserves ky = 0 exactly (see "Proposed fix").

## What the code does (`main` @ `e6b37b0e2`)

- `gyrokinetic/apps/gkyl_gyrokinetic.h`, `struct gkyl_gyrokinetic_core_parallel_bcs`:
  - `ts_upsample_factor`: default 4;
  - `ts_filter_half_width`: default 1, in cells of the simulation grid;
  - `ts_filter_cutoff_wavelength`: default 2*dx;
  - any of the three < 0 disables the filter (`gyrokinetic/zero/bc_twistshift.c`, `gkyl_bc_twistshift_inew`).
- `gyrokinetic/zero/bc_twistshift.c`: the ghost plane is supersampled along `shear_dir` (x) by
  `upsample_factor`, twist-shifted (`twistshift_dg.c`), filtered with `gkyl_dg_lowpass_filter`
  along `shear_dir`, then restricted back.
- `core/zero/dg_lowpass_filter.c`:
  - Blackman-windowed sinc kernel along one direction, normalized to preserve the zero mode.
  - At faces that coincide with the domain boundary, the stencil is reflected.
  - At interior edges of the update range (e.g. the LCFS cut), donors are dropped and the weights
    renormalized.
- Callers that pass the filter parameters:
  - species TS updaters in `gyrokinetic/apps/gk_species.c` (twist-shift-only and IWL/LCFS
    configurations);
  - field TS updaters in `gyrokinetic/apps/gk_field_2x3x.c` (phi parallel smoothing with TS ghost
    cells, twist-shift-only and IWL).
- `gyrokinetic/unit/ctest_bc_twistshift.c` checks that an identity filter reproduces the plain TS,
  that the restrict/prolong round trip is lossless, and that the filter smooths along the shear
  direction. **No test checks that the y-average is preserved.**

## Evidence

### 1. Minimal reproduction on `main` (CPU, local)

Test: `gyrokinetic/creg/rt_gk_tcv_iwl_adapt_source_3x2v_p1` (ES, IWL, Nx = 15 with the LCFS at
1/3 of the domain, Ny = Nz = 8, Nvpar = Nmu = 4, t_end = 1 us).
- The stock input uses the default filter.
- A copy with `.ts_upsample_factor = -1, .ts_filter_half_width = -1, .ts_filter_cutoff_wavelength = -1`
  in `core_parallel_bcs` serves as the control.

Diagnostic: per x cell, |<phi>_y(z_up edge) - <phi>_y(z_lo edge)| / max|<phi>_y|.
For 3x p1 serendipity, <f>_y at the z edges = (c0 +- sqrt(3) c3)/(2 sqrt(2)), averaged over y cells.

| frame | filter on, core cells (5) | filter off, core cells | SOL cells (both, free ends) |
|---|---|---|---|
| 0 (t = 0) | 1e-3 - 8.9e-2 | <= 1.6e-16 | 1e-2 - 2e-2 |
| 1 (t = 1 us) | 1.6e-3 - 2.9e-2 | <= 1.7e-16 | 1e-2 - 4e-2 |

After 1 us, filter on vs off: phi differs by up to 19 % (max relative), n_e by 1.1 %, ion density by
0.15 %. The jump is already present at t = 0, because the initial field solve also uses the
filtered TS ghost cells.

Conclusion: the defaults on `main` break the ky = 0 TS consistency, and disabling the filter restores
it to round-off.

### 2. TCV 0.5 MW coarse production runs (Perlmutter, GPU, ES)

Same input (3x2v p1, 24x16x12x12x8, x_LCFS = 4 cm out of Lx = 12 cm), differing only in the TS
filter:
- Reference, no filter: `/pscratch/sd/a/ah1032/gkeyll_em/tcv_pt_0.5MW_coarse/electrostatic/wk`
  (built before the filter existed; frames 0-49 come from an older input).
- Filter branch, `fc3c696ed+` (2026-09-09, branch `gk_tsbc_filter_fix`):
  `/pscratch/sd/a/ah1032/gkeyll_tsbc_filter/tcv_pt_coarse_0.5mw/wk`, with
  `ts_upsample_factor = 4`, `ts_filter_half_width = 4`, `ts_filter_cutoff_wavelength = 2*Lx/Nx`.
  The stencil then spans 9 coarse cells (4.5 cm), wider than the 4 cm core region.

y-averaged phi jump across z = +-pi inside the LCFS, filter run, frames 0-175 (0-350 us):

| core cell | jump |
|---|---|
| last core cell (next to the LCFS) | 9-16 % (from frame 1 on) |
| first cell (x = 0, reflected stencil) | 1-3.5 % |
| interior | 1e-4 - 5e-3 |
| reference run, all core cells | ~1e-15 |

Outboard-midplane, y-averaged profiles at t = 200-350 us (ranges over frames 100, 150, 175):

| | no filter | filter |
|---|---|---|
| phi at x = 0.75 cm | -180 to -235 V | -83 to -117 V |
| phi at x = 2.75 cm | -43 to +11 V | +113 to +256 V |
| n_e at x = 0.75 cm | 6.0-6.1e19 m^-3 | 5.4-5.5e19 m^-3 |
| T_e at x = 3.75 cm (inside LCFS) | 95-119 eV | 194-216 eV |
| rms dphi at x = 2.75 cm | 7-16 V | 13-39 V |
| rms dphi at x = 3.75 cm | 16-21 V | 41-71 V |

At t = 100 us the filter run already has rms dphi = 7.5 V at x = 2.75 cm, versus ~0 in the reference.

### 3. Controlled restart (same state, filter on vs off)

The reference ES run was restarted from frame 390 (t = 780 us) for 20-28 us in the debug QOS, in
`/pscratch/sd/a/ah1032/gkeyll_em/tcv_pt_0.5MW_coarse/electrostatic/wk_tsfilter_{on,off}`.
- Build: `get_em_back` `f62c83dfb`, i.e. `main` @ `e6b37b0e2` plus EM-only changes. The run is ES,
  so the EM code is inactive.
- Input: TS shift functions moved to `core_parallel_bcs` and `SHEATH` renamed to `SHEATH_CONDUCTING`.
  The only difference between the two restarts is the three `ts_*` = -1.

| | filter on | filter off |
|---|---|---|
| y-avg phi TS jump in core | 6-8 % from the restart frame | ~1e-15 |
| max\|phi\| at t = 782 us (core Er well) | 170 V (from 265 V) | 266 V |
| n_e,max over 20 us | 6.43 -> 5.90e19 (-8 %) | 6.43 -> 6.46e19 |
| rms dphi at t = 800 us | 12.6 V | 11.5 V |
| RK stage-2 dt reductions | 1714 of 4361 steps (dt -1 % to -6.8 %) | 0 of 6002 |

The filter-off restart matches the original reference run to the digits shown at the two
overlapping frames (391, 392). So the other changes merged into `main` since the reference run (the
sheath BC rewrite among them) do not change this simulation; the filter alone does.

## Minimal tests

**What a test needs.** The shift moves the field along y at each x, so it preserves the y-average of
every x cell exactly. The filter is linear and acts only along x. Therefore:

  <out>_y(x) = F_x[ <in>_y ](x)

- The error comes only from the ky = 0 part. A pure ky != 0 field keeps a zero y-average, so a
  ky != 0 initial condition does not show the problem.
- What exposes it is a y-averaged profile with radial structure.
- A y-independent field is the cleanest case, because the exact answer is known: the twist-shifted
  ghost plane must equal the skin plane.

### Unit test: `test_bc_twistshift_preserves_y_average_ho`

Patch against `main`: `tsbc_filter_y_average_test.patch`, which adds the test to
`gyrokinetic/unit/ctest_bc_twistshift.c`.

- **Donor.** f = 1 + tanh((x - 0.4)/0.25), with dx = 0.5 on the existing 8x8x4 test grid, optionally
  plus 0.5 cos(2 pi y/Ly) exp(-x^2) of ky != 0 content.
- **Shift.** offset 0.75, shear 2.9.
- **Check.** For every x cell, the y-sums of the y-independent DG coefficients (0, 1, 3, 5 of the 3x
  p1 serendipity basis) must agree between ghost and skin planes to 1e-12, for both edges.

| configuration | relative change of the y-average |
|---|---|
| no filter | 5-6e-16 (pass) |
| default filter (upsample 4, half-width 1, cutoff 2dx) | 3.3e-2 (fail) |
| filter (4, 4, 2dx) | 3.0e-2 (fail) |

- The result is the same with and without the ky != 0 content.
- All other tests in `ctest_bc_twistshift` still pass.
- The existing `test_bc_twistshift_conserves_particles_ho` checks only the integral over the whole
  ghost plane. The filter conserves that, so that test cannot catch the radial redistribution.

### One-step GK check (stock regression test, < 1 s)

```sh
./build/gyrokinetic/creg/rt_gk_tcv_iwl_adapt_source_3x2v_p1 -s 1 -x 15 -y 4 -z 4 -u 4 -v 2
```

- Evaluate the y-averaged phi z-edge mismatch of `field_0` in the core cells (x < x_LCFS, the first
  Nx/3 cells), with the snippet in "Reproduction" below.
- Frame 0 already shows it, because the initial field solve fills the phi ghost cells with the
  filtered twist-shift.

| 15x4x4x4x2, frame 0 | max core jump |
|---|---|
| stock (default filter) | 6.3e-2 |
| `ts_*` = -1 | 8e-16 |

At t = 0 the test is quasineutral and phi is round-off sized (~1e-11 V). The relative jump still
discriminates, but the physically meaningful check is frame 1 of the full 1 us test, where
|phi| ~ 700 V: 2.9e-2 with the default filter versus 1.7e-16 with the filter off.

Nx = 15 is the smallest usable value here. Nx = 12 trips the `ts_calc_mats` assertion: the shift
varies by more than Ly across one SOL cell at this shear.

## Workaround

Disable the filter in inputs with closed field lines:

```c
.core_parallel_bcs = {
  .lower_shift_func = ..., .upper_shift_func = ...,
  .lower_shift_ctx = ..., .upper_shift_ctx = ...,
  .ts_upsample_factor = -1,
  .ts_filter_half_width = -1,
  .ts_filter_cutoff_wavelength = -1,
},
```

## Proposed fix: filter only the y-fluctuating part

The filter is meant to de-alias the sheared, shifted ky != 0 content. The ky = 0 part is mapped
exactly by the shift and needs no de-aliasing. Keep it out of the filter.

With P the prolongation (supersampling in x), T the twist-shift on the fine grid, F the low-pass
filter along x, R the restriction, and <.>_y the average along `shift_dir`:

  current:   f_ghost = R F T P f_skin
  proposed:  f_ghost = R ( <g>_y + F[ g - <g>_y ] ),   g = T P f_skin

- **Exact y-average.** <T g>_y = <g>_y exactly, and P and R act in x only, so they commute with
  <.>_y. Hence <f_ghost>_y = R P <f_skin>_y = <f_skin>_y; the prolong/restrict round trip is
  lossless (`ctest_bc_twistshift`). The y-average is preserved to round-off for any filter
  parameters.
- **De-aliasing kept.** The ky != 0 content is filtered exactly as now.
- **Equivalent cheaper form.** F is linear and <g>_y is y-independent, so

    f_ghost = R ( F[g] + ( <g>_y - F[<g>_y] ) ).

  The correction <g>_y - F[<g>_y] lives on the y-reduced plane (2D for fields, 4D for f in 3x2v).
  It is computed with the existing filter on the reduced array and added back to every y cell.
- **Location.** All in `gyrokinetic/zero/bc_twistshift.c`, between the shift and the restriction. It
  then applies to f, phi and Apar alike, with no change to the callers or the boundary layout.
- **Ingredients.**
  - The y-average of a DG field keeps only the coefficients without y dependence (0, 1, 3, 5 for 3x
    p1 serendipity; the analogous subset for the 3x2v hybrid basis). `gkyl_array_average` can
    reduce along y.
  - Adding the reduced correction back to every y cell needs a broadcast; this may need a small
    kernel (CPU and GPU).
- **Validation.**
  - `test_bc_twistshift_preserves_y_average_ho` must pass with the filter on.
  - `rt_gk_tcv_iwl_adapt_source_3x2v_p1 -s 1 -x 15 -y 4 -z 4 -u 4 -v 2` must give a core phi jump of
    about 1e-16 with default parameters.
  - `test_bc_twistshift_dealiasing_smooths_shear_direction_ho` must still show the smoothing. Its
    donor has a ky = 0 part too, so its metric may need restricting to the fluctuating part.
  - The TCV ES restart (filter on vs off) should then give matching profiles and Er well.

**Implementation and validation** (branch `fix_of_the_tsbc_filter_fix_fix`, on top of `main` @ `e6b37b0e2`)

- **Code.**
  - `bc_twistshift.c`: `bc_twistshift_shift_dir_avg` computes `favg`, the shift_dir average of the
    shifted fine ghost plane, kept in the full basis and repeated in every shift_dir cell. The
    filtered path then subtracts `favg`, filters, and adds it back. The basis functions independent
    of shift_dir are detected once, at creation, with `basis->eval`, so any modal product basis
    works (3x serendipity, 3x2v gkhybrid).
  - `bc_twistshift_cu.cu`: the GPU kernel.
  - `gkyl_bc_twistshift_priv.h`: the new members.
- **Unit tests.** All 10 tests of `ctest_bc_twistshift` pass, including the new
  `test_bc_twistshift_preserves_y_average_ho` with the filter on and the existing de-aliasing test.
  `ctest_twistshift_dg` passes.
- **`rt_gk_tcv_iwl_adapt_source_3x2v_p1`** (default filter), max core jump of the y-averaged phi at
  frame 1 (1 us): 2.9e-2 on `main`, 2.6e-16 with the fix (filter off: 1.7e-16).
  - With the filter off, results are bitwise identical to `main`.
  - At t = 0 the fixed run matches the filter-off run to round-off: all f files identical, phi
    within 5e-13 V, and the same first two time steps. (On `main`, the first time step already
    differs by 11 %.)
  - The end states at 1 us still differ, because this test's strong early transient (over 110
    RK stage-2 retakes in 160 steps) amplifies round-off.
- **`rt_gk_cbc_3x2v_p1`** (twist-shift, no LCFS, ky != 0 initial condition), max jump of the
  y-averaged phi at 1 us: 7.4e-2 on `main`, 2.5e-16 with the fix.
  - phi differs from the filter-off run by up to 25x (max relative) on `main`, and by 1.1 % with
    the fix, which is the intended de-aliasing of the ky != 0 content.
  - A jump of about 1e-6 at t = 0 is also there with the filter off, so it is unrelated to the fix.
- **GPU (Perlmutter A100, branch `a569ddd8b`, install `/pscratch/sd/a/ah1032/installs/gkylsoft_tsbc_fix`).**
  - `ctest_bc_twistshift`: all 10 tests pass.
  - `rt_gk_tcv_iwl_adapt_source_3x2v_p1 -g`, max core jump of the y-averaged phi at frame 1:
    1.3e-15 with the default filter and 1.3e-15 with it off.
  - `rt_gk_cbc_3x2v_p1 -g`, frame 1: 1.3e-15 with the default filter and 2.3e-15 with it off.
  - TCV ES restart from frame 390, 26 us on 4 GPUs (`electrostatic/wk_tsfilter_on_fixed`), against
    the filter-off restart (`wk_tsfilter_off`):

    | | default filter, fixed | filter off | default filter, `main` |
    |---|---|---|---|
    | core phi twist-shift jump | ~1e-15 | ~1e-15 | 6-8 % |
    | max\|phi\| at 800 us | 274.5 V | 276.0 V | 212 V |
    | n_e,max at 800 us | 6.455e19 | 6.455e19 | 5.94e19 |
    | RK stage-2 retakes | 0 of 6058 | 0 of 6002 | 1714 of 4361 |

    - The zonal phi and densities stay within 0.1-0.2 % of the filter-off run after 4 us, and 2-4 %
      after 26 us. With the old filter, the zonal phi in the core differed by about 100 % within 4 us.
    - The fluctuations decorrelate from the filter-off run over the 26 us, as expected once the
      filter acts on the ky != 0 content.

**A more symmetric application does not help.** In `main`, f is already twist-shifted at both
z ends (`bc_ts_lo` and `bc_ts_up` in `gk_species.c`); phi and Apar are twist-shifted at the lower end
only, and the Dirichlet-ghost FEM projection closes the loop. Every filtered application gives
<f_ghost>_y = F_x[<f_skin>_y], with F_x a smoothing operator (gain <= 1):
- splitting the shift between both ends doubles the number of smoothings per transit;
- a weaker filter at each end (e.g. sqrt(F)) gives the same total damping per transit.

On closed field lines, any non-identity operator on the ky = 0 part compounds over repeated passes
into a radial diffusion; only an operator that is the identity on ky = 0 avoids it. Symmetry could
still restore reversibility: with the filter, the operator out through one end is not the inverse
of the one back in, so a round trip does not return f unchanged. That is secondary to the fix above.

## Open questions and suggested investigation on `main`

1. **Preserve ky = 0.** Implement and validate the fix above. An alternative to compare it with is
   to filter only the difference between the supersampled twist-shift and the plain one.
2. **Default off** until the above is validated. Regression references generated with the current
   defaults (e.g. `rt_gk_tcv_iwl_adapt_source_{2x2v,3x2v}_p1`) would need regenerating.
3. **Edge treatment.**
   - The largest errors are at the LCFS cut (dropped donors, renormalized weights) and at x = 0
     (reflected stencil). Check what the reflection does to the odd DG coefficients (`sign_mirror`)
     and whether the ghost-plane integral (particle content) is conserved at both edges.
   - `ts_filter_half_width` larger than the closed-field-line width should probably be rejected.
4. **Conservation.** Measure the change of the particle and energy content of the ghost plane per
   application (filter on/off), and the resulting radial flux at the z = +-pi boundary on closed
   field lines. That flux is the numerical diffusion estimated above.
5. **Tests.**
   - Add a unit test in `ctest_bc_twistshift.c`: the y-average of the TS ghost plane equals the
     y-average of the donor skin plane, for any shift and with the filter on. A draft is provided;
     see "Minimal tests" above.
   - Add a regression check: in an IWL test, the y-averaged phi is periodic in z inside the LCFS to
     round-off.
6. **Interaction with the parallel smoothing of phi.** The field TS updaters fill the ghost cells
   used by the Dirichlet-ghost FEM parallel projection. Even if f were left untouched, filtering
   phi's ghost plane alone breaks the zonal phi continuity.

## Reproduction

- **Local (CPU).** Build `main`, then compile `rt_gk_tcv_iwl_adapt_source_3x2v_p1` and a copy with the
  three `ts_*` = -1 added to `core_parallel_bcs`. Run both (~2 min each on a laptop) and evaluate
  the edge values of the y-averaged phi:

```python
import numpy as np, postgkyl as pg
d = pg.GData("<name>-field_1.gkyl"); v = d.get_values()   # 3x p1 serendipity, shape (Nx, Ny, Nz, 8)
lo = ((v[..., 0][:, :, 0]  - np.sqrt(3)*v[..., 3][:, :, 0]) /(2*np.sqrt(2))).mean(axis=1)
up = ((v[..., 0][:, :, -1] + np.sqrt(3)*v[..., 3][:, :, -1])/(2*np.sqrt(2))).mean(axis=1)
print(np.abs(up - lo)/np.abs(lo).max())                   # cells x < x_LCFS should be ~1e-16
```

- **Production data.** The Perlmutter directories listed above. The commit and build date of every
  frame are in the file metadata (`pgkyl <file> info`, or `pg.GData(f).ctx["changeset"]`).
