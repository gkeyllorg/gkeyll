// Tests for the PKPM distribution-function variables updater
// (gkyl_dg_calc_pkpm_dist_vars).
//
// The div(p_par b) operator is a recovery-based divergence: for a distribution
// function that is spatially *uniform* (only the cell-average component is
// non-zero, identical in every cell) and a uniform magnetic-field unit vector
// (b = x_hat), the parallel pressure p_par is constant and b is constant, so
//   div(p_par b) == 0
// identically.  This is a strong, physically-motivated check that exercises the
// full surface-recovery machinery and confirms it preserves a constant state.
//
// We also include a construction smoke test in 2x2v.
#include <acutest.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_dist_vars.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

void
test_div_ppar_uniform_1x1v_p1()
{
  int poly_order = 1;
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  double lower[] = {-2.0, -3.0}, upper[] = {2.0, 3.0};
  int cells[] = {8, 6};

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  struct gkyl_rect_grid grid, confGrid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_hybrid(&basis, cdim, vdim); // 6 components for 1x1v
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order); // 2 components

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_dg_calc_pkpm_dist_vars *up =
    gkyl_dg_calc_pkpm_dist_vars_new(&grid, &confBasis, false);
  TEST_CHECK( up != NULL );

  int nbp = basis.num_basis;     // 6
  int nbc = confBasis.num_basis; // 2

  // fIn: [F_0, T_perp/m G] -> 2 phase blocks.
  struct gkyl_array *fIn = mkarr(2*nbp, local_ext.volume);
  // bvar volume expansion: 9 conf blocks [bx,by,bz,bxbx,...,bzbz].
  struct gkyl_array *bvar = mkarr(9*nbc, confLocal_ext.volume);
  // bvar_surf: 2*cdim*4 * Nbasis_surf = 8 * 1 = 8 components.
  struct gkyl_array *bvar_surf = mkarr(8, confLocal_ext.volume);
  // max_b: 2*cdim*Nbasis_surf = 2 components.
  struct gkyl_array *max_b = mkarr(2, confLocal_ext.volume);
  struct gkyl_array *div_ppar = mkarr(nbc, confLocal_ext.volume);

  // Uniform F_0: only the phase cell-average component (index 0 of block 0).
  // Any uniform value gives div = 0; we pick a representative non-zero value.
  gkyl_array_clear(fIn, 0.0);
  gkyl_array_shiftc(fIn, 1.3, 0);      // F_0 cell-average component
  gkyl_array_shiftc(fIn, 0.5, nbp);    // G cell-average (unused by div_ppar but realistic)

  // Uniform b = x_hat: bx = 1, bxbx = 1, all else 0 (volume expansion).
  gkyl_array_clear(bvar, 0.0);
  gkyl_array_shiftc(bvar, sqrt(2.0), 0*nbc);     // bx cell-average physical value 1
  gkyl_array_shiftc(bvar, sqrt(2.0), 3*nbc);     // bxbx cell-average physical value 1

  // Surface b expansion: [bx_xl, bx_xr, bxbx_xl, bxbx_xr, ...]; set b = 1 on
  // both surfaces and the surface unit tensor likewise.
  gkyl_array_clear(bvar_surf, 0.0);
  gkyl_array_shiftc(bvar_surf, 1.0, 0); // bx_xl
  gkyl_array_shiftc(bvar_surf, 1.0, 1); // bx_xr
  gkyl_array_shiftc(bvar_surf, 1.0, 2); // bxbx_xl
  gkyl_array_shiftc(bvar_surf, 1.0, 3); // bxbx_xr

  // max_b = |b| = 1 on each edge.
  gkyl_array_clear(max_b, 0.0);
  gkyl_array_shiftc(max_b, 1.0, 0);
  gkyl_array_shiftc(max_b, 1.0, 1);

  gkyl_array_clear(div_ppar, 0.0);
  gkyl_dg_calc_pkpm_dist_vars_div_ppar(up, &confLocal, &local,
    bvar_surf, bvar, fIn, max_b, div_ppar);

  // For a uniform state div(p_par b) must be identically zero in every cell
  // and every basis component.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&confLocal, iter.idx);
    const double *d = gkyl_array_cfetch(div_ppar, loc);
    for (int k=0; k<nbc; ++k)
      TEST_CHECK( gkyl_compare(d[k], 0.0, 1e-12) );
  }

  gkyl_array_release(div_ppar);
  gkyl_array_release(max_b);
  gkyl_array_release(bvar_surf);
  gkyl_array_release(bvar);
  gkyl_array_release(fIn);
  gkyl_dg_calc_pkpm_dist_vars_release(up);
}

void
test_dist_vars_new_2x2v_p1()
{
  int poly_order = 1;
  int cdim = 2, vdim = 2, pdim = cdim+vdim;
  double lower[] = {-1.0, -1.0, -2.0, -2.0}, upper[] = {1.0, 1.0, 2.0, 2.0};
  int cells[] = {4, 4, 4, 4};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);

  struct gkyl_basis confBasis;
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_dg_calc_pkpm_dist_vars *up =
    gkyl_dg_calc_pkpm_dist_vars_new(&grid, &confBasis, false);
  TEST_CHECK( up != NULL );
  gkyl_dg_calc_pkpm_dist_vars_release(up);
}

TEST_LIST = {
  { "div_ppar_uniform_1x1v_p1", test_div_ppar_uniform_1x1v_p1 },
  { "dist_vars_new_2x2v_p1",    test_dist_vars_new_2x2v_p1 },
  { NULL, NULL },
};
