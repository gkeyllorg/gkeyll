// Test construction of the gyrokinetic boundary-flux updater
// (gkyl_boundary_flux_new). Verifies that the stored direction, edge,
// number of equations, copied grid, and skin/ghost range volumes match the
// inputs, and that the equation object is acquired (so it survives release of
// the caller's reference). A lightweight gyrokinetic-diffusion DG equation is
// used as the equation object.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_boundary_flux.h>
#include <gkyl_boundary_flux_priv.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion_gyrokinetic.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static struct gkyl_dg_eqn*
make_diffusion_eqn(struct gkyl_range *diff_range)
{
  int cdim = 1, vdim = 1, poly_order = 1;
  struct gkyl_basis basis, cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);

  gkyl_range_init(diff_range, cdim, (int[]) { 0 }, (int[]) { 7 });
  bool dir[] = { true };
  return gkyl_dg_diffusion_gyrokinetic_new(&basis, &cbasis, true, dir, 2, diff_range, false);
}

static void
check_boundary_flux(int dir, enum gkyl_edge_loc edge)
{
  // Phase-space grid: 1x1v.
  int cells[] = { 8, 8 };
  int ghost[] = { 1, 0 };
  double lower[] = { 0.0, -1.0 }, upper[] = { 1.0, 1.0 };

  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  struct gkyl_range diff_range;
  struct gkyl_dg_eqn *eqn = make_diffusion_eqn(&diff_range);

  const struct gkyl_dg_eqn *eqns[] = { eqn };
  struct gkyl_boundary_flux *bf = gkyl_boundary_flux_new(dir, edge, &grid,
    &skin_r, &ghost_r, 1, eqns, false);

  TEST_CHECK( bf != NULL );
  TEST_CHECK( bf->dir == dir );
  TEST_CHECK( bf->edge == edge );
  TEST_CHECK( bf->num_eqns == 1 );
  TEST_CHECK( bf->use_gpu == false );
  TEST_CHECK( bf->grid.ndim == 2 );
  TEST_CHECK( bf->grid.cells[0] == 8 );
  TEST_CHECK( bf->skin_r.volume == skin_r.volume );
  TEST_CHECK( bf->ghost_r.volume == ghost_r.volume );
  TEST_CHECK( bf->eqns[0] != NULL );

  // The updater acquired its own reference; release ours and confirm the
  // equation object is still alive (num_equations readable).
  gkyl_dg_eqn_release(eqn);
  TEST_CHECK( bf->eqns[0]->num_equations == 1 );

  gkyl_boundary_flux_release(bf);
}

void test_boundary_flux_lower() { check_boundary_flux(0, GKYL_LOWER_EDGE); }
void test_boundary_flux_upper() { check_boundary_flux(0, GKYL_UPPER_EDGE); }

TEST_LIST = {
  { "boundary_flux_lower", test_boundary_flux_lower },
  { "boundary_flux_upper", test_boundary_flux_upper },
  { NULL, NULL },
};
