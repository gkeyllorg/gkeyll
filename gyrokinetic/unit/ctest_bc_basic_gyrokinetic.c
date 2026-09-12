// Tests for the gyrokinetic basic-BC updater (gkyl_bc_basic_gyrokinetic).
//
// 1) Constructor field checks: dir/edge/cdim/bctype recorded, skin/ghost
//    range pointers stored.
// 2) COPY-BC advance compute check: fill the skin cell with known values,
//    run the updater, and verify the ghost cell ends up holding an exact
//    copy of the skin data (the COPY boundary condition).
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_bc_basic_gyrokinetic.h>
#include <gkyl_bc_basic_gyrokinetic_priv.h>
#include <gkyl_gk_bc_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

// Build a 1D conf grid + skin/ghost ranges at the given edge.
static void
setup_1x(enum gkyl_edge_loc edge, struct gkyl_basis *basis,
  struct gkyl_range *local, struct gkyl_range *local_ext,
  struct gkyl_range *skin_r, struct gkyl_range *ghost_r)
{
  int cells[] = { 8 };
  int ghost[] = { 1 };
  double lower[] = { 0.0 }, upper[] = { 1.0 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, local_ext, local);

  gkyl_cart_modal_serendip(basis, 1, 1);

  gkyl_skin_ghost_ranges(skin_r, ghost_r, 0, edge, local_ext, ghost);
}

void
test_bc_basic_ctor()
{
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext, skin_r, ghost_r;
  setup_1x(GKYL_LOWER_EDGE, &basis, &local, &local_ext, &skin_r, &ghost_r);

  struct gkyl_bc_basic_gyrokinetic *bc = gkyl_bc_basic_gyrokinetic_new(
    0, GKYL_LOWER_EDGE, GKYL_BC_GK_SPECIES_COPY, &basis, &skin_r, &ghost_r,
    basis.num_basis, 1, false);

  TEST_CHECK( bc != NULL );
  TEST_CHECK( bc->dir == 0 );
  TEST_CHECK( bc->cdim == 1 );
  TEST_CHECK( bc->edge == GKYL_LOWER_EDGE );
  TEST_CHECK( bc->bctype == GKYL_BC_GK_SPECIES_COPY );
  TEST_CHECK( bc->skin_r == &skin_r );
  TEST_CHECK( bc->ghost_r == &ghost_r );
  TEST_CHECK( bc->use_gpu == false );
  TEST_CHECK( bc->array_copy_func != NULL );

  gkyl_bc_basic_gyrokinetic_release(bc);
}

// Verify that applying a COPY BC fills the ghost cell with an exact copy of
// the skin cell.
static void
check_copy_advance(enum gkyl_edge_loc edge)
{
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext, skin_r, ghost_r;
  setup_1x(edge, &basis, &local, &local_ext, &skin_r, &ghost_r);

  int nc = basis.num_basis; // 2 for 1x p1

  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, nc, local_ext.volume);
  gkyl_array_clear(f, 0.0);

  // Buffer big enough to hold the skin cells.
  long buff_sz = skin_r.volume;
  struct gkyl_array *buff = gkyl_array_new(GKYL_DOUBLE, nc, buff_sz);
  gkyl_array_clear(buff, 0.0);

  // Fill the skin cell with known, distinct values.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skin_r);
  while (gkyl_range_iter_next(&iter)) {
    long sidx = gkyl_range_idx(&skin_r, iter.idx);
    double *fs = gkyl_array_fetch(f, sidx);
    for (int c=0; c<nc; c++)
      fs[c] = 10.0 + c + 0.5;
  }

  struct gkyl_bc_basic_gyrokinetic *bc = gkyl_bc_basic_gyrokinetic_new(
    0, edge, GKYL_BC_GK_SPECIES_COPY, &basis, &skin_r, &ghost_r, nc, 1, false);

  gkyl_bc_basic_gyrokinetic_advance(bc, buff, f);

  // Ghost cell should now equal the skin values.
  gkyl_range_iter_init(&iter, &ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    long gidx = gkyl_range_idx(&ghost_r, iter.idx);
    const double *fg = gkyl_array_cfetch(f, gidx);
    for (int c=0; c<nc; c++)
      TEST_CHECK( gkyl_compare(fg[c], 10.0 + c + 0.5, 1e-14) );
  }

  gkyl_bc_basic_gyrokinetic_release(bc);
  gkyl_array_release(f);
  gkyl_array_release(buff);
}

void test_bc_basic_copy_lower() { check_copy_advance(GKYL_LOWER_EDGE); }
void test_bc_basic_copy_upper() { check_copy_advance(GKYL_UPPER_EDGE); }

void
test_bc_basic_ctor_upper_reflect()
{
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext, skin_r, ghost_r;
  setup_1x(GKYL_UPPER_EDGE, &basis, &local, &local_ext, &skin_r, &ghost_r);

  struct gkyl_bc_basic_gyrokinetic *bc = gkyl_bc_basic_gyrokinetic_new(
    0, GKYL_UPPER_EDGE, GKYL_BC_GK_SPECIES_REFLECT, &basis, &skin_r, &ghost_r,
    basis.num_basis, 1, false);

  TEST_CHECK( bc != NULL );
  TEST_CHECK( bc->edge == GKYL_UPPER_EDGE );
  TEST_CHECK( bc->bctype == GKYL_BC_GK_SPECIES_REFLECT );

  gkyl_bc_basic_gyrokinetic_release(bc);
}

TEST_LIST = {
  { "bc_basic_ctor", test_bc_basic_ctor },
  { "bc_basic_copy_lower", test_bc_basic_copy_lower },
  { "bc_basic_copy_upper", test_bc_basic_copy_upper },
  { "bc_basic_ctor_upper_reflect", test_bc_basic_ctor_upper_reflect },
  { NULL, NULL },
};
