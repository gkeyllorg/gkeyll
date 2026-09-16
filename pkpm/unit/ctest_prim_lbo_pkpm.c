// Test construction of PKPM primitive-moment (LBO) objects.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_pkpm.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static void
mk_conf_range(int cdim, struct gkyl_range *local, struct gkyl_range *local_ext)
{
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) { lower[d] = -1.0; upper[d] = 1.0; cells[d] = 4; ghost[d] = 1; }
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, local_ext, local);
}

void
test_prim_lbo_pkpm_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(1, &local, &local_ext);

  struct gkyl_prim_lbo_type *prim =
    gkyl_prim_lbo_pkpm_new(&cbasis, &pbasis, &local, false);

  TEST_CHECK( prim->cdim == 1 );
  TEST_CHECK( prim->pdim == 2 );
  TEST_CHECK( prim->poly_order == 1 );
  TEST_CHECK( prim->num_config == cbasis.num_basis );
  TEST_CHECK( prim->num_phase == pbasis.num_basis );
  TEST_CHECK( prim->udim == 1 );           // pkpm has 1D momentum (parallel)
  TEST_CHECK( prim->self_prim != NULL );

  gkyl_prim_lbo_type_release(prim);
}

void
test_prim_lbo_pkpm_1x1v_p2()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(1, &local, &local_ext);

  struct gkyl_prim_lbo_type *prim =
    gkyl_prim_lbo_pkpm_new(&cbasis, &pbasis, &local, false);
  TEST_CHECK( prim->poly_order == 2 );
  TEST_CHECK( prim->udim == 1 );
  TEST_CHECK( prim->self_prim != NULL );
  gkyl_prim_lbo_type_release(prim);
}

void
test_prim_lbo_pkpm_2x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 3, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(2, &local, &local_ext);

  struct gkyl_prim_lbo_type *prim =
    gkyl_prim_lbo_pkpm_new(&cbasis, &pbasis, &local, false);
  TEST_CHECK( prim->cdim == 2 );
  TEST_CHECK( prim->pdim == 3 );
  TEST_CHECK( prim->udim == 1 );
  gkyl_prim_lbo_type_release(prim);
}

TEST_LIST = {
  { "prim_lbo_pkpm_1x1v_p1", test_prim_lbo_pkpm_1x1v_p1 },
  { "prim_lbo_pkpm_1x1v_p2", test_prim_lbo_pkpm_1x1v_p2 },
  { "prim_lbo_pkpm_2x1v_p1", test_prim_lbo_pkpm_2x1v_p1 },
  { NULL, NULL },
};
