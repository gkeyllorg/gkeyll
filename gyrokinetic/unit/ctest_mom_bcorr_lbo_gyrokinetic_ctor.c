// Test construction of the gyrokinetic LBO boundary-correction moment type.
// Verifies the returned gkyl_mom_type carries the expected dim/basis fields,
// that num_mom == 2 (parallel-velocity + energy corrections), and that
// acquire/release reference counting works. Uses an identity velocity map.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>

static struct gkyl_velocity_map*
make_identity_vmap(int cdim, int vdim, double *plower, double *pupper, int *pcells)
{
  int pdim = cdim + vdim;

  double vlower[3], vupper[3];
  int vcells[3];
  for (int d=0; d<vdim; d++) {
    vlower[d] = plower[cdim+d];
    vupper[d] = pupper[cdim+d];
    vcells[d] = pcells[cdim+d];
  }

  struct gkyl_rect_grid grid, vgrid;
  gkyl_rect_grid_init(&grid, pdim, plower, pupper, pcells);
  gkyl_rect_grid_init(&vgrid, vdim, vlower, vupper, vcells);

  int ghost[GKYL_MAX_DIM] = {0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vghost[3] = {0};
  struct gkyl_range vlocal, vlocal_ext;
  gkyl_create_grid_ranges(&vgrid, vghost, &vlocal_ext, &vlocal);

  struct gkyl_mapc2p_inp c2p_in = { };
  return gkyl_velocity_map_new(c2p_in, grid, vgrid, local, local_ext,
    vlocal, vlocal_ext, false);
}

static void
check_bcorr(int cdim, int vdim, int poly_order, double *plower, double *pupper, int *pcells)
{
  int pdim = cdim + vdim;

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&pbasis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_velocity_map *gvm = make_identity_vmap(cdim, vdim, plower, pupper, pcells);

  double mass = 1.25;
  struct gkyl_mom_type *bcorr = gkyl_mom_bcorr_lbo_gyrokinetic_new(&cbasis, &pbasis,
    mass, gvm, false);

  TEST_CHECK( bcorr != NULL );
  TEST_CHECK( bcorr->cdim == cdim );
  TEST_CHECK( bcorr->pdim == pdim );
  TEST_CHECK( bcorr->poly_order == poly_order );
  TEST_CHECK( bcorr->num_config == cbasis.num_basis );
  TEST_CHECK( bcorr->num_phase == pbasis.num_basis );
  // Two boundary-correction moments: parallel velocity and energy.
  TEST_CHECK( bcorr->num_mom == 2 );
  TEST_CHECK( gkyl_mom_type_num_mom(bcorr) == 2 );

  gkyl_mom_type_release(bcorr);
  gkyl_velocity_map_release(gvm);
}

void test_bcorr_1x1v_p1()
{
  double lo[] = {-1.0, -6.0}, up[] = {1.0, 6.0};
  int cells[] = {4, 8};
  check_bcorr(1, 1, 1, lo, up, cells);
}
void test_bcorr_1x2v_p1()
{
  double lo[] = {-1.0, -6.0, 0.0}, up[] = {1.0, 6.0, 36.0};
  int cells[] = {4, 8, 4};
  check_bcorr(1, 2, 1, lo, up, cells);
}
void test_bcorr_2x2v_p1()
{
  double lo[] = {-1.0, -1.0, -6.0, 0.0}, up[] = {1.0, 1.0, 6.0, 36.0};
  int cells[] = {4, 4, 8, 4};
  check_bcorr(2, 2, 1, lo, up, cells);
}
void test_bcorr_3x2v_p1()
{
  double lo[] = {-1.0, -1.0, -1.0, -6.0, 0.0}, up[] = {1.0, 1.0, 1.0, 6.0, 36.0};
  int cells[] = {2, 2, 2, 8, 4};
  check_bcorr(3, 2, 1, lo, up, cells);
}

void
test_bcorr_acquire()
{
  double lo[] = {-1.0, -6.0, 0.0}, up[] = {1.0, 6.0, 36.0};
  int cells[] = {4, 8, 4};

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_gkhybrid(&pbasis, 1, 2);

  struct gkyl_velocity_map *gvm = make_identity_vmap(1, 2, lo, up, cells);

  struct gkyl_mom_type *bcorr = gkyl_mom_bcorr_lbo_gyrokinetic_new(&cbasis, &pbasis,
    1.0, gvm, false);

  struct gkyl_mom_type *bcorr2 = gkyl_mom_type_acquire(bcorr);
  TEST_CHECK( bcorr2 == bcorr );

  gkyl_mom_type_release(bcorr);

  TEST_CHECK( bcorr2->num_mom == 2 );
  TEST_CHECK( bcorr2->cdim == 1 );

  gkyl_mom_type_release(bcorr2);
  gkyl_velocity_map_release(gvm);
}

TEST_LIST = {
  { "bcorr_1x1v_p1", test_bcorr_1x1v_p1 },
  { "bcorr_1x2v_p1", test_bcorr_1x2v_p1 },
  { "bcorr_2x2v_p1", test_bcorr_2x2v_p1 },
  { "bcorr_3x2v_p1", test_bcorr_3x2v_p1 },
  { "bcorr_acquire", test_bcorr_acquire },
  { NULL, NULL },
};
