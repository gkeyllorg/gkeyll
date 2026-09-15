#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

void
test_dg_max()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  // Identity position map providing the (required) conf-space Jacobian.
  double lower[] = { 0.0 }, upper[] = { 1.0 };
  int cells[] = { 8 }, ghost[] = { 1 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range crange, crange_ext;
  gkyl_create_grid_ranges(&grid, ghost, &crange_ext, &crange);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&grid,
    &crange, &crange_ext, &basis, inp_pmap, false);

  // Input structure for building the dg eqn object
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &basis,
    .crange = &crange,
    .jacob_pos = pos_map->jacob_pos,
    .conf_flux_surf = 0,
    .lightSpeed = 1.0,
    .field_id = GKYL_FIELD_E_B,
    .elcErrorSpeedFactor = 0.5,
    .mgnErrorSpeedFactor = 0.25,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn* eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);

  TEST_CHECK( eqn->num_equations == 8 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  TEST_CHECK( maxwell->maxwell_data.c == 1.0 );
  TEST_CHECK( maxwell->maxwell_data.chi == 0.5 );
  TEST_CHECK( maxwell->maxwell_data.gamma == 0.25 );

  struct dg_maxwell *m_on_dev = container_of(eqn->on_dev,
    struct dg_maxwell, eqn);

  TEST_CHECK( m_on_dev->maxwell_data.c == 1.0 );
  TEST_CHECK( m_on_dev->maxwell_data.chi == 0.5 );
  TEST_CHECK( m_on_dev->maxwell_data.gamma == 0.25 );

  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_position_map_release(pos_map);
}

#ifdef GKYL_HAVE_CUDA

int cu_maxwell_test(const struct gkyl_dg_eqn *eqn);

void
test_cu_dg_max()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  // Identity position map providing the (required) conf-space Jacobian.
  double lower[] = { 0.0 }, upper[] = { 1.0 };
  int cells[] = { 8 }, ghost[] = { 1 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range crange, crange_ext;
  gkyl_create_grid_ranges(&grid, ghost, &crange_ext, &crange);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&grid,
    &crange, &crange_ext, &basis, inp_pmap, true);

  // Input structure for building the dg eqn object
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &basis,
    .crange = &crange,
    .jacob_pos = pos_map->jacob_pos,
    .conf_flux_surf = 0,
    .lightSpeed = 1.0,
    .field_id = GKYL_FIELD_E_B,
    .elcErrorSpeedFactor = 0.5,
    .mgnErrorSpeedFactor = 0.25,
    .use_gpu = true,
  };
  struct gkyl_dg_eqn* eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  TEST_CHECK( maxwell->maxwell_data.c == 1.0 );
  TEST_CHECK( maxwell->maxwell_data.chi == 0.5 );
  TEST_CHECK( maxwell->maxwell_data.gamma == 0.25 );
  
  // call CUDA test
  /* int nfail = cu_maxwell_test(eqn->on_dev); */

  /* TEST_CHECK( nfail == 0 ); */

  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_position_map_release(pos_map);
}

#endif

TEST_LIST = {
  { "dg_max", test_dg_max },
#ifdef GKYL_HAVE_CUDA
  { "cu_dg_max", test_cu_dg_max },
#endif  
  { NULL, NULL },
};
