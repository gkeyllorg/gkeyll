#include <acutest.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
test_dg_vlasov()
{
  int cdim = 1, vdim = 1, pdim = 2;
  struct gkyl_basis cbasis, vbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_serendip(&vbasis, 1, 1);
  gkyl_cart_modal_serendip(&pbasis, 2, 1);

  int cells[] = {100, 12 };
  int ghost[] = {1, 0 };
  double lower[] = {0., -1. };
  double upper[] = {1., 1. };

  int velcells[] = {12};
  int velghost[] = { 0};
  double vellower[] = { -1.};
  double velupper[] = { 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid velGrid;
  struct gkyl_range velRange, velRange_ext;
  gkyl_rect_grid_init(&velGrid, vdim, vellower, velupper, velcells);
  gkyl_create_grid_ranges(&velGrid, velghost, &velRange_ext, &velRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize eqn
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;

    // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr1(false, vbasis.num_basis, velRange.volume);
  struct gkyl_array *gamma_inv = mkarr1(false, vbasis.num_basis, velRange.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velRange, &vbasis, inp_vmap, false, false);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confRange, &confRange_ext, &cbasis, inp_pmap, false);

  gkyl_dg_vlasov_calc_hamil(&velGrid, &vbasis, &velRange, 
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, false);

  // Sturcture pointers for input objects (but not used)
  int num_pt_indices[3] = { 1 , 6, 18 }; 
  struct gkyl_array *poisson_tensor_conf = mkarr1(false, cbasis.num_basis*num_pt_indices[vdim-1], confRange.volume );
  struct gkyl_array *pot_tot = mkarr1(false, cbasis.num_basis*4, confRange_ext.volume );
  struct gkyl_array *vel_flux_surf = mkarr1(false, pbasis.num_basis*vdim, phaseRange_ext.volume );
  struct gkyl_array *f_no_J = mkarr1(false, pbasis.num_basis, phaseRange_ext.volume); ;
  struct gkyl_array *rad = mkarr1(false, vdim*vbasis.num_basis, velRange.volume);
  struct gkyl_array *qmem = mkarr1(false, 8*cbasis.num_basis, confRange_ext.volume);

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .conf_range =  &confRange,
    .hamil_range = &velRange,
    .phase_range = &phaseRange,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .poisson_tensor_conf = poisson_tensor_conf,
    .hamil = hamil,
    .qmem = qmem, 
    .pot_tot = pot_tot, 
    .vel_flux_surf = vel_flux_surf, 
    .f_no_J = f_no_J, 
    .rad = rad, 
    .use_lo = false,
    .use_gpu = false,
  };  
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  TEST_CHECK( eqn->num_equations == 1 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  TEST_CHECK( vlasov->cdim == 1 );
  TEST_CHECK( vlasov->pdim == 2 );
  TEST_CHECK( vlasov->conf_range.volume == 100 );

  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_vlasov_position_map_release(pos_map);
  gkyl_array_release(qmem);
  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(poisson_tensor_conf);
  gkyl_array_release(pot_tot);
  gkyl_array_release(vel_flux_surf);
  gkyl_array_release(f_no_J);
  gkyl_array_release(rad);
}

#ifdef GKYL_HAVE_CUDA

void
test_cu_dg_vlasov()
{
  int cdim = 1, vdim = 1, pdim = 2;
  struct gkyl_basis cbasis, vbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_serendip(&vbasis, 1, 1);
  gkyl_cart_modal_serendip(&pbasis, 2, 1);

  int cells[] = {100, 12 };
  int ghost[] = {1, 0 };
  double lower[] = {0., -1. };
  double upper[] = {1., 1. };

  int velcells[] = {12};
  int velghost[] = { 0};
  double vellower[] = { -1.};
  double velupper[] = { 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid velGrid;
  struct gkyl_range velRange, velRange_ext;
  gkyl_rect_grid_init(&velGrid, vdim, vellower, velupper, velcells);
  gkyl_create_grid_ranges(&velGrid, velghost, &velRange_ext, &velRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize eqn
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;

    // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr1(true, vbasis.num_basis, velRange.volume);
  struct gkyl_array *gamma_inv = mkarr1(true, vbasis.num_basis, velRange.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velRange, &vbasis, inp_vmap, false, true);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confRange, &confRange_ext, &cbasis, inp_pmap, true);

  gkyl_dg_vlasov_calc_hamil(&velGrid, &vbasis, &velRange, 
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, true);

  // Sturcture pointers for input objects (but not used)
  int num_pt_indices[3] = { 1 , 6, 18 }; 
  struct gkyl_array *poisson_tensor_conf = mkarr1(true, cbasis.num_basis*num_pt_indices[vdim-1], confRange.volume );
  struct gkyl_array *pot_tot = mkarr1(true, cbasis.num_basis*4, confRange_ext.volume );
  struct gkyl_array *vel_flux_surf = mkarr1(true, pbasis.num_basis*vdim, phaseRange_ext.volume );
  struct gkyl_array *f_no_J = mkarr1(true, pbasis.num_basis, phaseRange_ext.volume); ;
  struct gkyl_array *rad = mkarr1(true, vdim*vbasis.num_basis, velRange.volume);
  struct gkyl_array *qmem = mkarr1(true, 8*cbasis.num_basis, confRange_ext.volume);

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .conf_range =  &confRange,
    .hamil_range = &velRange,
    .phase_range = &phaseRange,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .poisson_tensor_conf = poisson_tensor_conf,
    .hamil = hamil,
    .qmem = qmem, 
    .pot_tot = pot_tot, 
    .vel_flux_surf = vel_flux_surf, 
    .f_no_J = f_no_J, 
    .rad = rad, 
    .use_lo = false,
    .use_gpu = true,
  };  
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  TEST_CHECK( eqn->num_equations == 1 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  TEST_CHECK( vlasov->cdim == 1 );
  TEST_CHECK( vlasov->pdim == 2 );
  TEST_CHECK( vlasov->conf_range.volume == 100 );

  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_vlasov_position_map_release(pos_map);
  gkyl_array_release(qmem);
  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(poisson_tensor_conf);
  gkyl_array_release(pot_tot);
  gkyl_array_release(vel_flux_surf);
  gkyl_array_release(f_no_J);
  gkyl_array_release(rad);
}

#endif

TEST_LIST = {
  { "dg_vlasov", test_dg_vlasov },
#ifdef GKYL_HAVE_CUDA
  { "cu_dg_vlasov", test_cu_dg_vlasov },
#endif  
  { NULL, NULL },
};
